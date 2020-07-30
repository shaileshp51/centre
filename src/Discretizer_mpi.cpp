/*
 * Descretizer_mpi.cpp
 *
 *  Created on: Oct 30, 2019
 *      Author: shailesh
 */

#include "Discretizer.h"

#ifdef USE_OMPMPI

#include <mpi.h>
#include "MpiCommTags.h"
#include "macrologger/macrologger.h"

using namespace std;

void Discretizer::run_dtype_mpi(BAT_t dtype,
		std::vector<ull_int> &shuffle_index, const int rank, const int numprocs,
		const int n_thread_perproc) {

	if (rank == MASTER_PROC) {

		Netcdf_BAT traj_inp;
		std::string fname = inputs.getControl().getInfilepath()
				+ inputs.getBats().getFnames()[0];

		NetcdfFile::NCTYPE fltype = traj_inp.GetNetcdfConventions(
				fname.c_str());
		bool is_accepted_traj_format =
				(fltype == NetcdfFile::NCTYPE::NC_CENTRETRAJ) ? true : false;
		if (!is_accepted_traj_format) {
			LOG_ERROR(
					"Rank[%d]>> Provided file format is not a valid file format",
					rank);
			exit(0);
		}
		traj_inp.NC_openRead(fname);
		traj_inp.setupFrameDim();
		traj_inp.setupTime();
		CoordinateInfo cInfo;
		traj_inp.setupRead(cInfo);

		nframes_trjs = traj_inp.Ncframe();
		if (nfrm4rl2int > nframes_trjs) {
			LOG_ERROR(
					"Rank[%d]>> Trajectory frames=%llu is less than mentioned in input %llu for discretization",
					rank, nframes_trjs, nfrm4rl2int);
			exit(0);
		}

		std::vector<float> data(nframes_tot_eff);

		ProgTaskSet *dscrt_dtyp;
		std::string int_traj_file;
		std::string content;
		u_int n_dtyptasks = 0;
		switch (dtype) {
		case BAT_t::BOND:
			dscrt_dtyp = &(progressstate.dscrt_b);
			int_traj_file.assign(
					inputs.getControl().getOutfilepath() + "bin_bonds.nc");
			content.assign("B");
			n_dtyptasks = traj_inp.Ncbond();
			break;
		case BAT_t::ANGLE:
			dscrt_dtyp = &(progressstate.dscrt_a);
			int_traj_file.assign(
					inputs.getControl().getOutfilepath() + "bin_angles.nc");
			content.assign("A");
			n_dtyptasks = traj_inp.Ncangle();
			break;
		case BAT_t::DIHEDRAL:
			dscrt_dtyp = &(progressstate.dscrt_d);
			int_traj_file.assign(
					inputs.getControl().getOutfilepath() + "bin_torsions.nc");
			content.assign("T");
			n_dtyptasks = traj_inp.Ncdihedral();
			break;
		default:
			LOG_ERROR("Rank[%d]>> Invalid BAT_t encountered", rank);
		}
		dscrt_dtyp->start = Clock::now();
		dscrt_dtyp->current = dscrt_dtyp->start;
		dscrt_dtyp->strt_read = dscrt_dtyp->curr_read = dscrt_dtyp->start;
		dscrt_dtyp->strt_comp = dscrt_dtyp->curr_comp = dscrt_dtyp->start;
		dscrt_dtyp->strt_comm = dscrt_dtyp->curr_comm = dscrt_dtyp->start;
		dscrt_dtyp->strt_write = dscrt_dtyp->curr_write = dscrt_dtyp->start;
		dscrt_dtyp->cstt = ExecState::RUNNING;
		if (progressstate.dscrt.cstt != ExecState::RUNNING) {
			progressstate.dscrt.cstt = ExecState::RUNNING;
			progressstate.dscrt.start = dscrt_dtyp->start;
			progressstate.dscrt.current = dscrt_dtyp->start;
		}

		Netcdf_TrjInt *ptr_int_dtype_traj = new Netcdf_TrjInt(int_traj_file,
				content, n_dtyptasks, 0, 1, (ull_int) nfrm4rl2int,
				(hbin_t) bin_schemes.size(), bin_schemes);
		ptr_int_dtype_traj->NC_create(content + " Integer Vectors");

		int chache_dims_per_proc = cache_dscrt_dims_per_cpu * n_thread_perproc;
		int chached_dims = chache_dims_per_proc * numprocs;

		MPI_Barrier(MPI_COMM_WORLD);

		LOG_DEBUG("Rank[%d]>> bin schemes master: %ld, , scheme: %d", rank, bin_schemes.size(), (int) bin_schemes[0]);

		for (int bl_strat_id = 0; bl_strat_id < (int)n_dtyptasks; bl_strat_id +=
				chached_dims) {
			int block_tasks =
					(bl_strat_id + chached_dims <= (int)n_dtyptasks) ?
							chached_dims : (n_dtyptasks - bl_strat_id);
			int curr_cache_per_proc = chache_dims_per_proc;
			int max_slave_rank = numprocs;
			if (block_tasks < chached_dims) {
				curr_cache_per_proc = std::ceil(
						(double) block_tasks / numprocs);
				max_slave_rank = std::ceil(
						(double) block_tasks / curr_cache_per_proc);
			}
			int bl_details[4];
			bl_details[0] = bl_strat_id;
			bl_details[1] = block_tasks;
			bl_details[2] = curr_cache_per_proc;
			bl_details[3] = max_slave_rank;

			// broadcast from master to all other
			MPI_Bcast(&bl_details, 4, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);

			LOG_DEBUG("Rank[%d]>> Bcast [bl_start=%d, bl_tasks=%d, curr_cache_perporoc=%d, max_rank=%d]",
					rank, bl_details[0], bl_details[1], bl_details[2], bl_details[3]);

			auto _startr = Clock::now();

			std::vector<std::vector<data_t>> dtypes_v(curr_cache_per_proc,
					std::vector<data_t>(nfrm4rl2int));
			std::vector<std::vector<hbin_t>> dtypes_int_v(curr_cache_per_proc,
					std::vector<hbin_t>(bin_schemes.size() * nfrm4rl2int));
			std::vector<std::vector<double>> dtypes_extrm_v(curr_cache_per_proc,
					std::vector<double>(4, 0.0));

			std::vector<hbin_t> nfoundmins(curr_cache_per_proc, 0);
			std::vector<std::vector<double>> min_poss(curr_cache_per_proc,
					std::vector<double>(inputs.getVmkde().getNmaxconf(), 0.0));

			std::vector<CoordBAT> vecs_CoordBAT(chache_dims_per_proc,
					CoordBAT(-1, bin_schemes.size(), nfrm4rl2int));

			auto _endr = Clock::now();
			auto _dur_rd = std::chrono::duration_cast<ClockResolution>(
					_endr - _startr).count();
			dscrt_dtyp->curr_read += ClockResolution(_dur_rd);
			dscrt_dtyp->curr_read += ClockResolution(_dur_rd);
			size_t n_reads_curr_block = 0;
			for (size_t process_idx = 1; process_idx < (size_t)max_slave_rank;
					++process_idx) {
				int trjreadstatus = 0;
				auto read_start_id = bl_strat_id
						+ (process_idx - 1) * curr_cache_per_proc;
				switch (dtype) {
				case BAT_t::BOND:
					trjreadstatus = traj_inp.readMultiBondFrames(read_start_id,
							curr_cache_per_proc, 0, nfrm4rl2int, dtypes_v);
					break;
				case BAT_t::ANGLE:
					trjreadstatus = traj_inp.readMultiAngleFrames(read_start_id,
							curr_cache_per_proc, 0, nfrm4rl2int, dtypes_v);
					break;
				case BAT_t::DIHEDRAL:
					trjreadstatus = traj_inp.readMultiDihedralFrames(
							read_start_id, curr_cache_per_proc, 0, nfrm4rl2int,
							dtypes_v);
					break;
					//case BAT_t::NONE:
				default:
					LOG_ERROR("Rank[%d]>> Unrecognised BAT_t type", rank);
					exit(0);
				}
				if (trjreadstatus != 0) {
					LOG_ERROR(
							"Rank[%d]>> Writing binned data for bonds-id range (%d, %d)",
							rank, bl_strat_id, bl_strat_id + block_tasks);
					exit(0);
				}

				for (size_t n_reads_curr_proc = 0;
						n_reads_curr_proc < (size_t)curr_cache_per_proc;
						++n_reads_curr_proc) {
					MPI_Send(dtypes_v[n_reads_curr_proc].data(), nfrm4rl2int,
					MPI_FLOAT, process_idx,
					DIM_REAL_TAG, MPI_COMM_WORLD);

					LOG_DEBUG("Rank[%d]>> for dim-id(%ld) Sent(%lld) floats to slave(%ld)",
							rank, bl_strat_id+n_reads_curr_proc, nfrm4rl2int,
							process_idx);

				}
				n_reads_curr_block += curr_cache_per_proc;
			}
			MPI_Barrier(MPI_COMM_WORLD);

			LOG_DEBUG("Rank[%d]>> data scattering to slaves is complete", rank);
			int n_task4master = block_tasks - n_reads_curr_block;
			LOG_DEBUG("Rank[%d]>> has %d block_tasks to do", rank, n_task4master);
			fflush(stdout);
			std::vector<std::vector<hbin_t>> dtypes_int_v_master;
			std::vector<std::vector<double>> dtypes_extrm_v_master;

			std::vector<hbin_t> nfoundmins_master;
			std::vector<std::vector<double>> min_poss_master;

			std::vector<CoordBAT> vecs_CoordBAT_master;
			if (n_task4master > 0) {
				for (int it4master = 0; it4master < n_task4master;
						++it4master) {
					dtypes_int_v_master.push_back(
							std::vector<hbin_t>(
									bin_schemes.size() * nfrm4rl2int));
					dtypes_extrm_v_master.push_back(
							std::vector<double>(4, 0.0));
					nfoundmins_master.push_back(0);
					min_poss_master.push_back(
							std::vector<double>(inputs.getVmkde().getNmaxconf(),
									0.0));
					vecs_CoordBAT_master.push_back(
							CoordBAT(-1, bin_schemes.size(), nfrm4rl2int));
				}
				// read for master process
				int trjreadstatmaster = 0;
				switch (dtype) {
				case BAT_t::BOND:
					trjreadstatmaster = traj_inp.readMultiBondFrames(
							bl_strat_id + n_reads_curr_block, n_task4master, 0,
							nfrm4rl2int, dtypes_v);
					break;
				case BAT_t::ANGLE:
					trjreadstatmaster = traj_inp.readMultiAngleFrames(
							bl_strat_id + n_reads_curr_block, n_task4master, 0,
							nfrm4rl2int, dtypes_v);
					break;
				case BAT_t::DIHEDRAL:
					trjreadstatmaster = traj_inp.readMultiDihedralFrames(
							bl_strat_id + n_reads_curr_block, n_task4master, 0,
							nfrm4rl2int, dtypes_v);
					break;
					//case BAT_t::NONE:
				default:
					LOG_ERROR("Rank[%d]>> Unrecognised BAT_t type", rank);
					exit(0);
				}
				if (trjreadstatmaster != 0) {
					LOG_ERROR(
							"Rank[%d]>> Writing binned data for bond-ids range (%d, %d)",
							rank, bl_strat_id, bl_strat_id + block_tasks);
					exit(0);
				}
				auto _endcomm1 = Clock::now();
				auto _durcomm1 = std::chrono::duration_cast<
						ClockResolution>(_endcomm1 - _endr).count();
				dscrt_dtyp->curr_comm += ClockResolution(_durcomm1);
				dscrt_dtyp->curr_comm += ClockResolution(_durcomm1);

#pragma omp parallel for
				for (int tidx = 0; tidx < n_task4master; ++tidx) {
					double min_d = Constants::DEFAULT_BAT_VALUE;
					double max_d = Constants::DEFAULT_BAT_VALUE;
					double avg_d = Constants::DEFAULT_BAT_VALUE;
					double phase_d = Constants::DEFAULT_BAT_VALUE;
					size_t did = bl_strat_id + n_reads_curr_block + tidx;
					if (inputs.getBats().getPdfmethod()
							== PDFMethod::HISTOGRAM) {
						int returnvalue = 0;
						switch (dtype) {
						case BAT_t::BOND:
							returnvalue = Real2Int::bndReal2Int(did,
									bin_schemes, dtypes_v[tidx], &min_d, &max_d,
									&avg_d, dtypes_int_v_master[tidx]);

							break;
						case BAT_t::ANGLE:
							returnvalue = Real2Int::angReal2Int(did,
									bin_schemes, dtypes_v[tidx], &min_d, &max_d,
									&avg_d, dtypes_int_v_master[tidx]);
							break;
						case BAT_t::DIHEDRAL:
							returnvalue = Real2Int::dihReal2Int(did,
									bin_schemes, inputs.getHist().getNbinRef(),
									inputs.getBats().isOptimizedih(),
									dtypes_v[tidx], &min_d, &max_d, &avg_d,
									&phase_d, dtypes_int_v_master[tidx]);
							break;
            case BAT_t::NONE:
              break;
						}
						if (returnvalue != 0) {
							LOG_ERROR(
									"Rank[%d]>> Converting Real2Int using HISTOGRAM for %s id(%ld)",
									rank, content.c_str(), did);
							exit(0);
						}
					} else if (inputs.getBats().getPdfmethod()
							== PDFMethod::vonMisesKDE) {
						if (Real2Int::getMinimaPositions(did, false,
								inputs.getVmkde().getNmaxconf(),
								inputs.getVmkde().getKappaValue(),
								inputs.getVmkde().getSdoNstep(),
								inputs.getVmkde().getSdoConvLimit(),
								inputs.getVmkde().getSdoNiter(), dtypes_v[tidx],
								min_poss_master[tidx], nfoundmins_master[tidx],
								dtypes_int_v_master[tidx])) {
							LOG_ERROR(
									"ERROR>> Converting Real2Int using vonMisesKDE for %s id(%ld)",
									content.c_str(), did);
							exit(0);
						}
					}

					if (inputs.getBats().isShuffleFrames()) {
						Real2Int::shuffleArray(shuffle_index,
								dtypes_int_v_master[tidx]);
					} else if (inputs.getBats().isShuffleDofs()) {
						unsigned seed =
								Clock::now().time_since_epoch().count();
						std::shuffle(shuffle_index.begin(), shuffle_index.end(),
								std::default_random_engine(seed));
						Real2Int::shuffleArray(shuffle_index,
								dtypes_int_v_master[tidx]);
					}

					dtypes_extrm_v_master[tidx][0] = min_d;
					dtypes_extrm_v_master[tidx][1] = max_d;
					dtypes_extrm_v_master[tidx][2] = avg_d;
					if (dtype == BAT_t::DIHEDRAL) {
						dtypes_extrm_v_master[tidx][3] = phase_d;
					} else {
						// phase in not valid in bonds/angle hence assigned 0.0
						dtypes_extrm_v_master[tidx][3] = 0.0;
					}

					(vecs_CoordBAT_master[tidx]).setId(did);
					(vecs_CoordBAT_master[tidx]).setCoords(
							dtypes_extrm_v_master[tidx][0],
							dtypes_extrm_v_master[tidx][1],
							dtypes_extrm_v_master[tidx][2],
							dtypes_extrm_v_master[tidx][3], 0,
							(hbin_t) bin_schemes.size(), 0, nfrm4rl2int - 1,
							dtypes_int_v_master[tidx]);
				}

				auto _endcomm2 = Clock::now();
				auto _durcomm2 = std::chrono::duration_cast<
						ClockResolution>(_endcomm2 - _endr).count();

				dscrt_dtyp->curr_comp += ClockResolution(_durcomm2);
				progressstate.dscrt.curr_comp += ClockResolution(
						_durcomm2);
			}

			MPI_Barrier(MPI_COMM_WORLD);

			LOG_DEBUG("Rank[%d]>> Waiting for slave results", rank);

			for (size_t process_idx = 1; process_idx < (size_t)max_slave_rank;
					++process_idx) {
				auto _startcomm = Clock::now();
				for (auto recv_idx = 0; recv_idx < curr_cache_per_proc;
						++recv_idx) {
					size_t did = bl_strat_id
							+ (process_idx - 1) * curr_cache_per_proc
							+ recv_idx;
					MPI_Recv(dtypes_extrm_v[recv_idx].data(), 4,
					MPI_DOUBLE, process_idx, DIM_EXTRM_TAG,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);

					LOG_DEBUG("Rank[%d]>> received extrema of did(%ld) from rank(%ld)", rank, did, process_idx);

					MPI_Recv(dtypes_int_v[recv_idx].data(),
							n_schemes * nfrm4rl2int,
							MPI_UNSIGNED_CHAR, process_idx,
							DIM_INT_TAG,
							MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					(vecs_CoordBAT[recv_idx]).setId(did);
					(vecs_CoordBAT[recv_idx]).setCoords(
							dtypes_extrm_v[recv_idx][0],
							dtypes_extrm_v[recv_idx][1],
							dtypes_extrm_v[recv_idx][2],
							dtypes_extrm_v[recv_idx][3], 0,
							(hbin_t) bin_schemes.size(), 0, nfrm4rl2int - 1,
							dtypes_int_v[recv_idx]);
				}
				auto _endcomm = Clock::now();
				auto _durcomm =
						std::chrono::duration_cast<ClockResolution>(
								_endcomm - _startcomm).count();

				dscrt_dtyp->curr_comm += ClockResolution(_durcomm);
				progressstate.dscrt.curr_comm += ClockResolution(
						_durcomm);
				auto trjintwriteresultslave = ptr_int_dtype_traj->writeRecords(
						vecs_CoordBAT, (ull_int) (curr_cache_per_proc));
				if (trjintwriteresultslave) {
					LOG_ERROR(
							"Rank[%d]>> error(%d): Writing binned data for %s range id(%d, %d)",
							rank, trjintwriteresultslave, content.c_str(),
							bl_strat_id, bl_strat_id + block_tasks);
					MPI_Abort(MPI_COMM_WORLD, 0);
				}
				auto _endwrt = Clock::now();
				auto _durwrt = std::chrono::duration_cast<
						ClockResolution>(_endwrt - _endcomm).count();

				dscrt_dtyp->curr_write += ClockResolution(_durwrt);
				progressstate.dscrt.curr_write += ClockResolution(
						_durwrt);
				dscrt_dtyp->done_tasks += curr_cache_per_proc;
				progressstate.dscrt.done_tasks += curr_cache_per_proc;
				dscrt_tasks_done += curr_cache_per_proc;
			}
			auto _startwrt = Clock::now();
			if (n_task4master > 0) {
				auto trjintwriteresult = ptr_int_dtype_traj->writeRecords(
						vecs_CoordBAT_master, (ull_int) (n_task4master));
				if (trjintwriteresult) {
					LOG_ERROR(
							"Rank[%d]>> error(%s): Writing binned data for range id(%d, %d)",
							rank, content.c_str(), bl_strat_id, bl_strat_id + block_tasks);
					MPI_Abort(MPI_COMM_WORLD, 0);
				}
				dscrt_dtyp->done_tasks += n_task4master;
				progressstate.dscrt.done_tasks += n_task4master;
				dscrt_tasks_done += n_task4master;
			}
			auto _endwrt = Clock::now();
			auto _durwrt = std::chrono::duration_cast<ClockResolution>(
					_endwrt - _startwrt).count();
			dscrt_dtyp->curr_write += ClockResolution(_durwrt);
			progressstate.dscrt.curr_write += ClockResolution(_durwrt);

			dscrt_dtyp->current = Clock::now();
			progressstate.dscrt.current = dscrt_dtyp->current;

			if (dscrt_tasks_done >= dscrt_tasks_freq) {
				std::ofstream info_strm(
						inputs.getControl().getOutfilepath()
								+ inputs.getControl().getInfofile());
				info_strm << progressstate.toString();
				info_strm.close();
				dscrt_tasks_done %= dscrt_tasks_freq;
			}

			MPI_Barrier(MPI_COMM_WORLD);
		}

		dscrt_dtyp->current = Clock::now();
		dscrt_dtyp->cstt = ExecState::COMPLETED;
		progressstate.dscrt.current = dscrt_dtyp->current;
		std::ofstream info_strm(
				inputs.getControl().getOutfilepath()
						+ inputs.getControl().getInfofile());
		info_strm << progressstate.toString();
		info_strm.close();
		traj_inp.NC_close();

	} else if (rank != MASTER_PROC) { // start slave block
		/*********************** SLAVE PROCESS *******************
		 * Double to int conversion
		 *********************************************************/
		std::string content;
		u_int n_dtyptasks = 0;
		switch (dtype) {
		case BAT_t::BOND:
			content.assign("B");
			n_dtyptasks = inputs.getBats().getNbond();
			break;
		case BAT_t::ANGLE:
			content.assign("A");
			n_dtyptasks = inputs.getBats().getNangle();
			break;
		case BAT_t::DIHEDRAL:
			content.assign("T");
			n_dtyptasks = inputs.getBats().getNdihed();
			break;
		default:
			LOG_ERROR("Rank[%d]>> Invalid BAT_t encountered", rank);
		}
		int chache_dims_per_proc = cache_dscrt_dims_per_cpu * numprocs;
		int chached_dims = chache_dims_per_proc * n_thread_perproc;

		int curr_cache_per_proc, max_slave_rank; // block_tasks, block_start_id;
		int lc_bl_det[4];
		// Process local variables for data receiving and processing
		MPI_Barrier(MPI_COMM_WORLD);

		for (int bl_strat_id = 0; bl_strat_id < (int)n_dtyptasks; bl_strat_id +=
				chached_dims) {
			MPI_Bcast(&lc_bl_det, 4, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			//block_start_id = lc_bl_det[0];
			//block_tasks = lc_bl_det[1];
			curr_cache_per_proc = lc_bl_det[2];
			max_slave_rank = lc_bl_det[3];

			LOG_DEBUG("Rank[%d]>> received [bl_start=%d, bl_tasks=%d, curr_cache_perporoc=%d, max_rank=%d]",
					rank, lc_bl_det[0], lc_bl_det[1], lc_bl_det[2], lc_bl_det[3]);

			if (rank < max_slave_rank) {
				std::vector<std::vector<data_t>> dtypslv_v(curr_cache_per_proc,
						std::vector<data_t>(nfrm4rl2int));
				std::vector<std::vector<hbin_t>> dtypslv_int_v(
						curr_cache_per_proc,
						std::vector<hbin_t>(n_schemes * nfrm4rl2int));
				std::vector<std::vector<double>> dim_extremas(
						curr_cache_per_proc, std::vector<double>(4));
				std::vector<hbin_t> nfoundmins(curr_cache_per_proc, 0);
				std::vector<std::vector<double>> min_poss(curr_cache_per_proc,
						std::vector<double>(inputs.getVmkde().getNmaxconf(),
								0.0));

				LOG_DEBUG("Rank[%d]>> Waiting to receive %d sets of real nframes_tot: %llu", rank, curr_cache_per_proc, nfrm4rl2int);

				for (int tid = 0; tid < curr_cache_per_proc; ++tid) {
					MPI_Recv(dtypslv_v[tid].data(), nfrm4rl2int, MPI_FLOAT,
					MASTER_PROC, DIM_REAL_TAG,
					MPI_COMM_WORLD, MPI_STATUS_IGNORE);

					LOG_DEBUG("Rank[%d]>> received real-array(%d) of size(%llu) of expected(%d)", rank, tid, nfrm4rl2int, curr_cache_per_proc);

				}
				MPI_Barrier(MPI_COMM_WORLD);

				LOG_DEBUG("Rank[%d]>> received (%d) data-sets successfully, Going to start computation",
						rank, curr_cache_per_proc);

#pragma omp parallel for
				for (int tidx = 0; tidx < curr_cache_per_proc; ++tidx) {
					double min_d = 0.0, max_d = 0.0, avg_d = 0.0, phase_d = 0.0;
					size_t did = bl_strat_id + (rank - 1) * curr_cache_per_proc
							+ tidx;
					LOG_DEBUG("Rank[%d]>> started converting did(%lu)", rank, did);

					if (inputs.getBats().getPdfmethod()
							== PDFMethod::HISTOGRAM) {
						int returnvalue = 0;
						switch (dtype) {
						case BAT_t::BOND:
							returnvalue = Real2Int::bndReal2Int(did,
									bin_schemes, dtypslv_v[tidx], &min_d,
									&max_d, &avg_d, dtypslv_int_v[tidx]);

							break;
						case BAT_t::ANGLE:
							returnvalue = Real2Int::angReal2Int(did,
									bin_schemes, dtypslv_v[tidx], &min_d,
									&max_d, &avg_d, dtypslv_int_v[tidx]);
							break;
						case BAT_t::DIHEDRAL:
							returnvalue = Real2Int::dihReal2Int(did,
									bin_schemes, inputs.getHist().getNbinRef(),
									inputs.getBats().isOptimizedih(),
									dtypslv_v[tidx], &min_d, &max_d, &avg_d,
									&phase_d, dtypslv_int_v[tidx]);
							break;
            case BAT_t::NONE:
              break;
						}
						if (returnvalue != 0) {
							LOG_ERROR(
									"Rank[%d]>> Converting Real2Int using HISTOGRAM for %s id(%lu)",
									rank, content.c_str(), did);
							MPI_Abort(MPI_COMM_WORLD, returnvalue);
						}
					} else if (inputs.getBats().getPdfmethod()
							== PDFMethod::vonMisesKDE) {
						if (Real2Int::getMinimaPositions(did, false,
								inputs.getVmkde().getNmaxconf(),
								inputs.getVmkde().getKappaValue(),
								inputs.getVmkde().getSdoNstep(),
								inputs.getVmkde().getSdoConvLimit(),
								inputs.getVmkde().getSdoNiter(),
								dtypslv_v[tidx], min_poss[tidx],
								nfoundmins[tidx], dtypslv_int_v[tidx])) {
							LOG_ERROR(
									"Rank[%d]>> Converting Real2Int using vonMisesKDE for %s id(%lu)",
									rank, content.c_str(), did);
							exit(0);
						}
					}
					if (inputs.getBats().isShuffleFrames()) {
						Real2Int::shuffleArray(shuffle_index,
								dtypslv_int_v[tidx]);
					} else if (inputs.getBats().isShuffleDofs()) {
						unsigned seed =
								Clock::now().time_since_epoch().count();
						std::shuffle(shuffle_index.begin(), shuffle_index.end(),
								std::default_random_engine(seed));
						Real2Int::shuffleArray(shuffle_index,
								dtypslv_int_v[tidx]);
					}
					dim_extremas[tidx][0] = min_d;
					dim_extremas[tidx][1] = max_d;
					dim_extremas[tidx][2] = avg_d;
					dim_extremas[tidx][3] = phase_d;
					LOG_DEBUG("Rank[%d]>> Finished converting did(%lu)", rank, did);
				}

				MPI_Barrier(MPI_COMM_WORLD);
				LOG_DEBUG("Rank[%d]>> real2int computation completed, waiting to send back to master", rank);

				for (int tidx = 0; tidx < curr_cache_per_proc; ++tidx) {
					MPI_Send(dim_extremas[tidx].data(), 4, MPI_DOUBLE,
					MASTER_PROC, DIM_EXTRM_TAG, MPI_COMM_WORLD);
					LOG_DEBUG("Rank[%d]>> sent dim-exterma 4-doubles to master for did(%d)",
							rank, bl_strat_id + (rank - 1) * curr_cache_per_proc
							+ tidx);

					MPI_Send(dtypslv_int_v[tidx].data(),
							n_schemes * nfrm4rl2int,
							MPI_UNSIGNED_CHAR, MASTER_PROC, DIM_INT_TAG,
							MPI_COMM_WORLD);
					LOG_DEBUG("Rank[%d]>> Sending n_schemes * nframes_tot ints= %llu to master", rank, n_schemes * nframes_tot_eff);
				}
			} else {
				MPI_Barrier(MPI_COMM_WORLD); // Synchronize if this rank does'nt expect to receive data from master
				MPI_Barrier(MPI_COMM_WORLD); // Synchronize all when real2int completed on all mpi-pocesses
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	} // end slave block

}

void Discretizer::run_mpi(std::vector<ull_int> &shuffle_index, const int rank,
		const int numprocs, const int n_thread_perproc) {

	if (inputs.getControl().isDiscetize()) {

		Timer time_real2int;
		if (rank == MASTER_PROC) {
			progressstate.dscrt.start =
					Clock::now();

			time_real2int.start();
			mprintf("Trajectory Real to Integer conversion starts.\n");
			Netcdf_BAT trajIn;
			std::string fname = inputs.getControl().getInfilepath()
					+ inputs.getBats().getFnames()[0];

			NetcdfFile::NCTYPE fltype = trajIn.GetNetcdfConventions(
					fname.c_str());
			if (fltype != NetcdfFile::NCTYPE::NC_CENTRETRAJ) {
				LOG_ERROR("Input Netcdf File is not of a valid CENTRETRAJ");
				MPI_Abort(MPI_COMM_WORLD, 0);
			}
		}

		if ((inputs.getEntropy().getWorkSet() & BATSet::B1D)
				|| (inputs.getEntropy().getWorkSet() & BATSet::BB2D)) {
			run_dtype_mpi(BAT_t::BOND, shuffle_index, rank, numprocs,
					n_thread_perproc);
		}
		if ((inputs.getEntropy().getWorkSet() & BATSet::A1D)
				|| (inputs.getEntropy().getWorkSet() & BATSet::AA2D)) {
			run_dtype_mpi(BAT_t::ANGLE, shuffle_index, rank, numprocs,
					n_thread_perproc);
		}
		if ((inputs.getEntropy().getWorkSet() & BATSet::D1D)
				|| (inputs.getEntropy().getWorkSet() & BATSet::DD2D)) {
			run_dtype_mpi(BAT_t::DIHEDRAL, shuffle_index, rank, numprocs,
					n_thread_perproc);
		}

		if (rank == MASTER_PROC) {
			time_real2int.stop();
			mprintf("TIME: Discretization of trajectory time: %.4f seconds.\n",
					time_real2int.total());
		}
	}

}
#endif

