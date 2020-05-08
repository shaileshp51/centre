#include "EntropyCalculator.h"

#ifdef USE_OMPMPI
#include <mpi.h>
#include "MpiCommTags.h"
#include "macrologger/macrologger.h"

using namespace std;

void EntropyCalculator::run1d_mpi(BAT_t bat_type, const int rank,
		const int numprocs, const int n_thread_perproc) {
	if (rank == MASTER_PROC) {
		Timer time_d1;
		time_d1.start();

		ProgTaskSet *ptaskset;
		string str_dim_type;
		u_int n_dim_eff = 0;
		Netcdf_TrjInt *ptr_int_traj_d1d;
		string int_traj_filename_d1d;
		Netcdf_HistUtil *ptr_hist_d1d;
		string hist_filename_d1d;
		Netcdf_EntContri *ent_dim1d;
		string ent_filename_d1d;
		BATSet curr_batset = BATSet::NOTHING;

		switch (bat_type) {
		case BAT_t::BOND:
			ptaskset = &(progressstate.entc_b1d);
			str_dim_type.assign("B");
			n_dim_eff = n_bnd_eff;
			curr_batset = BATSet::B1D;

			int_traj_filename_d1d.assign("bin_bonds.nc");
			hist_filename_d1d.assign("hist_bnd-1d.nc");
			ent_filename_d1d.assign("entcontri_bnd-1d.nc");

			break;
		case BAT_t::ANGLE:
			ptaskset = &(progressstate.entc_a1d);
			str_dim_type.assign("A");
			n_dim_eff = n_ang_eff;
			curr_batset = BATSet::A1D;

			int_traj_filename_d1d.assign("bin_angles.nc");
			hist_filename_d1d.assign("hist_ang-1d.nc");
			ent_filename_d1d.assign("entcontri_ang-1d.nc");

			break;
		case BAT_t::DIHEDRAL:
			ptaskset = &(progressstate.entc_d1d);
			str_dim_type.assign("D");
			n_dim_eff = n_dih_eff;
			curr_batset = BATSet::D1D;

			int_traj_filename_d1d.assign("bin_torsions.nc");
			hist_filename_d1d.assign("hist_tor-1d.nc");
			ent_filename_d1d.assign("entcontri_tor-1d.nc");
			break;
		default:
			LOG_ERROR("Rank[%d]>> Invalid 1D type found", rank);
			exit(0);
		}

		vector<ull_int> dim_lens(1, n_dim_eff);
		ptr_int_traj_d1d = new Netcdf_TrjInt(
				inputs.getControl().getOutfilepath() + int_traj_filename_d1d,
				str_dim_type, n_dim_eff, startframe, strideframe,
				(ull_int) nframes_tot, n_schemes, bin_schemes);

		ptr_hist_d1d = new Netcdf_HistUtil(
				inputs.getControl().getOutfilepath() + hist_filename_d1d,
				str_dim_type, 1, 1, start_hist_step, nsteps, stride_hist_step,
				dim_lens, TensorType::FULL, n_schemes, bin_schemes);

		ent_dim1d = new Netcdf_EntContri(
				inputs.getControl().getOutfilepath() + ent_filename_d1d,
				str_dim_type, 1, 1, nsteps, dim_lens, TensorType::FULL,
				bin_schemes, nestimators);

		ptaskset->start = chrono::high_resolution_clock::now();
		ptaskset->current = ptaskset->start;
		ptaskset->strt_read = ptaskset->curr_read = ptaskset->start;
		ptaskset->strt_comp = ptaskset->curr_comp = ptaskset->start;
		ptaskset->strt_comm = ptaskset->curr_comm = ptaskset->start;
		ptaskset->strt_write = ptaskset->curr_write = ptaskset->start;
		ptaskset->cstt = ExecState::RUNNING;
		if (progressstate.entc.cstt != ExecState::RUNNING) {
			progressstate.entc.cstt = ExecState::RUNNING;
			progressstate.entc.start = ptaskset->start;
			progressstate.entc.current = ptaskset->start;
		}

		ull_int nfrm_eff_entropy = 0;

		if (isWritefreq && (frqWriteset & bat_type)) {
			ptr_hist_d1d->NC_create(
					(str_dim_type + string("-1D histograms")).c_str());
		}

		ent_dim1d->NC_create(
				(str_dim_type + string("-1D Entropy contributions")).c_str());

		int chache_dims_per_proc = cache_entc_dims_per_cpu;
		int n_cpus = numprocs * n_thread_perproc;
		int chached_dims = chache_dims_per_proc * n_cpus;

		MPI_Barrier(MPI_COMM_WORLD);
		/*
		 * This loop divides the total `n_dim_eff` task into blocks of `chached_dims`.
		 * For each task in the cached_dims work of computing histogram/state-probability
		 * followed by entropy computation from the frequency for requested set of
		 * entropy estimators and number of steps and binning schemes is carried.
		 *
		 * For each block of tasks
		 * 1. All the relevant data are read for each MPI process on master and scattered to it.
		 * 2. respective computation is performed on all of the slaves processes
		 * 3. Master process reads data for its own share of tasks and does computation
		 * 4. Results from master process are written to respective files.
		 * 5. Results data are gathered from slave processes on master
		 * 6. Results from slave processes are written to respective files.
		 */
		for (u_int bl_strat_id = 0; bl_strat_id < n_dim_eff; bl_strat_id +=
				chached_dims) {
			const int block_tasks =
					(bl_strat_id + chached_dims <= n_dim_eff) ?
							chached_dims : (n_dim_eff - bl_strat_id);

			int curr_cache_per_cpu = chache_dims_per_proc;
			int max_slave_rank = numprocs;
			if (block_tasks < chached_dims) {
				curr_cache_per_cpu = ceil((double) block_tasks / n_cpus);
				max_slave_rank = ceil(
						(double) block_tasks
								/ (n_thread_perproc * curr_cache_per_cpu));
			}
			int bl_details[4];
			bl_details[0] = bl_strat_id;
			bl_details[1] = block_tasks;
			bl_details[2] = curr_cache_per_cpu;
			bl_details[3] = max_slave_rank;
			int curr_cache_per_proc = n_thread_perproc * curr_cache_per_cpu;

			MPI_Bcast(&bl_details, 4, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);

			LOG_DEBUG("Rank[%d]>> Broadcast [bl_start=%d, bl_tasks=%d, cache_per_cpu=%d, max_rank=%d]", rank, bl_details[0], bl_details[1], bl_details[2], bl_details[3]);

			vector<hbin_t*> d1ds_int_v(curr_cache_per_proc);
			vector<vector<double>> d1ds_extrm_v(curr_cache_per_proc,
					vector<double>(4, 0.0));
			vector<CoordBAT> crd(curr_cache_per_proc,
					CoordBAT(0, n_schemes, nframes_tot_eff));
			vector<vector<ull_int>> freq_obs_v(curr_cache_per_proc,
					vector<ull_int>(nsteps * bin_schemes_sum));
			vector<vector<double>> bin_mids_v(curr_cache_per_proc,
					vector<double>(bin_schemes_sum, 0.0));
			vector<vector<double>> entropy_v(curr_cache_per_proc,
					vector<double>(nestimators * nsteps * n_schemes, 0.0));

			vector<BinGroup> bingrp_v(curr_cache_per_proc,
					BinGroup(-1, ptr_hist_d1d->getStepsEff(),
							ptr_hist_d1d->getDimXBinSchemesSum()));

			vector<DimEntropy> dimentropy_v(curr_cache_per_proc,
					DimEntropy(-1, (hbin_t) nsteps, n_schemes, nestimators));

			/*
			 * DATA READING FOR SLAVE PROCESSES ON MASTER AND SCATTERING FROM MASTER TO REQUIRED
			 * NUMBER OF SLAVES.
			 */
			size_t n_reads_curr_block = 0;
			for (size_t process_idx = 1; process_idx < max_slave_rank;
					++process_idx) {
				/*
				 * Read descretized data needed from file for a slave process
				 */
				auto _startr = chrono::high_resolution_clock::now();
				for (int idx = 0; idx < curr_cache_per_proc; ++idx) {
					const u_int dim_id = inputs.getSubset().getBATTypeId(
							bat_type, bl_strat_id + n_reads_curr_block + idx);
					crd[idx].setId(dim_id);
					if (ptr_int_traj_d1d->readCoords(dim_id, 0,
							(hbin_t) n_schemes, crd[idx]) != 0) {
						LOG_ERROR(
								"Rank[%d]>> for rank[%d] reading %s-int data for id(%d)",
								rank, process_idx, str_dim_type.c_str(),
								dim_id);
					}
					crd[idx].getCoords(&d1ds_extrm_v[idx][0],
							&d1ds_extrm_v[idx][1], &d1ds_extrm_v[idx][2],
							&d1ds_extrm_v[idx][3], &n_schemes,
							&nfrm_eff_entropy, &d1ds_int_v[idx]);
					LOG_DEBUG("Rank[%d]>> for rank[%d], successfully read %s-int data for id(%d)", rank, process_idx, str_dim_type.c_str(),
							dim_id);
				}

				auto _endr = chrono::high_resolution_clock::now();
				auto _dur_rd = chrono::duration_cast<chrono::nanoseconds>(
						_endr - _startr).count();
				ptaskset->curr_read += chrono::nanoseconds(_dur_rd);
				progressstate.entc.curr_read += chrono::nanoseconds(_dur_rd);

				/*
				 * Send data for bin frequency calculation and entropy estimation to the slave process
				 */
				for (int idx = 0; idx < curr_cache_per_proc; ++idx) {
					const u_int dim_id = inputs.getSubset().getBATTypeId(
							bat_type, bl_strat_id + n_reads_curr_block + idx);
					MPI_Send(d1ds_extrm_v[idx].data(), 4, MPI_DOUBLE,
							process_idx, DIM_EXTRM_TAG, MPI_COMM_WORLD);
					MPI_Send(d1ds_int_v[idx], n_schemes * nframes_tot_eff,
					MPI_UNSIGNED_CHAR, process_idx,
					DIM_INT_TAG, MPI_COMM_WORLD);

					LOG_DEBUG("Rank[%d]>> => Rank[%d] %s id(%d) sent(%d) ints", rank, process_idx, str_dim_type.c_str(), dim_id, n_schemes * nframes_tot_eff);

				}
				auto _endcomm1 = chrono::high_resolution_clock::now();
				auto _durcomm1 = chrono::duration_cast<chrono::nanoseconds>(
						_endcomm1 - _endr).count();
				ptaskset->curr_comm += chrono::nanoseconds(_durcomm1);
				progressstate.entc.curr_comm += chrono::nanoseconds(_durcomm1);
				n_reads_curr_block += curr_cache_per_proc;
			}
			/*
			 * DATA SCATTERING COMPLETED, WAIT TILL ALL(MASTER AND SLAVES) THE MPI PROCESSES
			 * HIT THIS POINT.
			 */
			MPI_Barrier(MPI_COMM_WORLD);
			int n_task4master = block_tasks - n_reads_curr_block;
			LOG_DEBUG("Rank[%d]>> %s data scattering to slaves is completed, will do %d tasks itself", rank, str_dim_type.c_str(), n_task4master);

			/*
			 * READ DATA FROM FILE AND DO ENTROPY COMPUTATION FOR SET OF TASK TO BE DONE ON
			 * MASTER PROCESS.
			 */
			vector<hbin_t*> d1ds_int_v_master;
			vector<vector<double>> d1ds_extrm_v_master;
			vector<CoordBAT> crd_master;
			vector<vector<ull_int>> freq_obs_v_master;
			vector<vector<double>> bin_mids_v_master;
			vector<vector<double>> entropy_v_master;
			vector<BinGroup> bingrp_v_master;
			vector<DimEntropy> dimentropy_v_master;

			if (n_task4master > 0) {
				for (int it4master = 0; it4master < n_task4master;
						++it4master) {
					hbin_t *hbptr;
					d1ds_int_v_master.push_back(hbptr);
					d1ds_extrm_v_master.push_back(vector<double>(4, 0.0));
					crd_master.push_back(
							CoordBAT(0, n_schemes, nframes_tot_eff));
					freq_obs_v_master.push_back(
							vector<ull_int>(nsteps * bin_schemes_sum));
					bin_mids_v_master.push_back(
							vector<double>(bin_schemes_sum, 0.0));
					entropy_v_master.push_back(
							vector<double>(nestimators * nsteps * n_schemes,
									0.0));

					bingrp_v_master.push_back(
							BinGroup(-1, ptr_hist_d1d->getStepsEff(),
									ptr_hist_d1d->getDimXBinSchemesSum()));

					dimentropy_v_master.push_back(
							DimEntropy(-1, (hbin_t) nsteps, n_schemes,
									nestimators));
				}
				/*
				 * Reading data from file for master process itself.
				 */
				auto _startr = chrono::high_resolution_clock::now();
				for (int idx = 0; idx < n_task4master; ++idx) {
					const u_int dim_id = inputs.getSubset().getBATTypeId(
							bat_type, bl_strat_id + n_reads_curr_block + idx);
					crd_master[idx].setId(dim_id);
					if (ptr_int_traj_d1d->readCoords(dim_id, 0,
							(hbin_t) n_schemes, crd_master[idx]) != 0) {
						LOG_ERROR("Rank[%d]>> Reading %s-int data for id(%d)",
								rank, str_dim_type.c_str(), dim_id);
					}
					crd_master[idx].getCoords(&d1ds_extrm_v_master[idx][0],
							&d1ds_extrm_v_master[idx][1],
							&d1ds_extrm_v_master[idx][2],
							&d1ds_extrm_v_master[idx][3], &n_schemes,
							&nfrm_eff_entropy, &d1ds_int_v_master[idx]);
					LOG_DEBUG("Rank[%d]>> Successfully read %s-int data for id(%d)", rank, str_dim_type.c_str(),
							dim_id);
				}
				auto _endr = chrono::high_resolution_clock::now();
				auto _dur_rd = chrono::duration_cast<chrono::nanoseconds>(
						_endr - _startr).count();
				ptaskset->curr_read += chrono::nanoseconds(_dur_rd);
				progressstate.entc.curr_read += chrono::nanoseconds(_dur_rd);
				auto _startcomp = chrono::high_resolution_clock::now();
				/*
				 * Computing bin frequency and entropy for task on master process.
				 */
#pragma omp parallel for
				for (int idx = 0; idx < n_task4master; ++idx) {
					size_t did = bl_strat_id + n_reads_curr_block + idx;
					if (binData1D(nsteps, bin_schemes, nframes_tot_eff,
							d1ds_extrm_v_master[idx][0],
							d1ds_extrm_v_master[idx][1], d1ds_int_v_master[idx],
							freq_obs_v_master[idx], bin_mids_v_master[idx])
							!= 0) {
						LOG_ERROR("Rank[%d]>> binning %s data id(%d)", rank,
								str_dim_type.c_str(), did);
					} else {
						entropy1D(bat_type, did, inputs.getEstimators(), nsteps,
								step_size, d1ds_extrm_v_master[idx][0],
								d1ds_extrm_v_master[idx][1],
								inputs.getEntropy().isJacobian(), isKDE,
								bin_schemes, freq_obs_v_master[idx],
								entropy_v_master[idx]);

						const u_int dim_id = inputs.getSubset().getBATTypeId(
								bat_type, did);
						if (isWritefreq && (frqWriteset & curr_batset)) {
							bingrp_v_master[idx].setId(dim_id);
							bingrp_v_master[idx].setExtremes(
									d1ds_extrm_v_master[idx]);
							bingrp_v_master[idx].setBinMids(0, bin_schemes_sum,
									bin_mids_v_master[idx]);
							bingrp_v_master[idx].setBinFreqs(
									ptr_hist_d1d->getFirstStep(),
									ptr_hist_d1d->getStepStride(),
									ptr_hist_d1d->getStepsEff(), 0,
									bin_schemes_sum, freq_obs_v_master[idx]);
						}
						dimentropy_v_master[idx].setId(dim_id);
						dimentropy_v_master[idx].setDimContri(0, nestimators, 0,
								nsteps, 0, n_schemes, entropy_v_master[idx]);
					}
				}
				auto _endcomp = chrono::high_resolution_clock::now();
				auto _durcomp = chrono::duration_cast<chrono::nanoseconds>(
						_endcomp - _startcomp).count();
				ptaskset->curr_comp += chrono::nanoseconds(_durcomp);
				progressstate.entc.curr_comp += chrono::nanoseconds(_durcomp);

			}

			MPI_Barrier(MPI_COMM_WORLD);
			int n_gather_curr_block = 0;
			int n_prepared2write_curr_block = 0;
			for (size_t process_idx = 1; process_idx < max_slave_rank;
					++process_idx) {
				auto _startcomm = chrono::high_resolution_clock::now();
				for (int idx = 0; idx < curr_cache_per_proc; ++idx) {
					MPI_Recv(bin_mids_v[idx].data(), bin_schemes_sum,
					MPI_DOUBLE, process_idx,
					DIM_EXTRM_TAG, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);

					if (isWritefreq && (frqWriteset & curr_batset)) {
						MPI_Recv(freq_obs_v[idx].data(),
								nsteps * bin_schemes_sum,
								MPI_UNSIGNED_LONG_LONG, process_idx,
								DIM_FREQ_TAG, MPI_COMM_WORLD,
								MPI_STATUS_IGNORE);
					}
					MPI_Recv(entropy_v[idx].data(),
							nestimators * nsteps * n_schemes,
							MPI_DOUBLE, process_idx, DIM_ENTR_TAG,
							MPI_COMM_WORLD, MPI_STATUS_IGNORE);

					n_gather_curr_block++;
					LOG_DEBUG("Rank[%d]>> <= Rank[%d] receiving entropy contribution id(%d)", rank, process_idx, bl_strat_id + idx);
				}
				auto _endcomm = chrono::high_resolution_clock::now();
				auto _durcomm = chrono::duration_cast<chrono::nanoseconds>(
						_startcomm - _endcomm).count();
				ptaskset->curr_comm += chrono::nanoseconds(_durcomm);
				progressstate.entc.curr_comm += chrono::nanoseconds(_durcomm);
				for (size_t idx = 0; idx < curr_cache_per_proc; ++idx) {
					const u_int dim_id = inputs.getSubset().getBATTypeId(
							bat_type,
							bl_strat_id + n_prepared2write_curr_block);
					if (isWritefreq && (frqWriteset & BATSet::B1D)) {
						bingrp_v[idx].setId(dim_id);
						bingrp_v[idx].setExtremes(d1ds_extrm_v[idx]);
						bingrp_v[idx].setBinMids(0, bin_schemes_sum,
								bin_mids_v[idx]);
						bingrp_v[idx].setBinFreqs(ptr_hist_d1d->getFirstStep(),
								ptr_hist_d1d->getStepStride(),
								ptr_hist_d1d->getStepsEff(), 0, bin_schemes_sum,
								freq_obs_v[idx]);
					}
					dimentropy_v[idx].setId(dim_id);
					dimentropy_v[idx].setDimContri(0, nestimators, 0, nsteps, 0,
							n_schemes, entropy_v[idx]);

					n_prepared2write_curr_block++;

				}
				if (isWritefreq && (frqWriteset & curr_batset)) {
					if (ptr_hist_d1d->writeRecords(bingrp_v) != 0) {
						LOG_ERROR(
								"Rank[%d]>> writing binned data for %s id range[%d, %d]",
								rank, str_dim_type.c_str(), bl_strat_id,
								bl_strat_id + n_prepared2write_curr_block);
					}
				}
				ent_dim1d->writeRecords(dimentropy_v);
				auto _endwrt = chrono::high_resolution_clock::now();
				auto _durwrt = chrono::duration_cast<chrono::nanoseconds>(
						_endwrt - _endcomm).count();
				ptaskset->curr_write += chrono::nanoseconds(_durwrt);
				progressstate.entc.curr_write += chrono::nanoseconds(_durwrt);
			}
			if (n_task4master > 0) {
				/*
				 * Write computed results bin_frequency and entropy to respective files.
				 */
				if (isWritefreq && (frqWriteset & curr_batset)) {
					if (ptr_hist_d1d->writeRecords(bingrp_v_master) != 0) {
						LOG_ERROR(
								"Rank[%d]>> writing binned data for %s id range[%d, %d]",
								rank, str_dim_type.c_str(),
								bl_strat_id + n_prepared2write_curr_block,
								bl_strat_id + n_prepared2write_curr_block
										+ n_task4master);
					}
				}
				ent_dim1d->writeRecords(dimentropy_v_master);
			}

			ptaskset->current = chrono::high_resolution_clock::now();
			progressstate.entc.current = ptaskset->current;
			ptaskset->done_tasks += block_tasks;
			progressstate.entc.done_tasks += block_tasks;

			entc_tasks_done += block_tasks;
			if (entc_tasks_done >= entc_tasks_freq) {
				ofstream info_strm(
						inputs.getControl().getOutfilepath()
								+ inputs.getControl().getInfofile());
				info_strm << progressstate.toString();
				info_strm.close();
				entc_tasks_done %= entc_tasks_freq;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		ptaskset->current = chrono::high_resolution_clock::now();
		ptaskset->cstt = ExecState::COMPLETED;
		progressstate.entc.current = ptaskset->current;

		int entc_active = static_cast<int>(inputs.getEntropy().getWorkSet());
		// The fall-through in case statement assume the order of execution BOND,ANGLE,DIHEDRAL
		switch (entc_active) {
		case 1:
			if (bat_type == BAT_t::BOND) {
				progressstate.entc.cstt = ExecState::COMPLETED;
			}
			break;
		case 2:
		case 3:
			if (bat_type == BAT_t::ANGLE) {
				progressstate.entc.cstt = ExecState::COMPLETED;
			}
			break;
		case 4:
		case 5:
		case 6:
		case 7:
			if (bat_type == BAT_t::DIHEDRAL) {
				progressstate.entc.cstt = ExecState::COMPLETED;
			}
			break;
		default:

			break;
		}
		ofstream info_strm(
				inputs.getControl().getOutfilepath()
						+ inputs.getControl().getInfofile());
		info_strm << progressstate.toString();
		info_strm.close();

		time_d1.stop();
		mprintf(
				"TIME: Histogram bin frequency calculation for %s: %.4f seconds.\n",
				str_dim_type.c_str(), time_d1.total());

	} else {
		/*
		 * Here goes the code which is run on slave processes
		 *
		 **/

		string str_dim_type;
		u_int n_dim_eff = 0;
		BATSet curr_batset = BATSet::NOTHING;

		switch (bat_type) {
		case BAT_t::BOND:
			str_dim_type.assign("B");
			n_dim_eff = n_bnd_eff;
			curr_batset = BATSet::B1D;
			break;
		case BAT_t::ANGLE:
			str_dim_type.assign("A");
			n_dim_eff = n_ang_eff;
			curr_batset = BATSet::A1D;
			break;
		case BAT_t::DIHEDRAL:
			str_dim_type.assign("D");
			n_dim_eff = n_dih_eff;
			curr_batset = BATSet::D1D;
			break;
		default:
			cout << "Error:: Invalid 1D type found" << endl;
			exit(0);
		}

		int chache_dims_per_proc = cache_entc_dims_per_cpu;
		int n_cpus = numprocs * n_thread_perproc;
		int chached_dims = chache_dims_per_proc * n_cpus;

		int block_tasks = 0, block_start_id = 0, curr_cache_per_proc = 0,
				max_slave_rank = 0;
		int bl_details[4];
		// Process local variables for data receiving and processing
		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);
		for (u_int bl_strat_id = 0; bl_strat_id < n_dim_eff; bl_strat_id +=
				chached_dims) {
			MPI_Bcast(&bl_details, 4, MPI_INT, MASTER_PROC, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			block_start_id = bl_details[0];
			block_tasks = bl_details[1];
			curr_cache_per_proc = bl_details[2];
			max_slave_rank = bl_details[3];

			LOG_DEBUG("Rank[%d]>> Broadcast [bl_start=%d, bl_tasks=%d, cache_per_cpu=%d, max_rank=%d]", rank, bl_details[0], bl_details[1], bl_details[2], bl_details[3]);

			vector<vector<hbin_t>> d1ds_int_v(curr_cache_per_proc,
					vector<hbin_t>(n_schemes * nframes_tot_eff, 0));
			vector<vector<double>> d1ds_extrm_v(curr_cache_per_proc,
					vector<double>(4, 0.0));
			vector<vector<ull_int>> freq_obs_v(curr_cache_per_proc,
					vector<ull_int>(nsteps * bin_schemes_sum));
			vector<vector<double>> bin_mids_v(curr_cache_per_proc,
					vector<double>(bin_schemes_sum, 0.0));
			vector<vector<double>> entropy_v(curr_cache_per_proc,
					vector<double>(nestimators * nsteps * n_schemes, 0.0));

			/*
			 * DATA READING FOR SLAVE PROCESSES ON MASTER AND SCATTERING FROM MASTER TO REQUIRED
			 * NUMBER OF SLAVES.
			 */
			size_t n_reads_curr_block = 0;
			if (rank < max_slave_rank) {
				/*
				 * Receive data from master for bin frequency calculation and entropy estimation on the slave process
				 */
				for (int idx = 0; idx < curr_cache_per_proc; ++idx) {
					MPI_Recv(d1ds_extrm_v[idx].data(), 4, MPI_DOUBLE,
					MASTER_PROC, DIM_EXTRM_TAG, MPI_COMM_WORLD,
							MPI_STATUS_IGNORE);
					MPI_Recv(d1ds_int_v[idx].data(),
							n_schemes * nframes_tot_eff,
							MPI_UNSIGNED_CHAR, MASTER_PROC,
							DIM_INT_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					LOG_DEBUG("Rank[%d]>> <= Rank[%d] %s id(%d) received %d ints", rank, MASTER_PROC, str_dim_type.c_str(), inputs.getSubset().getBATTypeId(
									bat_type, bl_strat_id + n_reads_curr_block + idx), n_schemes * nframes_tot_eff);
				}
				n_reads_curr_block += curr_cache_per_proc;
			}
			/*
			 * RECEIVINGG SCATTERED DATA COMPLETED, WAIT TILL ALL(MASTER AND SLAVES) THE MPI PROCESSES
			 * HIT THIS POINT.
			 */
			MPI_Barrier(MPI_COMM_WORLD);
			LOG_DEBUG("Rank[%d]>> %s data collecting from master is complete", rank, str_dim_type.c_str());

			if (rank < max_slave_rank) {
#pragma omp parallel for
				for (int idx = 0; idx < curr_cache_per_proc; ++idx) {
					size_t did = bl_strat_id + (rank - 1) * curr_cache_per_proc
							+ idx;
					if (binData1D(nsteps, bin_schemes, nframes_tot_eff,
							d1ds_extrm_v[idx][0], d1ds_extrm_v[idx][1],
							d1ds_int_v[idx].data(), freq_obs_v[idx],
							bin_mids_v[idx]) != 0) {
						LOG_ERROR("Rank[%d]>> %s binning data id(%d)", rank,
								str_dim_type.c_str(), did);
					} else {
						entropy1D(bat_type, did, inputs.getEstimators(), nsteps,
								step_size, d1ds_extrm_v[idx][0],
								d1ds_extrm_v[idx][1],
								inputs.getEntropy().isJacobian(), isKDE,
								bin_schemes, freq_obs_v[idx], entropy_v[idx]);
					}

				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			/*
			 * SENDING BACK COMPUTED BIN-MIDPOINTS, BIN-FREQUENCY AND ENTROPY RESULTS TO MASTER PROCESS
			 */
			if (rank < max_slave_rank) {
				// auto _startcomm = chrono::high_resolution_clock::now();
				for (int idx = 0; idx < curr_cache_per_proc; ++idx) {
					MPI_Send(bin_mids_v[idx].data(), bin_schemes_sum,
					MPI_DOUBLE, MASTER_PROC,
					DIM_EXTRM_TAG, MPI_COMM_WORLD);

					if (isWritefreq && (frqWriteset & curr_batset)) {
						MPI_Send(freq_obs_v[idx].data(),
								nsteps * bin_schemes_sum,
								MPI_UNSIGNED_LONG_LONG, MASTER_PROC,
								DIM_FREQ_TAG, MPI_COMM_WORLD);
					}
					MPI_Send(entropy_v[idx].data(),
							nestimators * nsteps * n_schemes,
							MPI_DOUBLE, MASTER_PROC, DIM_ENTR_TAG,
							MPI_COMM_WORLD);

					LOG_DEBUG("Rank[%d]>> => Rank[%d] sending entropy contribution id(%d)", rank, MASTER_PROC, bl_strat_id + idx);
				}

			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
}

void EntropyCalculator::run2d_xx_mpi(BATSet dim_type, const int rank,
		const int numprocs, const int n_thread_perproc) {
	if (rank == MASTER_PROC) {
		Timer time_ba;
		time_ba.start();

		ProgTaskSet *ptaskset;
		string str_dim_type;
		u_int n_dim_eff = 0;
		u_int num_xx2d;

		Netcdf_TrjInt *ptr_int_traj_xx2d;
		string int_traj_filename_xx2d;
		Netcdf_HistUtil *ptr_hist_xx2d;
		string hist_filename_xx2d;
		Netcdf_EntContri *ptr_ent_xx2d;
		string ent_filename_xx2d;

		vector<u_int> dim_neigh_keys;
		vector<u_int> dimtypes_v;
		BAT_t dtype1st = BAT_t::NONE;
		switch (dim_type) {
		case BATSet::BB2D:
			dtype1st = BAT_t::BOND;
			ptaskset = &(progressstate.entc_b2d);
			str_dim_type.assign("B/B");
			n_dim_eff = n_bnd_eff;

			int_traj_filename_xx2d.assign("bin_bonds.nc");
			hist_filename_xx2d.assign("hist_bnd-2d.nc");
			ent_filename_xx2d.assign("entcontri_bnd-2d.nc");

			inputs.getNeighbors().bondKeys(dim_neigh_keys);
			for (size_t d1 = 0; d1 < dim_neigh_keys.size(); ++d1) {
				num_xx2d += inputs.getNeighbors().bondNeighSize(
						dim_neigh_keys[d1]);
			}
			copy(inputs.getSubset().getBonds().begin(),
					inputs.getSubset().getBonds().end(),
					back_inserter(dimtypes_v));

			break;
		case BATSet::AA2D:
			dtype1st = BAT_t::ANGLE;
			ptaskset = &(progressstate.entc_a2d);
			str_dim_type.assign("A/A");
			n_dim_eff = n_ang_eff;

			int_traj_filename_xx2d.assign("bin_angles.nc");
			hist_filename_xx2d.assign("hist_ang-2d.nc");
			ent_filename_xx2d.assign("entcontri_ang-2d.nc");

			inputs.getNeighbors().angleKeys(dim_neigh_keys);
			for (size_t d1 = 0; d1 < dim_neigh_keys.size(); ++d1) {
				num_xx2d += inputs.getNeighbors().angleNeighSize(
						dim_neigh_keys[d1]);
			}

			copy(inputs.getSubset().getAngles().begin(),
					inputs.getSubset().getAngles().end(),
					back_inserter(dimtypes_v));
			break;
		case BATSet::DD2D:
			dtype1st = BAT_t::DIHEDRAL;
			ptaskset = &(progressstate.entc_d2d);
			str_dim_type.assign("D/D");
			n_dim_eff = n_dih_eff;

			int_traj_filename_xx2d.assign("bin_torsions.nc");
			hist_filename_xx2d.assign("hist_tor-2d.nc");
			ent_filename_xx2d.assign("entcontri_tor-2d.nc");

			inputs.getNeighbors().torsionKeys(dim_neigh_keys);
			for (size_t d1 = 0; d1 < dim_neigh_keys.size(); ++d1) {
				num_xx2d += inputs.getNeighbors().torsionNeighSize(
						dim_neigh_keys[d1]);
			}
			copy(inputs.getSubset().getTorsions().begin(),
					inputs.getSubset().getTorsions().end(),
					back_inserter(dimtypes_v));
			break;
		default:
			LOG_ERROR("Rank[%d]>> Invalid xx2D type found", rank);
			exit(0);
		}

		vector<ull_int> dim_lens(2, n_dim_eff);
		ptr_int_traj_xx2d = new Netcdf_TrjInt(
				inputs.getControl().getOutfilepath() + int_traj_filename_xx2d,
				str_dim_type, n_dim_eff, startframe, strideframe,
				(ull_int) nframes_tot, n_schemes, bin_schemes);

		ptr_hist_xx2d = new Netcdf_HistUtil(
				inputs.getControl().getOutfilepath() + hist_filename_xx2d,
				str_dim_type, 2, 2, start_hist_step, nsteps, stride_hist_step,
				dim_lens, TensorType::FULL, n_schemes, bin_schemes);

		if (isWritefreq && (frqWriteset & dim_type)) {
			ptr_hist_xx2d->NC_create(str_dim_type + " 2-D histograms");
		}

		ptr_ent_xx2d = new Netcdf_EntContri(
				inputs.getControl().getOutfilepath() + ent_filename_xx2d,
				str_dim_type, 2, 2, nsteps, dim_lens, TensorType::UPPER,
				bin_schemes, nestimators);

		ptaskset->start = chrono::high_resolution_clock::now();
		ptaskset->current = ptaskset->start;
		ptaskset->strt_read = ptaskset->curr_read = ptaskset->start;
		ptaskset->strt_comp = ptaskset->curr_comp = ptaskset->start;
		ptaskset->strt_comm = ptaskset->curr_comm = ptaskset->start;
		ptaskset->strt_write = ptaskset->curr_write = ptaskset->start;
		ptaskset->cstt = ExecState::RUNNING;
		if (progressstate.entc.cstt != ExecState::RUNNING) {
			progressstate.entc.cstt = ExecState::RUNNING;
			progressstate.entc.start = ptaskset->start;
			progressstate.entc.current = ptaskset->start;
		}

		ull_int nfrm_eff_entropy = 0;
		/*********************** MASTER PROCESS *******************
		 * 2-D ENTROPY estimation for: BOND/BOND
		 *********************************************************/

		ptr_ent_xx2d->NC_create(str_dim_type + " 2-D Entropy contributions");

		vector<u_int> block_boundry(5);
		int block_count = 0, block_tasks = 0;
		bool block_ready = false;

		const u_int num_dim_neigh_keys = dim_neigh_keys.size();

		int chache_dims_per_proc = cache_entc_dims_per_cpu * n_thread_perproc
				/ 2;
		int n_cpus = numprocs * n_thread_perproc;
		int chached_dims = chache_dims_per_proc * numprocs;

		MPI_Barrier(MPI_COMM_WORLD);
		/*
		 * This loop divides the total `n_dim_eff` task into blocks of `chached_dims`.
		 * For each task in the cached_dims work of computing histogram/state-probability
		 * followed by entropy computation from the frequency for requested set of
		 * entropy estimators and number of steps and binning schemes is carried.
		 *
		 * For each block of tasks
		 */
		for (u_int bi = 0; bi < num_dim_neigh_keys; bi +=
				chache_dims_per_proc) {
			const int bl_rows =
					(bi + chache_dims_per_proc <= num_dim_neigh_keys) ?
							chache_dims_per_proc : (num_dim_neigh_keys - bi);

			int row_cache_info[2];
			row_cache_info[0] = bi;
			row_cache_info[1] = bl_rows;

			LOG_DEBUG("Rank[%d]>> Bcast row-cache-info [bl_start=%d, bl_rows=%d] for %s",
					rank, row_cache_info[0], row_cache_info[1], str_dim_type.c_str());

			vector<hbin_t*> dtyps1_int_v(bl_rows);
			vector<vector<double>> dtyps1_extrm_v(bl_rows,
					vector<double>(4, 0.0));
			vector<CoordBAT> crd_rows(bl_rows,
					CoordBAT(0, n_schemes, nframes_tot_eff));
			std: map<u_int, u_int> row_id2index;
			vector<u_int> vec_row_id2index(2 * bl_rows);
			// Fill Cache with rows
			auto _startr = chrono::high_resolution_clock::now();
			for (auto rno = 0; rno < bl_rows; ++rno) {
				const u_int dtyp_id1 = dim_neigh_keys[bi + rno];
				row_id2index[dtyp_id1] = rno;
				vec_row_id2index[2 * rno] = dtyp_id1;
				vec_row_id2index[2 * rno + 1] = rno;
				crd_rows[rno].setId(dtyp_id1);
				if (ptr_int_traj_xx2d->readCoords(dtyp_id1, 0, n_schemes,
						crd_rows[rno]) != 0) {
					LOG_ERROR("Rank[%d]>> reading row-ints of id(%d) for %s",
							rank, dtyp_id1, str_dim_type.c_str());
				}
				crd_rows[rno].getCoords(&dtyps1_extrm_v[rno][0],
						&dtyps1_extrm_v[rno][1], &dtyps1_extrm_v[rno][2],
						&dtyps1_extrm_v[rno][3], &n_schemes, &nfrm_eff_entropy,
						&dtyps1_int_v[rno]);
			}
			auto _endr = chrono::high_resolution_clock::now();
			auto _dur_rd = chrono::duration_cast<chrono::nanoseconds>(
					_endr - _startr).count();
			ptaskset->curr_read += chrono::nanoseconds(_dur_rd);
			progressstate.entc.curr_read += chrono::nanoseconds(_dur_rd);

			/*
			 * Send data for bin frequency calculation and entropy estimation to the slave process
			 */
			for (int idx = 0; idx < bl_rows; ++idx) {
				const u_int dim_id = dim_neigh_keys[bi + idx];
				MPI_Bcast(dtyps1_extrm_v[idx].data(), 4, MPI_DOUBLE,
				MASTER_PROC, MPI_COMM_WORLD);
				MPI_Bcast(dtyps1_int_v[idx], n_schemes * nframes_tot_eff,
				MPI_UNSIGNED_CHAR, MASTER_PROC, MPI_COMM_WORLD);

				LOG_DEBUG("Rank[%d]>> Bcast %d ints of row id(%d) for %s", rank, n_schemes * nframes_tot_eff, dim_id, str_dim_type.c_str());
			}
			auto _endcomm1 = chrono::high_resolution_clock::now();
			auto _durcomm1 = chrono::duration_cast<chrono::nanoseconds>(
					_endcomm1 - _endr).count();
			ptaskset->curr_comm += chrono::nanoseconds(_durcomm1);
			progressstate.entc.curr_comm += chrono::nanoseconds(_durcomm1);

			MPI_Barrier(MPI_COMM_WORLD);
			for (u_int bj = 0; bj < n_dim_eff; bj += chached_dims) {
				const int bl_cols =
						(bj + chached_dims <= n_dim_eff) ?
								chached_dims : (n_dim_eff - bj);

				int curr_col_cache_per_proc = chache_dims_per_proc;
				int max_slave_rank = numprocs;
				if (bl_cols < chached_dims) {
					curr_col_cache_per_proc = ceil((double) bl_cols / numprocs);
					max_slave_rank = ceil(
							(double) bl_cols / (curr_col_cache_per_proc));
				}
				int col_cache_info[4];
				col_cache_info[0] = bj;
				col_cache_info[1] = bl_cols;
				col_cache_info[2] = curr_col_cache_per_proc;
				col_cache_info[3] = max_slave_rank;

				MPI_Barrier(MPI_COMM_WORLD);

				int n_cols_sent2slaves = 0;
				for (auto process_idx = 1; process_idx < max_slave_rank;
						process_idx++, n_cols_sent2slaves +=
								curr_col_cache_per_proc) {
					// Calc
					u_int block_tasks = 0, col_idx = 0;
					map<u_int, u_int> read_cols_dims;
					map<pair<u_int, u_int>, u_int> id2index;

					for (auto ri = bi; ri < bi + bl_rows; ++ri) {
						const u_int dtyp_id1 = dim_neigh_keys[ri];
						const vector<u_int> neigh =
								inputs.getNeighbors().getDimTypeNeighs(dim_type,
										dtyp_id1);
						auto min_cj = dimtypes_v[bj
								+ (process_idx - 1) * curr_col_cache_per_proc];
						auto max_cj = dimtypes_v[bj
								+ (process_idx * curr_col_cache_per_proc) - 1];
						for (auto const cj : neigh) {
							if (min_cj <= cj && cj <= max_cj) {
								if (!read_cols_dims.count(cj)) {
									read_cols_dims[cj] = col_idx;
									++col_idx;
								}
								id2index[make_pair(dtyp_id1, cj)] = block_tasks;
								++block_tasks;
							}
						}
					}
					u_int proc_cols_block_tasks[2];
					proc_cols_block_tasks[0] = block_tasks;
					proc_cols_block_tasks[1] = read_cols_dims.size();
					vector<u_int> vec_read_cols_dims(2 * read_cols_dims.size());
					vector<u_int> vec_id2index(3 * block_tasks);
					MPI_Send(&proc_cols_block_tasks, 2, MPI_UNSIGNED,
							process_idx, XnX_CBLK_TAG, MPI_COMM_WORLD);

					LOG_DEBUG("Rank[%d]>> => Rank[%d] sent for %s [block_tasks=%d, read_cols_dims=%d]",
							rank, process_idx, str_dim_type.c_str(), proc_cols_block_tasks[0],
							proc_cols_block_tasks[1]);

					if (block_tasks > 0) {
						auto itmp = 0;
						for (auto it = read_cols_dims.begin();
								it != read_cols_dims.end(); ++it) {
							vec_read_cols_dims[itmp] = it->first;
							vec_read_cols_dims[itmp + 1] = it->second;
							itmp += 2;
						}
						itmp = 0;
						for (map<pair<u_int, u_int>, u_int>::iterator it =
								id2index.begin(); it != id2index.end(); ++it) {
							const pair<u_int, u_int> &id_pair = it->first;
							vec_id2index[itmp] = id_pair.first;
							vec_id2index[itmp + 1] = id_pair.second;
							vec_id2index[itmp + 2] = it->second;
							itmp += 3;
						}
						MPI_Send(vec_read_cols_dims.data(),
								vec_read_cols_dims.size(), MPI_UNSIGNED,
								process_idx, XnX_CIDX_TAG, MPI_COMM_WORLD);
						MPI_Send(vec_id2index.data(), vec_id2index.size(),
						MPI_UNSIGNED, process_idx, XnX_CMAP_TAG,
						MPI_COMM_WORLD);

						u_int n_cache_cols = read_cols_dims.size();
						vector<hbin_t*> dtyps2_int_v(n_cache_cols);
						vector<vector<double>> bnds2_extrm_v(n_cache_cols,
								vector<double>(4, 0.0));
						vector<CoordBAT> crd_cols(n_cache_cols,
								CoordBAT(0, n_schemes, nframes_tot_eff));

						// read and cache cols
						auto _startr2 = chrono::high_resolution_clock::now();
						for (map<u_int, u_int>::iterator it =
								read_cols_dims.begin();
								it != read_cols_dims.end(); ++it) {
							crd_cols[it->second].setId(it->first);
							if (ptr_int_traj_xx2d->readCoords(it->first, 0,
									n_schemes, crd_cols[it->second]) != 0) {
								LOG_ERROR(
										"Rank[%d]>> reading %d ints of column id(%d) for %s",
										rank, n_schemes * nframes_tot_eff,
										it->first, str_dim_type.c_str());
							}
							crd_cols[it->second].getCoords(
									&bnds2_extrm_v[it->second][0],
									&bnds2_extrm_v[it->second][1],
									&bnds2_extrm_v[it->second][2],
									&bnds2_extrm_v[it->second][3], &n_schemes,
									&nfrm_eff_entropy,
									&dtyps2_int_v[it->second]);
						}
						auto _endr2 = chrono::high_resolution_clock::now();
						auto _dur2_rd = chrono::duration_cast<
								chrono::nanoseconds>(_endr2 - _startr2).count();
						ptaskset->curr_read += chrono::nanoseconds(_dur2_rd);
						progressstate.entc.curr_read += chrono::nanoseconds(
								_dur2_rd);

						for (map<u_int, u_int>::iterator it =
								read_cols_dims.begin();
								it != read_cols_dims.end(); ++it) {
							MPI_Send(bnds2_extrm_v[it->second].data(), 4,
									MPI_DOUBLE, process_idx, XnX_EXTRM_TAG,
									MPI_COMM_WORLD);
							MPI_Send(dtyps2_int_v[it->second],
									n_schemes * nframes_tot_eff,
									MPI_UNSIGNED_CHAR, process_idx,
									XnX_INT_TAG, MPI_COMM_WORLD);

							LOG_DEBUG("Rank[%d]>> => Rank[%d] sent %d ints of column-id(%d) for for %s",
									rank, process_idx, n_schemes * nframes_tot_eff, it->first, str_dim_type.c_str());

						}

						auto _endcomm1 = chrono::high_resolution_clock::now();
						auto _durcomm1 = chrono::duration_cast<
								chrono::nanoseconds>(_endcomm1 - _endr).count();
						ptaskset->curr_comm += chrono::nanoseconds(_durcomm1);
						progressstate.entc.curr_comm += chrono::nanoseconds(
								_durcomm1);
					}
				}

				LOG_DEBUG("Rank[%d]>> scattering columns ints to slaves completed", rank);

				int n_cols4master = bl_cols - n_cols_sent2slaves;

				u_int master_tasks = 0, mastercol_idx = 0;
				map<u_int, u_int> masterread_cols_dims;
				map<pair<u_int, u_int>, u_int> masterid2index;

				for (auto ri = bi; ri < bi + bl_rows; ++ri) {
					const u_int dtyp_id1 = dim_neigh_keys[ri];
					const vector<u_int> neigh =
							inputs.getNeighbors().getDimTypeNeighs(dim_type,
									dtyp_id1);
					auto min_cj = dimtypes_v[bj + n_cols_sent2slaves];
					auto max_cj = dimtypes_v[bj + bl_cols - 1];
					for (auto const cj : neigh) {
						if (cj >= min_cj && cj <= max_cj) {
							if (!masterread_cols_dims.count(cj)) {
								masterread_cols_dims[cj] = mastercol_idx;
								++mastercol_idx;
							}
							masterid2index[make_pair(dtyp_id1, cj)] =
									master_tasks;
							++master_tasks;
						} else if (cj > max_cj) {
							break;
						}
					}
				}

				u_int n_mastercache_cols = masterread_cols_dims.size();
				LOG_DEBUG("Rank[%d]>> data for %s [block_tasks=%d, read_cols_dims=%d]",
						rank, str_dim_type.c_str(), master_tasks, n_mastercache_cols);

				vector<hbin_t*> dtyps2_int_v_master;
				vector<vector<double>> bnds2_extrm_v_master;
				vector<CoordBAT> crd_cols_master;

				vector<vector<ull_int>> freq2_bb_obs_v_master;
				vector<vector<double>> bin2_bb_mids_v_master;

				vector<vector<double>> entropy_bb_v_master;

				vector<u_int> dummy_ids_master(2, -1);
				vector<BinGroup> bingrp_v_master;
				vector<DimEntropy> dimentropy_v_master;

				vector<u_int> vec_masterid2index;

				if (master_tasks > 0) {
					vec_masterid2index.resize(3 * master_tasks, 0);
					auto ii = 0;
					for (map<pair<u_int, u_int>, u_int>::iterator it =
							masterid2index.begin(); it != masterid2index.end();
							++it) {
						const pair<u_int, u_int> &id_pair = it->first;
						vec_masterid2index[ii] = id_pair.first;
						vec_masterid2index[ii + 1] = id_pair.second;
						vec_masterid2index[ii + 2] = it->second;
						ii += 3;
					}
					for (int it4master = 0; it4master < master_tasks;
							++it4master) {
						hbin_t *hbptr;
						dtyps2_int_v_master.push_back(hbptr);
						bnds2_extrm_v_master.push_back(vector<double>(4, 0.0));
						crd_cols_master.push_back(
								CoordBAT(0, n_schemes, nframes_tot_eff));
						freq2_bb_obs_v_master.push_back(
								vector<ull_int>(nsteps * nbins_sum_sq));
						bin2_bb_mids_v_master.push_back(
								vector<double>(2 * bin_schemes_sum, 0.0));
						entropy_bb_v_master.push_back(
								vector<double>(nestimators * nsteps * n_schemes,
										0.0));

						bingrp_v_master.push_back(
								BinGroup(dummy_ids_master, (hbin_t) nsteps,
										2 * bin_schemes_sum, nbins_sum_sq));

						dimentropy_v_master.push_back(
								DimEntropy(dummy_ids_master, (hbin_t) nsteps,
										n_schemes, nestimators));
					}

					// read and cache cols
					auto _startr2 = chrono::high_resolution_clock::now();
					for (map<u_int, u_int>::iterator it =
							masterread_cols_dims.begin();
							it != masterread_cols_dims.end(); ++it) {
						crd_cols_master[it->second].setId(it->first);
						if (ptr_int_traj_xx2d->readCoords(it->first, 0,
								n_schemes, crd_cols_master[it->second]) != 0) {
							LOG_ERROR(
									"Rank[%d]>> reading %d ints of column id(%d) for %s",
									rank, n_schemes * nframes_tot_eff,
									it->first, str_dim_type.c_str());

						}
						crd_cols_master[it->second].getCoords(
								&bnds2_extrm_v_master[it->second][0],
								&bnds2_extrm_v_master[it->second][1],
								&bnds2_extrm_v_master[it->second][2],
								&bnds2_extrm_v_master[it->second][3],
								&n_schemes, &nfrm_eff_entropy,
								&dtyps2_int_v_master[it->second]);
					}

					LOG_DEBUG("Rank[%d]>> Reading columns for master's tasks completed for %s", rank, str_dim_type.c_str());

					auto _endr2 = chrono::high_resolution_clock::now();
					auto _dur2_rd = chrono::duration_cast<chrono::nanoseconds>(
							_endr2 - _startr2).count();
					ptaskset->curr_read += chrono::nanoseconds(_dur2_rd);
					progressstate.entc.curr_read += chrono::nanoseconds(
							_dur2_rd);

					auto _endcomm1 = chrono::high_resolution_clock::now();

#pragma omp parallel for
					for (int idx = 0; idx < master_tasks; ++idx) {
						auto rid = vec_masterid2index[3 * idx];
						auto cid = vec_masterid2index[3 * idx + 1];
						auto rno = row_id2index[rid];
						auto cno = masterread_cols_dims[cid];

						LOG_DEBUG("Rank[%d]>> computing entropy of %s [tasks=%d idx=%d rid=%d, rno=%d, cid=%d, cno=%d]",
								rank, str_dim_type.c_str(), master_tasks, idx, rid, rno, cid, cno);

						if (binData2D(nsteps, bin_schemes, nframes_tot_eff,
								dtyps1_extrm_v[rno][0], dtyps1_extrm_v[rno][1],
								bnds2_extrm_v_master[cno][0],
								bnds2_extrm_v_master[cno][1], dtyps1_int_v[rno],
								dtyps2_int_v_master[cno],
								freq2_bb_obs_v_master[idx],
								bin2_bb_mids_v_master[idx]) != 0) {
							LOG_ERROR("Rank[%d]>> binning data %s id(%d, %d)",
									str_dim_type.c_str(), rid, cid);
						}
						entropy2D(dtype1st, dtype1st, rid, cid,
								inputs.getEstimators(), nsteps, step_size,
								dtyps1_extrm_v[rno][0], dtyps1_extrm_v[rno][1],
								bnds2_extrm_v_master[cno][0],
								bnds2_extrm_v_master[cno][1],
								inputs.getEntropy().isJacobian(), isKDE,
								bin_schemes, freq2_bb_obs_v_master[idx],
								entropy_bb_v_master[idx]);

						vector<u_int> dim_ids(2);
						dim_ids[0] = rid;
						dim_ids[1] = cid;
						vector<double> bb_extrm(8, 0.0);
						for (int xx = 0; xx < 4; ++xx) {
							bb_extrm[xx] = dtyps1_extrm_v[rno][xx];
							bb_extrm[xx + 4] = bnds2_extrm_v_master[cno][xx];
						}
						if (isWritefreq && (frqWriteset & dim_type)) {
							bingrp_v_master[idx].setId(0, 2, dim_ids);
							bingrp_v_master[idx].setExtremes(bb_extrm);
							bingrp_v_master[idx].setBinMids(0,
									2 * bin_schemes_sum,
									bin2_bb_mids_v_master[idx]);
							bingrp_v_master[idx].setBinFreqs(
									ptr_hist_xx2d->getFirstStep(),
									ptr_hist_xx2d->getStepStride(),
									ptr_hist_xx2d->getStepsEff(), 0,
									nbins_sum_sq, freq2_bb_obs_v_master[idx]);
						}
						dimentropy_v_master[idx].setId(0, 2, dim_ids);
						dimentropy_v_master[idx].setDimContri(0, nestimators, 0,
								nsteps, 0, n_schemes, entropy_bb_v_master[idx]);

					}
					auto _endcomp = chrono::high_resolution_clock::now();
					auto _durcomp = chrono::duration_cast<chrono::nanoseconds>(
							_endcomp - _endcomm1).count();
					ptaskset->curr_comp += chrono::nanoseconds(_durcomp);
					progressstate.entc.curr_comp += chrono::nanoseconds(
							_durcomp);
				} LOG_DEBUG("Rank[%d]>> Computing entropy for master completed..", rank);
				// start of receiving entropy and frequency data from slaves and writing to files
				int n_gather_curr_block = 0;
				int n_prepared2write_curr_block = 0;
				for (size_t process_idx = 1; process_idx < max_slave_rank;
						++process_idx) {
					u_int blcoks_tasks_slave = 0;
					u_int slave_cols_block_tasks[2];
					MPI_Recv(&slave_cols_block_tasks, 2, MPI_UNSIGNED,
							process_idx, XnX_CBLK_TAG, MPI_COMM_WORLD,
							MPI_STATUS_IGNORE);
					blcoks_tasks_slave = slave_cols_block_tasks[0];

					LOG_DEBUG("Rank[%d]>> <= Rank[%d] received for %s [block_tasks=%d, read_cols_dims=%d]",
							rank, process_idx, str_dim_type.c_str(), slave_cols_block_tasks[0],
							slave_cols_block_tasks[1]);

					if (blcoks_tasks_slave > 0) {
						vector<u_int> vec_read_cols_dims(
								2 * slave_cols_block_tasks[1]);
						vector<u_int> vec_id2index(3 * blcoks_tasks_slave);
						vector<double> dd_extrema_block(8 * blcoks_tasks_slave,
								0.0);

						MPI_Recv(vec_read_cols_dims.data(),
								2 * slave_cols_block_tasks[1],
								MPI_UNSIGNED, process_idx, XnX_CIDX_TAG,
								MPI_COMM_WORLD,
								MPI_STATUS_IGNORE);

						MPI_Recv(vec_id2index.data(), 3 * blcoks_tasks_slave,
						MPI_UNSIGNED, process_idx, XnX_CMAP_TAG,
						MPI_COMM_WORLD, MPI_STATUS_IGNORE);

						LOG_DEBUG("Rrank[%d]>> <= Rank[%d] received block_task=%d row/col index-id vectors for %s",
								rank, process_idx, blcoks_tasks_slave, str_dim_type.c_str());

						MPI_Recv(dd_extrema_block.data(),
								8 * blcoks_tasks_slave,
								MPI_DOUBLE, process_idx, XnX_EXTRM_TAG,
								MPI_COMM_WORLD,
								MPI_STATUS_IGNORE);

						LOG_DEBUG("Rank[%d]>> <= Rank[%d] received block_task=%d row/col index-id extremas=%d vectors for %s",
								rank, process_idx, blcoks_tasks_slave, 8*blcoks_tasks_slave, str_dim_type.c_str());

						vector<vector<ull_int>> freq2_bb_obs_v(
								blcoks_tasks_slave,
								vector<ull_int>(nsteps * nbins_sum_sq));
						vector<vector<double>> bin2_bb_mids_v(
								blcoks_tasks_slave,
								vector<double>(2 * bin_schemes_sum, 0.0));
						vector<vector<double>> entropy_bb_v(blcoks_tasks_slave,
								vector<double>(nestimators * nsteps * n_schemes,
										0.0));

						vector<double>(4, 0.0);
						vector<u_int> dim_ids(2);

						vector<BinGroup> bingrp_v(blcoks_tasks_slave,
								BinGroup(dummy_ids_master, (hbin_t) nsteps,
										2 * bin_schemes_sum, nbins_sum_sq));

						vector<DimEntropy> dimentropy_v(blcoks_tasks_slave,
								DimEntropy(dummy_ids_master, (hbin_t) nsteps,
										n_schemes, nestimators));

						auto _startcomm = chrono::high_resolution_clock::now();
						for (int idx = 0; idx < blcoks_tasks_slave; ++idx) {
							MPI_Recv(bin2_bb_mids_v[idx].data(),
									2 * bin_schemes_sum,
									MPI_DOUBLE, process_idx, XnX_EXTRM_TAG,
									MPI_COMM_WORLD,
									MPI_STATUS_IGNORE);

							if (isWritefreq && (frqWriteset & dim_type)) {
								MPI_Recv(freq2_bb_obs_v[idx].data(),
										nsteps * nbins_sum_sq,
										MPI_UNSIGNED_LONG_LONG, process_idx,
										XnX_FREQ_TAG, MPI_COMM_WORLD,
										MPI_STATUS_IGNORE);
							}
							MPI_Recv(entropy_bb_v[idx].data(),
									nestimators * nsteps * n_schemes,
									MPI_DOUBLE, process_idx, XnX_ENTR_TAG,
									MPI_COMM_WORLD, MPI_STATUS_IGNORE);

							n_gather_curr_block++;
							vector<u_int> dim_ids(2);
							dim_ids[0] = vec_id2index[3 * idx];
							dim_ids[1] = vec_id2index[3 * idx + 1];

							vector<double> bb_extrm(8, 0.0);
							for (int xx = 0; xx < 8; ++xx) {
								bb_extrm[xx] = dd_extrema_block[8 * idx + xx];
							}
							if (isWritefreq && (frqWriteset & dim_type)) {
								bingrp_v[idx].setId(0, 2, dim_ids);
								bingrp_v[idx].setExtremes(bb_extrm);
								bingrp_v[idx].setBinMids(0, 2 * bin_schemes_sum,
										bin2_bb_mids_v[idx]);
								bingrp_v[idx].setBinFreqs(
										ptr_hist_xx2d->getFirstStep(),
										ptr_hist_xx2d->getStepStride(),
										ptr_hist_xx2d->getStepsEff(), 0,
										nbins_sum_sq, freq2_bb_obs_v[idx]);
							}
							dimentropy_v[idx].setId(0, 2, dim_ids);
							dimentropy_v[idx].setDimContri(0, nestimators, 0,
									nsteps, 0, n_schemes, entropy_bb_v[idx]);
							LOG_DEBUG("Rank[%d]>> <= Rank[%d] received entropy contribution of [row-id=%d, col-id=%d] for %s",
									rank, process_idx, dim_ids[0], dim_ids[1], str_dim_type.c_str());
						}
						if (isWritefreq && (frqWriteset & dim_type)) {
							if (ptr_hist_xx2d->writeRecords(bingrp_v) != 0) {
								LOG_ERROR(
										"Rank[%d]>> writing histogram of block_tasks=%d for %s from Rank[%d]",
										rank, block_tasks, str_dim_type.c_str(),
										process_idx);
							}
						}
						ptr_ent_xx2d->writeRecords(dimentropy_v);
						auto _endcomm = chrono::high_resolution_clock::now();
						auto _durcomm =
								chrono::duration_cast<chrono::nanoseconds>(
										_startcomm - _endcomm).count();
						ptaskset->curr_comm += chrono::nanoseconds(_durcomm);
						progressstate.entc.curr_comm += chrono::nanoseconds(
								_durcomm);
					}

					ptaskset->done_tasks += blcoks_tasks_slave;
					progressstate.entc.done_tasks += blcoks_tasks_slave;
					entc_tasks_done += blcoks_tasks_slave;
				}
				// End of receiving entropy and frequency data from slaves and writing to files

				// Start writing entropy and frequency data for master_tasks
				if (master_tasks > 0) {
					auto _strtwrt = chrono::high_resolution_clock::now();
					if (isWritefreq && (frqWriteset & BATSet::BB2D)) {
						if (ptr_hist_xx2d->writeRecords(bingrp_v_master) != 0) {
							LOG_ERROR(
									"Rank[%d]>> writing histogram of block_tasks=%d for %s",
									rank, block_tasks, str_dim_type.c_str());
						}
					}
					ptr_ent_xx2d->writeRecords(dimentropy_v_master);

					auto _endwrt = chrono::high_resolution_clock::now();
					auto _durwrt = chrono::duration_cast<chrono::nanoseconds>(
							_endwrt - _strtwrt).count();
					ptaskset->curr_write += chrono::nanoseconds(_durwrt);
					progressstate.entc.curr_write += chrono::nanoseconds(
							_durwrt);
					ptaskset->current = chrono::high_resolution_clock::now();
					progressstate.entc.current = ptaskset->current;
					ptaskset->done_tasks += master_tasks;
					progressstate.entc.done_tasks += master_tasks;

					entc_tasks_done += master_tasks;
				}

				if (entc_tasks_done >= entc_tasks_freq) {
					ofstream info_strm(
							inputs.getControl().getOutfilepath()
									+ inputs.getControl().getInfofile());
					info_strm << progressstate.toString();
					info_strm.close();
					entc_tasks_done %= entc_tasks_freq;
				}
			}
		} // dih-dih block calc for d1 ends

		ptaskset->current = chrono::high_resolution_clock::now();
		ptaskset->cstt = ExecState::COMPLETED;
		progressstate.entc.current = ptaskset->current;
		int entc_workset_intval =
				static_cast<int>(inputs.getEntropy().getWorkSet());
		int entc_aa2d_intval = static_cast<int>(BATSet::AA2D);
		int entc_dd2d_intval = static_cast<int>(BATSet::DD2D);
		int entc_xx2d_intval = static_cast<int>(BATSet::XX2D);
		if ((entc_workset_intval < entc_aa2d_intval)
				&& (dim_type == BATSet::BB2D)) {
			progressstate.entc.cstt = ExecState::COMPLETED;
		} else if ((entc_workset_intval < entc_dd2d_intval)
				&& (dim_type == BATSet::AA2D)) {
			progressstate.entc.cstt = ExecState::COMPLETED;

		} else if ((entc_workset_intval <= entc_xx2d_intval)
				&& (dim_type == BATSet::DD2D)) {
			progressstate.entc.cstt = ExecState::COMPLETED;

		}
		ofstream info_strm(
				inputs.getControl().getOutfilepath()
						+ inputs.getControl().getInfofile());
		info_strm << progressstate.toString();
		info_strm.close();
		time_ba.stop();
		mprintf(
				"TIME: Histogram bin frequency calculation for %s-2D: %.4f seconds.\n",
				str_dim_type.c_str(), time_ba.total());

	} else if (rank != MASTER_PROC) {
		u_int n_dim_eff = 0;
		u_int num_xx2d;
		string str_dim_type;

		vector<u_int> dimTypeKeys;
		vector<u_int> dimtypes_v;
		BAT_t dtype1st = BAT_t::NONE;
		switch (dim_type) {
		case BATSet::BB2D:
			dtype1st = BAT_t::BOND;
			n_dim_eff = n_bnd_eff;
			str_dim_type.assign("B/B");
			inputs.getNeighbors().bondKeys(dimTypeKeys);
			for (size_t d1 = 0; d1 < dimTypeKeys.size(); ++d1) {
				num_xx2d += inputs.getNeighbors().bondNeighSize(
						dimTypeKeys[d1]);
			}
			copy(inputs.getSubset().getBonds().begin(),
					inputs.getSubset().getBonds().end(),
					back_inserter(dimtypes_v));

			break;
		case BATSet::AA2D:
			dtype1st = BAT_t::ANGLE;
			n_dim_eff = n_ang_eff;
			str_dim_type.assign("A/A");
			inputs.getNeighbors().angleKeys(dimTypeKeys);
			for (size_t d1 = 0; d1 < dimTypeKeys.size(); ++d1) {
				num_xx2d += inputs.getNeighbors().angleNeighSize(
						dimTypeKeys[d1]);
			}

			copy(inputs.getSubset().getAngles().begin(),
					inputs.getSubset().getAngles().end(),
					back_inserter(dimtypes_v));
			break;
		case BATSet::DD2D:
			dtype1st = BAT_t::DIHEDRAL;
			n_dim_eff = n_dih_eff;
			str_dim_type.assign("D/D");
			inputs.getNeighbors().torsionKeys(dimTypeKeys);
			for (size_t d1 = 0; d1 < dimTypeKeys.size(); ++d1) {
				num_xx2d += inputs.getNeighbors().torsionNeighSize(
						dimTypeKeys[d1]);
			}
			copy(inputs.getSubset().getTorsions().begin(),
					inputs.getSubset().getTorsions().end(),
					back_inserter(dimtypes_v));
			break;
		default:
			LOG_ERROR("Rank[%d]>> Invalid xx2D type found", rank);
			exit(0);
		}

		vector<ull_int> dim_lens(2, n_dim_eff);
		ull_int nfrm_eff_entropy = 0;
		/*********************** MASTER PROCESS *******************
		 * 2-D ENTROPY estimation for: BOND/ANGLE
		 *********************************************************/
		int chache_dims_per_proc = cache_entc_dims_per_cpu * n_thread_perproc
				/ 2;
		int n_cpus = numprocs * n_thread_perproc;
		int chached_dims = chache_dims_per_proc * numprocs;

		int block_count = 0, block_tasks = 0;
		bool block_ready = false;

		const u_int num_dimTypeKeys = dimTypeKeys.size();

		MPI_Barrier(MPI_COMM_WORLD);

		for (u_int bi = 0; bi < num_dimTypeKeys; bi += chache_dims_per_proc) {
			const int bl_rows =
					(bi + chache_dims_per_proc <= num_dimTypeKeys) ?
							chache_dims_per_proc : (num_dimTypeKeys - bi);

			int row_cache_info[2];
			row_cache_info[0] = bi;
			row_cache_info[1] = bl_rows;

			LOG_DEBUG("Rank[%d]>> Bcast row-cache-info [bl_start=%d, bl_rows=%d] for %s",
					rank, row_cache_info[0], row_cache_info[1], str_dim_type.c_str());

			vector<vector<hbin_t>> dtyps1_int_v(bl_rows,
					vector<hbin_t>(n_schemes * nframes_tot_eff));
			vector<vector<double>> dtyps1_extrm_v(bl_rows,
					vector<double>(4, 0.0));
			vector<CoordBAT> crd_rows(bl_rows,
					CoordBAT(0, n_schemes, nframes_tot_eff));

			map<u_int, u_int> row_id2index;
			vector<u_int> vec_row_id2index(2 * bl_rows);

			/*
			 * Send data for bin frequency calculation and entropy estimation to the slave process
			 */
			for (int idx = 0; idx < bl_rows; ++idx) {
				const u_int dim_id = dimTypeKeys[bi + idx];

				MPI_Bcast(dtyps1_extrm_v[idx].data(), 4, MPI_DOUBLE,
				MASTER_PROC, MPI_COMM_WORLD);
				MPI_Bcast(dtyps1_int_v[idx].data(), n_schemes * nframes_tot_eff,
				MPI_UNSIGNED_CHAR, MASTER_PROC, MPI_COMM_WORLD);

				row_id2index[dim_id] = idx;
				vec_row_id2index[2 * idx] = dim_id;
				vec_row_id2index[2 * idx + 1] = idx;

				LOG_DEBUG("Rank[%d]>> Bcast %d ints of row id(%d) for %s",
						rank, n_schemes * nframes_tot_eff, dim_id, str_dim_type.c_str());
			}

			MPI_Barrier(MPI_COMM_WORLD);
			for (u_int bj = 0; bj < n_dim_eff; bj += chached_dims) {
				const int bl_cols =
						(bj + chached_dims <= n_dim_eff) ?
								chached_dims : (n_dim_eff - bj);

				int curr_col_cache_per_proc = chache_dims_per_proc;
				int max_slave_rank = numprocs;
				if (bl_cols < chached_dims) {
					curr_col_cache_per_proc = ceil((double) bl_cols / numprocs);
					max_slave_rank = ceil(
							(double) bl_cols / (curr_col_cache_per_proc));
				}
				int col_cache_info[4];
				col_cache_info[0] = bj;
				col_cache_info[1] = bl_cols;
				col_cache_info[2] = curr_col_cache_per_proc;
				col_cache_info[3] = max_slave_rank;

				MPI_Barrier(MPI_COMM_WORLD);
				// Calc
				u_int block_tasks = 0, col_idx = 0;
				u_int slave_cols_block_tasks[2];
				map<u_int, u_int> read_cols_dims;
				map<pair<u_int, u_int>, u_int> id2index;

				if (rank < max_slave_rank) {
					MPI_Recv(&slave_cols_block_tasks, 2, MPI_UNSIGNED,
					MASTER_PROC, XnX_CBLK_TAG, MPI_COMM_WORLD,
							MPI_STATUS_IGNORE);
					block_tasks = slave_cols_block_tasks[0];
					u_int n_cache_cols = slave_cols_block_tasks[1];

					LOG_DEBUG("Rank[%d]>> <= Rank[%d] received for %s [block_tasks=%d, read_cols_dims=%d]",
							rank, MASTER_PROC, str_dim_type.c_str(), slave_cols_block_tasks[0],
							slave_cols_block_tasks[1]);

					if (block_tasks > 0) {
						vector<u_int> vec_read_cols_dims(2 * n_cache_cols);
						vector<u_int> vec_id2index(3 * block_tasks);
						MPI_Recv(vec_read_cols_dims.data(), 2 * n_cache_cols,
						MPI_UNSIGNED, MASTER_PROC,
						XnX_CIDX_TAG, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);
						MPI_Recv(vec_id2index.data(), 3 * block_tasks,
						MPI_UNSIGNED, MASTER_PROC,
						XnX_CMAP_TAG, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);

						for (int ia = 0; ia < 2 * n_cache_cols; ia += 2) {
							read_cols_dims[vec_read_cols_dims[ia]] =
									vec_read_cols_dims[ia + 1];
						}
						for (int ia = 0; ia < 3 * block_tasks; ia += 3) {
							id2index[make_pair(vec_id2index[ia],
									vec_id2index[ia + 1])] =
									vec_id2index[ia + 2];
						}

						vector<vector<hbin_t>> dtyps2_int_v(n_cache_cols,
								vector<hbin_t>(n_schemes * nframes_tot_eff));
						vector<vector<double>> bnds2_extrm_v(n_cache_cols,
								vector<double>(4, 0.0));
						vector<CoordBAT> crd_cols(n_cache_cols,
								CoordBAT(0, n_schemes, nframes_tot_eff));

						// read and cache cols
						auto _startr2 = chrono::high_resolution_clock::now();
						for (map<u_int, u_int>::iterator it =
								read_cols_dims.begin();
								it != read_cols_dims.end(); ++it) {
							MPI_Recv(bnds2_extrm_v[it->second].data(), 4,
							MPI_DOUBLE, MASTER_PROC,
							XnX_EXTRM_TAG, MPI_COMM_WORLD,
							MPI_STATUS_IGNORE);
							MPI_Recv(dtyps2_int_v[it->second].data(),
									n_schemes * nframes_tot_eff,
									MPI_UNSIGNED_CHAR, MASTER_PROC,
									XnX_INT_TAG, MPI_COMM_WORLD,
									MPI_STATUS_IGNORE);
							LOG_DEBUG("Rank[%d]>> <= Rank[%d] received %d ints of column-id(%d) for for %s",
									rank, MASTER_PROC, n_schemes * nframes_tot_eff, it->first, str_dim_type.c_str());
						}

						vector<vector<ull_int>> freq2_bb_obs_v(block_tasks,
								vector<ull_int>(nsteps * nbins_sum_sq));
						vector<vector<double>> bin2_bb_mids_v(block_tasks,
								vector<double>(2 * bin_schemes_sum, 0.0));

						vector<vector<double>> entropy_bb_v(block_tasks,
								vector<double>(nestimators * nsteps * n_schemes,
										0.0));

						vector<u_int> dummy_ids(2, -1);
						vector<BinGroup> bingrp_v(block_tasks,
								BinGroup(dummy_ids, (hbin_t) nsteps,
										2 * bin_schemes_sum, nbins_sum_sq));
						vector<DimEntropy> dimentropy_v(block_tasks,
								DimEntropy(dummy_ids, (hbin_t) nsteps,
										n_schemes, nestimators));

						auto _endcomm1 = chrono::high_resolution_clock::now();
						vector<double> bb_extrm(8 * block_tasks, 0.0);

#pragma omp parallel for
						for (int idx = 0; idx < block_tasks; ++idx) {
							auto rid = vec_id2index[3 * idx];
							auto cid = vec_id2index[3 * idx + 1];
							auto rno = row_id2index[rid];
							auto cno = read_cols_dims[cid];

							LOG_DEBUG("Rank[%d]>> computing entropy of %s [tasks=%d idx=%d rid=%d, rno=%d, cid=%d, cno=%d]",
									rank, str_dim_type.c_str(), block_tasks, idx, rid, rno, cid, cno);

							if (binData2D(nsteps, bin_schemes, nframes_tot_eff,
									dtyps1_extrm_v[rno][0],
									dtyps1_extrm_v[rno][1],
									bnds2_extrm_v[cno][0],
									bnds2_extrm_v[cno][1],
									dtyps1_int_v[rno].data(),
									dtyps2_int_v[cno].data(),
									freq2_bb_obs_v[idx], bin2_bb_mids_v[idx])
									!= 0) {
								LOG_ERROR(
										"Rank[%d]>> binning data %s id(%d, %d)",
										str_dim_type.c_str(), rid, cid);
							}
							entropy2D(dtype1st, dtype1st, rid, cid,
									inputs.getEstimators(), nsteps, step_size,
									dtyps1_extrm_v[rno][0],
									dtyps1_extrm_v[rno][1],
									bnds2_extrm_v[cno][0],
									bnds2_extrm_v[cno][1],
									inputs.getEntropy().isJacobian(), isKDE,
									bin_schemes, freq2_bb_obs_v[idx],
									entropy_bb_v[idx]);

							vector<u_int> dim_ids(2);
							dim_ids[0] = rid;
							dim_ids[1] = cid;

							for (int xx = 0; xx < 4; ++xx) {
								bb_extrm[8 * idx + xx] =
										dtyps1_extrm_v[rno][xx];
								bb_extrm[(8 * idx) + xx + 4] =
										bnds2_extrm_v[cno][xx];
							}
						} // #end-omp parallel for

						LOG_DEBUG("Rank[%d]>> Computing entropy for master completed..", rank);

						MPI_Send(&slave_cols_block_tasks, 2, MPI_UNSIGNED,
						MASTER_PROC, XnX_CBLK_TAG, MPI_COMM_WORLD);
						MPI_Send(vec_read_cols_dims.data(),
								vec_read_cols_dims.size(), MPI_UNSIGNED,
								MASTER_PROC, XnX_CIDX_TAG, MPI_COMM_WORLD);
						MPI_Send(vec_id2index.data(), vec_id2index.size(),
						MPI_UNSIGNED, MASTER_PROC, XnX_CMAP_TAG,
						MPI_COMM_WORLD);

						LOG_DEBUG("Rrank[%d]>> => Rank[%d] sent block_task=%d row/col index-id vectors for %s",
								rank, MASTER_PROC, block_tasks, str_dim_type.c_str());

						MPI_Send(bb_extrm.data(), bb_extrm.size(),
						MPI_DOUBLE, MASTER_PROC, XnX_EXTRM_TAG,
						MPI_COMM_WORLD);

						LOG_DEBUG("Rank[%d]>> => Rank[%d] sent block_task=%d row/col index-id extremas=%d vectors for %s",
								rank, MASTER_PROC, block_tasks, 8*block_tasks, str_dim_type.c_str());

						for (int idx = 0; idx < block_tasks; ++idx) {
							MPI_Send(bin2_bb_mids_v[idx].data(),
									2 * bin_schemes_sum,
									MPI_DOUBLE, MASTER_PROC,
									XnX_EXTRM_TAG, MPI_COMM_WORLD);

							if (isWritefreq && (frqWriteset & dim_type)) {
								MPI_Send(freq2_bb_obs_v[idx].data(),
										nsteps * nbins_sum_sq,
										MPI_UNSIGNED_LONG_LONG, MASTER_PROC,
										XnX_FREQ_TAG, MPI_COMM_WORLD);
							}
							MPI_Send(entropy_bb_v[idx].data(),
									nestimators * nsteps * n_schemes,
									MPI_DOUBLE, MASTER_PROC, XnX_ENTR_TAG,
									MPI_COMM_WORLD);

							LOG_DEBUG("Rank[%d]>> => Rank[%d] sending entropy contribution of [row-id=%d, col-id=%d] for %s",
									rank, MASTER_PROC, vec_id2index[2*idx], vec_id2index[2*idx+1], str_dim_type.c_str());

						}
					} else if (block_tasks == 0) {
						MPI_Send(&slave_cols_block_tasks, 2, MPI_UNSIGNED,
						MASTER_PROC, XnX_CBLK_TAG, MPI_COMM_WORLD);
					}
				}
				//MPI_Barrier(MPI_COMM_WORLD);
			}
		} // dih-dih block calc for d1 ends
	}
}

void EntropyCalculator::run2d_xy_mpi(BATSet dim_type, const int rank,
		const int numprocs, const int n_thread_perproc) {
	if (rank == MASTER_PROC) {
		Timer time_xy;
		time_xy.start();

		ProgTaskSet *ptaskset;
		string str_dim_type;
		u_int n_dim_x_eff = 0;
		u_int n_dim_y_eff = 0;
		u_int num_xy2d;
		Netcdf_TrjInt *intxofXY2DTraj;
		Netcdf_TrjInt *intyofXY2DTraj;
		string intxofXY2DTrajFileName;
		string intyofXY2DTrajFileName;
		Netcdf_HistUtil *xy2DHist;
		string xy2DHistFileName;
		Netcdf_EntContri *ent_xy2d;
		string ent_xy2dFileName;

		vector<u_int> dim_neigh_keys;
		vector<u_int> dimtypes_v;
		BAT_t dtype1st = BAT_t::NONE;
		BAT_t dtype2nd = BAT_t::NONE;
		switch (dim_type) {
		case BATSet::BA2D:
			dtype1st = BAT_t::BOND;
			dtype2nd = BAT_t::ANGLE;
			ptaskset = &(progressstate.entc_ba2d);
			str_dim_type.assign("B/A");
			n_dim_x_eff = n_bnd_eff;
			n_dim_y_eff = n_ang_eff;

			intxofXY2DTrajFileName.assign("bin_bonds.nc");
			intyofXY2DTrajFileName.assign("bin_angles.nc");
			xy2DHistFileName.assign("hist_ba-2d.nc");
			ent_xy2dFileName.assign("entcontri_ba-2d.nc");

			inputs.getNeighbors().bacrossKeys(dim_neigh_keys);
			for (size_t d1 = 0; d1 < dim_neigh_keys.size(); ++d1) {
				num_xy2d += inputs.getNeighbors().baNeighSize(
						dim_neigh_keys[d1]);
			}

			copy(inputs.getSubset().getAngles().begin(),
					inputs.getSubset().getAngles().end(),
					back_inserter(dimtypes_v));
			break;
		case BATSet::BD2D:
			dtype1st = BAT_t::BOND;
			dtype2nd = BAT_t::DIHEDRAL;
			ptaskset = &(progressstate.entc_bd2d);
			str_dim_type.assign("B/D");
			n_dim_x_eff = n_bnd_eff;
			n_dim_y_eff = n_dih_eff;

			intxofXY2DTrajFileName.assign("bin_bonds.nc");
			intyofXY2DTrajFileName.assign("bin_torsions.nc");
			xy2DHistFileName.assign("hist_bd-2d.nc");
			ent_xy2dFileName.assign("entcontri_bd-2d.nc");

			inputs.getNeighbors().bdcrossKeys(dim_neigh_keys);
			for (size_t d1 = 0; d1 < dim_neigh_keys.size(); ++d1) {
				num_xy2d += inputs.getNeighbors().bdNeighSize(
						dim_neigh_keys[d1]);
			}

			copy(inputs.getSubset().getTorsions().begin(),
					inputs.getSubset().getTorsions().end(),
					back_inserter(dimtypes_v));
			break;
		case BATSet::AD2D:
			dtype1st = BAT_t::ANGLE;
			dtype2nd = BAT_t::DIHEDRAL;
			ptaskset = &(progressstate.entc_ad2d);
			str_dim_type.assign("A/D");
			n_dim_x_eff = n_ang_eff;
			n_dim_y_eff = n_dih_eff;

			intxofXY2DTrajFileName.assign("bin_angles.nc");
			intyofXY2DTrajFileName.assign("bin_torsions.nc");
			xy2DHistFileName.assign("hist_ad-2d.nc");
			ent_xy2dFileName.assign("entcontri_ad-2d.nc");

			inputs.getNeighbors().adcrossKeys(dim_neigh_keys);
			for (size_t d1 = 0; d1 < dim_neigh_keys.size(); ++d1) {
				num_xy2d += inputs.getNeighbors().adNeighSize(
						dim_neigh_keys[d1]);
			}
			copy(inputs.getSubset().getTorsions().begin(),
					inputs.getSubset().getTorsions().end(),
					back_inserter(dimtypes_v));
			break;
		default:
			LOG_ERROR("Rank[%d]>> Invalid xy2D type found", rank);
			exit(0);
		}

		vector<ull_int> dim_lens(2, 0);
		dim_lens[0] = n_dim_x_eff;
		dim_lens[1] = n_dim_y_eff;

		intxofXY2DTraj = new Netcdf_TrjInt(
				inputs.getControl().getOutfilepath() + intxofXY2DTrajFileName,
				str_dim_type, n_dim_x_eff, startframe, strideframe,
				(ull_int) nframes_tot, n_schemes, bin_schemes);

		intyofXY2DTraj = new Netcdf_TrjInt(
				inputs.getControl().getOutfilepath() + intyofXY2DTrajFileName,
				str_dim_type, n_dim_y_eff, startframe, strideframe,
				(ull_int) nframes_tot, n_schemes, bin_schemes);

		xy2DHist = new Netcdf_HistUtil(
				inputs.getControl().getOutfilepath() + xy2DHistFileName,
				str_dim_type, 2, 2, start_hist_step, nsteps, stride_hist_step,
				dim_lens, TensorType::FULL, n_schemes, bin_schemes);

		if (isWritefreq && (frqWriteset & dim_type)) {
			xy2DHist->NC_create(str_dim_type + " 2-D histograms");
		}

		ent_xy2d = new Netcdf_EntContri(
				inputs.getControl().getOutfilepath() + ent_xy2dFileName,
				str_dim_type, 2, 2, nsteps, dim_lens, TensorType::FULL,
				bin_schemes, nestimators);

		ent_xy2d->NC_create(str_dim_type + "-2D Entropy contributions");

		ptaskset->start = chrono::high_resolution_clock::now();
		ptaskset->current = ptaskset->start;
		ptaskset->strt_read = ptaskset->curr_read = ptaskset->start;
		ptaskset->strt_comp = ptaskset->curr_comp = ptaskset->start;
		ptaskset->strt_comm = ptaskset->curr_comm = ptaskset->start;
		ptaskset->strt_write = ptaskset->curr_write = ptaskset->start;
		ptaskset->cstt = ExecState::RUNNING;
		if (progressstate.entc.cstt != ExecState::RUNNING) {
			progressstate.entc.cstt = ExecState::RUNNING;
			progressstate.entc.start = ptaskset->start;
			progressstate.entc.current = ptaskset->start;
		}

		ull_int nfrm_eff_entropy = 0;

		/*********************** MASTER PROCESS *******************
		 * 2-D ENTROPY estimation for: DimTypeX/DimTypeY
		 *********************************************************/

		vector<u_int> block_boundry(5);
		int block_count = 0, block_tasks = 0;
		bool block_ready = false;

		const u_int num_dim_neigh_keys = dim_neigh_keys.size();

		int chache_dims_per_proc = cache_entc_dims_per_cpu * n_thread_perproc
				/ 2;
		int n_cpus = numprocs * n_thread_perproc;
		int chached_dims = chache_dims_per_proc * numprocs;

		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);
		/*
		 * This loop divides the total `n_dim_eff` task into blocks of `chached_dims`.
		 * For each task in the cached_dims work of computing histogram/state-probability
		 * followed by entropy computation from the frequency for requested set of
		 * entropy estimators and number of steps and binning schemes is carried.
		 *
		 * For each block of tasks
		 */
		for (u_int bi = 0; bi < num_dim_neigh_keys; bi +=
				chache_dims_per_proc) {
			const int bl_rows =
					(bi + chache_dims_per_proc <= num_dim_neigh_keys) ?
							chache_dims_per_proc : (num_dim_neigh_keys - bi);

			int row_cache_info[2];
			row_cache_info[0] = bi;
			row_cache_info[1] = bl_rows;

			LOG_DEBUG("Rank[%d]>> row-cache-info [bl_start=%d, bl_rows=%d] for %s",
					rank, row_cache_info[0], row_cache_info[1], str_dim_type.c_str());

			vector<hbin_t*> dtyps1_int_v(bl_rows);
			vector<vector<double>> dtyps1_extrm_v(bl_rows,
					vector<double>(4, 0.0));
			vector<CoordBAT> crd_rows(bl_rows,
					CoordBAT(0, n_schemes, nframes_tot_eff));
			std: map<u_int, u_int> row_id2index;
			vector<u_int> vec_row_id2index(2 * bl_rows);
			// Fill Cache with rows
			auto _startr = chrono::high_resolution_clock::now();
			for (auto rno = 0; rno < bl_rows; ++rno) {
				const u_int dtyp_id1 = dim_neigh_keys[bi + rno];
				row_id2index[dtyp_id1] = rno;
				vec_row_id2index[2 * rno] = dtyp_id1;
				vec_row_id2index[2 * rno + 1] = rno;
				crd_rows[rno].setId(dtyp_id1);
				if (intxofXY2DTraj->readCoords(dtyp_id1, 0, n_schemes,
						crd_rows[rno]) != 0) {
					LOG_ERROR("Rank[%d]>> reading row-ints of id(%d) for %s",
							rank, dtyp_id1, str_dim_type.c_str());
				}
				crd_rows[rno].getCoords(&dtyps1_extrm_v[rno][0],
						&dtyps1_extrm_v[rno][1], &dtyps1_extrm_v[rno][2],
						&dtyps1_extrm_v[rno][3], &n_schemes, &nfrm_eff_entropy,
						&dtyps1_int_v[rno]);
			}
			auto _endr = chrono::high_resolution_clock::now();
			auto _dur_rd = chrono::duration_cast<chrono::nanoseconds>(
					_endr - _startr).count();
			ptaskset->curr_read += chrono::nanoseconds(_dur_rd);
			progressstate.entc.curr_read += chrono::nanoseconds(_dur_rd);

			fflush(stdout);

			/*
			 * Send data for bin frequency calculation and entropy estimation to the slave process
			 */
			for (int idx = 0; idx < bl_rows; ++idx) {
				const u_int dim_id = dim_neigh_keys[bi + idx];
				MPI_Bcast(dtyps1_extrm_v[idx].data(), 4, MPI_DOUBLE,
				MASTER_PROC, MPI_COMM_WORLD);
				MPI_Bcast(dtyps1_int_v[idx], n_schemes * nframes_tot_eff,
				MPI_UNSIGNED_CHAR, MASTER_PROC, MPI_COMM_WORLD);

				LOG_DEBUG("Rank[%d]>> Bcast %d ints of row id(%d) for %s", rank, n_schemes * nframes_tot_eff, dim_id, str_dim_type.c_str());
			}
			auto _endcomm1 = chrono::high_resolution_clock::now();
			auto _durcomm1 = chrono::duration_cast<chrono::nanoseconds>(
					_endcomm1 - _endr).count();
			ptaskset->curr_comm += chrono::nanoseconds(_durcomm1);
			progressstate.entc.curr_comm += chrono::nanoseconds(_durcomm1);
			fflush(stdout);
			MPI_Barrier(MPI_COMM_WORLD);
			for (u_int bj = 0; bj < dim_lens[1]; bj += chached_dims) {
				const int bl_cols =
						(bj + chached_dims <= dim_lens[1]) ?
								chached_dims : (dim_lens[1] - bj);

				int curr_col_cache_per_proc = chache_dims_per_proc;
				int max_slave_rank = numprocs;
				if (bl_cols < chached_dims) {
					curr_col_cache_per_proc = ceil((double) bl_cols / numprocs);
					max_slave_rank = ceil(
							(double) bl_cols / (curr_col_cache_per_proc));
				}
				int col_cache_info[4];
				col_cache_info[0] = bj;
				col_cache_info[1] = bl_cols;
				col_cache_info[2] = curr_col_cache_per_proc;
				col_cache_info[3] = max_slave_rank;

				MPI_Barrier(MPI_COMM_WORLD);

				int n_cols_sent2slaves = 0;
				for (auto process_idx = 1; process_idx < max_slave_rank;
						process_idx++, n_cols_sent2slaves +=
								curr_col_cache_per_proc) {
					// Calc
					u_int block_tasks = 0, col_idx = 0;
					map<u_int, u_int> read_cols_dims;
					map<pair<u_int, u_int>, u_int> id2index;

					for (auto ri = bi; ri < bi + bl_rows; ++ri) {
						const u_int dtyp_id1 = dim_neigh_keys[ri];
						const vector<u_int> neigh =
								inputs.getNeighbors().getDimTypeNeighs(dim_type,
										dtyp_id1);
						auto min_cj = dimtypes_v[bj
								+ (process_idx - 1) * curr_col_cache_per_proc];
						auto max_cj = dimtypes_v[bj
								+ (process_idx * curr_col_cache_per_proc) - 1];
						for (auto const cj : neigh) {
							if (min_cj <= cj && cj <= max_cj) {
								if (!read_cols_dims.count(cj)) {
									read_cols_dims[cj] = col_idx;
									++col_idx;
								}
								id2index[make_pair(dtyp_id1, cj)] = block_tasks;
								++block_tasks;
							}
						}
					}
					u_int proc_cols_block_tasks[2];
					proc_cols_block_tasks[0] = block_tasks;
					proc_cols_block_tasks[1] = read_cols_dims.size();
					vector<u_int> vec_read_cols_dims(2 * read_cols_dims.size());
					vector<u_int> vec_id2index(3 * block_tasks);
					MPI_Send(&proc_cols_block_tasks, 2, MPI_UNSIGNED,
							process_idx, XnY_CBLK_TAG, MPI_COMM_WORLD);

					LOG_DEBUG("Rank[%d]>> => Rank[%d] sent for %s [block_tasks=%d, read_cols_dims=%d]",
							rank, process_idx, str_dim_type.c_str(), proc_cols_block_tasks[0],
							proc_cols_block_tasks[1]);

					if (block_tasks > 0) {
						auto itmp = 0;
						for (auto it = read_cols_dims.begin();
								it != read_cols_dims.end(); ++it) {
							vec_read_cols_dims[itmp] = it->first;
							vec_read_cols_dims[itmp + 1] = it->second;
							itmp += 2;
						}
						itmp = 0;
						for (map<pair<u_int, u_int>, u_int>::iterator it =
								id2index.begin(); it != id2index.end(); ++it) {
							const pair<u_int, u_int> &id_pair = it->first;
							vec_id2index[itmp] = id_pair.first;
							vec_id2index[itmp + 1] = id_pair.second;
							vec_id2index[itmp + 2] = it->second;
							itmp += 3;
						}
						MPI_Send(vec_read_cols_dims.data(),
								vec_read_cols_dims.size(), MPI_UNSIGNED,
								process_idx, XnY_CIDX_TAG, MPI_COMM_WORLD);
						MPI_Send(vec_id2index.data(), vec_id2index.size(),
						MPI_UNSIGNED, process_idx, XnY_CMAP_TAG,
						MPI_COMM_WORLD);

						u_int n_cache_cols = read_cols_dims.size();
						vector<hbin_t*> dtyps2_int_v(n_cache_cols);
						vector<vector<double>> bnds2_extrm_v(n_cache_cols,
								vector<double>(4, 0.0));
						vector<CoordBAT> crd_cols(n_cache_cols,
								CoordBAT(0, n_schemes, nframes_tot_eff));

						// read and cache cols
						auto _startr2 = chrono::high_resolution_clock::now();
						for (map<u_int, u_int>::iterator it =
								read_cols_dims.begin();
								it != read_cols_dims.end(); ++it) {
							//const u_int dtyp_id2 = it->first;
							crd_cols[it->second].setId(it->first);
							if (intyofXY2DTraj->readCoords(it->first, 0,
									n_schemes, crd_cols[it->second]) != 0) {
								LOG_ERROR(
										"Rank[%d]>> reading %d ints of column id(%d) for %s",
										rank, n_schemes * nframes_tot_eff,
										it->first, str_dim_type.c_str());
							}
							crd_cols[it->second].getCoords(
									&bnds2_extrm_v[it->second][0],
									&bnds2_extrm_v[it->second][1],
									&bnds2_extrm_v[it->second][2],
									&bnds2_extrm_v[it->second][3], &n_schemes,
									&nfrm_eff_entropy,
									&dtyps2_int_v[it->second]);
						}
						auto _endr2 = chrono::high_resolution_clock::now();
						auto _dur2_rd = chrono::duration_cast<
								chrono::nanoseconds>(_endr2 - _startr2).count();
						ptaskset->curr_read += chrono::nanoseconds(_dur2_rd);
						progressstate.entc.curr_read += chrono::nanoseconds(
								_dur2_rd);

						for (map<u_int, u_int>::iterator it =
								read_cols_dims.begin();
								it != read_cols_dims.end(); ++it) {

							MPI_Send(bnds2_extrm_v[it->second].data(), 4,
									MPI_DOUBLE, process_idx, XnY_EXTRM_TAG,
									MPI_COMM_WORLD);
							MPI_Send(dtyps2_int_v[it->second],
									n_schemes * nframes_tot_eff,
									MPI_UNSIGNED_CHAR, process_idx,
									XnY_INT_TAG, MPI_COMM_WORLD);
							LOG_DEBUG("Rank[%d]>> => Rank[%d] sent %d ints of column-id(%d) for for %s",
									rank, process_idx, n_schemes * nframes_tot_eff, it->first, str_dim_type.c_str());
						}

						auto _endcomm1 = chrono::high_resolution_clock::now();
						auto _durcomm1 = chrono::duration_cast<
								chrono::nanoseconds>(_endcomm1 - _endr).count();
						ptaskset->curr_comm += chrono::nanoseconds(_durcomm1);
						progressstate.entc.curr_comm += chrono::nanoseconds(
								_durcomm1);
					}
				}

				LOG_DEBUG("Rank[%d]>> scattering columns ints to slaves completed", rank);

				int n_cols4master = bl_cols - n_cols_sent2slaves;

				u_int master_tasks = 0, mastercol_idx = 0;
				map<u_int, u_int> masterread_cols_dims;
				map<pair<u_int, u_int>, u_int> masterid2index;

				for (auto ri = bi; ri < bi + bl_rows; ++ri) {
					const u_int dtyp_id1 = dim_neigh_keys[ri];
					const vector<u_int> neigh =
							inputs.getNeighbors().getDimTypeNeighs(dim_type,
									dtyp_id1);
					auto min_cj = dimtypes_v[bj + n_cols_sent2slaves];
					auto max_cj = dimtypes_v[bj + bl_cols - 1];
					for (auto const cj : neigh) {
						if (cj >= min_cj && cj <= max_cj) {
							if (!masterread_cols_dims.count(cj)) {
								masterread_cols_dims[cj] = mastercol_idx;
								++mastercol_idx;
							}
							masterid2index[make_pair(dtyp_id1, cj)] =
									master_tasks;
							++master_tasks;
						} else if (cj > max_cj) {
							break;
						}
					}
				}

				u_int n_mastercache_cols = masterread_cols_dims.size();
				LOG_DEBUG("Rank[%d]>> data for %s [block_tasks=%d, read_cols_dims=%d]",
						rank, str_dim_type.c_str(), master_tasks, n_mastercache_cols);

				vector<hbin_t*> dtyps2_int_v_master;
				vector<vector<double>> dtyps2_extrm_v_master;
				vector<CoordBAT> crd_cols_master;

				vector<vector<ull_int>> freq2_xy_obs_v_master;
				vector<vector<double>> bin2_xy_mids_v_master;

				vector<vector<double>> entropy_xy_v_master;

				vector<u_int> dummy_ids_master(2, -1);
				vector<BinGroup> bingrp_v_master;
				vector<DimEntropy> dimentropy_v_master;

				vector<u_int> vec_masterid2index;

				if (master_tasks > 0) {
					vec_masterid2index.resize(3 * master_tasks, 0);
					auto ii = 0;
					for (map<pair<u_int, u_int>, u_int>::iterator it =
							masterid2index.begin(); it != masterid2index.end();
							++it) {
						const pair<u_int, u_int> &id_pair = it->first;
						vec_masterid2index[ii] = id_pair.first;
						vec_masterid2index[ii + 1] = id_pair.second;
						vec_masterid2index[ii + 2] = it->second;
						ii += 3;
					}
					for (int it4master = 0; it4master < master_tasks;
							++it4master) {
						hbin_t *hbptr;
						dtyps2_int_v_master.push_back(hbptr);
						dtyps2_extrm_v_master.push_back(vector<double>(4, 0.0));
						crd_cols_master.push_back(
								CoordBAT(0, n_schemes, nframes_tot_eff));
						freq2_xy_obs_v_master.push_back(
								vector<ull_int>(nsteps * nbins_sum_sq));
						bin2_xy_mids_v_master.push_back(
								vector<double>(2 * bin_schemes_sum, 0.0));
						entropy_xy_v_master.push_back(
								vector<double>(nestimators * nsteps * n_schemes,
										0.0));

						bingrp_v_master.push_back(
								BinGroup(dummy_ids_master, (hbin_t) nsteps,
										2 * bin_schemes_sum, nbins_sum_sq));

						dimentropy_v_master.push_back(
								DimEntropy(dummy_ids_master, (hbin_t) nsteps,
										n_schemes, nestimators));
					}

					// read and cache cols
					auto _startr2 = chrono::high_resolution_clock::now();
					for (map<u_int, u_int>::iterator it =
							masterread_cols_dims.begin();
							it != masterread_cols_dims.end(); ++it) {
						crd_cols_master[it->second].setId(it->first);
						if (intyofXY2DTraj->readCoords(it->first, 0, n_schemes,
								crd_cols_master[it->second]) != 0) {
							LOG_ERROR(
									"Rank[%d]>> reading %d ints of column id(%d) for %s",
									rank, n_schemes * nframes_tot_eff,
									it->first, str_dim_type.c_str());

						}
						crd_cols_master[it->second].getCoords(
								&dtyps2_extrm_v_master[it->second][0],
								&dtyps2_extrm_v_master[it->second][1],
								&dtyps2_extrm_v_master[it->second][2],
								&dtyps2_extrm_v_master[it->second][3],
								&n_schemes, &nfrm_eff_entropy,
								&dtyps2_int_v_master[it->second]);
					} LOG_DEBUG("Rank[%d]>> Reading columns for master's tasks completed for %s", rank, str_dim_type.c_str());

					auto _endr2 = chrono::high_resolution_clock::now();
					auto _dur2_rd = chrono::duration_cast<chrono::nanoseconds>(
							_endr2 - _startr2).count();
					ptaskset->curr_read += chrono::nanoseconds(_dur2_rd);
					progressstate.entc.curr_read += chrono::nanoseconds(
							_dur2_rd);
					//

					auto _endcomm1 = chrono::high_resolution_clock::now();

#pragma omp parallel for
					for (int idx = 0; idx < master_tasks; ++idx) {
						auto rid = vec_masterid2index[3 * idx];
						auto cid = vec_masterid2index[3 * idx + 1];
						auto rno = row_id2index[rid];
						auto cno = masterread_cols_dims[cid];

						LOG_DEBUG("Rank[%d]>> computing entropy of %s [tasks=%d idx=%d rid=%d, rno=%d, cid=%d, cno=%d]",
								rank, str_dim_type.c_str(), master_tasks, idx, rid, rno, cid, cno);

						if (binData2D(nsteps, bin_schemes, nframes_tot_eff,
								dtyps1_extrm_v[rno][0], dtyps1_extrm_v[rno][1],
								dtyps2_extrm_v_master[cno][0],
								dtyps2_extrm_v_master[cno][1],
								dtyps1_int_v[rno], dtyps2_int_v_master[cno],
								freq2_xy_obs_v_master[idx],
								bin2_xy_mids_v_master[idx]) != 0) {
							LOG_ERROR("Rank[%d]>> binning data %s id(%d, %d)",
									str_dim_type.c_str(), rid, cid);
						}
						entropy2D(dtype1st, dtype2nd, rid, cid,
								inputs.getEstimators(), nsteps, step_size,
								dtyps1_extrm_v[rno][0], dtyps1_extrm_v[rno][1],
								dtyps2_extrm_v_master[cno][0],
								dtyps2_extrm_v_master[cno][1],
								inputs.getEntropy().isJacobian(), isKDE,
								bin_schemes, freq2_xy_obs_v_master[idx],
								entropy_xy_v_master[idx]);

						vector<u_int> dim_ids(2);
						dim_ids[0] = rid;
						dim_ids[1] = cid;
						vector<double> bb_extrm(8, 0.0);
						for (int xx = 0; xx < 4; ++xx) {
							bb_extrm[xx] = dtyps1_extrm_v[rno][xx];
							bb_extrm[xx + 4] = dtyps2_extrm_v_master[cno][xx];
						}
						if (isWritefreq && (frqWriteset & dim_type)) {
							bingrp_v_master[idx].setId(0, 2, dim_ids);
							bingrp_v_master[idx].setExtremes(bb_extrm);
							bingrp_v_master[idx].setBinMids(0,
									2 * bin_schemes_sum,
									bin2_xy_mids_v_master[idx]);
							bingrp_v_master[idx].setBinFreqs(
									xy2DHist->getFirstStep(),
									xy2DHist->getStepStride(),
									xy2DHist->getStepsEff(), 0, nbins_sum_sq,
									freq2_xy_obs_v_master[idx]);
						}
						dimentropy_v_master[idx].setId(0, 2, dim_ids);
						dimentropy_v_master[idx].setDimContri(0, nestimators, 0,
								nsteps, 0, n_schemes, entropy_xy_v_master[idx]);

					}
					auto _endcomp = chrono::high_resolution_clock::now();
					auto _durcomp = chrono::duration_cast<chrono::nanoseconds>(
							_endcomp - _endcomm1).count();
					ptaskset->curr_comp += chrono::nanoseconds(_durcomp);
					progressstate.entc.curr_comp += chrono::nanoseconds(
							_durcomp);

				} LOG_DEBUG("Rank[%d]>> Computing entropy for master completed..", rank);
				// start of receiving entropy and frequency data from slaves and writing to files
				int n_gather_curr_block = 0;
				int n_prepared2write_curr_block = 0;
				for (size_t process_idx = 1; process_idx < max_slave_rank;
						++process_idx) {
					u_int blcoks_tasks_slave = 0;
					u_int slave_cols_block_tasks[2];
					MPI_Recv(&slave_cols_block_tasks, 2, MPI_UNSIGNED,
							process_idx,
							XnY_CBLK_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					blcoks_tasks_slave = slave_cols_block_tasks[0];

					LOG_DEBUG("Rank[%d]>> <= Rank[%d] received for %s [block_tasks=%d, read_cols_dims=%d]",
							rank, process_idx, str_dim_type.c_str(), slave_cols_block_tasks[0],
							slave_cols_block_tasks[1]);
					if (blcoks_tasks_slave > 0) {

						vector<u_int> vec_read_cols_dims(
								2 * slave_cols_block_tasks[1]);
						vector<u_int> vec_id2index(3 * blcoks_tasks_slave);
						vector<double> dd_extrema_block(8 * blcoks_tasks_slave,
								0.0);

						MPI_Recv(vec_read_cols_dims.data(),
								2 * slave_cols_block_tasks[1],
								MPI_UNSIGNED, process_idx, XnY_CIDX_TAG,
								MPI_COMM_WORLD,
								MPI_STATUS_IGNORE);

						MPI_Recv(vec_id2index.data(), 3 * blcoks_tasks_slave,
						MPI_UNSIGNED, process_idx, XnY_CMAP_TAG,
						MPI_COMM_WORLD, MPI_STATUS_IGNORE);

						LOG_DEBUG("Rrank[%d]>> <= Rank[%d] received block_task=%d row/col index-id vectors for %s",
								rank, process_idx, blcoks_tasks_slave, str_dim_type.c_str());
						MPI_Recv(dd_extrema_block.data(),
								8 * blcoks_tasks_slave,
								MPI_DOUBLE, process_idx, XnY_EXTRM_TAG,
								MPI_COMM_WORLD,
								MPI_STATUS_IGNORE);
						LOG_DEBUG("Rank[%d]>> <= Rank[%d] received block_task=%d row/col index-id extremas=%d vectors for %s",
								rank, process_idx, blcoks_tasks_slave, 8*blcoks_tasks_slave, str_dim_type.c_str());

						vector<vector<ull_int>> freq2_bb_obs_v(
								blcoks_tasks_slave,
								vector<ull_int>(nsteps * nbins_sum_sq));
						vector<vector<double>> bin2_bb_mids_v(
								blcoks_tasks_slave,
								vector<double>(2 * bin_schemes_sum, 0.0));
						vector<vector<double>> entropy_bb_v(blcoks_tasks_slave,
								vector<double>(nestimators * nsteps * n_schemes,
										0.0));

						vector<double>(4, 0.0);
						vector<u_int> dim_ids(2);

						vector<BinGroup> bingrp_v(blcoks_tasks_slave,
								BinGroup(dummy_ids_master, (hbin_t) nsteps,
										2 * bin_schemes_sum, nbins_sum_sq));

						vector<DimEntropy> dimentropy_v(blcoks_tasks_slave,
								DimEntropy(dummy_ids_master, (hbin_t) nsteps,
										n_schemes, nestimators));

						auto _startcomm = chrono::high_resolution_clock::now();
						for (int idx = 0; idx < blcoks_tasks_slave; ++idx) {
							MPI_Recv(bin2_bb_mids_v[idx].data(),
									2 * bin_schemes_sum,
									MPI_DOUBLE, process_idx, XnY_EXTRM_TAG,
									MPI_COMM_WORLD,
									MPI_STATUS_IGNORE);

							if (isWritefreq && (frqWriteset & dim_type)) {
								MPI_Recv(freq2_bb_obs_v[idx].data(),
										nsteps * nbins_sum_sq,
										MPI_UNSIGNED_LONG_LONG, process_idx,
										XnY_FREQ_TAG, MPI_COMM_WORLD,
										MPI_STATUS_IGNORE);
							}
							MPI_Recv(entropy_bb_v[idx].data(),
									nestimators * nsteps * n_schemes,
									MPI_DOUBLE, process_idx, XnY_ENTR_TAG,
									MPI_COMM_WORLD, MPI_STATUS_IGNORE);

							n_gather_curr_block++;
							vector<u_int> dim_ids(2);
							dim_ids[0] = vec_id2index[3 * idx];
							dim_ids[1] = vec_id2index[3 * idx + 1];

							vector<double> bb_extrm(8, 0.0);
							for (int xx = 0; xx < 8; ++xx) {
								bb_extrm[xx] = dd_extrema_block[8 * idx + xx];
							}
							if (isWritefreq && (frqWriteset & dim_type)) {
								bingrp_v[idx].setId(0, 2, dim_ids);
								bingrp_v[idx].setExtremes(bb_extrm);
								bingrp_v[idx].setBinMids(0, 2 * bin_schemes_sum,
										bin2_bb_mids_v[idx]);
								bingrp_v[idx].setBinFreqs(
										xy2DHist->getFirstStep(),
										xy2DHist->getStepStride(),
										xy2DHist->getStepsEff(), 0,
										nbins_sum_sq, freq2_bb_obs_v[idx]);
							}
							dimentropy_v[idx].setId(0, 2, dim_ids);
							dimentropy_v[idx].setDimContri(0, nestimators, 0,
									nsteps, 0, n_schemes, entropy_bb_v[idx]);

							LOG_DEBUG("Rank[%d]>> <= Rank[%d] received entropy contribution of [row-id=%d, col-id=%d] for %s",
									rank, process_idx, dim_ids[0], dim_ids[1], str_dim_type.c_str());
						}
						if (isWritefreq && (frqWriteset & dim_type)) {
							if (xy2DHist->writeRecords(bingrp_v) != 0) {
								LOG_ERROR(
										"Rank[%d]>> writing histogram of block_tasks=%d for %s from Rank[%d]",
										rank, block_tasks, str_dim_type.c_str(),
										process_idx);
							}
						}
						ent_xy2d->writeRecords(dimentropy_v);
						auto _endcomm = chrono::high_resolution_clock::now();
						auto _durcomm =
								chrono::duration_cast<chrono::nanoseconds>(
										_startcomm - _endcomm).count();
						ptaskset->curr_comm += chrono::nanoseconds(_durcomm);
						progressstate.entc.curr_comm += chrono::nanoseconds(
								_durcomm);
					}

					ptaskset->done_tasks += blcoks_tasks_slave;
					progressstate.entc.done_tasks += blcoks_tasks_slave;
					entc_tasks_done += blcoks_tasks_slave;
				}
				// End of receiving entropy and frequency data from slaves and writing to files

				// Start writing entropy and frequency data for master_tasks
				if (master_tasks > 0) {
					auto _strtwrt = chrono::high_resolution_clock::now();
					if (isWritefreq && (frqWriteset & BATSet::BB2D)) {
						if (xy2DHist->writeRecords(bingrp_v_master) != 0) {
							LOG_ERROR(
									"Rank[%d]>> writing histogram of block_tasks=%d for %s",
									rank, block_tasks, str_dim_type.c_str());
						}
					}
					ent_xy2d->writeRecords(dimentropy_v_master);

					auto _endwrt = chrono::high_resolution_clock::now();
					auto _durwrt = chrono::duration_cast<chrono::nanoseconds>(
							_endwrt - _strtwrt).count();
					ptaskset->curr_write += chrono::nanoseconds(_durwrt);
					progressstate.entc.curr_write += chrono::nanoseconds(
							_durwrt);
					ptaskset->current = chrono::high_resolution_clock::now();
					progressstate.entc.current = ptaskset->current;
					ptaskset->done_tasks += master_tasks;
					progressstate.entc.done_tasks += master_tasks;

					entc_tasks_done += master_tasks;
				}

				if (entc_tasks_done >= entc_tasks_freq) {
					ofstream info_strm(
							inputs.getControl().getOutfilepath()
									+ inputs.getControl().getInfofile());
					info_strm << progressstate.toString();
					info_strm.close();
					entc_tasks_done %= entc_tasks_freq;
				}
			}
		} // dih-dih block calc for d1 ends

		ptaskset->current = chrono::high_resolution_clock::now();
		ptaskset->cstt = ExecState::COMPLETED;
		progressstate.entc.current = ptaskset->current;
		int entc_workset_intval =
				static_cast<int>(inputs.getEntropy().getWorkSet());
		int entc_bd2d_intval = static_cast<int>(BATSet::BD2D);
		int entc_ad2d_intval = static_cast<int>(BATSet::AD2D);
		int entc_xy2d_intval = static_cast<int>(BATSet::XY2D);
		if ((entc_workset_intval < entc_bd2d_intval)
				&& (dim_type == BATSet::BA2D)) {
			progressstate.entc.cstt = ExecState::COMPLETED;
		} else if ((entc_workset_intval < entc_ad2d_intval)
				&& (dim_type == BATSet::BD2D)) {
			progressstate.entc.cstt = ExecState::COMPLETED;

		} else if ((entc_workset_intval <= entc_xy2d_intval)
				&& (dim_type == BATSet::AD2D)) {
			progressstate.entc.cstt = ExecState::COMPLETED;

		}
		ofstream info_strm(
				inputs.getControl().getOutfilepath()
						+ inputs.getControl().getInfofile());
		info_strm << progressstate.toString();
		info_strm.close();
		time_xy.stop();
		mprintf(
				"TIME: Histogram bin frequency calculation for %s-2D: %.4f seconds.\n",
				str_dim_type.c_str(), time_xy.total());

	} else if (rank != MASTER_PROC) {
		ProgTaskSet *ptaskset;
		string str_dim_type;
		u_int n_dim_x_eff = 0;
		u_int n_dim_y_eff = 0;
		u_int num_xy2d;
		Netcdf_TrjInt *intxofXY2DTraj;
		Netcdf_TrjInt *intyofXY2DTraj;
		string intxofXY2DTrajFileName;
		string intyofXY2DTrajFileName;
		Netcdf_HistUtil *xy2DHist;
		string xy2DHistFileName;
		Netcdf_EntContri *ent_xy2d;
		string ent_xy2dFileName;

		vector<u_int> dim_neigh_keys;
		vector<u_int> dimtypes_v;
		BAT_t dtype1st = BAT_t::NONE;
		BAT_t dtype2nd = BAT_t::NONE;
		switch (dim_type) {
		case BATSet::BA2D:
			dtype1st = BAT_t::BOND;
			dtype2nd = BAT_t::ANGLE;
			ptaskset = &(progressstate.entc_ba2d);
			str_dim_type.assign("B/A");
			n_dim_x_eff = n_bnd_eff;
			n_dim_y_eff = n_ang_eff;

			intxofXY2DTrajFileName.assign("bin_bonds.nc");
			intyofXY2DTrajFileName.assign("bin_angles.nc");
			xy2DHistFileName.assign("hist_ba-2d.nc");
			ent_xy2dFileName.assign("entcontri_ba-2d.nc");

			inputs.getNeighbors().bacrossKeys(dim_neigh_keys);
			for (size_t d1 = 0; d1 < dim_neigh_keys.size(); ++d1) {
				num_xy2d += inputs.getNeighbors().baNeighSize(
						dim_neigh_keys[d1]);
			}

			copy(inputs.getSubset().getAngles().begin(),
					inputs.getSubset().getAngles().end(),
					back_inserter(dimtypes_v));
			break;
		case BATSet::BD2D:
			dtype1st = BAT_t::BOND;
			dtype2nd = BAT_t::DIHEDRAL;
			ptaskset = &(progressstate.entc_bd2d);
			str_dim_type.assign("B/D");
			n_dim_x_eff = n_bnd_eff;
			n_dim_y_eff = n_dih_eff;

			intxofXY2DTrajFileName.assign("bin_bonds.nc");
			intyofXY2DTrajFileName.assign("bin_torsions.nc");
			xy2DHistFileName.assign("hist_bd-2d.nc");
			ent_xy2dFileName.assign("entcontri_bd-2d.nc");

			inputs.getNeighbors().bdcrossKeys(dim_neigh_keys);
			for (size_t d1 = 0; d1 < dim_neigh_keys.size(); ++d1) {
				num_xy2d += inputs.getNeighbors().bdNeighSize(
						dim_neigh_keys[d1]);
			}

			copy(inputs.getSubset().getTorsions().begin(),
					inputs.getSubset().getTorsions().end(),
					back_inserter(dimtypes_v));
			break;
		case BATSet::AD2D:
			dtype1st = BAT_t::ANGLE;
			dtype2nd = BAT_t::DIHEDRAL;
			ptaskset = &(progressstate.entc_ad2d);
			str_dim_type.assign("A/D");
			n_dim_x_eff = n_ang_eff;
			n_dim_y_eff = n_dih_eff;

			intxofXY2DTrajFileName.assign("bin_angles.nc");
			intyofXY2DTrajFileName.assign("bin_torsions.nc");
			xy2DHistFileName.assign("hist_ad-2d.nc");
			ent_xy2dFileName.assign("entcontri_ad-2d.nc");

			inputs.getNeighbors().adcrossKeys(dim_neigh_keys);
			for (size_t d1 = 0; d1 < dim_neigh_keys.size(); ++d1) {
				num_xy2d += inputs.getNeighbors().adNeighSize(
						dim_neigh_keys[d1]);
			}
			copy(inputs.getSubset().getTorsions().begin(),
					inputs.getSubset().getTorsions().end(),
					back_inserter(dimtypes_v));
			break;
		default:
			LOG_ERROR("Rank[%d]>> Invalid xy2D type found", rank);
			exit(0);
		}

		vector<ull_int> dim_lens(2, 0);
		dim_lens[0] = n_dim_x_eff;
		dim_lens[1] = n_dim_y_eff;
		ull_int nfrm_eff_entropy = 0;
		/*********************** MASTER PROCESS *******************
		 * 2-D ENTROPY estimation for: BOND/ANGLE
		 *********************************************************/
		int chache_dims_per_proc = cache_entc_dims_per_cpu * n_thread_perproc
				/ 2;
		int n_cpus = numprocs * n_thread_perproc;
		int chached_dims = chache_dims_per_proc * numprocs;

		int block_count = 0, block_tasks = 0;
		bool block_ready = false;

		const u_int n_dim_neigh_keys = dim_neigh_keys.size();

		//u_int d1 = 0;
		// Process local variables for data receiving and processing
		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);
		//int chached_dims = cache_entc_dims_per_cpu * numprocs / 2;
		for (u_int bi = 0; bi < n_dim_neigh_keys; bi += chache_dims_per_proc) {
			const int bl_rows =
					(bi + chache_dims_per_proc <= n_dim_neigh_keys) ?
							chache_dims_per_proc : (n_dim_neigh_keys - bi);

			int row_cache_info[2];
			row_cache_info[0] = bi;
			row_cache_info[1] = bl_rows;

			LOG_DEBUG("Rank[%d]>> row-cache-info [bl_start=%d, bl_rows=%d] for %s",
					rank, row_cache_info[0], row_cache_info[1], str_dim_type.c_str());

			vector<vector<hbin_t>> dtyps1_int_v(bl_rows,
					vector<hbin_t>(n_schemes * nframes_tot_eff));
			vector<vector<double>> dtyps1_extrm_v(bl_rows,
					vector<double>(4, 0.0));
			vector<CoordBAT> crd_rows(bl_rows,
					CoordBAT(0, n_schemes, nframes_tot_eff));

			map<u_int, u_int> row_id2index;
			vector<u_int> vec_row_id2index(2 * bl_rows);

			// Fill Cache with rows
			fflush(stdout);
			/*
			 * Send data for bin frequency calculation and entropy estimation to the slave process
			 */
			for (int idx = 0; idx < bl_rows; ++idx) {
				const u_int dim_id = dim_neigh_keys[bi + idx];

				MPI_Bcast(dtyps1_extrm_v[idx].data(), 4, MPI_DOUBLE,
				MASTER_PROC, MPI_COMM_WORLD);
				MPI_Bcast(dtyps1_int_v[idx].data(), n_schemes * nframes_tot_eff,
				MPI_UNSIGNED_CHAR, MASTER_PROC, MPI_COMM_WORLD);

				row_id2index[dim_id] = idx;
				vec_row_id2index[2 * idx] = dim_id;
				vec_row_id2index[2 * idx + 1] = idx;

				LOG_DEBUG("Rank[%d]>> Bcast %d ints of row id(%d) for %s",
						rank, n_schemes * nframes_tot_eff, dim_id, str_dim_type.c_str());
			}

			MPI_Barrier(MPI_COMM_WORLD);
			for (u_int bj = 0; bj < dim_lens[1]; bj += chached_dims) {
				const int bl_cols =
						(bj + chached_dims <= dim_lens[1]) ?
								chached_dims : (dim_lens[1] - bj);

				int curr_col_cache_per_proc = chache_dims_per_proc;
				int max_slave_rank = numprocs;
				if (bl_cols < chached_dims) {
					curr_col_cache_per_proc = ceil((double) bl_cols / numprocs);
					max_slave_rank = ceil(
							(double) bl_cols / (curr_col_cache_per_proc));
				}
				int col_cache_info[4];
				col_cache_info[0] = bj;
				col_cache_info[1] = bl_cols;
				col_cache_info[2] = curr_col_cache_per_proc;
				col_cache_info[3] = max_slave_rank;

				MPI_Barrier(MPI_COMM_WORLD);
				// Calc
				u_int block_tasks = 0, col_idx = 0;
				u_int slave_cols_block_tasks[2];
				map<u_int, u_int> read_cols_dims;
				map<pair<u_int, u_int>, u_int> id2index;

				if (rank < max_slave_rank) {
					MPI_Recv(&slave_cols_block_tasks, 2, MPI_UNSIGNED,
					MASTER_PROC, XnY_CBLK_TAG, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);
					block_tasks = slave_cols_block_tasks[0];
					u_int n_cache_cols = slave_cols_block_tasks[1];

					LOG_DEBUG("Rank[%d]>> <= Rank[%d] received for %s [block_tasks=%d, read_cols_dims=%d]",
							rank, MASTER_PROC, str_dim_type.c_str(), slave_cols_block_tasks[0],
							slave_cols_block_tasks[1]);

					if (block_tasks > 0) {
						vector<u_int> vec_read_cols_dims(2 * n_cache_cols);
						vector<u_int> vec_id2index(3 * block_tasks);
						MPI_Recv(vec_read_cols_dims.data(), 2 * n_cache_cols,
						MPI_UNSIGNED, MASTER_PROC,
						XnY_CIDX_TAG, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);
						MPI_Recv(vec_id2index.data(), 3 * block_tasks,
						MPI_UNSIGNED, MASTER_PROC,
						XnY_CMAP_TAG, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE);

						for (int ia = 0; ia < 2 * n_cache_cols; ia += 2) {
							read_cols_dims[vec_read_cols_dims[ia]] =
									vec_read_cols_dims[ia + 1];
						}
						for (int ia = 0; ia < 3 * block_tasks; ia += 3) {
							id2index[make_pair(vec_id2index[ia],
									vec_id2index[ia + 1])] =
									vec_id2index[ia + 2];
						}

						vector<vector<hbin_t>> dtyps2_int_v(n_cache_cols,
								vector<hbin_t>(n_schemes * nframes_tot_eff));
						vector<vector<double>> bnds2_extrm_v(n_cache_cols,
								vector<double>(4, 0.0));
						vector<CoordBAT> crd_cols(n_cache_cols,
								CoordBAT(0, n_schemes, nframes_tot_eff));

						// read and cache cols
						auto _startr2 = chrono::high_resolution_clock::now();
						for (map<u_int, u_int>::iterator it =
								read_cols_dims.begin();
								it != read_cols_dims.end(); ++it) {
							MPI_Recv(bnds2_extrm_v[it->second].data(), 4,
							MPI_DOUBLE, MASTER_PROC,
							XnY_EXTRM_TAG, MPI_COMM_WORLD,
							MPI_STATUS_IGNORE);
							MPI_Recv(dtyps2_int_v[it->second].data(),
									n_schemes * nframes_tot_eff,
									MPI_UNSIGNED_CHAR, MASTER_PROC,
									XnY_INT_TAG, MPI_COMM_WORLD,
									MPI_STATUS_IGNORE);
							LOG_DEBUG("Rank[%d]>> <= Rank[%d] received %d ints of column-id(%d) for for %s",
									rank, MASTER_PROC, n_schemes * nframes_tot_eff, it->first, str_dim_type.c_str());
						}

						vector<vector<ull_int>> freq2_bb_obs_v(block_tasks,
								vector<ull_int>(nsteps * nbins_sum_sq));
						vector<vector<double>> bin2_bb_mids_v(block_tasks,
								vector<double>(2 * bin_schemes_sum, 0.0));

						vector<vector<double>> entropy_bb_v(block_tasks,
								vector<double>(nestimators * nsteps * n_schemes,
										0.0));

						vector<u_int> dummy_ids(2, -1);
						vector<BinGroup> bingrp_v(block_tasks,
								BinGroup(dummy_ids, (hbin_t) nsteps,
										2 * bin_schemes_sum, nbins_sum_sq));
						vector<DimEntropy> dimentropy_v(block_tasks,
								DimEntropy(dummy_ids, (hbin_t) nsteps,
										n_schemes, nestimators));

						auto _endcomm1 = chrono::high_resolution_clock::now();
						vector<double> bb_extrm(8 * block_tasks, 0.0);

#pragma omp parallel for
						for (int idx = 0; idx < block_tasks; ++idx) {
							auto rid = vec_id2index[3 * idx];
							auto cid = vec_id2index[3 * idx + 1];
							auto rno = row_id2index[rid];
							auto cno = read_cols_dims[cid];

							LOG_DEBUG("Rank[%d]>> computing entropy of %s [tasks=%d idx=%d rid=%d, rno=%d, cid=%d, cno=%d]",
									rank, str_dim_type.c_str(), block_tasks, idx, rid, rno, cid, cno);

							if (binData2D(nsteps, bin_schemes, nframes_tot_eff,
									dtyps1_extrm_v[rno][0],
									dtyps1_extrm_v[rno][1],
									bnds2_extrm_v[cno][0],
									bnds2_extrm_v[cno][1],
									dtyps1_int_v[rno].data(),
									dtyps2_int_v[cno].data(),
									freq2_bb_obs_v[idx], bin2_bb_mids_v[idx])
									!= 0) {
								LOG_ERROR(
										"Rank[%d]>> binning data %s id(%d, %d)",
										str_dim_type.c_str(), rid, cid);
							}
							entropy2D(dtype1st, dtype2nd, rid, cid,
									inputs.getEstimators(), nsteps, step_size,
									dtyps1_extrm_v[rno][0],
									dtyps1_extrm_v[rno][1],
									bnds2_extrm_v[cno][0],
									bnds2_extrm_v[cno][1],
									inputs.getEntropy().isJacobian(), isKDE,
									bin_schemes, freq2_bb_obs_v[idx],
									entropy_bb_v[idx]);

							vector<u_int> dim_ids(2);
							dim_ids[0] = rid;
							dim_ids[1] = cid;

							for (int xx = 0; xx < 4; ++xx) {
								bb_extrm[8 * idx + xx] =
										dtyps1_extrm_v[rno][xx];
								bb_extrm[(8 * idx) + xx + 4] =
										bnds2_extrm_v[cno][xx];
							}
						} // #end-omp parallel for

						LOG_DEBUG("Rank[%d]>> Computing entropy for master completed..", rank);

						MPI_Send(&slave_cols_block_tasks, 2, MPI_UNSIGNED,
						MASTER_PROC,
						XnY_CBLK_TAG, MPI_COMM_WORLD);
						MPI_Send(vec_read_cols_dims.data(),
								vec_read_cols_dims.size(), MPI_UNSIGNED,
								MASTER_PROC, XnY_CIDX_TAG, MPI_COMM_WORLD);
						MPI_Send(vec_id2index.data(), vec_id2index.size(),
						MPI_UNSIGNED, MASTER_PROC, XnY_CMAP_TAG,
						MPI_COMM_WORLD);

						LOG_DEBUG("Rrank[%d]>> => Rank[%d] sent block_task=%d row/col index-id vectors for %s",
								rank, MASTER_PROC, block_tasks, str_dim_type.c_str());

						MPI_Send(bb_extrm.data(), bb_extrm.size(),
						MPI_DOUBLE, MASTER_PROC, XnY_EXTRM_TAG,
						MPI_COMM_WORLD);

						LOG_DEBUG("Rank[%d]>> => Rank[%d] sent block_task=%d row/col index-id extremas=%d vectors for %s",
								rank, MASTER_PROC, block_tasks, 8*block_tasks, str_dim_type.c_str());

						for (int idx = 0; idx < block_tasks; ++idx) {
							MPI_Send(bin2_bb_mids_v[idx].data(),
									2 * bin_schemes_sum,
									MPI_DOUBLE, MASTER_PROC,
									XnY_EXTRM_TAG, MPI_COMM_WORLD);

							if (isWritefreq && (frqWriteset & dim_type)) {
								MPI_Send(freq2_bb_obs_v[idx].data(),
										nsteps * nbins_sum_sq,
										MPI_UNSIGNED_LONG_LONG, MASTER_PROC,
										XnY_FREQ_TAG, MPI_COMM_WORLD);
							}
							MPI_Send(entropy_bb_v[idx].data(),
									nestimators * nsteps * n_schemes,
									MPI_DOUBLE, MASTER_PROC, XnY_ENTR_TAG,
									MPI_COMM_WORLD);

							LOG_DEBUG("Rank[%d]>> => Rank[%d] sending entropy contribution of [row-id=%d, col-id=%d] for %s",
									rank, MASTER_PROC, vec_id2index[2*idx], vec_id2index[2*idx+1], str_dim_type.c_str());

						}
					} else if (block_tasks == 0) {
						MPI_Send(&slave_cols_block_tasks, 2, MPI_UNSIGNED,
						MASTER_PROC,
						XnY_CBLK_TAG, MPI_COMM_WORLD);

					}
				}
				//MPI_Barrier(MPI_COMM_WORLD);
			}
		} // dih-dih block calc for d1 ends
	}
}

void EntropyCalculator::run2d_xy(BATSet dim_type, int num_threads) {

	ProgTaskSet *ptaskset;
	string str_dim_type;
	u_int n_dim_x_eff = 0;
	u_int n_dim_y_eff = 0;
	u_int num_xy2d;
	Netcdf_TrjInt *intxofXY2DTraj;
	Netcdf_TrjInt *intyofXY2DTraj;
	string intxofXY2DTrajFileName;
	string intyofXY2DTrajFileName;
	Netcdf_HistUtil *xy2DHist;
	string xy2DHistFileName;
	Netcdf_EntContri *ent_xy2d;
	string ent_xy2dFileName;

	vector<u_int> dimTypeKeys;
	vector<u_int> dimtypes_v;
	BAT_t dtype1st = BAT_t::NONE;
	BAT_t dtype2nd = BAT_t::NONE;
	switch (dim_type) {
	case BATSet::BA2D:
		dtype1st = BAT_t::BOND;
		dtype2nd = BAT_t::ANGLE;
		ptaskset = &(progressstate.entc_ba2d);
		str_dim_type.assign("B/A");
		n_dim_x_eff = n_bnd_eff;
		n_dim_y_eff = n_ang_eff;

		intxofXY2DTrajFileName.assign("bin_bonds.nc");
		intyofXY2DTrajFileName.assign("bin_angles.nc");
		xy2DHistFileName.assign("hist_ba-2d.nc");
		ent_xy2dFileName.assign("entcontri_ba-2d.nc");

		inputs.getNeighbors().bacrossKeys(dimTypeKeys);
		for (size_t d1 = 0; d1 < dimTypeKeys.size(); ++d1) {
			num_xy2d += inputs.getNeighbors().baNeighSize(dimTypeKeys[d1]);
		}

		copy(inputs.getSubset().getAngles().begin(),
				inputs.getSubset().getAngles().end(),
				back_inserter(dimtypes_v));
		break;
	case BATSet::BD2D:
		dtype1st = BAT_t::BOND;
		dtype2nd = BAT_t::DIHEDRAL;
		ptaskset = &(progressstate.entc_bd2d);
		str_dim_type.assign("B/D");
		n_dim_x_eff = n_bnd_eff;
		n_dim_y_eff = n_dih_eff;

		intxofXY2DTrajFileName.assign("bin_bonds.nc");
		intyofXY2DTrajFileName.assign("bin_torsions.nc");
		xy2DHistFileName.assign("hist_bd-2d.nc");
		ent_xy2dFileName.assign("entcontri_bd-2d.nc");

		inputs.getNeighbors().bdcrossKeys(dimTypeKeys);
		for (size_t d1 = 0; d1 < dimTypeKeys.size(); ++d1) {
			num_xy2d += inputs.getNeighbors().bdNeighSize(dimTypeKeys[d1]);
		}

		copy(inputs.getSubset().getTorsions().begin(),
				inputs.getSubset().getTorsions().end(),
				back_inserter(dimtypes_v));
		break;
	case BATSet::AD2D:
		dtype1st = BAT_t::ANGLE;
		dtype2nd = BAT_t::DIHEDRAL;
		ptaskset = &(progressstate.entc_ad2d);
		str_dim_type.assign("A/D");
		n_dim_x_eff = n_ang_eff;
		n_dim_y_eff = n_dih_eff;

		intxofXY2DTrajFileName.assign("bin_angles.nc");
		intyofXY2DTrajFileName.assign("bin_torsions.nc");
		xy2DHistFileName.assign("hist_ad-2d.nc");
		ent_xy2dFileName.assign("entcontri_ad-2d.nc");

		inputs.getNeighbors().adcrossKeys(dimTypeKeys);
		for (size_t d1 = 0; d1 < dimTypeKeys.size(); ++d1) {
			num_xy2d += inputs.getNeighbors().adNeighSize(dimTypeKeys[d1]);
		}
		copy(inputs.getSubset().getTorsions().begin(),
				inputs.getSubset().getTorsions().end(),
				back_inserter(dimtypes_v));
		break;
	default:
		cout << "Error:: Invalid xy2D type found" << endl;
		exit(0);
	}

	vector<ull_int> dim_lens(2, 0);
	dim_lens[0] = n_dim_x_eff;
	dim_lens[1] = n_dim_y_eff;

	intxofXY2DTraj = new Netcdf_TrjInt(
			inputs.getControl().getOutfilepath() + intxofXY2DTrajFileName,
			str_dim_type, n_dim_x_eff, startframe, strideframe,
			(ull_int) nframes_tot, n_schemes, bin_schemes);

	intyofXY2DTraj = new Netcdf_TrjInt(
			inputs.getControl().getOutfilepath() + intyofXY2DTrajFileName,
			str_dim_type, n_dim_y_eff, startframe, strideframe,
			(ull_int) nframes_tot, n_schemes, bin_schemes);

	xy2DHist = new Netcdf_HistUtil(
			inputs.getControl().getOutfilepath() + xy2DHistFileName,
			str_dim_type, 2, 2, start_hist_step, nsteps, stride_hist_step,
			dim_lens, TensorType::UPPER, n_schemes, bin_schemes);

	if (isWritefreq && (frqWriteset & dim_type)) {
		xy2DHist->NC_create(str_dim_type + " 2-D histograms");
	}

	ent_xy2d = new Netcdf_EntContri(
			inputs.getControl().getOutfilepath() + ent_xy2dFileName,
			str_dim_type, 2, 2, nsteps, dim_lens, TensorType::UPPER,
			bin_schemes, nestimators);

	ent_xy2d->NC_create(str_dim_type + "-2D Entropy contributions");

	ptaskset->start = chrono::high_resolution_clock::now();
	ptaskset->current = ptaskset->start;
	ptaskset->strt_read = ptaskset->curr_read = ptaskset->start;
	ptaskset->strt_comp = ptaskset->curr_comp = ptaskset->start;
	ptaskset->strt_comm = ptaskset->curr_comm = ptaskset->start;
	ptaskset->strt_write = ptaskset->curr_write = ptaskset->start;
	ptaskset->cstt = ExecState::RUNNING;
	if (progressstate.entc.cstt != ExecState::RUNNING) {
		progressstate.entc.cstt = ExecState::RUNNING;
		progressstate.entc.start = ptaskset->start;
		progressstate.entc.current = ptaskset->start;
	}
	Timer time_xy;
	time_xy.start();
	ull_int nfrm_eff_entropy = 0;
	/*********************** MASTER PROCESS *******************
	 * 2-D ENTROPY estimation for: BOND/ANGLE
	 *********************************************************/

//vector<u_int> proc_rids(numprocs,
//		numeric_limits<u_int>::max());
// 5 elements: (start-d1, end-d1, start-d2, end-d2, block_tasks) for every block
	vector<u_int> block_boundry(5);
	int block_count = 0, block_tasks = 0;
	bool block_ready = false;

	const u_int num_xy2DKeys = dimTypeKeys.size();
//u_int d1 = 0;

	int chached_dims = cache_entc_dims_per_cpu * num_threads / 2;
	for (u_int bi = 0; bi < num_xy2DKeys; bi += chached_dims) {
		const int bl_rows =
				(bi + chached_dims <= num_xy2DKeys) ?
						chached_dims : (num_xy2DKeys - bi);

		vector<hbin_t*> dtyps1_int_v(bl_rows);
		vector<vector<double>> dtyps1_extrm_v(bl_rows, vector<double>(4, 0.0));
		vector<CoordBAT> crd_rows(bl_rows,
				CoordBAT(0, n_schemes, nframes_tot_eff));

		// Fill Cache with rows
		auto _startr = chrono::high_resolution_clock::now();
		for (auto rno = 0; rno < bl_rows; ++rno) {
			const u_int dtyp_id1 = dimTypeKeys[bi + rno];
			crd_rows[rno].setId(dtyp_id1);
			if (intxofXY2DTraj->readCoords(dtyp_id1, 0, n_schemes,
					crd_rows[rno]) != 0) {
				mprinterr("Error: Reading int ang id(%d)", dtyp_id1);
			}
			crd_rows[rno].getCoords(&dtyps1_extrm_v[rno][0],
					&dtyps1_extrm_v[rno][1], &dtyps1_extrm_v[rno][2],
					&dtyps1_extrm_v[rno][3], &n_schemes, &nfrm_eff_entropy,
					&dtyps1_int_v[rno]);
		}
		auto _endr = chrono::high_resolution_clock::now();
		auto _dur_rd = chrono::duration_cast<chrono::nanoseconds>(
				_endr - _startr).count();
		ptaskset->curr_read += chrono::nanoseconds(_dur_rd);
		progressstate.entc.curr_read += chrono::nanoseconds(_dur_rd);

		for (u_int bj = 0; bj < n_dim_y_eff; bj += chached_dims) {
			const int bl_cols =
					(bj + chached_dims <= n_dim_y_eff) ?
							chached_dims : (n_dim_y_eff - bj);

			// Calc
			u_int block_tasks = 0, col_idx = 0;
			map<u_int, u_int> read_cols_dims;
			map<pair<u_int, u_int>, u_int> id2index;

			for (auto ri = bi; ri < bi + bl_rows; ++ri) {
				const u_int dtyp_id1 = dimTypeKeys[ri];
				const vector<u_int> neigh =
						inputs.getNeighbors().getDimTypeNeighs(dim_type,
								dtyp_id1);
				for (auto const cj : neigh) {
					if (cj >= dimtypes_v[bj]
							&& cj <= dimtypes_v[bj + bl_cols - 1]) {
						if (!read_cols_dims.count(cj)) {
							read_cols_dims[cj] = col_idx;
							++col_idx;
						}
						id2index[make_pair(dtyp_id1, cj)] = block_tasks;
						++block_tasks;
					}
				}
			}
			if (block_tasks) {
				u_int n_cache_cols = read_cols_dims.size();
				vector<hbin_t*> dtyps2_int_v(n_cache_cols);
				vector<vector<double>> dtyps2_extrm_v(n_cache_cols,
						vector<double>(4, 0.0));
				vector<CoordBAT> crd_cols(n_cache_cols,
						CoordBAT(0, n_schemes, nframes_tot_eff));

				// read and cache cols
				auto _startr2 = chrono::high_resolution_clock::now();
				for (map<u_int, u_int>::iterator it = read_cols_dims.begin();
						it != read_cols_dims.end(); ++it) {
					//const u_int ang_id2 = it->first;
					crd_cols[it->second].setId(it->first);
					if (intyofXY2DTraj->readCoords(it->first, 0, n_schemes,
							crd_cols[it->second]) != 0) {
						mprinterr("Error: Reading int ang id(%d)", it->first);
					}
					crd_cols[it->second].getCoords(
							&dtyps2_extrm_v[it->second][0],
							&dtyps2_extrm_v[it->second][1],
							&dtyps2_extrm_v[it->second][2],
							&dtyps2_extrm_v[it->second][3], &n_schemes,
							&nfrm_eff_entropy, &dtyps2_int_v[it->second]);
				}
				auto _endr2 = chrono::high_resolution_clock::now();
				auto _dur2_rd = chrono::duration_cast<chrono::nanoseconds>(
						_endr2 - _startr2).count();
				ptaskset->curr_read += chrono::nanoseconds(_dur2_rd);
				progressstate.entc.curr_read += chrono::nanoseconds(_dur2_rd);
				//

				vector<vector<ull_int>> freq2_ad_obs_v(block_tasks,
						vector<ull_int>(nsteps * nbins_sum_sq));
				vector<vector<double>> bin2_ad_mids_v(block_tasks,
						vector<double>(2 * bin_schemes_sum, 0.0));

				vector<vector<double>> entropy_xy2d_v(block_tasks,
						vector<double>(nestimators * nsteps * n_schemes, 0.0));

				vector<u_int> dummy_ids(2, -1);
				vector<BinGroup> bingrp_v(block_tasks,
						BinGroup(dummy_ids, (hbin_t) nsteps,
								2 * bin_schemes_sum, nbins_sum_sq));
				vector<DimEntropy> dimentropy_v(block_tasks,
						DimEntropy(dummy_ids, (hbin_t) nsteps, n_schemes,
								nestimators));

				auto _endcomm1 = chrono::high_resolution_clock::now();

#pragma omp parallel for
				for (auto ri = bi; ri < bi + bl_rows; ++ri) {
					auto rno = ri - bi;
					const u_int dtyp_id1 = dimTypeKeys[ri];
					const vector<u_int> neigh =
							inputs.getNeighbors().getDimTypeNeighs(dim_type,
									dtyp_id1);
					for (auto const cj : neigh) {
						if (cj >= dimtypes_v[bj]
								&& cj <= dimtypes_v[bj + bl_cols - 1]) {
							u_int idx = id2index[make_pair(dtyp_id1, cj)];
							auto cno = read_cols_dims[cj];
							if (binData2D(nsteps, bin_schemes, nframes_tot_eff,
									dtyps1_extrm_v[rno][0],
									dtyps1_extrm_v[rno][1],
									dtyps2_extrm_v[cno][0],
									dtyps2_extrm_v[cno][1], dtyps1_int_v[rno],
									dtyps2_int_v[cno], freq2_ad_obs_v[idx],
									bin2_ad_mids_v[idx]) != 0) {
								mprinterr("Error: Binning data %s id(%d, %d)",
										str_dim_type.c_str(), dtyp_id1, cj);
							}
							entropy2D(dtype1st, dtype2nd, dtyp_id1, cj,
									inputs.getEstimators(), nsteps, step_size,
									dtyps1_extrm_v[rno][0],
									dtyps1_extrm_v[rno][1],
									dtyps2_extrm_v[cno][0],
									dtyps2_extrm_v[cno][1],
									inputs.getEntropy().isJacobian(), isKDE,
									bin_schemes, freq2_ad_obs_v[idx],
									entropy_xy2d_v[idx]);

							vector<u_int> dim_ids(2);
							dim_ids[0] = dtyp_id1;
							dim_ids[1] = cj;
							vector<double> ad_extrm(8, 0.0);
							for (int xx = 0; xx < 4; ++xx) {
								ad_extrm[xx] = dtyps1_extrm_v[rno][xx];
								ad_extrm[xx + 4] = dtyps2_extrm_v[cno][xx];
							}
							if (isWritefreq && (frqWriteset & dim_type)) {
								bingrp_v[idx].setId(0, 2, dim_ids);
								bingrp_v[idx].setExtremes(ad_extrm);
								bingrp_v[idx].setBinMids(0, 2 * bin_schemes_sum,
										bin2_ad_mids_v[idx]);
								bingrp_v[idx].setBinFreqs(
										xy2DHist->getFirstStep(),
										xy2DHist->getStepStride(),
										xy2DHist->getStepsEff(), 0,
										nbins_sum_sq, freq2_ad_obs_v[idx]);
							}
							dimentropy_v[idx].setId(0, 2, dim_ids);
							dimentropy_v[idx].setDimContri(0, nestimators, 0,
									nsteps, 0, n_schemes, entropy_xy2d_v[idx]);

						}
					}
				}
				auto _endcomp = chrono::high_resolution_clock::now();
				auto _durcomp = chrono::duration_cast<chrono::nanoseconds>(
						_endcomp - _endcomm1).count();
				ptaskset->curr_comp += chrono::nanoseconds(_durcomp);
				progressstate.entc.curr_comp += chrono::nanoseconds(_durcomp);

				auto _endcomm2 = chrono::high_resolution_clock::now();
				auto _durcomm2 = chrono::duration_cast<chrono::nanoseconds>(
						_endcomm2 - _endcomp).count();
				ptaskset->curr_comm += chrono::nanoseconds(_durcomm2);
				progressstate.entc.curr_comm += chrono::nanoseconds(_durcomm2);

				if (isWritefreq && (frqWriteset & dim_type)) {
					if (xy2DHist->writeRecords(bingrp_v) != 0) {
						mprinterr(
								"Error: Writing binned data for %s id(%d, %d)",
								str_dim_type.c_str(), 0, 0);
					}
				}
				ent_xy2d->writeRecords(dimentropy_v);

				auto _endwrt = chrono::high_resolution_clock::now();
				auto _durwrt = chrono::duration_cast<chrono::nanoseconds>(
						_endwrt - _endcomm2).count();
				ptaskset->curr_write += chrono::nanoseconds(_durwrt);
				progressstate.entc.curr_write += chrono::nanoseconds(_durwrt);
				ptaskset->current = chrono::high_resolution_clock::now();
				progressstate.entc.current = ptaskset->current;
				ptaskset->done_tasks += block_tasks;
				progressstate.entc.done_tasks += block_tasks;

				entc_tasks_done += block_tasks;
				if (entc_tasks_done >= entc_tasks_freq) {
					ofstream info_strm(
							inputs.getControl().getOutfilepath()
									+ inputs.getControl().getInfofile());
					info_strm << progressstate.toString();
					info_strm.close();
					entc_tasks_done %= entc_tasks_freq;
				}
			}
		}
	} // ang-dih block calc for d1 ends
	ptaskset->current = chrono::high_resolution_clock::now();
	ptaskset->cstt = ExecState::COMPLETED;
	progressstate.entc.current = ptaskset->current;
	int entc_workset_intval = static_cast<int>(inputs.getEntropy().getWorkSet());
	int entc_bd2d_intval = static_cast<int>(BATSet::BD2D);
	int entc_ad2d_intval = static_cast<int>(BATSet::AD2D);
	int entc_xy2d_intval = static_cast<int>(BATSet::XY2D);
	if ((entc_workset_intval < entc_bd2d_intval)
			&& (dim_type == BATSet::BA2D)) {
		progressstate.entc.cstt = ExecState::COMPLETED;
	} else if ((entc_workset_intval < entc_ad2d_intval)
			&& (dim_type == BATSet::BD2D)) {
		progressstate.entc.cstt = ExecState::COMPLETED;

	} else if ((entc_workset_intval <= entc_xy2d_intval)
			&& (dim_type == BATSet::AD2D)) {
		progressstate.entc.cstt = ExecState::COMPLETED;

	}
	ofstream info_strm(
			inputs.getControl().getOutfilepath()
					+ inputs.getControl().getInfofile());
	info_strm << progressstate.toString();
	info_strm.close();
	time_xy.stop();
	mprintf(
			"TIME: Histogram bin frequency calculation for %s-2D: %.4f seconds.\n",
			str_dim_type.c_str(), time_xy.total());
}

void EntropyCalculator::run_mpi(const int rank, const int numprocs,
		const int n_thread_perproc) {
	if (inputs.getControl().isCalcEntropy()) {
		BATSet workset = inputs.getEntropy().getWorkSet();

		if ((workset & BATSet::B1D) || (workset & BATSet::BB2D)) {
			run1d_mpi(BAT_t::BOND, rank, numprocs, n_thread_perproc);
		}
		if ((workset & BATSet::A1D) || (workset & BATSet::AA2D)) {
			run1d_mpi(BAT_t::ANGLE, rank, numprocs, n_thread_perproc);
		}
		if ((workset & BATSet::D1D) || (workset & BATSet::DD2D)) {
			run1d_mpi(BAT_t::DIHEDRAL, rank, numprocs, n_thread_perproc);
		}
		if ((workset & BATSet::BB2D)) {
			run2d_xx_mpi(BATSet::BB2D, rank, numprocs, n_thread_perproc);
		}
		if ((workset & BATSet::AA2D)) {
			run2d_xx_mpi(BATSet::AA2D, rank, numprocs, n_thread_perproc);
		}
		if ((workset & BATSet::DD2D)) {
			run2d_xx_mpi(BATSet::DD2D, rank, numprocs, n_thread_perproc);
		}
		if ((workset & BATSet::BA2D)) {
			run2d_xy_mpi(BATSet::BA2D, rank, numprocs, n_thread_perproc);
		}
		if ((workset & BATSet::BD2D)) {
			run2d_xy_mpi(BATSet::BD2D, rank, numprocs, n_thread_perproc);
		}
		if ((workset & BATSet::AD2D)) {
			run2d_xy_mpi(BATSet::AD2D, rank, numprocs, n_thread_perproc);
		}

	}
}

#endif
