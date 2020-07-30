/*
 * Discretizer.cpp
 *
 *  Created on: Sep 5, 2019
 *      Author: shailesh
 */

#include "Discretizer.h"

#ifndef USE_OMPMPI

void Discretizer::run_dtype(BAT_t dtype, std::vector<ull_int> &shuffle_index,
		const int num_threads) {

	Netcdf_BAT trajIn;
	NetcdfFile::NCTYPE fltype = trajIn.GetNetcdfConventions(fname.c_str());

	if (fltype == NetcdfFile::NCTYPE::NC_CENTRETRAJ) {
		trajIn.NC_openRead(fname);
		trajIn.setupFrameDim();
		trajIn.setupTime();
		CoordinateInfo cInfo;
		trajIn.setupRead(cInfo);

		nframes_trjs = trajIn.Ncframe();
		if (nfrm4rl2int > nframes_trjs) {
			mprinterr(
					"frames in traj is less than mentioned in input for discretization");
			exit(0);
		}

		std::vector<float> data(nframes_tot_eff);

		ProgTaskSet *dscrt_dtyp;
		std::string intTrajFile;
		std::string content;
		u_int n_dtyptasks = 0;

		std::stringstream vmkde_minima_log;
		switch (dtype) {
		case BAT_t::BOND:
			dscrt_dtyp = &(progressstate.dscrt_b);
			intTrajFile.assign(
					inputs.getControl().getOutfilepath() + "bin_bonds.nc");
			content.assign("B");
			n_dtyptasks = trajIn.Ncbond();
			break;
		case BAT_t::ANGLE:
			dscrt_dtyp = &(progressstate.dscrt_a);
			intTrajFile.assign(
					inputs.getControl().getOutfilepath() + "bin_angles.nc");
			content.assign("A");
			n_dtyptasks = trajIn.Ncangle();
			break;
		case BAT_t::DIHEDRAL:
			dscrt_dtyp = &(progressstate.dscrt_d);
			intTrajFile.assign(
					inputs.getControl().getOutfilepath() + "bin_torsions.nc");
			content.assign("T");
			n_dtyptasks = trajIn.Ncdihedral();
			break;
		default:
			cout << "ERROR:: Invalid BAT_t encountered" << endl;
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

		Netcdf_TrjInt *intDtypeTraj = new Netcdf_TrjInt(intTrajFile, content,
				n_dtyptasks, 0, 1, (ull_int) nfrm4rl2int,
				(hbin_t) bin_schemes.size(), bin_schemes);
		intDtypeTraj->NC_create(content + " Integer Vectors");

		int chached_dims = cache_dscrt_dims_per_cpu * num_threads;

		for (int bl_strat_id = 0; bl_strat_id < (int) n_dtyptasks;
				bl_strat_id += chached_dims) {
			int block_tasks =
					(bl_strat_id + chached_dims <= (int) n_dtyptasks) ?
							chached_dims : (n_dtyptasks - bl_strat_id);

			auto _startr = Clock::now();
			std::vector<std::vector<data_t>> dtypes_v(block_tasks,
					std::vector<data_t>(nfrm4rl2int));
			std::vector<std::vector<hbin_t>> dtypes_int_v(block_tasks,
					std::vector<hbin_t>(bin_schemes.size() * nfrm4rl2int));
			std::vector<std::vector<double>> dtypes_extrm_v(block_tasks,
					std::vector<double>(4, 0.0));

			std::vector<hbin_t> nfoundmins(block_tasks, 0);
			std::vector<std::vector<double>> min_poss(block_tasks,
					std::vector<double>(inputs.getVmkde().getNmaxconf(), 0.0));

			std::vector<CoordBAT> vecs_CoordBAT(block_tasks,
					CoordBAT(-1, bin_schemes.size(), nfrm4rl2int));

			int trjreadstatus = 0;

			switch (dtype) {
			case BAT_t::BOND:
				trjreadstatus = trajIn.readMultiBondFrames(bl_strat_id,
						block_tasks, 0, nfrm4rl2int, dtypes_v);
				break;
			case BAT_t::ANGLE:
				trjreadstatus = trajIn.readMultiAngleFrames(bl_strat_id,
						block_tasks, 0, nfrm4rl2int, dtypes_v);
				break;
			case BAT_t::DIHEDRAL:
				trjreadstatus = trajIn.readMultiDihedralFrames(bl_strat_id,
						block_tasks, 0, nfrm4rl2int, dtypes_v);
				break;
				//case BAT_t::NONE:
			default:
				std::cout << "Unrecognised BAT_t type" << std::endl;
				exit(0);
			}
			if (trjreadstatus != 0) {
				mprinterr(
						"Error: Writing binned data for bonds range id(%d, %d)",
						bl_strat_id, bl_strat_id + block_tasks);
				exit(0);
			}
			auto _endr = Clock::now();
			auto _dur_rd = std::chrono::duration_cast<ClockResolution>(
					_endr - _startr).count();
			dscrt_dtyp->curr_read += ClockResolution(_dur_rd);
			progressstate.dscrt.curr_read += ClockResolution(_dur_rd);

#pragma omp parallel for
			for (int idx = 0; idx < block_tasks; ++idx) {
				double min_d = Constants::DEFAULT_BAT_VALUE;
				double max_d = Constants::DEFAULT_BAT_VALUE;
				double avg_d = Constants::DEFAULT_BAT_VALUE;
				double phase_d = Constants::DEFAULT_BAT_VALUE;
				size_t did = bl_strat_id + idx;
				if (inputs.getBats().getPdfmethod() == PDFMethod::HISTOGRAM) {
					int returnvalue = 0;
					switch (dtype) {
					case BAT_t::BOND:
						returnvalue = Real2Int::bndReal2Int(did, bin_schemes,
								dtypes_v[idx], &min_d, &max_d, &avg_d,
								dtypes_int_v[idx]);

						break;
					case BAT_t::ANGLE:
						returnvalue = Real2Int::angReal2Int(did, bin_schemes,
								dtypes_v[idx], &min_d, &max_d, &avg_d,
								dtypes_int_v[idx]);
						break;
					case BAT_t::DIHEDRAL:
						returnvalue = Real2Int::dihReal2Int(did, bin_schemes,
								inputs.getHist().getNbinRef(),
								inputs.getBats().isOptimizedih(), dtypes_v[idx],
								&min_d, &max_d, &avg_d, &phase_d,
								dtypes_int_v[idx]);
						break;
					case BAT_t::NONE:
						break;
					}
					if (returnvalue != 0) {
						mprinterr(
								"Error: Converting Real2Int using HISTOGRAM for %s id(%d)",
								content.c_str(), did);
						exit(0);
					}
				} else if (inputs.getBats().getPdfmethod()
						== PDFMethod::vonMisesKDE) {
					if (Real2Int::getMinimaPositions(did, false,
							inputs.getVmkde().getNmaxconf(),
							inputs.getVmkde().getKappaValue(),
							inputs.getVmkde().getSdoNstep(),
							inputs.getVmkde().getSdoConvLimit(),
							inputs.getVmkde().getSdoNiter(), dtypes_v[idx],
							min_poss[idx], nfoundmins[idx],
							dtypes_int_v[idx])) {
						mprinterr(
								"Error: Converting Real2Int using vonMisesKDE for %s id(%d)",
								content.c_str(), did);
						exit(0);
					}
				}

				if (inputs.getBats().isShuffleFrames()) {
					Real2Int::shuffleArray(shuffle_index, dtypes_int_v[idx]);
				} else if (inputs.getBats().isShuffleDofs()) {
					unsigned seed =
							Clock::now().time_since_epoch().count();
					std::shuffle(shuffle_index.begin(), shuffle_index.end(),
							std::default_random_engine(seed));
					Real2Int::shuffleArray(shuffle_index, dtypes_int_v[idx]);
				}

				dtypes_extrm_v[idx][0] = min_d;
				dtypes_extrm_v[idx][1] = max_d;
				dtypes_extrm_v[idx][2] = avg_d;
				if (dtype == BAT_t::DIHEDRAL) {
					dtypes_extrm_v[idx][3] = phase_d;
				} else {
					// phase in not valid in bonds/angle hence assigned 0.0
					dtypes_extrm_v[idx][3] = 0.0;
				}

				(vecs_CoordBAT[idx]).setId(bl_strat_id + idx);
				(vecs_CoordBAT[idx]).setCoords(dtypes_extrm_v[idx][0],
						dtypes_extrm_v[idx][1], dtypes_extrm_v[idx][2],
						dtypes_extrm_v[idx][3], 0, (hbin_t) bin_schemes.size(),
						0, nfrm4rl2int - 1, dtypes_int_v[idx]);
			}
			auto _endcomm2 = Clock::now();
			auto _durcomm2 =
					std::chrono::duration_cast<ClockResolution>(
							_endcomm2 - _endr).count();

			dscrt_dtyp->curr_comp += ClockResolution(_durcomm2);
			progressstate.dscrt.curr_comp += ClockResolution(
					_durcomm2);

			if (intDtypeTraj->writeRecords(vecs_CoordBAT)) {
				mprinterr("Error: Writing binned data for %s range id(%d, %d)",
						content.c_str(), bl_strat_id,
						bl_strat_id + block_tasks);
				exit(0);
			}
			auto _endwrt = Clock::now();
			auto _durwrt = std::chrono::duration_cast<ClockResolution>(
					_endwrt - _endcomm2).count();
			dscrt_dtyp->curr_write += ClockResolution(_durwrt);
			progressstate.dscrt.curr_write += ClockResolution(_durwrt);

			dscrt_dtyp->current = Clock::now();
			progressstate.dscrt.current = dscrt_dtyp->current;
			dscrt_dtyp->done_tasks += block_tasks;
			progressstate.dscrt.done_tasks += block_tasks;
			dscrt_tasks_done += block_tasks;
			if ((dtype == BAT_t::DIHEDRAL)
					&& (inputs.getBats().getPdfmethod()
							== PDFMethod::vonMisesKDE)) {
				for (int idx = 0; idx < block_tasks; ++idx) {
					size_t did = bl_strat_id + idx;
					vmkde_minima_log << "torsion " << std::setw(6) << std::right
							<< did << " " << (int) nfoundmins[idx];
					for (auto itmpl = 0; itmpl < nfoundmins[idx]; ++itmpl)
						vmkde_minima_log << " " << std::fixed
								<< std::setprecision(8) << min_poss[idx][itmpl];
					vmkde_minima_log << ";" << std::endl;
				}
			}
			if (dscrt_tasks_done >= dscrt_tasks_freq) {
				std::ofstream vmkdelogfile(
						inputs.getControl().getOutfilepath()
								+ "/vmkde_torsions_minimas.out");
				vmkdelogfile << vmkde_minima_log.str();
				std::ofstream info_strm(
						inputs.getControl().getOutfilepath()
								+ inputs.getControl().getInfofile());
				info_strm << progressstate.toString();
				info_strm.close();
				dscrt_tasks_done %= dscrt_tasks_freq;
			}
		}
		dscrt_dtyp->current = Clock::now();
		dscrt_dtyp->cstt = ExecState::COMPLETED;
		int dscrt_active = 0;
		if (progressstate.dscrt_b.active)
			dscrt_active += 1;
		if (progressstate.dscrt_a.active)
			dscrt_active += 2;
		if (progressstate.dscrt_d.active)
			dscrt_active += 4;
		progressstate.dscrt.current = dscrt_dtyp->current;
		// The fall-through in case statement assume the order of execution BOND,ANGLE,DIHEDRAL
		switch (dscrt_active) {
		case 1:
			if (dtype == BAT_t::BOND) {
				progressstate.dscrt.cstt = ExecState::COMPLETED;
			}
			break;
		case 2:
		case 3:
			if (dtype == BAT_t::ANGLE) {
				progressstate.dscrt.cstt = ExecState::COMPLETED;
			}
			break;
		case 4:
		case 5:
		case 6:
		case 7:
			if (dtype == BAT_t::DIHEDRAL) {
				progressstate.dscrt.cstt = ExecState::COMPLETED;
			}
			break;
		default:

			break;
		}
		std::ofstream info_strm(
				inputs.getControl().getOutfilepath()
						+ inputs.getControl().getInfofile());
		info_strm << progressstate.toString();
		info_strm.close();
		trajIn.NC_close();
	}
}

void Discretizer::run(std::vector<ull_int> &shuffle_index,
		const int num_threads) {

	if (inputs.getControl().isDiscetize()) {
		{
			progressstate.dscrt.start =
					Clock::now();
			Timer time_real2int;
			time_real2int.start();
			mprintf("Trajectory Real to Integer conversion starts.\n");
			Netcdf_BAT trajIn;
			std::string fname(inputs.getControl().getInfilepath());
			fname.append(inputs.getBats().getFnames()[0]);

			NetcdfFile::NCTYPE fltype = trajIn.GetNetcdfConventions(
					fname.c_str());
			if (fltype == NetcdfFile::NCTYPE::NC_CENTRETRAJ) {
				if ((inputs.getEntropy().getWorkSet() & BATSet::B1D)
						|| (inputs.getEntropy().getWorkSet() & BATSet::BB2D)
						|| (inputs.getEntropy().getWorkSet() & BATSet::BA2D)
						|| (inputs.getEntropy().getWorkSet() & BATSet::BD2D)) {
					run_dtype(BAT_t::BOND, shuffle_index, num_threads);
				}
				if ((inputs.getEntropy().getWorkSet() & BATSet::A1D)
						|| (inputs.getEntropy().getWorkSet() & BATSet::BA2D)
						|| (inputs.getEntropy().getWorkSet() & BATSet::AA2D)
						|| (inputs.getEntropy().getWorkSet() & BATSet::AD2D)) {
					run_dtype(BAT_t::ANGLE, shuffle_index, num_threads);
				}
				if ((inputs.getEntropy().getWorkSet() & BATSet::D1D)
						|| (inputs.getEntropy().getWorkSet() & BATSet::BD2D)
						|| (inputs.getEntropy().getWorkSet() & BATSet::AD2D)
						|| (inputs.getEntropy().getWorkSet() & BATSet::DD2D)) {
					run_dtype(BAT_t::DIHEDRAL, shuffle_index, num_threads);
				}
			} else {
				std::cout << "Input Netcdf File is not of a valid CENTRETRAJ"
						<< std::endl;
				exit(0);
			}

			progressstate.dscrt.cstt = ExecState::COMPLETED;
			time_real2int.stop();
			mprintf("TIME: Discretization of trajectory time: %.4f seconds.\n",
					time_real2int.total());
		}
	}
}

#endif
