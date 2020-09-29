#include "EntropyCalculator.h"
#include "macrologger/macrologger.h"

#ifndef USE_OMPMPI

void EntropyCalculator::run1d(BAT_t bat_type, int num_threads)
{
	auto _inittime = Clock::now();
	Timer time_d1;
	time_d1.start();
	ProgTaskSet *ptaskset;
	std::string str_dim_type;
	u_int n_dim_eff = 0;
	Netcdf_TrjInt *intD1DTraj;
	std::string intD1DTrajFileName;
	Netcdf_HistUtil *d1DHist;
	std::string d1DHistFileName;
	Netcdf_EntContri *ent_dim1d;
	std::string ent_dim1dFileName;

	switch (bat_type)
	{
	case BAT_t::BOND:
		ptaskset = &(progressstate.entc_b1d);
		str_dim_type.assign("B");
		n_dim_eff = n_bnd_eff;
		intD1DTrajFileName.assign("bin_bonds.nc");
		d1DHistFileName.assign("hist_bnd-1d.nc");
		ent_dim1dFileName.assign("entcontri_bnd-1d.nc");
		break;
	case BAT_t::ANGLE:
		ptaskset = &(progressstate.entc_a1d);
		str_dim_type.assign("A");
		n_dim_eff = n_ang_eff;
		intD1DTrajFileName.assign("bin_angles.nc");
		d1DHistFileName.assign("hist_ang-1d.nc");
		ent_dim1dFileName.assign("entcontri_ang-1d.nc");
		break;
	case BAT_t::DIHEDRAL:
		ptaskset = &(progressstate.entc_d1d);
		str_dim_type.assign("D");
		n_dim_eff = n_dih_eff;
		intD1DTrajFileName.assign("bin_torsions.nc");
		d1DHistFileName.assign("hist_tor-1d.nc");
		ent_dim1dFileName.assign("entcontri_tor-1d.nc");
		break;
	default:
		cout << "Error:: Invalid 1D type found" << std::endl;
		exit(0);
	}

	std::vector<ull_int> dim_lens(1, n_dim_eff);
	intD1DTraj = new Netcdf_TrjInt(
		inputs.getControl().getOutfilepath() + intD1DTrajFileName,
		str_dim_type, n_dim_eff, startframe, strideframe,
		(ull_int)nframes_tot, n_schemes, bin_schemes);

	d1DHist = new Netcdf_HistUtil(
		inputs.getControl().getOutfilepath() + d1DHistFileName,
		str_dim_type, 1, 1, start_hist_step, nsteps, stride_hist_step,
		dim_lens, TensorType::FULL, n_schemes, bin_schemes);

	ent_dim1d = new Netcdf_EntContri(
		inputs.getControl().getOutfilepath() + ent_dim1dFileName,
		str_dim_type, 1, 1, nsteps, dim_lens, TensorType::FULL, bin_schemes,
		nestimators);

	ptaskset->start = _inittime;
	ptaskset->current = ptaskset->start;
	ptaskset->strt_read = ptaskset->curr_read = ptaskset->start;
	ptaskset->strt_comp = ptaskset->curr_comp = ptaskset->start;
	ptaskset->strt_comm = ptaskset->curr_comm = ptaskset->start;
	ptaskset->strt_write = ptaskset->curr_write = ptaskset->start;
	ptaskset->cstt = ExecState::RUNNING;
	if (progressstate.entc.cstt != ExecState::RUNNING)
	{
		progressstate.entc.cstt = ExecState::RUNNING;
		progressstate.entc.start = ptaskset->start;
		progressstate.entc.current = ptaskset->start;
	}

	ull_int nfrm_eff_entropy = 0;

	if (isWritefreq && (frqWriteset & bat_type))
	{
		d1DHist->NC_create((str_dim_type + string("-1D histograms")).c_str());
	}

	ent_dim1d->NC_create(
		(str_dim_type + string("-1D Entropy contributions")).c_str());

	auto _endinittime = Clock::now();
	auto _durinit = std::chrono::duration_cast<ClockResolution>(
						_endinittime - _inittime)
						.count();
	ptaskset->curr_comp += ClockResolution(_durinit);
	progressstate.entc.curr_comp += ClockResolution(_durinit);

	int chached_dims = cache_entc_dims_per_cpu * num_threads;
	for (u_int bl_strat_id = 0; bl_strat_id < n_dim_eff; bl_strat_id +=
														 chached_dims)
	{
		const int block_tasks =
			(bl_strat_id + chached_dims <= n_dim_eff) ? chached_dims : (n_dim_eff - bl_strat_id);
		auto _startr = Clock::now();
		std::vector<hbin_t *> d1ds_int_v(block_tasks);
		std::vector<std::vector<double>> d1ds_extrm_v(block_tasks,
													  std::vector<double>(4, 0.0));
		std::vector<CoordBAT> crd(block_tasks,
								  CoordBAT(0, n_schemes, nframes_tot_eff));
		std::vector<std::vector<ull_int>> freq_obs_v(block_tasks,
													 std::vector<ull_int>(nsteps * bin_schemes_sum));
		std::vector<std::vector<double>> bin_mids_v(block_tasks,
													std::vector<double>(bin_schemes_sum, 0.0));
		std::vector<std::vector<double>> entropy_v(block_tasks,
												   std::vector<double>(nestimators * nsteps * n_schemes, 0.0));

		std::vector<BinGroup> bingrp_v(block_tasks,
									   BinGroup(-1, d1DHist->getStepsEff(),
												d1DHist->getDimXBinSchemesSum()));

		std::vector<DimEntropy> dimentropy_v(block_tasks,
											 DimEntropy(-1, (hbin_t)nsteps, n_schemes, nestimators));

		for (int idx = 0; idx < block_tasks; ++idx)
		{
			const u_int dim_id = inputs.getSubset().getBATTypeId(bat_type,
																 bl_strat_id + idx);
			crd[idx].setId(dim_id);
			if (intD1DTraj->readCoords(dim_id, 0, (hbin_t)n_schemes, crd[idx]) != 0)
			{
				LOG_ERROR("Reading int %s id(%u)", str_dim_type.c_str(), dim_id);
			}
			crd[idx].getCoords(&d1ds_extrm_v[idx][0], &d1ds_extrm_v[idx][1],
							   &d1ds_extrm_v[idx][2], &d1ds_extrm_v[idx][3], &n_schemes,
							   &nfrm_eff_entropy, &d1ds_int_v[idx]);
		}
		auto _endr = Clock::now();
		auto _dur_rd = std::chrono::duration_cast<ClockResolution>(
						   _endr - _startr)
						   .count();
		ptaskset->curr_read += ClockResolution(_dur_rd);
		progressstate.entc.curr_read += ClockResolution(_dur_rd);

		auto _endcomm1 = Clock::now();
		auto _durcomm1 = std::chrono::duration_cast<ClockResolution>(
							 _endcomm1 - _endr)
							 .count();
		ptaskset->curr_comm += ClockResolution(_durcomm1);
		progressstate.entc.curr_comm += ClockResolution(_durcomm1);

#pragma omp parallel for
		for (int idx = 0; idx < block_tasks; ++idx)
		{
			size_t did = bl_strat_id + idx;
			if (binData1D(nsteps, bin_schemes, nframes_tot_eff,
						  d1ds_extrm_v[idx][0], d1ds_extrm_v[idx][1], d1ds_int_v[idx],
						  freq_obs_v[idx], bin_mids_v[idx]) != 0)
			{
				LOG_ERROR("Binning %s data bond id(%u)",
						  str_dim_type.c_str(), bl_strat_id + idx);
			}
			else
			{
				entropy1D(bat_type, did, inputs.getEstimators(), nsteps,
						  step_size, d1ds_extrm_v[idx][0], d1ds_extrm_v[idx][1],
						  inputs.getEntropy().isJacobian(), isKDE, bin_schemes,
						  freq_obs_v[idx], entropy_v[idx]);

				const u_int dim_id = inputs.getSubset().getBATTypeId(bat_type,
																	 did);
				if (isWritefreq && (frqWriteset & BATSet::B1D))
				{
					bingrp_v[idx].setId(dim_id);
					bingrp_v[idx].setExtremes(d1ds_extrm_v[idx]);
					bingrp_v[idx].setBinMids(0, bin_schemes_sum,
											 bin_mids_v[idx]);
					bingrp_v[idx].setBinFreqs(d1DHist->getFirstStep(),
											  d1DHist->getStepStride(), d1DHist->getStepsEff(), 0,
											  bin_schemes_sum, freq_obs_v[idx]);
				}
				dimentropy_v[idx].setId(dim_id);
				dimentropy_v[idx].setDimContri(0, nestimators, 0, nsteps, 0,
											   n_schemes, entropy_v[idx]);
			}
		}

		auto _endcomp = Clock::now();
		auto _durcomp = std::chrono::duration_cast<ClockResolution>(
							_endcomp - _endcomm1)
							.count();
		ptaskset->curr_comp += ClockResolution(_durcomp);
		progressstate.entc.curr_comp += ClockResolution(_durcomp);

		auto _endcomm2 = Clock::now();
		auto _durcomm2 = std::chrono::duration_cast<ClockResolution>(
							 _endcomm2 - _endcomp)
							 .count();
		ptaskset->curr_comm += ClockResolution(_durcomm2);
		progressstate.entc.curr_comm += ClockResolution(_durcomm2);

		if (isWritefreq && (frqWriteset & bat_type))
		{
			if (d1DHist->writeRecords(bingrp_v) != 0)
			{
				LOG_ERROR("Writing binned data for %s id(%u, %u)",
						  str_dim_type.c_str(), bl_strat_id,
						  bl_strat_id + block_tasks);
			}
		}
		ent_dim1d->writeRecords(dimentropy_v);

		auto _endwrt = Clock::now();
		auto _durwrt = std::chrono::duration_cast<ClockResolution>(
						   _endwrt - _endcomm2)
						   .count();
		ptaskset->curr_write += ClockResolution(_durwrt);
		progressstate.entc.curr_write += ClockResolution(_durwrt);
		ptaskset->current = Clock::now();
		progressstate.entc.current = ptaskset->current;
		ptaskset->done_tasks += block_tasks;
		progressstate.entc.done_tasks += block_tasks;

		entc_tasks_done += block_tasks;
		if (entc_tasks_done >= entc_tasks_freq)
		{
			std::ofstream info_strm(
				inputs.getControl().getOutfilepath() + inputs.getControl().getInfofile());
			info_strm << progressstate.toString();
			info_strm.close();
			entc_tasks_done %= entc_tasks_freq;
		}
	}
	ptaskset->current = Clock::now();
	ptaskset->cstt = ExecState::COMPLETED;
	progressstate.entc.current = ptaskset->current;
	time_d1.stop();
	mprintf("TIME: Histogram bin frequency calculation for %s: %.4f seconds.\n",
			str_dim_type.c_str(), time_d1.total());

	int entc_active = static_cast<int>(inputs.getEntropy().getWorkSet());
	// The fall-through in case statement assume the order of execution BOND,ANGLE,DIHEDRAL
	switch (entc_active)
	{
	case 1:
		if (bat_type == BAT_t::BOND)
		{
			progressstate.entc.cstt = ExecState::COMPLETED;
		}
		break;
	case 2:
	case 3:
		if (bat_type == BAT_t::ANGLE)
		{
			progressstate.entc.cstt = ExecState::COMPLETED;
		}
		break;
	case 4:
	case 5:
	case 6:
	case 7:
		if (bat_type == BAT_t::DIHEDRAL)
		{
			progressstate.entc.cstt = ExecState::COMPLETED;
		}
		break;
	default:

		break;
	}
	std::ofstream info_strm(
		inputs.getControl().getOutfilepath() + inputs.getControl().getInfofile());
	info_strm << progressstate.toString();
	info_strm.close();
}

void EntropyCalculator::run2d_xx(BATSet dim_type, int num_threads)
{
	auto _inittime = Clock::now();
	Timer time_ba;
	time_ba.start();
	ProgTaskSet *ptaskset;
	std::string str_dim_type;
	u_int n_dim_eff = 0;
	u_int num_xx2d;
	Netcdf_TrjInt *intXX2DTraj;
	std::string intXX1DTrajFileName;
	Netcdf_HistUtil *xx2DHist;
	std::string xx2DHistFileName;
	Netcdf_EntContri *ent_xx2d;
	std::string ent_xx2dFileName;

	std::vector<u_int> dimTypeKeys;
	std::vector<u_int> dimtypes_v;
	BAT_t dtype1st = BAT_t::NONE;
	switch (dim_type)
	{
	case BATSet::BB2D:
		dtype1st = BAT_t::BOND;
		ptaskset = &(progressstate.entc_b2d);
		str_dim_type.assign("B/B");
		n_dim_eff = n_bnd_eff;
		intXX1DTrajFileName.assign("bin_bonds.nc");
		xx2DHistFileName.assign("hist_bnd-2d.nc");
		ent_xx2dFileName.assign("entcontri_bnd-2d.nc");
		inputs.getNeighbors().bondKeys(dimTypeKeys);
		for (size_t d1 = 0; d1 < dimTypeKeys.size(); ++d1)
		{
			num_xx2d += inputs.getNeighbors().bondNeighSize(dimTypeKeys[d1]);
		}
		std::copy(inputs.getSubset().getBonds().begin(),
				  inputs.getSubset().getBonds().end(),
				  std::back_inserter(dimtypes_v));
		break;
	case BATSet::AA2D:
		dtype1st = BAT_t::ANGLE;
		ptaskset = &(progressstate.entc_a2d);
		str_dim_type.assign("A/A");
		n_dim_eff = n_ang_eff;
		intXX1DTrajFileName.assign("bin_angles.nc");
		xx2DHistFileName.assign("hist_ang-2d.nc");
		ent_xx2dFileName.assign("entcontri_ang-2d.nc");
		inputs.getNeighbors().angleKeys(dimTypeKeys);
		for (size_t d1 = 0; d1 < dimTypeKeys.size(); ++d1)
		{
			num_xx2d += inputs.getNeighbors().angleNeighSize(dimTypeKeys[d1]);
		}
		std::copy(inputs.getSubset().getAngles().begin(),
				  inputs.getSubset().getAngles().end(),
				  std::back_inserter(dimtypes_v));
		break;
	case BATSet::DD2D:
		dtype1st = BAT_t::DIHEDRAL;
		ptaskset = &(progressstate.entc_d2d);
		str_dim_type.assign("D/D");
		n_dim_eff = n_dih_eff;
		intXX1DTrajFileName.assign("bin_torsions.nc");
		xx2DHistFileName.assign("hist_tor-2d.nc");
		ent_xx2dFileName.assign("entcontri_tor-2d.nc");
		inputs.getNeighbors().torsionKeys(dimTypeKeys);
		for (size_t d1 = 0; d1 < dimTypeKeys.size(); ++d1)
		{
			num_xx2d += inputs.getNeighbors().torsionNeighSize(dimTypeKeys[d1]);
		}
		std::copy(inputs.getSubset().getTorsions().begin(),
				  inputs.getSubset().getTorsions().end(),
				  std::back_inserter(dimtypes_v));
		break;
	default:
		cout << "Error:: Invalid xx2D type found" << std::endl;
		exit(0);
	}

	std::vector<ull_int> dim_lens(2, n_dim_eff);
	intXX2DTraj = new Netcdf_TrjInt(
		inputs.getControl().getOutfilepath() + intXX1DTrajFileName,
		str_dim_type, n_dim_eff, startframe, strideframe,
		(ull_int)nframes_tot, n_schemes, bin_schemes);

	xx2DHist = new Netcdf_HistUtil(
		inputs.getControl().getOutfilepath() + xx2DHistFileName,
		str_dim_type, 2, 2, start_hist_step, nsteps, stride_hist_step,
		dim_lens, TensorType::FULL, n_schemes, bin_schemes);

	if (isWritefreq && (frqWriteset & dim_type))
	{
		xx2DHist->NC_create(str_dim_type + " 2-D histograms");
	}

	ent_xx2d = new Netcdf_EntContri(
		inputs.getControl().getOutfilepath() + ent_xx2dFileName,
		str_dim_type, 2, 2, nsteps, dim_lens, TensorType::UPPER,
		bin_schemes, nestimators);

	ptaskset->start = _inittime;
	ptaskset->current = ptaskset->start;
	ptaskset->strt_read = ptaskset->curr_read = ptaskset->start;
	ptaskset->strt_comp = ptaskset->curr_comp = ptaskset->start;
	ptaskset->strt_comm = ptaskset->curr_comm = ptaskset->start;
	ptaskset->strt_write = ptaskset->curr_write = ptaskset->start;
	ptaskset->cstt = ExecState::RUNNING;
	if (progressstate.entc.cstt != ExecState::RUNNING)
	{
		progressstate.entc.cstt = ExecState::RUNNING;
		progressstate.entc.start = ptaskset->start;
		progressstate.entc.current = ptaskset->start;
	}

	ull_int nfrm_eff_entropy = 0;
	/*********************** MASTER PROCESS *******************
	 * 2-D ENTROPY estimation for: BOND/ANGLE
	 *********************************************************/

	ent_xx2d->NC_create(str_dim_type + " 2-D Entropy contributions");
	auto _endinittime = Clock::now();
	auto _durinit = std::chrono::duration_cast<ClockResolution>(
						_endinittime - _inittime)
						.count();
	ptaskset->curr_comp += ClockResolution(_durinit);
	progressstate.entc.curr_comp += ClockResolution(_durinit);

	std::vector<u_int> block_boundry(5);
	const u_int num_dimTypeKeys = dimTypeKeys.size();

	int chached_dims = cache_entc_dims_per_cpu * num_threads / 2;
	for (u_int bi = 0; bi < num_dimTypeKeys; bi += chached_dims)
	{
		const int bl_rows =
			(bi + chached_dims <= num_dimTypeKeys) ? chached_dims : (num_dimTypeKeys - bi);

		std::vector<hbin_t *> dtyps1_int_v(bl_rows);
		std::vector<std::vector<double>> dtyps1_extrm_v(bl_rows,
														std::vector<double>(4, 0.0));
		std::vector<CoordBAT> crd_rows(bl_rows,
									   CoordBAT(0, n_schemes, nframes_tot_eff));

		// Fill Cache with rows
		auto _startr = Clock::now();
		for (auto rno = 0; rno < bl_rows; ++rno)
		{
			const u_int dtyp_id1 = dimTypeKeys[bi + rno];
			crd_rows[rno].setId(dtyp_id1);
			if (intXX2DTraj->readCoords(dtyp_id1, 0, n_schemes, crd_rows[rno]) != 0)
			{
				LOG_ERROR("Error: Reading int %s for row id(%u)",
						  str_dim_type.c_str(), dtyp_id1);
			}
			crd_rows[rno].getCoords(&dtyps1_extrm_v[rno][0],
									&dtyps1_extrm_v[rno][1], &dtyps1_extrm_v[rno][2],
									&dtyps1_extrm_v[rno][3], &n_schemes, &nfrm_eff_entropy,
									&dtyps1_int_v[rno]);
		}
		auto _endr = Clock::now();
		auto _dur_rd = std::chrono::duration_cast<ClockResolution>(
						   _endr - _startr)
						   .count();
		ptaskset->curr_read += ClockResolution(_dur_rd);
		progressstate.entc.curr_read += ClockResolution(_dur_rd);

		for (u_int bj = 0; bj < n_dim_eff; bj += chached_dims)
		{
			const int bl_cols =
				(bj + chached_dims <= n_dim_eff) ? chached_dims : (n_dim_eff - bj);

			// Calc
			u_int block_tasks = 0, col_idx = 0;
			std::map<u_int, u_int> read_cols_dims;
			std::map<std::pair<u_int, u_int>, u_int> id2index;

			for (auto ri = bi; ri < bi + bl_rows; ++ri)
			{
				const u_int dtyp_id1 = dimTypeKeys[ri];
				const std::vector<u_int> neigh =
					inputs.getNeighbors().getDimTypeNeighs(dim_type,
														   dtyp_id1);
				for (auto const cj : neigh)
				{
					if (cj >= dimtypes_v[bj] && cj <= dimtypes_v[bj + bl_cols - 1])
					{
						if (!read_cols_dims.count(cj))
						{
							read_cols_dims[cj] = col_idx;
							++col_idx;
						}
						id2index[std::make_pair(dtyp_id1, cj)] = block_tasks;
						++block_tasks;
					}
				}
			}
			if (block_tasks)
			{
				u_int n_cache_cols = read_cols_dims.size();
				std::vector<hbin_t *> dtyps2_int_v(n_cache_cols);
				std::vector<std::vector<double>> bnds2_extrm_v(n_cache_cols,
															   std::vector<double>(4, 0.0));
				std::vector<CoordBAT> crd_cols(n_cache_cols,
											   CoordBAT(0, n_schemes, nframes_tot_eff));

				// read and cache cols
				auto _startr2 = Clock::now();
				for (std::map<u_int, u_int>::iterator it =
						 read_cols_dims.begin();
					 it != read_cols_dims.end();
					 ++it)
				{
					//const u_int dtyp_id2 = it->first;
					crd_cols[it->second].setId(it->first);
					if (intXX2DTraj->readCoords(it->first, 0, n_schemes,
												crd_cols[it->second]) != 0)
					{
						LOG_ERROR("Reading int of %s for column id(%u)",
								  str_dim_type.c_str(), it->first);
					}
					crd_cols[it->second].getCoords(
						&bnds2_extrm_v[it->second][0],
						&bnds2_extrm_v[it->second][1],
						&bnds2_extrm_v[it->second][2],
						&bnds2_extrm_v[it->second][3], &n_schemes,
						&nfrm_eff_entropy, &dtyps2_int_v[it->second]);
				}
				auto _endr2 = Clock::now();
				auto _dur2_rd = std::chrono::duration_cast<
									ClockResolution>(_endr2 - _startr2)
									.count();
				ptaskset->curr_read += ClockResolution(_dur2_rd);
				progressstate.entc.curr_read += ClockResolution(
					_dur2_rd);

				std::vector<std::vector<ull_int>> freq2_bb_obs_v(block_tasks,
																 std::vector<ull_int>(nsteps * nbins_sum_sq));
				std::vector<std::vector<double>> bin2_bb_mids_v(block_tasks,
																std::vector<double>(2 * bin_schemes_sum, 0.0));

				std::vector<std::vector<double>> entropy_bb_v(block_tasks,
															  std::vector<double>(nestimators * nsteps * n_schemes,
																				  0.0));

				std::vector<u_int> dummy_ids(2, -1);
				std::vector<BinGroup> bingrp_v(block_tasks,
											   BinGroup(dummy_ids, (hbin_t)nsteps,
														2 * bin_schemes_sum, nbins_sum_sq));
				std::vector<DimEntropy> dimentropy_v(block_tasks,
													 DimEntropy(dummy_ids, (hbin_t)nsteps, n_schemes,
																nestimators));

#pragma omp parallel for
				for (auto ri = bi; ri < bi + bl_rows; ++ri)
				{
					auto rno = ri - bi;
					const u_int bnd_id1 = dimTypeKeys[ri];
					const std::vector<u_int> neigh =
						inputs.getNeighbors().getDimTypeNeighs(dim_type,
															   bnd_id1);

					for (auto const cj : neigh)
					{
						if (cj >= dimtypes_v[bj] && cj <= dimtypes_v[bj + bl_cols - 1])
						{
							u_int idx = id2index[std::make_pair(bnd_id1, cj)];
							auto cno = read_cols_dims[cj];
							if (binData2D(nsteps, bin_schemes, nframes_tot_eff,
										  dtyps1_extrm_v[rno][0],
										  dtyps1_extrm_v[rno][1],
										  bnds2_extrm_v[cno][0],
										  bnds2_extrm_v[cno][1], dtyps1_int_v[rno],
										  dtyps2_int_v[cno], freq2_bb_obs_v[idx],
										  bin2_bb_mids_v[idx]) != 0)
							{
								LOG_ERROR("Bining data %s row=%u, col=%u]",
										  str_dim_type.c_str(), bnd_id1, cj);
							}
							entropy2D(dtype1st, dtype1st, bnd_id1, cj,
									  inputs.getEstimators(), nsteps, step_size,
									  dtyps1_extrm_v[rno][0],
									  dtyps1_extrm_v[rno][1],
									  bnds2_extrm_v[cno][0],
									  bnds2_extrm_v[cno][1],
									  inputs.getEntropy().isJacobian(), isKDE,
									  bin_schemes, freq2_bb_obs_v[idx],
									  entropy_bb_v[idx]);

							std::vector<u_int> dim_ids(2);
							dim_ids[0] = bnd_id1;
							dim_ids[1] = cj;
							std::vector<double> bb_extrm(8, 0.0);
							for (int xx = 0; xx < 4; ++xx)
							{
								bb_extrm[xx] = dtyps1_extrm_v[rno][xx];
								bb_extrm[xx + 4] = bnds2_extrm_v[cno][xx];
							}
							if (isWritefreq && (frqWriteset & BATSet::BB2D))
							{
								bingrp_v[idx].setId(0, 2, dim_ids);
								bingrp_v[idx].setExtremes(bb_extrm);
								bingrp_v[idx].setBinMids(0, 2 * bin_schemes_sum,
														 bin2_bb_mids_v[idx]);
								bingrp_v[idx].setBinFreqs(
									xx2DHist->getFirstStep(),
									xx2DHist->getStepStride(),
									xx2DHist->getStepsEff(), 0,
									nbins_sum_sq, freq2_bb_obs_v[idx]);
							}
							dimentropy_v[idx].setId(0, 2, dim_ids);
							dimentropy_v[idx].setDimContri(0, nestimators, 0,
														   nsteps, 0, n_schemes, entropy_bb_v[idx]);
						}
					}
				}
				auto _endcomp = Clock::now();
				auto _durcomp = std::chrono::duration_cast<
									ClockResolution>(_endcomp - _endr2)
									.count();
				ptaskset->curr_comp += ClockResolution(_durcomp);
				progressstate.entc.curr_comp += ClockResolution(
					_durcomp);

				auto _endcomm2 = Clock::now();
				auto _durcomm2 = std::chrono::duration_cast<
									 ClockResolution>(_endcomm2 - _endcomp)
									 .count();
				ptaskset->curr_comm += ClockResolution(_durcomm2);
				progressstate.entc.curr_comm += ClockResolution(
					_durcomm2);

				if (isWritefreq && (frqWriteset & BATSet::BB2D))
				{
					if (xx2DHist->writeRecords(bingrp_v) != 0)
					{
						LOG_ERROR(
							"Writing histogram bin data for %s id rows[%u, %u] cols[%u, %u]",
							str_dim_type.c_str(), bi, bi + bl_rows, bj, bj + bl_cols);
					}
				}
				ent_xx2d->writeRecords(dimentropy_v);

				auto _endwrt = Clock::now();
				auto _durwrt = std::chrono::duration_cast<
								   ClockResolution>(_endwrt - _endcomm2)
								   .count();
				ptaskset->curr_write += ClockResolution(_durwrt);
				progressstate.entc.curr_write += ClockResolution(
					_durwrt);
				ptaskset->current = Clock::now();
				progressstate.entc.current = ptaskset->current;
				ptaskset->done_tasks += block_tasks;
				progressstate.entc.done_tasks += block_tasks;

				entc_tasks_done += block_tasks;
				if (entc_tasks_done >= entc_tasks_freq)
				{
					std::ofstream info_strm(
						inputs.getControl().getOutfilepath() + inputs.getControl().getInfofile());
					info_strm << progressstate.toString();
					info_strm.close();
					entc_tasks_done %= entc_tasks_freq;
				}
			}
		}
	} // dih-dih block calc for d1 ends

	ptaskset->current = Clock::now();
	ptaskset->cstt = ExecState::COMPLETED;
	progressstate.entc.current = ptaskset->current;
	int entc_workset_intval = static_cast<int>(inputs.getEntropy().getWorkSet());
	int entc_aa2d_intval = static_cast<int>(BATSet::AA2D);
	int entc_dd2d_intval = static_cast<int>(BATSet::DD2D);
	int entc_xx2d_intval = static_cast<int>(BATSet::XX2D);
	if ((entc_workset_intval < entc_aa2d_intval) && (dim_type == BATSet::BB2D))
	{
		progressstate.entc.cstt = ExecState::COMPLETED;
	}
	else if ((entc_workset_intval < entc_dd2d_intval) && (dim_type == BATSet::AA2D))
	{
		progressstate.entc.cstt = ExecState::COMPLETED;
	}
	else if ((entc_workset_intval <= entc_xx2d_intval) && (dim_type == BATSet::DD2D))
	{
		progressstate.entc.cstt = ExecState::COMPLETED;
	}
	std::ofstream info_strm(
		inputs.getControl().getOutfilepath() + inputs.getControl().getInfofile());
	info_strm << progressstate.toString();
	info_strm.close();
	time_ba.stop();
	mprintf(
		"TIME: Histogram bin frequency calculation for %s-2D: %.4f seconds.\n",
		str_dim_type.c_str(), time_ba.total());
}

void EntropyCalculator::run2d_xy(BATSet dim_type, int num_threads)
{
	auto _inittime = Clock::now();
	Timer time_xy;
	time_xy.start();
	ProgTaskSet *ptaskset;
	std::string str_dim_type;
	u_int n_dim_x_eff = 0;
	u_int n_dim_y_eff = 0;
	u_int num_xy2d;
	Netcdf_TrjInt *intxofXY2DTraj;
	Netcdf_TrjInt *intyofXY2DTraj;
	std::string intxofXY2DTrajFileName;
	std::string intyofXY2DTrajFileName;
	Netcdf_HistUtil *xy2DHist;
	std::string xy2DHistFileName;
	Netcdf_EntContri *ent_xy2d;
	std::string ent_xy2dFileName;

	std::vector<u_int> dimTypeKeys;
	std::vector<u_int> dimtypes_v;
	BAT_t dtype1st = BAT_t::NONE;
	BAT_t dtype2nd = BAT_t::NONE;
	switch (dim_type)
	{
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
		for (size_t d1 = 0; d1 < dimTypeKeys.size(); ++d1)
		{
			num_xy2d += inputs.getNeighbors().baNeighSize(dimTypeKeys[d1]);
		}

		std::copy(inputs.getSubset().getAngles().begin(),
				  inputs.getSubset().getAngles().end(),
				  std::back_inserter(dimtypes_v));
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
		for (size_t d1 = 0; d1 < dimTypeKeys.size(); ++d1)
		{
			num_xy2d += inputs.getNeighbors().bdNeighSize(dimTypeKeys[d1]);
		}

		std::copy(inputs.getSubset().getTorsions().begin(),
				  inputs.getSubset().getTorsions().end(),
				  std::back_inserter(dimtypes_v));
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
		for (size_t d1 = 0; d1 < dimTypeKeys.size(); ++d1)
		{
			num_xy2d += inputs.getNeighbors().adNeighSize(dimTypeKeys[d1]);
		}
		std::copy(inputs.getSubset().getTorsions().begin(),
				  inputs.getSubset().getTorsions().end(),
				  std::back_inserter(dimtypes_v));
		break;
	default:
		cout << "Error:: Invalid xy2D type found" << std::endl;
		exit(0);
	}

	std::vector<ull_int> dim_lens(2, 0);
	dim_lens[0] = n_dim_x_eff;
	dim_lens[1] = n_dim_y_eff;

	intxofXY2DTraj = new Netcdf_TrjInt(
		inputs.getControl().getOutfilepath() + intxofXY2DTrajFileName,
		str_dim_type, n_dim_x_eff, startframe, strideframe,
		(ull_int)nframes_tot, n_schemes, bin_schemes);

	intyofXY2DTraj = new Netcdf_TrjInt(
		inputs.getControl().getOutfilepath() + intyofXY2DTrajFileName,
		str_dim_type, n_dim_y_eff, startframe, strideframe,
		(ull_int)nframes_tot, n_schemes, bin_schemes);

	xy2DHist = new Netcdf_HistUtil(
		inputs.getControl().getOutfilepath() + xy2DHistFileName,
		str_dim_type, 2, 2, start_hist_step, nsteps, stride_hist_step,
		dim_lens, TensorType::UPPER, n_schemes, bin_schemes);

	if (isWritefreq && (frqWriteset & dim_type))
	{
		xy2DHist->NC_create(str_dim_type + " 2-D histograms");
	}

	ent_xy2d = new Netcdf_EntContri(
		inputs.getControl().getOutfilepath() + ent_xy2dFileName,
		str_dim_type, 2, 2, nsteps, dim_lens, TensorType::UPPER,
		bin_schemes, nestimators);

	ent_xy2d->NC_create(str_dim_type + "-2D Entropy contributions");

	ptaskset->start = _inittime;
	ptaskset->current = ptaskset->start;
	ptaskset->strt_read = ptaskset->curr_read = ptaskset->start;
	ptaskset->strt_comp = ptaskset->curr_comp = ptaskset->start;
	ptaskset->strt_comm = ptaskset->curr_comm = ptaskset->start;
	ptaskset->strt_write = ptaskset->curr_write = ptaskset->start;
	ptaskset->cstt = ExecState::RUNNING;
	if (progressstate.entc.cstt != ExecState::RUNNING)
	{
		progressstate.entc.cstt = ExecState::RUNNING;
		progressstate.entc.start = ptaskset->start;
		progressstate.entc.current = ptaskset->start;
	}
	auto _endinittime = Clock::now();
	auto _durinit = std::chrono::duration_cast<ClockResolution>(
						_endinittime - _inittime)
						.count();
	ptaskset->curr_comp += ClockResolution(_durinit);
	progressstate.entc.curr_comp += ClockResolution(_durinit);
	ull_int nfrm_eff_entropy = 0;
	/*********************** MASTER PROCESS *******************
	 * 2-D ENTROPY estimation for: BOND/ANGLE
	 *********************************************************/

	std::vector<u_int> block_boundry(5);
	const u_int num_xy2DKeys = dimTypeKeys.size();

	int chached_dims = cache_entc_dims_per_cpu * num_threads / 2;
	for (u_int bi = 0; bi < num_xy2DKeys; bi += chached_dims)
	{
		const int bl_rows =
			(bi + chached_dims <= num_xy2DKeys) ? chached_dims : (num_xy2DKeys - bi);

		std::vector<hbin_t *> dtyps1_int_v(bl_rows);
		std::vector<std::vector<double>> dtyps1_extrm_v(bl_rows,
														std::vector<double>(4, 0.0));
		std::vector<CoordBAT> crd_rows(bl_rows,
									   CoordBAT(0, n_schemes, nframes_tot_eff));

		// Fill Cache with rows
		auto _startr = Clock::now();
		for (auto rno = 0; rno < bl_rows; ++rno)
		{
			const u_int dtyp_id1 = dimTypeKeys[bi + rno];
			crd_rows[rno].setId(dtyp_id1);
			if (intxofXY2DTraj->readCoords(dtyp_id1, 0, n_schemes,
										   crd_rows[rno]) != 0)
			{
				LOG_ERROR("Reading %s int vector id(%u)", str_dim_type.c_str(), dtyp_id1);
			}
			crd_rows[rno].getCoords(&dtyps1_extrm_v[rno][0],
									&dtyps1_extrm_v[rno][1], &dtyps1_extrm_v[rno][2],
									&dtyps1_extrm_v[rno][3], &n_schemes, &nfrm_eff_entropy,
									&dtyps1_int_v[rno]);
		}
		auto _endr = Clock::now();
		auto _dur_rd = std::chrono::duration_cast<ClockResolution>(
						   _endr - _startr)
						   .count();
		ptaskset->curr_read += ClockResolution(_dur_rd);
		progressstate.entc.curr_read += ClockResolution(_dur_rd);

		for (u_int bj = 0; bj < n_dim_y_eff; bj += chached_dims)
		{
			const int bl_cols =
				(bj + chached_dims <= n_dim_y_eff) ? chached_dims : (n_dim_y_eff - bj);

			// Calc
			u_int block_tasks = 0, col_idx = 0;
			std::map<u_int, u_int> read_cols_dims;
			std::map<std::pair<u_int, u_int>, u_int> id2index;

			for (auto ri = bi; ri < bi + bl_rows; ++ri)
			{
				const u_int dtyp_id1 = dimTypeKeys[ri];
				const std::vector<u_int> neigh =
					inputs.getNeighbors().getDimTypeNeighs(dim_type,
														   dtyp_id1);
				for (auto const cj : neigh)
				{
					if (cj >= dimtypes_v[bj] && cj <= dimtypes_v[bj + bl_cols - 1])
					{
						if (!read_cols_dims.count(cj))
						{
							read_cols_dims[cj] = col_idx;
							++col_idx;
						}
						id2index[std::make_pair(dtyp_id1, cj)] = block_tasks;
						++block_tasks;
					}
				}
			}
			if (block_tasks)
			{
				u_int n_cache_cols = read_cols_dims.size();
				std::vector<hbin_t *> dtyps2_int_v(n_cache_cols);
				std::vector<std::vector<double>> dtyps2_extrm_v(n_cache_cols,
																std::vector<double>(4, 0.0));
				std::vector<CoordBAT> crd_cols(n_cache_cols,
											   CoordBAT(0, n_schemes, nframes_tot_eff));

				// read and cache cols
				auto _startr2 = Clock::now();
				for (std::map<u_int, u_int>::iterator it =
						 read_cols_dims.begin();
					 it != read_cols_dims.end();
					 ++it)
				{
					crd_cols[it->second].setId(it->first);
					if (intyofXY2DTraj->readCoords(it->first, 0, n_schemes,
												   crd_cols[it->second]) != 0)
					{
						LOG_ERROR("Reading %s int vector  id(%u)", str_dim_type.c_str(), it->first);
					}
					crd_cols[it->second].getCoords(
						&dtyps2_extrm_v[it->second][0],
						&dtyps2_extrm_v[it->second][1],
						&dtyps2_extrm_v[it->second][2],
						&dtyps2_extrm_v[it->second][3], &n_schemes,
						&nfrm_eff_entropy, &dtyps2_int_v[it->second]);
				}
				auto _endr2 = Clock::now();
				auto _dur2_rd = std::chrono::duration_cast<
									ClockResolution>(_endr2 - _startr2)
									.count();
				ptaskset->curr_read += ClockResolution(_dur2_rd);
				progressstate.entc.curr_read += ClockResolution(
					_dur2_rd);

				std::vector<std::vector<ull_int>> freq2_ad_obs_v(block_tasks,
																 std::vector<ull_int>(nsteps * nbins_sum_sq));
				std::vector<std::vector<double>> bin2_ad_mids_v(block_tasks,
																std::vector<double>(2 * bin_schemes_sum, 0.0));

				std::vector<std::vector<double>> entropy_xy2d_v(block_tasks,
																std::vector<double>(nestimators * nsteps * n_schemes,
																					0.0));

				std::vector<u_int> dummy_ids(2, -1);
				std::vector<BinGroup> bingrp_v(block_tasks,
											   BinGroup(dummy_ids, (hbin_t)nsteps,
														2 * bin_schemes_sum, nbins_sum_sq));
				std::vector<DimEntropy> dimentropy_v(block_tasks,
													 DimEntropy(dummy_ids, (hbin_t)nsteps, n_schemes,
																nestimators));

#pragma omp parallel for
				for (auto ri = bi; ri < bi + bl_rows; ++ri)
				{
					auto rno = ri - bi;
					const u_int dtyp_id1 = dimTypeKeys[ri];
					const std::vector<u_int> neigh =
						inputs.getNeighbors().getDimTypeNeighs(dim_type,
															   dtyp_id1);
					for (auto const cj : neigh)
					{
						if (cj >= dimtypes_v[bj] && cj <= dimtypes_v[bj + bl_cols - 1])
						{
							u_int idx = id2index[std::make_pair(dtyp_id1, cj)];
							auto cno = read_cols_dims[cj];
							if (binData2D(nsteps, bin_schemes, nframes_tot_eff,
										  dtyps1_extrm_v[rno][0],
										  dtyps1_extrm_v[rno][1],
										  dtyps2_extrm_v[cno][0],
										  dtyps2_extrm_v[cno][1], dtyps1_int_v[rno],
										  dtyps2_int_v[cno], freq2_ad_obs_v[idx],
										  bin2_ad_mids_v[idx]) != 0)
							{
								LOG_ERROR("Error: Binning data %s id(%u, %u)",
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

							std::vector<u_int> dim_ids(2);
							dim_ids[0] = dtyp_id1;
							dim_ids[1] = cj;
							std::vector<double> ad_extrm(8, 0.0);
							for (int xx = 0; xx < 4; ++xx)
							{
								ad_extrm[xx] = dtyps1_extrm_v[rno][xx];
								ad_extrm[xx + 4] = dtyps2_extrm_v[cno][xx];
							}
							if (isWritefreq && (frqWriteset & BATSet::AD2D))
							{
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
				auto _endcomp = Clock::now();
				auto _durcomp = std::chrono::duration_cast<
									ClockResolution>(_endcomp - _endr2)
									.count();
				ptaskset->curr_comp += ClockResolution(_durcomp);
				progressstate.entc.curr_comp += ClockResolution(
					_durcomp);

				auto _endcomm2 = Clock::now();
				auto _durcomm2 = std::chrono::duration_cast<
									 ClockResolution>(_endcomm2 - _endcomp)
									 .count();
				ptaskset->curr_comm += ClockResolution(_durcomm2);
				progressstate.entc.curr_comm += ClockResolution(
					_durcomm2);

				if (isWritefreq && (frqWriteset & BATSet::AD2D))
				{
					if (xy2DHist->writeRecords(bingrp_v) != 0)
					{
						LOG_ERROR("Writing binned data for %s ids rows[%u, %u] ids cols[%u, %u]",
								  str_dim_type.c_str(), bi, bi + bl_rows, bj, bj + bl_cols);
					}
				}
				ent_xy2d->writeRecords(dimentropy_v);

				auto _endwrt = Clock::now();
				auto _durwrt = std::chrono::duration_cast<
								   ClockResolution>(_endwrt - _endcomm2)
								   .count();
				ptaskset->curr_write += ClockResolution(_durwrt);
				progressstate.entc.curr_write += ClockResolution(
					_durwrt);
				ptaskset->current = Clock::now();
				progressstate.entc.current = ptaskset->current;
				ptaskset->done_tasks += block_tasks;
				progressstate.entc.done_tasks += block_tasks;

				entc_tasks_done += block_tasks;
				if (entc_tasks_done >= entc_tasks_freq)
				{
					std::ofstream info_strm(
						inputs.getControl().getOutfilepath() + inputs.getControl().getInfofile());
					info_strm << progressstate.toString();
					info_strm.close();
					entc_tasks_done %= entc_tasks_freq;
				}
			}
		}
	} // ang-dih block calc for d1 ends
	ptaskset->current = Clock::now();
	ptaskset->cstt = ExecState::COMPLETED;
	progressstate.entc.current = ptaskset->current;
	int entc_workset_intval = static_cast<int>(inputs.getEntropy().getWorkSet());
	int entc_bd2d_intval = static_cast<int>(BATSet::BD2D);
	int entc_ad2d_intval = static_cast<int>(BATSet::AD2D);
	int entc_xy2d_intval = static_cast<int>(BATSet::XY2D);
	if ((entc_workset_intval < entc_bd2d_intval) && (dim_type == BATSet::BA2D))
	{
		progressstate.entc.cstt = ExecState::COMPLETED;
	}
	else if ((entc_workset_intval < entc_ad2d_intval) && (dim_type == BATSet::BD2D))
	{
		progressstate.entc.cstt = ExecState::COMPLETED;
	}
	else if ((entc_workset_intval <= entc_xy2d_intval) && (dim_type == BATSet::AD2D))
	{
		progressstate.entc.cstt = ExecState::COMPLETED;
	}
	std::ofstream info_strm(
		inputs.getControl().getOutfilepath() + inputs.getControl().getInfofile());
	info_strm << progressstate.toString();
	info_strm.close();
	time_xy.stop();
	mprintf(
		"TIME: Histogram bin frequency calculation for %s-2D: %.4f seconds.\n",
		str_dim_type.c_str(), time_xy.total());
}

void EntropyCalculator::run(int num_threads)
{
	if (inputs.getControl().isCalcEntropy())
	{

		if ((inputs.getEntropy().getWorkSet() & BATSet::B1D) || (inputs.getEntropy().getWorkSet() & BATSet::BB2D))
		{
			run1d(BAT_t::BOND, num_threads);
		}

		if ((inputs.getEntropy().getWorkSet() & BATSet::A1D) || (inputs.getEntropy().getWorkSet() & BATSet::AA2D))
		{
			run1d(BAT_t::ANGLE, num_threads);
		}

		if ((inputs.getEntropy().getWorkSet() & BATSet::D1D) || (inputs.getEntropy().getWorkSet() & BATSet::DD2D))
		{
			run1d(BAT_t::DIHEDRAL, num_threads);
		}

		if ((inputs.getEntropy().getWorkSet() & BATSet::BB2D))
		{
			run2d_xx(BATSet::BB2D, num_threads);
		}

		if ((inputs.getEntropy().getWorkSet() & BATSet::AA2D))
		{
			run2d_xx(BATSet::AA2D, num_threads);
		}

		if ((inputs.getEntropy().getWorkSet() & BATSet::DD2D))
		{
			run2d_xx(BATSet::DD2D, num_threads);
		}

		if ((inputs.getEntropy().getWorkSet() & BATSet::BA2D))
		{
			run2d_xy(BATSet::BA2D, num_threads);
		}

		if ((inputs.getEntropy().getWorkSet() & BATSet::BD2D))
		{
			run2d_xy(BATSet::BD2D, num_threads);
		}

		if ((inputs.getEntropy().getWorkSet() & BATSet::AD2D))
		{
			run2d_xy(BATSet::AD2D, num_threads);
		}
	}
}

#endif
