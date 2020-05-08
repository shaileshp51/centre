#include "main_common.h"
#include "EntropyScorer.h"
#include "macrologger/macrologger.h"

using namespace std;

void EntropyScorer::run() {
	auto bin_schemes = inputs.getHist().getBinSchemes();
	auto nestimators = inputs.getEstimators().size();

	hbin_t n_schemes = (hbin_t) bin_schemes.size();
	u_int bin_schemes_sum = 0;
	ull_int nbins_sum_sq = 0;
	for (size_t i = 0; i < bin_schemes.size(); ++i) {
		bin_schemes_sum += bin_schemes[i];
		nbins_sum_sq += bin_schemes[i] * bin_schemes[i];
	}

	const hbin_t nsteps = inputs.getControl().getNsteps();
	const ull_int step_size = inputs.getEntropy().getNframesEff() / nsteps;
	const u_int n_bnd_eff = inputs.getSubset().bondsSize();
	const u_int n_ang_eff = inputs.getSubset().anglesSize();
	const u_int n_dih_eff = inputs.getSubset().torsionsSize();

	u_int n_bnd2d = 0;
	u_int n_ang2d = 0;
	u_int n_dih2d = 0;
	u_int n_ba2d = 0;
	u_int n_bd2d = 0;
	u_int n_ad2d = 0;
	if (inputs.getEntropy().getWorkSet() & BATSet::BB2D) {
		vector<u_int> bndKeys;
		inputs.getNeighbors().bondKeys(bndKeys);
		const u_int num_bndKeys = bndKeys.size();
		for (size_t d1 = 0; d1 < num_bndKeys; ++d1) {
			n_bnd2d += inputs.getNeighbors().bondNeighSize(bndKeys[d1]);
		}
	}
	if (inputs.getEntropy().getWorkSet() & BATSet::AA2D) {
		vector<u_int> angKeys;
		inputs.getNeighbors().angleKeys(angKeys);
		const u_int num_angKeys = angKeys.size();
		for (size_t d1 = 0; d1 < num_angKeys; ++d1) {
			n_ang2d += inputs.getNeighbors().angleNeighSize(angKeys[d1]);
		}
	}
	if (inputs.getEntropy().getWorkSet() & BATSet::DD2D) {
		vector<u_int> dihKeys;
		inputs.getNeighbors().torsionKeys(dihKeys);
		const u_int num_dihKeys = dihKeys.size();
		for (size_t d1 = 0; d1 < num_dihKeys; ++d1) {
			n_dih2d += inputs.getNeighbors().torsionNeighSize(dihKeys[d1]);
		}
	}
	if (inputs.getEntropy().getWorkSet() & BATSet::BA2D) {
		vector<u_int> baKeys;
		inputs.getNeighbors().bacrossKeys(baKeys);
		const u_int num_baKeys = baKeys.size();
		for (size_t d1 = 0; d1 < num_baKeys; ++d1) {
			n_ba2d += inputs.getNeighbors().baNeighSize(baKeys[d1]);
		}
	}
	if (inputs.getEntropy().getWorkSet() & BATSet::BD2D) {
		vector<u_int> bdKeys;
		inputs.getNeighbors().bdcrossKeys(bdKeys);
		const u_int num_bdKeys = bdKeys.size();
		for (size_t d1 = 0; d1 < num_bdKeys; ++d1) {
			n_bd2d += inputs.getNeighbors().bdNeighSize(bdKeys[d1]);
		}
	}
	if (inputs.getEntropy().getWorkSet() & BATSet::AD2D) {
		vector<u_int> adKeys;
		inputs.getNeighbors().adcrossKeys(adKeys);
		const u_int num_adKeys = adKeys.size();
		for (size_t d1 = 0; d1 < num_adKeys; ++d1) {
			n_ad2d += inputs.getNeighbors().adNeighSize(adKeys[d1]);
		}
	}
	if (inputs.getControl().isGenConvgdata()) {
		//
		// Convergence build code goes here
		//
		Netcdf_EntContri *ent_bnd1d;
		Netcdf_EntContri *ent_ang1d;
		Netcdf_EntContri *ent_dih1d;
		Netcdf_EntContri *ent_bnd2d;
		Netcdf_EntContri *ent_ang2d;
		Netcdf_EntContri *ent_dih2d;
		Netcdf_EntContri *ent_ba2d;
		Netcdf_EntContri *ent_bd2d;
		Netcdf_EntContri *ent_ad2d;

		Timer time_conv;

		time_conv.start();
		vector<ull_int> dim_len(1, n_bnd_eff);
		ent_bnd1d = new Netcdf_EntContri(
				inputs.getControl().getOutfilepath() + "entcontri_bnd-1d.nc",
				"B", 1, 1, nsteps, dim_len, TensorType::FULL, bin_schemes,
				nestimators);

		dim_len[0] = n_ang_eff;
		ent_ang1d = new Netcdf_EntContri(
				inputs.getControl().getOutfilepath() + "entcontri_ang-1d.nc",
				"A", 1, 1, nsteps, dim_len, TensorType::FULL, bin_schemes,
				nestimators);

		dim_len[0] = n_dih_eff;
		ent_dih1d = new Netcdf_EntContri(
				inputs.getControl().getOutfilepath() + "entcontri_tor-1d.nc",
				"T", 1, 1, nsteps, dim_len, TensorType::FULL, bin_schemes,
				nestimators);

		vector<ull_int> dim_lens(1, n_bnd_eff);
		dim_lens[0] = n_bnd_eff;
		dim_lens[1] = n_bnd_eff;
		ent_bnd2d = new Netcdf_EntContri(
				inputs.getControl().getOutfilepath() + "entcontri_bnd-2d.nc",
				"B/B", 1, 2, nsteps, dim_lens, TensorType::FULL, bin_schemes,
				nestimators);

		dim_lens[0] = n_ang_eff;
		dim_lens[1] = n_ang_eff;
		ent_ang2d = new Netcdf_EntContri(
				inputs.getControl().getOutfilepath() + "entcontri_ang-2d.nc",
				"A/A", 1, 2, nsteps, dim_lens, TensorType::FULL, bin_schemes,
				nestimators);

		dim_lens[0] = n_dih_eff;
		dim_lens[1] = n_dih_eff;
		ent_dih2d = new Netcdf_EntContri(
				inputs.getControl().getOutfilepath() + "entcontri_tor-2d.nc",
				"T/T", 1, 2, nsteps, dim_lens, TensorType::FULL, bin_schemes,
				nestimators);

		dim_lens[0] = n_bnd_eff;
		dim_lens[1] = n_ang_eff;
		ent_ba2d = new Netcdf_EntContri(
				inputs.getControl().getOutfilepath() + "entcontri_ba-2d.nc",
				"B/A", 1, 2, nsteps, dim_lens, TensorType::FULL, bin_schemes,
				nestimators);

		dim_lens[0] = n_bnd_eff;
		dim_lens[1] = n_dih_eff;
		ent_bd2d = new Netcdf_EntContri(
				inputs.getControl().getOutfilepath() + "entcontri_bd-2d.nc",
				"B/T", 1, 2, nsteps, dim_lens, TensorType::FULL, bin_schemes,
				nestimators);

		dim_lens[0] = n_ang_eff;
		dim_lens[1] = n_dih_eff;
		ent_ad2d = new Netcdf_EntContri(
				inputs.getControl().getOutfilepath() + "entcontri_ad-2d.nc",
				"A/D", 1, 2, nsteps, dim_lens, TensorType::FULL, bin_schemes,
				nestimators);

		map<u_int, u_int> id2index_b;
		map<u_int, u_int> id2index_a;
		map<u_int, u_int> id2index_d;

		map<u_int, pair<u_int, u_int>> index2id_bb;
		map<u_int, pair<u_int, u_int>> index2id_ba;
		map<u_int, pair<u_int, u_int>> index2id_bd;
		map<u_int, pair<u_int, u_int>> index2id_aa;
		map<u_int, pair<u_int, u_int>> index2id_ad;
		map<u_int, pair<u_int, u_int>> index2id_dd;

		if (inputs.getEntropy().getWorkSet() & BATSet::B1D) {
			vector<u_int> ids_bnd(n_bnd_eff, 0);
			ent_bnd1d->readIds(0, n_bnd_eff, ids_bnd);
			for (size_t itmp = 0u; itmp < ids_bnd.size(); ++itmp) {
				id2index_b[ids_bnd[itmp]] = itmp;
			}
		}

		if (inputs.getEntropy().getWorkSet() & BATSet::A1D) {
			vector<u_int> ids_ang(n_ang_eff, 0);
			ent_ang1d->readIds(0, n_ang_eff, ids_ang);
			for (size_t itmp = 0; itmp < ids_ang.size(); ++itmp) {
				id2index_a[ids_ang[itmp]] = itmp;
			}
		}

		if (inputs.getEntropy().getWorkSet() & BATSet::D1D) {
			vector<u_int> ids_dih(n_dih_eff, 0);
			ent_dih1d->readIds(0, n_dih_eff, ids_dih);
			for (size_t itmp = 0; itmp < ids_dih.size(); ++itmp) {
				id2index_d[ids_dih[itmp]] = itmp;
			}
		}

		if (inputs.getEntropy().getWorkSet() & BATSet::BB2D) {
			vector<u_int> ids_bb2(2 * n_bnd2d, 0);
			ent_bnd2d->readIds(0, n_bnd2d, ids_bb2);
			for (size_t itmp = 0; itmp < n_bnd2d; ++itmp) {
				index2id_bb[itmp] = make_pair(ids_bb2[2 * itmp],
						ids_bb2[2 * itmp + 1]);
			}
		}

		if (inputs.getEntropy().getWorkSet() & BATSet::BA2D) {
			vector<u_int> ids_ba2(2 * n_ba2d, 0);

			ent_ba2d->readIds(0, n_ba2d, ids_ba2);
			for (size_t itmp = 0; itmp < n_ba2d; ++itmp) {
				index2id_ba[itmp] = make_pair(ids_ba2[2 * itmp],
						ids_ba2[2 * itmp + 1]);
			}
		}

		if (inputs.getEntropy().getWorkSet() & BATSet::BD2D) {
			vector<u_int> ids_bd2(2 * n_bd2d, 0);

			ent_bd2d->readIds(0, n_bd2d, ids_bd2);
			for (size_t itmp = 0; itmp < n_bd2d; ++itmp) {
				index2id_bd[itmp] = make_pair(ids_bd2[2 * itmp],
						ids_bd2[2 * itmp + 1]);
			}
		}

		if (inputs.getEntropy().getWorkSet() & BATSet::AA2D) {
			vector<u_int> ids_aa2(2 * n_ang2d, 0);

			ent_ang2d->readIds(0, n_ang2d, ids_aa2);
			for (size_t itmp = 0; itmp < n_ang2d; ++itmp) {
				index2id_aa[itmp] = make_pair(ids_aa2[2 * itmp],
						ids_aa2[2 * itmp + 1]);
			}
		}

		if (inputs.getEntropy().getWorkSet() & BATSet::AD2D) {
			vector<u_int> ids_ad2(2 * n_ad2d, 0);

			ent_ad2d->readIds(0, n_ad2d, ids_ad2);
			for (size_t itmp = 0; itmp < n_ad2d; ++itmp) {
				index2id_ad[itmp] = make_pair(ids_ad2[2 * itmp],
						ids_ad2[2 * itmp + 1]);
			}
		}

		if (inputs.getEntropy().getWorkSet() & BATSet::DD2D) {
			vector<u_int> ids_dd2(2 * n_dih2d, 0);
			ent_dih2d->readIds(0, n_dih2d, ids_dd2);
			for (size_t itmp = 0; itmp < n_dih2d; ++itmp) {
				index2id_dd[itmp] = make_pair(ids_dd2[2 * itmp],
						ids_dd2[2 * itmp + 1]);
			}
		}

		auto &estimators = inputs.getEstimators();
		for (hbin_t estid = 0; estid < nestimators; ++estid) {
			EntropyEstimators est = estimators[estid];
			for (hbin_t b_sch = 0; b_sch < n_schemes; ++b_sch) {
				double sum_b = 0.0, sum_a = 0.0, sum_d = 0.0, singlet = 0.0;
				double sum_bb = 0.0, sum_aa = 0.0, sum_dd = 0.0, sum_ba = 0.0,
						sum_bd = 0.0, sum_ad = 0.0;
				double mi_bb = 0.0, mi_aa = 0.0, mi_dd = 0.0, mi_ba = 0.0,
						mi_bd = 0.0, mi_ad = 0.0, doublet = 0.0;
				double mist_bb = 0.0, mist_aa = 0.0, mist_dd = 0.0, mist_ba =
						0.0, mist_bd = 0.0, mist_ad = 0.0, mist_doublet;

				ofstream miefile;
				ofstream mistfile;
				for (hbin_t stp_id = 0; stp_id < nsteps; ++stp_id) {

					sum_b = 0.0;
					sum_a = 0.0;
					sum_d = 0.0;
					singlet = 0.0;
					sum_bb = 0.0;
					sum_aa = 0.0;
					sum_dd = 0.0;
					sum_ba = 0.0;
					sum_bd = 0.0;
					sum_ad = 0.0;
					mist_bb = 0.0;
					mist_aa = 0.0;
					mist_dd = 0.0;
					mist_ba = 0.0;
					mist_bd = 0.0;
					mist_ad = 0.0;
					mist_doublet = 0.0;

					vector<double> entropy_b1dv(n_bnd_eff, 0.0);
					vector<double> entropy_a1dv(n_ang_eff, 0.0);
					vector<double> entropy_d1dv(n_dih_eff, 0.0);

					vector<double> entropy_bbv(n_bnd2d, 0.0);
					vector<double> mi_bbv(n_bnd2d, 0.0);

					vector<double> entropy_aav(n_ang2d, 0.0);
					vector<double> mi_aav(n_ang2d, 0.0);

					vector<double> entropy_ddv(n_dih2d, 0.0);
					vector<double> mi_ddv(n_dih2d, 0.0);

					vector<double> entropy_bav(n_ba2d, 0.0);
					vector<double> mi_bav(n_ba2d, 0.0);

					vector<double> entropy_bdv(n_bd2d, 0.0);
					vector<double> mi_bdv(n_bd2d, 0.0);

					vector<double> entropy_adv(n_ad2d, 0.0);
					vector<double> mi_adv(n_ad2d, 0.0);

					vector<nodeMIST> miSpanTree;

					if ((inputs.getEntropy().getWorkSet() & BATSet::B1D)
							|| (inputs.getEntropy().getWorkSet() & BATSet::BB2D)) {
						if (stp_id == 0 && b_sch == 0) {
							ent_bnd1d->setupRead();
						}
						if (ent_bnd1d->readEntrContrib(0, n_bnd_eff, estid, 1,
								b_sch, 1, stp_id, 1, entropy_b1dv) != 0) {
							LOG_ERROR(
									"Reading entropy contributions for bonds step-id(%d)",
									stp_id);
						}
						sumContrib1D(n_bnd_eff, entropy_b1dv, sum_b);
					}

					if ((inputs.getEntropy().getWorkSet() & BATSet::A1D)
							|| (inputs.getEntropy().getWorkSet() & BATSet::AA2D)) {
						if (stp_id == 0 && b_sch == 0) {
							ent_ang1d->setupRead();
						}
						if (ent_ang1d->readEntrContrib(0, n_ang_eff, estid, 1,
								b_sch, 1, stp_id, 1, entropy_a1dv) != 0) {
							LOG_ERROR(
									"Reading entropy contributions for angles step-id(%d)",
									stp_id);
						}
						sumContrib1D(n_ang_eff, entropy_a1dv, sum_a);

					}

					if ((inputs.getEntropy().getWorkSet() & BATSet::D1D)
							|| (inputs.getEntropy().getWorkSet() & BATSet::DD2D)) {
						if (stp_id == 0 && b_sch == 0) {
							ent_dih1d->setupRead();
						}
						if (ent_dih1d->readEntrContrib(0, n_dih_eff, estid, 1,
								b_sch, 1, stp_id, 1, entropy_d1dv) != 0) {
							LOG_ERROR(
									"Reading entropy contributions for torsions step-id(%d)",
									stp_id);
						}
						sumContrib1D(n_dih_eff, entropy_d1dv, sum_d);
					}

					if (inputs.getEntropy().getWorkSet() & BATSet::BB2D) {

						if (stp_id == 0 && b_sch == 0) {
							ent_bnd2d->setupRead();
						}

						u_int nei_bb2d = inputs.getNeighbors().bondsSize();
						vector<u_int> bndKeys;
						inputs.getNeighbors().bondKeys(bndKeys);

						if (ent_bnd2d->readEntrContrib(0, n_bnd2d, estid, 1,
								b_sch, 1, stp_id, 1, entropy_bbv) != 0) {
							LOG_ERROR(
									"Reading entropy contributions for B/B step-id(%d)",
									stp_id);
						}

						double sum_bb_mst;
						sumContrib2D(0, n_bnd2d, entropy_bbv, sum_bb_mst);

						getMaxMI2D(id2index_b, id2index_b, entropy_b1dv,
								entropy_b1dv, index2id_bb, entropy_bbv, mi_bbv);

						sum_bb = sum_bb_mst;
					}

					if (inputs.getEntropy().getWorkSet() & BATSet::AA2D) {

						if (stp_id == 0 && b_sch == 0) {
							ent_ang2d->setupRead();
						}

						if (ent_ang2d->readEntrContrib(0, n_ang2d, estid, 1,
								b_sch, 1, stp_id, 1, entropy_aav) != 0) {
							LOG_ERROR(
									"Reading entropy contributions for A/A step-id(%d)",
									stp_id);
						}

						double sum_aa_mst;
						sumContrib2D(0, n_ang2d, entropy_aav, sum_aa_mst);

						getMaxMI2D(id2index_a, id2index_a, entropy_a1dv,
								entropy_a1dv, index2id_aa, entropy_aav, mi_aav);

						sum_aa = sum_aa_mst;

						double sum_aa_slv = 0.0;
					}

					if (inputs.getEntropy().getWorkSet() & BATSet::DD2D) {

						if (stp_id == 0 && b_sch == 0) {
							ent_dih2d->setupRead();
						}

						if (ent_dih2d->readEntrContrib(0, n_dih2d, estid, 1,
								b_sch, 1, stp_id, 1, entropy_ddv) != 0) {
							LOG_ERROR(
									"Reading entropy contributions for D/D step-id(%d)",
									stp_id);
						}

						double sum_dd_mst;
						sumContrib2D(0, n_dih2d, entropy_ddv, sum_dd_mst);

						getMaxMI2D(id2index_d, id2index_d, entropy_d1dv,
								entropy_d1dv, index2id_dd, entropy_ddv, mi_ddv);

						sum_dd = sum_dd_mst;
					}

					if (inputs.getEntropy().getWorkSet() & BATSet::BA2D) {

						if (stp_id == 0 && b_sch == 0) {
							ent_ba2d->setupRead();
						}

						if (ent_ba2d->readEntrContrib(0, n_ba2d, estid, 1,
								b_sch, 1, stp_id, 1, entropy_bav) != 0) {
							LOG_ERROR(
									"Reading entropy contributions for B/A step-id(%d)",
									stp_id);
						}

						double sum_ba_mst;
						sumContrib2D(0, n_ba2d, entropy_bav, sum_ba_mst);

						getMaxMI2D(id2index_b, id2index_a, entropy_b1dv,
								entropy_a1dv, index2id_ba, entropy_bav, mi_bav);

						sum_ba = sum_ba_mst;

					}

					if (inputs.getEntropy().getWorkSet() & BATSet::BD2D) {

						if (stp_id == 0 && b_sch == 0) {
							ent_bd2d->setupRead();
						}

						if (ent_bd2d->readEntrContrib(0, n_bd2d, estid, 1,
								b_sch, 1, stp_id, 1, entropy_bdv) != 0) {
							LOG_ERROR(
									"Reading entropy contributions for B/D step-id(%d)",
									stp_id);
						}

						double sum_bd_mst;
						sumContrib2D(0, n_bd2d, entropy_bdv, sum_bd_mst);

						getMaxMI2D(id2index_b, id2index_d, entropy_b1dv,
								entropy_d1dv, index2id_bd, entropy_bdv, mi_bdv);

						sum_bd = sum_bd_mst;
					}

					if (inputs.getEntropy().getWorkSet() & BATSet::AD2D) {

						if (stp_id == 0 && b_sch == 0) {
							ent_ad2d->setupRead();
						}

						if (ent_ad2d->readEntrContrib(0, n_ad2d, estid, 1,
								b_sch, 1, stp_id, 1, entropy_adv) != 0) {
							mprinterr(
									"Reading entropy contributions for A/D step-id(%d)",
									stp_id);
						}

						double sum_ad_mst;
						sumContrib2D(0, n_ad2d, entropy_adv, sum_ad_mst);

						getMaxMI2D(id2index_a, id2index_d, entropy_a1dv,
								entropy_d1dv, index2id_ad, entropy_adv, mi_adv);

						sum_ad = sum_ad_mst;
					}

					singlet = sum_b + sum_a + sum_d;

					mi_bb = accumulate(mi_bbv.begin(), mi_bbv.end(), 0.0);
					mi_aa = accumulate(mi_aav.begin(), mi_aav.end(), 0.0);
					mi_dd = accumulate(mi_ddv.begin(), mi_ddv.end(), 0.0);

					mi_ba = accumulate(mi_bav.begin(), mi_bav.end(), 0.0);
					mi_bd = accumulate(mi_bdv.begin(), mi_bdv.end(), 0.0);
					mi_ad = accumulate(mi_adv.begin(), mi_adv.end(), 0.0);

					doublet = singlet - mi_bb - mi_aa - mi_dd - mi_ba - mi_bd
							- mi_ad;
					int fw = 13;
					if (stp_id == 0) {

						string fname = inputs.getControl().getOutfilepath()
								+ "MIE" + "-" + getEstimatorName(est)
								+ "_convergence_b-"
								+ to_string((int) bin_schemes[b_sch]) + "."
								+ inputs.getEntropy().getConvergenceFileExt();
						miefile.open(fname.c_str());

						miefile << setw(fw) << fixed << left << "TIME" << ", "
								<< setw(fw) << left << "BOND-S1" << ", "
								<< setw(fw) << fixed << left << "ANGLE-S1"
								<< ", " << setw(fw) << fixed << left
								<< "DIHED-S1" << ", " << setw(fw) << fixed
								<< left << "B/B-S2" << ", " << setw(fw) << left
								<< "B/B-I2" << ", " << setw(fw) << fixed << left
								<< "A/A-S2" << ", " << setw(fw) << fixed << left
								<< "A/A-I2" << ", " << setw(fw) << fixed << left
								<< "D/D-S2" << ", " << setw(fw) << left
								<< "D/D-I2" << ", " << setw(fw) << fixed << left
								<< "B/A-S2" << ", " << setw(fw) << fixed << left
								<< "B/A-I2" << ", " << setw(fw) << fixed << left
								<< "B/D-S2" << ", " << setw(fw) << left
								<< "B/D-I2" << ", " << setw(fw) << fixed << left
								<< "A/D-S2" << ", " << setw(fw) << fixed << left
								<< "A/D-I2" << ", " << setw(fw) << fixed << left
								<< "S1" << ", " << setw(fw) << left << "S2"
								<< endl;
					}
					ull_int time = step_size * (stp_id + 1);
					miefile.precision(4);
					miefile << setw(fw) << fixed << left << time << ", "
							<< setw(fw) << fixed << left << sum_b << ", "
							<< setw(fw) << fixed << left << sum_a << ", "
							<< setw(fw) << fixed << left << sum_d << ", "
							<< setw(fw) << fixed << left << sum_bb << ", "
							<< setw(fw) << fixed << left << mi_bb << ", "
							<< setw(fw) << fixed << left << sum_aa << ", "
							<< setw(fw) << left << mi_aa << ", " << setw(fw)
							<< fixed << left << sum_dd << ", " << setw(fw)
							<< fixed << left << mi_dd << ", " << setw(fw)
							<< left << sum_ba << ", " << setw(fw) << fixed
							<< left << mi_ba << ", " << setw(fw) << fixed
							<< left << sum_bd << ", " << setw(fw) << left
							<< mi_bd << ", " << setw(fw) << fixed << left
							<< sum_ad << ", " << setw(fw) << fixed << left
							<< mi_ad << ", " << setw(fw) << left << singlet
							<< ", " << setw(fw) << fixed << left << doublet
							<< endl;
					time += step_size;

					if (inputs.getMethod() == EntropyMethods::MIST) {
						Graph graph;
						u_int ang_node_id_offset = 0;
						u_int dih_node_id_offset = 0;
						// add nodes to graph, here bond nodes are added first and their indices start
						// from 0
						if ((inputs.getEntropy().getWorkSet() & BATSet::BB2D)
								|| (inputs.getEntropy().getWorkSet()
										& BATSet::BA2D)
								|| (inputs.getEntropy().getWorkSet()
										& BATSet::BD2D)) {
							for (u_int i1 = 0; i1 < n_bnd_eff; ++i1) {
								graph.nodes.push_back(
										Node(BAT_t::BOND,
												inputs.getSubset().getBonds()[i1]));
							}
							ang_node_id_offset = n_bnd_eff;
							dih_node_id_offset = n_bnd_eff;
						}
						// add nodes to graph, here angle nodes are added after bonds and their indices start
						// from n_bnd_eff
						if ((inputs.getEntropy().getWorkSet() & BATSet::AA2D)
								|| (inputs.getEntropy().getWorkSet()
										& BATSet::AD2D)) {
							for (u_int i1 = 0; i1 < n_ang_eff; ++i1) {
								graph.nodes.push_back(
										Node(BAT_t::ANGLE,
												inputs.getSubset().getAngles()[i1]));
							}
							dih_node_id_offset += n_ang_eff;
						}
						// add nodes to graph, here torsion nodes are added after angles and their indices start
						// from n_bnd_eff + n_ang_eff
						if ((inputs.getEntropy().getWorkSet() & BATSet::DD2D)) {
							for (u_int i1 = 0; i1 < n_dih_eff; ++i1) {
								graph.nodes.push_back(
										Node(BAT_t::DIHEDRAL,
												inputs.getSubset().getTorsions()[i1]));
							}
						}

						// Now mi terms are being converted to graph edges and being added to the graph
						if ((inputs.getEntropy().getWorkSet() & BATSet::BB2D)) {
							for (auto const &pp : index2id_bb) {
								auto pv = pp.second;
								Edge edg(id2index_b[pv.first],
										id2index_b[pv.second],
										mi_bbv[pp.first]);
								graph.edges.push_back(edg);
							}
						}

						if ((inputs.getEntropy().getWorkSet() & BATSet::BA2D)) {
							for (auto const &pp : index2id_ba) {
								auto pv = pp.second;
								Edge edg(id2index_b[pv.first],
										id2index_a[pv.second]
												+ ang_node_id_offset,
										mi_bav[pp.first]);
								graph.edges.push_back(edg);
							}
						}

						if ((inputs.getEntropy().getWorkSet() & BATSet::BD2D)) {
							for (auto const &pp : index2id_bd) {
								auto pv = pp.second;
								Edge edg(id2index_b[pv.first],
										id2index_d[pv.second]
												+ dih_node_id_offset,
										mi_bdv[pp.first]);
								graph.edges.push_back(edg);
							}
						}
						if ((inputs.getEntropy().getWorkSet() & BATSet::AA2D)) {
							for (auto const &pp : index2id_aa) {
								auto pv = pp.second;
								Edge edg(
										id2index_a[pv.first]
												+ ang_node_id_offset,
										id2index_a[pv.second]
												+ ang_node_id_offset,
										mi_aav[pp.first]);
								graph.edges.push_back(edg);
							}
						}

						if ((inputs.getEntropy().getWorkSet() & BATSet::AD2D)) {
							for (auto const &pp : index2id_ad) {
								auto pv = pp.second;
								Edge edg(
										id2index_a[pv.first]
												+ ang_node_id_offset,
										id2index_d[pv.second]
												+ dih_node_id_offset,
										mi_adv[pp.first]);
								graph.edges.push_back(edg);
							}
						}

						if ((inputs.getEntropy().getWorkSet() & BATSet::DD2D)) {
							for (auto const &pp : index2id_dd) {
								auto pv = pp.second;
								Edge edg(
										id2index_d[pv.first]
												+ dih_node_id_offset,
										id2index_d[pv.second]
												+ dih_node_id_offset,
										mi_ddv[pp.first]);
								graph.edges.push_back(edg);
							}
						}

						bool isTreeBuild = false;
						vector<Edge> mistree = boruvkaMST(graph);

						ofstream treefile;
						string fmisttree = inputs.getControl().getOutfilepath()
								+ "Tree-MIST" + "-" + getEstimatorName(est)
								+ "_b-" + to_string((int) bin_schemes[b_sch])
								+ "_s-" + to_string((int) stp_id) + "."
								+ inputs.getEntropy().getConvergenceFileExt();
						treefile.open(fmisttree.c_str());
						treefile << "TD1, TD2, IndexDim1, IndexDim2, MI"
								<< endl;
						char type_c[5] = "-BAT";

						sort(mistree.begin(), mistree.end(), cmpEdgeMIST);
						for (size_t tmp = 0; tmp < mistree.size(); ++tmp) {
							treefile << setw(3) << left
									<< type_c[(int) graph.nodes[mistree[tmp].source].type]
									<< ", " << setw(3) << left
									<< type_c[(int) graph.nodes[mistree[tmp].dest].type]
									<< ", " << setw(9) << left
									<< graph.nodes[mistree[tmp].source].id
									<< ", " << setw(9) << left
									<< graph.nodes[mistree[tmp].dest].id << ", "
									<< setw(3) << left << setprecision(6)
									<< fixed << mistree[tmp].weight << endl;
						}

						treefile.close();
						mist_doublet = 0.0, mist_bb = 0.0, mist_ba = 0.0, mist_bd =
								0.0;
						mist_aa = 0.0, mist_ad = 0.0, mist_dd = 0.0;
						for (u_int t = 0; t < mistree.size(); t++) {
							Edge e = mistree[t];
							if (graph.nodes[e.source].type == BAT_t::BOND
									&& graph.nodes[e.dest].type
											== BAT_t::BOND) {
								mist_bb += e.weight;
							} else if ((graph.nodes[e.source].type
									== BAT_t::BOND
									&& graph.nodes[e.dest].type == BAT_t::ANGLE)
									|| (graph.nodes[e.source].type
											== BAT_t::ANGLE
											&& graph.nodes[e.dest].type
													== BAT_t::BOND)) {
								mist_ba += e.weight;
							} else if ((graph.nodes[e.source].type
									== BAT_t::BOND
									&& graph.nodes[e.dest].type
											== BAT_t::DIHEDRAL)
									|| (graph.nodes[e.source].type
											== BAT_t::DIHEDRAL
											&& graph.nodes[e.dest].type
													== BAT_t::BOND)) {
								mist_bd += e.weight;
							} else if (graph.nodes[e.source].type
									== BAT_t::ANGLE
									&& graph.nodes[e.dest].type
											== BAT_t::ANGLE) {
								mist_aa += e.weight;
							} else if ((graph.nodes[e.source].type
									== BAT_t::ANGLE
									&& graph.nodes[e.dest].type
											== BAT_t::DIHEDRAL)
									|| (graph.nodes[e.source].type
											== BAT_t::DIHEDRAL
											&& graph.nodes[e.dest].type
													== BAT_t::ANGLE)) {
								mist_ad += e.weight;
							} else if (graph.nodes[e.source].type
									== BAT_t::DIHEDRAL
									&& graph.nodes[e.dest].type
											== BAT_t::DIHEDRAL) {
								mist_dd += e.weight;
							}
						}
						mist_doublet = singlet - mist_bb - mist_aa - mist_dd
								- mist_ba - mist_bd - mist_ad;
						int fw = 13;
						if (stp_id == 0) {
							string fname2 =
									inputs.getControl().getOutfilepath()
											+ "MIST" + "-"
											+ getEstimatorName(est)
											+ "_convergence_b-"
											+ to_string(
													(int) bin_schemes[b_sch])
											+ "."
											+ inputs.getEntropy().getConvergenceFileExt();
							mistfile.open(fname2.c_str());

							mistfile << setw(fw) << fixed << left << "TIME"
									<< ", " << setw(fw) << left << "BOND-S1"
									<< ", " << setw(fw) << fixed << left
									<< "ANGLE-S1" << ", " << setw(fw) << fixed
									<< left << "DIHED-S1" << ", " << setw(fw)
									<< fixed << left << "B/B-S2" << ", "
									<< setw(fw) << fixed << left << "B/B-I2"
									<< ", " << setw(fw) << left << "A/A-S2"
									<< ", " << setw(fw) << fixed << left
									<< "A/A-I2" << ", " << setw(fw) << fixed
									<< left << "D/D-S2" << ", " << setw(fw)
									<< fixed << left << "D/D-I2" << ", "
									<< setw(fw) << fixed << left << "B/A-S2"
									<< ", " << setw(fw) << left << "B/A-I2"
									<< ", " << setw(fw) << fixed << left
									<< "B/D-S2" << ", " << setw(fw) << fixed
									<< left << "B/D-I2" << ", " << setw(fw)
									<< fixed << left << "A/D-S2" << ", "
									<< setw(fw) << left << "A/D-I2" << ", "
									<< setw(fw) << fixed << left << "S1" << ", "
									<< setw(fw) << fixed << left << "S2"
									<< endl;
						}
						ull_int time = step_size * (stp_id + 1);
						mistfile.precision(4);
						mistfile << setw(fw) << fixed << left << time << ", "
								<< setw(fw) << fixed << left << sum_b << ", "
								<< setw(fw) << fixed << left << sum_a << ", "
								<< setw(fw) << left << sum_d << ", " << setw(fw)
								<< fixed << left << sum_bb << ", " << setw(fw)
								<< fixed << left << mist_bb << ", " << setw(fw)
								<< left << sum_aa << ", " << setw(fw) << fixed
								<< left << mist_aa << ", " << setw(fw) << fixed
								<< left << sum_dd << ", " << setw(fw) << left
								<< mist_dd << ", " << setw(fw) << fixed << left
								<< sum_ba << ", " << setw(fw) << fixed << left
								<< mist_ba << ", " << setw(fw) << left << sum_bd
								<< ", " << setw(fw) << fixed << left << mist_bd
								<< ", " << setw(fw) << fixed << left << sum_ad
								<< ", " << setw(fw) << left << mist_ad << ", "
								<< setw(fw) << fixed << left << singlet << ", "
								<< setw(fw) << fixed << left << mist_doublet
								<< endl;
						time += step_size;
						if (stp_id == nsteps - 1) {
							mistfile.close();
						}
					}

				}
			}
		}

		time_conv.stop();
		mprintf("TIME: MIE convergence data building time: %.4f seconds.\n",
				time_conv.total());

	}
}
