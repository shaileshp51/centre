/*
 * EntropyScorer.h
 *
 *  Created on: Sep 5, 2019
 *      Author: shailesh
 */

#ifndef SRC_ENTROPYCALCULATOR_H_
#define SRC_ENTROPYCALCULATOR_H_

#include "EntropyCalculator.h"
#include "configcentre.h"
#include <algorithm>
#include "main_common.h"

class EntropyCalculator {
	Inputs &inputs;
	ProgState &progressstate;
	bool isWritefreq;
	BATSet frqWriteset;

	std::vector<hbin_t> bin_schemes;
	hbin_t n_schemes;
	u_int bin_schemes_sum;
	ull_int nbins_sum_sq;

	ull_int nframes_tot;
	ull_int startframe;
	ull_int strideframe;
	ull_int nframes_tot_eff;
	hbin_t nsteps;
	hbin_t start_hist_step;
	hbin_t stride_hist_step;
	hbin_t eff_hist_step;
	hbin_t nestimators;

	ull_int nfrm4rl2int;
	ull_int step_size;
	u_int n_bonds;
	u_int n_angles;
	u_int n_diheds;
	u_int n_bnd_eff;
	u_int n_ang_eff;
	u_int n_dih_eff;
	bool isKDE;

	ull_int nframes_trjs;

	ull_int entc_tasks_done;

	ull_int entc_tasks_freq;
	u_int cache_entc_dims_per_cpu;
private:
	void run1d(BAT_t bat_type, int num_threads);

	void run1d_mpi(BAT_t bat_type, const int rank, const int numprocs,
			const int n_thread_perproc);

	void run2d_xx(BATSet dim_type, int num_threads);

	void run2d_xx_mpi(BATSet dim_type, const int rank, const int numprocs,
			const int n_thread_perproc);

	void run2d_xy_mpi(BATSet dim_type, const int rank, const int numprocs,
			const int n_thread_perproc);

	void run2d_xy(BATSet dim_type, int num_threads);

	void run2d_mpi(int num_threads);

public:
	EntropyCalculator(Inputs &inps, ProgState &prg_state) :
			inputs(inps), progressstate(prg_state) {

		isWritefreq = inputs.getHist().isWritefreq();
		frqWriteset = inputs.getHist().getWriteSet();

		if (inputs.getBats().getPdfmethod() == PDFMethod::vonMisesKDE) {
			isWritefreq = inputs.getVmkde().isWritefreq();
			frqWriteset = inputs.getVmkde().getWriteSet();
		}

		bin_schemes = inputs.getHist().getBinSchemes();
		n_schemes = (hbin_t) bin_schemes.size();
		bin_schemes_sum = 0;
		nbins_sum_sq = 0;

		for (size_t i = 0; i < bin_schemes.size(); ++i) {
			bin_schemes_sum += bin_schemes[i];
			nbins_sum_sq += bin_schemes[i] * bin_schemes[i];
		}

		nframes_tot = inputs.getEntropy().getNumframe();
		startframe = inputs.getEntropy().getStartframe();
		strideframe = inputs.getEntropy().getStrideframe();
		nframes_tot_eff = inputs.getEntropy().getNframesEff();
		nsteps = inputs.getControl().getNsteps();
		start_hist_step = inputs.getHist().getWritestepstart();
		stride_hist_step = inputs.getHist().getWritestepstride();
		eff_hist_step = (inputs.getControl().getNsteps() - start_hist_step
				+ stride_hist_step - 1) / stride_hist_step;
		nestimators = inputs.getEstimators().size();
		ull_int tmp_frm = 0;
		for (auto v : inputs.getBats().getNframes()) {
			tmp_frm += v;
		}
		nfrm4rl2int = tmp_frm; // inp->getEntropy().getNumframe();
		step_size = nframes_tot_eff / nsteps;
		n_bonds = inputs.getBats().getNbond();
		n_angles = inputs.getBats().getNangle();
		n_diheds = inputs.getBats().getNdihed();
		n_bnd_eff = inputs.getSubset().bondsSize();
		n_ang_eff = inputs.getSubset().anglesSize();
		n_dih_eff = inputs.getSubset().torsionsSize();
		isKDE = inputs.getBats().getPdfmethod() == PDFMethod::vonMisesKDE; // false if histogram
		nframes_trjs = 0;

		entc_tasks_done = 0;
		entc_tasks_freq = inputs.getControl().getEntcinfofreq();
		cache_entc_dims_per_cpu = inputs.getControl().getCacheentcdimspercpu();
	}
	void run(int num_threads);

	void run_mpi(const int rank, const int numprocs,
			const int n_thread_perproc);
};

#endif /* SRC_ENTROPYCALCULATOR_H_ */
