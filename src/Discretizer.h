/*
 * Discretizer.h
 *
 *  Created on: Sep 5, 2019
 *      Author: shailesh
 */

#ifndef SRC_DISCRETIZER_H_
#define SRC_DISCRETIZER_H_

#include "configcentre.h"
#include <algorithm>
#include <iostream>
#include "main_common.h"

class Discretizer {
	Inputs &inputs;
	ProgState &progressstate;
	ull_int nframes_tot_eff;
	ull_int nfrm4rl2int;
	ull_int dscrt_tasks_done;
	ull_int dscrt_tasks_freq;
	u_int cache_dscrt_dims_per_cpu;
	ull_int nframes_trjs = 0;
	std::vector<hbin_t> bin_schemes;
	hbin_t n_schemes;
	std::string fname;
private:
	void run_dtype(BAT_t dtype, std::vector<ull_int> &shuffle_index,
			const int num_threads);
#ifdef USE_OMPMPI
	void run_dtype_mpi(BAT_t dtype, std::vector<ull_int> &shuffle_index,
			const int rank, const int numprocs, const int n_thread_perproc);
#endif
public:
	Discretizer(Inputs &inps, ProgState &prg_state) :
			inputs(inps), progressstate(prg_state) {
		nfrm4rl2int = 0;
		dscrt_tasks_done = 0;
		nframes_tot_eff = inputs.getEntropy().getNframesEff();

		for (auto v : inputs.getBats().getNframes()) {
			nfrm4rl2int += v;
		}
		dscrt_tasks_freq = inputs.getControl().getDscrinfofreq();
		cache_dscrt_dims_per_cpu =
				inputs.getControl().getCachedscrtdimspercpu();

		//ull_int nframes_trjs = 0;
		bin_schemes = inputs.getHist().getBinSchemes();
		n_schemes = bin_schemes.size();
		fname.assign(
				inputs.getControl().getInfilepath()
						+ inputs.getBats().getFnames()[0]);
	}
	void run(std::vector<ull_int> &shuffle_index, const int num_threads);
#ifdef USE_OMPMPI
	void run_mpi(std::vector<ull_int> &shuffle_index, const int rank,
			const int numprocs, const int n_thread_perproc);
#endif
};

#endif /* SRC_DISCRETIZER_H_ */
