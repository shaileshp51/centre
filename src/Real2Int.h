/*
 * Real2Int.h
 *
 *  Created on: 28-Sep-2016
 *      Author: shailesh
 */

#ifndef REAL2INT_H_
#define REAL2INT_H_

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <vector>
#include <cmath>
#include <cstddef>

#include "configcentre.h"
#include "Constants.h"
#include "Utils.h"

class Real2Int {
public:
	static int shuffleArray(const std::vector<ull_int> &rand_indices,
			std::vector<hbin_t> &data_int);

	static int bndReal2Int(const u_int id,
			const std::vector<hbin_t> &bin_schemes,
			const std::vector<data_t> &data, double *min_d, double *max_d,
			double *avg_d, std::vector<hbin_t> &data_int);

	static int bndReal2Int(const u_int id,
			const std::vector<hbin_t> &bin_schemes,
			const std::vector<ull_int> &rand_indices,
			const std::vector<data_t> &data, double *min_d, double *max_d,
			double *avg_d, std::vector<hbin_t> &data_int);

	static int angReal2Int(const u_int id,
			const std::vector<hbin_t> &bin_schemes,
			const std::vector<data_t> &data, double *min_d, double *max_d,
			double *avg_d, std::vector<hbin_t> &data_int);

	static int angReal2Int(const u_int id,
			const std::vector<hbin_t> &bin_schemes,
			const std::vector<ull_int> &rand_indices,
			const std::vector<data_t> &data, double *min_d, double *max_d,
			double *avg_d, std::vector<hbin_t> &data_int);

	static int dihReal2Int(const u_int id,
			const std::vector<hbin_t> &bin_schemes, const hbin_t bins_ref,
			const bool optimize, const std::vector<data_t> &data, double *min_d,
			double *max_d, double *avg_d, double *phase_d,
			std::vector<hbin_t> &data_int);

	static int getMinimaPositions(const u_int id, bool analyticGrad,
			hbin_t n_conf_max, double k_value, u_int step_sdo,
			double conv_criterion, u_int max_iterations,
			std::vector<data_t> &data, std::vector<double> &minimum,
			hbin_t &found_mins, std::vector<hbin_t> &data_int);
private:
	static void modfitReal2Int(const u_int id, const hbin_t nddens,
			const std::vector<data_t> &data, std::vector<hbin_t> *data_int);

	static double calculateVonMisesGradient(double SmoothingParam, double X,
			std::vector<data_t> &ang);

	static double vonMises(double SmoothingParam, double X,
			std::vector<data_t> &ang, double &vMdensity, double &vMgradient,
			double &vMsec_deriv);

};

#endif /* REAL2INT_H_ */
