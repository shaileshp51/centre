/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <vector>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <limits>
#include <omp.h>

#include "Real2Int.h"
#include "configcentre.h"

int Real2Int::shuffleArray(const std::vector<ull_int> &rand_indices,
		std::vector<hbin_t> &data_int) {
	if (data_int.size() != rand_indices.size()) {
		return -1;
	} else {
		hbin_t tmp;
		for (size_t i = 0; i < data_int.size(); ++i) {
			tmp = data_int[i];
			data_int[i] = data_int[rand_indices[i]];
			data_int[rand_indices[i]] = tmp;
		}
	}
	return 0;
}

int Real2Int::bndReal2Int(const u_int id,
		const std::vector<hbin_t> &bin_schemes, const std::vector<data_t> &data,
		double *min_d, double *max_d, double *avg_d,
		std::vector<hbin_t> &data_int) {
	/**
	 * Converts the double precision bond lengths to an integer array
	 * indicating bin which containes bond snapshot.
	 */

	if (data.size() * bin_schemes.size() != data_int.size()) {
		return -1;
	}
	const double dx = DELTA_X_RANGE_EXTEND;
	double min_data = 1.0e100, max_data = -1.0e100, avg_data = 0.0, sum_data =
			0.0;
	double bin_width = 0;
	/**
	 * Start by computing the sum of each bond degree of freedom for the
	 * average, For each bond, compute the bin width based on the number of bins and
	 * the bond extrema .
	 */

	for (size_t i = 0; i < data.size(); ++i) {
		sum_data += data[i];
		if (min_data > data[i]) {
			min_data = data[i];
		}
		if (max_data < data[i]) {
			max_data = data[i];
		}
	}
	min_data -= dx;
	if (min_data < 0.0) {
		min_data = 0.0;
	}
	max_data += dx;

	if ((max_data - min_data) < 1.0e-4) {
		max_data += 0.05;
	}
	avg_data = sum_data / data.size();

	for (size_t bidx = 0; bidx < bin_schemes.size(); ++bidx) {
		bin_width = (max_data - min_data) / (int) bin_schemes[bidx];
		ull_int offset = bidx * data.size();
		hbin_t ti;
		for (size_t i = 0; i < data.size(); ++i) {
			ti = (hbin_t) floor((data[i] - min_data) / bin_width);
			data_int[offset + i] = (ti == bin_schemes[bidx]) ? ti - 1 : ti;
		}
	}
	*min_d = min_data;
	*max_d = max_data;
	*avg_d = avg_data;
	return 0;
}

int Real2Int::bndReal2Int(const u_int id,
		const std::vector<hbin_t> &bin_schemes,
		const std::vector<ull_int> &rand_indices,
		const std::vector<data_t> &data, double *min_d, double *max_d,
		double *avg_d, std::vector<hbin_t> &data_int) {
	/**
	 * Converts the double precision bond lengths to an integer array
	 * indicating correct bin location.
	 */

	if (data.size() * bin_schemes.size() != data_int.size()
			|| data_int.size() != rand_indices.size()) {
		return -1;
	}

	const double dx = DELTA_X_RANGE_EXTEND;
	double min_data = 1.0e100, max_data = -1.0e100, avg_data = 0.0, sum_data =
			0.0;
	double bin_width = 0;
	/**
	 * Start by computing the sum of each bond degree of freedom for the
	 * average For each bond, compute the bin width based on the number of bins and
	 * the bond extrema .
	 */

	for (size_t i = 0; i < data.size(); ++i) {
		sum_data += data[i];
		if (min_data > data[i]) {
			min_data = data[i];
		}
		if (max_data < data[i]) {
			max_data = data[i];
		}
	}

	min_data -= dx;
	if (min_data < 0.0) {
		min_data = 0.0;
	}
	max_data += dx;

	if ((max_data - min_data) < 1.0e-4) {
		max_data += 0.05;
	}

	avg_data = sum_data / data.size();

	for (size_t bidx = 0; bidx < bin_schemes.size(); ++bidx) {
		bin_width = (max_data - min_data) / (int) bin_schemes[bidx];
		ull_int offset = bidx * data.size();
		hbin_t ti;
		for (size_t i = 0; i < data.size(); ++i) {
			ti = (hbin_t) floor((data[i] - min_data) / bin_width);
			data_int[offset + i] = (ti == bin_schemes[bidx]) ? ti - 1 : ti;
		}
	}
	*min_d = min_data;
	*max_d = max_data;
	*avg_d = avg_data;

	shuffleArray(rand_indices, data_int);

	return 0;
}

int Real2Int::angReal2Int(const u_int id,
		const std::vector<hbin_t> &bin_schemes, const std::vector<data_t> &data,
		double *min_d, double *max_d, double *avg_d,
		std::vector<hbin_t> &data_int) {
	/**
	 * Converts the double precision bond lengths to an integer array
	 * indicating bin which containes bond snapshot.
	 */

	if (data.size() * bin_schemes.size() != data_int.size()) {
		return -1;
	}

	const double _PI = Constants::PI;
	const double dx = DELTA_X_RANGE_EXTEND;

	double min_data = 1.0e100, max_data = -1.0e100, avg_data = 0.0, sum_data =
			0.0;
	double bin_width = 0;
	/**
	 * Start by computing the sum of each bond degree of freedom for the
	 * average, For each bond, compute the bin width based on the number of bins and
	 * the bond extrema .
	 */

	for (size_t i = 0; i < data.size(); ++i) {
		sum_data += data[i];
		if (min_data > data[i]) {
			min_data = data[i];
		}
		if (max_data < data[i]) {
			max_data = data[i];
		}
	}

	min_data -= dx;
	if (min_data < 0.0) {
		min_data = 0.0;
	}
	max_data += dx;
	if (max_data > _PI) {
		max_data = _PI;
	}

	if ((max_data - min_data) < 1.0e-4) {
		max_data += 0.05;
	}
	avg_data = sum_data / data.size();

	for (size_t bidx = 0; bidx < bin_schemes.size(); ++bidx) {
		bin_width = (max_data - min_data) / (int) bin_schemes[bidx];
		ull_int offset = bidx * data.size();
		hbin_t ti;
		for (size_t i = 0; i < data.size(); ++i) {
			ti = (hbin_t) floor((data[i] - min_data) / bin_width);
			data_int[offset + i] = (ti == bin_schemes[bidx]) ? ti - 1 : ti;
		}
	}

	*min_d = min_data;
	*max_d = max_data;
	*avg_d = avg_data;
	return 0;
}

int Real2Int::angReal2Int(const u_int id,
		const std::vector<hbin_t> &bin_schemes,
		const std::vector<ull_int> &rand_indices,
		const std::vector<data_t> &data, double *min_d, double *max_d,
		double *avg_d, std::vector<hbin_t> &data_int) {

	if (data.size() * bin_schemes.size() != data_int.size()
			|| data_int.size() != rand_indices.size()) {
		return -1;
	}
	double min_data = 1.0e100, max_data = -1.0e100, avg_data = 0.0, sum_data =
			0.0;
	double bin_width = 0;
	const double _PI = Constants::PI;
	const double dx = DELTA_X_RANGE_EXTEND;
	/**
	 * Start by computing the sum of each bond degree of freedom for the
	 * average For each bond, compute the bin width based on the number of bins and
	 * the bond extrema .
	 */

	for (size_t i = 0; i < data.size(); ++i) {
		sum_data += data[i];
		if (min_data > data[i]) {
			min_data = data[i];
		}
		if (max_data < data[i]) {
			max_data = data[i];
		}
	}

	min_data -= dx;
	if (min_data < 0.0) {
		min_data = 0.0;
	}
	max_data += dx;
	if (max_data > _PI) {
		max_data = _PI;
	}
	if ((max_data - min_data) < 1.0e-4) {
		max_data += 0.05;
	}
	avg_data = sum_data / data.size();

	for (size_t bidx = 0; bidx < bin_schemes.size(); ++bidx) {
		bin_width = (max_data - min_data) / (int) bin_schemes[bidx];
		ull_int offset = bidx * data.size();
		hbin_t ti;
		for (size_t i = 0; i < data.size(); ++i) {
			ti = (hbin_t) floor((data[i] - min_data) / bin_width);
			data_int[offset + i] = (ti == bin_schemes[bidx]) ? ti - 1 : ti;
		}
	}

	*min_d = min_data;
	*max_d = max_data;
	*avg_d = avg_data;

	shuffleArray(rand_indices, data_int);

	return 0;
}

int Real2Int::dihReal2Int(const u_int id,
		const std::vector<hbin_t> &bin_schemes, const hbin_t bins_ref,
		const bool optimize, const std::vector<data_t> &data, double *min_d,
		double *max_d, double *avg_d, double *phase_d,
		std::vector<hbin_t> &data_int) {

	if (data.size() * bin_schemes.size() != data_int.size()) {
		return -1;
	}
	hbin_t nbins_max = bins_ref;
//   for (size_t i = 0; i < bin_schemes.size(); ++i) {
//      if (bin_schemes[i] > nbins_max) {
//         nbins_max = bin_schemes[i];
//      }
//   }
	const double _PI = Constants::PI;
	const double _2PI = 2 * _PI;
	const double dx = DELTA_X_RANGE_EXTEND;
	double modfit_min, modfit_max;
	double modfit_displac;
	hbin_t modfitnbins = (hbin_t) 0;
	int dpind;

	double min_data = 1.0e100, max_data = -1.0e100, avg_data = 0.0;
	double min_data1 = 1.0e100, max_data1 = -1.0e100;
	double min_data2 = 1.0e100, max_data2 = -1.0e100;
	double min_data3 = 1.0e100, max_data3 = -1.0e100;
	double sum_data = 0.0, bin_width = 0.0, binsize = 0.0;

	modfitnbins = ((3 * (int) nbins_max) > 180) ? 180 : (3 * (int) nbins_max);

	max_data1 = 0.0;
	min_data1 = 6.5;
	max_data2 = 0.0;
	min_data2 = 6.5;
	max_data3 = 0.0;
	min_data1 = 6.5;
	avg_data = 0.0;
	dpind = 0;

	if (optimize) {
		modfit_min = 0.0;
		modfit_max = _2PI + dx;

		binsize = (modfit_max - modfit_min) / modfitnbins;

		for (size_t i = 0; i < data.size(); ++i) {
			int t1 = (int) fabs(floor((data[i] - modfit_min) / binsize));
			data_int[i] = t1;
		}

		/***************************************************************
		 * ModfitDih1D routines included here
		 **************************************************************/
		std::vector<ull_int> prob(modfitnbins, 0);
		std::vector<double> bincenter(modfitnbins);
		std::vector<bool> zerohist(2 * modfitnbins);

		u_int sumzerohist;
		double probrecord;
		u_int minidx;
		u_int countzeros;
		u_int countzerosrecord;
		u_int firstzero;

		// !
		// ! Loop over the number of convergence steps.
		// !
		//
		// assigning each dihedral to its bin.

		for (ull_int i = 0; i < data.size(); i++) {
			prob[data_int[i]]++;
		}

		//
		//  Having determined the number of data points included in each bin, now we can estimate the
		//  probability density.
		//
		//  Loop through the number of bins to "integrate" the probability density.
		//

		sumzerohist = 0;
		for (int i = 0; i < modfitnbins; ++i) {
			bincenter[i] = modfit_min + (i + 0.5) * binsize;
			// zerohist: Put false for POPULATED bin, true if it is EMPTY
			if (prob[i] > 0) {
				zerohist[i] = false;
			} else {
				zerohist[i] = true;
				sumzerohist++;
			}
			zerohist[modfitnbins + i] = false; // initialize duplicate of zerohist
		}
		//
		// For the case when all histogram bins have some population:
		// Set edge at minimum density
		if (sumzerohist == 0) {
			// Get minidx: index of bin with minimum density
			probrecord = prob[0];
			minidx = 0;
			for (int j = 1; j < modfitnbins; j++) {
				if (prob[j] < probrecord) {
					probrecord = prob[j];
					minidx = j;
				}
			}
		} else {
			// For the case where there are unoccupied bins: Set edge at
			// the beginning of the largest unoccupied span
			// Find first occurrence of an unoccupied bin
			// (prob(j)==0) or (zerohist[j] == false)
			int j = 0;
			while ((zerohist[j] == false) && (j < modfitnbins)) {
				++j;
			}
			minidx = j;
			firstzero = j;
			// Duplicate to cover edges in circular variable
			for (int k = 0; k < modfitnbins; k++) {
				zerohist[k + modfitnbins] = zerohist[k];
			}
			// Find the largest stretch of unoccupied bins
			countzeros = 0;
			countzerosrecord = 0;
			for (j = firstzero; j < 2 * modfitnbins; j++) {
				if (zerohist[j] == true) {
					countzeros++;
				} else {
					if (countzeros > countzerosrecord) {
						minidx = j - countzeros;
						countzerosrecord = countzeros;
					}
					countzeros = 0;
				}
			}
		}
		// modfitdisplac holds the center of the histogram where the
		// largest continuous unoccupied stretch begins, or where the
		// minimum density occurs (for the case that all bins are occupied)
		modfit_displac = _2PI - bincenter[minidx];
	}

	/**
	 * Loop through all individual dihedral array binary files and read
	 * each. First allocate the size of dvec read the dihedral arrays
	 * find the current dihedral extrema mina/maxa holds normal range,
	 * minb/maxb holds range of shifted distribution (by 180 deg) find
	 * out the smaller of the 2 ranges
	 */

	/***************************************************************
	 * DihExtrema calculation starts here
	 **************************************************************/

	/**
	 * Find the smallest and largest value of each degree of freedom Then
	 * the range is shifted by pi and recomputed. This allows for narrow
	 * distributions centered near 2pi to be better represented.
	 */

	for (size_t i = 0; i < data.size(); ++i) {
		sum_data += data[i];
		if (min_data1 > data[i]) {
			min_data1 = data[i];
		}
		if (max_data1 < data[i]) {
			max_data1 = data[i];
		}

		double data_PI_shifted = fmod(data[i] + _PI, _2PI);
		if (data_PI_shifted < min_data2) {
			min_data2 = data_PI_shifted;
		}
		if (data_PI_shifted > max_data2) {
			max_data2 = data_PI_shifted;
		}

		if (optimize) {
			double data_optimized = fmod(data[i] + modfit_displac, _2PI);
			if (data_optimized < min_data3) {
				min_data3 = data_optimized;
			}
			if (data_optimized > max_data3) {
				max_data3 = data_optimized;
			}
		}
	}

	avg_data = sum_data / data.size();

	if (!optimize) {
		if (fabs((max_data1 - min_data1) - (max_data2 - min_data2)) <= 1.0e-1) {
			min_data = min_data1;
			max_data = max_data1;
			dpind = 0;
		} else {
			if ((max_data1 - min_data1) < (max_data2 - min_data2)) {
				min_data = min_data1;
				max_data = max_data1;
				dpind = 0;
			} else {
				min_data = min_data2;
				max_data = max_data2;
				dpind = 1;
			}
		}
	} else {
		min_data = min_data3;
		max_data = max_data3;
		dpind = 2;
#ifdef DEBUG_CODE
		std::cout.precision(10);
		std::cout << std::fixed << min_data3 << " " << std::fixed << max_data3 << std::endl;
#endif
	}
	/**
	 * Extend the dihedral extrema by 0.05 degrees and finish computing the
	 * average
	 */
	min_data -= dx;
	if (min_data < 0.0) {
		min_data = 0.0;
	}
	max_data += dx;
	if (max_data > _2PI) {
		max_data = _2PI;
	}

	if ((max_data - min_data) < 1.0e-4) {
		max_data += 0.05;
	}
	if (dpind == 0) {
		avg_data = sum_data / data.size();
	} else if (dpind == 1) {
		avg_data = fmod((sum_data / data.size()) + _PI, _2PI);
	} else if (dpind == 2) {
		avg_data = fmod((sum_data / data.size()) + modfit_displac, _2PI);
	}

	/**
	 * Write the extrema for each dihdral. write the extrema information to
	 * file convert dihedrals to integer vector Deallocate Array
	 */

	if (!optimize) {
		*min_d = min_data;
		*max_d = max_data;
		*avg_d = avg_data;
		*phase_d = dpind * _PI;

	} else {
		*min_d = min_data;
		*max_d = max_data;
		*avg_d = avg_data;
		*phase_d = modfit_displac;
	}

	/***************************************************************
	 * DihExtrema calculation ends here
	 **************************************************************/
	/**
	 * Converts the double precision dihedral values to an integer array
	 * indicating bin which containers diherdal snapshot.
	 */
	for (size_t bidx = 0; bidx < bin_schemes.size(); ++bidx) {
		hbin_t nbins_loc = bin_schemes[bidx];
		ull_int offset = bidx * data.size();
		double bin_width1 = (max_data - min_data) / (int) nbins_loc;

		for (size_t jj = 0; jj < data.size(); ++jj) {
			hbin_t tt;
			if (dpind == 0) {
				tt = (hbin_t) floor((data[jj] - min_data) / bin_width1);
				data_int[offset + jj] = (tt == nbins_loc) ? tt - 1 : tt;
			}
			if (dpind == 1) {
				tt = (hbin_t) floor(
						(fmod(data[jj] + _PI, _2PI) - min_data) / bin_width1);
				data_int[offset + jj] = (tt == nbins_loc) ? tt - 1 : tt;
			}
			if (dpind == 2) {
				tt = (hbin_t) floor(
						(fmod(data[jj] + modfit_displac, _2PI) - min_data)
								/ bin_width1);
				data_int[offset + jj] = (tt == nbins_loc) ? tt - 1 : tt;
			}
		}
	}
	return 0;
}

int Real2Int::getMinimaPositions(const u_int id, bool analyticGrad,
		hbin_t n_conf_max, double k_value, u_int step_sdo,
		double conv_criterion, u_int max_iterations, std::vector<data_t> &data,
		std::vector<double> &minimum, hbin_t &found_mins,
		std::vector<hbin_t> &data_int) {

	double coordmaxs[9];    // Positions of the maximums in the PDF
	double secderivmaxs[9] = { 0.0 }; // Positions of the maximums of the 2nd derivate of the (PDF)

	double density[361];    // Probability Density Function (PDF)
	double gradient[361];   // Gradient of the PDF
	double secderiv[361];   // Second derivate of the PDF
	double smoothingparam;  // Smoothing parameter (see ref. 2)
	double vMgradient; // Gradient in a general position (not only in {0,...359})
	double vMdensity;       // Density in a general position
	double vMsec_deriv;     // Second derivate in a general position

	double step_times_grad;   //  Step x gradient
	int n_obsmax;   //  Observed Number of Maximums

	int contador;
	int Itrs, m;
	double temp;
	bool Converged;
	double grad;

	int result = 0;
	ull_int NumSnap = data.size();

	// Computation of distributions Bandwith
	smoothingparam = pow(
			(3 * NumSnap * (k_value * k_value)
					* Utils::modifiedBesselI(2, 2 * k_value)
					/ (4 * Constants::SQRTPI
							* pow(Utils::modifiedBesselI(0, k_value), 2.0))),
			(2.0 / 5.0));

//
// Calculating densities, gradient, second derivate
// and looking for local maximums
	n_obsmax = 0;   //    Observed Number of Maximums
	for (int i = 0; i <= 360; ++i) {
		vonMises(smoothingparam, (i * Constants::DEGRAD), data, vMdensity,
				vMgradient, vMsec_deriv);
		gradient[i] = vMgradient;
		density[i] = vMdensity;
		secderiv[i] = vMsec_deriv;
		if (i > 0) {
			if (((gradient[i - 1] * gradient[i]) <= 0.0)
					&& (vMsec_deriv < (-1.0e-4)) && (n_obsmax < 9)) {
				coordmaxs[n_obsmax] = i / 2.0 + (i - 1) / 2.0; // Approximate coord of the maximum
				secderivmaxs[n_obsmax] = secderiv[i] / 2.0
						+ secderiv[i - 1] / 2.0; //   Approximate second derivate
				n_obsmax += 1; //
			}
		}
	}
	//
	// Sorting to give preference to those maximums with
	// lower(more negative) second derivate
	for (int i = 0; i < n_obsmax - 1; ++i) {
		for (int j = i + 1; j < n_obsmax; ++j) {
			if (secderivmaxs[i] > secderivmaxs[j]) {
				temp = secderivmaxs[i];
				secderivmaxs[i] = secderivmaxs[j];
				secderivmaxs[j] = temp;
				temp = coordmaxs[i];
				coordmaxs[i] = coordmaxs[j];
				coordmaxs[j] = temp;
			}
		}
	}

// In a circular variable, the Max. number of minimums
// most be equal to the Max. number of maximums, we get
// the lowest value between the predefined MaxNumConf
// and the observed number of maximums/minimums (ObsNmax)
	if (n_conf_max < n_obsmax) {
		found_mins = n_conf_max;
	} else if (n_conf_max >= n_obsmax) {
		found_mins = n_obsmax;
	}
//
// Sorting only the first FNumMins using the coordinates
	for (int i = 0; i < found_mins - 1; ++i) {
		for (int j = i + 1; j < found_mins; ++j) {
			if (coordmaxs[i] > coordmaxs[j]) {
				temp = coordmaxs[i];
				coordmaxs[i] = coordmaxs[j];
				coordmaxs[j] = temp;
			}
		}
	}
//
// Guessing where are the minimums
	if (found_mins == 1) {
		minimum[0] = (coordmaxs[0] - 180.0);
		if (minimum[0] < 0.0) {
			minimum[0] = minimum[0] + 360.0;
		}
	} else {
		for (int i = 0; i < found_mins; ++i) {
			if (i == 0) {
				minimum[i] = (coordmaxs[i] + coordmaxs[found_mins - 1] - 360.0)
						/ 2.0;
				if (minimum[i] < 0.0) {
					minimum[i] = minimum[i] + 360;
				}
			} else {
				minimum[i] = (coordmaxs[i] + coordmaxs[i - 1]) / 2.0;
			}
		}
	}
//
// Optimizing Minimums
	if (found_mins > 1) {
		for (int m = 0; m < found_mins; ++m) {
			Itrs = 0;
			Converged = false;
			while ((!Converged) && (Itrs < max_iterations)) {
				Itrs = Itrs + 1;
				if (analyticGrad) {
					grad = calculateVonMisesGradient(smoothingparam,
							Constants::DEGRAD * minimum[m], data);
				} else {
					//  !Gradient Linear interpolation
					grad = (fabs(minimum[m] - 1 - int(minimum[m])))
							* gradient[int(minimum[m])]
							+ (fabs(minimum[m] - int(minimum[m]))
									* gradient[1 + int(minimum[m])]);
				}
				if ((fabs(grad) < conv_criterion)) {
					Converged = true;
				} else {
					if (fabs(step_sdo * grad) > 5.0) {
						step_times_grad = ((step_sdo * grad)
								/ fabs(step_sdo * grad)) * 5.0;
					} else {
						step_times_grad = step_sdo * grad;
					}
					minimum[m] = minimum[m] - step_times_grad;
					if (minimum[m] < 0.0) {
						minimum[m] = minimum[m] + 360.0;
					}
					if (minimum[m] >= 360.0) {
						minimum[m] = minimum[m] - 360.0;
					}
				}
			}
		}
	}

	for (int i = 0; i < found_mins - 1; ++i) {
		for (int j = i + 1; j < found_mins; ++j) {
			if (minimum[i] > minimum[j]) {
				temp = minimum[i];
				minimum[i] = minimum[j];
				minimum[j] = temp;
			}
		}
	}

	for (int i = 0; i < found_mins; ++i) {
		minimum[i] *= Constants::DEGRAD;
	}

	//-codify
	for (int k = 0; k < NumSnap; ++k) {
		if (found_mins == 1) {
			data_int[k] = 0;
		} else {
			for (int j = 0; j < found_mins; ++j) {
				if (j == 0) {
					if ((data[k] < minimum[j])
							|| (data[k] > minimum[found_mins - 1])) {
						data_int[k] = j;
					}
				} else {
					if ((data[k] < minimum[j]) && (data[k] > minimum[j - 1])) {
						data_int[k] = j;
					}
				}
			}
		}
	}

	return result;

}

double Real2Int::vonMises(double SmoothingParam, double X,
		std::vector<data_t> &ang, double &vMdensity, double &vMgradient,
		double &vMsec_deriv) {

	double aj;                   // Data Angles in radians
	double expo;                 // Auxiliary real variable
	ull_int n_data = ang.size();

	vMgradient = 0.0;
	vMgradient = 0.0;
	vMsec_deriv = 0.0;
	double X_minus_aj = 0.0;
	double cos_of_X_minus_aj = 0.0;
	double sin_of_X_minus_aj = 0.0;
	double sin_sqr_of_X_minus_aj = 0.0;
	//double piover2 = std::asin(1.0);
	for (size_t j = 0; j < n_data; ++j) {
		aj = ang[j];
		X_minus_aj = X - aj;
		cos_of_X_minus_aj = std::cos(X_minus_aj);
		sin_of_X_minus_aj = std::sin(X_minus_aj);
		sin_sqr_of_X_minus_aj = sin_of_X_minus_aj * sin_of_X_minus_aj;
		//cos_of_X_minus_aj = sqrt(1 - sin_sqr_of_X_minus_aj);
		//cos_of_X_minus_aj= (abs(X_minus_aj) > piover2) ? -cos_of_X_minus_aj: cos_of_X_minus_aj;

		expo = exp(SmoothingParam * cos_of_X_minus_aj);

		vMdensity = vMdensity + expo;
		vMgradient = vMgradient - sin_of_X_minus_aj * expo;
		vMsec_deriv = vMsec_deriv
				+ (SmoothingParam * sin_sqr_of_X_minus_aj - cos_of_X_minus_aj)
						* expo;
	}
	double denom = (Constants::TWOPI * n_data
			* Utils::modifiedBesselI(0, SmoothingParam));
	vMdensity = vMdensity / denom;
	vMgradient = vMgradient * SmoothingParam / denom;
	vMsec_deriv = vMsec_deriv * SmoothingParam / denom;
	return vMsec_deriv;
}

double Real2Int::calculateVonMisesGradient(double SmoothingParam, double X,
		std::vector<data_t> &ang) {
	double result = 0.0;
	ull_int n_data = ang.size();
	double cos_of_X_minus_aj = 0.0;
	double sin_of_X_minus_aj = 0.0;
	for (int j = 0; j < n_data; ++j) {
		cos_of_X_minus_aj = cos(X - ang[j]);
		sin_of_X_minus_aj = sqrt(1.0 - (cos_of_X_minus_aj * cos_of_X_minus_aj));

		result -= SmoothingParam * sin_of_X_minus_aj
				* exp(SmoothingParam * cos_of_X_minus_aj)
				/ (Constants::TWOPI * n_data
						* Utils::modifiedBesselI(0, SmoothingParam));
	}
	return result;
}

