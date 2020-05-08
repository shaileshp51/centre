/*!\brief This file defines basic functions used for calculating bin frequency and entropy.
 *
 * File:   Estimators.h
 * Author: shailesh
 *
 * Created on Created on: 30-Sep-2016
 */

#ifndef ESTIMATORS_H_
#define ESTIMATORS_H_

#include <string>

/*!
 \fn std::string getEstimatorName(EntropyEstimators est)

 A function go get the name[string] of EntropyEstimator from enum value.
 \param[in] est the entropy estimator

 Return the name string of entropy estimator.
 */
std::string getEstimatorName(EntropyEstimators est);

/*!
 \fn std::string getMethodName(EntropyMethods mthd)

 A function go get the name[string] of EntropyMethods from enum value.
 \param[in] mthd the EntropyMethods value

 Return the name string of entropy method.
 */
std::string getMethodName(EntropyMethods mthd);

size_t get2D21DIndex(const hbin_t nD, const hbin_t *indices,
		const hbin_t dim_len);

/*!
 \fn int binData1D(const hbin_t nsteps, const std::vector<hbin_t>& bin_shemes,
 const ull_int nframes, const double min_d, const double max_d,
 hbin_t* data, std::vector<ull_int> &freqs_obs,
 std::vector<double> &bin_mids)

 A function go calculate observed bin frequency.
 \param[in] nsteps the number of steps for entropy estimation
 \param[in] bin_shemes a vector of bining schemes used for discretization
 \param[in] min_d the minimum value observed  in the data
 \param[in] max_d the maximum value observed  in the data
 \param[in] *data a pointer to the data array
 \param[out] freq_obs a vector of observed bin frequencies
 \param[out] bin_mids a vector of mind points of bins

 Return \c 0 on success; non-zero otherwise.
 */
int binData1D(const hbin_t nsteps, const std::vector<hbin_t> &bin_shemes,
		const ull_int nframes, const double min_d, const double max_d,
		hbin_t *data, std::vector<ull_int> &freqs_obs,
		std::vector<double> &bin_mids);

int binData1DS(hbin_t nsteps, hbin_t nbins, ull_int nframes, double min_d,
		double max_d, hbin_t *data, std::vector<ull_int> *freqs_obs);

/*!
 \fn int binData2D(const hbin_t nsteps, const std::vector<hbin_t>& bin_shemes,
 const ull_int nframes, const double min_d1, const double max_d1,
 const double min_d2, const double max_d2, hbin_t* data1, hbin_t* data2,
 std::vector<ull_int>& freqs_obs, std::vector<double>& bin_mids)

 A function go calculate observed bin frequency.
 \param[in] nsteps the number of steps for entropy estimation
 \param[in] bin_shemes a vector of bining schemes used for discretization
 \param[in] min_d1 the minimum value observed  in the data1
 \param[in] max_d1 the maximum value observed  in the data1
 \param[in] min_d2 the minimum value observed  in the data2
 \param[in] max_d2 the maximum value observed  in the data2
 \param[in] *data1 a pointer to the data1 array
 \param[in] *data2 a pointer to the data2 array
 \param[out] freq_obs a vector of observed bin frequencies
 \param[out] bin_mids a vector of mind points of bins

 Return \c 0 on success; non-zero otherwise.
 */
int binData2D(const hbin_t nsteps, const std::vector<hbin_t> &bin_shemes,
		const ull_int nframes, const double min_d1, const double max_d1,
		const double min_d2, const double max_d2, hbin_t *data1, hbin_t *data2,
		std::vector<ull_int> &freqs_obs, std::vector<double> &bin_mids);

/*!
 \fn int entropy1D(const BAT_t type_d, const u_int id,
 const std::vector<EntropyEstimators>& ests, const hbin_t nsteps,
 const ull_int step_size, const double min_d, const double max_d,
 const bool jacobian, const bool kde,
 const std::vector<hbin_t> bin_schemes,
 const std::vector<ull_int>& freqs_obs, std::vector<double>& entropy)

 A function go calculate observed bin frequency.
 \param[in] type_d the BAT type of DOF
 \param[in] id the id of DOF
 \param[in] ests a veco of entropy estimators to be used
 \param[in] nsteps the number of steps for entropy estimation
 \param[in] step_size the number of points in each step
 \param[in] min_d the minimum value observed  in the data
 \param[in] max_d the maximum value observed  in the data
 \param[in] jacobian a bool to instruct whether to use jacobian in calculation
 \param[in] kde a bool true means kde method is used otherwise histogram
 \param[in] bin_shemes a vector of bining schemes used for discretization
 \param[out] freq_obs a vector of observed bin frequencies
 \param[out] entropy a vector of entropy values

 Return \c 0 on success; non-zero otherwise.
 */
int entropy1D(const BAT_t type_d, const u_int id,
		const std::vector<EntropyEstimators> &ests, const hbin_t nsteps,
		const ull_int step_size, const double min_d, const double max_d,
		const bool jacobian, const bool kde,
		const std::vector<hbin_t> bin_schemes,
		const std::vector<ull_int> &freqs_obs, std::vector<double> &entropy);

/*!
 \fn int binData2D(const hbin_t nsteps, const std::vector<hbin_t>& bin_shemes,
 const ull_int nframes, const double min_d1, const double max_d1,
 const double min_d2, const double max_d2, hbin_t* data1, hbin_t* data2,
 std::vector<ull_int>& freqs_obs, std::vector<double>& bin_mids)

 A function go calculate observed bin frequency.
 \param[in] type_d1 the BAT type of first DOF in the pair of DOFs
 \param[in] type_d2 the BAT type of second DOF in the pair of DOFs
 \param[in] id1 the id of first DOF in the pair of DOFs
 \param[in] id2 the id of second DOF in the pair of DOFs
 \param[in] ests a veco of entropy estimators to be used
 \param[in] nsteps the number of steps for entropy estimation
 \param[in] step_size the number of points in each step
 \param[in] min_d1 the minimum value observed  in the data1
 \param[in] max_d1 the maximum value observed  in the data1
 \param[in] min_d2 the minimum value observed  in the data2
 \param[in] max_d2 the maximum value observed  in the data2
 \param[in] jacobian a bool to instruct whether to use jacobian in calculation
 \param[in] kde a bool true means kde method is used otherwise histogram
 \param[in] bin_shemes a vector of bining schemes used for discretization
 \param[out] freq_obs a vector of observed bin frequencies
 \param[out] entropy a vector of entropy values

 Return \c 0 on success; non-zero otherwise.
 */
int entropy2D(const BAT_t type_d1, const BAT_t type_d2, const u_int id1,
		const u_int id2, const std::vector<EntropyEstimators> &ests,
		const hbin_t nsteps, const ull_int step_size, const double min_d1,
		const double max_d1, const double min_d2, const double max_d2,
		const bool jacobian, const bool kde,
		const std::vector<hbin_t> bin_schemes,
		const std::vector<ull_int> &freqs_obs, std::vector<double> &entropy);

#endif /* ESTIMATORS_H_ */
