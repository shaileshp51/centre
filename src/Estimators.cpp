/*
 * Estimators.cpp
 *
 *  Created on: 26-Sep-2016
 *      Author: shailesh
 */
#include <cmath>
#include <vector>

#include <omp.h>

#include "Centre.h"
#include "Estimators.h"

string getEstimatorName(EntropyEstimators est) {
	std::string res = "ML";
	switch (est) {
	case EntropyEstimators::ML:
		res = "ML";
		break;
	case EntropyEstimators::MM:
		res = "MM";
		break;
	case EntropyEstimators::CS:
		res = "CS";
		break;
	case EntropyEstimators::JS:
		res = "JS";
		break;
	}
	return res;
}

string getMethodName(EntropyMethods mthd) {
	std::string res = "MIE";
	switch (mthd) {
	case EntropyMethods::MIE:
		res = "MIE";
		break;
	case EntropyMethods::MIST:
		res = "MIST";
		break;
	case EntropyMethods::AMIE:
		res = "AMIE";
		break;
	case EntropyMethods::AMIST:
		res = "AMIST";
		break;
	case EntropyMethods::QH:
		res = "QH";
		break;
	}
	return res;
}

int binData1D(const hbin_t nsteps, const std::vector<hbin_t> &bin_shemes,
		const ull_int nframes, const double min_d, const double max_d,
		hbin_t *data, std::vector<ull_int> &freqs_obs,
		std::vector<double> &bin_mids) {
	size_t nbins_sum = 0;
	for (size_t i = 0; i < bin_shemes.size(); ++i) {
		nbins_sum += bin_shemes[i];
	}
	if (nsteps * nbins_sum != freqs_obs.size()) {
		return -1;
	} else {
		int binid_off = 0;
		for (size_t tmp = 0; tmp < bin_shemes.size(); ++tmp) {
			int nbins = bin_shemes[tmp];
			double bin_width = (max_d - min_d) / nbins;
			double half_bin_width = bin_width / 2.0;
			double bin_offset_dbl = min_d;
			for (int i = 0; i < nbins; ++i) {
				bin_mids[binid_off + i] = bin_offset_dbl + half_bin_width;
				bin_offset_dbl += bin_width;
			}
			binid_off += nbins;
		}
#ifdef USE_OMPMPI
		int num_threads = omp_get_max_threads();
#endif
		ull_int step_size = nframes / nsteps;
		ull_int frameid_shift = 0;
		u_int binid_shift = 0;
		for (u_int schemes_id = 0; schemes_id < bin_shemes.size();
				++schemes_id) {
			hbin_t nbins = bin_shemes[schemes_id];
			for (size_t i = 0; i < nsteps; ++i) {
				ull_int offset = i * step_size;
				int offset_freqs = (nbins_sum * i) + binid_shift;

#ifdef USE_OMPMPI
				std::vector<std::vector<int>> local_freqs(num_threads + 1,
						std::vector<int>(nbins));
#pragma omp parallel
				{
					int tid = omp_get_thread_num();
					// __declspec (align(64)) int local_freqs[num_threads + 1][(int)nbins];

#pragma omp for nowait
					for (size_t j = 0; j < step_size; ++j) {
						hbin_t bin = (*(data + frameid_shift + offset + j));
						// std::cout << (int) bin << " ";
						local_freqs[tid][bin]++;
					}
					// std::cout << std::endl << "hist freq thread step(" << i << ") ";
				}
//#pragma omp parallel
//			{
//#pragma omp for
				for (size_t j = 0; j < nbins; ++j) {
					for (size_t t = 0; t < num_threads; ++t) {
						// std::cout << local_freqs[t][j] << " ";
						freqs_obs[offset_freqs + j] += local_freqs[t][j];
					}
				}
//			}

#else
				hbin_t bin;
				for (size_t j = 0; j < step_size; ++j) {
					bin = (*(data + frameid_shift + offset + j));
					freqs_obs[offset_freqs + bin] += 1;
				}
#endif
			}

			for (size_t i = 1; i < nsteps; ++i) {
				int off = (nbins_sum * i) + binid_shift;
				for (size_t j = 0; j < nbins; ++j) {
					freqs_obs[off + j] += freqs_obs[off - nbins_sum + j];
				}
			}
#ifdef DEBUG_CODE_2
			std::cout << "Scheme: (" << schemes_id << ", " << (int) nbins << ")" << std::endl;
			for (size_t i = 0; i < nsteps; ++i) {
				int off = (nbins_sum * i) + binid_shift;
				std::cout << std::endl << "hist freq: final step  {";
				for (size_t j = 0; j < nbins; ++j) {
					std::cout << "(" << off + j << ", " << freqs_obs[off + j] << "), ";
				}
				std::cout << "} " << std::endl;
			}
#endif
			frameid_shift += nframes;
			binid_shift += nbins;
		}
	}
	return 0;
}

size_t get2D21DIndex(const hbin_t nD, const hbin_t *indices,
		const hbin_t dim_len) {
	size_t result = 0;
	for (int i = 0; i < nD - 1; ++i) {
		result += *(indices + i) * pow(dim_len, nD - i - 1);
	}
	result += *(indices + nD - 1);
	return result;
}

int binData2D(const hbin_t nsteps, const std::vector<hbin_t> &bin_shemes,
		const ull_int nframes, const double min_d1, const double max_d1,
		const double min_d2, const double max_d2, hbin_t *data1, hbin_t *data2,
		std::vector<ull_int> &freqs_obs, std::vector<double> &bin_mids) {
	size_t nbins_sum_sq = 0;
	for (size_t i = 0; i < bin_shemes.size(); ++i) {
		nbins_sum_sq += bin_shemes[i] * bin_shemes[i];
	}
	if (nsteps * nbins_sum_sq != freqs_obs.size()) {
		return -1;
	} else {
		int binid_off = 0;
		for (size_t tmp = 0; tmp < bin_shemes.size(); ++tmp) {
			int nbins = bin_shemes[tmp];
			double bin_width1 = (max_d1 - min_d1) / nbins;
			double half_bin_width1 = bin_width1 / 2;
			double bin_offset_dbl1 = min_d1;
			double bin_width2 = (max_d2 - min_d2) / nbins;
			double half_bin_width2 = bin_width2 / 2;
			double bin_offset_dbl2 = min_d2;
			int binid_off2 = binid_off + nbins;
			for (int i = 0; i < nbins; ++i) {
				bin_mids[binid_off + i] = bin_offset_dbl1 + half_bin_width1;
				bin_offset_dbl1 += bin_width1;

				// bin mids for second data: data2
				bin_mids[binid_off2 + i] = bin_offset_dbl2 + half_bin_width2;
				bin_offset_dbl2 += bin_width2;
			}
			binid_off += 2 * nbins;
		}
#ifdef USE_OMPMPI
		int num_threads = omp_get_max_threads();
#endif
		ull_int step_size = nframes / nsteps;
		ull_int frameid_shift = 0;
		u_int binid_shift = 0;
		for (u_int schemes_id = 0; schemes_id < bin_shemes.size();
				++schemes_id) {
			hbin_t nbins = bin_shemes[schemes_id];
			size_t nbins_sq = nbins * nbins;
			for (size_t i = 0; i < nsteps; ++i) {
				ull_int offset = i * step_size;
				int offset_freqs = (nbins_sum_sq * i) + binid_shift;

				// int local_freqs[num_threads + 1][(int) nbins][(int) nbins] { };
#ifdef USE_OMPMPI
				std::vector<std::vector<int>> local_freqs(num_threads + 1,
						std::vector<int>(nbins_sq));
#pragma omp parallel
				{
					int tid = omp_get_thread_num();
					// __declspec (align(64)) int local_freqs[num_threads + 1][(int)nbins];
					hbin_t bin_ids[2];
#pragma omp for nowait
					for (size_t j = 0; j < step_size; ++j) {
						bin_ids[0] = (*(data1 + frameid_shift + offset + j));
						bin_ids[1] = (*(data2 + frameid_shift + offset + j));
						// local_freqs[tid][bin1][bin2]++;
						local_freqs[tid][get2D21DIndex(2, &bin_ids[0], nbins)]++;
					}
					// std::cout << std::endl << "hist freq thread step(" << i << ") ";
				}
//#pragma omp parallel
//			{
//#pragma omp for
				for (size_t j = 0; j < nbins_sq; ++j) {
					for (size_t t = 0; t < num_threads; ++t) {
						// std::cout << local_freqs[t][j] << " ";
						freqs_obs[offset_freqs + j] += local_freqs[t][j];
					}
				}
//			}
#else
				hbin_t bin_ids[2];
				size_t d1_index;
				for (size_t j = 0; j < step_size; ++j) {
					bin_ids[0] = (*(data1 + frameid_shift + offset + j));
					bin_ids[1] = (*(data2 + frameid_shift + offset + j));
					d1_index = get2D21DIndex(2, &bin_ids[0], nbins);
					freqs_obs[offset_freqs + d1_index] += 1;
				}
#endif
			}

			for (size_t i = 1; i < nsteps; ++i) {
				int off = (nbins_sum_sq * i) + binid_shift;
				for (size_t j = 0; j < nbins_sq; ++j) {
					freqs_obs[off + j] += freqs_obs[off - nbins_sum_sq + j];
				}
			}
#ifdef DEBUG_CODE_2
			std::cout << "Scheme: (" << schemes_id << ", " << (int) nbins << ")" << std::endl;
			for (size_t i = 0; i < nsteps; ++i) {
				int off = (nbins_sum_sq * i) + binid_shift;
				std::cout << std::endl << "hist freq: final step  {";
				for (size_t j = 0; j < nbins_sq; ++j) {
					std::cout << "(" << off + j << ", " << freqs_obs[off + j] << "), ";
				}
				std::cout << "} " << std::endl;
			}
#endif
			frameid_shift += nframes;
			binid_shift += nbins_sq;
		}
	}
	return 0;
}

size_t getND21DIndex(const hbin_t nD, const hbin_t *indices,
		const hbin_t dim_len) {
	size_t result = 0;
	for (int i = 0; i < nD - 1; ++i) {
		result += *(indices + i) * pow(dim_len, nD - i - 1);
	}
	result += *(indices + nD - 1);
	return result;
}

int binDataND(const hbin_t nD, const hbin_t nsteps,
		const std::vector<hbin_t> &bin_schemes, const ull_int nframes,
		const std::vector<double> min_d, const std::vector<double> max_d,
		std::vector<hbin_t*> data_pv, std::vector<ull_int> *freqs_obs,
		std::vector<double> *bin_mids) {
	size_t nbins_powN_sum = 0;
	for (size_t i = 0; i < bin_schemes.size(); ++i) {
		nbins_powN_sum += pow(bin_schemes[i], nD);
	}

	if (nsteps * nbins_powN_sum != freqs_obs->size()) {
		return -1;
	} else {
		int binid_off = 0;
		for (size_t tmp = 0; tmp < bin_schemes.size(); ++tmp) {
			hbin_t nbins = bin_schemes[tmp];
			std::vector<double> freqs(nbins_powN_sum, 0.0);
			std::vector<double> bin_width(nD, 0.0);
			std::vector<double> bin_width_half(nD, 0.0);

			for (int i = 0; i < nD; ++i) {
				bin_width[i] = (max_d[i] - min_d[i]) / nbins;
				bin_width_half[i] = bin_width[i] / 2.0;
			}
			std::vector<double> bin_offset_dbl(min_d);
			for (int d = 0; d < nD; ++d) {
				int binid_off2 = binid_off + d * nbins;
				for (int i = 0; i < nbins; ++i) {
					(*bin_mids)[binid_off2 + i] = bin_offset_dbl[i]
							+ bin_width_half[i];
					bin_offset_dbl[i] += bin_width[i];
				}
			}
			binid_off += nD * nbins;
		}

		int num_threads = omp_get_max_threads();
		ull_int step_size = nframes / nsteps;
		ull_int frameid_shift = 0;
		u_int binid_shift = 0;
		for (u_int schemes_id = 0; schemes_id < bin_schemes.size();
				++schemes_id) {
			hbin_t nbins = bin_schemes[schemes_id];
			size_t nbins_powN = pow(nbins, nD);
			for (size_t i = 0; i < nsteps; ++i) {
				// ull_int offset = i * step_size;
				int offset_freqs = (nbins_powN_sum * i) + binid_shift;

				std::vector<std::vector<int>> local_freqs(num_threads + 1,
						std::vector<int>(nbins_powN));

#pragma omp parallel
				{
					int tid = omp_get_thread_num();
					// __declspec (align(64)) int local_freqs[num_threads + 1][(int)nbins];

#pragma omp for nowait
					for (size_t j = 0; j < step_size; ++j) {
						hbin_t bin_ids[nD];
						for (int did = 0; did < nD; ++did) {
							hbin_t bin1 = (*(data_pv[did] + frameid_shift + j));
							bin_ids[did] = bin1;
						}

						// std::cout << (int) bin << " ";
						local_freqs[tid][getND21DIndex(nD, &bin_ids[0], nbins)]++;
					}
					// std::cout << std::endl << "hist freq thread step(" << i << ") ";
				}
//#pragma omp parallel
//          {
//#pragma omp for
				for (size_t j = 0; j < nbins_powN; ++j) {
					for (size_t t = 0; t < num_threads; ++t) {
						// std::cout << local_freqs[t][j] << " ";
						(*freqs_obs)[offset_freqs + j] += local_freqs[t][j];
					}
				}
//          }
			}
			for (size_t i = 1; i < nsteps; ++i) {
				int off = (nbins_powN_sum * i) + binid_shift;
				for (size_t j = 0; j < nbins_powN; ++j) {
					(*freqs_obs)[off + j] += (*freqs_obs)[off - nbins_powN_sum
							+ j];
				}
			}
#ifdef DEBUG_CODE
			std::cout << "Scheme: (" << schemes_id << ", " << (int) nbins << ")" << std::endl;
			for (size_t i = 0; i < nsteps; ++i) {
				int off = (nbins_powN_sum * i) + binid_shift;
				std::cout << std::endl << "hist freq: final step  {";
				for (size_t j = 0; j < nbins_powN; ++j) {
					std::cout << "(" << off + j << ", " << (*freqs_obs)[off + j] << "), ";
				}
				std::cout << "} " << std::endl;
			}
#endif
			frameid_shift += nframes;
			binid_shift += nbins_powN;
		}
	}
	return 0;
}

int entropy1D(const BAT_t type_d, const u_int id,
		const std::vector<EntropyEstimators> &ests, const hbin_t nsteps,
		const ull_int step_size, const double min_d, const double max_d,
		const bool jacobian, const bool kde,
		const std::vector<hbin_t> bin_schemes,
		const std::vector<ull_int> &freqs_obs, std::vector<double> &entropy) {
	u_int bin_schemes_sum = 0;
	hbin_t nestm = ests.size();
	for (size_t tmp = 0; tmp < bin_schemes.size(); ++tmp) {
		bin_schemes_sum += bin_schemes[tmp];
	}
	if (nsteps * bin_schemes_sum != freqs_obs.size()) {
		return -1;
	} else if (nestm * nsteps * bin_schemes.size() != entropy.size()) {
		return -2;
	} else {

#ifdef DEBUD_CODE_2
		std::cout << "calc entropy: " << type_d << ", " << id << ", " << (int) nsteps << ", " << step_size << ", " << min_d << ", " << max_d << std::endl;
#endif
		for (hbin_t estid = 0; estid < nestm; ++estid) {
			EntropyEstimators est = ests[estid];
			u_int estm_offset = estid * nsteps * bin_schemes.size();

			u_int frq_off = 0;
			for (size_t tmp = 0, ent_off = 0; tmp < bin_schemes.size();
					ent_off += nsteps, ++tmp) {
				hbin_t nbins = bin_schemes[tmp];
				std::vector<double> freqs(nbins, 0.0);
				double bin_width = (max_d - min_d) / nbins;
				double bin_width_half = bin_width / 2;

				if (kde) {
					bin_width = 1.0;
				}
				std::vector<double> bin_mid(nbins, 0.0);
				double sum_plnp = 0.0;
				double probdens = 0.0;
				double jacobian_factor = 1.0;
				bool has_jacob_fact = false;
				hbin_t blank_bins = (hbin_t) 0;
				hbin_t singleton_bins = (hbin_t) 0;
				if (type_d == BAT_t::BOND || type_d == BAT_t::ANGLE) {
					has_jacob_fact = true;
					for (size_t j = 0; j < nbins; ++j) {
						bin_mid[j] = min_d + j * bin_width + bin_width_half;
					}
				}

				for (size_t i = 0; i < nsteps; ++i) {
					ull_int offset = (i * bin_schemes_sum) + frq_off;
					double n_data = (i + 1.0) * step_size;
					for (size_t j = 0; j < nbins; ++j) {
						freqs[j] = freqs_obs[offset + j] / (n_data);
					}
					sum_plnp = 0.0;
					blank_bins = (hbin_t) 0;
					singleton_bins = (hbin_t) 0;
					switch (est) {
					case EntropyEstimators::ML:
						break;
					case EntropyEstimators::MM:
						for (size_t j = 0; j < nbins; ++j) {
							if (freqs_obs[offset + j] == 0) {
								++blank_bins;
							}
						}
						break;
					case EntropyEstimators::CS:
						for (size_t j = 0; j < nbins; ++j) {
							if (freqs_obs[offset + j] == 1) {
								++singleton_bins;
							}
						}
						// FREQUENCY_CORRECTION:: for Chao-Shen, Good-Turing correction
						for (size_t j = 0; j < nbins; ++j) {
							freqs[j] = (1 - (singleton_bins / n_data))
									* (freqs_obs[offset + j] / n_data);
							// freqs[j] /= bin_width;
						}
						break;
					case EntropyEstimators::JS:
						double lambda_hat = 0.0;
						double tk = 1.0 / (double) nbins;
						double sum_theta_ML2 = 0;
						double tk_sum_theta_ML2 = 0;

						for (size_t j = 0; j < nbins; ++j) {
							sum_theta_ML2 += freqs[j] * freqs[j];
							tk_sum_theta_ML2 += ((tk - freqs[j])
									* (tk - freqs[j]));
						}
						lambda_hat = (1 - sum_theta_ML2)
								/ ((n_data - 1) * tk_sum_theta_ML2);

						for (size_t j = 0; j < nbins; ++j) {
							freqs[j] = (lambda_hat * tk)
									+ ((1 - lambda_hat) * freqs[j]);
						}
						break;
					}
					if (est == EntropyEstimators::CS) {
						for (size_t j = 0; j < nbins; ++j) {
							probdens = freqs[j] / bin_width;
							if (probdens > 0.0) {
								// HT ESTIMATOR: Horowitz-Thompson estimators, denomirator
								double denom = 1 - pow((1 - freqs[j]), n_data);
								if (jacobian && has_jacob_fact) {
									switch (type_d) {
									case BAT_t::BOND:
										jacobian_factor = pow(bin_mid[j], 2);
										break;
									case BAT_t::ANGLE:
										jacobian_factor = sin(bin_mid[j]);
										break;
									default: // When it is a dihedral, associated jacobian is 1.0
										break;
									}

									probdens = probdens / jacobian_factor;
									sum_plnp += (jacobian_factor * probdens
											* log(probdens) / denom);
								} else {
									sum_plnp += (probdens * log(probdens)
											/ denom);
								}
							}
						}
						sum_plnp *= -bin_width;
					} else {
						for (size_t j = 0; j < nbins; ++j) {
							probdens = freqs[j] / bin_width;
							if (probdens > 0.0) {
								if (jacobian && has_jacob_fact) {
									switch (type_d) {
									case BAT_t::BOND:
										jacobian_factor = pow(bin_mid[j], 2);
										break;
									case BAT_t::ANGLE:
										jacobian_factor = sin(bin_mid[j]);
										break;
									default: // When it is a dihedral, associated jacobian is 1.0
										break;
									}
									probdens = probdens / jacobian_factor;
									sum_plnp += (jacobian_factor * probdens
											* log(probdens));

								} else {
									sum_plnp += (probdens * log(probdens));
								}
							}
						}
						sum_plnp *= -bin_width;
						// Now apply bias corrections for `sum_plnp`
						switch (est) {
						case EntropyEstimators::MM:
							// Miller-Madow estimator corrector
							sum_plnp += ((nbins - blank_bins) - 1)
									/ (2.0 * n_data);
							break;
						default:
							break;
						}
					}
					entropy[estm_offset + ent_off + i] = sum_plnp;
				}
				frq_off += nbins;
			}
		}
	}
	return 0;
}

int entropy2D(const BAT_t type_d1, const BAT_t type_d2, const u_int id1,
		const u_int id2, const std::vector<EntropyEstimators> &ests,
		const hbin_t nsteps, const ull_int step_size, const double min_d1,
		const double max_d1, const double min_d2, const double max_d2,
		const bool jacobian, const bool kde,
		const std::vector<hbin_t> bin_schemes,
		const std::vector<ull_int> &freqs_obs, std::vector<double> &entropy) {
	hbin_t nestm = ests.size();
	size_t nbins_sum_sq = 0;
	for (size_t i = 0; i < bin_schemes.size(); ++i) {
		nbins_sum_sq += bin_schemes[i] * bin_schemes[i];
	}
	if (nsteps * nbins_sum_sq != freqs_obs.size()) {
		return -1;
	} else if (nestm * nsteps * bin_schemes.size() != entropy.size()) {
		return -2;
	} else {
		for (hbin_t estid = 0; estid < nestm; ++estid) {
			EntropyEstimators est = ests[estid];
			u_int estm_offset = estid * nsteps * bin_schemes.size();

			u_int frq_off = 0;
			for (size_t tmp = 0, ent_off = 0; tmp < bin_schemes.size();
					ent_off += nsteps, ++tmp) {
				hbin_t nbins = bin_schemes[tmp];
				u_int nbins_sq = nbins * nbins;
				std::vector<std::vector<double>> freqs(nbins,
						std::vector<double>(nbins, 0.0)); // Normalized frequency
				double bin_width1 = (max_d1 - min_d1) / nbins;
				double bin_width_half1 = bin_width1 / 2;
				double bin_width2 = (max_d2 - min_d2) / nbins;
				double bin_width_half2 = bin_width2 / 2;

				if (kde) {
					bin_width1 = 1.0;
					bin_width2 = 1.0;
				}

				std::vector<double> bin_mid1(nbins, 0.0);
				std::vector<double> bin_mid2(nbins, 0.0);
				double sum_plnp = 0.0;
				double probdens = 0.0;
				double jacobian_factor1 = 1.0;
				double jacobian_factor2 = 1.0;
				bool has_jacob_fact1 = false;
				bool has_jacob_fact2 = false;
				u_int blank_bins = 0;
				u_int singleton_bins = 0;
				if (type_d1 == BAT_t::BOND || type_d1 == BAT_t::ANGLE) {
					has_jacob_fact1 = true;
					for (size_t j = 0; j < nbins; ++j) {
						bin_mid1[j] = min_d1 + (j * bin_width1)
								+ bin_width_half1;
					}
				}

				if (type_d2 == BAT_t::BOND || type_d2 == BAT_t::ANGLE) {
					has_jacob_fact2 = true;
					for (size_t j = 0; j < nbins; ++j) {
						bin_mid2[j] = min_d2 + (j * bin_width2)
								+ bin_width_half2;
					}
				}
				double bin_area = bin_width1 * bin_width2;
				for (size_t i = 0; i < nsteps; ++i) {
					ull_int offset = (i * nbins_sum_sq) + frq_off;
					double n_data = (i + 1.0) * step_size;
					double n_data_normlized = 1.0;
					n_data_normlized =
							(est == EntropyEstimators::CS) ? n_data : (n_data);

					// Normalize observed frequency,
					for (size_t j = 0; j < nbins; ++j) {
						u_int off_j = j * nbins;
						for (size_t k = 0; k < nbins; ++k) {
							freqs[j][k] = freqs_obs[offset + off_j + k]
									/ n_data_normlized;
						}
					}
					sum_plnp = 0.0;
					blank_bins = 0;
					singleton_bins = 0;

					// Apply frequency correction if any applicable
					switch (est) {
					case EntropyEstimators::ML:
						break;
					case EntropyEstimators::MM:
						for (size_t j = 0; j < nbins; ++j) {
							u_int off_j = j * nbins;
							for (size_t k = 0; k < nbins; ++k) {
								if (freqs_obs[offset + off_j + k] == 0) {
									++blank_bins;
								}
							}
						}
						break;
					case EntropyEstimators::CS:
						// Count: Singleton bins in distribution
						for (size_t j = 0; j < nbins; ++j) {
							u_int off_j = j * nbins;
							for (size_t k = 0; k < nbins; ++k) {
								if (freqs_obs[offset + off_j + k] == 1) {
									++singleton_bins;
								}
							}
						}
						// FREQUENCY_CORRECTION:: Good-Turing correction for Chao-Shen,
						for (size_t j = 0; j < nbins; ++j) {
							u_int off_j = j * nbins;
							for (size_t k = 0; k < nbins; ++k) {
								freqs[j][k] = (1 - (singleton_bins / n_data))
										* (freqs_obs[offset + off_j + k]
												/ n_data);
							}
						}
						break;
					case EntropyEstimators::JS:
						double lambda_hat = 0.0;
						double tk = 1.0 / (double) nbins;
						double sum_theta_ML2 = 0;
						double tk_sum_theta_ML2 = 0;

						for (size_t j = 0; j < nbins; ++j) {
							for (size_t k = 0; k < nbins; ++k) {
								sum_theta_ML2 += freqs[j][k] * freqs[j][k];
								tk_sum_theta_ML2 += ((tk - freqs[j][k])
										* (tk - freqs[j][k]));
							}
						}
						lambda_hat = (1 - sum_theta_ML2)
								/ ((n_data - 1) * tk_sum_theta_ML2);

						for (size_t j = 0; j < nbins; ++j) {
							for (size_t k = 0; k < nbins; ++k) {
								freqs[j][k] = (lambda_hat * tk)
										+ ((1 - lambda_hat) * freqs[j][k]);
							}
						}
						break;
					}

					// Estimate associated entropy from normalized frequencies
					if (est == EntropyEstimators::CS) {
						for (size_t j = 0; j < nbins; ++j) {
							if (has_jacob_fact1) {
								switch (type_d1) {
								case BAT_t::BOND:
									jacobian_factor1 = pow(bin_mid1[j], 2);
									break;
								case BAT_t::ANGLE:
									jacobian_factor1 = sin(bin_mid1[j]);
									break;
								default: // When it is a dihedral, associated jacobian is 1.0
									break;
								}
							}
							double plnp_tmp = 0.0;
							for (size_t k = 0; k < nbins; ++k) {
								probdens = freqs[j][k] / bin_area;
								if (probdens > 0.0) {
									// HT ESTIMATOR: Horowitz-Thompson estimators, denomirator
									double denom = 1
											- pow((1 - freqs[j][k]), n_data);
									if (has_jacob_fact2) {
										switch (type_d2) {
										case BAT_t::BOND:
											jacobian_factor2 = pow(bin_mid2[j],
													2);
											break;
										case BAT_t::ANGLE:
											jacobian_factor2 = sin(bin_mid2[j]);
											break;
										default: // When it is a dihedral, associated jacobian is 1.0
											break;
										}
									}
									if (jacobian) {
										double jacobian_factor =
												(jacobian_factor1
														* jacobian_factor2);
										probdens = probdens / jacobian_factor;
										plnp_tmp += (jacobian_factor * probdens
												* log(probdens) / denom);
									} else {
										plnp_tmp += (probdens * log(probdens)
												/ denom);
									}
								}
							}
							plnp_tmp *= -bin_area;
							sum_plnp += plnp_tmp;
						}

					} else {
						for (size_t j = 0; j < nbins; ++j) {
							if (has_jacob_fact1) {
								switch (type_d1) {
								case BAT_t::BOND:
									jacobian_factor1 = pow(bin_mid1[j], 2);
									break;
								case BAT_t::ANGLE:
									jacobian_factor1 = sin(bin_mid1[j]);
									break;
								default: // When it is a dihedral, associated jacobian is 1.0
									break;
								}
							}
							double plnp_tmp = 0.0;
							for (size_t k = 0; k < nbins; ++k) {
								probdens = freqs[j][k] / bin_area;
								if (probdens > 0.0) {
									if (has_jacob_fact2) {
										switch (type_d2) {
										case BAT_t::BOND:
											jacobian_factor2 = pow(bin_mid2[k],
													2);
											break;
										case BAT_t::ANGLE:
											jacobian_factor2 = sin(bin_mid2[k]);
											break;
										default: // When it is a dihedral, associated jacobian is 1.0
											break;
										}
									}
									if (jacobian) {
										double jacobian_factor =
												(jacobian_factor1
														* jacobian_factor2);
										probdens = probdens / jacobian_factor;
										plnp_tmp += (jacobian_factor * probdens
												* log(probdens));

									} else {
										plnp_tmp += (probdens * log(probdens));
									}
								}
							}
							plnp_tmp *= -bin_area;
							sum_plnp += plnp_tmp;
						}

						// Now apply bias corrections for `sum_plnp`
						switch (est) {
						case EntropyEstimators::MM:
							// Miller-Madow estimator corrector
							sum_plnp += ((nbins_sq - blank_bins) - 1)
									/ (2.0 * n_data);
							break;
						default:
							break;
						}
					}
					entropy[estm_offset + ent_off + i] = sum_plnp;
				}
				frq_off += nbins_sq;
			}
		}
	}
	return 0;
}

int entropyND(const hbin_t nD, const EntropyEstimators est, const hbin_t nsteps,
		const ull_int step_size, const std::vector<double> min_d,
		const std::vector<double> max_d, const std::vector<hbin_t> bin_schemes,
		const std::vector<ull_int> &freqs_obs, std::vector<double> *entropy) {

	size_t nbins_powN_sum = 0;
	for (size_t i = 0; i < bin_schemes.size(); ++i) {
		nbins_powN_sum += pow(bin_schemes[i], nD);
	}
	if (nsteps * nbins_powN_sum != freqs_obs.size()) {
		return -1;
	} else if (nsteps * bin_schemes.size() != (*entropy).size()) {
		return -2;
	} else {
		u_int frq_off = 0;
		for (size_t tmp = 0, ent_off = 0; tmp < bin_schemes.size(); ent_off +=
				nsteps, ++tmp) {
			hbin_t nbins = bin_schemes[tmp];
			u_int nbins_powN = pow(nbins, nD);
			std::vector<double> freqs(nbins_powN_sum, 0.0);
			std::vector<double> bin_width(nD, 0.0);

			double bin_hypvolm = 1.0;
			for (int i = 0; i < nD; ++i) {
				bin_width[i] = (max_d[i] - min_d[i]) / nbins;
				bin_hypvolm *= bin_width[i];
			}
			double sum_plnp = 0.0;
			double probdens = 0.0;
			u_int blank_bins = 0;
			u_int singleton_bins = 0;

			for (size_t i = 0; i < nsteps; ++i) {
				sum_plnp = 0.0;
				ull_int offset = (i * nbins_powN_sum) + frq_off;
				double n_data = 1.0 * (i + 1) * step_size;
				for (size_t j = 0; j < nbins_powN_sum; ++j) {
					freqs[j] = freqs_obs[offset + j] / (n_data * bin_hypvolm);
				}
				sum_plnp = 0.0;
				blank_bins = 0;
				singleton_bins = 0;
				switch (est) {
				case EntropyEstimators::ML:
					break;
				case EntropyEstimators::MM:
					for (size_t j = 0; j < nbins_powN_sum; ++j) {
						if (freqs_obs[offset + j] == 0) {
							++blank_bins;
						}
					}
					break;
				case EntropyEstimators::CS:
					for (size_t j = 0; j < nbins_powN_sum; ++j) {
						if (freqs_obs[offset + j] == 1) {
							++singleton_bins;
						}
					}
					// FREQUENCY_CORRECTION:: for Chao-Shen, Good-Turing correction
					for (size_t j = 0; j < nbins_powN_sum; ++j) {
						freqs[j] = (1 - (singleton_bins * freqs[j] / n_data));
					}
					break;
				case EntropyEstimators::JS:
					double lambda_hat = 0.0;
					double tk = 1.0 / (double) nbins;
					double sum_theta_ML2 = 0;
					double tk_sum_theta_ML2 = 0;

					for (size_t j = 0; j < nbins_powN_sum; ++j) {
						sum_theta_ML2 += freqs[j] * freqs[j];
						tk_sum_theta_ML2 += ((tk - freqs[j]) * (tk - freqs[j]));
					}
					lambda_hat = (1 - sum_theta_ML2)
							/ ((n_data - 1) * tk_sum_theta_ML2);

					for (size_t j = 0; j < nbins_powN_sum; ++j) {
						freqs[j] = (lambda_hat * tk)
								+ ((1 - lambda_hat) * freqs[j]);
					}
					break;
				}
				if (est == EntropyEstimators::CS) {
					double plnp_tmp = 0.0;
					for (size_t j = 0; j < nbins_powN_sum; ++j) {
						probdens = freqs[j];
						if (probdens > 0.0) {
							// HT ESTIMATOR: Horowitz-Thompson estimators, denomirator
							double denom = 1 - pow((1 - freqs[j]), n_data);
							plnp_tmp += (probdens * log(probdens) / denom);
						}
					}
					plnp_tmp *= -bin_hypvolm;
					sum_plnp += plnp_tmp;
				} else {
					double plnp_tmp = 0.0;
					for (size_t j = 0; j < nbins_powN_sum; ++j) {
						probdens = freqs[j];
						if (probdens > 0.0) {
							plnp_tmp += (probdens * log(probdens));
						}
					}
					plnp_tmp *= -bin_hypvolm;
					sum_plnp += plnp_tmp;

					// Now apply bias corrections for `sum_plnp`
					switch (est) {
					case EntropyEstimators::MM:
						// Miller-Madow estimator corrector
						sum_plnp += ((nbins_powN - blank_bins) - 1)
								/ (2.0 * n_data);
						break;
					default:
						break;
					}
				}
				(*entropy)[ent_off + i] = sum_plnp;
			}
			frq_off += nbins_powN;
		}
	}
	return 0;
}

