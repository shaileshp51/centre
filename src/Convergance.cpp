/*
 * Convergance.cpp
 *
 *  Created on: 21-Oct-2016
 *      Author: shailesh
 */

#include "Convergence.h"
#include "Centre.h"

void sumContrib1D(const u_int dim_len, const std::vector<double> &dataIn,
		double &entSum) {
	ull_int n_data = dataIn.size();

	double sum = 0.0;
#ifdef USE_OMPMPI
#pragma omp parallel for reduction(+:sum)
#endif
	for (u_int n = 0; n < n_data; ++n)
		sum += dataIn[n];
	entSum = sum;
}

void sumContrib2D(const u_int start_indx, const u_int count,
		const std::vector<double> &dataIn, double &entSum) {
	ull_int n_data = dataIn.size();
	u_int n_end = start_indx + count;
	if (n_end <= n_data) {
		double sum = 0.0;
#ifdef USE_OMPMPI
#pragma omp parallel for reduction(+:sum)
#endif
		for (u_int n = start_indx; n < n_end; ++n)
			sum += dataIn[n];
		entSum = sum;
	}
}

void sumContrib2D(const hbin_t n_schemes, const u_int dim_len,
		const std::vector<double> &dataIn, std::vector<double> &entSum) {
	ull_int n_data = dataIn.size();
	for (hbin_t s = 0; s < n_schemes; ++s) {
		double sum = 0.0;
#ifdef USE_OMPMPI
#pragma omp parallel for reduction(+:sum)
#endif
		for (u_int n = s; n < n_data; n += n_schemes)
			sum += dataIn[n];
		entSum[s] = sum;
	}
}

void getMaxMI2D(const SubsetData &subset, const Neighbors &neigh,
		const char d1type, const char d2type, const std::vector<double> S1,
		const std::vector<int> &block, const std::vector<double> &dataIn,
		const u_int dataOffset, std::vector<double> &dataMI) {

	u_int dataIdx = dataOffset;
	u_int maxMIsIdx = 0;
	std::map<u_int, u_int> id2indx;
	std::vector<u_int> sub_dims;
	std::vector<u_int> d1_vec;
	switch (d1type) {
	case 'B':
		neigh.bondKeys(d1_vec);
		sub_dims = subset.getBonds();
		break;
	case 'A':
		neigh.angleKeys(d1_vec);
		sub_dims = subset.getAngles();
		break;
	case 'T':
		neigh.torsionKeys(d1_vec);
		sub_dims = subset.getTorsions();
		break;
	}
	u_int dim1_len = d1_vec.size();
	for (u_int i = 0; i < sub_dims.size(); ++i) {
		id2indx[sub_dims[i]] = i;
	}
	for (u_int b1id = block[0]; b1id <= block[1]; ++b1id) {
		double mi_max = -1.0e200;
		u_int tmp_d1 = d1_vec[b1id];
		std::vector<u_int> b_neigh;
		switch (d2type) {
		case 'B':
			b_neigh = neigh.getBond(d1_vec[b1id]);
			break;
		case 'A':
			b_neigh = neigh.getAngle(d1_vec[b1id]);
			break;
		case 'T':
			b_neigh = neigh.getTorsion(d1_vec[b1id]);
			break;
		}

		u_int n_b1neigh = b_neigh.size();
		u_int b2_s = (block[0] == b1id) ? block[2] : 0;
		u_int b2_e = (block[1] == b1id) ? block[3] : n_b1neigh;
		for (u_int b2id = b2_s; b2id < b2_e; ++b2id) {
			u_int tmp_d2 = b_neigh[b2id];
			dataMI[dataIdx] = S1[id2indx[tmp_d1]] + S1[id2indx[tmp_d2]]
					- dataIn[dataIdx];
			++dataIdx;
		}
	}
}

void getMaxMI2D(const SubsetData &subset, const Neighbors &neigh,
		const char d1type, const char d2type, const std::vector<double> S1,
		const std::vector<double> &dataIn, std::vector<double> &dataMI) {

	u_int dataIdx = 0;
	u_int maxMIsIdx = 0;
	std::map<u_int, u_int> id2indx;
	std::vector<u_int> sub_dims;
	std::vector<u_int> d1_vec;
	switch (d1type) {
	case 'B':
		neigh.bondKeys(d1_vec);
		sub_dims = subset.getBonds();
		break;
	case 'A':
		neigh.angleKeys(d1_vec);
		sub_dims = subset.getAngles();
		break;
	case 'T':
		neigh.torsionKeys(d1_vec);
		sub_dims = subset.getTorsions();
		break;
	}
	u_int dim1_len = d1_vec.size();
	for (u_int i = 0; i < sub_dims.size(); ++i) {
		id2indx[sub_dims[i]] = i;
	}
	for (u_int b1id = 0; b1id < sub_dims.size(); ++b1id) {
		double mi_max = -1.0e200;
		u_int tmp_d1 = d1_vec[b1id];
		std::vector<u_int> b_neigh;
		switch (d2type) {
		case 'B':
			b_neigh = neigh.getBond(d1_vec[b1id]);
			break;
		case 'A':
			b_neigh = neigh.getAngle(d1_vec[b1id]);
			break;
		case 'T':
			b_neigh = neigh.getTorsion(d1_vec[b1id]);
			break;
		}

		u_int n_b1neigh = b_neigh.size();
		u_int b2_s = 0;
		u_int b2_e = n_b1neigh;
		for (u_int b2id = b2_s; b2id < b2_e; ++b2id) {
			u_int tmp_d2 = b_neigh[b2id];
			dataMI[dataIdx] = S1[id2indx[tmp_d1]] + S1[id2indx[tmp_d2]]
					- dataIn[dataIdx];
			++dataIdx;
		}
	}
}

void getMaxMI2D(const std::map<u_int, u_int> &id2index_1d1,
		const std::map<u_int, u_int> &id2index_1d2,
		const std::vector<double> &S1, const std::vector<double> &S2,
		std::map<u_int, std::pair<u_int, u_int>> &index2id_2d,
		const std::vector<double> &dataIn, std::vector<double> &dataMI) {

	/*for (u_int idx1 = 0; idx1 < index2id_2d.size(); ++idx1) {
	 auto p = index2id_2d[idx1];

	 dataMI[idx1] = S1[id2index_1d1.at(p.first)]
	 + S2[id2index_1d2.at(p.second)] - dataIn[idx1];

	 }
	 */
	for (auto const &pp : index2id_2d) {
		dataMI[pp.first] = S1[id2index_1d1.at(pp.second.first)]
				+ S2[id2index_1d2.at(pp.second.second)] - dataIn[pp.first];
	}
}

void getMaxCrossMI2D(const SubsetData &subset, const Neighbors &neigh,
		const char d1type, const char d2type, const std::vector<double> S1d1,
		const std::vector<double> S1d2, const std::vector<int> &block,
		const std::vector<double> &dataIn, const u_int dataOffset,
		std::vector<double> &mi_data) {

	u_int dataIdx = dataOffset;
	u_int maxMIsIdx = 0;
	std::map<u_int, u_int> id2indx1;
	std::map<u_int, u_int> id2indx2;

	std::vector<u_int> d1_vec;
	std::vector<u_int> d2_vec;

	switch (d1type) {
	case 'B':
		d1_vec = subset.getBonds();
		break;
	case 'A':
		d1_vec = subset.getAngles();
		break;
	case 'T':
		d1_vec = subset.getTorsions();
		break;
	}
	switch (d2type) {
	case 'B':
		d2_vec = subset.getBonds();
		break;
	case 'A':
		d2_vec = subset.getAngles();
		break;
	case 'T':
		d2_vec = subset.getTorsions();
		break;
	}

	u_int dim1_len = d1_vec.size();
	u_int dim2_len = d2_vec.size();
	for (u_int i = 0; i < dim1_len; ++i) {
		id2indx1[d1_vec[i]] = i;
	}
	for (u_int i = 0; i < dim2_len; ++i) {
		id2indx2[d2_vec[i]] = i;
	}
	std::vector<u_int> dim1_Keys;
	if (d1type == 'B' && d2type == 'B') {
		neigh.bondKeys(dim1_Keys);
	} else if (d1type == 'B' && d2type == 'A') {
		neigh.bacrossKeys(dim1_Keys);
	} else if (d1type == 'B' && d2type == 'T') {
		neigh.bdcrossKeys(dim1_Keys);
	} else if (d1type == 'A' && d2type == 'A') {
		neigh.angleKeys(dim1_Keys);
	} else if (d1type == 'A' && d2type == 'T') {
		neigh.adcrossKeys(dim1_Keys);
	} else if (d1type == 'T' && d2type == 'T') {
		neigh.torsionKeys(dim1_Keys);
	}
	for (u_int b1id = block[0]; b1id <= block[1]; ++b1id) {
		double mi_max = -1.0e200;
		u_int tmp_d1 = dim1_Keys[b1id];
		std::vector<u_int> b_neigh;
		if (d1type == 'B' && d2type == 'B') {
			b_neigh = neigh.getBond(tmp_d1);
		} else if (d1type == 'B' && d2type == 'A') {
			b_neigh = neigh.getBacross(tmp_d1);
		} else if (d1type == 'B' && d2type == 'T') {
			b_neigh = neigh.getBdcross(tmp_d1);
		} else if (d1type == 'A' && d2type == 'A') {
			b_neigh = neigh.getAngle(tmp_d1);
		} else if (d1type == 'A' && d2type == 'T') {
			b_neigh = neigh.getAdcross(tmp_d1);
		} else if (d1type == 'T' && d2type == 'T') {
			b_neigh = neigh.getTorsion(tmp_d1);
		}

		u_int b2_s = (block[0] == b1id) ? block[2] : 0;
		u_int b2_e = (block[1] == b1id) ? block[3] : b_neigh.size();
		for (u_int b2id = b2_s; b2id < b2_e; ++b2id) {
			u_int tmp_d2 = b_neigh[b2id];
			mi_data[dataIdx] = S1d1[id2indx1[tmp_d1]] + S1d2[id2indx2[tmp_d2]]
					- dataIn[dataIdx];
			++dataIdx;
		}
	}
}

void getMaxCrossMI2D(const SubsetData &subset, const Neighbors &neigh,
		const char d1type, const char d2type, const std::vector<double> S1d1,
		const std::vector<double> S1d2, const std::vector<double> &dataIn,
		std::vector<double> &mi_data) {

	u_int dataIdx = 0;
	u_int maxMIsIdx = 0;
	std::map<u_int, u_int> id2indx1;
	std::map<u_int, u_int> id2indx2;

	std::vector<u_int> d1_vec;
	std::vector<u_int> d2_vec;

	switch (d1type) {
	case 'B':
		d1_vec = subset.getBonds();
		break;
	case 'A':
		d1_vec = subset.getAngles();
		break;
	case 'T':
		d1_vec = subset.getTorsions();
		break;
	}
	switch (d2type) {
	case 'B':
		d2_vec = subset.getBonds();
		break;
	case 'A':
		d2_vec = subset.getAngles();
		break;
	case 'T':
		d2_vec = subset.getTorsions();
		break;
	}

	u_int dim1_len = d1_vec.size();
	u_int dim2_len = d2_vec.size();
	for (u_int i = 0; i < dim1_len; ++i) {
		id2indx1[d1_vec[i]] = i;
	}
	for (u_int i = 0; i < dim2_len; ++i) {
		id2indx2[d2_vec[i]] = i;
	}
	std::vector<u_int> dim1_Keys;
	if (d1type == 'B' && d2type == 'B') {
		neigh.bondKeys(dim1_Keys);
	} else if (d1type == 'B' && d2type == 'A') {
		neigh.bacrossKeys(dim1_Keys);
	} else if (d1type == 'B' && d2type == 'T') {
		neigh.bdcrossKeys(dim1_Keys);
	} else if (d1type == 'A' && d2type == 'A') {
		neigh.angleKeys(dim1_Keys);
	} else if (d1type == 'A' && d2type == 'T') {
		neigh.adcrossKeys(dim1_Keys);
	} else if (d1type == 'T' && d2type == 'T') {
		neigh.torsionKeys(dim1_Keys);
	}
	for (u_int b1id = 0; b1id < dim1_Keys.size(); ++b1id) {
		double mi_max = -1.0e200;
		u_int tmp_d1 = dim1_Keys[b1id];
		std::vector<u_int> b_neigh;
		if (d1type == 'B' && d2type == 'B') {
			b_neigh = neigh.getBond(tmp_d1);
		} else if (d1type == 'B' && d2type == 'A') {
			b_neigh = neigh.getBacross(tmp_d1);
		} else if (d1type == 'B' && d2type == 'T') {
			b_neigh = neigh.getBdcross(tmp_d1);
		} else if (d1type == 'A' && d2type == 'A') {
			b_neigh = neigh.getAngle(tmp_d1);
		} else if (d1type == 'A' && d2type == 'T') {
			b_neigh = neigh.getAdcross(tmp_d1);
		} else if (d1type == 'T' && d2type == 'T') {
			b_neigh = neigh.getTorsion(tmp_d1);
		}

		u_int b2_s = 0;
		u_int b2_e = b_neigh.size();
		for (u_int b2id = b2_s; b2id < b2_e; ++b2id) {
			u_int tmp_d2 = b_neigh[b2id];
			mi_data[dataIdx] = S1d1[id2indx1[tmp_d1]] + S1d2[id2indx2[tmp_d2]]
					- dataIn[dataIdx];
			++dataIdx;
		}
	}
}

