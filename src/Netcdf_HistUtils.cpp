#include <cstdlib>

#include <iostream>
#include <cmath>

#include "Netcdf_HistUtils.h"
#include "NetcdfFile.h"

BinGroup::BinGroup(const u_int id, const hbin_t nstep_eff,
		const ull_int bin_schemes_sum) {
	record_id.resize(1);
	record_id[0] = id;
	nstep_eff_ = nstep_eff;
	bin_dims_sum_mid_ = bin_schemes_sum;
	bin_pow_dims_sum_ = bin_schemes_sum;
	bin_mid.resize(bin_dims_sum_mid_);
	bin_mid.assign(bin_dims_sum_mid_, 0.0);
	ull_int frq_sz = nstep_eff * bin_pow_dims_sum_;
	dim_extremes.resize(record_id.size() * 4);
	dim_extremes.assign(record_id.size() * 4, 0.0);
	count.resize(frq_sz);
	count.assign(frq_sz, 0);
}

BinGroup::BinGroup(const std::vector<u_int> &id, const hbin_t nstep_eff,
		const ull_int bin_mids_dim, const ull_int bin_pow_dim_sum) {
	record_id.assign(id.begin(), id.end());
	nstep_eff_ = nstep_eff;
	bin_dims_sum_mid_ = bin_mids_dim;
	bin_pow_dims_sum_ = bin_pow_dim_sum;
	bin_mid.assign(bin_dims_sum_mid_, 0.0);
	ull_int frq_sz = nstep_eff * bin_pow_dims_sum_;
	dim_extremes.resize(record_id.size() * 4);
	dim_extremes.assign(record_id.size() * 4, 0.0);
	count.resize(frq_sz);
	count.assign(frq_sz, 0);
}

void BinGroup::getId(ull_int *size, u_int **p) {
	*size = (ull_int) record_id.size();
	*p = record_id.data();
}

void BinGroup::getExtremes(ull_int *size, double **p) {
	*size = (ull_int) dim_extremes.size();
	*p = dim_extremes.data();
}

void BinGroup::getBinMids(ull_int *sz, double **data) {
	*sz = (ull_int) bin_mid.size();
	*data = bin_mid.data();
}

void BinGroup::getBinFreqs(ull_int *nsteps, ull_int *nbins, ull_int **data) {
	*nsteps = nstep_eff_;
	*nbins = (ull_int) count.size() / nstep_eff_;
	*data = count.data();
}

void BinGroup::setId(const ull_int start_index_in, const ull_int in_count,
		const std::vector<u_int> &in) {
	std::vector<u_int>::const_iterator it;
	it = in.begin() + start_index_in;
	record_id.assign(it, it + in_count);
}

void BinGroup::setId(const u_int &in) {
	record_id[0] = in;
}

void BinGroup::setExtremes(const std::vector<double> &in) {
	std::vector<double>::const_iterator it;
	it = in.begin();
	dim_extremes.assign(it, it + in.size());
}

int BinGroup::setBinMids(const ull_int start_bin_id, const ull_int bins_count,
		const std::vector<double> &in) {
	ull_int end_ = start_bin_id + bins_count;
	for (ull_int i = start_bin_id, j = 0; i < end_; ++i, ++j) {
		bin_mid[i] = in[j];
	}
	return 0;
}

int BinGroup::setBinMids(const ull_int start_bin_id, const ull_int bins_count,
		ull_int start_in, const std::vector<double> &in) {
	ull_int end_ = start_bin_id + bins_count;
	for (ull_int i = start_bin_id, j = start_in; i < end_; ++i, ++j) {
		bin_mid[i] = in[j];
	}
	return 0;
}

int BinGroup::setBinFreqs(const hbin_t step_id, const ull_int start_bin_id,
		const ull_int bin_count, const std::vector<ull_int> &in) {
	ull_int end_ = start_bin_id + bin_count;
	for (ull_int i = start_bin_id, j = 0; i < end_; ++i, ++j) {
		count[((int) step_id * bin_pow_dims_sum_) + i] = in[j];
	}
	return 0;
}

int BinGroup::setBinFreqs(const hbin_t start_step_id, const hbin_t steps_count,
		const ull_int start_bin_id, const ull_int bins_count,
		const std::vector<ull_int> &in) {
	ull_int steps_end_ = (int) start_step_id + (int) steps_count;
	ull_int bins_end_ = start_bin_id + bins_count;
	ull_int b = 0;
	for (ull_int a = start_step_id; a < steps_end_; ++a) {
		for (ull_int i = start_bin_id; i < bins_end_; ++i) {
			count[(a * bin_pow_dims_sum_) + i] = in[b];
			++b;
		}
	}
	return 0;
}

int BinGroup::setBinFreqs(const hbin_t start_step_in,
		const hbin_t step_stride_in, const hbin_t steps_cnt_eff_in,
		const ull_int start_bin_id, const ull_int bins_count,
		const std::vector<ull_int> &in) {
	ull_int steps_end_ = (int) start_step_in
			+ (int) step_stride_in * steps_cnt_eff_in;
	ull_int bins_end_ = start_bin_id + bins_count;
	ull_int b = 0;
	ull_int offset_in = 0;
	ull_int offset_bgrp = 0;
	for (ull_int a = start_step_in; a < steps_end_; a += step_stride_in) {
		offset_in = a * bin_pow_dims_sum_;
		offset_bgrp = b * bin_pow_dims_sum_;
		for (ull_int i = start_bin_id; i < bins_end_; ++i) {
			count[offset_bgrp + i] = in[offset_in + i];
		}
		++b;
	}
	return 0;
}

int BinGroup::setBinFreqs(const hbin_t start_step_id, const hbin_t steps_count,
		const ull_int start_bin_id, const ull_int bins_count,
		const ull_int start_in, const ull_int *in) {
	ull_int steps_end_ = (int) start_step_id + (int) steps_count;
	ull_int bins_end_ = start_bin_id + bins_count;
	ull_int b = start_in;
	for (ull_int a = start_step_id; a < steps_end_; ++a) {
		ull_int offset = (a * bin_pow_dims_sum_);
		for (ull_int i = start_bin_id; i < bins_end_; ++i) {
			count[offset + i] = *(in + b);
			++b;
		}
	}
	return 0;
}

int BinGroup::setBinFreqs(const hbin_t start_step_id, const hbin_t steps_count,
		const ull_int start_bin_id, const ull_int bins_count,
		const ull_int start_in, const std::vector<ull_int> &in) {
	ull_int steps_end_ = (int) start_step_id + (int) steps_count;
	ull_int bins_end_ = start_bin_id + bins_count;
	ull_int b = start_in;
	for (ull_int a = start_step_id; a < steps_end_; ++a) {
		ull_int offset = (a * bin_pow_dims_sum_);
		for (ull_int i = start_bin_id; i < bins_end_; ++i) {
			count[offset + i] = in[b];
			++b;
		}
	}
	return 0;
}

int BinGroup::setBinFreqs(hbin_t step_id, ull_int bin_id, ull_int value) {
	count[(step_id * bin_pow_dims_sum_) + bin_id] = value;
	return 0;
}

std::ostream& operator<<(std::ostream &os, const BinGroup &bg) {
	os
			<< "++++++++++++++++++++++++++++++++++++ BinGroup Object +++++++++++++++++++++++++++++++++++"
			<< std::endl;
	os << "Record Id: {";
	for (u_int i = 0; i < bg.record_id.size(); ++i) {
		os << (u_int) bg.record_id[i];
		if (i + 1 != bg.record_id.size())
			os << ", ";
	}
	os << "}" << std::endl;
	os << "N-Step Effective: " << (int) bg.nstep_eff_ << std::endl;
	os << "Bin mids dim: " << (int) bg.bin_dims_sum_mid_ << std::endl;
	os << "Bin Mids: ";
	for (u_int i = 0; i < bg.bin_mid.size(); ++i) {
		os << bg.bin_mid[i];
		if (i + 1 != bg.bin_mid.size())
			os << ", ";
	}
	os << std::endl;
	os << "Bin Power Dim Sum: " << (int) bg.bin_pow_dims_sum_ << std::endl;
	os << "Bin Frequences (nstep x bin_pow_dims_sum): {"
			<< (u_int) bg.nstep_eff_ << ", " << (u_int) bg.bin_pow_dims_sum_
			<< "}" << std::endl;
	for (u_int i = 0; i < bg.nstep_eff_; ++i) {
		std::cout << "Step [" << (i + 1) << "]: ";
		for (u_int j = 0; j < bg.bin_pow_dims_sum_; ++j) {
			ull_int offset = (i * bg.bin_pow_dims_sum_);
			os << (ull_int) bg.count[offset + j];
			if ((j + 1) % bg.bin_pow_dims_sum_ != 0)
				os << ", ";
		}
		std::cout << std::endl;
	}
	os
			<< "====================================== BinGroup Object End ======================================="
			<< std::endl;
	return os;
}

Netcdf_HistUtil::Netcdf_HistUtil(const std::string &filename,
		const std::string &c, const hbin_t ndims, const hbin_t order,
		const hbin_t first_step, const hbin_t nsteps, const hbin_t step_stride,
		const std::vector<ull_int> &dim_lens, const TensorType &typ,
		const hbin_t bin_schemes_count, const std::vector<hbin_t> &bin_schemes) :
		NetcdfFile(filename) {
	filename_ = std::string(filename.c_str());
	contents_ = std::string(c.c_str());
	ndims_ = (ndims < 1) ? 1 : ndims;
	order_ = (order < 1) ? 1 : order;
	nsteps_ = (nsteps < 1) ? 1 : nsteps;
	first_step_ = (first_step < 0) ? nsteps_ - 1 :
					(first_step >= nsteps_) ? nsteps_ - 1 : first_step;
	step_stride_ = (step_stride < 1) ? 1 :
					(step_stride > nsteps_) ? nsteps_ : step_stride;
	nstep_eff_ = (nsteps_ - first_step_ + step_stride_ - 1) / step_stride_;
	bin_schemes_count_ = bin_schemes_count;
	storage_type_ = typ;
	if (storage_type_ != TensorType::FULL) {
		ndiag_ = 1;
	} else {
		ndiag_ = 0;
	}

	nrecords = 0;
	std::vector<ull_int>::const_iterator it;
	for (it = dim_lens.begin(); it != dim_lens.end(); ++it) {
		dim_lens_.push_back(*it);
	}
	u_int bin_schemes_sum_ = 0;
	bin_pow_dims_sum_ = 0;
	bin_schemes_.reserve(bin_schemes_count_);
	for (size_t i = 0; i < bin_schemes_count_; ++i) {
		bin_schemes_.push_back(bin_schemes[i]);
		bin_schemes_sum_ += bin_schemes[i];
		bin_pow_dims_sum_ += pow(bin_schemes[i], ndims);
	}
	dimsXbin_schemes_sum_ = bin_schemes_sum_ * ndims;
	bin_dims_sum_mid_ = dimsXbin_schemes_sum_;
	nstepsEffDID_ = -1;
	recordIdDID_ = -1;
	binDimsDID_ = -1;
	binPowDimDID_ = -1;
	recordIdDID_ = -1;
	recordsDID_ = -1;
	extrmDID_ = -1;

	binFreqVID_ = -1;
	recordIdVID_ = -1;
	binMidsVID_ = -1;
	extrmVID_ = -1;
}

std::ostream& operator<<(std::ostream &os, const Netcdf_HistUtil &ht) {
	os
			<< "++++++++++++++++++++++++++++++++++++ HistUtil Object +++++++++++++++++++++++++++++++++++"
			<< std::endl;
	os << "Filename: " << ht.filename_ << std::endl;
	os << "Content: " << ht.contents_ << std::endl;
	os << "N-Dims: " << (int) ht.ndims_ << std::endl;
	os << "Order: " << (int) ht.order_ << std::endl;
	os << "N-Step: " << (int) ht.nsteps_ << std::endl;
	os << "Step Stride: " << (int) ht.step_stride_ << std::endl;
	os << "N-Step Effective: " << (int) ht.nstep_eff_ << std::endl;
	os << "Bin Schemes Count: " << (int) ht.bin_schemes_count_ << std::endl;
	os << "Bin Schemes: ";
	for (int i = 0; i < ht.bin_schemes_count_; ++i) {
		os << (int) ht.bin_schemes_[i];
		if (i + 1 != ht.bin_schemes_count_)
			os << ", ";
	}
	os << std::endl;
	os << "Bin Schemes Sum: " << (int) ht.dimsXbin_schemes_sum_ << std::endl;
	os << "Storage Type: " << ht.storage_type_ << std::endl;
	os << "N-Diag: " << (int) ht.ndiag_ << std::endl;
	os
			<< "====================================== HistUtil Object End ======================================="
			<< std::endl;
	return os;
}

int Netcdf_HistUtil::removeRecord(const BinGroup &rec) {
	return 0;
}

int Netcdf_HistUtil::removeRecords(const std::vector<BinGroup> &recs) {
	return 0;
}

void Netcdf_HistUtil::clearRecords() {
	vec_records.clear();
}

int Netcdf_HistUtil::setupRead() {
	if (ncmode_ != NC_NOWRITE && ncid_ != -1) {
		NC_close();
	}
	if (ncid_ == -1) {
		NC_openRead(filename_);

		int temp;
		recordsDID_ = getDimInfo(HU_D_NC_RECORD, &temp);
		nrecords = temp;
		recordIdDID_ = getDimInfo(HU_D_NC_RECORDID, &temp);
		order_ = temp;
		extrmDID_ = getDimInfo(HU_D_NC_DIMEXTRM, &temp);
		binDimsDID_ = getDimInfo(HU_D_NC_BINMIDS, &temp);
		bin_dims_sum_mid_ = temp;
		binPowDimDID_ = getDimInfo(HU_D_NC_BINFREQS, &temp);
		bin_pow_dims_sum_ = temp;
		nstepsEffDID_ = getDimInfo(HU_D_NC_NSTEPEFF, &temp);
		nstep_eff_ = temp;

		// Get coord info
		recordIdVID_ = -1;
		if (nc_inq_varid(ncid_, HU_V_NC_RECORDID, &recordIdVID_) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file has recordIds.\n");
		}
		extrmVID_ = -1;
		if (nc_inq_varid(ncid_, HU_V_NC_DIMEXTRM, &extrmVID_) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file has recordIds.\n");
		}
		binMidsVID_ = -1;
		if (nc_inq_varid(ncid_, HU_V_NC_BINMID, &binMidsVID_) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file has bin mids.\n");
		}
		binFreqVID_ = -1;
		if (nc_inq_varid(ncid_, HU_V_NC_BINFREQ, &binFreqVID_) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file has frequencies.\n");
		}
	}
	return 0;
}

int Netcdf_HistUtil::setupWrite() {
	if (ncmode_ != NC_WRITE && ncid_ != -1) {
		NC_close();
	}
	if (ncid_ == -1) {
		NC_openWrite(filename_);

		int temp;
		recordsDID_ = getDimInfo(HU_D_NC_RECORD, &temp);
		nrecords = temp;
		recordIdDID_ = getDimInfo(HU_D_NC_RECORDID, &temp);
		order_ = temp;
		extrmDID_ = getDimInfo(HU_D_NC_DIMEXTRM, &temp);
		binDimsDID_ = getDimInfo(HU_D_NC_BINMIDS, &temp);
		bin_dims_sum_mid_ = temp;
		binPowDimDID_ = getDimInfo(HU_D_NC_BINFREQS, &temp);
		bin_pow_dims_sum_ = temp;
		nstepsEffDID_ = getDimInfo(HU_D_NC_NSTEPEFF, &temp);
		nstep_eff_ = temp;

		// Get coord info
		recordIdVID_ = -1;
		if (nc_inq_varid(ncid_, HU_V_NC_RECORDID, &recordIdVID_) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file has recordIds.\n");
		}

		extrmVID_ = -1;
		if (nc_inq_varid(ncid_, HU_V_NC_DIMEXTRM, &extrmVID_) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file has recordIds.\n");
		}

		binMidsVID_ = -1;
		if (nc_inq_varid(ncid_, HU_V_NC_BINMID, &binMidsVID_) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file has bin mids.\n");
		}
		binFreqVID_ = -1;
		if (nc_inq_varid(ncid_, HU_V_NC_BINFREQ, &binFreqVID_) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file has frequencies.\n");
		}
	}
	return 0;
}

int Netcdf_HistUtil::readRecord(const size_t dim_index,
		std::vector<double> &extrm_v, std::vector<ull_int> &freqs_v) {
	int ret_code = 0;
	if (extrm_v.size() != order_ * 4) {
		ret_code = -2;
	} else if (freqs_v.size() != nstep_eff_ * bin_pow_dims_sum_) {
		ret_code = -4;
	}
	if (ret_code != 0) {
		return ret_code;
	}

	setupRead();

	size_t start[3], count[3];

	start[0] = dim_index;
	start[1] = 0;
	start[2] = 0;
	count[0] = 1;
	count[1] = order_;
	count[2] = 4;

	if (checkNCerr(
			nc_get_vara_double(ncid_, extrmVID_, start, count,
					extrm_v.data()))) {
		mprinterr("Error: Reading record extrema data.\n");
		return 1;
	}

	start[1] = 0;
	count[1] = nstep_eff_;

	start[2] = 0;
	count[2] = dimsXbin_schemes_sum_;

	if (checkNCerr(
			nc_get_vara_ulonglong(ncid_, binFreqVID_, start, count,
					freqs_v.data()))) {
		mprinterr("Error: Reading record frequency data.\n");
		return 1;
	}

	NC_close();
	return 0;
}

int Netcdf_HistUtil::readRecord(const size_t dim_index,
		std::vector<u_int> &ids_v, std::vector<double> &extrm_v,
		std::vector<ull_int> &freqs_v) {
	int ret_code = 0;
	if (ids_v.size() != order_) {
		ret_code = -1;
	} else if (extrm_v.size() != order_ * 4) {
		ret_code = -2;
	} else if (freqs_v.size() != nstep_eff_ * bin_pow_dims_sum_) {
		ret_code = -4;
	}
	if (ret_code != 0) {
		return ret_code;
	}

	setupRead();

	size_t start[3], count[3];

	start[0] = dim_index;
	start[1] = 0;
	count[0] = 1;
	count[1] = order_;

	if (checkNCerr(
			nc_get_vara_uint(ncid_, recordIdVID_, start, count,
					ids_v.data()))) {
		mprinterr("Error: Reading record_id data.\n");
		return 1;
	}

	start[0] = dim_index;
	start[1] = 0;
	start[2] = 0;
	count[0] = 1;
	count[1] = order_;
	count[2] = 4;

	if (checkNCerr(
			nc_get_vara_double(ncid_, extrmVID_, start, count,
					extrm_v.data()))) {
		mprinterr("Error: Reading record_id data.\n");
		return 1;
	}

	start[1] = 0;
	count[1] = nstep_eff_;

	start[2] = 0;
	count[2] = bin_pow_dims_sum_;

	if (checkNCerr(
			nc_get_vara_ulonglong(ncid_, binFreqVID_, start, count,
					freqs_v.data()))) {
		mprinterr("Error: Reading record frequency data.\n");
		return 1;
	}

	NC_close();
	return 0;
}

int Netcdf_HistUtil::readRecord(const size_t dim_index,
		std::vector<u_int> &ids_v, std::vector<double> &extrm_v,
		std::vector<double> &mids_v, std::vector<ull_int> &freqs_v) {
	int ret_code = 0;
	if (ids_v.size() != order_) {
		ret_code = -1;
	} else if (extrm_v.size() != order_ * 4) {
		ret_code = -2;
	} else if (mids_v.size() != ndims_ * dimsXbin_schemes_sum_) {
		ret_code = -3;
	} else if (freqs_v.size() != nstep_eff_ * bin_pow_dims_sum_) {
		ret_code = -4;
	}
	if (ret_code != 0) {
		return ret_code;
	}
	if (ncid_ == -1) {
		NC_openRead(filename_);
	}
	size_t start[3], count[3];

	start[0] = dim_index;
	start[1] = 0;
	count[0] = 1;
	count[1] = order_;

	if (checkNCerr(
			nc_get_vara_uint(ncid_, recordIdVID_, start, count,
					ids_v.data()))) {
		mprinterr("Error: Reading record_id data.\n");
		return 1;
	}

	start[0] = dim_index;
	start[1] = 0;
	start[2] = 0;
	count[0] = 1;
	count[1] = order_;
	count[2] = 4;

	if (checkNCerr(
			nc_get_vara_double(ncid_, extrmVID_, start, count,
					extrm_v.data()))) {
		mprinterr("Error: Reading record_id data.\n");
		return 1;
	}

	start[1] = 0;
	count[1] = ndims_ * dimsXbin_schemes_sum_;
	if (checkNCerr(
			nc_get_vara_double(ncid_, binMidsVID_, start, count,
					mids_v.data()))) {
		mprinterr("Error: Reading bin mid data.\n");
		return 1;
	}

	start[1] = 0;
	count[1] = nstep_eff_;

	start[2] = 0;
	count[2] = bin_pow_dims_sum_;

	if (checkNCerr(
			nc_get_vara_ulonglong(ncid_, binFreqVID_, start, count,
					freqs_v.data()))) {
		mprinterr("Error: Reading record frequency data.\n");
		return 1;
	}

	NC_close();
	return 0;
}

int Netcdf_HistUtil::readRecord(const size_t start_index, BinGroup rd) {
	std::vector<u_int> ids_v(order_);
	std::vector<double> extrm_v(order_ * 4);
	std::vector<double> mids_v(ndims_ * dimsXbin_schemes_sum_);
	std::vector<ull_int> freqs_v(nstep_eff_ * bin_pow_dims_sum_);

	setupRead();
	size_t start[3], count[3];
	start[0] = start_index;
	start[1] = 0;
	count[0] = 1;
	count[1] = order_;

	if (checkNCerr(
			nc_get_vara_uint(ncid_, recordIdVID_, start, count,
					ids_v.data()))) {
		mprinterr("Error: Reading record_id data.\n");
		return 1;
	}

	start[0] = start_index;
	start[1] = 0;
	start[2] = 0;
	count[0] = 1;
	count[1] = order_;
	count[2] = 4;

	if (checkNCerr(
			nc_get_vara_double(ncid_, extrmVID_, start, count,
					extrm_v.data()))) {
		mprinterr("Error: Reading record_id data.\n");
		return 1;
	}

	start[1] = 0;
	count[1] = ndims_ * dimsXbin_schemes_sum_;
	if (checkNCerr(
			nc_get_vara_double(ncid_, binMidsVID_, start, count,
					mids_v.data()))) {
		mprinterr("Error: Reading bin mid data.\n");
		return 1;
	}

	start[1] = 0;
	count[1] = nstep_eff_;

	start[2] = 0;
	count[2] = bin_pow_dims_sum_;

	if (checkNCerr(
			nc_get_vara_ulonglong(ncid_, binFreqVID_, start, count,
					freqs_v.data()))) {
		mprinterr("Error: Reading record frequency data.\n");
		return 1;
	}

	rd.setId(0, order_, ids_v);

	rd.setExtremes(extrm_v);

	rd.setBinMids(0, ndims_ * dimsXbin_schemes_sum_, 0, mids_v);

	rd.setBinFreqs(0, 1, nstep_eff_, 0, bin_pow_dims_sum_, freqs_v);

	NC_close();
	return 0;
}

int Netcdf_HistUtil::writeRecord(BinGroup &rec) {
	setupWrite();
	ull_int sz_id, *freqs, nsteps_eff, nbins_mid, nbins_frq;
	u_int *id;

	rec.getId(&sz_id, &id);
	size_t start[3], count[3];

	start[0] = nrecords;
	start[1] = 0;
	count[0] = 1;
	count[1] = sz_id;
	unsigned int id_u = (*id);
	if (checkNCerr(nc_put_vara_uint(ncid_, recordIdVID_, start, count, id))) {
		mprinterr("Error: Writing record data with Id=%d\n", id_u);
		return 1;
	}

	ull_int sz_ext;
	double *ext;
	rec.getExtremes(&sz_ext, &ext);
	start[0] = nrecords;
	start[1] = 0;
	start[2] = 0;
	count[0] = 1;
	count[1] = sz_id;
	count[2] = sz_ext / sz_id;
	if (checkNCerr(nc_put_vara_double(ncid_, extrmVID_, start, count, ext))) {
		mprinterr("Error: Writing record extrema data.\n");
		return 1;
	}

	double *mids;
	rec.getBinMids(&nbins_mid, &mids);
	start[1] = 0;
	count[1] = nbins_mid;
	if (checkNCerr(
			nc_put_vara_double(ncid_, binMidsVID_, start, count, mids))) {
		mprinterr("Error: Writing bin mid data.\n");
		return 1;
	}

	rec.getBinFreqs(&nsteps_eff, &nbins_frq, &freqs);
	start[1] = 0;
	count[1] = nstep_eff_;

	start[2] = 0;
	count[2] = nbins_frq;

	if (checkNCerr(
			nc_put_vara_ulonglong(ncid_, binFreqVID_, start, count, freqs))) {
		mprinterr("Error: Writing record Id data.\n");
		return 1;
	}

	NC_close();
	return 0;
}

int Netcdf_HistUtil::writeRecords(std::vector<BinGroup> &vec_rec) {

	ull_int sz_id, *freqs, nsteps_eff, nbins_mid, nbins_frq;
	u_int *id;
	for (BinGroup &rec : vec_rec) {
		setupWrite();
		rec.getId(&sz_id, &id);
		size_t start[3], count[3];

		start[0] = nrecords;
		start[1] = 0;
		count[0] = 1;
		count[1] = sz_id;
		unsigned int id_u = (*id);
		if (checkNCerr(
				nc_put_vara_uint(ncid_, recordIdVID_, start, count, id))) {
			mprinterr("Error: Writing record data with Id = %d\n", id_u);
			return 1;
		}

		ull_int sz_ext;
		double *ext;
		rec.getExtremes(&sz_ext, &ext);
		start[0] = nrecords;
		start[1] = 0;
		start[2] = 0;
		count[0] = 1;
		count[1] = sz_id;
		count[2] = sz_ext / sz_id;
		if (checkNCerr(
				nc_put_vara_double(ncid_, extrmVID_, start, count, ext))) {
			mprinterr("Error: Writing record Id data.\n");
			return 1;
		}

		double *mids;
		rec.getBinMids(&nbins_mid, &mids);
		start[1] = 0;
		count[1] = nbins_mid;
		if (checkNCerr(
				nc_put_vara_double(ncid_, binMidsVID_, start, count, mids))) {
			mprinterr("Error: Writing bin mid data.\n");
			return 1;
		}

		rec.getBinFreqs(&nsteps_eff, &nbins_frq, &freqs);
		start[1] = 0;
		count[1] = nstep_eff_;

		start[2] = 0;
		count[2] = nbins_frq;

		if (checkNCerr(
				nc_put_vara_ulonglong(ncid_, binFreqVID_, start, count,
						freqs))) {
			mprinterr("Error: Writing record Id data.\n");
			return 1;
		}

		NC_close();
	}

	return 0;
}

void Netcdf_HistUtil::printRecords() {
	size_t n_records = vec_records.size();
	for (size_t i = 0; i < n_records; ++i) {
		std::cout << vec_records[i];
	}
}

int Netcdf_HistUtil::addRecords(const std::vector<BinGroup> &recs) {
	std::vector<BinGroup>::const_iterator it;
	for (it = recs.begin(); it != recs.end(); ++it) {
		vec_records.push_back(*it);
	}
	return 0;
}

int Netcdf_HistUtil::addRecords(std::vector<BinGroup> &recs) {
	std::vector<BinGroup>::iterator it;
	for (it = recs.begin(); it != recs.end(); ++it) {
		vec_records.push_back(*it);
	}
	return 0;
}

int Netcdf_HistUtil::addRecord(BinGroup rec) {
	vec_records.push_back(rec);
	return 0;
}

int Netcdf_HistUtil::NC_create(std::string const &title) {
	return NC_create(filename_, title);
}

int Netcdf_HistUtil::NC_create(std::string const &filename,
		std::string const &title) {
	if (filename.empty())
		return 1;

	if (checkNCerr(
			nc_create(filename.c_str(), NC_NOCLOBBER | NC_64BIT_DATA, &ncid_)))
		return 1;

	// Attributes
	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, HU_GA_NC_TITLE, title.size(),
					title.c_str()))) {
		mprinterr("Error: Writing title.\n");
		return 1;
	}
	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, HU_GA_NC_APPNAME, 6, "CENTRE"))) {
		mprinterr("Error: Writing application.\n");
		return 1;
	}
	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, HU_GA_NC_PROG, 8, "Real2Int"))) {
		mprinterr("Error: Writing program.\n");
		return 1;
	}
	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, HU_GA_NC_PROGVER, 5, "1.0.0"))) {
		mprinterr("Error: Writing program version.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, HU_GA_NC_CONV, 11,
					"CENTRE_HIST"))) {
		mprinterr("Error: Writing conventions.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, HU_GA_NC_CONVVER, 3, "1.0"))) {
		mprinterr("Error: Writing conventions version.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, HU_GA_NC_CONTENT,
					contents_.size(), contents_.c_str()))) {
		mprinterr("Error: Writing contents.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, HU_GA_NC_NDIG, NC_UBYTE, 1,
					&ndiag_))) {
		mprinterr("Error: Writing ndiag.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, HU_GA_NC_ORDER, NC_UBYTE, 1,
					&order_))) {
		mprinterr("Error: Writing order .\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, HU_GA_NC_STEP_1, NC_UBYTE, 1,
					&first_step_))) {
		mprinterr("Error: Writing order .\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, HU_GA_NC_NSTEPS, NC_UBYTE, 1,
					&nsteps_))) {
		mprinterr("Error: Writing nsteps.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, HU_GA_NC_STEPSTRIDE, NC_UBYTE, 1,
					&step_stride_))) {
		mprinterr("Error: Writing step_stride.\n");
		return 1;
	}
	hbin_t store_type = (hbin_t) storage_type_;
	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, HU_GA_NC_STORAGETYPE, NC_UBYTE,
					1, &store_type))) {
		mprinterr("Error: Writing store type.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, HU_GA_NC_NDIMS, NC_UBYTE, 1,
					&ndims_))) {
		mprinterr("Error: Writing ndims.\n");
		return 1;
	}

	if (checkNCerr(nc_put_att_ulonglong(ncid_, NC_GLOBAL, HU_GA_NC_DIMLENGHS,
	NC_UINT64, dim_lens_.size(), &dim_lens_[0]))) {
		mprinterr("Error: Writing dim lengths.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, HU_GA_NC_NBINSCHEMES, NC_UBYTE,
					1, &bin_schemes_count_))) {
		mprinterr("Error: Writing binning schemes count.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, HU_GA_NC_BINSCHEMES, NC_UBYTE,
					bin_schemes_.size(), &bin_schemes_[0]))) {
		mprinterr("Error: Writing binning schemes.\n");
		return 1;
	}
	// Attribute Definitions ends here

	// Dimension definition starts
	if (checkNCerr(
			nc_def_dim(ncid_, HU_D_NC_RECORD, NC_UNLIMITED, &recordsDID_))) {
		mprinterr("Error: Writing records.\n");
		return 1;
	}

	if (checkNCerr(
			nc_def_dim(ncid_, HU_D_NC_BINMIDS, bin_dims_sum_mid_,
					&binDimsDID_))) {
		mprinterr("Error: Writing bin schemes sum.\n");
		return 1;
	}

	if (checkNCerr(
			nc_def_dim(ncid_, HU_D_NC_BINFREQS, bin_pow_dims_sum_,
					&binPowDimDID_))) {
		mprinterr("Error: Writing bin schemes sum.\n");
		return 1;
	}

	if (checkNCerr(
			nc_def_dim(ncid_, HU_D_NC_NSTEPEFF, nstep_eff_, &nstepsEffDID_))) {
		mprinterr(
				"Error (HistUtils): Defining dimenstion effective_number_of_steps.\n");
		return 1;
	}

	if (checkNCerr(
			nc_def_dim(ncid_, HU_D_NC_RECORDID, ndims_, &recordIdDID_))) {
		mprinterr("Error: Writing id dim.\n");
		return 1;
	}

	if (checkNCerr(nc_def_dim(ncid_, HU_D_NC_DIMEXTRM, 4, &extrmDID_))) {
		mprinterr("Error: Writing id dim.\n");
		return 1;
	}

	int dimRecIds[2];
	dimRecIds[0] = recordsDID_;
	dimRecIds[1] = recordIdDID_;

	if (checkNCerr(
			nc_def_var(ncid_, HU_V_NC_RECORDID, NC_UINT, 2, &dimRecIds[0],
					&recordIdVID_))) {
		mprinterr("Error: Writing mid values for bins.\n");
		return 1;
	}

	int dimExtrem[3];
	dimExtrem[0] = recordsDID_;
	dimExtrem[1] = recordIdDID_;
	dimExtrem[2] = extrmDID_;

	if (checkNCerr(
			nc_def_var(ncid_, HU_V_NC_DIMEXTRM, NC_DOUBLE, 3, &dimExtrem[0],
					&extrmVID_))) {
		mprinterr("Error: Writing mid values for bins.\n");
		return 1;
	}

	int dimsBinMids[2];
	dimsBinMids[0] = recordsDID_;
	dimsBinMids[1] = binDimsDID_;

	if (checkNCerr(
			nc_def_var(ncid_, HU_V_NC_BINMID, NC_DOUBLE, 2, &dimsBinMids[0],
					&binMidsVID_))) {
		mprinterr("Error: Writing mid values for bins.\n");
		return 1;
	}

	int dimsBinFreq[3];
	dimsBinFreq[0] = recordsDID_;
	dimsBinFreq[1] = nstepsEffDID_;
	dimsBinFreq[2] = binPowDimDID_;

	if (checkNCerr(
			nc_def_var(ncid_, HU_V_NC_BINFREQ, NC_UINT64, 3, &dimsBinFreq[0],
					&binFreqVID_))) {
		mprinterr("Error: Writing bin frequencies.\n");
		return 1;
	}

	// Set fill mode
	if (checkNCerr(nc_set_fill(ncid_, NC_NOFILL, dimsBinFreq))) {
		mprinterr("Error: NetCDF setting fill value.\n");
		return 1;
	}

	// End netcdf definitions
	if (checkNCerr(nc_enddef(ncid_))) {
		mprinterr("NetCDF error on ending definitions.");
		return 1;
	}

	NC_close();
	return 0;
}
