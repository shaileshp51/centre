#include <cstdlib>

#include <iostream>

#include "Netcdf_EntContri.h"
#include "NetcdfFile.h"

DimEntropy::DimEntropy(const u_int id, const hbin_t nsteps,
		const hbin_t nbin_schemes, const hbin_t nestimators) {
	dimid_.resize(1);
	dimid_[0] = id;
	nsteps_ = nsteps;
	nbin_schemes_ = nbin_schemes;
	nestimators_ = nestimators;
	ull_int rec_sz = nestimators_ * nsteps_ * nbin_schemes_;
	contri_S_.resize(rec_sz);
	contri_S_.assign(rec_sz, 0);
}

DimEntropy::DimEntropy(const std::vector<u_int> &id, const hbin_t nsteps,
		const hbin_t nbin_schemes, const hbin_t nestimators) {
	dimid_.assign(id.begin(), id.end());
	nsteps_ = nsteps;
	nbin_schemes_ = nbin_schemes;
	nestimators_ = nestimators;
	ull_int rec_sz = nestimators_ * nsteps_ * nbin_schemes_;
	contri_S_.resize(rec_sz);
	contri_S_.assign(rec_sz, 0);
}

void DimEntropy::getId(hbin_t *size, u_int **p) {
	*size = (hbin_t) dimid_.size();
	*p = dimid_.data();
}

void DimEntropy::getDimContri(hbin_t *nestm, hbin_t *nsteps,
		hbin_t *nbin_schemes, double **data) {
	*nestm = nestimators_;
	*nsteps = nsteps_;
	*nbin_schemes = (hbin_t) (nbin_schemes_);
	*data = contri_S_.data();
}

void DimEntropy::setId(const ull_int start_index_in, const ull_int in_count,
		const std::vector<u_int> &in) {
	std::vector<u_int>::const_iterator it;
	it = in.begin() + start_index_in;
	dimid_.assign(it, it + in_count);
}

void DimEntropy::setId(const u_int &in) {
	dimid_[0] = in;
}

//void DimEntropy::setDimContri(const hbin_t estm_id, const hbin_t step_id,
//      const hbin_t start_bin_scheme_id, const hbin_t bin_scheme_count,
//      const std::vector<double>& in) {
//   hbin_t end_ = start_bin_scheme_id + bin_scheme_count;
//   for (hbin_t i = start_bin_scheme_id, j = 0; i < end_; ++i, ++j) {
//      contri_S_[(estm_id * nbin_schemes_ * nsteps_) + ((int) i * nsteps_) + (int) step_id] = in[j];
//   }
//}

void DimEntropy::setDimContri(const hbin_t start_estm, const hbin_t count_estm,
		const hbin_t start_step_id, const hbin_t steps_count,
		const hbin_t start_bin_scheme_id, const hbin_t bin_scheme_count,
		const std::vector<double> &in) {
	hbin_t end_estm = start_estm + count_estm;
	ull_int steps_end_ = (int) start_step_id + (int) steps_count;
	ull_int bins_end_ = start_bin_scheme_id + bin_scheme_count;
	ull_int b = 0;
	for (hbin_t u = start_estm; u < end_estm; ++u) {
		ull_int est_offset = (ull_int) (u * nsteps_ * nbin_schemes_);
		for (ull_int i = start_bin_scheme_id; i < bins_end_; ++i) {
			ull_int offset = est_offset + (int) i * (int) nsteps_;
			for (ull_int a = start_step_id; a < steps_end_; ++a) {
				contri_S_[offset + a] = in[b];
				++b;
			}
		}
	}
}

void DimEntropy::setDimContri(const hbin_t start_estm, const hbin_t count_estm,
		const hbin_t start_step_id, const hbin_t steps_count,
		const hbin_t start_bin_scheme_id, const hbin_t bin_scheme_count,
		const ull_int start_in, const std::vector<double> &in) {
	hbin_t end_estm = start_estm + count_estm;
	ull_int steps_end_ = (int) start_step_id + (int) steps_count;
	ull_int bins_end_ = start_bin_scheme_id + bin_scheme_count;
	ull_int b = start_in;
	for (hbin_t u = start_estm; u < end_estm; ++u) {
		ull_int est_offset = (ull_int) (u * nsteps_ * nbin_schemes_);
		for (ull_int i = start_bin_scheme_id; i < bins_end_; ++i) {
			ull_int offset = est_offset + (int) i * (int) nsteps_;
			for (ull_int a = start_step_id; a < steps_end_; ++a) {
				contri_S_[offset + a] = in[b];
				++b;
			}
		}
	}
}

//void DimEntropy::setDimContri(hbin_t step_id, hbin_t bin_id, double value) {
//   contri_S_[(bin_id * nsteps_) + step_id] = value;
//}

std::ostream& operator<<(std::ostream &os, const DimEntropy &bg) {
	os
			<< "++++++++++++++++++++++++++++++++++++ DimEntContri Object +++++++++++++++++++++++++++++++++++"
			<< std::endl;
	os << "Record Id: {";
	for (int i = 0; i < bg.dimid_.size(); ++i) {
		os << (int) bg.dimid_[i];
		if (i + 1 != bg.dimid_.size())
			os << ", ";
	}
	os << "}" << std::endl;
	os << "N-Step: " << (int) bg.nsteps_ << std::endl;
	os << "nBin Schemes Count: " << (int) bg.nbin_schemes_ << std::endl;

	os << "Ent Contribution (nstep x nbin): {" << (int) bg.nsteps_ << ", "
			<< (int) bg.nbin_schemes_ << "}" << std::endl;
	for (int i = 0; i < bg.nbin_schemes_; ++i) {
		std::cout << "Scheme [" << (i + 1) << "]: ";
		for (int j = 0; j < bg.nsteps_; ++j) {
			os << (double) bg.contri_S_[i * bg.nsteps_ + j];
			if ((j + 1) % bg.nsteps_ != 0)
				os << ", ";
		}
		std::cout << std::endl;
	}
	os
			<< "====================================== DimEntContri Object End ======================================="
			<< std::endl;
	return os;
}

Netcdf_EntContri::Netcdf_EntContri(const std::string &filename,
		const std::string &c, const hbin_t ndims, const hbin_t order,
		const hbin_t nsteps, const std::vector<ull_int> &dim_lens,
		const TensorType &typ, const std::vector<hbin_t> &bin_schemes,
		const hbin_t n_estimators) {
	filename_ = std::string(filename.c_str());
	contents_ = std::string(c.c_str());
	ndims_ = (ndims < 1) ? 1 : ndims;
	order_ = (order < 1) ? 1 : order;
	nsteps_ = (nsteps < 1) ? 1 : nsteps;
	bin_schemes_count_ = bin_schemes.size();
	nestimators_ = n_estimators;
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

	bin_schemes_.reserve(bin_schemes_count_);
	for (size_t i = 0; i < bin_schemes_count_; ++i) {
		bin_schemes_.push_back(bin_schemes[i]);
	}
	recordsDID_ = -1;
	estimDID_ = -1;
	binSchemesDID_ = -1;
	nstepsDID_ = -1;
	recordIdDID_ = -1;

	// variables for storing variable ids
	contri_SVID_ = -1;
	recordIdVID_ = -1;
}

std::ostream& operator<<(std::ostream &os, const Netcdf_EntContri &ht) {
	os
			<< "++++++++++++++++++++++++++++++++++++ Netcdf_EntContri Object +++++++++++++++++++++++++++++++++++"
			<< std::endl;
	os << "Filename: " << ht.filename_ << std::endl;
	os << "Content: " << ht.contents_ << std::endl;
	os << "N-Dims: " << (int) ht.ndims_ << std::endl;
	os << "Order: " << (int) ht.order_ << std::endl;
	os << "N-Step: " << (int) ht.nsteps_ << std::endl;

	os << "Bin Schemes Count: " << (int) ht.bin_schemes_count_ << std::endl;
	os << "Bin Schemes: ";
	for (int i = 0; i < ht.bin_schemes_count_; ++i) {
		os << (int) ht.bin_schemes_[i];
		if (i + 1 != ht.bin_schemes_count_)
			os << ", ";
	}
	os << std::endl;

	os << "Storage Type: " << ht.storage_type_ << std::endl;
	os << "N-Diag: " << (int) ht.ndiag_ << std::endl;
	os
			<< "====================================== Netcdf_EntContri Object End ======================================="
			<< std::endl;
	return os;
}

int Netcdf_EntContri::removeRecord(const DimEntropy &rec) {
	return 0;
}

int Netcdf_EntContri::removeRecords(const std::vector<DimEntropy> &recs) {
	return 0;
}

void Netcdf_EntContri::clearRecords() {
	vec_records.clear();
}

int Netcdf_EntContri::setupRead() {
	if (ncmode_ != NC_NOWRITE && ncid_ != -1) {
		NC_close();
	}
	if (ncid_ == -1) {
		NC_openRead(filename_);
		int temp;
		recordsDID_ = getDimInfo(EC_D_NC_RECORD, &temp);
		nrecords = temp;
		recordIdDID_ = getDimInfo(EC_D_NC_RECORDID, &temp);
		order_ = temp;
		estimDID_ = getDimInfo("estimators", &temp);
		nestimators_ = temp;
		binSchemesDID_ = getDimInfo(EC_D_NC_NBINSCHEMES, &temp);
		bin_schemes_count_ = temp;
		nstepsDID_ = getDimInfo(EC_D_NC_NSTEPS, &temp);
		nsteps_ = temp;

		// Get coord info
		recordIdVID_ = -1;
		if (nc_inq_varid(ncid_, EC_V_NC_RECORDID, &recordIdVID_) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file has recordIds.\n");
		}

		contri_SVID_ = -1;
		if (nc_inq_varid(ncid_, EC_V_NC_ENTCONTRI, &contri_SVID_) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file has frequencies.\n");
		}
	}

	return 0;
}

int Netcdf_EntContri::setupWrite() {
	if (ncmode_ != NC_WRITE && ncid_ != -1) {
		NC_close();
	}
	if (ncid_ == -1) {
		NC_openWrite(filename_);
		int temp;
		recordsDID_ = getDimInfo(EC_D_NC_RECORD, &temp);
		nrecords = temp;
		recordIdDID_ = getDimInfo(EC_D_NC_RECORDID, &temp);
		order_ = temp;
		estimDID_ = getDimInfo("estimators", &temp);
		nestimators_ = temp;
		binSchemesDID_ = getDimInfo(EC_D_NC_NBINSCHEMES, &temp);
		bin_schemes_count_ = temp;
		nstepsDID_ = getDimInfo(EC_D_NC_NSTEPS, &temp);
		nsteps_ = temp;

		// Get coord info
		recordIdVID_ = -1;
		if (nc_inq_varid(ncid_, EC_V_NC_RECORDID, &recordIdVID_) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file has recordIds.\n");
		}

		contri_SVID_ = -1;
		if (nc_inq_varid(ncid_, EC_V_NC_ENTCONTRI, &contri_SVID_) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file has frequencies.\n");
		}
	}
	return 0;
}

int Netcdf_EntContri::readRecords(const size_t start_index,
		const size_t rec_counts, bool is_append) {
	std::vector<DimEntropy> rd(rec_counts,
			DimEntropy(0, nsteps_, bin_schemes_count_, nestimators_));
	std::vector<u_int> ids_v(rec_counts * order_);
	std::vector<double> contri_S_v(rec_counts * nsteps_ * bin_schemes_count_);
#ifdef DEBUG_CODE
	std::cout << "Freq: (" << contri_S_v.size() << ", " << contri_S_v.size()
			<< ") " << std::endl;
#endif

	setupRead();

	size_t start[4], count[4];

	start[0] = start_index;
	count[0] = rec_counts;

	start[1] = 0;
	count[1] = order_;

	if (checkNCerr(
			nc_get_vara_uint(ncid_, recordIdVID_, start, count,
					ids_v.data()))) {
		mprinterr("Error: Reading record_id data.\n");
		return 1;
	}

	//std::cout << "Record ID: start: [" << start[0] << ", " << start[1] << "] count: [" << count[0] << ", " << count[1] << "]" << std::endl;
	start[1] = 0;
	count[1] = nestimators_;

	start[2] = 0;
	count[2] = bin_schemes_count_;

	start[3] = 0;
	count[3] = nsteps_;

	if (checkNCerr(
			nc_get_vara_double(ncid_, contri_SVID_, start, count,
					contri_S_v.data()))) {
		mprinterr("Error: Reading record frequency data.\n");
		return 1;
	}
	size_t k = 0;
	for (size_t i = 0; i < rec_counts; ++i) {
		rd[i].setId(k, order_, ids_v);
		k += order_;
	}
	k = 0;
	for (size_t i = 0; i < rec_counts; ++i) {
		rd[i].setDimContri(0, nestimators_, 0, nsteps_, 0, bin_schemes_count_,
				k, contri_S_v);
		k += (nestimators_ * nsteps_ * bin_schemes_count_);
	}
#ifdef DEBUG_CODE_2
	std::cout << "Record Bin Contri: start: [" << start[0] << ", " << start[1] << ", " << start[2];
	std::cout << "] count: [" << count[0] << ", " << count[1] << ", " << count[2] << "]"
	<< std::endl;
#endif
	if (is_append) {
		addRecords(rd);
	} else {
		clearRecords();
		addRecords(rd);
	}
	return 0;
}

int Netcdf_EntContri::readEntrContrib(const u_int start_index,
		const u_int rec_counts, const hbin_t start_estm, const hbin_t n_estm,
		const hbin_t start_schemes, const hbin_t n_schemes,
		const hbin_t start_step, const hbin_t nsteps,
		std::vector<double> &dataOut) {

	setupRead();
	size_t start[4], count[4];

	start[0] = start_index;
	count[0] = rec_counts;

	start[1] = start_estm;
	count[1] = n_estm;

	start[2] = start_schemes;
	count[2] = n_schemes;

	start[3] = start_step;
	count[3] = nsteps;

	if (checkNCerr(
			nc_get_vara_double(ncid_, contri_SVID_, start, count,
					dataOut.data()))) {
		mprinterr("Error: Reading entropy contribution data.\n");
		return 1;
	}

	return 0;
}

int Netcdf_EntContri::readIds(const u_int start_index, const u_int rec_counts,
		std::vector<u_int> &ids_v) {

	setupRead();
	size_t start[2], count[2];

	start[0] = start_index;
	count[0] = rec_counts;

	start[1] = 0;
	count[1] = order_;

	if (checkNCerr(
			nc_get_vara_uint(ncid_, recordIdVID_, start, count,
					ids_v.data()))) {
		mprinterr("Error: Reading record_id data.\n");
		return 1;
	}

	return 0;
}

int Netcdf_EntContri::readEntrContrib(const u_int start_index,
		const u_int rec_counts, const hbin_t start_estm, const hbin_t n_estm,
		const hbin_t start_schemes, const hbin_t n_schemes,
		const hbin_t start_step, const hbin_t nsteps, const hbin_t order,
		std::vector<u_int> &ids_v, std::vector<double> &dataOut) {

	setupRead();
	size_t start[4], count[4];

	start[0] = start_index;
	count[0] = rec_counts;

	start[1] = 0;
	count[1] = order;

	if (checkNCerr(
			nc_get_vara_uint(ncid_, recordIdVID_, start, count,
					ids_v.data()))) {
		mprinterr("Error: Reading record_id data.\n");
		return 1;
	}

	start[1] = start_estm;
	count[1] = n_estm;

	start[2] = start_schemes;
	count[2] = n_schemes;

	start[3] = start_step;
	count[3] = nsteps;

	if (checkNCerr(
			nc_get_vara_double(ncid_, contri_SVID_, start, count,
					dataOut.data()))) {
		mprinterr("Error: Reading record frequency data.\n");
		return 1;
	}

	return 0;
}

int Netcdf_EntContri::writeRecord(DimEntropy &rec) {
	setupWrite();
	hbin_t szid, nestm, nsteps, nbin_schemes, nbins_frq;
	double *contri;
	u_int *id;
	rec.getId(&szid, &id);
	size_t start[4], count[4];
#ifdef DEBUG_CODE_2
	std::cout << "Writing EntContri Records in file: " << nrecords << ", order: " << (int) szid << " Id: " << *id
	<< std::endl;
#endif
	start[0] = nrecords;
	count[0] = 1;

	start[1] = 0;
	count[1] = szid;

	if (checkNCerr(nc_put_vara_uint(ncid_, recordIdVID_, start, count, id))) {
		mprinterr("Error: Writing record Id data.\n");
		return 1;
	} else {
#ifdef DEBUG_CODE_2
		std::cout << " Variable ID written successfully...";
#endif
	}

	rec.getDimContri(&nestm, &nsteps, &nbin_schemes, &contri);
	start[1] = 0;
	count[1] = nestimators_;

	start[2] = 0;
	count[2] = bin_schemes_count_;

	start[3] = 0;
	count[3] = nsteps_;
#ifdef DEBUG_CODE_2
	std::cout << "Record ID: start: [" << start[0] << ", " << start[1] << ", " << start[2];
	std::cout << "] count: [" << count[0] << ", " << count[1] << ", " << count[2] << "]"
	<< std::endl;
#endif
	if (checkNCerr(
			nc_put_vara_double(ncid_, contri_SVID_, start, count, contri))) {
		mprinterr("Error: Writing entropy contribution data.\n");
		return 1;
	} else {
#ifdef DEBUG_CODE_2
		std::cout << " Contribution written successfully...";
#endif
	}
#ifdef DEBUG_CODE_2
	std::cout << "Record Id: " << *id;
#endif
	NC_close();
	return 0;
}

int Netcdf_EntContri::writeRecords(std::vector<DimEntropy> &recs) {
	ull_int nrecs = recs.size();
	if (nrecs > 0) {
		setupWrite();

		u_int *id;
		hbin_t szid, nestm, nsteps, nbin_schemes;
		double *contri_S;
		(recs[0]).getId(&szid, &id);
		size_t start[4], count[4];
#ifdef DEBUG_CODE
		std::cout << "Writing EntContriRecords in file: nrecords x recdim: ["
				<< nrecs << ", " << (int) szid << "]" << std::endl;
#endif
		std::vector<u_int> ids_v(szid * nrecs);
		size_t k = 0;
		for (size_t i = 0; i < nrecs; ++i) {
			(recs[i]).getId(&szid, &id);
			for (size_t j = 0; j < szid; ++j, ++k) {
				ids_v[k] = *(id + j);
				// ++id;
			}
		}
		start[0] = nrecords;
		count[0] = nrecs;

		start[1] = 0;
		count[1] = szid;

#ifdef DEBUG_CODE
		std::cout << "Record ID: start: [" << start[0] << ", " << start[1]
				<< "] count: [" << count[0] << ", " << count[1] << "]: "
				<< ids_v.size() << std::endl;
#endif

		if (checkNCerr(
				nc_put_vara_uint(ncid_, recordIdVID_, start, count,
						ids_v.data()))) {
			mprinterr("Error: Writing record Id data.\n");
			return 1;
		}

		recs[0].getDimContri(&nestm, &nsteps, &nbin_schemes, &contri_S);
		std::vector<double> contri_S_v(nrecs * nestm * nbin_schemes * nsteps);
		k = 0;
		for (size_t i = 0; i < nrecs; ++i) {
			recs[i].getDimContri(&nestm, &nsteps, &nbin_schemes, &contri_S);
			int contri_index = 0;
			for (size_t m = 0; m < nestm; ++m) {
				for (size_t j = 0; j < nbin_schemes; ++j) {
					for (size_t m = 0; m < nsteps; ++m) {
						contri_S_v[k] = *(contri_S + contri_index);
						//++contri_S;
						contri_index++;
						k++;
					}
				}
			}
		}

		start[1] = 0;
		count[1] = nestm;

		start[2] = 0;
		count[2] = nbin_schemes;

		start[3] = 0;
		count[3] = nsteps;

#ifdef DEBUG_CODE
		std::cout << "Entropy Record Data: start: [" << start[0] << ", "
				<< start[1] << ", " << start[2] << ", " << start[3];
		std::cout << "] count: [" << count[0] << ", " << count[1] << ", "
				<< count[2] << ", " << count[3] << "] " << contri_S_v.size()
				<< std::endl;
#endif

		if (checkNCerr(
				nc_put_vara_double(ncid_, contri_SVID_, start, count,
						contri_S_v.data()))) {
			mprinterr("Error: Writing record data.\n");
			return 1;
		}
#ifdef DEBUG_CODE
		std::cout << "Record Bin Freqs: start: [" << start[0] << ", "
				<< start[1] << ", " << start[2];
		std::cout << "] count: [" << count[0] << ", " << count[1] << ", "
				<< count[2] << "]" << std::endl;
		// std::cout << "Record Id: " << *id << std::endl;
#endif
		NC_close();
	}
	return 0;
}

void Netcdf_EntContri::printRecords() {
	size_t n_records = vec_records.size();
	for (size_t i = 0; i < n_records; ++i) {
		std::cout << vec_records[i];
	}
}

int Netcdf_EntContri::writeRecords() {
	int ret = writeRecords(vec_records);
	return ret;
}

int Netcdf_EntContri::addRecords(const std::vector<DimEntropy> &recs) {
	std::vector<DimEntropy>::const_iterator it;
	for (it = recs.begin(); it != recs.end(); ++it) {
		vec_records.push_back(*it);
	}
	return 0;
}

int Netcdf_EntContri::addRecords(std::vector<DimEntropy> &recs) {
	std::vector<DimEntropy>::iterator it;
	for (it = recs.begin(); it != recs.end(); ++it) {
		vec_records.push_back(*it);
	}
	return 0;
}

int Netcdf_EntContri::addRecord(DimEntropy rec) {
	vec_records.push_back(rec);
	return 0;
}

int Netcdf_EntContri::NC_create(std::string const &title) {
	return NC_create(filename_, title);
}

int Netcdf_EntContri::NC_create(std::string &filename,
		std::string const &title) {
	if (filename.empty())
		return 1;
	int dimensionID[NC_MAX_VAR_DIMS];

	nc_type dataType;
	NCTYPE type = NCTYPE::NC_CENTREHIST;
	/*
	 if (ncdebug_>1)
	 mprintf("DEBUG: NC_create: %s  natom=%i V=%i  box=%i  temp=%i  time=%i\n",
	 Name.c_str(),natomIn,(int)coordInfo.HasVel(),(int)coordInfo.HasBox(),
	 (int)coordInfo.HasTemp(),(int)coordInfo.HasTime()); NETCDF3_64BIT_DATA
	 * */
	if (checkNCerr(
			nc_create(filename.c_str(), NC_NOCLOBBER | NC_NETCDF4, &ncid_)))
		return 1;

	// Attributes
	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, EC_GA_NC_TITLE, title.size(),
					title.c_str()))) {
		mprinterr("Error: Writing title.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, EC_GA_NC_APPNAME, 6, "CENTRE"))) {
		mprinterr("Error: Writing application.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, EC_GA_NC_PROG, 8, "Real2Int"))) {
		mprinterr("Error: Writing program.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, EC_GA_NC_PROGVER, 5, "1.0.0"))) {
		mprinterr("Error: Writing program version.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, EC_GA_NC_CONV, 11,
					"EntropyContibution"))) {
		mprinterr("Error: Writing conventions.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, EC_GA_NC_CONVVER, 3, "1.0"))) {
		mprinterr("Error: Writing conventions version.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, EC_GA_NC_CONTENT,
					contents_.size(), contents_.c_str()))) {
		mprinterr("Error: Writing contents.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, EC_GA_NC_NDIG, NC_UBYTE, 1,
					&ndiag_))) {
		mprinterr("Error: Writing ndiag.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, EC_GA_NC_ORDER, NC_UBYTE, 1,
					&order_))) {
		mprinterr("Error: Writing order .\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, EC_GA_NC_NSTEPS, NC_UBYTE, 1,
					&nsteps_))) {
		mprinterr("Error: Writing nsteps.\n");
		return 1;
	}

	hbin_t store_type = (hbin_t) storage_type_;
	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, EC_GA_NC_STORAGETYPE, NC_UBYTE,
					1, &store_type))) {
		mprinterr("Error: Writing store type.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, EC_GA_NC_NDIMS, NC_UBYTE, 1,
					&ndims_))) {
		mprinterr("Error: Writing ndims.\n");
		return 1;
	}

	if (checkNCerr(nc_put_att_ulonglong(ncid_, NC_GLOBAL, EC_GA_NC_DIMLENGHS,
	NC_UINT64, dim_lens_.size(), &dim_lens_[0]))) {
		mprinterr("Error: Writing dim lengths.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, EC_GA_NC_NBINSCHEMES, NC_UBYTE,
					1, &bin_schemes_count_))) {
		mprinterr("Error: Writing binning schemes count.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, EC_GA_NC_BINSCHEMES, NC_UBYTE,
					bin_schemes_.size(), &bin_schemes_[0]))) {
		mprinterr("Error: Writing binning schemes.\n");
		return 1;
	}
	// Attribute Definitions ends here

	// Dimension definition starts
	if (checkNCerr(
			nc_def_dim(ncid_, EC_D_NC_RECORD, NC_UNLIMITED, &recordsDID_))) {
		mprinterr("Error: Writing records.\n");
		return 1;
	}

	if (checkNCerr(nc_def_dim(ncid_, "estimators", nestimators_, &estimDID_))) {
		mprinterr("Error: Writing bin schemes sum.\n");
		return 1;
	}

	if (checkNCerr(
			nc_def_dim(ncid_, EC_D_NC_NBINSCHEMES, bin_schemes_count_,
					&binSchemesDID_))) {
		mprinterr("Error: Writing bin schemes sum.\n");
		return 1;
	}

	if (checkNCerr(nc_def_dim(ncid_, EC_D_NC_NSTEPS, nsteps_, &nstepsDID_))) {
		mprinterr("Error: Defining dimension number_of_steps.\n");
		return 1;
	}

	if (checkNCerr(
			nc_def_dim(ncid_, EC_D_NC_RECORDID, order_, &recordIdDID_))) {
		mprinterr("Error: Writing id dim.\n");
		return 1;
	}

	int dimRecIds[2];
	dimRecIds[0] = recordsDID_;
	dimRecIds[1] = recordIdDID_;
	if (checkNCerr(
			nc_def_var(ncid_, EC_V_NC_RECORDID, NC_UINT, 2, &dimRecIds[0],
					&recordIdVID_))) {
		mprinterr("Error: Writing mid values for bins.\n");
		return 1;
	}

	int dimsContriS[4];
	dimsContriS[0] = recordsDID_;
	dimsContriS[1] = estimDID_;
	dimsContriS[2] = binSchemesDID_;
	dimsContriS[3] = nstepsDID_;

	unsigned long chunkSizes[4];
	chunkSizes[0] = 16;
	chunkSizes[1] = nestimators_;
	chunkSizes[2] = bin_schemes_count_;
	chunkSizes[3] = nsteps_;

	if (checkNCerr(
			nc_def_var(ncid_, EC_V_NC_ENTCONTRI, NC_DOUBLE, 4, &dimsContriS[0],
					&contri_SVID_))) {
		mprinterr("Error: Writing bin frequencies.\n");
		return 1;
	}

	if (checkNCerr(
			nc_def_var_chunking(ncid_, contri_SVID_, NC_CHUNKED,
					&chunkSizes[0]))) {
		mprinterr("Error: Writing bin frequencies.\n");
		return 1;
	}

	// Set fill mode
	if (checkNCerr(nc_set_fill(ncid_, NC_NOFILL, dimsContriS))) {
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
