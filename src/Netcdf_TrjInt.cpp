#include <cstdlib>
#include <sstream>
#include <iostream>

#include "Netcdf_TrjInt.h"
#include "Exceptions.h"
#include "NetcdfFile.h"

using namespace std;

CoordBAT::CoordBAT(const u_int id, const hbin_t bin_schemes_count) {
	dimid_ = id;
	bin_schemes_count_ = bin_schemes_count;
	n_frame_ = 0;
}

CoordBAT::CoordBAT(const u_int id, const hbin_t bin_schemes_count,
		const ull_int nframes) {
	dimid_ = id;
	bin_schemes_count_ = bin_schemes_count;
	n_frame_ = nframes;
	extremas_.assign(4, 0.0); // for min, mav, avg & phase
	coord_.assign(bin_schemes_count * n_frame_, 0);
}

u_int CoordBAT::getId() const {
	return dimid_;
}

void CoordBAT::setId(const u_int id) {
	dimid_ = id;
}

void CoordBAT::getCoords(double *min_v, double *max_v, double *avg_v,
		double *phase_v, hbin_t *bin_schemes_count, ull_int *nframes,
		hbin_t **data) {
	*min_v = extremas_[0]; // minimum value
	*max_v = extremas_[1]; // maximum value
	*avg_v = extremas_[2]; // average value
	*phase_v = extremas_[3]; // phase value, if any otherwise 0.0
	*bin_schemes_count = bin_schemes_count_;
	*nframes = n_frame_;
	*data = coord_.data();
}

int CoordBAT::setCoords(const double min_v, const double max_v,
		const double avg_v, const double phase_v, const ull_int start_binscheme,
		const ull_int schemes_count, const ull_int start_frame,
		const ull_int frame_end, const std::vector<hbin_t> &in) {

	extremas_[0] = min_v;
	extremas_[1] = max_v;
	extremas_[2] = avg_v;
	extremas_[3] = phase_v;
	ull_int bins_end_ = start_binscheme + schemes_count;
	ull_int b = 0;
	for (ull_int i = start_binscheme; i < bins_end_; ++i) {
		ull_int offset = i * n_frame_;
		for (ull_int a = start_frame; a <= frame_end; ++a) {
			coord_[offset + a] = in[b];
			++b;
		}
	}
	return 0;
}

int CoordBAT::setCoords(const double min_v, const double max_v,
		const double avg_v, const double phase_v, const ull_int start_binscheme,
		const ull_int schemes_count, const ull_int start_frame,
		const ull_int frame_stride, const ull_int frame_end,
		const std::vector<hbin_t> &in) {
	extremas_[0] = min_v;
	extremas_[1] = max_v;
	extremas_[2] = avg_v;
	extremas_[3] = phase_v;
	ull_int bins_end_ = start_binscheme + schemes_count;
	ull_int b = 0;
	for (ull_int i = start_binscheme; i < bins_end_; ++i) {
		ull_int offset = i * n_frame_;
		for (ull_int a = start_frame; a <= frame_end; a += frame_stride) {
			coord_[offset + a] = in[b];
			++b;
		}
	}
	return 0;
}

std::ostream& operator<<(std::ostream &os, const CoordBAT &bg) {
	os
			<< "++++++++++++++++++++++++++++++++++++ BinGroup Object +++++++++++++++++++++++++++++++++++"
			<< std::endl;
	os << "Dim Id: " << bg.getId();
	os << "Bin Schemes Count: " << (int) bg.bin_schemes_count_ << std::endl;
	os << "Bin Data (schemes x nframe): {" << (int) bg.bin_schemes_count_
			<< ", " << (int) bg.n_frame_ << "}" << std::endl;
	for (auto i = 0; i < bg.bin_schemes_count_; ++i) {
		std::cout << "Scheme [" << (i + 1) << "]: ";
		for (auto j = 0; j < bg.n_frame_; ++j) {
			os << (ull_int) bg.coord_[i * bg.bin_schemes_count_ + j];
			if ((j + 1) % bg.n_frame_ != 0)
				os << ", ";
		}
		std::cout << std::endl;
	}
	os
			<< "====================================== BinGroup Object End ======================================="
			<< std::endl;
	return os;
}

Netcdf_TrjInt::Netcdf_TrjInt(const std::string &filename, const std::string &c,
		const u_int ndims, const ull_int start_frame,
		const ull_int frame_stride, const ull_int nframes,
		const hbin_t bin_schemes_count, const std::vector<hbin_t> &bin_schemes) :
		NetcdfFile(filename) {
	filename_ = std::string(filename.c_str());
	contents_ = std::string(c.c_str());
	ndims_ = (ndims < 1) ? 1 : ndims;
	n_frame_ = nframes;
	start_frame_ = (n_frame_ <= start_frame) ? 0 : start_frame;
	frame_stride_ = (n_frame_ < frame_stride) ? 1 : frame_stride;
	n_frame_eff_ = (n_frame_ - start_frame_ + frame_stride_ - 1)
			/ frame_stride_;
	bin_schemes_count_ = bin_schemes_count;
	nrecords = 0;
	nwrittendims_ = 0;
	bin_schemes_.reserve(bin_schemes_count_);
	for (size_t i = 0; i < bin_schemes_count_; ++i) {
		bin_schemes_.push_back(bin_schemes[i]);
	}
	binSchemesDID_ = -1;
	frameDID_ = -1;
	dimDID_ = -1;
	extrDID_ = -1;
	dimExtrVID_ = -1;
	dimIdVID_ = -1;
	coordVID_ = -1;
	ncid_ = -1;
}

std::ostream& operator<<(std::ostream &os, const Netcdf_TrjInt &ht) {
	os
			<< "++++++++++++++++++++++++++++++++++++ Netcdf_TrjInt Object +++++++++++++++++++++++++++++++++++"
			<< std::endl;
	os << "Filename: " << ht.filename_ << std::endl;
	os << "Content: " << ht.contents_ << std::endl;
	os << "N-Dims: " << (int) ht.ndims_ << std::endl;
	os << "N-Frames: " << ht.n_frame_ << std::endl;
	os << "Bin Schemes Count: " << (int) ht.bin_schemes_count_ << std::endl;
	os << "Bin Schemes: ";
	for (int i = 0; i < ht.bin_schemes_count_; ++i) {
		os << (int) ht.bin_schemes_[i];
		if (i + 1 != ht.bin_schemes_count_)
			os << ", ";
	}
	os << std::endl;
	os
			<< "====================================== Netcdf_TrjInt Object End ======================================="
			<< std::endl;
	return os;
}

int Netcdf_TrjInt::removeRecord(const CoordBAT &rec) {
	return 0;
}

int Netcdf_TrjInt::removeRecords(const std::vector<CoordBAT> &recs) {
	return 0;
}

void Netcdf_TrjInt::clearRecords() {
	vec_records.clear();
}

int Netcdf_TrjInt::setupRead() {
	if (ncmode_ != NC_NOWRITE && (ncid_ != -1)) {
		NC_close();
	}
	if (ncid_ == -1) {
		NC_openRead(filename_);

		int temp;
		dimDID_ = getDimInfo(TI_D_NC_DIMS, &temp);
		ndims_ = temp;
		extrDID_ = getDimInfo(TI_D_NC_DIMEXTM, &temp);
		frameDID_ = getDimInfo(TI_D_NC_FRAME, &temp);
		n_frame_ = temp;
		binSchemesDID_ = getDimInfo(TI_D_NC_BINSCHEMESCOUNT, &temp);
		bin_schemes_count_ = temp;

		// Get coord info
		dimIdVID_ = -1;
		if (nc_inq_varid(ncid_, TI_V_NC_DIMID, &dimIdVID_) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file has bin mids.\n");
		}
		dimExtrVID_ = -1;
		if (nc_inq_varid(ncid_, TI_V_NC_DIMEXTM, &dimExtrVID_) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file has bin mids.\n");
		}
		coordVID_ = -1;
		if (nc_inq_varid(ncid_, TI_V_NC_COORD, &coordVID_) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file has frequencies.\n");
		}
	}
	return 0;
}

int Netcdf_TrjInt::setupWrite() {
	if (ncmode_ != NC_WRITE && ncid_ != -1) {
		NC_close();
	}
	if (ncid_ == -1) {
		NC_openWrite(filename_);
	}
	int temp;
	dimDID_ = getDimInfo(TI_D_NC_DIMS, &temp);
	nwrittendims_ = temp;
	extrDID_ = getDimInfo(TI_D_NC_DIMEXTM, &temp);
	frameDID_ = getDimInfo(TI_D_NC_FRAME, &temp);
	n_frame_ = temp;
	binSchemesDID_ = getDimInfo(TI_D_NC_BINSCHEMESCOUNT, &temp);
	bin_schemes_count_ = temp;

	// Get coord info
	dimIdVID_ = -1;
	if (nc_inq_varid(ncid_, TI_V_NC_DIMID, &dimIdVID_) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf file has dimids.\n");
	}
	dimExtrVID_ = -1;
	if (nc_inq_varid(ncid_, TI_V_NC_DIMEXTM, &dimExtrVID_) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf file has bin mids.\n");
	}
	coordVID_ = -1;
	if (nc_inq_varid(ncid_, TI_V_NC_COORD, &coordVID_) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf file has frequencies.\n");
	}

	return 0;
}

int Netcdf_TrjInt::readCoords(const u_int dim_idx, const hbin_t bin_id,
		CoordBAT &rd) {

	setupRead();

	size_t start[3], count[3];
	ptrdiff_t stride[3];
	std::vector<double> extrm_v(4);
	std::vector<hbin_t> coords_v(1 * n_frame_eff_);
	u_int dimid_v;
	start[0] = dim_idx;
	count[0] = 1;
	stride[0] = 1;

	if (checkNCerr(
			nc_get_vara_uint(ncid_, dimIdVID_, start, count, &dimid_v))) {
		ostringstream oss;
		oss << "NC_ERR>> Reading DimensionID '" << dimIdVID_
				<< "' data from file: '" << ncfilename_ << "'";
		throw NCDataReadError(oss.str());

		return 1;
	}

	start[1] = 0;
	count[1] = 4;

	if (checkNCerr(
			nc_get_vara_double(ncid_, dimExtrVID_, start, count,
					extrm_v.data()))) {
		ostringstream oss;
		oss << "NC_ERR>> Reading Extremas data from file: '" << ncfilename_
				<< "'";
		throw NCDataReadError(oss.str());

		return 1;
	}

	start[1] = bin_id;
	count[1] = 1;
	stride[1] = 1;

	start[2] = start_frame_;
	count[2] = n_frame_;
	stride[2] = frame_stride_;

	if (checkNCerr(
			nc_get_vars_uchar(ncid_, coordVID_, start, count, stride,
					coords_v.data()))) {
		ostringstream oss;
		oss << "NC_ERR>> Reading coords data from file: '" << ncfilename_
				<< "'";
		throw NCDataReadError(oss.str());

		return 1;
	}
	//size_t k = 0;
	rd.setId(dimid_v);
	rd.setCoords(extrm_v[0], extrm_v[1], extrm_v[2], extrm_v[3], 0, 1, 0,
			n_frame_eff_ - 1, coords_v);

	return 0;
}

int Netcdf_TrjInt::readCoords(const u_int dim_idx,
		const hbin_t binscheme_start_id, const hbin_t binscheme_count,
		CoordBAT &rd) {

	setupRead();

	size_t start[3], count[3];
	ptrdiff_t stride[3];
	std::vector<double> extrm_v(4);
	std::vector<hbin_t> coords_v(binscheme_count * n_frame_eff_);
	u_int dimid_v;
	start[0] = dim_idx;
	count[0] = 1;
	stride[0] = 1;

	if (checkNCerr(
			nc_get_vara_uint(ncid_, dimIdVID_, start, count, &dimid_v))) {
		ostringstream oss;
		oss << "NC_ERR>> Reading DimensionID '" << dimIdVID_ << "' from file: '"
				<< ncfilename_ << "'";
		throw NCDataReadError(oss.str());

		return 1;
	}

	start[1] = 0;
	count[1] = 4;

	if (checkNCerr(
			nc_get_vara_double(ncid_, dimExtrVID_, start, count,
					extrm_v.data()))) {
		ostringstream oss;
		oss << "NC_ERR>> Reading extremas from file: '" << ncfilename_ << "'";
		throw NCDataReadError(oss.str());

		return 1;
	}

	start[1] = binscheme_start_id;
	count[1] = binscheme_count;
	stride[1] = 1;

	start[2] = start_frame_;
	count[2] = n_frame_eff_;
	stride[2] = frame_stride_;

	if (checkNCerr(
			nc_get_vars_uchar(ncid_, coordVID_, start, count, stride,
					coords_v.data()))) {
		ostringstream oss;
		oss << "NC_ERR>> Reading coords from file: '" << ncfilename_ << "'";
		throw NCDataReadError(oss.str());

		return 1;
	}

	rd.setId(dimid_v);
	rd.setCoords(extrm_v[0], extrm_v[1], extrm_v[2], extrm_v[3],
			binscheme_start_id, binscheme_count, 0, n_frame_eff_ - 1, coords_v);

	// NC_close();

	return 0;
}

int Netcdf_TrjInt::writeRecord(CoordBAT &rec) {

	setupWrite();

	ull_int nframes;
	u_int id;
	hbin_t *data, bin_schemes_count;
	std::vector<double> extrm_v(4);

	id = rec.getId();
	rec.getCoords(&extrm_v[0], &extrm_v[1], &extrm_v[2], &extrm_v[3],
			&bin_schemes_count, &nframes, &data);
	size_t start[3], count[3];

	start[0] = nwrittendims_;
	count[0] = 1;

	if (checkNCerr(nc_put_vara_uint(ncid_, dimIdVID_, start, count, &id))) {
		ostringstream oss;
		oss << "NC_ERR>> Writing record-id to file: '" << ncfilename_ << "'";
		throw NCDataWriteError(oss.str());

		return 1;
	}

	start[1] = 0;
	count[1] = 4;

	if (checkNCerr(
			nc_put_vara_double(ncid_, dimExtrVID_, start, count,
					extrm_v.data()))) {
		ostringstream oss;
		oss << "NC_ERR>> Writing extremas data to file: '" << ncfilename_
				<< "'";
		throw NCDataWriteError(oss.str());

		return 1;
	}

	start[1] = 0;
	count[1] = bin_schemes_count;

	start[2] = 0;
	count[2] = nframes;

	if (checkNCerr(nc_put_vara_uchar(ncid_, coordVID_, start, count, data))) {
		ostringstream oss;
		oss << "NC_ERR>> Writing coords data to file: '" << ncfilename_ << "'";
		throw NCDataWriteError(oss.str());

		return 1;
	}

	NC_close();
	return 0;
}

int Netcdf_TrjInt::writeRecords(std::vector<CoordBAT> &vec_recs,
		ull_int n_recs /* = 0 */) {
	ull_int nrecs = n_recs;
	if (n_recs == 0) {
		nrecs = vec_recs.size();
	}
	if (nrecs > 0) {
		setupWrite();

		ull_int nframes;
		std::vector<u_int> ids_v(nrecs, -1);
		hbin_t bin_schemes_count;
		hbin_t *data;
		std::vector<double> extrm_v(4 * nrecs);

		(vec_recs[0]).getCoords(&extrm_v[0], &extrm_v[1], &extrm_v[2],
				&extrm_v[3], &bin_schemes_count, &nframes, &data);

		std::vector<hbin_t> data_v(nrecs * bin_schemes_count * nframes);

		for (auto i = 0u; i < nrecs; ++i) {
			ids_v[i] = (vec_recs[i]).getId();
			(vec_recs[i]).getCoords(&extrm_v[i * 4], &extrm_v[i * 4 + 1],
					&extrm_v[i * 4 + 2], &extrm_v[i * 4 + 3],
					&bin_schemes_count, &nframes, &data);
			ull_int offset_j = i * bin_schemes_count * nframes;
			ull_int data_sz = bin_schemes_count * nframes;
			for (auto j = 0u; j < data_sz; ++j) {
				data_v[offset_j + j] = *(data + j);
			}
		}

		size_t start[3], count[3];

		start[0] = nwrittendims_;
		count[0] = nrecs;

		if (checkNCerr(
				nc_put_vara_uint(ncid_, dimIdVID_, start, count,
						ids_v.data()))) {
			ostringstream oss;
			oss << "NC_ERR>> Writing dimension-id data to file: '"
					<< ncfilename_ << "'";
			throw NCDataWriteError(oss.str());

			return 1;
		}

		start[1] = 0;
		count[1] = 4;

		if (checkNCerr(
				nc_put_vara_double(ncid_, dimExtrVID_, start, count,
						extrm_v.data()))) {
			ostringstream oss;
			oss << "NC_ERR>> Writing extremas data to file: '" << ncfilename_
					<< "'";
			throw NCDataWriteError(oss.str());

			return 1;
		}

		start[1] = 0;
		count[1] = bin_schemes_count;

		start[2] = 0;
		count[2] = nframes;

		if (checkNCerr(
				nc_put_vara_uchar(ncid_, coordVID_, start, count,
						data_v.data()))) {
			ostringstream oss;
			oss << "NC_ERR>> Writing coords data to file: '" << ncfilename_
					<< "'";
			throw NCDataWriteError(oss.str());

			return 1;
		}

		NC_close();
	}
	return 0;
}

int Netcdf_TrjInt::addRecords(const std::vector<CoordBAT> &recs) {
	std::vector<CoordBAT>::const_iterator it;
	for (it = recs.begin(); it != recs.end(); ++it) {
		vec_records.push_back(*it);
	}
	return 0;
}

int Netcdf_TrjInt::addRecords(std::vector<CoordBAT> &recs) {
	std::vector<CoordBAT>::iterator it;
	for (it = recs.begin(); it != recs.end(); ++it) {
		vec_records.push_back(*it);
	}
	return 0;
}

int Netcdf_TrjInt::addRecord(CoordBAT rec) {
	vec_records.push_back(rec);
	return 0;
}

int Netcdf_TrjInt::NC_create(std::string const &title) {
	return NC_create(filename_, title);
}

int Netcdf_TrjInt::NC_create(std::string const &filename,
		std::string const &title) {
	if (filename.empty())
		return 1;
	// int dimensionID[NC_MAX_VAR_DIMS];

	//nc_type dataType;
	NCTYPE type = NCTYPE::NC_CENTREHIST;

	// NC_64BIT_DATA NC_NETCDF4
	if (checkNCerr(
			nc_create(filename.c_str(), NC_NOCLOBBER | NC_64BIT_DATA, &ncid_)))
		return 1;

	// Attributes
	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, TI_GA_NC_TITLE, title.size(),
					title.c_str()))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'title' to file: '" << ncfilename_ << "'";
		throw NCAttributeCreateError(oss.str());

		return 1;
	}
	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, TI_GA_NC_APPNAME, 6, "CENTRE"))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'application' to file: '" << ncfilename_
				<< "'";
		throw NCAttributeCreateError(oss.str());

		return 1;
	}
	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, TI_GA_NC_PROG, 8, "Real2Int"))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'program' to file: '" << ncfilename_ << "'";
		throw NCAttributeCreateError(oss.str());
		return 1;
	}
	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, TI_GA_NC_PROGVER, 5, "1.0.0"))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'program version' to file: '" << ncfilename_
				<< "'";
		throw NCAttributeCreateError(oss.str());

		return 1;
	}

	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, TI_GA_NC_CONV, 11,
					"CENTRE_Real2Int"))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'conventions' to file: '" << ncfilename_
				<< "'";
		throw NCAttributeCreateError(oss.str());

		return 1;
	}

	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, TI_GA_NC_CONVVER, 3, "1.0"))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'conventions version' to file: '"
				<< ncfilename_ << "'";
		throw NCAttributeCreateError(oss.str());

		return 1;
	}

	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, TI_GA_NC_CONTENT,
					contents_.size(), contents_.c_str()))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'contents' to file: '" << ncfilename_ << "'";
		throw NCAttributeCreateError(oss.str());

		return 1;
	}

	if (checkNCerr(
			nc_put_att_uint(ncid_, NC_GLOBAL, TI_GA_NC_COORDIMS, NC_UINT, 1,
					&ndims_))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'ndims' to file: '" << ncfilename_ << "'";
		throw NCAttributeCreateError(oss.str());

		return 1;
	}

	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, TI_GA_NC_NBINSCHEMES, NC_UBYTE,
					1, &bin_schemes_count_))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'binning schemes count' to file: '"
				<< ncfilename_ << "'";
		throw NCAttributeCreateError(oss.str());

		return 1;
	}

	if (checkNCerr(
			nc_put_att_uchar(ncid_, NC_GLOBAL, TI_GA_NC_BINSCHEMES, NC_UBYTE,
					bin_schemes_.size(), &bin_schemes_[0]))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'binning schemes' to file: '" << ncfilename_
				<< "'";
		throw NCAttributeCreateError(oss.str());

		return 1;
	}
	// Attribute Definitions ends here

	// Dimension definition starts
	if (checkNCerr(nc_def_dim(ncid_, TI_D_NC_DIMS, NC_UNLIMITED, &dimDID_))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'records' to file: '" << ncfilename_ << "'";
		throw NCDimensionCreateError(oss.str());

		return 1;
	}

	if (checkNCerr(nc_def_dim(ncid_, TI_D_NC_DIMEXTM, 4, &extrDID_))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'extremas' to file: '" << ncfilename_ << "'";
		throw NCDimensionCreateError(oss.str());

		return 1;
	}

	if (checkNCerr(
			nc_def_dim(ncid_, TI_D_NC_BINSCHEMESCOUNT, bin_schemes_count_,
					&binSchemesDID_))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'bin schemes sum' to file: '" << ncfilename_
				<< "'";
		throw NCDimensionCreateError(oss.str());

		return 1;
	}

	if (checkNCerr(nc_def_dim(ncid_, TI_D_NC_FRAME, n_frame_, &frameDID_))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'frame' to file: '" << ncfilename_ << "'";
		throw NCDimensionCreateError(oss.str());

		return 1;
	}

	int dimRecIds[1];
	dimRecIds[0] = dimDID_;

	if (checkNCerr(
			nc_def_var(ncid_, TI_V_NC_DIMID, NC_UINT, 1, &dimRecIds[0],
					&dimIdVID_))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'bin-mid-points' to file: '" << ncfilename_
				<< "'";
		throw NCVariableCreateError(oss.str());

		return 1;
	}

	int dimRecExtrms[2];
	dimRecExtrms[0] = dimDID_;
	dimRecExtrms[1] = extrDID_;

	if (checkNCerr(
			nc_def_var(ncid_, TI_V_NC_DIMEXTM, NC_DOUBLE, 2, &dimRecExtrms[0],
					&dimExtrVID_))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'dim extremas' to file: '" << ncfilename_
				<< "'";
		throw NCVariableCreateError(oss.str());

		return 1;
	}

	int dimsBinFreq[3];
	dimsBinFreq[0] = dimDID_;
	dimsBinFreq[1] = binSchemesDID_;
	dimsBinFreq[2] = frameDID_;

	unsigned long chunkSizes[4];
	chunkSizes[0] = 1;
	chunkSizes[1] = bin_schemes_count_; // bin_schemes_count_ : 1
	chunkSizes[2] = n_frame_;

	if (checkNCerr(
			nc_def_var(ncid_, TI_V_NC_COORD, NC_UBYTE, 3, &dimsBinFreq[0],
					&coordVID_))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'frequencies' to file: '" << ncfilename_
				<< "'";
		throw NCVariableCreateError(oss.str());

		return 1;
	}

	// Set fill mode
	if (checkNCerr(nc_set_fill(ncid_, NC_NOFILL, dimsBinFreq))) {
		ostringstream oss;
		oss << "NC_ERR>> Creating 'fill value' to file: '" << ncfilename_
				<< "'";
		throw NCVariableCreateError(oss.str());

		return 1;
	}

	// End netcdf definitions
	if (checkNCerr(nc_enddef(ncid_))) {
		ostringstream oss;
		oss << "NC_ERR>> Finalizing 'definitions' to file: '" << ncfilename_
				<< "'";
		throw NCVariableCreateError(oss.str());

		return 1;
	}

	NC_close();
	return 0;
}
