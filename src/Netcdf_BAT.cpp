#include "Netcdf_BAT.h"

/** Get the frame dimension ID and # of frames (ncframe). */
int Netcdf_BAT::setupFrameDim() {
	frameDID_ = getDimInfo(NCFRAME, &ncframe_);
	if (frameDID_ == -1)
		return 1;
	return 0;
}

/** Determine if Netcdf file contains time; set up timeVID and check units. */
int Netcdf_BAT::setupTime() {
	if (nc_inq_varid(ncid_, NC_TIME, &timeVID_) == NC_NOERR) {
		std::string attrText = getAttrText(timeVID_, "units");
		if (attrText != "picosecond")
			mprintf(
					"Warning: NetCDF file has time units of '%s' - expected picosecond.\n",
					attrText.c_str());
		return 0;
	}
	return 1;
}

int Netcdf_BAT::setupRead(CoordinateInfo &cInfo) {

	if (ncmode_ != NC_NOWRITE) {
		NC_close();
	}
	if (ncid_ == -1) {
		NC_openRead(filename_);
	}

	if ((setupFrameDim() && setupTime()) != 0)
		return -1;

	pseudoDID_ = getDimInfo(NCPSEUDOS, &ncpseudos_);

	pseudoBondsVID_ = -1;

	if (nc_inq_varid(ncid_, NC_PSEUDO, &pseudoBondsVID_) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf file has pseudos.\n");
	} else {
		return -2;
	}

	bondDID_ = getDimInfo(NCBOND, &ncbonds_);
	angleDID_ = getDimInfo(NCANGLE, &ncangles_);
	dihedralDID_ = getDimInfo(NCDIHEDRAL, &ncdihedrals_);

	// Get coord info
	bondVID_ = -1;
	if (nc_inq_varid(ncid_, NC_BND, &bondVID_) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf file has bonds.\n");
		std::string attrText = getAttrText(bondVID_, "units");
		if (attrText != "angstrom")
			mprintf(
					"Warning: Netcdf file has length units of '%s' - expected angstrom.\n",
					attrText.c_str());

	}
	angleVID_ = -1;
	if (nc_inq_varid(ncid_, NC_ANG, &angleVID_) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf file has angles.\n");
		std::string attrText = getAttrText(angleVID_, "units");
		if (!(attrText == "degree" || attrText == "radian"))
			mprintf(
					"Warning: Netcdf file has angle units of '%s' - expected degree/radian.\n",
					attrText.c_str());
	}
	dihedralVID_ = -1;
	if (nc_inq_varid(ncid_, NC_DIH, &dihedralVID_) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf file has dihedral.\n");
		std::string attrText = getAttrText(dihedralVID_, "units");
		if (!(attrText == "degree" || attrText == "radian"))
			mprintf(
					"Warning: Netcdf file has dihedral-angle units of '%s' - expected degree/radian.\n",
					attrText.c_str());
	}

	std::string atmStart = getAttrText(NC_GLOBAL, "atomStartIndex");
	std::istringstream stream1(atmStart);
	long atmStartInd;
	stream1 >> atmStartInd;
	cInfo.SetAtomStartIndex(atmStartInd);
	std::string rootIndices = getAttrText(NC_GLOBAL, "rootIndices");
	std::istringstream stream(rootIndices);

	int r1, r2, r3;
	stream >> r1;
	if (!stream.fail()) {
		stream >> r2;
		if (!stream.fail()) {
			stream >> r3;
			if (!stream.fail()) {
				cInfo.SetRoots(r1, r2, r3);
			}
		}
	}
	std::string hasPseudos = getAttrText(NC_GLOBAL, "hasPseudo");
	if (hasPseudos != "N") {
		cInfo.SetPseudo(true);
		std::vector<size_t> start_(2);
		std::vector<size_t> count_(2);
		start_[0] = 0;
		start_[1] = 0;
		count_[0] = ncpseudos_;
		count_[1] = 2;
		std::vector<long> pseudo(2 * ncpseudos_);
		size_t num_pseudo = 2 * ncpseudos_;
		if (nc_get_vara_long(ncid_, pseudoBondsVID_, &start_[0], &count_[0],
				&pseudo[0]) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file has %d pseudos.\n", num_pseudo/2);
		}
		cInfo.SetPseudo(ncpseudos_);
		cInfo.SetPseudoBonds(pseudo);
	} else {
		cInfo.SetPseudo(false);
	}
	if (timeVID_ != -1) {
		cInfo.SetTime(true);
	}
	cInfo.SetBond(ncbonds_);
	cInfo.SetAngle(ncangles_);
	cInfo.SetDihedral(ncdihedrals_);

	return 0;
}

int Netcdf_BAT::readFrames(size_t frm_frame, size_t frm_count,
		std::vector<float> &times, std::vector<float> &bonds,
		std::vector<float> &angles, std::vector<float> &diheds) {
	std::vector<size_t> start_(2);
	std::vector<size_t> count_(2);
	start_[0] = frm_frame;
	count_[0] = frm_count;

	if (nc_get_vara_float(ncid_, timeVID_, &start_[0], &count_[0],
			times.data()) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf file reading bonds.\n");
	} else {
		return -1;
	}

	start_[1] = 0;
	count_[1] = Ncbond();
	if (nc_get_vara_float(ncid_, bondVID_, &start_[0], &count_[0],
			bonds.data()) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf file reading bonds.\n");
	} else {
		return -1;
	}

	start_[1] = 0;
	count_[1] = Ncangle();
	if (nc_get_vara_float(ncid_, angleVID_, &start_[0], &count_[0],
			angles.data()) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf file reading bonds.\n");
	} else {
		return -1;
	}

	start_[1] = 0;
	count_[1] = Ncdihedral();
	if (nc_get_vara_float(ncid_, dihedralVID_, &start_[0], &count_[0],
			diheds.data()) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf file reading bonds.\n");
	} else {
		return -1;
	}
	return 0;
}

int Netcdf_BAT::readBondFrames(size_t bnd_id, size_t frm_frame,
		size_t frm_count, std::vector<float> &data) {
	std::vector<size_t> start_(2);
	std::vector<size_t> count_(2);
	start_[0] = frm_frame;
	start_[1] = bnd_id;
	count_[0] = frm_count;
	count_[1] = 1;
	if (nc_get_vara_float(ncid_, bondVID_, &start_[0], &count_[0],
			&data[0]) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf file reading bonds.\n");
	} else {
		return -1;
	}
	return 0;
}

int Netcdf_BAT::readAngleFrames(size_t bnd_id, size_t frm_frame,
		size_t frm_count, std::vector<float> &data) {
	std::vector<size_t> start_(2);
	std::vector<size_t> count_(2);
	start_[0] = frm_frame;
	start_[1] = bnd_id;
	count_[0] = frm_count;
	count_[1] = 1;
	if (nc_get_vara_float(ncid_, angleVID_, &start_[0], &count_[0],
			&data[0]) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf file reading angles.\n");
	} else {
		return -1;
	}
	return 0;
}

int Netcdf_BAT::readDihedralFrames(size_t bnd_id, size_t frm_frame,
		size_t frm_count, std::vector<float> &data) {
	std::vector<size_t> start_(2);
	std::vector<size_t> count_(2);
	start_[0] = frm_frame;
	start_[1] = bnd_id;
	count_[0] = frm_count;
	count_[1] = 1;
	if (nc_get_vara_float(ncid_, dihedralVID_, &start_[0], &count_[0],
			&data[0]) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf file reading dihedrals.\n");
	} else {
		return -1;
	}
	return 0;
}

int Netcdf_BAT::readMultiBondFrames(size_t from_bnd_id, size_t n_bnds,
		size_t frm_frame, size_t frm_count,
		std::vector<std::vector<float>> &data) {
	if (data.size()) {
		std::vector<float> tmp(data[0].size() * data.size());
		std::vector<size_t> start_(2);
		std::vector<size_t> count_(2);
		start_[0] = frm_frame;
		start_[1] = from_bnd_id;
		count_[0] = frm_count;
		count_[1] = n_bnds;
		if (nc_get_vara_float(ncid_, bondVID_, &start_[0], &count_[0],
				tmp.data()) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file reading bonds.\n");
		} else {
			return -1;
		}
		// transpose tmp to get data, every dim in inner vector (column wise)
		size_t tt = 0;
		for (auto ti = 0U; ti < data[0].size(); ++ti) {
			for (auto tj = 0U; tj < data.size(); ++tj) {
				data[tj][ti] = tmp[tt];
				++tt;
			}
		}
	} else {
		return -2;
	}
	return 0;
}

int Netcdf_BAT::readMultiAngleFrames(size_t from_ang_id, size_t n_angs,
		size_t frm_frame, size_t frm_count,
		std::vector<std::vector<float>> &data) {
	if (data.size()) {
		std::vector<float> tmp(data[0].size() * data.size());
		std::vector<size_t> start_(2);
		std::vector<size_t> count_(2);
		start_[0] = frm_frame;
		start_[1] = from_ang_id;
		count_[0] = frm_count;
		count_[1] = n_angs;
		if (nc_get_vara_float(ncid_, angleVID_, &start_[0], &count_[0],
				tmp.data()) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file reading angles.\n");
		} else {
			return -1;
		}
		// transpose tmp to get data, every dim in inner vector (column wise)
		size_t tt = 0;
		for (auto ti = 0U; ti < data[0].size(); ++ti) {
			for (auto tj = 0U; tj < data.size(); ++tj) {
				data[tj][ti] = tmp[tt];
				++tt;
			}
		}
	} else {
		return -2;
	}
	return 0;
}

int Netcdf_BAT::readMultiDihedralFrames(size_t from_dih_id, size_t n_tors,
		size_t frm_frame, size_t frm_count,
		std::vector<std::vector<float>> &data) {
	if (data.size()) {
		std::vector<float> tmp(data[0].size() * data.size());
		std::vector<size_t> start_(2);
		std::vector<size_t> count_(2);
		start_[0] = frm_frame;
		start_[1] = from_dih_id;
		count_[0] = frm_count;
		count_[1] = n_tors;
		if (nc_get_vara_float(ncid_, dihedralVID_, &start_[0], &count_[0],
				tmp.data()) == NC_NOERR) {
			if (ncdebug_ > 0)
				mprintf("\tNetcdf file reading dihedrals.\n");
		} else {
			return -1;
		}
// transpose tmp to get data, every dim in inner vector (column wise)
		size_t tt = 0;
		for (auto ti = 0U; ti < data[0].size(); ++ti) {
			for (auto tj = 0U; tj < data.size(); ++tj) {
				data[tj][ti] = tmp[tt];
				++tt;
			}
		}
	} else {
		return -2;
	}
	return 0;
}

int Netcdf_BAT::setupWrite() {
	if (ncmode_ != NC_WRITE && ncid_ != -1) {
		NC_close();
	}
	if (ncid_ == -1) {
		NC_openWrite(filename_);
	}
	//int temp;
	tuple2DID_ = getDimInfo(NCTUPLE2, &nctuple2_);
	pseudoDID_ = getDimInfo(NCPSEUDOS, &ncpseudos_);
	frameDID_ = getDimInfo(NCFRAME, &ncframe_);
	bondDID_ = getDimInfo(NCBOND, &ncbonds_);
	angleDID_ = getDimInfo(NCANGLE, &ncangles_);
	dihedralDID_ = getDimInfo(NCDIHEDRAL, &ncdihedrals_);

	pseudoVID_ = -1;
	if (nc_inq_varid(ncid_, NC_PSEUDO, &pseudoVID_) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf error reading pseudo variable.\n");
	}

	// Get coord info
	timeVID_ = -1;
	if (nc_inq_varid(ncid_, NC_TIME, &timeVID_) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf error reading time variable.\n");
	}

	bondVID_ = -1;
	if (nc_inq_varid(ncid_, NC_BND, &bondVID_) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf error reading bond variable.\n");
	}

	angleVID_ = -1;
	if (nc_inq_varid(ncid_, NC_ANG, &angleVID_) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf error reading angle variable.\n");
	}

	dihedralVID_ = -1;
	if (nc_inq_varid(ncid_, NC_DIH, &dihedralVID_) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf error reading torsion variable.\n");
	}

	return 0;
}

int Netcdf_BAT::appendFrames(const std::vector<data_t> &time,
		const std::vector<data_t> &bnd, const std::vector<data_t> &ang,
		const std::vector<data_t> &dih) {
	if (ncid_ != -1) {
		NC_close();
	}
	if (NC_openWrite(filename_) != 0) {
		mprinterr("Error: Opening file for writing.\n");
		return -1;
	}
//  std::cout << "Writing record..." << std::endl;
	setupWrite();
	//ull_int nframes;

	size_t start[2], count[2];

	start[0] = ncframe_;
	count[0] = time.size();

	if (checkNCerr(
			nc_put_vara_float(ncid_, timeVID_, start, count, time.data()))) {
		mprinterr("Error: Writing time variable data.\n");
		return 1;
	}

	start[1] = 0;
	count[1] = Ncbond();
	if (checkNCerr(
			nc_put_vara_float(ncid_, bondVID_, start, count, bnd.data()))) {
		mprinterr("Error: Writing bond variable data.\n");
		return 1;
	}

	start[1] = 0;
	count[1] = Ncangle();
	if (checkNCerr(
			nc_put_vara_float(ncid_, angleVID_, start, count, ang.data()))) {
		mprinterr("Error: Writing angle variable data.\n");
		return 1;
	}

	start[1] = 0;
	count[1] = Ncdihedral();
	if (checkNCerr(
			nc_put_vara_float(ncid_, dihedralVID_, start, count, dih.data()))) {
		mprinterr("Error: Writing dihedral variable data.\n");
		return 1;
	}
	NC_close();
	return 0;
}

int Netcdf_BAT::NC_create(std::string const &Name, NCTYPE type,
		CoordinateInfo &coordInfo, std::string const &title) {
	if (Name.empty())
		return 1;
	int dimensionID[NC_MAX_VAR_DIMS];
	// chunksizes=(4096, 4,)

	unsigned long chunkSizes[4];
	chunkSizes[0] = 4096;
	chunkSizes[1] = 4;

	//int NDIM;
	int numpseduos;
	nc_type dataType;

	if (checkNCerr(nc_create(Name.c_str(), NC_NOCLOBBER | NC_NETCDF4, &ncid_)))
		return 1;

	// Set number of dimensions based on file type
	switch (type) {
	case NC_CENTRETRAJ:
		//NDIM = 3;
		dataType = NC_FLOAT;
		break;
	default:
		mprinterr("Error: NC_create (%s): Unrecognized type (%i)\n",
				Name.c_str(), (int) type);
		return 1;
	}

	ncframe_ = 0;
	if (type == NC_CENTRETRAJ) {
		if (checkNCerr(nc_def_dim(ncid_, NCTUPLE2, 2, &tuple2DID_))) {
			mprinterr("Error: Defining dimension tuple2.\n");
			return 1;
		}

		numpseduos = coordInfo.nPseudos() > 1 ? coordInfo.nPseudos() : 1;
		if (checkNCerr(nc_def_dim(ncid_, NCPSEUDOS, numpseduos, &pseudoDID_))) {
			mprinterr("Error: Defining dimension pseudo.\n");
			return 1;
		}
		// Frame dimension for traj
		if (checkNCerr(nc_def_dim(ncid_, NCFRAME, NC_UNLIMITED, &frameDID_))) {
			mprinterr("Error: Defining dimension frame.\n");
			return 1;
		}
		// Since frame is UNLIMITED, it must be lowest dim.
		dimensionID[0] = frameDID_;
		if (checkNCerr(
				nc_def_dim(ncid_, NC_BND, coordInfo.nBonds(), &bondDID_))) {
			mprinterr("Error: Defining dimension bonds.\n");
			return 1;
		}
		if (checkNCerr(
				nc_def_dim(ncid_, NC_ANG, coordInfo.nAngles(), &angleDID_))) {
			mprinterr("Error: Defining dimension angles.\n");
			return 1;
		}
		if (checkNCerr(
				nc_def_dim(ncid_, NC_DIH, coordInfo.nDihedrals(),
						&dihedralDID_))) {
			mprinterr("Error: Defining dimension dihedrals.\n");
			return 1;
		}
	}

	dimensionID[0] = pseudoDID_;
	dimensionID[1] = tuple2DID_;
	if (checkNCerr(
			nc_def_var(ncid_, NC_PSEUDO, NC_INT, 2, dimensionID,
					&pseudoBondsVID_))) {
		mprinterr("Error: Defining variable pseudoBonds.\n");
		return 1;
	}

	// Time variable and units
	if (coordInfo.HasTime()) {
		dimensionID[0] = frameDID_;
		if (checkNCerr(
				nc_def_var(ncid_, NC_TIME, dataType, 1, dimensionID,
						&timeVID_))) {
			mprinterr("Error: Defining variable time.\n");
			return 1;
		}
		if (checkNCerr(
				nc_put_att_text(ncid_, timeVID_, "units", 10, "picosecond"))) {
			mprinterr("Error: Defining unit attribute for variable time.\n");
			return 1;
		}
	}

	// Bond variable
	if (coordInfo.HasAngle()) {
		dimensionID[0] = frameDID_;
		dimensionID[1] = bondDID_;
		if (checkNCerr(
				nc_def_var(ncid_, NC_BND, dataType, 2, dimensionID,
						&bondVID_))) {
			mprinterr("Error: Defining varible bonds.\n");
			return 1;
		}
		if (checkNCerr(
				nc_def_var_chunking(ncid_, bondVID_, NC_CHUNKED,
						&chunkSizes[0]))) {
			mprinterr("Error: Writing bin frequencies.\n");
			return 1;
		}
		if (checkNCerr(
				nc_put_att_text(ncid_, bondVID_, "units", 8, "angstrom"))) {
			mprinterr("Error: Defining unit attribute for variable bonds.\n");
			return 1;
		}
	}
	// Angle variable

	if (coordInfo.HasAngle()) {
		dimensionID[0] = frameDID_;
		dimensionID[1] = angleDID_;
		if (checkNCerr(
				nc_def_var(ncid_, NC_ANG, dataType, 2, dimensionID,
						&angleVID_))) {
			mprinterr("Error: Defining variable angles.\n");
			return 1;
		}
		if (checkNCerr(
				nc_def_var_chunking(ncid_, bondVID_, NC_CHUNKED,
						&chunkSizes[0]))) {
			mprinterr("Error: Writing bin frequencies.\n");
			return 1;
		}
		if (checkNCerr(
				nc_put_att_text(ncid_, angleVID_, "units", 6, "degree"))) {
			mprinterr("Error: Defining unit attribute for variable angle.\n");
			return 1;
		}
	}
	// Dihedral variable
	if (coordInfo.HasDihedral()) {
		dimensionID[0] = frameDID_;
		dimensionID[1] = dihedralDID_;
		if (checkNCerr(
				nc_def_var(ncid_, NC_DIH, dataType, 2, dimensionID,
						&dihedralVID_))) {
			mprinterr("Error: Defining varible dihedral\n");
			return 1;
		}
		if (checkNCerr(
				nc_def_var_chunking(ncid_, bondVID_, NC_CHUNKED,
						&chunkSizes[0]))) {
			mprinterr("Error: Writing bin frequencies.\n");
			return 1;
		}
		if (checkNCerr(
				nc_put_att_text(ncid_, dihedralVID_, "units", 6, "degree"))) {
			mprinterr(
					"Error: Defining unit attribute for variable dihedral.\n");
			return 1;
		}
	}

	// Attributes
	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, "title", title.size(),
					title.c_str()))) {
		mprinterr("Error: Writing attribute title.\n");
		return 1;
	}
	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, "application", 6, "CENTRE"))) {
		mprinterr("Error: Writing attribute application.\n");
		return 1;
	}
	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, "program", 8, "Cart2BAT"))) {
		mprinterr("Error: Writing attribute program.\n");
		return 1;
	}
	if (checkNCerr(nc_put_att_text(ncid_, NC_GLOBAL, "programVersion",
	NETCDF_VERSION_STRLEN,
	NETCDF_VERSION_STRING))) {
		mprinterr("Error: Writing attribute program version.\n");
		return 1;
	}
	// TODO: Make conventions a static string
	bool errOccurred = false;
	if (type == NC_CENTRETRAJ)
		errOccurred = checkNCerr(
				nc_put_att_text(ncid_, NC_GLOBAL, "conventions", 6, "CENTRE"));
	if (errOccurred) {
		mprinterr("Error: Writing conventions.\n");
		return 1;
	}

	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, "conventionVersion", 3, "1.0"))) {
		mprinterr("Error: Writing conventions version.\n");
		return 1;
	}

	std::string tmp = std::to_string(coordInfo.atomStartIdx());
	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, "atomStartIndex", tmp.length(),
					tmp.c_str()))) {
		mprinterr("Error: Writing conventions version.\n");
		return 1;
	}

	std::ostringstream ss1;
	int r1, r2, r3;
	coordInfo.rootAtoms(r1, r2, r3);
	ss1 << r1 << " " << r2 << " " << r3;
	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, "rootIndices", ss1.str().length(),
					ss1.str().c_str()))) {
		mprinterr("Error: Writing conventions version.\n");
		return 1;
	}

	std::string tmp2 = coordInfo.HasPseudo() ? "Y" : "N";
	if (checkNCerr(
			nc_put_att_text(ncid_, NC_GLOBAL, "hasPseudo", tmp2.length(),
					tmp2.c_str()))) {
		mprinterr("Error: Writing conventions version.\n");
		return 1;
	}

	// Set fill mode
	if (checkNCerr(nc_set_fill(ncid_, NC_NOFILL, dimensionID))) {
		mprinterr("Error: NetCDF setting fill value.\n");
		return 1;
	}

	// End netcdf definitions
	if (checkNCerr(nc_enddef(ncid_))) {
		mprinterr("Error: NetCDF error on ending definitions.\n");
		return 1;
	}

	std::vector<long> pseudo_v;
	if (coordInfo.HasPseudo()) {
		coordInfo.getPseudoBonds(pseudo_v);
	} else {
		pseudo_v.push_back(-1);
		pseudo_v.push_back(-1);
	}
	if (ncid_ != -1) {
		NC_close();
	}
	if (NC_openWrite(filename_) != 0) {
		mprinterr("Error: Opening file for writing.\n");
		return -1;
	}

	tuple2DID_ = getDimInfo(NCTUPLE2, &nctuple2_);
	pseudoDID_ = getDimInfo(NCPSEUDOS, &ncpseudos_);

	pseudoVID_ = -1;
	if (nc_inq_varid(ncid_, NC_PSEUDO, &pseudoVID_) == NC_NOERR) {
		if (ncdebug_ > 0)
			mprintf("\tNetcdf file has dimids.\n");
	}
	size_t start[2], count[2];
	start[0] = 0;
	count[0] = 2;
	start[1] = 0;
	count[1] = numpseduos;
	if (checkNCerr(
			nc_put_vara_long(ncid_, pseudoVID_, start, count,
					pseudo_v.data()))) {
		mprinterr("Error: Writing pseudoBond data.\n");
		return 1;
	}
	if (ncid_ != -1) {
		NC_close();
	}
	return 0;
}
