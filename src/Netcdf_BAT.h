/* 
 * File:   Netcdf_BAT.h
 * Author: shailesh
 *
 * Created on 19 September, 2016, 3:22 PM
 */

#include "NetcdfFile.h"

// #include "netcdf.h"

#include <string>
#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <sstream>

#include "configcentre.h"
#include "CoordinateInfo.h"
#include "NetcdfFile.h"
#include "CentreStdio.h"
#include "Version.h"

#define BINTRAJ 1

#ifdef BINTRAJ
// DEFINES
#define NCPSEUDOS    "pseudos"
#define NCTUPLE2     "tuple2"
#define NCFRAME      "frame"
#define NCBOND       "bond"
#define NCANGLE      "angle"
#define NCDIHEDRAL   "dihedral"
#define NC_TIME      "time"
#define NC_BND       "bond"
#define NC_ANG       "angle"
#define NC_DIH       "dihedral"
#define NC_TUP2      "tuple2"
#define NC_PSEUDO    "pseudoBonds"

/// The base interface to NetCDF trajectory files.
class Netcdf_BAT: public NetcdfFile {
public:
	Netcdf_BAT() :
			NetcdfFile() {
		ncframe_ = -1;
		timeVID_ = -1;
		frameDID_ = -1;
		pseudoDID_ = -1;
		pseudoVID_ = -1;
		ncid_ = -1;
		dihedralDID_ = -1;
		dihedralVID_ = -1;
		bondDID_ = -1;
		bondVID_ = -1;
		angleDID_ = -1;
		angleVID_ = -1;
		ncdebug_ = -1;
		ncpseudos_ = -1;
		tuple2DID_ = -1;
		pseudoBondsVID_ = -1;
		labelDID_ = -1;
		ncdihedrals_ = -1;
		ncbonds_ = -1;
		ncangles_ = -1;
		nctuple2_ = -1;
	}
	Netcdf_BAT(std::string fname) :
			NetcdfFile(fname) {
		ncframe_ = -1;
		timeVID_ = -1;
		frameDID_ = -1;
		filename_ = fname;
		pseudoDID_ = -1;
		pseudoVID_ = -1;
		ncid_ = -1;
		dihedralDID_ = -1;
		dihedralVID_ = -1;
		bondDID_ = -1;
		bondVID_ = -1;
		angleDID_ = -1;
		angleVID_ = -1;
		ncdebug_ = -1;
		ncpseudos_ = -1;
		tuple2DID_ = -1;
		pseudoBondsVID_ = -1;
		labelDID_ = -1;
		ncdihedrals_ = -1;
		ncbonds_ = -1;
		ncangles_ = -1;
		nctuple2_ = -1;
	}
	int NC_create(std::string const&, NCTYPE, CoordinateInfo&,
			std::string const&);

	int setupFrameDim();
	int SetupCoordsBAT();
	int setupTime();
	int setupRead(CoordinateInfo &c);
	int setupWrite();
	int appendFrames(const std::vector<data_t> &time,
			const std::vector<data_t> &bnd, const std::vector<data_t> &ang,
			const std::vector<data_t> &dih);

	int readFrames(size_t frm_frame, size_t frm_count,
			std::vector<float> &times, std::vector<float> &bonds,
			std::vector<float> &angles, std::vector<float> &diheds);

	int readBondFrames(size_t bnd_id, size_t frm_frame, size_t frm_count,
			std::vector<float> &data);
	int readAngleFrames(size_t ang_id, size_t frm_frame, size_t frm_count,
			std::vector<float> &data);
	int readDihedralFrames(size_t dih_id, size_t frm_frame, size_t frm_count,
			std::vector<float> &data);

	int readMultiBondFrames(size_t from_bnd_id, size_t n_bnds, size_t frm_frame,
			size_t frm_count, std::vector<std::vector<float>> &data);

	int readMultiAngleFrames(size_t from_ang_id, size_t n_angs,
			size_t frm_frame, size_t frm_count,
			std::vector<std::vector<float>> &data);

	int readMultiDihedralFrames(size_t from_dih_id, size_t n_tors,
			size_t frm_frame, size_t frm_count,
			std::vector<std::vector<float>> &data);

	void floatToDouble(double*, const float*, int);
	void doubleToFloat(float*, const double*, int);

	inline int Ncbond() const {
		return ncbonds_;
	}

	inline int Ncangle() const {
		return ncangles_;
	}

	inline int Ncdihedral() const {
		return ncdihedrals_;
	}

	inline int Ncpseudo() const {
		return ncpseudos_;
	}

	inline int Ncframe() const {
		return ncframe_;
	}

	inline int BondVID() const {
		return bondVID_;
	}

	inline int AngleVID() const {
		return angleVID_;
	}

	inline int DihedralVID() const {
		return dihedralVID_;
	}

	bool HasPseudos() {
		return (pseudoBondsVID_ != -1);
	}

	bool HasBonds() {
		return (bondVID_ != -1);
	}

	bool HasAngles() {
		return (angleVID_ != -1);
	}

	bool HasDihedrals() {
		return (dihedralVID_ != -1);
	}

	bool HasTimes() {
		return (timeVID_ != -1);
	}

protected:
	// TODO: Make all private
	size_t start_[4];
	size_t count_[4];
	int ncdebug_;
	int ncframe_;
	int ncbonds_;
	int ncangles_;
	int ncdihedrals_;
	int ncpseudos_;
	int nctuple2_;
	int pseudoVID_;
	int bondVID_;
	int angleVID_;
	int dihedralVID_;
	int timeVID_; ///< Time variable ID.

private:
	std::string filename_;
	int frameDID_;
	int bondDID_;
	int angleDID_;
	int dihedralDID_;
	int tuple2DID_;
	int pseudoDID_;
	int labelDID_;
	int pseudoBondsVID_;
	// #   endif
};

#endif /* NETCDF_BAT_H */

