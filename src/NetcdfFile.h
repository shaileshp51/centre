#ifndef INC_NETCDFFILE_H
#define INC_NETCDFFILE_H

#include <netcdf.h>

#include <string>
#include <cstddef>
#include "CoordinateInfo.h"
/// The base interface to NetCDF trajectory files.

class NetcdfFile {
public:

	/// For determining netcdf file type
	enum NCTYPE {
		NC_UNKNOWN = 0, NC_CENTRETRAJ, NC_CENTREHIST, NC_CENTREENTCONTRI
	};
	NCTYPE GetNetcdfConventions(const char*);

	NetcdfFile() :
			ncid_(-1), ncdebug_(-1), ncmode_(-1) {
		ncfilename_ = "";
	}

	std::string getAttrText(const char*);
	NCTYPE getNetcdfConventions();
	int NC_openRead(std::string&);
	int NC_openWrite(std::string&);

	void NC_close();

	inline int Ncid() const {
		return ncid_;
	}
	void writeVIDs() const;
protected:
	// TODO: Make all private
	int ncdebug_;
	int ncid_;
	int ncmode_;
	std::string ncfilename_;

	NetcdfFile(std::string filename) :
			ncid_(-1), ncdebug_(-1), ncmode_(-1), ncfilename_(filename) {
	}

	bool checkNCerr(int);
	std::string getAttrText(int, const char*);
	int getDimInfo(const char*, int*);
	void sync();
};
#endif
