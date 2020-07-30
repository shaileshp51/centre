#include "NetcdfFile.h"
//#include "CentreStdio.h"
#include "Constants.h"
#include "Exceptions.h"
#include "Version.h"

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sstream>

using namespace std;

#define NCTRAJ 1

// NetcdfFile::GetNetcdfConventions()

NetcdfFile::NCTYPE NetcdfFile::GetNetcdfConventions(const char *fname) {
	NCTYPE nctype = NC_UNKNOWN;
	// NOTE: Do not use checkNCerr so this fails silently. Allows routine to
	//       be used in file autodetection.
	if (nc_open(fname, NC_NOWRITE, &ncid_) != NC_NOERR)
		return NC_UNKNOWN;
	nctype = getNetcdfConventions();
	NC_close();

	return nctype;
}

#ifdef NCTRAJ
// DEFINES
#define NCPSEUDOS  "pseudos"
#define NCTUPLE2   "tuple2"
#define NCFRAME    "frame"
#define NCBOND     "bond"
#define NCANGLE    "angle"
#define NCDIHEDRAL "dihedral"
#define NC_TIME    "time"
#define NC_BND     "bond"
#define NC_ANG     "angle"
#define NC_DIH     "dihedral"
#define NC_TUP2    "tuple2"
#define NC_PSEUDO  "pseudoBonds"

/** Get the information about a netcdf attribute with given vid and 
 * attribute text.
 * Since there is no guarantee that null char at the end of retrieved string
 * append one.
 */
std::string NetcdfFile::getAttrText(int vid, const char *attribute) {
	size_t attlen;
	std::string attrOut;
	// Get attr length
	if (checkNCerr(nc_inq_attlen(ncid_, vid, attribute, &attlen))) {
		ostringstream oss;
		oss << "NC_ERR>> Getting length for attribute '" << attribute
				<< "' from file: '" << ncfilename_ << "'";
		throw NCAttibuteTextError(oss.str());
		return attrOut;
	}
	// Allocate space for attr text, plus one for null char
	char *attrText = new char[(attlen + 1)];
	// Get attr text
	if (checkNCerr(nc_get_att_text(ncid_, vid, attribute, attrText))) {
		ostringstream oss;
		oss << "NC_ERR>> Getting attribute text for '" << attribute
				<< "' from file: '" << ncfilename_ << "'";
		throw NCAttibuteTextError(oss.str());
		delete[] attrText;
		return attrOut;
	}
	// Append null char - NECESSARY?
	attrText[attlen] = '\0';
	attrOut.assign(attrText);
	delete[] attrText;

	return attrOut;
}

// NetcdfFile::GetAttrText()

/** Get information about a netcdf global attribute. */
std::string NetcdfFile::getAttrText(const char *attribute) {
	return getAttrText(NC_GLOBAL, attribute);
}

// NetcdfFile::GetNetcdfConventions()

NetcdfFile::NCTYPE NetcdfFile::getNetcdfConventions() {
	NCTYPE nctype = NCTYPE::NC_UNKNOWN;
	std::string attrText = getAttrText(NC_GLOBAL, "convention");
	if (attrText == "CENTRE")
		nctype = NCTYPE::NC_CENTRETRAJ;
	else if (attrText == "CENTRE_HIST")
		nctype = NCTYPE::NC_CENTREHIST;
	else if (attrText == "CENTRE_ENT_CONTRI")
		nctype = NCTYPE::NC_CENTREENTCONTRI;
	else if (attrText.empty()) {
		ostringstream oss;
		oss << "NC_ERR>> Could not get conventions from file: '" << ncfilename_
				<< "'";
		throw NCConventionReadError(oss.str());

	} else {
		ostringstream oss;
		oss << "NC_ERR>> Unrecognized convention '" << attrText
				<< "' from file: '" << ncfilename_ << "'" << endl;
		oss << "Expected convention 'CENTRE'" << endl;
		throw NCUnknownConventionError(oss.str());
	}
	return nctype;
}

/** Return the dimension ID of a given attribute in netcdf file ncid.
 * Also set dimension length.
 */
int NetcdfFile::getDimInfo(const char *attribute, int *length) {
	int dimID;
	size_t slength = 0;

	*length = 0;
	// Get dimid 
	if (checkNCerr(nc_inq_dimid(ncid_, attribute, &dimID))) {
		ostringstream oss;
		oss << "NC_ERR>> Getting DimensionID for attribute '" << attribute
				<< "' from file: '" << ncfilename_ << "'";
		throw NCDimensionIDReadError(oss.str());
		return -1;
	}
	// get Dim length 
	if (checkNCerr(nc_inq_dimlen(ncid_, dimID, &slength))) {
		ostringstream oss;
		oss << "NC_ERR>> Getting Dimension length for attribute '" << attribute
				<< "' from file: '" << ncfilename_ << "'";
		throw NCDimensionLengthReadError(oss.str());
		return -1;
	}
	*length = (int) slength;
	return dimID;
}

// NetcdfFile::checkNCerr()
bool NetcdfFile::checkNCerr(int ncerr) {
	if (ncerr != NC_NOERR) {
		ostringstream oss;
		oss << "NC_ERR>> NETCDF GENERIC ERROR: '" << nc_strerror(ncerr) << "'"
				<< endl;
		throw NCGenericLibraryError(oss.str());
		return true;
	}
	return false;
}

int NetcdfFile::NC_openRead(std::string &Name) {
	if (Name.empty()) {
		ostringstream oss;
		oss << "Empty filename '" << Name << "'" << endl;
		throw NCOpenForReadError(oss.str());
		return 1;
	}
	if (checkNCerr(nc_open(Name.c_str(), NC_NOWRITE, &ncid_))) {
		ostringstream oss;
		oss << "NC_ERR>> Opening in READMODE filename '" << Name << "'" << endl;
		throw NCOpenForReadError(oss.str());
		return 1;
	}
	ncmode_ = NC_NOWRITE;
	return 0;
}

// NetcdfFile::NC_close()

void NetcdfFile::NC_close() {
	if (ncid_ == -1)
		return;
	bool err = checkNCerr(nc_close(ncid_));
	if (ncdebug_ > 0 && !err)
		cout << "Successfully closed ncid " << ncid_ << endl;
	ncid_ = -1;
	ncmode_ = -1;
}

int NetcdfFile::NC_openWrite(std::string &Name) {
	if (Name.empty())
		return 1;
	if (checkNCerr(nc_open(Name.c_str(), NC_WRITE, &ncid_))) {
		ostringstream oss;
		oss << "Opening in NC_WRITE filename '" << Name << "'" << endl;
		throw NCOpenForWriteError(oss.str());
		return 1;
	}
	ncmode_ = NC_WRITE;
	return 0;
}

#endif
