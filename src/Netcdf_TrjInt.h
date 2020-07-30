/* 
 * File:   Netcdf_TrjInt.h
 * Author: shailesh
 *
 * Created on 21 September, 2016, 11:19 AM
 * This class is used for input/oputput of histogram bin trajectory in NetCDF
 * file format.
 * This class is designed to handle all input output of bond or angle or dihedral
 * bin membership represented in integer, it also defines a file format for the
 * purpose and uses NetCDF file format and library as requisite.
 *
 * Dependencies: netcdf-c
 */

#ifndef NETCDF_TRJINT_H
#define NETCDF_TRJINT_H

// for ostream overload
#include <iostream>

#include <iosfwd>
#include <cstdint>
#include <vector>
#include <string>
#include <netcdf.h>

#include "configcentre.h"
#include "NetcdfFile.h"
#include "CentreStdio.h"
#include "Constants.h"
#include "Version.h"

#define TI_GA_NC_TITLE             "title"
#define TI_GA_NC_APPNAME           "application"
#define TI_GA_NC_PROG              "program"
#define TI_GA_NC_PROGVER           "programVersion"
#define TI_GA_NC_CONV              "convention"
#define TI_GA_NC_CONVVER           "conventionVersion"
#define TI_GA_NC_CONTENT           "content"
#define TI_GA_NC_COORDIMS          "ncoordims"
#define TI_GA_NC_NBINSCHEMES       "nbinSchemes"
#define TI_GA_NC_BINSCHEMES        "binSchemes"

#define TI_D_NC_DIMS               "dims"
#define TI_D_NC_DIMEXTM            "extremes"
#define TI_D_NC_FRAME              "frames"
#define TI_D_NC_BINSCHEMESCOUNT    "bin_scheme_count"

#define TI_V_NC_DIMID              "dimid"
#define TI_V_NC_DIMEXTM            "extreme"
#define TI_V_NC_COORD              "coord"

/**
 *  @brief Class that represents a single dimension of one of {bonds, angles, dihedral}.
 */
class CoordBAT {
protected:
	u_int dimid_; /*!< an unsigned integer for dimension Id of the coordinate */
	hbin_t bin_schemes_count_; /*!< number of bins schemes used for discretization of dimension  */
	ull_int n_frame_; /*!< number of frames in dataset used for entropy estimation */
	std::vector<double> extremas_; /*!< variable for min, max, avg, and phase of dimension */
	std::vector<hbin_t> coord_; /*!< integer representing bin membership for dimension for all frames */
public:
	/*! \fn constructor
	 *  \brief A constructor taking id and bin_scheme_count to initialize the object.
	 *  \param id[in] an unsigned int dimension id.
	 *  \param bin_schemes_count[in] a unsigned int > 0 for bin_schemes_count.
	 */
	CoordBAT(const u_int id, const hbin_t bin_schemes_count);

	/*! \fn constructor
	 *  \brief A constructor taking id and bin_scheme_count to initialize the object.
	 *  \param id[in] an unsigned int dimension id.
	 *  \param bin_schemes_count[in] a unsigned int > 0 for bin_schemes_count.
	 *  \param nframes[in] a unsigned int > 0 for number of frames.
	 */
	CoordBAT(const u_int id, const hbin_t bin_schemes_count,
			const ull_int nframes);

	// Output the coordinate set to output stream
	friend std::ostream& operator<<(std::ostream &os, const CoordBAT &bp);

	// Set DimensionID for the coordinate set
	void setId(const u_int in);

	// Get DimensionID for the coordinate set
	u_int getId() const;

	/*! \fn const char *Fn_Test::member(char c,int n)
	 *  \brief A member function for getting coordinate values to the object.
	 *  \param *min_v[out]             a pointer to double to hold minimum value of the dimension.
	 *  \param *max_v[out]             a pointer to double to hold maximum value of the dimension.
	 *  \param *avg_v[out]             a pointer to double to hold average value of the dimension.
	 *  \param *phase_v[out]           a pointer to double to hold phase value of the dimension.
	 *  \param *bin_schemes_count[out] an integer > 0 for count of bin schemes of coord set.
	 *  \param *nframes[out]           an unsigned integer for number of frame in coord set.
	 *  \param **data[out]             a pointer to hbin_t vector of size (bin_schemes_count x nframes)
	 *                                 holding bin membership of frames for given schemes and frames of
	 *                                 the dimension
	 *                                 |_schm-1..nfrm elm/frm stride_|...|schm-k..nfrm elm/frm stride_|.
	 */
	// TODO: Make all Get??? members 'const' without adding any performance overhead
	void getCoords(double *min_v, double *max_v, double *avg_v, double *phase_v,
			hbin_t *bin_schemes_count, ull_int *nframes, hbin_t **data);

	/*! \fn const char *Fn_Test::member(char c,int n)
	 *  \brief A member function for setting coordinate values to the object.
	 *  \param min_v[in]            a double for minimum value of the dimension.
	 *  \param max_v[in]            a double for maximum value of the dimension.
	 *  \param avg_v[in]            a double for average value of the dimension.
	 *  \param phase_v[in]          a double for phase value of the dimension.
	 *  \param start_binscheme[in]  an integer 0-based index for start bin_scheme to set coord.
	 *  \param binschemes_count[in] an integer for count of bin schemes to set coord.
	 *  \param start_frame[in]      an unsigned integer 0-based index of start frame to set coord.
	 *  \param start_end[in]        an unsigned integer 0-based index of end frame to set coord.
	 *  \param in[in]               a vector of size () of integers holding bin membership of frames for
	 *                              given schemes and frames the dimension.
	 *                              |_schm-1..nfrm elm/frm stride_|...|schm-k..nfrm elm/frm stride_|.
	 *  \return                     an integer 0: success, -ve: failure.
	 */
	int setCoords(const double min_v, const double max_v, const double avg_v,
			const double phase_v, const ull_int start_binscheme,
			const ull_int binschemes_count, const ull_int start_frame,
			const ull_int frame_end, const std::vector<hbin_t> &in);

	/*! \fn const char *Fn_Test::member(char c,int n)
	 *  \brief A member function for setting coordinate values to the object.
	 *  \param min_v[in]            a double for minimum value of the dimension.
	 *  \param max_v[in]            a double for maximum value of the dimension.
	 *  \param avg_v[in]            a double for average value of the dimension.
	 *  \param phase_v[in]          a double for phase value of the dimension.
	 *  \param start_binscheme[in]  an integer 0-based index for start bin_scheme to set coord.
	 *  \param binschemes_count[in] an integer for count of bin schemes to set coord.
	 *  \param start_frame[in]      an unsigned integer 0-based index of start frame to set coord.
	 *  \param frame_stride[in]     an unsigned integer > 0, setting every frame_stride th frame
	 *                              to set coord.
	 *  \param start_end[in]        an unsigned integer 0-based index of end frame to set coord.
	 *  \param in[in]               a vector of size () of integers holding bin membership of frames for
	 *                              given schemes and frames the dimension.
	 *                              |_schm-1..nfrm elm/frm stride_|...|schm-k..nfrm elm/frm stride_|.
	 *  \return                     an integer 0: success, -ve: failure.
	 */
	int setCoords(const double min_v, const double max_v, const double avg_v,
			const double phase_v, const ull_int start_binscheme,
			const ull_int binschemes_count, const ull_int start_frame,
			const ull_int frame_stride, const ull_int frame_end,
			const std::vector<hbin_t> &in);

};

class Netcdf_TrjInt: public NetcdfFile {
private:
	std::string filename_;
	std::string contents_;
	u_int ndims_;
	ull_int start_frame_;
	ull_int frame_stride_;
	ull_int n_frame_;
	ull_int n_frame_eff_;
	hbin_t bin_schemes_count_;
	std::vector<hbin_t> bin_schemes_;
	ull_int nrecords, nwrittendims_;

	std::vector<CoordBAT> vec_records;

	// variables for storing dim ids
	int dimDID_;
	int extrDID_;
	int binSchemesDID_;
	int frameDID_;

	// variables for storing variable ids
	int dimIdVID_;
	int dimExtrVID_;
	int coordVID_;

public:

	Netcdf_TrjInt();
	Netcdf_TrjInt(const std::string &filename, const std::string &c,
			const u_int ndims, const ull_int start_frame,
			const ull_int frame_stride, const ull_int nframes,
			const hbin_t bin_schemes_count,
			const std::vector<hbin_t> &bin_schemes);

	friend std::ostream& operator<<(std::ostream &os, const Netcdf_TrjInt &ht);

	// Clear all the record_vec<CoordBAT> members of the trajectory object
	void clearRecords();

	// Delete a record from record_vec<CoordBAT> members of the trajectory object
	int removeRecord(const CoordBAT &rec);

	// Delete record from record_vec<CoordBAT> members of the trajectory object, supplied in input vector
	int removeRecords(const std::vector<CoordBAT> &recs);

	int readCoords(const u_int dim_idx, const hbin_t bin_id, CoordBAT &rd);

	int readCoords(const u_int dim_idx, const hbin_t binscheme_start_id,
			const hbin_t binscheme_count, CoordBAT &rd);

	int readRecords(const u_int dim_idx, const hbin_t bin_id, bool is_append);

	int writeRecord(CoordBAT &rec);

	int writeRecords(std::vector<CoordBAT> &recs, ull_int n_recs = 0);

	// int WriteRecords();

	int addRecords(const std::vector<CoordBAT> &recs);

	int addRecords(std::vector<CoordBAT> &recs);

	int addRecord(const CoordBAT rec);

	/* Setup object for the reading CoordBAT record/records from the corresponding
	 * NetCDF file.
	 */
	int setupRead();

	/* Setup object for the writing CoordBAT record/records to the corresponding
	 * NetCDF file.
	 */
	int setupWrite();

	/* Create a NetCDF file for handing I/O of CoordBAT records with the filename
	 * as set with the constructor.
	 */
	int NC_create(std::string const &title);

private:
	int NC_create(std::string const &Name, std::string const &title);
};

#endif /* NETCDF_TRJINT_H */

