/* 
 * File:   HistUtils.h
 * Author: shailesh
 *
 * Created on 16 September, 2016, 11:42 PM
 */

#ifndef HISTUTILS_H
#define HISTUTILS_H
// for ostream overload
#include <iostream>

#include <iosfwd>
#include <cstdint>
#include <vector>
#include <string>
#include <cmath>
#include <netcdf.h>

#include "configcentre.h"
#include "NetcdfFile.h"
#include "CentreStdio.h"
#include "Constants.h"
#include "Version.h"

#define HU_GA_NC_TITLE             "title"
#define HU_GA_NC_APPNAME           "appication"
#define HU_GA_NC_PROG              "program"
#define HU_GA_NC_PROGVER           "programVersion"
#define HU_GA_NC_CONV              "convention"
#define HU_GA_NC_CONVVER           "conventionVersion"
#define HU_GA_NC_NDIG              "title"
#define HU_GA_NC_CONTENT           "content"
#define HU_GA_NC_ORDER             "order"
#define HU_GA_NC_STEP_1            "firstStep"
#define HU_GA_NC_NSTEPS            "nsteps"
#define HU_GA_NC_STEPSTRIDE        "stepStride"
#define HU_GA_NC_STORAGETYPE       "storagetype"
#define HU_GA_NC_NDIMS             "ndims"
#define HU_GA_NC_DIMLENGHS         "dimLengths"
#define HU_GA_NC_NBINSCHEMES       "nbinSchemes"
#define HU_GA_NC_BINSCHEMES        "binSchemes"

#define HU_D_NC_RECORD             "records"
#define HU_D_NC_RECORDID           "ids"
#define HU_D_NC_DIMEXTRM           "extremes"
#define HU_D_NC_BINMIDS      	   "dimsXbin_sum"
#define HU_D_NC_BINFREQS      	   "bin_pow_dims_sum"
#define HU_D_NC_NSTEPEFF           "nsteps_eff"

#define HU_V_NC_RECORDID           "id"
#define HU_V_NC_BINMID             "mid"
#define HU_V_NC_BINFREQ            "freq"
#define HU_V_NC_DIMEXTRM           "extrema"

class BinGroup {
	hbin_t nstep_eff_;
	ull_int bin_dims_sum_mid_;
	ull_int bin_pow_dims_sum_;
	std::vector<u_int> record_id;
	std::vector<double> dim_extremes;
	std::vector<double> bin_mid;
	std::vector<ull_int> count;
public:
	BinGroup(const u_int id, const hbin_t nstep_eff,
			const ull_int bin_schemes_sum);

	BinGroup(const std::vector<u_int> &id, const hbin_t nstep_eff,
			const ull_int bin_mids_dim, const ull_int bin_pow_dim_sum);

	friend std::ostream& operator<<(std::ostream &os, const BinGroup &bp);

	void setId(const u_int &in);

	void setId(const ull_int start_index_in, const ull_int in_count,
			const std::vector<u_int> &in);

	// TODO: Make all Get??? members 'const' without adding any performance overhead
	void getId(ull_int *size, u_int **p);

	void getExtremes(ull_int *size, double **p);

	void setExtremes(const std::vector<double> &in);

	void getBinMids(ull_int *size, double **data);

	void getBinFreqs(ull_int *nsteps, ull_int *nbins, ull_int **data);

	int setBinMids(const ull_int start_bin_id, const ull_int bin_count,
			const std::vector<double> &in);

	int setBinMids(const ull_int start_bin_id, const ull_int bin_count,
			const ull_int start_in, const std::vector<double> &in);

	int setBinFreqs(const hbin_t step_id, const ull_int bin_id,
			const ull_int value);

	int setBinFreqs(const hbin_t step_id, const ull_int start_bin_id,
			const ull_int bin_count, const std::vector<ull_int> &in);

	int setBinFreqs(const hbin_t start_step_id, const hbin_t steps_count,
			const ull_int start_bin_id, const ull_int bins_count,
			const std::vector<ull_int> &in);

	int setBinFreqs(const hbin_t start_step_in, const hbin_t step_stride_in,
			const hbin_t steps_cnt_eff_in, const ull_int start_bin_id,
			const ull_int bins_count, const std::vector<ull_int> &in);

	int setBinFreqs(const hbin_t start_step_id, const hbin_t steps_count,
			const ull_int start_bin_id, const ull_int bins_count,
			const ull_int start_in, const ull_int *in);

	int setBinFreqs(const hbin_t start_step_id, const hbin_t steps_count,
			const ull_int start_bin_id, const ull_int bins_count,
			const ull_int start_in, const std::vector<ull_int> &in);

};

class Netcdf_HistUtil: public NetcdfFile {
private:
	std::string filename_;
	std::string contents_;
	hbin_t ndiag_;
	hbin_t ndims_;
	std::vector<ull_int> dim_lens_;
	hbin_t order_;
	hbin_t first_step_;
	hbin_t nsteps_;
	hbin_t step_stride_;
	hbin_t nstep_eff_;
	TensorType storage_type_;
	hbin_t bin_schemes_count_;
	std::vector<hbin_t> bin_schemes_;
	ull_int dimsXbin_schemes_sum_;
	ull_int bin_dims_sum_mid_;
	ull_int bin_pow_dims_sum_;
	ull_int nrecords;
	std::vector<BinGroup> vec_records;

	// variables for storing dim ids
	int recordsDID_;
	int binDimsDID_;
	int binPowDimDID_;
	int nstepsEffDID_;
	int recordIdDID_;
	int extrmDID_;

	// variables for storing variable ids
	int binMidsVID_;
	int binFreqVID_;
	int recordIdVID_;
	int extrmVID_;

public:

	Netcdf_HistUtil();
	Netcdf_HistUtil(const std::string &filename, const std::string &c,
			const hbin_t ndims, const hbin_t order, const hbin_t first_step,
			const hbin_t nsteps, const hbin_t step_stride,
			const std::vector<ull_int> &dim_lens, const TensorType &typ,
			const hbin_t bin_schemes_count,
			const std::vector<hbin_t> &bin_schemes);

	inline hbin_t getFirstStep() {
		return first_step_;
	}

	inline hbin_t getStepStride() {
		return step_stride_;
	}

	inline hbin_t getNSteps() {
		return nsteps_;
	}

	inline hbin_t getStepsEff() {
		return nstep_eff_;
	}

	inline ull_int getDimXBinSchemesSum() {
		return dimsXbin_schemes_sum_;
	}

	friend std::ostream& operator<<(std::ostream &os,
			const Netcdf_HistUtil &ht);

	void clearRecords();

	int removeRecord(const BinGroup &rec);

	int removeRecords(const std::vector<BinGroup> &recs);

	int readRecord(const size_t dim_index, std::vector<double> &extrm_v,
			std::vector<ull_int> &freqs_v);

	int readRecord(const size_t dim_index, std::vector<u_int> &ids_v,
			std::vector<double> &extrm_v, std::vector<ull_int> &freqs_v);

	int readRecord(const size_t dim_index, std::vector<u_int> &ids_v,
			std::vector<double> &extrm_v, std::vector<double> &mids_v,
			std::vector<ull_int> &freqs_v);

	int readRecord(const size_t start_index, BinGroup rd);

	int writeRecord(BinGroup &rec);
	int writeRecords(std::vector<BinGroup> &vec_rec);

	int addRecords(const std::vector<BinGroup> &recs);

	int addRecords(std::vector<BinGroup> &recs);

	int addRecord(const BinGroup rec);

	void printRecords();

	int NC_create(std::string const &title);

private:
	int setupRead();

	int setupWrite();

	int NC_create(std::string const &Name, std::string const &title);
};

#endif /* HISTUTILS_H */
