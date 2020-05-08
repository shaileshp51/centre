/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Netcdf_EntConti.h
 * Author: shailesh
 *
 * Created on 23 September, 2016, 7:25 PM
 */

#ifndef NETCDF_ENTCONTI_H
#define NETCDF_ENTCONTI_H

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

#define EC_GA_NC_TITLE             "title"
#define EC_GA_NC_APPNAME           "application"
#define EC_GA_NC_PROG              "program"
#define EC_GA_NC_PROGVER           "programVersion"
#define EC_GA_NC_CONV              "convention"
#define EC_GA_NC_CONVVER           "conventionVersion"
#define EC_GA_NC_NDIG              "title"
#define EC_GA_NC_CONTENT           "content"
#define EC_GA_NC_ORDER             "order"
#define EC_GA_NC_NSTEPS            "nsteps"
#define EC_GA_NC_STORAGETYPE       "storagetype"
#define EC_GA_NC_NDIMS             "ndims"
#define EC_GA_NC_DIMLENGHS         "dimLengths"
#define EC_GA_NC_NBINSCHEMES       "nbinSchemes"
#define EC_GA_NC_BINSCHEMES        "binSchemes"

#define EC_D_NC_RECORD             "records"
#define EC_D_NC_RECORDID           "ids"
#define EC_D_NC_NBINSCHEMES        "nbin_schemes"
#define EC_D_NC_NSTEPS             "nsteps"

#define EC_V_NC_RECORDID           "id"
#define EC_V_NC_ENTCONTRI          "contri_S"

class DimEntropy {
private:
	hbin_t nsteps_;
	hbin_t nbin_schemes_;
	hbin_t nestimators_;
	std::vector<u_int> dimid_;
	std::vector<double> contri_S_;
public:
	DimEntropy(const u_int id, const hbin_t nsteps, const hbin_t nbin_schemes,
			const hbin_t nestimators);

	DimEntropy(const std::vector<u_int> &id, const hbin_t nsteps,
			const hbin_t nbin_schemes, const hbin_t nestimators);

	friend std::ostream& operator<<(std::ostream &os, const DimEntropy &bp);

	void setId(const u_int &in);

	void setId(const ull_int start_index_in, const ull_int in_count,
			const std::vector<u_int> &in);

	// TODO: Make all Get??? members 'const' without adding any performance overhead
	void getId(hbin_t *size, u_int **p);

	void getDimContri(hbin_t *nestm, hbin_t *nsteps, hbin_t *nbin_schemes,
			double **data);

	void setDimContri(const hbin_t start_estm, const hbin_t count_estm,
			const hbin_t start_step_id, const hbin_t steps_count,
			const hbin_t start_bin_scheme_id, const hbin_t bin_scheme_count,
			const std::vector<double> &in);

	void setDimContri(const hbin_t start_estm, const hbin_t count_estm,
			const hbin_t start_step_id, const hbin_t steps_count,
			const hbin_t start_bin_scheme_id, const hbin_t bin_scheme_count,
			const ull_int start_in, const std::vector<double> &in);

	void setDimContri(hbin_t step_id, hbin_t bin_id, double value);

};

class Netcdf_EntContri: public NetcdfFile {
private:
	std::string filename_;
	std::string contents_;
	hbin_t ndiag_;
	hbin_t ndims_;
	std::vector<ull_int> dim_lens_;
	hbin_t order_;
	hbin_t nsteps_;
	TensorType storage_type_;
	hbin_t bin_schemes_count_;
	std::vector<hbin_t> bin_schemes_;
	hbin_t nestimators_;
	ull_int nrecords;
	std::vector<DimEntropy> vec_records;

	// variables for storing dim ids
	int recordsDID_;
	int estimDID_;
	int binSchemesDID_;
	int nstepsDID_;
	int recordIdDID_;

	// variables for storing variable ids
	int contri_SVID_;
	int recordIdVID_;

public:

	Netcdf_EntContri();
	Netcdf_EntContri(const std::string &filename, const std::string &c,
			const hbin_t ndims, const hbin_t order, const hbin_t nsteps,
			const std::vector<ull_int> &dim_lens, const TensorType &typ,
			const std::vector<hbin_t> &bin_schemes, const hbin_t nestimators);

	friend std::ostream& operator<<(std::ostream &os,
			const Netcdf_EntContri &ht);

	void clearRecords();

	int removeRecord(const DimEntropy &rec);

	int removeRecords(const std::vector<DimEntropy> &recs);

	int readRecords(const size_t start_index, const size_t count,
			bool is_append);

	int readIds(const u_int start_index, const u_int rec_counts,
			std::vector<u_int> &ids_v);

	int readEntrContrib(const u_int start_index, const u_int rec_counts,
			const hbin_t start_estm, const hbin_t n_estm,
			const hbin_t start_schemes, const hbin_t n_schemes,
			const hbin_t start_step, const hbin_t nsteps,
			std::vector<double> &dataOut);

	int readEntrContrib(const u_int start_index, const u_int rec_counts,
			const hbin_t start_estm, const hbin_t n_estm,
			const hbin_t start_schemes, const hbin_t n_schemes,
			const hbin_t start_step, const hbin_t nsteps, const hbin_t order,
			std::vector<u_int> &ids_v, std::vector<double> &dataOut);

	int writeRecord(DimEntropy &rec);

	int writeRecords(std::vector<DimEntropy> &vec_recs);

	int writeRecords();

	int addRecords(const std::vector<DimEntropy> &recs);

	int addRecords(std::vector<DimEntropy> &recs);

	int addRecord(const DimEntropy rec);

	void printRecords();

	int setupRead();

	int setupWrite();

	int NC_create(std::string const &title);
private:
	int NC_create(std::string &Name, std::string const &title);
};

#endif /* NETCDF_ENTCONTI_H */

