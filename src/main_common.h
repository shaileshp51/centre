/*
 * main.h
 *
 *  Created on: 01-Dec-2016
 *      Author: shailesh
 */

#ifndef SRC_MAIN_H_
#define SRC_MAIN_H_

#include <cstdlib>
#include <chrono>
#include <ctime>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <fstream>
#include <cstddef>
#include <utility>
#include <memory>
#include <limits>
#include <string>
#include <cmath>
#include <map>
#include <random>
#include <algorithm>
#include <numeric>

#include <omp.h>

#include "configcentre.h"
#include "Timer.h"
#include "Centre.h"
#include "NetcdfFile.h"
#include "Netcdf_HistUtils.h"
#include "Netcdf_BAT.h"
#include "Netcdf_TrjInt.h"
#include "Netcdf_EntContri.h"
#include "Inputs.h"
#include "Real2Int.h"
#include "Estimators.h"
#include "Convergence.h"

#define BINTRAJ 1

enum ExecState {
	DISABLED = 0, WAITING = 1, RUNNING = 2, COMPLETED = 4
};

struct ProgTaskSet {
	typedef Clock::time_point CTimePoint;
	std::string name;
	bool active;
	ExecState cstt;
	CTimePoint start;
	CTimePoint current;
	CTimePoint strt_read;
	CTimePoint curr_read;
	CTimePoint strt_comp;
	CTimePoint curr_comp;
	CTimePoint strt_comm;
	CTimePoint curr_comm;
	CTimePoint strt_write;
	CTimePoint curr_write;
	ull_int done_tasks;
	ull_int tot_tasks;

	ProgTaskSet(std::string const &name_) {
		name = std::string(name_.c_str());
		active = false;
		cstt = ExecState::DISABLED;
		auto NOW = Clock::now();
		start = NOW;
		current = NOW;
		strt_read = NOW;
		curr_read = NOW;
		strt_comp = NOW;
		curr_comp = NOW;
		strt_comm = NOW;
		curr_comm = NOW;
		strt_write = NOW;
		curr_write = NOW;
		done_tasks = 0;
		tot_tasks = 0;
	}

	ProgTaskSet() {
		name = std::string("UN-NAMED");
		active = false;
		cstt = ExecState::DISABLED;
		auto NOW = Clock::now();
		start = NOW;
		current = NOW;
		strt_read = NOW;
		curr_read = NOW;
		strt_comp = NOW;
		curr_comp = NOW;
		strt_comm = NOW;
		curr_comm = NOW;
		strt_write = NOW;
		curr_write = NOW;
		done_tasks = 0;
		tot_tasks = 0;
	}
	friend std::ostream& operator<<(std::ostream &os, const ProgTaskSet &tsk);

	std::string toString();
};

std::ostream& operator<<(std::ostream &os, const ProgTaskSet &elm);

char* getCmdOption(char **begin, char **end, const std::string &option);

bool cmdOptionExists(char **begin, char **end, const std::string &option);

std::string to_durration_string(int64_t dur, bool show_ms = true);

struct ProgState {
	ProgTaskSet dscrt;

	ProgTaskSet dscrt_b;
	ProgTaskSet dscrt_a;
	ProgTaskSet dscrt_d;

	ProgTaskSet entc;
	ProgTaskSet entc_b1d;
	ProgTaskSet entc_a1d;
	ProgTaskSet entc_d1d;

	ProgTaskSet entc_b2d;
	ProgTaskSet entc_a2d;
	ProgTaskSet entc_d2d;

	ProgTaskSet entc_ba2d;
	ProgTaskSet entc_bd2d;
	ProgTaskSet entc_ad2d;

	ProgState() {
		dscrt = ProgTaskSet("Discrete");
		dscrt_b = ProgTaskSet("BND");
		dscrt_a = ProgTaskSet("ANG");
		dscrt_d = ProgTaskSet("DIH");
		entc = ProgTaskSet("Entropy");

		entc_b1d = ProgTaskSet("BND-1D");
		entc_a1d = ProgTaskSet("ANG-1D");
		entc_d1d = ProgTaskSet("DIH-1D");

		entc_b2d = ProgTaskSet("B/B-2D");
		entc_a2d = ProgTaskSet("A/A-2D");
		entc_d2d = ProgTaskSet("D/D-2D");

		entc_ba2d = ProgTaskSet("B/A-2D");
		entc_bd2d = ProgTaskSet("B/D-2D");
		entc_ad2d = ProgTaskSet("A/D-2D");
	}

	friend std::ostream& operator<<(std::ostream &os, const ProgState &stt);

	std::string toString();
};

std::ostream& operator<<(std::ostream &os, const ProgState &stt);

bool cmpNodeMIST(nodeMIST &lh, nodeMIST &rh);

bool cmpEdgeMIST(const Edge &lh, const Edge &rh);

void printHeader();

void printUsage();

void printSample();

void printFooter();

#endif /* SRC_MAIN_H_ */
