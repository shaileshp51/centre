/* 
 * File:   config.h
 * Author: shailesh
 *
 * Created on 15 September, 2016, 3:44 PM
 */

#ifndef INC_CENTRE_CONFIG_H
#define INC_CENTRE_CONFIG_H

#ifndef LOG_LEVEL
#define LOG_LEVEL   ERROR_LEVEL
#endif

#include "../config.h"
#include <chrono>

//#ifndef USE_OMPMPI
//#define USE_OMPMPI
//#endif

#define DELTA_X_RANGE_EXTEND 5.0e-8

using hbin_t = unsigned char;
using data_t = float;
using ull_int = unsigned long long;
using u_int = unsigned int;


#ifdef OS_IS_OSX
using Clock = std::chrono::system_clock;
using ClockResolution = std::chrono::microseconds;
#endif

#ifdef OS_IS_LINUX
using Clock = std::chrono::system_clock;
using ClockResolution = std::chrono::nanoseconds;
#endif


enum TensorType {
	FULL = 0, LOWER = 1, UPPER = 3
};

/*! \mainpage CENTRE
 *
 * \section intro_sec Introduction
 * The CENTRE is a parallel and scalable tool for estimating Configurational entropy from MD simulation data using Information
 * Theoretic methods: Mutual Information Expansion (MIE), Maximum Information Spanning Tree (MIST), Approximate Mutual Information Expansion
 * (or Neighbor Approximated MIE i.e. A-MIE), and Neighbor Approximated Mutual Information Expansion (A-MIST). This program accepts
 * BOND/ANGLE/TORSION (BAT) format trajectory of the system for which the entropy is to be computed. An auxilary program \c cart2bat
 * shipped with \c pyceneter python package can be used to generate a BAT format file from topology and trajectory files generated
 * using popular Molecular Dynamics programs e.g. AMBER, NAMD, CHARMM, ACEMD, GROMACS.
 *
 * \section install_sec Installation
 *
 * \subsection requirements Step 0: Checking and Installing Required Libraries
 * This program requires NetCDF Library version 4.6.0 or higher. Whether NetCDF library exists can be checked
 * on Linux Distributions using command \c nc-config as below
 * \code nc-config --version \endcode
 * If output of above command is 4.6.0 or higher then we are ready to go for installation of CENTRE.
 * Otherwise, NetCDF library has to be installed before proceeding.
 * If Linux Distribution is \c Ubuntu \c >= \c 18.04, then these libraries and dependencies can be installed using
 * \code
 * sudo apt install netcdf-bin libnetcdf11 libnetcdf-dev libhdf5-10 libhdf5-dev
 * \endcode
 * For other distributions NetCDF source code can be downloaded and compiled from
 * http://github.com/Unidata/netcdf-c
 *
 * \subsection step1 Step 1: Install pycentre python3.6 package
 *
 *
 * etc...
 */

#endif

