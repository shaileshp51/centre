ABOUT
=====
The CENTRE is a parallel and scalable tool for estimating Configurational 
entropy from MD simulation data using Information Theoretic methods: 
     * Mutual Information Expansion (MIE), 
     * Maximum Information Spanning Tree (MIST), 
     * Approximate Mutual Information Expansion (or Neighbor Approximated MIE i.e. A-MIE), 
     * Neighbor Approximated Maximum Information Spanning Tree (A-MIST). 

This program accepts BOND/ANGLE/TORSION (BAT) format trajectory of the system 
for which the entropy is to be computed. An auxilary program cart2bat shipped 
with pycentre (https://github.com/shaileshp51/pycentre) python package can be used to 
generate a BAT format file from topology and trajectory files generated using popular 
Molecular Dynamics programs e.g. AMBER, NAMD, CHARMM, ACEMD, GROMACS.


DEPENDENCIES
============
NetCDF C Library version >= 4.6, with support for NETCDF4, NC_64BIT_OFFSET and 
NC_UBYTE.

NetCDF library can be installed using following commands on ubuntu based Linux
Distributions based on Ubuntu 18.04 or higher.

sudo apt install netcdf-bin libnetcdf11 libnetcdf-dev libhdf5-10 libhdf5-dev

NetCDF library source-code distribution can be ontained from following 
http://github.com/Unidata/netcdf-c 
and build following instructions provided on 

INSTALL
=======

Extract the distribution tar file using

tar -xvf centre-0.1.0.tar.bz
change directory to the extracted directory centre-0.1.0
cd centre-0.1.0

On LINUX
--------
Do an OpenMP thread based parallel build

./configure --prefix=`pwd`
make install

This step will create a file centre_OMP in bin folder, check with 

ls bin/centre_OMP

Do an MPI based proecess+thread based hybrid parallel build.
This build requires mpi headers and mpic++ processor capable to compile c++11 
code

./configure --enable-mpi --prefix=`pwd`
make install

This step will create a file centre_OMPMPI in bin folder, check with 
ls bin/centre_OMPMPI

On macOS
--------
On macOS using clang++ compiler
Do an OpenMP thread based parallel build
./configure --prefix=`pwd` CXX="clang++"

make install LIBS=-lomp

This step will create a file centre_OMP in bin folder, check with 

ls bin/centre_OMP

Do an MPI based proecess+thread based hybrid parallel build.
This build requires mpi headers and mpic++ processor capable to compile c++11 
code

export OMPI_CXX=clang++
./configure --enable-mpi --prefix=`pwd`
make install LIBS=-lomp

This step will create a file centre_OMPMPI in bin folder, check with 
ls bin/centre_OMPMPI
USAGE
=====
USAGE: 
centre_OMP [OPTION] -i \<input-file\>
Optional arguments for the CENTRE V0.1.0
  -h --help   print this help message and exit
  -s --sample print a sample input file with description of parameters and exit
  -O 		  overwite any already existing files in output directory
    
    
Running MPI+OMP hybrid build   
USAGE: 
mpirun -n $NUM_PROC centre_OMPMPI [OPTION] -i \<input-file\>
Optional arguments for the CENTRE V0.1.0
  -h --help    print this help message and exit
  -s --sample  print a sample input file with description of parameters and exit
  -O 		      overwite any already existing files in output directory


TERMS & CONDITIONS
==================
CENTRE  Copyright (C) 2020  Shailesh Kumar Panday
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License  version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
