ABOUT
=====
The CENTRE is a parallel and scalable tool for estimating Configurational 
entropy from MD simulation data using Information Theoretic methods: 
- Mutual Information Expansion (MIE), 
- Maximum Information Spanning Tree (MIST), 
- Approximate Mutual Information Expansion (or Neighbor Approximated MIE i.e. A-MIE), 
- Neighbor Approximated Maximum Information Spanning Tree (A-MIST). 

This program accepts BOND/ANGLE/TORSION (BAT) format trajectory of the system 
for which the entropy is to be computed. An auxilary program cart2bat.py shipped 
with pycentre (https://github.com/shaileshp51/pycentre) python package can be used to 
generate a BAT format file from topology and trajectory files generated using popular 
Molecular Dynamics programs e.g. AMBER, NAMD, CHARMM, ACEMD, GROMACS.


DEPENDENCIES
============
NetCDF C Library version >= 4.6, with support for NETCDF4, NC_64BIT_OFFSET and 
NC_UBYTE.

NetCDF library can be installed using following commands on ubuntu based Linux
Distributions based on Ubuntu 18.04 or higher.

```
        sudo apt install netcdf-bin libnetcdf11 libnetcdf-dev libhdf5-10 libhdf5-dev
```

NetCDF library source-code distribution can be ontained from following 
http://github.com/Unidata/netcdf-c 
and build following instructions provided on 

INSTALL
=======

Extract the distribution tar file using

```
tar -xvf centre-v0.1.2-b.tar.bz
```
change directory to the extracted directory centre-0.1.0
cd centre-0.1.0

On LINUX
--------
Check the location of netcdf include and library file. In some of the Linux
Distributions, these files are no longer installed in default locations i.e.
/usr/local/include and /usr/local/bin respectively. However, the configuration 
searches for them in default location and fails. In such a case, try finding
the locations of these files using command:

```
nc-config --all
```
It should print the detials of netcdf configuration on your system as below:
```
This netCDF 4.7.3 has been built with the following features: 

  --cc            -> /usr/bin/cc
  --cflags        -> -I/usr/include -I/usr/include/hdf5/serial
  --libs          -> -L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lnetcdf
  --static        -> -lhdf5_hl -lhdf5 -lpthread -lsz -lz -ldl -lm -lcurl
  
  ... many more lines
```
As above we see the include files are are in /usr/include and lib file in
/usr/lib/x86_64-linux-gnu. If so is the case please add

```
LDFLAGS="-L/usr/lib/x86_64-linux-gnu"
```
after --prefix option on on the ./configure command below.

Do an OpenMP thread based parallel build
```
./configure --prefix=`pwd`
make install
```
This step will create a file centre_OMP in bin folder, check with 
```
ls bin/centre_OMP
```
Do an MPI based proecess+thread based hybrid parallel build.
This build requires mpi headers and mpic++ processor capable to compile c++11 
code
```
./configure --enable-mpi --prefix=`pwd`
make install
```
This step will create a file centre_OMPMPI in bin folder, check with 
```
ls bin/centre_OMPMPI
```
On macOS
--------
On macOS using clang++ compiler
Do an OpenMP thread based parallel build
```
./configure --prefix=`pwd` CXX="clang++"
make install LIBS=-lomp
```
This step will create a file centre_OMP in bin folder, check with 
```
ls bin/centre_OMP
```
Do an MPI based proecess+thread based hybrid parallel build.
This build requires mpi headers and mpic++ processor capable to compile c++11 
code
```
export OMPI_CXX=clang++
./configure --enable-mpi --prefix=`pwd`
make install LIBS=-lomp
```
This step will create a file centre_OMPMPI in bin folder, check with 
ls bin/centre_OMPMPI

USAGE
=====
```
USAGE: 
centre_OMP [OPTION] -i <input-file>
Optional arguments for the CENTRE v0.1.2-b
    -h  --help       print this help message and exit
    -s  --sample     print a sample input file with description of parameters and exit
    -c  --cite       print citation details and exit
    -O               overwite any already existing files in output directory
```
    
    
Running MPI+OMP hybrid build   
```
USAGE: 
mpirun -n $NUM_PROC centre_OMPMPI [OPTION] -i <input-file>
Optional arguments for the CENTRE v0.1.2-b
    -h  --help       print this help message and exit
    -s  --sample     print a sample input file with description of parameters and exit
    -c  --cite       print citation details and exit
    -O               overwite any already existing files in output directory
```

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
