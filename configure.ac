dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([centre], [0.1.0], [shaileshp51@gmail.com])


AC_ARG_ENABLE([full-info],
[  --enable-full-info    Turn on full progress information details in info file],
[case "${enableval}" in
  yes)  full_info=true ;; 
  no)  full_info=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-full-info]) ;;
esac],[full_info=false])

AS_IF([test x$full_info = xtrue], [AC_DEFINE([DETAILED_TIMING], [1], [Enable full progress information.])])

AC_ARG_ENABLE([debug-level],
[  --enable-debug-level    Turn on debugging],
[case "${enableval}" in
  0)  debug_level=0 ;; 
  1)  debug_level=1 ;;
  2)  debug_level=2 ;;
  3)  debug_level=3 ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-debug-level]) ;;
esac],[debug_level=0])

AC_ARG_ENABLE([mpi],
[  --enable-mpi    Turn on openmp-mpi hydrid compilation],
[case "${enableval}" in
  yes)  use_ompmpi=true ;; 
  no)  use_ompmpi=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-openmp-mpi]) ;;
esac],[use_ompmpi=false])

AM_CONDITIONAL([USE_OMPMPI], [test x$use_ompmpi = xtrue])

AS_IF([test x$use_ompmpi = xtrue], [AC_DEFINE([USE_OMPMPI], [1], [Enable omp-mpi hydrid build])])

AC_PROG_CXX
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_11
AC_OPENMP



if test "$use_ompmpi" = true; then
 CXX="mpic++"
 LDFLAGS+=" -L/usr/local/lib" 
 CXXFLAGS+=" -D__cplusplus=201103L -std=c++11 -O3 -fopenmp -fPIC -MMD "
 prog_suffix="_OMPMPI"
else 
 CXX="g++"
 LDFLAGS+=" -L/usr/local/lib"
 AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS -std=c++11 -O3 -fmessage-length=0 -fPIC")
 prog_suffix="_OMP"
fi

AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_AUX_DIR([build-aux])

AM_INIT_AUTOMAKE([-Wall -Werror foreign 1.13])

AC_CONFIG_HEADERS([config.h])

AS_IF([test x$use_ompmpi = xtrue], [AC_CHECK_HEADERS([mpi.h], [
                                        echo "Header mpi.h was found"; has_mpi=yes
                                    ], [
                                        echo "Header mpi.h was not found\n Please install a package providing mpi support e.g. openmpi or mpich"
                                       ]
                                    )
                                   ])

AC_CHECK_HEADERS([netcdf.h], [ AC_CHECK_LIB([netcdf], [nc_create], [
                                ], [
                                    echo "Error! You need to have netcdf around."
                                    exit -1
                                ])
                              has_netcdf=true
                             ], [ 
                              echo "Error! Header netcdf.h was not found\n Please install netcdf package" 
                              exit -1
                             ])


AC_CONFIG_FILES([Makefile src/Makefile])

AC_PROG_CXX

AC_OUTPUT

