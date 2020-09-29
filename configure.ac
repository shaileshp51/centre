dnl Process this file with autoconf to produce a configure script.
AC_PREREQ(2.62)
AC_INIT([centre], [0.1.1], [shaileshp51@gmail.com])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build-aux])

###############################################################################
# Let's get the host OS details and proceed accordingly.
###############################################################################
AC_CANONICAL_HOST

build_linux=no
build_mac=no

# Detect the target system
case "${host_os}" in
    linux*)
        build_linux=true
        ;;
    darwin*)
        build_mac=true
        ;;
    *)
    AC_MSG_ERROR(["Host OS $host_os is not supported"])
    ;;
esac

# Pass the conditionals to automake
AM_CONDITIONAL([LINUX], [test "$build_linux" = true])
AM_CONDITIONAL([OSX], [test "$build_mac" = true])

AS_IF([test x$build_mac = xtrue], 
      [AC_DEFINE([OS_IS_OSX], [1], [Current OS is macOS.])]
)
AS_IF([test x$build_linux = xtrue], 
      [AC_DEFINE([OS_IS_LINUX], [1], [Current OS is Linux.])]
)

##########################################################################


##########################################################################
# adding the netcdf library (e.g. with static name 'libnetcdf.a')
##########################################################################

# adding the lib to the files to link
LIBS="-lnetcdf $LIBS"

# introduce the optional configure parameter for a non-standard install prefix 
# of netcdf
AC_ARG_WITH([netcdf],
    [AS_HELP_STRING([--with-netcdf=prefix],
        [try this for a non-standard install prefix of the netcdf library])],
    [NETCDFPATHSET=1],
    [NETCDFPATHSET=0])

# if optional parameter used, extend path flags for compliler and linker
if test $NETCDFPATHSET = 1 ; then
    # extend the compiler and linker flags according to the path set
    AM_CXXFLAGS="$AM_CXXFLAGS -I$with_netcdf/include"
    AM_LDFLAGS="$AM_LDFLAGS -L$with_netcdf/lib"
else
    AM_CXXFLAGS="$AM_CXXFLAGS -I/usr/local/include"
    AM_LDFLAGS="$AM_LDFLAGS -L/usr/local/lib"
fi

##########################################################################


##########################################################################
# configuring for compilation with debug support
##########################################################################

AC_MSG_CHECKING([whether to build with debug information])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [enable debug data generation (def=no)])],
    [debugmode="$enableval"],
    [debugmode=no])
AC_MSG_RESULT([$debugmode])

if test x"$debugmode" = x"yes"; then
    AC_DEFINE([DEBUG],[],[Debug Mode])
    AM_CXXFLAGS="$AM_CXXFLAGS -g -Wall -Werror -Wno-uninitialized -O0"
else
    AC_DEFINE([NDEBUG],[],[No-debug Mode])
    AM_CXXFLAGS="$AM_CXXFLAGS -Ofast"
fi

##########################################################################


##########################################################################
# Whether pedentic progress info is requested
##########################################################################
AC_ARG_ENABLE([full-info],
    [AS_HELP_STRING([--enable-full-info],
        [Turn on full progress information details in info file])],
[case "${enableval}" in
  yes)  full_info=true ;; 
  no)  full_info=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-full-info]) ;;
esac],[full_info=false])

AS_IF([test x$full_info = xtrue], 
      [AC_DEFINE([DETAILED_TIMING], [1], 
          [Enable extra-detailed progress information.])])

##########################################################################


##########################################################################
# Lets check whether MPI build is requested
##########################################################################
AC_ARG_ENABLE([mpi],
[AS_HELP_STRING([--enable-mpi],
               [ Turn on openmp-mpi hydrid compilation])],
[case "${enableval}" in
  yes)  use_ompmpi=true ;; 
  no)  use_ompmpi=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-openmp-mpi]) ;;
esac],[use_ompmpi=false])

AM_CONDITIONAL([USE_OMPMPI], [test x$use_ompmpi = xtrue])

AS_IF([test x$use_ompmpi = xtrue], [AC_DEFINE([USE_OMPMPI], [1], 
      [Enable omp-mpi hydrid build])])

if test "$use_ompmpi" = true; then
 MPOMP=""
 if test x"$CXX" = x"clang++"; then
   MPOMP=""
   # echo "CXX=$CXX: MPOMP=$MPOMP"
 else
   MPOMP="-fopenmp -MMD"
   # echo "else CXX=$CXX: MPOMP=$MPOMP"
 fi
 CXX="mpic++"
 AM_CXXFLAGS="$AM_CXXFLAGS -fmessage-length=0 $MPOMP -fPIC "
 prog_suffix="_OMPMPI"
else 
 AM_CXXFLAGS="$AM_CXXFLAGS -fmessage-length=0 -fPIC "
 prog_suffix="_OMP"
fi

##########################################################################


CXXFLAGS="$AM_CXXFLAGS"
LDFLAGS="$AM_LDFLAGS"

AC_PROG_CXX
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_11
AC_OPENMP

AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS $AM_CXXFLAGS")

AM_INIT_AUTOMAKE([foreign 1.13])

AC_CONFIG_HEADERS([config.h])

AS_IF([test x$use_ompmpi = xtrue], 
      [AC_CHECK_HEADERS([mpi.h], 
          [ echo "Header mpi.h was found"; has_mpi=yes ], 
          [ echo "Header mpi.h was not found\n Please install a package providing mpi support e.g. openmpi or mpich" ]
          )
      ]
)

AC_CHECK_HEADERS([netcdf.h], 
		 [AC_CHECK_LIB([netcdf], 
			    [nc_create], 
			    [], 
			    [echo "ERROR! You need to have netcdf around."; exit -1 ])
     has_netcdf=true ], 
		 [echo "ERROR! netcdf.h was not found\n Please install netcdf package" 
      exit -1 ]
)

AC_SUBST([LIBS])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])

AC_CONFIG_FILES([Makefile src/Makefile])

AC_PROG_CXX

AC_OUTPUT

