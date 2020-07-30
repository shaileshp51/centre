
if USE_OMPMPI
bin_PROGRAMS=centre_OMPMPI
LDADD = -lnetcdf
centre_OMPMPI_SOURCES= \
				Centre.cpp \
				main_common.h \
				macrologger/macrologger.h \
				CentreStdio.h \
				Convergence.h \
				Netcdf_EntContri.h \
				Netcdf_HistUtils.h \
				Real2Int.h \
				Utils.h \
				configcentre.h \
				CoordinateInfo.h \
				Netcdf_BAT.cpp \
				NetcdfFile.cpp \
				Netcdf_TrjInt.cpp \
				Version.h \
				Centre.h \
				Constants.h \
				Inputs.h \
				MpiCommTags.h \
				Estimators.cpp \
				Inputs.cpp \
				Netcdf_BAT.h \
				NetcdfFile.h \
				Netcdf_TrjInt.h \
				Timer.h \
				CentreStdio.cpp \
				Convergance.cpp \
				Estimators.h \
				CentreIException.h \
				CentreIException.cpp \
				Discretizer.h  \
				Exceptions.h \
				EntropyCalculator.h  \
				EntropyScorer.h  \
				Discretizer_mpi.cpp \
				EntropyCalculator_mpi.cpp \
				EntropyScorer.cpp \
				EntropyScorer.h \
				Inputs.h \
				Netcdf_EntContri.cpp \
				Netcdf_HistUtils.cpp \
				Real2Int.cpp \
				Utils.cpp \
				main_common.cpp 
centre_OMPMPI_SOURCES+= main_mpi.cpp
else
bin_PROGRAMS=centre_OMP
LDADD = -lnetcdf
centre_OMP_SOURCES= \
				Centre.cpp \
				macrologger/macrologger.h \
				main_common.h \
				CentreStdio.h \
				Convergence.h \
				Netcdf_EntContri.h \
				Netcdf_HistUtils.h \
				Real2Int.h \
				Utils.h \
				configcentre.h \
				CoordinateInfo.h \
				Netcdf_BAT.cpp \
				NetcdfFile.cpp \
				Netcdf_TrjInt.cpp \
				Version.h \
				Centre.h \
				Inputs.h \
				Constants.h \
				Estimators.cpp \
				Inputs.cpp \
				Netcdf_BAT.h \
				NetcdfFile.h \
				Netcdf_TrjInt.h \
				Timer.h \
				CentreStdio.cpp \
				Convergance.cpp \
				Estimators.h \
				Inputs.h \
				CentreIException.h \
				Discretizer.h  \
				EntropyCalculator.h  \
				EntropyScorer.h  \
				Discretizer.cpp \
				CentreIException.cpp \
				Exceptions.h \
				EntropyCalculator.cpp \
				EntropyScorer.cpp \
				EntropyScorer.h \
				Netcdf_EntContri.cpp \
				Netcdf_HistUtils.cpp \
				Real2Int.cpp \
				Utils.cpp \
				main_common.cpp
centre_OMP_SOURCES+= main_smp.cpp 
endif

