cmake_minimum_required(VERSION 2.6)

# If we need specific header file search paths add them here
#include_directories()

# Glob together files for each library.
## libciomr
file(GLOB CIOMR_SOURCES libciomr/*.cc)
file(GLOB CIOMR_HEADERS libciomr/*.h)
add_library(PSI_ciomr ${CIOMR_SOURCES})

## libchkpt
file(GLOB CHKPT_SOURCES libchkpt/*.cc)
file(GLOB CHKPT_HEADERS libchkpt/*.h libchkpt/*.hpp)
add_library(PSI_chkpt ${CHKPT_SOURCES})

## libqt
file(GLOB QT_SOURCES libqt/*.cc)
file(GLOB QT_HEADERS libqt/*.h)
add_library(PSI_qt ${QT_SOURCES})

## libpsio
file(GLOB PSIO_SOURCES libpsio/*.cc)
file(GLOB PSIO_HEADERS libpsio/*.h)
add_library(PSI_psio ${PSIO_SOURCES})

## libiwl
file(GLOB IWL_SOURCES libiwl/*.cc)
file(GLOB IWL_HEADERS libiwl/*.h libiwl/*.hpp)
add_library(PSI_iwl ${IWL_SOURCES})

## libdpd
file(GLOB DPD_SOURCES libdpd/*.cc)
file(GLOB DPD_HEADERS libdpd/*.h libdpd/*.hpp)
add_library(PSI_dpd ${DPD_SOURCES})

## libutil
file(GLOB UTIL_SOURCES libutil/*.cc)
file(GLOB UTIL_HEADERS libutil/*.h libutil/*.hpp)
add_library(PSI_util ${UTIL_SOURCES})

## liboptions
file(GLOB OPTIONS_SOURCES liboptions/*.cc)
file(GLOB OPTIONS_HEADERS liboptions/*.h liboptions/*.hpp)
add_library(PSI_options ${OPTIONS_SOURCES})

## lib3index
file(GLOB THREEINDEX_SOURCES lib3index/*.cc)
file(GLOB THREEINDEX_HEADERS lib3index/*.h lib3index/*.hpp)
add_library(PSI_3index ${THREEINDEX_SOURCES})

## libfunctional
file(GLOB FUNCTIONAL_SOURCES libfunctional/*.cc)
file(GLOB FUNCTIONAL_HEADERS libfunctional/*.h libfunctional/*.hpp)
add_library(PSI_functional ${FUNCTIONAL_SOURCES})

## libmints
file(GLOB MINTS_SOURCES libmints/*.cc)
file(GLOB MINTS_HEADERS libmints/*.h libmints/*.hpp)
add_library(PSI_mints ${MINTS_SOURCES})

## libtrans
file(GLOB TRANS_SOURCES libtrans/*.cc)
file(GLOB TRANS_HEADERS libtrans/*.h libtrans/*.hpp)
add_library(PSI_trans ${TRANS_SOURCES})

## libparallel
file(GLOB PARALLEL_SOURCES libparallel/*.cc)
file(GLOB PARALLEL_HEADERS libparallel/*.h libparallel/*.hpp)
add_library(PSI_parallel ${PARALLEL_SOURCES})

## libscf_solver
file(GLOB SCF_SOLVER_SOURCES libscf_solver/*.cc)
file(GLOB SCF_SOLVER_HEADERS libscf_solver/*.h libscf_solver/*.hpp)
add_library(PSI_scf_solver ${SCF_SOLVER_SOURCES})

## libsapt_solver
file(GLOB SAPT_SOLVER_SOURCES libsapt_solver/*.cc)
file(GLOB SAPT_SOLVER_HEADERS libsapt_solver/*.h libsapt_solver/*.hpp)
add_library(PSI_sapt_solver ${SAPT_SOLVER_SOURCES})

## libdfcc_solver
file(GLOB DFCC_SOLVER_SOURCES libdfcc_solver/*.cc)
file(GLOB DFCC_SOLVER_HEADERS libdfcc_solver/*.h libdfcc_solver/*.hpp)
add_library(PSI_dfcc_solver ${DFCC_SOLVER_SOURCES})

## libboost_python
file(GLOB BOOST_PYTHON_SOURCES libboost_python/*.cc)
file(GLOB BOOST_PYTHON_HEADERS libboost_python/*.h libboost_python/*.hpp)
add_library(boost_python ${BOOST_PYTHON_SOURCES})

## libboost_regex
file(GLOB BOOST_REGEX_SOURCES libboost_regex/*.cc)
file(GLOB BOOST_REGEX_HEADERS libboost_regex/*.h libboost_regex/*.hpp)
add_library(boost_regex ${BOOST_REGEX_SOURCES})

## libboost_thread
file(GLOB BOOST_THREAD_SOURCES libboost_thread/*.cc)
file(GLOB BOOST_THREAD_HEADERS libboost_thread/*.h libboost_thread/*.hpp)
add_library(boost_thread ${BOOST_THREAD_SOURCES})

## libboost_filesystem
file(GLOB BOOST_FILESYSTEM_SOURCES libboost_filesystem/*.cc)
file(GLOB BOOST_FILESYSTEM_HEADERS libboost_filesystem/*.h libboost_filesystem/*.hpp)
add_library(boost_filesystem ${BOOST_FILESYSTEM_SOURCES})

## libboost_system
file(GLOB BOOST_SYSTEM_SOURCES libboost_system/*.cc)
file(GLOB BOOST_SYSTEM_HEADERS libboost_system/*.h libboost_system/*.hpp)
add_library(boost_system ${BOOST_SYSTEM_SOURCES})

## libdiis
file(GLOB DIIS_SOURCES libdiis/*.cc)
file(GLOB DIIS_HEADERS libdiis/*.h libdiis/*.hpp)
add_library(PSI_diis ${DIIS_SOURCES})

## libplugin
file(GLOB PLUGIN_SOURCES libplugin/*.cc)
file(GLOB PLUGIN_HEADERS libplugin/*.h libplugin/*.hpp)
add_library(PSI_plugin ${PLUGIN_SOURCES})

## libmoinfo
file(GLOB MOINFO_SOURCES libmoinfo/*.cc)
file(GLOB MOINFO_HEADERS libmoinfo/*.h libmoinfo/*.hpp)
add_library(PSI_moinfo ${MOINFO_SOURCES})

## libint
file(GLOB LIBINT_SOURCES libint/libint/*.cc libint/libint/vrr_build.c)
file(GLOB LIBINT_HEADERS libint/libint/libint.h)
add_library(PSI_int ${LIBINT_SOURCES})

## libderiv
file(GLOB LIBDERIV_SOURCES libderiv/libderiv/*.cc)
file(GLOB LIBDERIV_HEADERS libderiv/libderiv/libderiv.h)
add_library(PSI_deriv ${LIBDERIV_SOURCES})

