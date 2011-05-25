cmake_minimum_required(VERSION 2.6)

## transqt2
file(GLOB TRANSQT2_SOURCES transqt2/*.cc)
add_library(transqt2 ${TRANSQT2_SOURCES})

## dcft
file(GLOB DCFT_SOURCES dcft/*.cc)
add_library(dcft ${DCFT_SOURCES})

## ccenergy
file(GLOB CCENERGY_SOURCES ccenergy/*.cc)
add_library(ccenergy ${CCENERGY_SOURCES})

## cctriples
file(GLOB CCTRIPLES_SOURCES cctriples/*.cc)
add_library(cctriples ${CCTRIPLES_SOURCES})

## ccsort
file(GLOB CCSORT_SOURCES ccsort/*.cc)
add_library(ccsort ${CCSORT_SOURCES})

## scf
file(GLOB SCF_SOURCES scf/*.cc)
add_library(scf ${SCF_SOURCES})

## mp2
file(GLOB MP2_SOURCES mp2/*.cc)
add_library(mp2 ${MP2_SOURCES})

## mcscf
file(GLOB MCSCF_SOURCES mcscf/*.cc)
add_library(mcscf ${MCSCF_SOURCES})

## df-mp2
file(GLOB DF_MP2_SOURCES df-mp2/*.cc)
add_library(df-mp2 ${DF_MP2_SOURCES})

## df-cc
file(GLOB DF_CC_SOURCES df-cc/*.cc)
add_library(df-cc ${DF_CC_SOURCES})

## sapt
file(GLOB SAPT_SOURCES sapt/*.cc)
add_library(sapt ${SAPT_SOURCES})

## cchbar
file(GLOB CCHBAR_SOURCES cchbar/*.cc)
add_library(cchbar ${CCHBAR_SOURCES})

## cclambda
file(GLOB CCLAMBDA_SOURCES cclambda/*.cc)
add_library(cclambda ${CCLAMBDA_SOURCES})

## ccdensity
file(GLOB CCDENSITY_SOURCES ccdensity/*.cc)
add_library(ccdensity ${CCDENSITY_SOURCES})

## transqt
file(GLOB TRANSQT_SOURCES transqt/*.cc)
add_library(transqt ${TRANSQT_SOURCES})

## mints
file(GLOB MINTS_SOURCES mints/*.cc)
add_library(mints ${MINTS_SOURCES})

## deriv
file(GLOB DERIV_SOURCES deriv/*.cc)
add_library(deriv ${DERIV_SOURCES})

## optking
file(GLOB OPTKING_SOURCES optking/*.cc)
add_library(optking ${OPTKING_SOURCES})

## psimrcc
file(GLOB PSIMRCC_SOURCES psimrcc/*.cc)
add_library(psimrcc ${PSIMRCC_SOURCES})

## TODO: Handle psi4
add_subdirectory(psi4)

