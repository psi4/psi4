cmake_minimum_required(VERSION 2.6)

file(GLOB PSI4_SOURCES *.cc)
add_executable(psi4 ${PSI4_SOURCES})
target_link_libraries(psi4 
    mints
    df-mp2
    scf
    ccenergy
    ccsort
    psimrcc
    transqt2
    cctriples
    mp2
    dcft
    mcscf
    sapt
    cchbar
    cclambda
    ccdensity
    cceom
    ccresponse
    transqt
    df-cc
    deriv
    optking
    PSI_sapt_solver
    PSI_scf_solver
    PSI_dfcc_solver
    PSI_functional
    PSI_3index
    PSI_mints
    PSI_trans
    PSI_dpd
    PSI_chkpt
    PSI_iwl
    PSI_psio
    PSI_qt
    PSI_ciomr
    PSI_options
    PSI_moinfo
    PSI_util
    PSI_deriv
    PSI_int
    PSI_util
    PSI_diis
    PSI_plugin
    PSI_parallel
    boost_python
    boost_regex
    boost_thread
    boost_filesystem
    boost_system
    ${PYTHON_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}
)

