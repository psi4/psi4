#These tests aren't set-up the same way as the others
set(TestsNotSetUpRight cfour optking testgrendel)

#In order to get CTest to run tests in three seperate categories:
#
#  LongTests=  The (as of right now) 4 really long tests
#  Tests=      Every test in the test folder less LongTests
#  QuickTests= Proper subset of tests deemed to run quickly
#
#  We use ctest -I [start,end] where start is the test number we start with
#  and end is the test number we end with.  We will then use math to back out
#  the appropriate starts and ends...
#
#
#  In order to ensure that this works properly, please add your tests to the
#  appropriate variables given below




set(ADCLong adc1 adc2)

set(CastUpShort castup1)
set(CastUpLong castup2 castup3)
set(CCShort cc1 cc2 cc3 cc4 cc8 cc8a cc12 cc15 cc16 cc17 cc18 cc19 cc23 cc24 cc28 cc32 cc49 )
set(CCLong cc10 cc11 cc13 cc13a cc14 cc21 cc22 cc25 cc26 cc27 cc29 cc30 cc31 cc33 cc34 cc35 cc36 cc37 cc38 cc39 cc40 cc41 cc42 cc43 cc44 cc45 cc46 cc47 cc48 cc4a cc5 cc50 cc51 cc52 cc53 cc54 cc55 cc5a cc6 cc8b cc8c cc9 cc9a)

set(CDShort cdomp2-1 cdomp2-2)
#No long CD jobs
#No CEPAShort
set(CEPALong cepa0-grad1 cepa0-grad2 cepa1 cepa2 cepa3)
set(CISDShort cisd-h2o+-0 cisd-opt-fd)
set(CISDLong cisd-h2o+-1 cisd-h2o+-2 cisd-h2o-clpse cisd-sp cisd-sp-2)

set(DCFTShort dcft-grad1 dcft1 dcft2 dcft3 dcft4)
set(DCFTLong dcft5 dcft6 dcft7)
set(DFShort dfmp2-1 dfomp2-1 dfomp2-2)
set(DFLong dfmp2-2 dfmp2-3 dfmp2-4 dfomp2-3 dfomp2-4 dfomp2-grad1 dfomp2-grad2 dfscf-bz2)
set(DFTShort dft1)
set(DFTLong dft1-alt dft2 dft3 dft-b2plyp dft-dldf dft-freq dft-grad dft-pbe0-2 dft-psivar)
#No Docs Short
set(DocsLong docs-bases docs-dft docs-psimod)

set(FCIShort fci-dipole)
set(FCILong fci-h2o fci-h2o-2 fci-h2o-fzcv fci-tdm fci-tdm-2)
set(FDLong fd-freq-energy fd-freq-energy-large fd-freq-gradient fd-gradient)
set(FNOCCShort fnocc2 fnocc4)
set(FNOCCLong fnocc1 fnocc3) 

set(MCSCFShort mcscf1)
set(MCSCFLong mcscf2 mcscf3)
set(MintsShort mints1 mints2)
set(MintsLong mints3 mints4 mints5 mints6 mints8)
set(MP2Short mp2_5-grad1 mpn-bh)
set(MP2Long mp2-1 mp2_5-grad2 mp2-grad1 mp2-grad2 mp3-grad1 mp3-grad2)
#No MRCCShort
set(MRCCLong mrcc1 mrcc2 mrcc3 mrcc4)

set(OCEPAShort ocepa1 ocepa-grad2)
set(OCEPALong ocepa2 ocepa3 ocepa-freq1 ocepa-grad1)
set(OMPShort omp2-2 omp3-1 omp2_5-2 omp2-grad1)
set(OMPLong omp2-1 omp2-3 omp2-4 omp2-5 omp2_5-1 omp2_5-grad1 omp2_5-grad2 omp2-grad2 omp3-2 omp3-3 omp3-4 omp3-5 omp3-grad1 omp3-grad2) 
set(OPTShort opt1 opt1-fd)
set(OPTLong opt2 opt2-fd opt3 opt4 opt5 opt6 opt7 opt-lindep-change)

#No props short
set(PropsLong props1 props2 props3)
set(PSIMRCCShort psimrcc-sp1)
set(PSIMRCCLong psimrcc-ccsd_t-1 psimrcc-ccsd_t-2 psimrcc-ccsd_t-3 psimrcc-ccsd_t-4 psimrcc-fd-freq1 psimrcc-fd-freq2 psimrcc-pt2)
set(PYWrapShort pywrap-db1 pywrap-db2 pywrap-cbs1 pywrap-all pywrap-alias pywrap-molecule pywrap-freq-e-sowreap)
set(PYWrapLong pywrap-basis pywrap-checkrun-convcrit pywrap-checkrun-rhf pywrap-checkrun-rohf pywrap-checkrun-uhf pywrap-db3 pywrap-opt-sowreap)

set(RASCIShort rasci-h2o)
set(RASCILong rasci-c2-active rasci-ne)

set(SAPTShort sapt1)
set(SAPTLong  sapt3 sapt5)
set(SCFShort scf1 scf2 scf3 scf5 scf-guess-read)
set(SCFLong scf11-freq-from-energies scf4 scf6 scf-bz2)

set(TutorialShort tu1-h2o-energy tu2-ch2-energy tu3-h2o-opt tu4-h2o-freq tu5-sapt )
set(TutorialLong tu6-cp-ne2)
set(MISCShort ci-multi sad1 matrix1 gibbs)
set(MISCLong frac ghosts large_atoms min_input mom psithon1 psithon2 pubchem1 zaptn-nh2)

set(AutoTestSubjects cc8 cc9 cc11 cc18 cc19 cc28 cc29 cc30 cc31 cc34 cc36 cc37 cc38 cc39 cc40 cc41 cc42 cc43 cc46 cc47 cc48 cc49 cc50 cc51 cc52 cc55)

###### END test definitions ################

#These are the truly long tests
set(LongTests mp2-def2 sapt2 sapt4 fd-freq-gradient-large)

#These are the tests that will become the quicktest target
set(QuickTests ${CastUpShort} ${CCShort} ${CDShort} ${CISDShort} ${DCFTShort}
    ${DFShort} ${DFTShort} ${FCIShort} ${FNOCCShort} ${MCSCFShort}
${MintsShort} ${MP2Short} ${OCEPAShort} ${OMPShort} ${OPTShort}
${PSIMRCCShort} ${PYWrapShort} ${RASCIShort} ${SAPTShort} ${SCFShort}
${TutorialShort} ${MISCShort}
)

#And these are the rest of the tests
set(LongerTests ${ADCLong} ${CastUpLong} ${CCLong} ${CEPALong} ${CISDLong}
${DCFTLong} ${DFLong} ${DFTLong} ${DocsLong} ${FCILong} ${FDLong} ${FNOCCLong}
${MCSCFLong} ${MintsLong} ${MP2Long} ${MRCCLong} ${OCEPALong} ${OMPLong}
${OPTLong} ${PropsLong} ${PSIMRCCLong} ${PYWrapLong} ${RASCILong} ${SAPTLong}
${SCFLong} ${TutorialLong} ${MISCLong})

###### START external-executable-dependent tests ################

# Set up groups of test cases analogous to dftd3 below
#   In particular, new ext-exe-dep tests should 
#   * get their own directory in tests (e.g., tests/dftd3)
#   * have a name that doesn't repeat the exe (e.g., test/dftd3/energy)
#   * be listed in the exeShort and exeLong variable name as exe-name (e.g., dftd3-energy)
#   * be sure the exe is listed in the python code below (search 'partition')

# <<<  DFTD3  >>>

set(DFTD3Short dftd3-energy)
set(DFTD3Long dftd3-grad dftd3-psithon2)

find_program (DFTD3_EXECUTABLE
             NAMES dftd3
             DOC "Path to Grimme's dftd3 executable")

if(DFTD3_EXECUTABLE)
    list(APPEND QuickTests ${DFTD3Short})
    list(APPEND LongerTests ${DFTD3Long})
    message("-- dftd3 interfaced executable found. Test cases added.")
endif()

###### END external-executable-dependent tests ################

set(TestSet ${LongTests} ${QuickTests} ${LongerTests})

list(LENGTH TestSet NumberOfTests)
list(LENGTH LongTests NLong)
list(LENGTH QuickTests NShort)
MATH(EXPR Start "${NLong} +1")
MATH(EXPR Stop "${NLong} + ${NShort}")

add_custom_target(longtests ${CMAKE_CTEST_COMMAND} -I 1,${NLong})
add_custom_target(quicktests ${CMAKE_CTEST_COMMAND} -I ${Start},${Stop})
add_custom_target(tests ${CMAKE_CTEST_COMMAND} -I ${Start},${NumberOfTests})


#The commands below this point are used to set up cmake's test suite with
#the info given up to this point

#This is the psi command to actually run the tests
set(PSIEXE ${CMAKE_BINARY_DIR}/bin/psi4)
#This is where the test directories live
set(TESTDIR ${CMAKE_CURRENT_SOURCE_DIR})
#This is the perl script that we call, to call the psi, to run the tests
set(TESTEXE ${TESTDIR}/runtest.py)

#Finally run these beasts
set(LOGFILE ${CMAKE_BINARY_DIR}/testresults.log)

find_package(Perl QUIET)
if(NOT PERL_FOUND)
   message(WARNING "Perl was not found. Some Coupled Cluster test cases will
not be thoroughly checked.")
endif()

# Declare each test as a possible test target even though no one will likely
#   want to run individual tests. In the end, makes our job easier.
set(index 1)
foreach(Tests2Run ${TestSet})

    # Check if a test case from subdirectory (ext-exec-dep)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c 
        "pieces = '${Tests2Run}'.partition('-'); print '%s;%s' % (pieces[0], pieces[2]) if pieces[0] in ['dftd3', 'cfour'] else ''"
        OUTPUT_VARIABLE TEST_PIECES OUTPUT_STRIP_TRAILING_WHITESPACE)
    list(LENGTH TEST_PIECES LEN_PIECES)

    # Directory the test will be run from and in
    if(${LEN_PIECES} EQUAL 2)
        list(GET TEST_PIECES 0 SubDir)
        list(GET TEST_PIECES 1 TestName)
        set(TEST_SRC_DIR ${TESTDIR}/${SubDir}/${TestName})
        set(TEST_RUN_DIR ${CMAKE_BINARY_DIR}/tests/${SubDir}/${TestName})
    else() 
        set(TEST_SRC_DIR ${TESTDIR}/${Tests2Run})
        set(TEST_RUN_DIR ${CMAKE_BINARY_DIR}/tests/${Tests2Run})
    endif()

    # Generic setup
    file(MAKE_DIRECTORY ${TEST_RUN_DIR})
    set(INPUTFILE ${TEST_SRC_DIR}/input.dat)
    set(OUTFILE ${TEST_RUN_DIR}/output.dat)

    # Turn on psitest.pl if eligible 
    set(AUTOTEST false)
    list(FIND AutoTestSubjects ${Tests2Run} InThrall)
    if(NOT(${InThrall} EQUAL -1) AND ${PERL_FOUND})
        set(AUTOTEST true)
    endif()

    # Add (serial) target for each test case
    add_test(
        NAME ${Tests2Run}
        WORKING_DIRECTORY ${TEST_RUN_DIR}
        COMMAND ${PYTHON_EXECUTABLE} ${TESTEXE} ${INPUTFILE} ${LOGFILE} ${AUTOTEST} ${PROJECT_SOURCE_DIR} ${OUTFILE} ${PSIEXE})
    add_custom_target(${Tests2Run} ${CMAKE_CTEST_COMMAND} -I ${index},${index})
    set_tests_properties(${Tests2Run} PROPERTIES DEPENDS ${PSIEXE})

    #In theory would like to make sure all tests run in parallel as well
    #but....most do not right now...
    if( MPI_FOUND AND FALSE)
        add_test(
            NAME MPI-${Tests2Run}
            WORKING_DIRECTORY ${RUNDIR}
            COMMAND mpirun -n 2 ${PYTHON_EXECUTABLE} ${TESTEXE} ${INPUTFILE} ${LOGFILE} ${AUTOTEST} ${PROJECT_SOURCE_DIR} ${OUTFILE} ${PSIEXE})
    endif()
    math(EXPR index "${index}+1")
endforeach(Tests2Run)
