find_package(Perl QUIET)

# In order to get CTest to run tests in separate categories,
# labels are used. All tests must carry the label "psi".
#
# The following labels are defined:
#
#  smoketests  = Does Psi4 work? A few basic operational tests plus
#                one from each addon and addin
#  quicktests  = Best coverage possible in a 20 min envelope (for Travis-CI).
#                This subset is chosen to provide one or more tests for each
#                major module, but it does not necessarily cover all features
#                of every module.
#  longtests   = The really long (> ~15 min) tests
#
# Common test runs:
#
#  ctest -L quicktests      run all quicktests
#  ctest -LE longtests      run all but longtests
#
# It is possible to specify more than one label per test. CMake
# expects a ;-separated list, e.g. "psi;quicktests;scf"
# It is suggested to add at least a second label to specify the method tested.
#
#  We use ctest -I [start,end] where start is the test number we start with
#  and end is the test number we end with.  We will then use math to back out
#  the appropriate starts and ends...
#
#
#  In order to ensure that this works properly, please add your tests to the
#  appropriate variables given below

set(py3_fail_list pywrap-freq-e-sowreap
                  pywrap-freq-g-sowreap pywrap-opt-sowreap
                  pywrap-db2) 
#set(py36_fail_list extern1 extern2)
foreach(test_name adc1 adc2 casscf-fzc-sp casscf-sa-sp ao-casscf-sp casscf-sp castup1 
                  castup2 castup3 cbs-delta-energy cbs-xtpl-energy 
                  cbs-xtpl-freq cbs-xtpl-gradient cbs-xtpl-opt cbs-xtpl-func 
                  cbs-xtpl-wrapper cc1 cc10 cc11 cc12 cc13 cc13a cc13b cc13c cc13d cc14 cc15 cc16 
                  cc17 cc18 cc19 cc2 cc21 cc22 cc23 cc24 cc25 cc26 cc27 cc28 
                  cc29 cc3 cc30 cc31 cc32 cc33 cc34 cc35 cc36 cc37 cc38 cc39 
                  cc4 cc40 cc41 cc42 cc43 cc44 cc45 cc46 cc47 cc48 cc49 cc4a 
                  cc50 cc51 cc52 cc53 cc54 cc55 cc5a cc6 cc8 cc8a cc8b cc8c 
                  cc9 cc9a cdomp2-1 cdomp2-2 cepa0-grad1 cepa0-grad2 cepa1 
                  cepa2 cepa3 cepa-module ci-multi cisd-h2o+-0 cisd-h2o+-1 
                  cisd-h2o+-2 cisd-h2o-clpse cisd-opt-fd cisd-sp cisd-sp-2 
                  ci-property cubeprop decontract dcft-grad1 dcft-grad2 
                  dcft-grad3 dcft-grad4 dcft1 dcft2 dcft3 dcft4 dcft5 dcft6 
                  dcft7 dcft8 dcft9 ao-dfcasscf-sp dfcasscf-sa-sp dfcasscf-fzc-sp dfcasscf-sp 
                  dfccd1 dfccdl1 dfccd-grad1 dfccsd1 dfccsdl1 dfccsd-grad1 
                  dfccsdt1 dfccsdat1 dfmp2-1 dfmp2-2 dfmp2-3 dfmp2-4 dfmp2-grad1
                  dfmp2-grad2 dfmp2-grad3 dfmp2-grad4 dfmp2-grad5 dfomp2-1 dfomp2-2 dfomp2-3
                  dfomp2-4 dfomp2-grad1 dfomp2-grad2 dfomp3-1 dfomp3-2 
                  dfomp3-grad1 dfomp3-grad2 dfomp2p5-1 dfomp2p5-2 dfomp2p5-grad1
                  dfomp2p5-grad2 dfrasscf-sp dfscf-bz2 dft-b2plyp dft-dldf 
                  dft-freq dft-grad dft-pbe0-2 dft-psivar dft-b3lyp dft1 
                  dft1-alt dft2 dft3 docs-bases docs-dft extern1 extern2
                  fsapt1 fsapt2 isapt1 isapt2
                  fci-dipole fci-h2o fci-h2o-2 fci-h2o-fzcv fci-tdm fci-tdm-2 
                  fd-freq-energy fd-freq-energy-large fd-freq-gradient 
                  fd-freq-gradient-large fd-gradient freq-isotope fnocc1 fnocc2 
                  fnocc3 fnocc4 frac ghosts gibbs matrix1 mcscf1 mcscf2 mcscf3 
                  mints1 mints2 mints3 mints4 mints5 mints6 mints8 
                  mints9 molden1 molden2 mom mp2-1 mp2-def2 mp2-grad1 mp2-grad2 
                  mp2-module mp2p5-grad1 mp2p5-grad2 mp3-grad1 mp3-grad2 
                  mp2-property mpn-bh nbody-he-cluster numpy-array-interface 
                  ocepa-freq1 ocepa-grad1 ocepa-grad2 ocepa1 ocepa2 ocepa3 
                  omp2-1 omp2-2 omp2-3 omp2-4 omp2-5 omp2-grad1 omp2-grad2 
                  omp2p5-1 omp2p5-2 omp2p5-grad1 omp2p5-grad2 omp3-1 omp3-2 
                  omp3-3 omp3-4 omp3-5 omp3-grad1 omp3-grad2 opt-lindep-change 
                  opt1 opt1-fd opt2 opt2-fd opt3 opt4 opt5 opt6 opt7 opt8 opt9 
                  opt11 opt12 opt13 opt14 opt-irc-1 opt-irc-2 opt-freeze-coords 
                  props1 props2 props3 psimrcc-ccsd_t-1 psimrcc-ccsd_t-2 
                  psimrcc-ccsd_t-3 psimrcc-ccsd_t-4 psimrcc-fd-freq1 
                  psimrcc-fd-freq2 psimrcc-pt2 psimrcc-sp1 psithon1 psithon2 
                  pubchem1 pubchem2 pywrap-alias pywrap-all pywrap-basis 
                  pywrap-cbs1 pywrap-checkrun-convcrit pywrap-checkrun-rhf 
                  pywrap-checkrun-rohf pywrap-checkrun-uhf pywrap-db1 pywrap-db2
                  pywrap-db3 pywrap-freq-e-sowreap pywrap-freq-g-sowreap 
                  pywrap-molecule pywrap-opt-sowreap rasci-c2-active rasci-h2o 
                  rasci-ne rasscf-sp sad1 sapt1 sapt2 sapt3 sapt4 sapt5 sapt6 
                  sapt7 sapt8 scf-bz2 scf-guess-read1 scf-guess-read2 scf-bs scf1 scf-occ
                  scf2 scf3 scf4 scf5 scf6 scf7 scf-property soscf1 soscf2 stability1
                  stability2 tu1-h2o-energy tu2-ch2-energy tu3-h2o-opt 
                  tu4-h2o-freq tu5-sapt tu6-cp-ne2 x2c1 x2c2 x2c3 zaptn-nh2 
                  options1 cubeprop-esp
)
    # Some tests fail in python 3
    if("${PYTHON_VERSION_MAJOR}" VERSION_EQUAL "3")
        if (${test_name} IN_LIST py3_fail_list)
            message("Python3 incompatible test: ${test_name}... skipping")
        elseif(("${PYTHON_VERSION_MINOR}" VERSION_EQUAL "6") AND
               (${test_name} IN_LIST py36_fail_list))
            message("Python3.6 recalcitrant test: ${test_name}... skipping")
        else()
            add_subdirectory(${test_name})
        endif()
    else()
       add_subdirectory(${test_name})
    endif()
endforeach()

# Add internal projects
add_subdirectory(psi4numpy)
add_subdirectory(cookbook)
add_subdirectory(python)
add_subdirectory(json)
if(ENABLE_pasture)
  add_subdirectory(pasture-ccsorttransqt2)
  message(STATUS "${Cyan}Found Pasture${ColourReset}")
  message(STATUS "Adding test cases: Psi4 + Pasture")
endif()

message(STATUS "Adding test cases: Psi4")
if(NOT PERL_FOUND)
   message(WARNING "Perl was not found. Some Coupled Cluster test cases will not be thoroughly checked.")
endif()

###########################  External Project Tests  ###########################

# For code that (1) is built or detected in the Psi4 superbuild but isn't
#   maintained in the Psi4 codebase and (2) has Psi4-dependent tests _not_ in
#   the external codebase, tests should
#   * get their own directory in tests (e.g., tests/dftd3)
#   * have a name that doesn't repeat the exe (e.g., tests/dftd3/energy, not test/dftd3/dftd3-energy)
#   * be listed in the test case CMakeLists.txt as addon-testname (e.g., `add_regression_test(dftd3-energy ...`)
#   * alphabetize the addons

# <<<  CFOUR  >>>

find_program(CFOUR_EXECUTABLE
             NAMES xcfour
             DOC "Path to Stanton and Gauss's xcfour executable")

if(CFOUR_EXECUTABLE)
    add_subdirectory(cfour)
    message(STATUS "${Cyan}Found CFOUR${ColourReset}: ${CFOUR_EXECUTABLE}")
    message(STATUS "Adding test cases: Psi4 + CFOUR")
endif()

# <<<  CHEMPS2  >>>

if(ENABLE_CheMPS2)
    add_subdirectory(chemps2)
    message(STATUS "Adding test cases: Psi4 + CheMPS2")
endif()

# <<<  DFTD3  >>>

find_program(DFTD3_EXECUTABLE
             NAMES dftd3
             DOC "Path to Grimme's dftd3 executable")

if(DFTD3_EXECUTABLE)
    add_subdirectory(dftd3)
    message(STATUS "${Cyan}Found DFTD3${ColourReset}: ${DFTD3_EXECUTABLE}")
    message(STATUS "Adding test cases: Psi4 + DFTD3")
endif()

# <<<  DKH  >>>

if(ENABLE_dkh)
    add_subdirectory(dkh)
    #message(STATUS "Adding test cases: Psi4 + dkh")  # uncomment when there's a test to add
endif()

# <<<  LIBEFP  >>>

if(ENABLE_libefp)
    add_subdirectory(libefp)
    message(STATUS "Adding test cases: Psi4 + libefp")
endif()

# <<<  ERD  >>>

if(ENABLE_erd)
    add_subdirectory(erd)
    message(STATUS "Adding test cases: Psi4 + erd")
endif()

# <<<  gCP  >>>

find_program(gCP_EXECUTABLE
             NAMES gcp
             DOC "Path to Grimme's gCP executable")

if(gCP_EXECUTABLE)
    add_subdirectory(gcp)
    message(STATUS "${Cyan}Found gCP${ColourReset}: ${gCP_EXECUTABLE}")
    message(STATUS "Adding test cases: Psi4 + gCP")
endif()

# <<<  GDMA  >>>

if(ENABLE_gdma)
    add_subdirectory(gdma)
    message(STATUS "Adding test cases: Psi4 + gdma")
endif()

# <<<  MRCC  >>>

find_program(MRCC_EXECUTABLE
             NAMES dmrcc
             DOC "Path to Kallay's mrcc executable")

if(MRCC_EXECUTABLE)
    add_subdirectory(mrcc)
    message(STATUS "${Cyan}Found MRCC${ColourReset}: ${MRCC_EXECUTABLE}")
    message(STATUS "Adding test cases: Psi4 + MRCC")
endif()

# <<<  PCMSOLVER  >>>

if(ENABLE_PCMSolver)
    add_subdirectory(pcmsolver)
    message(STATUS "Adding test cases: Psi4 + PCMSolver")
endif()

# <<<  SIMINT  >>>

if(ENABLE_simint)
    add_subdirectory(simint)
    message(STATUS "Adding test cases: Psi4 + simint")
endif()

# plugin regime
if(ENABLE_PLUGIN_TESTING)
    add_subdirectory(plugins)
    message(STATUS "Adding test cases: Psi4 plugins")
endif()

# pytest regime
add_subdirectory(pytest)
