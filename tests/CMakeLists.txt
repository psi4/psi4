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
#  ctest -L quick      run all quicktests
#  ctest -LE long      run all but longtests
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

foreach(test_name aediis-1
                  casscf-fzc-sp casscf-semi casscf-sa-sp ao-casscf-sp casscf-sp castup1
                  castup2 castup3 cbs-delta-energy cbs-parser cbs-xtpl-alpha cbs-xtpl-energy
                  cbs-xtpl-freq cbs-xtpl-gradient cbs-xtpl-opt cbs-xtpl-func cbs-xtpl-nbody
                  cbs-xtpl-wrapper cbs-xtpl-dict cc1 cc10 cc11 cc12 cc13 cc13a cc13b cc13c
                  cc13d cc14 cc15 cc16 cc17 cc18 cc19 cc2 cc21 cc22 cc23 cc24 cc25 cc26 cc27 cc28
                  cc29 cc3 cc30 cc31 cc32 cc33 cc34 cc35 cc36 cc37 cc38 cc39
                  cc4 cc40 cc41 cc42 cc43 cc44 cc45 cc46 cc47 cc48 cc49 cc4a
                  cc50 cc51 cc52 cc53 cc54 cc55 cc56 cc5 cc6 cc7 cc8 cc8a cc8b cc8c
                  cc9 cc9a cdomp2-1 cdomp2-2 cdoremp-energy1 cdoremp-energy2 cdremp-1 cdremp-2 cepa1
                  cepa2 cepa3 cepa-module ci-multi cisd-h2o+-0 cisd-h2o+-1
                  cisd-h2o+-2 cisd-h2o-clpse cisd-opt-fd cisd-sp cisd-sp-2
                  ci-property cubeprop cubeprop-frontier decontract dct-grad1 dct-grad2
                  dct-grad3 dct-grad4 dct1 dct2 dct3 dct4 dct5 dct6 dct7 dct8 dct9
                  dct10 dct11 dct12 ao-dfcasscf-sp density-screen-1 density-screen-2 dfcasscf-sa-sp
                  dfcasscf-fzc-sp dfcasscf-sp dfccd1 dfccdl1 dfccd-grad1 dfccsd1 dfccsdl1 dfccsd-grad1
                  dfccsd-t-grad1
                  dfccsdt1 dfccsdat1 dfmp2-1 dfmp2-2 dfmp2-3 dfmp2-4 dfmp2-5 dfmp2-fc dfmp2-freq1 dfmp2-freq2
                  dfccsd-grad2 dfccsd-t-grad2 dfccsdat2 dfccsdt2
                  dfmp2-grad1 dfmp2-grad2 dfmp2-grad3 dfmp2-grad4 dfmp2-grad5 dfomp2-1 dfomp2-2 dfomp2-3
                  dfomp2-4 dfomp2-grad1 dfomp2-grad2 dfomp2-grad3 dfomp3-1 dfomp3-2
                  dfomp3-grad1 dfomp3-grad2 dfomp2p5-1 dfomp2p5-2 dfomp2p5-grad1
                  dforemp-grad1 dforemp-grad2 dfremp-1 dfremp-2
                  dft-grad-lr1 dft-grad-lr2 dft-grad-lr3 dft-grad-disk
                  dfomp2p5-grad2 dfrasscf-sp dfscf-bz2 dft-b2plyp dft-grac dft-ghost dft-grad-meta
                  dft-freq dft-freq-analytic1 dft-freq-analytic2 dft-grad1 dft-grad2 dft-psivar dft-b3lyp dft1 dft-vv10
                  dft1-alt dft2 dft3 dft-omega dft-dens-cut dlpnomp2-1 dlpnomp2-2 dlpnomp2-3
                  docs-bases docs-dft embpot1 explicit-am-basis extern1 extern2 extern3 extern4 extern5
                  fsapt1 fsapt2 fsapt-terms fsapt-allterms fsapt-ext fsapt-ext-abc fsapt-ext-abc2
                  fsapt-ext-abc-au isapt1 isapt2 isapt-siao1 fisapt-siao1 isapt-charged
                  fci-dipole fci-h2o fci-h2o-2 fci-h2o-fzcv fci-tdm fci-tdm-2
                  fci-coverage
                  fcidump
                  fd-freq-energy fd-freq-energy-large fd-freq-gradient
                  fd-freq-gradient-large fd-gradient freq-isotope1 freq-isotope2 fnocc1 fnocc2
                  fnocc3 fnocc4 fnocc5 fnocc6 fnocc7 frac frac-ip-fitting frac-sym frac-traverse ghosts gibbs
                  lccd-grad1 lccd-grad2 matrix1 matrix2
                  mbis-1 mbis-2 mbis-3 mbis-4 mbis-5 mbis-6 mcscf1 mcscf2 mcscf3
                  mints1 mints2 mints3 mints4 mints5 mints6 mints8 mints-benchmark mints-helper
                  mints9 mints10 mints15 molden1 molden2 mom mom-h2o-3 mom-h2o-4
                  mp2-1  mp2-def2 mp2-grad1 mp2-grad2 mp2-h mp2p5-grad1 mp2p5-grad2 mp3-grad1 mp3-grad2
                  mp2-property mpn-bh nbody-he-cluster nbody-intermediates nbody-nocp-gradient
                  nbo nbody-cp-gradient nbody-vmfc-gradient nbody-vmfc-hessian nbody-hessian nbody-convergence
                  nbody-freq nbody-multi-level nbody-multi-level-2 numpy-array-interface
                  olccd-freq1 olccd-grad1 olccd-grad2 olccd1 olccd2 olccd3
                  omp2-1 omp2-2 omp2-3 omp2-4 omp2-5 omp2-grad1 omp2-grad2
                  omp2p5-1 omp2p5-2 omp2p5-grad1 omp2p5-grad2 omp3-1 omp3-2
                  omp3-3 omp3-4 omp3-5 omp3-grad1 omp3-grad2 opt-lindep-change
                  opt1 opt1-fd opt2 opt2-fd opt3 opt4 opt5 opt6 opt7 opt8 opt9
                  opt11 opt12 opt13 opt14 opt15 opt16 opt-irc-1 opt-irc-2 opt-irc-3
                  opt-freeze-coords opt-full-hess-every
                  oremp-grad1 oremp-grad2
                  phi-ao
                  props1 props2 props3 props4 psimrcc-ccsd_t-1 psimrcc-ccsd_t-2
                  psimrcc-ccsd_t-3 psimrcc-ccsd_t-4 psimrcc-fd-freq1
                  psimrcc-fd-freq2 psimrcc-pt2 psimrcc-sp1 psithon1 psithon2
                  pubchem1 pubchem2 pywrap-alias pywrap-all pywrap-basis
                  pywrap-cbs1 pywrap-checkrun-convcrit pywrap-checkrun-rhf
                  pywrap-checkrun-rohf pywrap-checkrun-uhf pywrap-db1
                  pywrap-db3
                  pywrap-molecule rasci-c2-active rasci-h2o
                  rasci-ne rasscf-sp sad-scf-type sad1 sapt1 sapt2 sapt3 sapt4 sapt5 sapt6 sapt-dft-api sapt-dft-lrc
                  remp-energy1 remp-energy2
                  sapt-exch-disp-inf sapt-exch-ind-inf sapt-exch-ind30-inf
                  sapt7 sapt8 scf-bz2 scf-dipder scf-guess scf-guess-read1 scf-upcast-custom-basis
                  scf-guess-read2 scf-guess-read3 scf1 scf-occ scf2 scf3 scf4 scf5 scf6
                  scf7 scf-level-shift-rks scf-level-shift-uhf scf-level-shift-cuhf scf-level-shift-rohf
                  scf-property soscf-large soscf-ref
                  soscf-dft stability1 dfep2-1 dfep2-2 sapt-dft1 sapt-dft2 sapt-compare sapt-sf1 dft-custom dft-reference
                  stability2 stability3 tu1-h2o-energy tu2-ch2-energy tu3-h2o-opt scf-response1 scf-response2 scf-response3
                  scf-cholesky-basis scf-auto-cholesky
                  tu4-h2o-freq tu5-sapt tu6-cp-ne2 x2c1 x2c2 x2c3 x2c-perturb-h zaptn-nh2
                  options1 cubeprop-esp dft-smoke scf-hess1 scf-hess2 scf-hess3 scf-hess4 scf-hess5 scf-freq1 dft-jk scf-coverage
                  dft-custom-dhdf dft-custom-hybrid dft-custom-mgga dft-custom-gga
                  pywrap-bfs pywrap-align pywrap-align-chiral mints12 cc-module
                  tdscf-1 tdscf-2 tdscf-3 tdscf-4 tdscf-5 tdscf-6 tdscf-7
                  dft-pruning freq-masses sapt9 sapt10 sapt11 scf-uhf-grad-nobeta
                  linK-1 linK-2 linK-3
                  cbs-xtpl-energy-conv ddd-deriv nbody-he-4b ddd-function-kwargs
                  )
    add_subdirectory(${test_name})
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

###########################  External Project Tests  ###########################

# For code that (1) is built or detected in the Psi4 superbuild but isn't
#   maintained in the Psi4 codebase and (2) has Psi4-dependent tests _not_ in
#   the external codebase, tests should
#   * get their own directory in tests (e.g., tests/dftd3)
#   * have a name that doesn't repeat the exe (e.g., tests/dftd3/energy, not tests/dftd3/dftd3-energy)
#   * be listed in the test case CMakeLists.txt as addon-testname (e.g., `add_regression_test(dftd3-energy ...`)
#   * alphabetize the addons

# <<<  BSE  >>>

if(ENABLE_bse)
    add_subdirectory(mints2-bse)
    message(STATUS "Adding test cases: Psi4 + bse")
endif()

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
             NAMES s-dftd3
             DOC "Path to Grimme's dftd3 executable")

if(DFTD3_EXECUTABLE)
    add_subdirectory(dftd3)
    add_subdirectory(sapt0-d)
    add_subdirectory(fsapt-d)
    add_subdirectory(fsaptd-terms)
    message(STATUS "${Cyan}Found DFTD3${ColourReset}: ${DFTD3_EXECUTABLE}")
    message(STATUS "Adding test cases: Psi4 + DFTD3")
endif()

# <<<  DFTD4  >>>

find_program(DFTD4_EXECUTABLE
             NAMES dftd4
             DOC "Path to Grimme's dftd4 executable")

if(DFTD4_EXECUTABLE)
    add_subdirectory(dftd4)
    message(STATUS "${Cyan}Found DFTD4${ColourReset}: ${DFTD4_EXECUTABLE}")
    message(STATUS "Adding test cases: Psi4 + DFTD4")
endif()

# <<<  DKH  >>>

if(ENABLE_dkh)
    add_subdirectory(dkh)
    message(STATUS "Adding test cases: Psi4 + dkh")
endif()

# <<<  LIBEFP  >>>

if(ENABLE_libefp)
    add_subdirectory(libefp)
    message(STATUS "Adding test cases: Psi4 + libefp")
endif()

# <<<  gCP  >>>

find_program(gCP_EXECUTABLE
             NAMES mctc-gcp
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

# <<<  V2RDM_CASSCF  >>>

if(ENABLE_v2rdm_casscf)
    add_subdirectory(v2rdm_casscf)
    message(STATUS "Adding test cases: Psi4 + v2rdm_casscf")
endif()

# <<<  GPU_DFCC  >>>

if(ENABLE_gpu_dfcc)
    add_subdirectory(gpu_dfcc)
    message(STATUS "Adding test cases: Psi4 + gpu_dfcc")
endif()

# <<<  SNSMP2  >>>

if(ENABLE_snsmp2)
    add_subdirectory(snsmp2)
    message(STATUS "Adding test cases: Psi4 + snsmp2")
endif()

# <<<  BrianQC  >>>
if(ENABLE_BrianQC)
    add_subdirectory(brianqc)
    message(STATUS "Adding test cases: Psi4 + BrianQC")
endif()

# <<<  ECPINT  >>>
if(ENABLE_ecpint)
    add_subdirectory(aediis-2)
    add_subdirectory(basis-ecp)
    add_subdirectory(dfmp2-ecp)
    add_subdirectory(sapt-ecp)
    add_subdirectory(scf-ecp)
    add_subdirectory(scf-ecp2)
    add_subdirectory(scf-ecp3)
    message(STATUS "Adding test cases: Psi4 + ecpint")
endif()

# plugin regime
if(ENABLE_PLUGIN_TESTING)
    add_subdirectory(plugins)
    message(STATUS "Adding test cases: Psi4 plugins")
endif()

# pytest regime
add_subdirectory(pytests)
