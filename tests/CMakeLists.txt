#These tests aren't set-up the same way as the others
set(TestsNotSetUpRight cfour optking testgrendel)

# In order to get CTest to run tests in separate categories,
# labels are used. All tests must carry the label "psi".
# 
# The following labels are defined:
#
#  longtests   = The (as of right now) 4 really long tests
#  longertests = Every test in the test folder less longtests
#  quicktests  = Proper subset of tests deemed to run quickly
#
# It is possible to specify more than one label per test. CMake
# expects a ;-separated list, e.g. "psi;quicktests;scf"
# It is mandatory that any test carries a length label!
# It is suggested to add at least a second label to specify the method tested.
#
#  We use ctest -I [start,end] where start is the test number we start with
#  and end is the test number we end with.  We will then use math to back out
#  the appropriate starts and ends...
#
#
#  In order to ensure that this works properly, please add your tests to the
#  appropriate variables given below

add_regression_test(adc1 "psi;longertests;adc")
add_regression_test(adc2 "psi;longertests;adc")

add_regression_test(castup1 "psi;quicktests;scf") 
add_regression_test(castup2 "psi;longertests;scf") 
add_regression_test(castup2 "psi;longertests;scf") 

add_regression_test(cc1 "psi;quicktests;cc")
add_regression_test(cc2 "psi;quicktests;cc")
add_regression_test(cc3 "psi;quicktests;cc")
add_regression_test(cc4 "psi;quicktests;cc")
add_regression_test(cc8 "psi;quicktests;cc")
add_regression_test(cc8a "psi;quicktests;cc")
add_regression_test(cc12 "psi;quicktests;cc")
add_regression_test(cc15 "psi;quicktests;cc")
add_regression_test(cc16 "psi;quicktests;cc")
add_regression_test(cc17 "psi;quicktests;cc")
add_regression_test(cc18 "psi;quicktests;cc")
add_regression_test(cc19 "psi;quicktests;cc")
add_regression_test(cc23 "psi;quicktests;cc")
add_regression_test(cc24 "psi;quicktests;cc")
add_regression_test(cc28 "psi;quicktests;cc")
add_regression_test(cc32 "psi;quicktests;cc")
add_regression_test(cc49 "psi;quicktests;cc")

add_regression_test(cc4a "psi;longertests;cc")
add_regression_test(cc5a "psi;longertests;cc")
add_regression_test(cc6 "psi;longertests;cc")
add_regression_test(cc8b "psi;longertests;cc")
add_regression_test(cc8c "psi;longertests;cc")
add_regression_test(cc8b "psi;longertests;cc")
add_regression_test(cc9 "psi;longertests;cc")
add_regression_test(cc9a "psi;longertests;cc")
add_regression_test(cc10 "psi;longertests;cc")
add_regression_test(cc11 "psi;longertests;cc")
add_regression_test(cc13 "psi;longertests;cc")
add_regression_test(cc13a "psi;longertests;cc")
add_regression_test(cc14 "psi;longertests;cc")
add_regression_test(cc21 "psi;longertests;cc")
add_regression_test(cc22 "psi;longertests;cc")
add_regression_test(cc25 "psi;longertests;cc")
add_regression_test(cc26 "psi;longertests;cc")
add_regression_test(cc27 "psi;longertests;cc")
add_regression_test(cc29 "psi;longertests;cc")
add_regression_test(cc30 "psi;longertests;cc")
add_regression_test(cc31 "psi;longertests;cc")
add_regression_test(cc33 "psi;longertests;cc")
add_regression_test(cc34 "psi;longertests;cc")
add_regression_test(cc35 "psi;longertests;cc")
add_regression_test(cc36 "psi;longertests;cc")
add_regression_test(cc37 "psi;longertests;cc")
add_regression_test(cc38 "psi;longertests;cc")
add_regression_test(cc39 "psi;longertests;cc")
add_regression_test(cc40 "psi;longertests;cc")
add_regression_test(cc41 "psi;longertests;cc")
add_regression_test(cc42 "psi;longertests;cc")
add_regression_test(cc43 "psi;longertests;cc")
add_regression_test(cc44 "psi;longertests;cc")
add_regression_test(cc45 "psi;longertests;cc")
add_regression_test(cc46 "psi;longertests;cc")
add_regression_test(cc47 "psi;longertests;cc")
add_regression_test(cc48 "psi;longertests;cc")
add_regression_test(cc49 "psi;longertests;cc")
add_regression_test(cc50 "psi;longertests;cc")
add_regression_test(cc51 "psi;longertests;cc")
add_regression_test(cc52 "psi;longertests;cc")
add_regression_test(cc53 "psi;longertests;cc")
add_regression_test(cc54 "psi;longertests;cc")
add_regression_test(cc55 "psi;longertests;cc")

add_regression_test(cdomp2-1 "psi;quicktests;cdomp2")
add_regression_test(cdomp2-2 "psi;quicktests;cdomp2")

add_regression_test(cepa0-grad1 "psi;longertests;cepa")
add_regression_test(cepa0-grad2 "psi;longertests;cepa")
add_regression_test(cepa1 "psi;longertests;cepa")
add_regression_test(cepa2 "psi;longertests;cepa")
add_regression_test(cepa3 "psi;longertests;cepa")

add_regression_test(cisd-h20+-0 "psi;quicktests;cisd")
add_regression_test(cisd-opt-fd "psi;quicktests;cisd")

add_regression_test(cisd-h20+-1 "psi;longertests;cisd")
add_regression_test(cisd-h20+-2 "psi;longertests;cisd")
add_regression_test(cisd-h20-clpse "psi;longertests;cisd")
add_regression_test(cisd-h20-sp "psi;longertests;cisd")
add_regression_test(cisd-h20-sp-2 "psi;longertests;cisd")

add_regression_test(dcft-grad1 "psi;quicktests;dcft")
add_regression_test(dcft1 "psi;quicktests;dcft")
add_regression_test(dcft2 "psi;quicktests;dcft")
add_regression_test(dcft3 "psi;quicktests;dcft")
add_regression_test(dcft4 "psi;quicktests;dcft")

add_regression_test(dcft5 "psi;longertests;dcft")
add_regression_test(dcft6 "psi;longertests;dcft")
add_regression_test(dcft7 "psi;longertests;dcft")

add_regression_test(dfmp2-1 "psi;quicktests;df;dfmp2")

add_regression_test(dfmp2-2 "psi;longertests;df;dfmp2")
add_regression_test(dfmp2-3 "psi;longertests;df;dfmp2")
add_regression_test(dfmp2-4 "psi;longertests;df;dfmp2")

add_regression_test(dfomp2-1 "psi;quicktests;df;dfomp2")
add_regression_test(dfomp2-2 "psi;quicktests;df;dfomp2")

add_regression_test(dfomp2-3 "psi;longertests;df;dfomp2")
add_regression_test(dfomp2-4 "psi;longertests;df;dfomp2")
add_regression_test(dfomp2-grad1 "psi;longertests;df;dfomp2")
add_regression_test(dfomp2-grad2 "psi;longertests;df;dfomp2")

add_regression_test(dfscf-bz2 "psi;longertests;df;dfscf")

add_regression_test(dft1 "psi;quicktests;dft;scf")

add_regression_test(dft1-alt "psi;longertests;dft;scf")
add_regression_test(dft2 "psi;longertests;dft;scf")
add_regression_test(dft3 "psi;longertests;dft;scf")
add_regression_test(dft-b2plyp "psi;longertests;dft;scf")
add_regression_test(dft-dldf "psi;longertests;dft;scf")
add_regression_test(dft-freq "psi;longertests;dft;scf")
add_regression_test(dft-grad "psi;longertests;dft;scf")
add_regression_test(dft-pbe0-2 "psi;longertests;dft;scf")
add_regression_test(dft-psivar "psi;longertests;dft;scf")

add_regression_test(docs-bases "psi;longertests;docs")
add_regression_test(docs-dft "psi;longertests;docs")
add_regression_test(docs-psimod "psi;longertests;docs")

add_regression_test(fci-dipole "psi;quicktests;fci")

add_regression_test(fci-h2o "psi;longertests;fci")
add_regression_test(fci-h2o-2 "psi;longertests;fci")
add_regression_test(fci-h2o-fzcv "psi;longertests;fci")
add_regression_test(fci-tdm "psi;longertests;fci")
add_regression_test(fci-tdm-2 "psi;longertests;fci")

add_regression_test(fd-freq-energy "psi;longertests;findif")
add_regression_test(fd-freq-energy-large "psi;longertests;findif")
add_regression_test(fd-freq-gradient "psi;longertests;findif")
add_regression_test(fd-gradient "psi;longertests;findif")

add_regression_test(fd-freq-gradient-large "psi;longtests;findif")

add_regression_test(fnocc2 "psi;quicktests;fnocc")
add_regression_test(fnocc4 "psi;quicktests;fnocc")

add_regression_test(fnocc1 "psi;longertests;fnocc")
add_regression_test(fnocc3 "psi;longertests;fnocc")

add_regression_test(mcscf1 "psi;quicktests;mcscf")

add_regression_test(mcscf2 "psi;longertests;mcscf")
add_regression_test(mcscf3 "psi;longertests;mcscf")

add_regression_test(mints1 "psi;quicktests;mints")
add_regression_test(mints2 "psi;quicktests;mints")

add_regression_test(mints3 "psi;longertests;mints")
add_regression_test(mints4 "psi;longertests;mints")
add_regression_test(mints5 "psi;longertests;mints")
add_regression_test(mints6 "psi;longertests;mints")
add_regression_test(mints8 "psi;longertests;mints")

add_regression_test(mp2_5-grad1 "psi;quicktests;mp2")
add_regression_test(mpn "psi;quicktests;mp2")

add_regression_test(mp2-1 "psi;longertests;mp2")
add_regression_test(mp2_5-grad2 "psi;longertests;mp2")
add_regression_test(mp2-grad1 "psi;longertests;mp2")
add_regression_test(mp2-grad2 "psi;longertests;mp2")
add_regression_test(mp2-def2 "psi;longtests;mp2")
add_regression_test(mp3-grad1 "psi;longertests;mp3")
add_regression_test(mp3-grad2 "psi;longertests;mp3")

add_regression_test(ocepa1 "psi;quicktests;ocepa")
add_regression_test(ocepa-grad2 "psi;quicktests;ocepa")

add_regression_test(ocepa2 "psi;longertests;ocepa")
add_regression_test(ocepa3 "psi;longertests;ocepa")
add_regression_test(ocepa-freq1 "psi;longertests;ocepa")
add_regression_test(ocepa-grad1 "psi;longertests;ocepa")

add_regression_test(omp2-2 "psi;quicktests;omp")
add_regression_test(omp3-1 "psi;quicktests;omp")
add_regression_test(omp2_5-2 "psi;quicktests;omp")
add_regression_test(omp2-grad1 "psi;quicktests;omp")

add_regression_test(omp2-1 "psi;longertests;omp")
add_regression_test(omp2-3 "psi;longertests;omp")
add_regression_test(omp2-4 "psi;longertests;omp")
add_regression_test(omp2_5-1 "psi;longertests;omp")
add_regression_test(omp2_5-grad1 "psi;longertests;omp")
add_regression_test(omp2_5-grad2 "psi;longertests;omp")
add_regression_test(omp2-grad2 "psi;longertests;omp")
add_regression_test(omp3-2 "psi;longertests;omp")
add_regression_test(omp3-3 "psi;longertests;omp")
add_regression_test(omp3-4 "psi;longertests;omp")
add_regression_test(omp3-5 "psi;longertests;omp")
add_regression_test(omp3-grad1 "psi;longertests;omp")
add_regression_test(omp3-grad2 "psi;longertests;omp")

add_regression_test(opt1 "psi;quicktests;opt")
add_regression_test(opt1-fd "psi;quicktests;opt")

add_regression_test(opt2 "psi;longertests;opt")
add_regression_test(opt2-fd "psi;longertests;opt")
add_regression_test(opt3 "psi;longertests;opt")
add_regression_test(opt4 "psi;longertests;opt")
add_regression_test(opt5 "psi;longertests;opt")
add_regression_test(opt6 "psi;longertests;opt")
add_regression_test(opt7 "psi;longertests;opt")
add_regression_test(opt-lindep-change "psi;longertests;opt")

add_regression_test(props1 "psi;longertests;prop")
add_regression_test(props2 "psi;longertests;prop")
add_regression_test(props3 "psi;longertests;prop")

add_regression_test(psimrcc-sp1 "psi;quicktests;psimrcc")

add_regression_test(psimrcc-ccsd_t-1 "psi;longertests;psimrcc")
add_regression_test(psimrcc-ccsd_t-2 "psi;longertests;psimrcc")
add_regression_test(psimrcc-ccsd_t-3 "psi;longertests;psimrcc")
add_regression_test(psimrcc-ccsd_t-4 "psi;longertests;psimrcc")
add_regression_test(psimrcc-fd-freq1 "psi;longertests;psimrcc")
add_regression_test(psimrcc-fd-freq2 "psi;longertests;psimrcc")
add_regression_test(psimrcc-pt2 "psi;longertests;psimrcc")

add_regression_test(pywrap-db1 "psi;quicktests;pywrap")
add_regression_test(pywrap-db2 "psi;quicktests;pywrap")
add_regression_test(pywrap-cbs1 "psi;quicktests;pywrap")
add_regression_test(pywrap-all "psi;quicktests;pywrap")
add_regression_test(pywrap-alias "psi;quicktests;pywrap")
add_regression_test(pywrap-molecule "psi;quicktests;pywrap")
add_regression_test(pywrap-freq-e-sowreap "psi;quicktests;pywrap")

add_regression_test(pywrap-basis "psi;longertests;pywrap")
add_regression_test(pywrap-checkrun-convcrit "psi;longertests;pywrap")
add_regression_test(pywrap-checkrun-rhf "psi;longertests;pywrap")
add_regression_test(pywrap-checkrun-rohf "psi;longertests;pywrap")
add_regression_test(pywrap-checkrun-uhf "psi;longertests;pywrap")
add_regression_test(pywrap-db2 "psi;longertests;pywrap")
add_regression_test(pywrap-opt-sowreap "psi;longertests;pywrap")

add_regression_test(rasci-h2o "psi;quicktests;rasci")

add_regression_test(rasci-c2-active "psi;longertests;rasci")
add_regression_test(rasci-ne "psi;longertests;rasci")

add_regression_test(sapt1 "psi;quicktests;sapt")

add_regression_test(sapt3 "psi;longertests;sapt")
add_regression_test(sapt5 "psi;longertests;sapt")

add_regression_test(sapt2 "psi;longtests;sapt")

add_regression_test(scf1 "psi;quicktests;scf")
add_regression_test(scf2 "psi;quicktests;scf")
add_regression_test(scf3 "psi;quicktests;scf")
add_regression_test(scf5 "psi;quicktests;scf")
add_regression_test(scf-guess-read "psi;quicktests;scf")

add_regression_test(scf11-freq-from-energies "psi;longertests;scf")
add_regression_test(scf4 "psi;longertests;scf")
add_regression_test(scf6 "psi;longertests;scf")
add_regression_test(scf-bz2 "psi;longertests;scf")

add_regression_test(tu1-h2o-energy "psi;quicktests;tutorial")
add_regression_test(tu2-ch2-energy "psi;quicktests;tutorial")
add_regression_test(tu3-h2o-opt "psi;quicktests;tutorial;opt")
add_regression_test(tu4-h2o-freq "psi;quicktests;tutorial;freq")
add_regression_test(tu5-sapt "psi;quicktests;tutorial;sapt")

add_regression_test(tu6-cp-ne2 "psi;longertests;tutorial")

add_regression_test(ci-multi "psi;quicktests;misc")
add_regression_test(sad1 "psi;quicktests;misc")
add_regression_test(matrix1 "psi;quicktests;misc")
add_regression_test(gibbs "psi;quicktests;misc")

add_regression_test(frac "psi;longertests;misc")
add_regression_test(ghosts "psi;longertests;misc")
add_regression_test(large_atoms "psi;longertests;misc")
add_regression_test(min_input "psi;longertests;misc")
add_regression_test(mom "psi;longertests;misc")
add_regression_test(psithon1 "psi;longertests;misc")
add_regression_test(psithon2 "psi;longertests;misc")
add_regression_test(pubchem1 "psi;longertests;misc")
add_regression_test(zaptn-nh2 "psi;longertests;misc")

add_regression_test(casscf-fzc-sp "psi;quicktests;casscf")
add_regression_test(casscf-sp "psi;quicktests;casscf")

set(AutoTestSubjects cc8 cc9 cc11 cc18 cc19 cc28 cc29 cc30 cc31 cc34 cc36 cc37 cc38 cc39 cc40 cc41 cc42 cc43 cc46 cc47 cc48 cc49 cc50 cc51 cc52 cc55)

###### END test definitions ################

###### START external-executable-dependent tests ################

# Set up groups of test cases analogous to dftd3 below
#   In particular, new ext-exe-dep tests should 
#   * get their own directory in tests (e.g., tests/dftd3)
#   * have a name that doesn't repeat the exe (e.g., test/dftd3/energy)
#   * be listed in the exeShort and exeLong variable name as exe-name (e.g., dftd3-energy)
#   * be sure the exe is listed in the python code below (search 'partition')

# <<<  DFTD3  >>>

find_program(DFTD3_EXECUTABLE
             NAMES dftd3
             DOC "Path to Grimme's dftd3 executable")

if(DFTD3_EXECUTABLE)
   add_regression_test(dftd3-energy "psi;quicktests;dftd;ext-exe-dep")
   add_regression_test(dftd3-grad "psi;longertests;dftd;ext-exe-dep")
   add_regression_test(dftd3-psithon2 "psi;longertests;dftd;ext-exe-dep")
   message(STATUS "dftd3 executable for interfaced DFTD3 found. Test cases added.")
endif()

# <<<  MRCC  >>>

set(MRCCShort mrcc-ccsdt mrcc-ccsd_t_)
set(MRCCLong mrcc-ccsdt_q_ mrcc-optfreq)

find_program(MRCC_EXECUTABLE
             NAMES dmrcc
             DOC "Path to Kallay's mrcc executable")

if(MRCC_EXECUTABLE)
   add_regression_test(mrcc-ccsdt "psi;quicktests;mrcc;cc;ext-exe-dep")
   add_regression_test(mrcc-ccsd_t_ "psi;quicktests;mrcc;cc;ext-exe-dep")
   add_regression_test(mrcc-ccsdt_q_ "psi;longertests;mrcc;cc;ext-exe-dep")
   add_regression_test(mrcc-optfreq "psi;longertests;mrcc;cc;ext-exe-dep")
   message(STATUS "dmrcc executable for interfaced MRCC found. Test cases added.")
endif()

# <<<  CFOUR  >>>

find_program(CFOUR_EXECUTABLE
             NAMES xcfour
             DOC "Path to Stanton and Gauss's xcfour executable")

if(CFOUR_EXECUTABLE)
   add_regression_test(cfour-sp-uhf-mp2 "psi;quicktests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-scf "psi;quicktests;cfour;scf;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-ccsd_t_-ecc "psi;quicktests;cfour;ext-exe-dep")
   add_regression_test(cfour-opt-rhf-scf "psi;quicktests;cfour;ext-exe-dep")
   add_regression_test(cfour-psi-rhf-scsmp2 "psi;quicktests;cfour;ext-exe-dep")
   add_regression_test(cfour-pywrap-cbs1 "psi;quicktests;cfour;ext-exe-dep")

   add_regression_test(cfour-puream "psi;longertests;cfour;ext-exe-dep")            
   add_regression_test(cfour-mints5 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-mints5-grad "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-kw-1 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-kw-2 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-kw-3 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-kw-4 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-kw-5 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-kw-6 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-kw-7 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-kw-8 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-psi-ghost-grad "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-scf4 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-mp2-1 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-pywrap-db1 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-pywrap-basis "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-dfmp2-1 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-scf "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-mp2 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-ccsd "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-ccsd-ao "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-ccsd_t_ "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-ccsd_t_-ao "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-ccsd_t_-ao-ecc "psi;longertests;cfour;ext-exe-dep") 
   add_regression_test(cfour-sp-rhf-ccsd_t_-fc "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-cc3 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-ccsdt "psi;longertests;cfour;ext-exe-dep") 
   add_regression_test(cfour-sp-uhf-scf "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-uhf-ccsd "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-uhf-ccsd_t_ "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-uhf-ccsd_t_-ao "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-uhf-ccsd_t_-ao-ecc "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-uhf-ccsd_t_-ecc "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-uhf-cc3 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-uhf-ccsdt "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-mp2-sc "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-mp4-sc "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-ccsd "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-ccsd_t_ "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-ccsdt "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-ccsd-ao "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-ccsd_t_-ao "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-ccsd_t_-ao-ecc "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-ccsd_t_-fc "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-opt-rhf-mp2 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-opt-rhf-ccsd_t_ "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-opt-rhf-ccsd_t_-ecc "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-psi-a24-grad "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-psi-uhf-scsmp2 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-psi-rohf-scsmp2 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-psi-rhf-mp3 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-psi-uhf-mp3 "psi;longertests;cfour;ext-exe-dep")
   add_regression_test(cfour-psi-mp4 "psi;longertests;cfour;ext-exe-dep")
   message(STATUS "xcfour executable for interfaced CFOUR found. Test cases added.")
endif()

###### END external-executable-dependent tests ################

#set(TestSet ${longtests} ${quicktests} ${longertests})
#
#list(LENGTH TestSet NumberOfTests)
#list(LENGTH longtests NLong)
#list(LENGTH quicktests NShort)
#math(EXPR Start "${NLong} +1")
#math(EXPR Stop "${NLong} + ${NShort}")
#
#add_custom_target(longtests ${CMAKE_CTEST_COMMAND} -I 1,${NLong})
#add_custom_target(quicktests ${CMAKE_CTEST_COMMAND} -I ${Start},${Stop})
#add_custom_target(tests ${CMAKE_CTEST_COMMAND} -I ${Start},${NumberOfTests})


if(NOT PERL_FOUND)
   message(WARNING "Perl was not found. Some Coupled Cluster test cases will not be thoroughly checked.")
endif()

# Declare each test as a possible test target even though no one will likely
#   want to run individual tests. In the end, makes our job easier.
set(index 1)
foreach(Tests2Run ${TestSet})

    # Check if a test case from subdirectory (ext-exec-dep)
    execute_process(COMMAND ${PYTHON_EXECUTABLE} -c 
        "pieces = '${Tests2Run}'.partition('-'); print '%s;%s' % (pieces[0], pieces[2]) if pieces[0] in ['dftd3', 'mrcc', 'cfour'] else ''"
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
endforeach(Tests2Run)
