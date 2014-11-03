include(TestingMacros)
#These tests aren't set-up the same way as the others
set(TestsNotSetUpRight cfour optking testgrendel)

# In order to get CTest to run tests in separate categories,
# labels are used. All tests must carry the label "psi".
# 
# The following labels are defined:
#
#  LongTests   = The (as of right now) 4 really long tests
#  LongerTests = Every test in the test folder less LongTests
#  QuickTests  = Proper subset of tests deemed to run quickly
#
# It is possible to specify more than one label per test. CMake
# expects a ;-separated list, e.g. "psi;QuickTests;scf"
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

add_regression_test(adc1 "psi;LongerTests;adc")
add_regression_test(adc2 "psi;LongerTests;adc")

add_regression_test(castup1 "psi;QuickTests;scf") 
add_regression_test(castup2 "psi;LongerTests;scf") 
add_regression_test(castup2 "psi;LongerTests;scf") 

add_regression_test(cc1 "psi;QuickTests;cc")
add_regression_test(cc2 "psi;QuickTests;cc")
add_regression_test(cc3 "psi;QuickTests;cc")
add_regression_test(cc4 "psi;QuickTests;cc")
add_regression_test(cc8 "psi;QuickTests;cc")
add_regression_test(cc8a "psi;QuickTests;cc")
add_regression_test(cc12 "psi;QuickTests;cc")
add_regression_test(cc15 "psi;QuickTests;cc")
add_regression_test(cc16 "psi;QuickTests;cc")
add_regression_test(cc17 "psi;QuickTests;cc")
add_regression_test(cc18 "psi;QuickTests;cc")
add_regression_test(cc19 "psi;QuickTests;cc")
add_regression_test(cc23 "psi;QuickTests;cc")
add_regression_test(cc24 "psi;QuickTests;cc")
add_regression_test(cc28 "psi;QuickTests;cc")
add_regression_test(cc32 "psi;QuickTests;cc")
add_regression_test(cc49 "psi;QuickTests;cc")

add_regression_test(cc4a "psi;LongerTests;cc")
add_regression_test(cc5a "psi;LongerTests;cc")
add_regression_test(cc6 "psi;LongerTests;cc")
add_regression_test(cc8b "psi;LongerTests;cc")
add_regression_test(cc8c "psi;LongerTests;cc")
add_regression_test(cc8b "psi;LongerTests;cc")
add_regression_test(cc9 "psi;LongerTests;cc")
add_regression_test(cc9a "psi;LongerTests;cc")
add_regression_test(cc10 "psi;LongerTests;cc")
add_regression_test(cc11 "psi;LongerTests;cc")
add_regression_test(cc13 "psi;LongerTests;cc")
add_regression_test(cc13a "psi;LongerTests;cc")
add_regression_test(cc14 "psi;LongerTests;cc")
add_regression_test(cc21 "psi;LongerTests;cc")
add_regression_test(cc22 "psi;LongerTests;cc")
add_regression_test(cc25 "psi;LongerTests;cc")
add_regression_test(cc26 "psi;LongerTests;cc")
add_regression_test(cc27 "psi;LongerTests;cc")
add_regression_test(cc29 "psi;LongerTests;cc")
add_regression_test(cc30 "psi;LongerTests;cc")
add_regression_test(cc31 "psi;LongerTests;cc")
add_regression_test(cc33 "psi;LongerTests;cc")
add_regression_test(cc34 "psi;LongerTests;cc")
add_regression_test(cc35 "psi;LongerTests;cc")
add_regression_test(cc36 "psi;LongerTests;cc")
add_regression_test(cc37 "psi;LongerTests;cc")
add_regression_test(cc38 "psi;LongerTests;cc")
add_regression_test(cc39 "psi;LongerTests;cc")
add_regression_test(cc40 "psi;LongerTests;cc")
add_regression_test(cc41 "psi;LongerTests;cc")
add_regression_test(cc42 "psi;LongerTests;cc")
add_regression_test(cc43 "psi;LongerTests;cc")
add_regression_test(cc44 "psi;LongerTests;cc")
add_regression_test(cc45 "psi;LongerTests;cc")
add_regression_test(cc46 "psi;LongerTests;cc")
add_regression_test(cc47 "psi;LongerTests;cc")
add_regression_test(cc48 "psi;LongerTests;cc")
add_regression_test(cc49 "psi;LongerTests;cc")
add_regression_test(cc50 "psi;LongerTests;cc")
add_regression_test(cc51 "psi;LongerTests;cc")
add_regression_test(cc52 "psi;LongerTests;cc")
add_regression_test(cc53 "psi;LongerTests;cc")
add_regression_test(cc54 "psi;LongerTests;cc")
add_regression_test(cc55 "psi;LongerTests;cc")

add_regression_test(cdomp2-1 "psi;QuickTests;cdomp2")
add_regression_test(cdomp2-2 "psi;QuickTests;cdomp2")

add_regression_test(cepa0-grad1 "psi;LongerTests;cepa")
add_regression_test(cepa0-grad2 "psi;LongerTests;cepa")
add_regression_test(cepa1 "psi;LongerTests;cepa")
add_regression_test(cepa2 "psi;LongerTests;cepa")
add_regression_test(cepa3 "psi;LongerTests;cepa")

add_regression_test(cisd-h20+-0 "psi;QuickTests;cisd")
add_regression_test(cisd-opt-fd "psi;QuickTests;cisd")

add_regression_test(cisd-h20+-1 "psi;LongerTests;cisd")
add_regression_test(cisd-h20+-2 "psi;LongerTests;cisd")
add_regression_test(cisd-h20-clpse "psi;LongerTests;cisd")
add_regression_test(cisd-h20-sp "psi;LongerTests;cisd")
add_regression_test(cisd-h20-sp-2 "psi;LongerTests;cisd")

add_regression_test(dcft-grad1 "psi;QuickTests;dcft")
add_regression_test(dcft1 "psi;QuickTests;dcft")
add_regression_test(dcft2 "psi;QuickTests;dcft")
add_regression_test(dcft3 "psi;QuickTests;dcft")
add_regression_test(dcft4 "psi;QuickTests;dcft")

add_regression_test(dcft5 "psi;LongerTests;dcft")
add_regression_test(dcft6 "psi;LongerTests;dcft")
add_regression_test(dcft7 "psi;LongerTests;dcft")

add_regression_test(dfmp2-1 "psi;QuickTests;df;dfmp2")

add_regression_test(dfmp2-2 "psi;LongerTests;df;dfmp2")
add_regression_test(dfmp2-3 "psi;LongerTests;df;dfmp2")
add_regression_test(dfmp2-4 "psi;LongerTests;df;dfmp2")

add_regression_test(dfomp2-1 "psi;QuickTests;df;dfomp2")
add_regression_test(dfomp2-2 "psi;QuickTests;df;dfomp2")

add_regression_test(dfomp2-3 "psi;LongerTests;df;dfomp2")
add_regression_test(dfomp2-4 "psi;LongerTests;df;dfomp2")
add_regression_test(dfomp2-grad1 "psi;LongerTests;df;dfomp2")
add_regression_test(dfomp2-grad2 "psi;LongerTests;df;dfomp2")

add_regression_test(dfscf-bz2 "psi;LongerTests;df;dfscf")

add_regression_test(dft1 "psi;QuickTests;dft;scf")

add_regression_test(dft1-alt "psi;LongerTests;dft;scf")
add_regression_test(dft2 "psi;LongerTests;dft;scf")
add_regression_test(dft3 "psi;LongerTests;dft;scf")
add_regression_test(dft-b2plyp "psi;LongerTests;dft;scf")
add_regression_test(dft-dldf "psi;LongerTests;dft;scf")
add_regression_test(dft-freq "psi;LongerTests;dft;scf")
add_regression_test(dft-grad "psi;LongerTests;dft;scf")
add_regression_test(dft-pbe0-2 "psi;LongerTests;dft;scf")
add_regression_test(dft-psivar "psi;LongerTests;dft;scf")

add_regression_test(docs-bases "psi;LongerTests;docs")
add_regression_test(docs-dft "psi;LongerTests;docs")
add_regression_test(docs-psimod "psi;LongerTests;docs")

add_regression_test(fci-dipole "psi;QuickTests;fci")

add_regression_test(fci-h2o "psi;LongerTests;fci")
add_regression_test(fci-h2o-2 "psi;LongerTests;fci")
add_regression_test(fci-h2o-fzcv "psi;LongerTests;fci")
add_regression_test(fci-tdm "psi;LongerTests;fci")
add_regression_test(fci-tdm-2 "psi;LongerTests;fci")

add_regression_test(fd-freq-energy "psi;LongerTests;findif")
add_regression_test(fd-freq-energy-large "psi;LongerTests;findif")
add_regression_test(fd-freq-gradient "psi;LongerTests;findif")
add_regression_test(fd-gradient "psi;LongerTests;findif")

add_regression_test(fd-freq-gradient-large "psi;LongTests;findif")

add_regression_test(fnocc2 "psi;QuickTests;fnocc")
add_regression_test(fnocc4 "psi;QuickTests;fnocc")

add_regression_test(fnocc1 "psi;LongerTests;fnocc")
add_regression_test(fnocc3 "psi;LongerTests;fnocc")

add_regression_test(mcscf1 "psi;QuickTests;mcscf")

add_regression_test(mcscf2 "psi;LongerTests;mcscf")
add_regression_test(mcscf3 "psi;LongerTests;mcscf")

add_regression_test(mints1 "psi;QuickTests;mints")
add_regression_test(mints2 "psi;QuickTests;mints")

add_regression_test(mints3 "psi;LongerTests;mints")
add_regression_test(mints4 "psi;LongerTests;mints")
add_regression_test(mints5 "psi;LongerTests;mints")
add_regression_test(mints6 "psi;LongerTests;mints")
add_regression_test(mints8 "psi;LongerTests;mints")

add_regression_test(mp2_5-grad1 "psi;QuickTests;mp2")
add_regression_test(mpn "psi;QuickTests;mp2")

add_regression_test(mp2-1 "psi;LongerTests;mp2")
add_regression_test(mp2_5-grad2 "psi;LongerTests;mp2")
add_regression_test(mp2-grad1 "psi;LongerTests;mp2")
add_regression_test(mp2-grad2 "psi;LongerTests;mp2")
add_regression_test(mp2-def2 "psi;LongTests;mp2")
add_regression_test(mp3-grad1 "psi;LongerTests;mp3")
add_regression_test(mp3-grad2 "psi;LongerTests;mp3")

add_regression_test(ocepa1 "psi;QuickTests;ocepa")
add_regression_test(ocepa-grad2 "psi;QuickTests;ocepa")

add_regression_test(ocepa2 "psi;LongerTests;ocepa")
add_regression_test(ocepa3 "psi;LongerTests;ocepa")
add_regression_test(ocepa-freq1 "psi;LongerTests;ocepa")
add_regression_test(ocepa-grad1 "psi;LongerTests;ocepa")

add_regression_test(omp2-2 "psi;QuickTests;omp")
add_regression_test(omp3-1 "psi;QuickTests;omp")
add_regression_test(omp2_5-2 "psi;QuickTests;omp")
add_regression_test(omp2-grad1 "psi;QuickTests;omp")

add_regression_test(omp2-1 "psi;LongerTests;omp")
add_regression_test(omp2-3 "psi;LongerTests;omp")
add_regression_test(omp2-4 "psi;LongerTests;omp")
add_regression_test(omp2_5-1 "psi;LongerTests;omp")
add_regression_test(omp2_5-grad1 "psi;LongerTests;omp")
add_regression_test(omp2_5-grad2 "psi;LongerTests;omp")
add_regression_test(omp2-grad2 "psi;LongerTests;omp")
add_regression_test(omp3-2 "psi;LongerTests;omp")
add_regression_test(omp3-3 "psi;LongerTests;omp")
add_regression_test(omp3-4 "psi;LongerTests;omp")
add_regression_test(omp3-5 "psi;LongerTests;omp")
add_regression_test(omp3-grad1 "psi;LongerTests;omp")
add_regression_test(omp3-grad2 "psi;LongerTests;omp")

add_regression_test(opt1 "psi;QuickTests;opt")
add_regression_test(opt1-fd "psi;QuickTests;opt")

add_regression_test(opt2 "psi;LongerTests;opt")
add_regression_test(opt2-fd "psi;LongerTests;opt")
add_regression_test(opt3 "psi;LongerTests;opt")
add_regression_test(opt4 "psi;LongerTests;opt")
add_regression_test(opt5 "psi;LongerTests;opt")
add_regression_test(opt6 "psi;LongerTests;opt")
add_regression_test(opt7 "psi;LongerTests;opt")
add_regression_test(opt-lindep-change "psi;LongerTests;opt")

add_regression_test(props1 "psi;LongerTests;prop")
add_regression_test(props2 "psi;LongerTests;prop")
add_regression_test(props3 "psi;LongerTests;prop")

add_regression_test(psimrcc-sp1 "psi;QuickTests;psimrcc")

add_regression_test(psimrcc-ccsd_t-1 "psi;LongerTests;psimrcc")
add_regression_test(psimrcc-ccsd_t-2 "psi;LongerTests;psimrcc")
add_regression_test(psimrcc-ccsd_t-3 "psi;LongerTests;psimrcc")
add_regression_test(psimrcc-ccsd_t-4 "psi;LongerTests;psimrcc")
add_regression_test(psimrcc-fd-freq1 "psi;LongerTests;psimrcc")
add_regression_test(psimrcc-fd-freq2 "psi;LongerTests;psimrcc")
add_regression_test(psimrcc-pt2 "psi;LongerTests;psimrcc")

add_regression_test(pywrap-db1 "psi;QuickTests;pywrap")
add_regression_test(pywrap-db2 "psi;QuickTests;pywrap")
add_regression_test(pywrap-cbs1 "psi;QuickTests;pywrap")
add_regression_test(pywrap-all "psi;QuickTests;pywrap")
add_regression_test(pywrap-alias "psi;QuickTests;pywrap")
add_regression_test(pywrap-molecule "psi;QuickTests;pywrap")
add_regression_test(pywrap-freq-e-sowreap "psi;QuickTests;pywrap")

add_regression_test(pywrap-basis "psi;LongerTests;pywrap")
add_regression_test(pywrap-checkrun-convcrit "psi;LongerTests;pywrap")
add_regression_test(pywrap-checkrun-rhf "psi;LongerTests;pywrap")
add_regression_test(pywrap-checkrun-rohf "psi;LongerTests;pywrap")
add_regression_test(pywrap-checkrun-uhf "psi;LongerTests;pywrap")
add_regression_test(pywrap-db2 "psi;LongerTests;pywrap")
add_regression_test(pywrap-opt-sowreap "psi;LongerTests;pywrap")

add_regression_test(rasci-h2o "psi;QuickTests;rasci")

add_regression_test(rasci-c2-active "psi;LongerTests;rasci")
add_regression_test(rasci-ne "psi;LongerTests;rasci")

add_regression_test(sapt1 "psi;QuickTests;sapt")

add_regression_test(sapt3 "psi;LongerTests;sapt")
add_regression_test(sapt5 "psi;LongerTests;sapt")

add_regression_test(sapt2 "psi;LongTests;sapt")

add_regression_test(scf1 "psi;QuickTests;scf")
add_regression_test(scf2 "psi;QuickTests;scf")
add_regression_test(scf3 "psi;QuickTests;scf")
add_regression_test(scf5 "psi;QuickTests;scf")
add_regression_test(scf-guess-read "psi;QuickTests;scf")

add_regression_test(scf11-freq-from-energies "psi;LongerTests;scf")
add_regression_test(scf4 "psi;LongerTests;scf")
add_regression_test(scf6 "psi;LongerTests;scf")
add_regression_test(scf-bz2 "psi;LongerTests;scf")

add_regression_test(tu1-h2o-energy "psi;QuickTests;tutorial")
add_regression_test(tu2-ch2-energy "psi;QuickTests;tutorial")
add_regression_test(tu3-h2o-opt "psi;QuickTests;tutorial;opt")
add_regression_test(tu4-h2o-freq "psi;QuickTests;tutorial;freq")
add_regression_test(tu5-sapt "psi;QuickTests;tutorial;sapt")

add_regression_test(tu6-cp-ne2 "psi;LongerTests;tutorial")

add_regression_test(ci-multi "psi;QuickTests;misc")
add_regression_test(sad1 "psi;QuickTests;misc")
add_regression_test(matrix1 "psi;QuickTests;misc")
add_regression_test(gibbs "psi;QuickTests;misc")

add_regression_test(frac "psi;LongerTests;misc")
add_regression_test(ghosts "psi;LongerTests;misc")
add_regression_test(large_atoms "psi;LongerTests;misc")
add_regression_test(min_input "psi;LongerTests;misc")
add_regression_test(mom "psi;LongerTests;misc")
add_regression_test(psithon1 "psi;LongerTests;misc")
add_regression_test(psithon2 "psi;LongerTests;misc")
add_regression_test(pubchem1 "psi;LongerTests;misc")
add_regression_test(zaptn-nh2 "psi;LongerTests;misc")

add_regression_test(casscf-fzc-sp "psi;QuickTests;casscf")
add_regression_test(casscf-sp "psi;QuickTests;casscf")

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
   add_regression_test(dftd3-energy "psi;QuickTests;dftd;ext-exe-dep")
   add_regression_test(dftd3-grad "psi;LongerTests;dftd;ext-exe-dep")
   add_regression_test(dftd3-psithon2 "psi;LongerTests;dftd;ext-exe-dep")
   message(STATUS "dftd3 executable for interfaced DFTD3 found. Test cases added.")
endif()

# <<<  MRCC  >>>

set(MRCCShort mrcc-ccsdt mrcc-ccsd_t_)
set(MRCCLong mrcc-ccsdt_q_ mrcc-optfreq)

find_program(MRCC_EXECUTABLE
             NAMES dmrcc
             DOC "Path to Kallay's mrcc executable")

if(MRCC_EXECUTABLE)
   add_regression_test(mrcc-ccsdt "psi;QuickTests;mrcc;cc;ext-exe-dep")
   add_regression_test(mrcc-ccsd_t_ "psi;QuickTests;mrcc;cc;ext-exe-dep")
   add_regression_test(mrcc-ccsdt_q_ "psi;LongerTests;mrcc;cc;ext-exe-dep")
   add_regression_test(mrcc-optfreq "psi;LongerTests;mrcc;cc;ext-exe-dep")
   message(STATUS "dmrcc executable for interfaced MRCC found. Test cases added.")
endif()

# <<<  CFOUR  >>>

find_program(CFOUR_EXECUTABLE
             NAMES xcfour
             DOC "Path to Stanton and Gauss's xcfour executable")

if(CFOUR_EXECUTABLE)
   add_regression_test(cfour-sp-uhf-mp2 "psi;QuickTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-scf "psi;QuickTests;cfour;scf;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-ccsd_t_-ecc "psi;QuickTests;cfour;ext-exe-dep")
   add_regression_test(cfour-opt-rhf-scf "psi;QuickTests;cfour;ext-exe-dep")
   add_regression_test(cfour-psi-rhf-scsmp2 "psi;QuickTests;cfour;ext-exe-dep")
   add_regression_test(cfour-pywrap-cbs1 "psi;QuickTests;cfour;ext-exe-dep")

   add_regression_test(cfour-puream "psi;LongerTests;cfour;ext-exe-dep")            
   add_regression_test(cfour-mints5 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-mints5-grad "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-kw-1 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-kw-2 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-kw-3 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-kw-4 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-kw-5 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-kw-6 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-kw-7 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-kw-8 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-psi-ghost-grad "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-scf4 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-mp2-1 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-pywrap-db1 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-pywrap-basis "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-dfmp2-1 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-scf "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-mp2 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-ccsd "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-ccsd-ao "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-ccsd_t_ "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-ccsd_t_-ao "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-ccsd_t_-ao-ecc "psi;LongerTests;cfour;ext-exe-dep") 
   add_regression_test(cfour-sp-rhf-ccsd_t_-fc "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-cc3 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rhf-ccsdt "psi;LongerTests;cfour;ext-exe-dep") 
   add_regression_test(cfour-sp-uhf-scf "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-uhf-ccsd "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-uhf-ccsd_t_ "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-uhf-ccsd_t_-ao "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-uhf-ccsd_t_-ao-ecc "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-uhf-ccsd_t_-ecc "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-uhf-cc3 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-uhf-ccsdt "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-mp2-sc "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-mp4-sc "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-ccsd "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-ccsd_t_ "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-ccsdt "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-ccsd-ao "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-ccsd_t_-ao "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-ccsd_t_-ao-ecc "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-sp-rohf-ccsd_t_-fc "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-opt-rhf-mp2 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-opt-rhf-ccsd_t_ "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-opt-rhf-ccsd_t_-ecc "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-psi-a24-grad "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-psi-uhf-scsmp2 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-psi-rohf-scsmp2 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-psi-rhf-mp3 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-psi-uhf-mp3 "psi;LongerTests;cfour;ext-exe-dep")
   add_regression_test(cfour-psi-mp4 "psi;LongerTests;cfour;ext-exe-dep")
   message(STATUS "xcfour executable for interfaced CFOUR found. Test cases added.")
endif()

###### END external-executable-dependent tests ################

#set(TestSet ${LongTests} ${QuickTests} ${LongerTests})
#
#list(LENGTH TestSet NumberOfTests)
#list(LENGTH LongTests NLong)
#list(LENGTH QuickTests NShort)
#math(EXPR Start "${NLong} +1")
#math(EXPR Stop "${NLong} + ${NShort}")
#
#add_custom_target(longtests ${CMAKE_CTEST_COMMAND} -I 1,${NLong})
#add_custom_target(quicktests ${CMAKE_CTEST_COMMAND} -I ${Start},${Stop})
#add_custom_target(tests ${CMAKE_CTEST_COMMAND} -I ${Start},${NumberOfTests})


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
