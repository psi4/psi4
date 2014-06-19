#These tests aren't set-up the same way as the others
set(TestsNotSetUpRight dftd3 cfour optking testgrendel)

#No ADC Short
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
set(FDLong fd-freq-energy fd-freq-energy-large fd-freq-gradient fd-freq-gradient-large fd-gradient)
set(FNOCCShort fnocc2 fnocc4)
set(FNOCCLong fnocc1 fnocc3) 

set(MCSCFShort mcscf1)
set(MCSCFLong mcscf2 mcscf3)
set(MintsShort mints1 mints2)
set(MintsLong mints3 mints4 mints5 mints6 mints8)
set(MP2Short mp2_5-grad1 mpn-bh)
set(MP2Long mp2-1 mp2_5-grad2 mp2-def2 mp2-grad1 mp2-grad2 mp3-grad1 mp3-grad2)
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
set(SAPTLong sapt2 sapt3 sapt4 sapt5)
set(SCFShort scf1 scf2 scf3 scf5 scf-guess-read)
set(SCFLong scf11-freq-from-energies scf4 scf6 scf-bz2)

set(TutorialShort tu1-h2o-energy tu2-ch2-energy tu3-h2o-opt tu4-h2o-freq tu5-sapt )
set(TutorialLong tu6-cp-ne2)
set(MISCShort ci-multi sad1 matrix1 gibbs)
set(MISCLong frac ghosts large_atoms min_input mom psithon1 psithon2 pubchem1 zaptn-nh2)

set(QuickTests ${CastUpShort} ${CCShort} ${CDShort} ${CISDShort} ${DCFTShort}
    ${DFShort} ${DFTShort} ${FCIShort} ${FNOCCShort} ${MCSCFShort}
${MintsShort} ${MP2Short} ${OCEPAShort} ${OMPShort} ${OPTShort}
${PSIMRCCShort} ${PYWrapShort} ${RASCIShort} ${SAPTShort} ${SCFShort}
${TutorialShort} ${MISCShort}
)

set(LongTests ${ADCLong} ${CastUpLong} ${CCLong} ${CEPALong} ${CISDLong}
${DCFTLong} ${DFLong} ${DFTLong} ${DocsLong} ${FCILong} ${FDLong} ${FNOCCLong}
${MCSCFLong} ${MintsLong} ${MP2Long} ${MRCCLong} ${OCEPALong} ${OMPLong}
${OPTLong} ${PropsLong} ${PSIMRCCLong} ${PYWrapLong} ${RASCILong} ${SAPTLong}
${SCFLong} ${TutorialLong} ${MISCLong})

set(AllTests ${QuickTests} ${LongTests})
set(TestSet ${QuickTests})

#This is the command to actually run the tests
set(PSIEXE ${CMAKE_BINARY_DIR}/bin/psi4)
set(TESTDIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TESTEXE ${TESTDIR}/runtest.pl)

#CMake version of the Make runtest.pl.in
file(WRITE ${TESTEXE}
"use strict;
use warnings;

die \"Usage\\n\\t\\t\$0 input_file logfile [output_file]\\n\" unless (@ARGV == 3) || (@ARGV == 4) || (@ARGV == 5);

my \$psi = \"../../bin/psi4\";
my \$input = \$ARGV[0];
my \$logfile = \$ARGV[1];
my \$psiautotest = \$ARGV[2];
my \$outfile = \"output.dat\";
if (@ARGV == 4) { \$outfile = $ARGV[3]; }
if (@ARGV == 5) { \$psi = $ARGV[4]; }

open(LOGFILE,\">\$logfile\") or die \"I can't write to \$logfile\\n\";

sub backtick(@)
{
    my \$pid = open(KID, '-|');
    die \"fork: \$!\" unless defined(\$pid);
    if (\$pid) {
        my \$output;
        while (<KID>) {
            print STDOUT;
            print LOGFILE;
        }
        close(KID);
        return \$output;
    } else {
        exec @_;
    }
}

my @cmd = (\$psi, \$input, \$outfile);
backtick(@cmd);

my \$return = \$? ? 1 : 0;
close LOGFILE;

if (\$psiautotest eq \"true\") {
  require(\"${TESTDIR}/psitest.pl\");
  do_tests(@ARGV);
}

exit \$return;\n")
#Finally run these beasts
set(LOGFILE ${CMAKE_BINARY_DIR}/testresults.log)
foreach(Tests2Run ${TestSet})
set(RUNDIR ${TESTDIR}/${Tests2Run})
set(INPUTFILE ${RUNDIR}/input.dat)
set(OUTFILE ${RUNDIR}/output.dat)
set(AUTOTEST 'true')
add_test(
    NAME ${Tests2Run}
    WORKING_DIRECTORY ${RUNDIR}
    COMMAND perl ${TESTEXE} ${INPUTFILE} ${LOGFILE} ${AUTOTEST} ${OUTFILE} ${PSIEXE})
endforeach(Tests2Run)
