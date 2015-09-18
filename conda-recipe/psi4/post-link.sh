set +x off
echo ""
echo ""
echo "  Thank you for installing psi4. Additional resources:"
echo "    Website: www.psicode.org"
echo "    Inputs:  ${PREFIX}/share/psi/samples"
echo "    Manual:  http://psicode.org/psi4manual/master/index.html"
echo "    GitHub:  https://github.com/psi4/psi4public/wiki"
echo "    Binary:  https://anaconda.org/psi4"
echo "    Runtime Environment Diagnostic: ${PREFIX}/share/psi/scripts/setenv.py"
echo ""
echo "  For csh/tcsh command-line use, add to shell or ~/.tcshrc file:"
echo "    setenv PATH ${PREFIX}/bin:\$PATH"
echo "    setenv PSI_SCRATCH /path/to/existing/writable/local-not-network/disk/for/scratch/files"
echo ""
echo "  For sh/bash command-line use, add to shell or ~/.bashrc file:"
echo "    export PATH=${PREFIX}/bin:\$PATH"
echo "    export PSI_SCRATCH=/path/to/existing/writable/local-not-network/disk/for/scratch/files"
echo ""

cat > linktest.in << EOL
memory 250 mb

molecule dimer {
0 1
C   0.000000  -0.667578  -2.124659
C   0.000000   0.667578  -2.124659
H   0.923621  -1.232253  -2.126185
H  -0.923621  -1.232253  -2.126185
H  -0.923621   1.232253  -2.126185
H   0.923621   1.232253  -2.126185
--
0 1
C   0.000000   0.000000   2.900503
C   0.000000   0.000000   1.693240
H   0.000000   0.000000   0.627352
H   0.000000   0.000000   3.963929
}

set {
  BASIS jun-cc-pVDZ
  SCF_TYPE DF
  FREEZE_CORE True
}

energy('sapt0')

compare_values(85.189064196429, dimer.nuclear_repulsion_energy(), 9, "Nuclear Repulsion Energy")
compare_values(-0.003378388762, psi4.get_variable("SAPT ELST ENERGY"), 6, "SAPT0 Eelst")
compare_values( 0.003704416103, psi4.get_variable("SAPT EXCH ENERGY"), 6, "SAPT0 Eexch")
compare_values(-0.000889316601, psi4.get_variable("SAPT IND ENERGY"), 6, "SAPT0 Eind")
compare_values(-0.001672292164, psi4.get_variable("SAPT DISP ENERGY"), 6, "SAPT0 Edisp")
compare_values(-0.002235581423, psi4.get_variable("SAPT SAPT0 ENERGY"), 6, "SAPT0 Etotal")
EOL

PSIOUT=`PSI_SCRATCH=/tmp ${PREFIX}/bin/psi4 linktest.in`
echo $PSIOUT
echo ""
