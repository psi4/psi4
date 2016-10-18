# modify objdir here. assumes running from test dir
OBJDIR="objdir-pn"
DVRDIR="" #"-l ../../share/"

# clean
rm -f psi*clean
rm -f *intco
rm -f output.dat
rm -f timer.dat
rm -f FREQ-*in
rm -f FREQ-*out
rm -f psi.*.1
rm -f OPT-*.grad
rm -f OPT-*.in
rm -f OPT-*.out
rm -f BASIC-*in
rm -f BASIC-*out

# echo to screen
set -x

# run calculation sequence
../../${OBJDIR}/bin/psi4 ${DVRDIR}

../../${OBJDIR}/bin/psi4 ${DVRDIR} -i FREQ-1.in                   -o FREQ-1.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i FREQ-2.in                   -o FREQ-2.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i FREQ-3.in                   -o FREQ-3.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i FREQ-4.in                   -o FREQ-4.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i FREQ-5.in                   -o FREQ-5.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i FREQ-6.in                   -o FREQ-6.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i FREQ-7.in                   -o FREQ-7.out

cat tests >> FREQ-master.in
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i FREQ-master.in              -o FREQ-master.out

