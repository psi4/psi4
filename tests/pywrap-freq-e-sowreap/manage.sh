# modify objdir here. assumes running from test dir
OBJDIR="objdir-pn/stage/usr/local/psi4"
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
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i FREQ-8.in                   -o FREQ-8.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i FREQ-9.in                   -o FREQ-9.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i FREQ-10.in                  -o FREQ-10.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i FREQ-11.in                  -o FREQ-11.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i FREQ-12.in                  -o FREQ-12.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i FREQ-13.in                  -o FREQ-13.out

cat tests >> FREQ-master.in
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i FREQ-master.in              -o FREQ-master.out

