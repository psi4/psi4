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
rm -f OPT.h2s.grad

# echo to screen
set -x

# run calculation sequence
../../${OBJDIR}/bin/psi4 ${DVRDIR}

../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-1-1.in                  -o OPT-1-1.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-1-2.in                  -o OPT-1-2.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-1-3.in                  -o OPT-1-3.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-1-4.in                  -o OPT-1-4.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-1-5.in                  -o OPT-1-5.out

../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-master.in               -o OPT-master.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-2-1.in                  -o OPT-2-1.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-2-2.in                  -o OPT-2-2.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-2-3.in                  -o OPT-2-3.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-2-4.in                  -o OPT-2-4.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-2-5.in                  -o OPT-2-5.out

../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-master.in               -o OPT-master.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-3-1.in                  -o OPT-3-1.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-3-2.in                  -o OPT-3-2.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-3-3.in                  -o OPT-3-3.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-3-4.in                  -o OPT-3-4.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-3-5.in                  -o OPT-3-5.out

../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-master.in               -o OPT-master.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-4-1.in                  -o OPT-4-1.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-4-2.in                  -o OPT-4-2.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-4-3.in                  -o OPT-4-3.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-4-4.in                  -o OPT-4-4.out
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-4-5.in                  -o OPT-4-5.out

cat tests >> OPT-master.in
../../${OBJDIR}/bin/psi4 ${DVRDIR} -i OPT-master.in               -o OPT-master.out

