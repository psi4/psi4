# modify objdir here. assumes running from test dir
OBJDIR="objdir-pn"

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
../../${OBJDIR}/bin/psi4 -l ../../share/

../../${OBJDIR}/bin/psi4 -l ../../share/ -i BASIC-ch4-reagent.in        -o BASIC-ch4-reagent.out
../../${OBJDIR}/bin/psi4 -l ../../share/ -i BASIC-nh3-reagent.in        -o BASIC-nh3-reagent.out

cat tests >> BASIC-master.in
../../${OBJDIR}/bin/psi4 -l ../../share/ -i BASIC-master.in             -o BASIC-master.out

