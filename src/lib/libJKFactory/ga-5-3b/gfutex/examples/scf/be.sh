#!/bin/sh
## Generate input file for SCF

if [ $# -lt 1 ]
then
    echo "Usage: $0 <number of Be atoms>"
    exit 1
fi

file=be.inpt
/bin/rm -rf $file

echo $1 >> $file
echo >> $file
iters=`expr $1 - 1`

for i in `seq 0 $iters`
do
    echo "4  `expr 4 "*" $i`.000 0.000 0.000" >> $file
done
