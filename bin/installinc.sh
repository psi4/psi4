#!/bin/sh -f
#
# 
#

DEST=0
updates=0

for ii in $*
do
	if [ -d $ii ] ; then DEST=$ii; fi
done

if [ $DEST = 0 ]
then
	cmp $1 $2 > /dev/null 2>&1
	if [ $? != 0 ]
	then
		if [ -f $2 ] ; then chmod 644 $2; fi
		/bin/cp $1 $2
		if [ $updates = 0 ] ; then echo -n " updating"; updates=1; fi
		echo -n " "$1
	fi
else
	
for ii in $*
do
	if [ ! -d $ii ] 
	then

		cmp $ii $DEST/`basename $ii` > /dev/null 2>&1
		if [ $? != 0 ]
		then
			if [ -f $DEST/`basename $ii` ] ; then chmod 644 `basename $DEST/$ii`; fi
			/bin/cp $ii $DEST/`basename $ii`
			if [ $updates = 0 ] ; then echo -n " updating"; updates=1; fi
			echo -n " "$ii
		fi
	fi
done
fi

if [ $updates = 0 ] ; then echo "no updates"; else echo ""; fi

exit
