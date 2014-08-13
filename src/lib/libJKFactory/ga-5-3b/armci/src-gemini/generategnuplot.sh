#!/bin/sh
echo "#Range number process" > profile_3d.dat
for ((  i = 0 ;  i <= 128;  i++  ))
do
  affile="armci_profile.${i}"
  if test -s $affile
  then 
	  head -n 28 $affile | tail -n 22 | awk '{print $7" "$0+$1+$2" '$i'"}' | awk -F- '{print $2}' | awk -F")" '{print $1" "$2}' >> profile_3d.dat

  fi
done
