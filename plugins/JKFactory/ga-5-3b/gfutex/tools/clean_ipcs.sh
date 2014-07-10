#!/bin/bash

for d in `ipcs -q | awk '{ print $2 }' | grep "[0-9]"`
do
  ipcrm -q $d
done

for d in `ipcs -m | awk '{ print $2 }' | grep "[0-9]"`
do
  ipcrm -m $d
done
