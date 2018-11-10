#!/bin/csh

cd displacements

foreach n (*.in)
   psi4 -i "$n"
end

