#!/bin/csh

set n = 0

cd displacements

while ($n < 19)
  psi4 -i $n-disp.in -o $n-disp.out
  @ n++
end

