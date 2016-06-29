C  Copyright (c) 2003-2010 University of Florida
C
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  The GNU General Public License is included in this distribution
C  in the file COPYRIGHT.
         SUBROUTINE  ERD__CTR_TU_EXPAND
     +
     +                    ( NXYZRS,NTU,
     +                      NT,NU,
     +                      X,
     +
     +                             Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__CTR_TU_EXPAND
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine expands the tu contraction indices:
C
C                      x (nxyzt,rs,t>=u) --> y (nxyzt,rs,tu)
C
C
C                  Input:
C
C                      NXYZRS     =  total # of cartesian monomial
C                                    quadruplets times # of contractions
C                                    for contraction shell pair rs
C                                    (invariant indices)
C                      NTU        =  total # of contraction index pairs
C                                    for contraction shells tu before
C                                    expansion
C                      NT(U)      =  # of contractions for contraction
C                                    shells t(u)
C                      X          =  original integral batch with
C                                    contraction index restriction
C                                    t>=u
C
C                  Output:
C
C                      Y          =  new integral batch with complete
C                                    set of contraction indices tu
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         INTEGER     N,T,U
         INTEGER     NT,NU
         INTEGER     NTU
         INTEGER     NXYZRS
         INTEGER     TU

         DOUBLE PRECISION  X (1:NXYZRS,1:NTU)
         DOUBLE PRECISION  Y (1:NXYZRS,1:NT,1:NU)
C
C
C------------------------------------------------------------------------
C
C
C             ...do the expansion.
C
C
         TU = NTU + 1

         DO 100 U = NU,1,-1

            DO 200 T = NT,U+1,-1

               TU = TU - 1

               DO N = 1,NXYZRS
                  Y (N,T,U) = X (N,TU)
               END DO

               DO N = 1,NXYZRS
                  Y (N,U,T) = X (N,TU)
               END DO

  200       CONTINUE

            TU = TU - 1
            DO N = 1,NXYZRS
               Y (N,U,U) = X (N,TU)
            END DO

  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
