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
         SUBROUTINE  OED__CTR_2INDEX_REORDER
     +
     +                    ( NXYZT,NCTR,
     +                      NCGTOR,NCGTOS,
     +                      IXOFFI,IXOFFJ,
     +                      X,
     +
     +                             Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__CTR_2INDEX_REORDER
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine reorders the 2 contraction indices r,s:
C
C                        x (nxyzt,rs) --> y (nxyzt,ij)
C
C                such that the contraction indices correspond to the
C                final ordering of the batch.
C
C
C                  Input:
C
C                      NXYZT      =  total # of cartesian monomial
C                                    quadruplets (invariant indices)
C                      NCTR       =  total # of contraction indices
C                      NCGTOx     =  # of contractions for contraction
C                                    shells x=R,S
C                      IXOFFx     =  index offsets for contraction
C                                    shells x=I,J
C                      X          =  original integral batch with
C                                    contraction shells in RS order
C
C                  Output:
C
C                      Y          =  new integral batch with contraction
C                                    shells in IJ order
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         INTEGER     J,IJ
         INTEGER     IXOFFI,IXOFFJ
         INTEGER     N,R,S
         INTEGER     NCGTOR,NCGTOS
         INTEGER     NCTR
         INTEGER     NXYZT
         INTEGER     RS

         DOUBLE PRECISION  X (1:NXYZT,1:NCTR)
         DOUBLE PRECISION  Y (1:NXYZT,1:NCTR)
C
C
C------------------------------------------------------------------------
C
C
C             ...do the reordering.
C
C
         RS = 0

         J = - IXOFFJ + 1
         DO 100 S = 1,NCGTOS
            J = J + IXOFFJ
            IJ = J - IXOFFI
            DO 200 R = 1,NCGTOR
               IJ = IJ + IXOFFI
               RS = RS + 1

               DO N = 1,NXYZT
                  Y (N,IJ) = X (N,RS)
               END DO

  200       CONTINUE
  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
