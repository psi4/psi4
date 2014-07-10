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
         SUBROUTINE  OED__CTR_3INDEX_REORDER
     +
     +                    ( NXYZT,NCTR,
     +                      NCGTOR,NCGTOS,NCGTOT,
     +                      IXOFFI,IXOFFJ,IXOFFK,
     +                      X,
     +
     +                             Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__CTR_3INDEX_REORDER
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine reorders the 3 contraction indices r,s,t:
C
C                        x (nxyzt,rst) --> y (nxyzt,ijk)
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
C                                    shells x=R,S,T
C                      IXOFFx     =  index offsets for contraction
C                                    shells x=I,J,K
C                      X          =  original integral batch with
C                                    contraction shells in RST order
C
C                  Output:
C
C                      Y          =  new integral batch with contraction
C                                    shells in IJK order
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

         INTEGER     IJK,JK,K
         INTEGER     IXOFFI,IXOFFJ,IXOFFK
         INTEGER     N,R,S,T
         INTEGER     NCGTOR,NCGTOS,NCGTOT
         INTEGER     NCTR
         INTEGER     NXYZT
         INTEGER     RST

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
         RST = 0

         K = - IXOFFK + 1
         DO 100 T = 1,NCGTOT
            K = K + IXOFFK
            JK = K - IXOFFJ
            DO 200 S = 1,NCGTOS
               JK = JK + IXOFFJ
               IJK = JK - IXOFFI
               DO 300 R = 1,NCGTOR
                  IJK = IJK + IXOFFI
                  RST = RST + 1

                  DO N = 1,NXYZT
                     Y (N,IJK) = X (N,RST)
                  END DO

  300          CONTINUE
  200       CONTINUE
  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
