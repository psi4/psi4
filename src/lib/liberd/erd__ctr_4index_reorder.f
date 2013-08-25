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
         SUBROUTINE  ERD__CTR_4INDEX_REORDER
     +
     +                    ( NXYZT,NCTR,
     +                      NCGTOR,NCGTOS,NCGTOT,NCGTOU,
     +                      IXOFFI,IXOFFJ,IXOFFK,IXOFFL,
     +                      X,
     +
     +                             Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__CTR_4INDEX_REORDER
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine reorders the contraction indices r,s,t,u:
C
C                        x (nxyzt,rstu) --> y (nxyzt,ijkl)
C
C                such that the contraction indices correspond to the
C                final ordering of the integral batch.
C
C                The trick used to do a very general permutational
C                reordering without the need to write specific routines
C                for all possible 4!=24 cases, is to use predefined
C                index offsets, which contain implicitly the permutation
C                information. These index offsets allow one to
C                calculate the target address for each rstu index
C                quadruple using a simple multiplication formula,
C                which of course is here transformed into a repeated
C                summation form to save cpu time.
C
C
C                  Input:
C
C                      NXYZT      =  total # of cartesian monomial
C                                    quadruplets (invariant indices)
C                      NCTR       =  total # of contraction indices
C                      NCGTOx     =  # of contractions for contraction
C                                    shells x=R,S,T,U
C                      IXOFFx     =  index offsets for contraction
C                                    shells x=I,J,K,L
C                      X          =  original integral batch with
C                                    contraction shells in RSTU order
C
C                  Output:
C
C                      Y          =  new integral batch with contraction
C                                    shells in IJKL order
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

         INTEGER     IJKL,JKL,KL,L,N
         INTEGER     IXOFFI,IXOFFJ,IXOFFK,IXOFFL
         INTEGER     NCGTOR,NCGTOS,NCGTOT,NCGTOU
         INTEGER     NCTR
         INTEGER     NXYZT
         INTEGER     R,S,T,U
         INTEGER     RSTU

         DOUBLE PRECISION  X (1:NXYZT,1:NCTR)
         DOUBLE PRECISION  Y (1:NXYZT,1:NCTR)
C
C
C------------------------------------------------------------------------
C
C
C             ...do the reordering. Target index IJKL is evaluated
C                from each R,S,T,U index quadruple using summations
C                only!
C
C
         RSTU = 0

         L = - IXOFFL + 1
         DO 100 U = 1,NCGTOU
            L = L + IXOFFL
            KL = L - IXOFFK
            DO 200 T = 1,NCGTOT
               KL = KL + IXOFFK
               JKL = KL - IXOFFJ
               DO 300 S = 1,NCGTOS
                  JKL = JKL + IXOFFJ
                  IJKL = JKL - IXOFFI
                  DO 400 R = 1,NCGTOR
                     IJKL = IJKL + IXOFFI
                     RSTU = RSTU + 1

                     DO N = 1,NXYZT
                        Y (N,IJKL) = X (N,RSTU)
                     END DO

  400             CONTINUE
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
