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
         SUBROUTINE  OED__NAI_DERV_INT1D_TO_00
     +
     +                    ( NEXP,NGQPCEN,NGEXCEN,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      DIFFY,DIFFZ,
     +                      TEMP1,
     +                      SCALE,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_DERV_INT1D_TO_00
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of cartesian
C                derivative nuclear attraction integrals [0|0], adding
C                up the contributions from all the 1D 00 integrals for
C                all nuclear centers.
C
C                Simplified version of the general AB routine to reduce
C                loop overheads for those cases where only s-shells
C                are involved. For comments and details see the general
C                AB routine.
C
C
C                  Input:
C
C                    NEXP        =  current # of exponent pairs
C                    NGQPCEN     =  product of # of gaussian quadrature
C                                   points (roots) times # of nuclear
C                                   attraction centers
C                    NGEXCEN     =  product of # of gaussian quadrature
C                                   points times # of exponent pairs
C                                   times # of nuclear attraction
C                                   centers
C                    INT1Dx      =  all current 1D 00 derivative
C                                   nuclear attraction integrals for
C                                   each cartesian component (x = X,Y,Z)
C                    DIFFx       =  is true, if differentiation was
C                                   performed along the x=Y,Z direction
C                    TEMP1       =  scratch array holding intermediate
C                                   1D 00 derivative integral products
C                    SCALE       =  the NGEXCEN scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive cartesian
C                                   [0|0] derivative nuclear attraction
C                                   integrals corresponding to all
C                                   current exponent pairs
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     DIFFY,DIFFZ

         INTEGER     K,M,N
         INTEGER     NEXP,NGQPCEN,NGEXCEN

         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  SCALE (1:NGEXCEN)
         DOUBLE PRECISION  TEMP1 (1:NGEXCEN)

         DOUBLE PRECISION  BATCH (1:NEXP)

         DOUBLE PRECISION  INT1DX (1:NGEXCEN)
         DOUBLE PRECISION  INT1DY (1:NGEXCEN)
         DOUBLE PRECISION  INT1DZ (1:NGEXCEN)

         PARAMETER  (ZERO  =  0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...no skipping of 1DX integrals can be done here,
C                since the 1DX integrals carry the Rys weight!
C
C
         DO N = 1,NGEXCEN
            TEMP1 (N) = SCALE (N) * INT1DX (N)
         END DO
C
C
C             ...skip the multiplication of 1DY integrals, if no
C                y-coordinate derivative was formed, as then the
C                1DY integrals are equal to 1. Same with the 1DZ
C                1DZ integrals.
C
C
         IF (.NOT.DIFFY .AND. .NOT.DIFFZ) THEN
             K = 0
             DO M = 1,NEXP
                SUM = ZERO
                DO N = 1,NGQPCEN
                   SUM = SUM + TEMP1 (K+N)
                END DO
                K = K + NGQPCEN
                BATCH (M) = SUM
             END DO
         ELSE IF (.NOT.DIFFY) THEN
             K = 0
             DO M = 1,NEXP
                SUM = ZERO
                DO N = 1,NGQPCEN
                   SUM = SUM + TEMP1 (K+N) * INT1DZ (K+N)
                END DO
                K = K + NGQPCEN
                BATCH (M) = SUM
             END DO
         ELSE IF (.NOT.DIFFZ) THEN
             K = 0
             DO M = 1,NEXP
                SUM = ZERO
                DO N = 1,NGQPCEN
                   SUM = SUM + TEMP1 (K+N) * INT1DY (K+N)
                END DO
                K = K + NGQPCEN
                BATCH (M) = SUM
             END DO
         ELSE
             DO N = 1,NGEXCEN
                TEMP1 (N) = TEMP1 (N) * INT1DY (N)
             END DO
             K = 0
             DO M = 1,NEXP
                SUM = ZERO
                DO N = 1,NGQPCEN
                   SUM = SUM + TEMP1 (K+N) * INT1DZ (K+N)
                END DO
                K = K + NGQPCEN
                BATCH (M) = SUM
             END DO
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
