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
         SUBROUTINE  OED__NAI_1D_P_INTEGRALS
     +
     +                    ( SHELLP,
     +                      NGEXCEN,
     +                      WTS,
     +                      R1X,R1Y,R1Z,
     +                      R2,
     +
     +                               INT1DX,
     +                               INT1DY,
     +                               INT1DZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_1D_P_INTEGRALS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation calculates a full table of 1D P X,Y,Z
C                nuclear attraction integrals using the Rys vertical
C                recurrence scheme explained below.
C
C                The recurrence scheme is due to Rys, Dupuis and King,
C                J. Comp. Chem. 4, p.154-157 (1983).
C
C
C                   INT1D (0,0) = 1.D0    (* WEIGHT for the 1DX case)
C                   INT1D (1,0) = R1      (* WEIGHT for the 1DX case)
C
C                   For I = 1,...,SHELLP-1
C                       INT1D (I+1,0) = I * R2 * INT1D (I-1,0)
C                                         + R1 * INT1D (I,0)
C
C
C                The 1D P X,Y,Z integrals are calculated for all nuclear
C                centers, all exponent pairs and all quadrature points
C                simultaneously and placed into a 2-dimensional array.
C
C
C                  Input:
C
C                    SHELLP      =  maximum shell type A+B
C                    NGEXCEN     =  # of roots times # of primitive
C                                   exponent pairs times # of nuclear
C                                   attraction centers
C                    WTS         =  all quadrature weights
C                    R1x         =  the VRR R1-coefficients (individual
C                                   cartesian components x=X,Y,Z) for
C                                   shell expansion on center P
C                    R2          =  the coordinate independent VRR
C                                   R2-coefficients
C
C
C                  Output:
C
C                    INT1Dx      =  all 1D P nuclear attraction
C                                   integrals for each cartesian
C                                   component (x = X,Y,Z)
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT  NONE

         INTEGER   I,N
         INTEGER   I1,I2
         INTEGER   NGEXCEN
         INTEGER   SHELLP

         DOUBLE PRECISION  F
         DOUBLE PRECISION  ONE
         DOUBLE PRECISION  R2N

         DOUBLE PRECISION  R1X  (1:NGEXCEN)
         DOUBLE PRECISION  R1Y  (1:NGEXCEN)
         DOUBLE PRECISION  R1Z  (1:NGEXCEN)
         DOUBLE PRECISION  R2   (1:NGEXCEN)
         DOUBLE PRECISION  WTS  (1:NGEXCEN)

         DOUBLE PRECISION  INT1DX (1:NGEXCEN,0:SHELLP)
         DOUBLE PRECISION  INT1DY (1:NGEXCEN,0:SHELLP)
         DOUBLE PRECISION  INT1DZ (1:NGEXCEN,0:SHELLP)

         PARAMETER  (ONE   =  1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...initialize 1D nuclear attraction integrals.
C
C
         DO 100 N = 1,NGEXCEN
            INT1DX (N,0) = WTS (N)
            INT1DY (N,0) = ONE
            INT1DZ (N,0) = ONE
  100    CONTINUE

         IF (SHELLP.EQ.0) RETURN
C
C
C             ...proceed, if total shell is > s-shell.
C
C
         DO 200 N = 1,NGEXCEN
            INT1DX (N,1) = WTS (N) * R1X (N)
            INT1DY (N,1) = R1Y (N)
            INT1DZ (N,1) = R1Z (N)
  200    CONTINUE

         F = ONE
         DO 220 I = 2,SHELLP
            I1 = I - 1
            I2 = I - 2
            DO 230 N = 1,NGEXCEN
               R2N = F * R2 (N)
               INT1DX (N,I) =       R2N * INT1DX (N,I2)
     +                        + R1X (N) * INT1DX (N,I1)
               INT1DY (N,I) =       R2N * INT1DY (N,I2)
     +                        + R1Y (N) * INT1DY (N,I1)
               INT1DZ (N,I) =       R2N * INT1DZ (N,I2)
     +                        + R1Z (N) * INT1DZ (N,I1)
  230       CONTINUE
            F = F + ONE
  220    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
