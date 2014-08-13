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
         SUBROUTINE  OED__NAI_INT1D_TO_E0
     +
     +                    ( SHELLA,SHELLP,
     +                      NEXP,NGQPCEN,NGEXCEN,
     +                      NXYZET,NXYZP,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      TEMP1,TEMP2,
     +                      SCALE,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_INT1D_TO_E0
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of cartesian
C                nuclear attraction integrals:
C
C                             [E|0] , E = A to P = A + B,
C
C                adding up the contributions from all the 1D integrals
C                for all nuclear centers.
C
C                The routine uses the reduced Rys multiplication scheme
C                as suggested in R.Lindh, U.Ryu and B.Liu, JCP 95, 5889.
C                This scheme reuses intermediate products between
C                1DX and 1DY integrals, which can be achieved by having
C                the outer loops run over all possible x and y monomial
C                parts and the inner loop over all allowed E shells.
C
C                For details of the algorithm used, especially how the
C                monomial basis is organized and how to find the batch
C                address for a specific x,y,z,E combination, see the
C                corresponding routine for the ERI evaluations.
C
C
C                  Input:
C
C                    SHELLx      =  shell types for individual csh
C                                   x=A and csh sum P=A+B
C                    NEXP        =  current # of exponent pairs
C                    NGQPCEN     =  product of # of gaussian quadrature
C                                   points (roots) times # of nuclear
C                                   attraction centers
C                    NGEXCEN     =  product of # of gaussian quadrature
C                                   points times # of exponent pairs
C                                   times # of nuclear attraction
C                                   centers
C                    NXYZET      =  sum of # of cartesian monomials
C                                   for all shells in the range
C                                   E = A,...,P=A+B
C                    NXYZP       =  # of cartesian monomials for the
C                                   P=A+B shell
C                    INT1Dx      =  all current 1D nuclear attraction
C                                   integrals for each cartesian
C                                   component (x = X,Y,Z)
C                    TEMP1(2)    =  scratch arrays holding intermediate
C                                   1D integral products
C                    SCALE       =  the NGEXCEN scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive cartesian
C                                   nuclear attraction integrals
C                                   corresponding to all current
C                                   exponent pairs
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         INTEGER     I,K,M,N
         INTEGER     NEXP,NGQPCEN,NGEXCEN
         INTEGER     NXYZE,NXYZET,NXYZP
         INTEGER     SE,SEEND
         INTEGER     SHELLA,SHELLP
         INTEGER     XE,YE,ZE
         INTEGER     XEMAX
         INTEGER     XEP,XYEP
         INTEGER     XYE
         INTEGER     YEEND

         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  SCALE (1:NGEXCEN)
         DOUBLE PRECISION  TEMP1 (1:NGEXCEN)
         DOUBLE PRECISION  TEMP2 (1:NGEXCEN)

         DOUBLE PRECISION  BATCH (1:NEXP,1:NXYZET)

         DOUBLE PRECISION  INT1DX (1:NGEXCEN,0:SHELLP)
         DOUBLE PRECISION  INT1DY (1:NGEXCEN,0:SHELLP)
         DOUBLE PRECISION  INT1DZ (1:NGEXCEN,0:SHELLP)

         PARAMETER  (ZERO  =  0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...outer loop over x contributions. No skipping of the
C                x contribution of 0-type can be done here, since
C                the 1DX integrals carry the Rys weight!
C
C
         XEP = NXYZET + 3
         DO 100 XE = 0,SHELLP
            XEP = XEP + XE - 2
            XEMAX = XE * SHELLP
            YEEND = SHELLP - XE

            DO N = 1,NGEXCEN
               TEMP1 (N) = SCALE (N) * INT1DX (N,XE)
            END DO
C
C
C             ...middle loop over y contributions. Skip multiplication
C                of y contributions, if we have y = 0, as then the
C                1DY integral is equal to 1.
C
C
            XYEP = XEP - XEMAX
            DO 200 YE = 0,YEEND
               XYE = XE + YE
               XYEP = XYEP - 1
               SEEND = MAX (SHELLA,XYE)

               IF (YE.EQ.0) THEN
                   DO N = 1,NGEXCEN
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGEXCEN
                      TEMP2 (N) = TEMP1 (N) * INT1DY (N,YE)
                   END DO
               END IF
C
C
C             ...inner loop over z,E contributions. Skip multiplication
C                of z,E contributions, if we have z = 0, as then the
C                1DY integral is equal to 1.
C
C
               I = XYEP
               NXYZE = NXYZP
               DO 300 SE = SHELLP,SEEND,-1
                  ZE = SE - XYE
C
C
C             ...all info concerning all three x,y and z contributions
C                have been collected for all exponent pairs, nuclear
C                centers and quadrature points at once. Sum up the
C                1D X,Y,Z integral products over the nuclear centers
C                and the quadrature points to the appropriate place of
C                the [E|0] batch.
C
C
                  IF (ZE.EQ.0) THEN
                      K = 0
                      DO M = 1,NEXP
                         SUM = ZERO
                         DO N = 1,NGQPCEN
                            SUM = SUM + TEMP2 (K+N)
                         END DO
                         K = K + NGQPCEN
                         BATCH (M,I) = SUM
                      END DO
                  ELSE
                      K = 0
                      DO M = 1,NEXP
                         SUM = ZERO
                         DO N = 1,NGQPCEN
                            SUM = SUM + TEMP2 (K+N) * INT1DZ (K+N,ZE)
                         END DO
                         K = K + NGQPCEN
                         BATCH (M,I) = SUM
                      END DO
                  END IF
C
C
C             ...next z,E contribution.
C
C
                  I = I - NXYZE + XE
                  NXYZE = NXYZE - SE - 1
  300          CONTINUE
C
C
C             ...next y and x contribution.
C
C
  200       CONTINUE
  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
