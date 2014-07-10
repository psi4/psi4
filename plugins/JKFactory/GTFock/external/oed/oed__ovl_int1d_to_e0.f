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
         SUBROUTINE  OED__OVL_INT1D_TO_E0
     +
     +                    ( SHELLA,SHELLP,
     +                      ATOMIC,
     +                      NEXP,
     +                      NXYZET,NXYZP,
     +                      INT1DX,INT1DY,INT1DZ,
     +                      TEMP1,TEMP2,
     +                      SCALE,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__OVL_INT1D_TO_E0
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of cartesian
C                overlap integrals:
C
C                             [E|0] , E = A to P = A + B,
C
C                adding up all the contributions from all the 1D
C                integrals.
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
C                    ATOMIC      =  indicates, if purely atomic [E|0]
C                                   integrals will be evaluated
C                    NEXP        =  current # of exponent pairs
C                    NXYZET      =  sum of # of cartesian monomials
C                                   for all shells in the range
C                                   E = A,...,P=A+B
C                    NXYZP       =  # of cartesian monomials for the
C                                   P=A+B shell
C                    INT1Dx      =  all current 1D overlap integrals
C                                   for each cartesian component
C                                   x = X,Y,Z
C                    TEMP1(2)    =  scratch arrays holding intermediate
C                                   1D overlap integral products
C                    SCALE       =  the NEXP scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive cartesian [E|0]
C                                   overlap integrals corresponding
C                                   to all current exponent pairs
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

         LOGICAL     ATOMIC

         INTEGER     I,N
         INTEGER     NEXP
         INTEGER     NXYZE,NXYZET,NXYZP
         INTEGER     SE,SEEND
         INTEGER     SHELLA,SHELLP
         INTEGER     XE,YE,ZE
         INTEGER     XEMAX
         INTEGER     XEP,XYEP
         INTEGER     XYE
         INTEGER     YEEND

         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  SCALE (1:NEXP)
         DOUBLE PRECISION  TEMP1 (1:NEXP)
         DOUBLE PRECISION  TEMP2 (1:NEXP)

         DOUBLE PRECISION  BATCH (1:NEXP,1:NXYZET)

         DOUBLE PRECISION  INT1DX (1:NEXP,0:SHELLP)
         DOUBLE PRECISION  INT1DY (1:NEXP,0:SHELLP)
         DOUBLE PRECISION  INT1DZ (1:NEXP,0:SHELLP)

         PARAMETER  (ZERO  =  0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...the atomic case. Outer loop is over x contributions,
C                middle loop over y contributions and inner loop over
C                z contributions for different E shells. Skip
C                multiplications for x,y,z = 0, since the corresponding
C                1D integrals are equal to 1. Skip contributions if
C                x,y or z are odd, since the corresponding 1D integrals
C                are equal to 0.
C
C
         IF (ATOMIC) THEN

             DO I = 1,NXYZET
             DO N = 1,NEXP
                BATCH (N,I) = ZERO
             END DO
             END DO

             XEP = NXYZET + 3
             DO 100 XE = 0,SHELLP
                XEP = XEP + XE - 2
                XEMAX = XE * SHELLP
                YEEND = SHELLP - XE

                IF (MOD (XE,2).EQ.1) GOTO 100

                IF (XE.EQ.0) THEN
                    DO N = 1,NEXP
                       TEMP1 (N) = SCALE (N)
                    END DO
                ELSE
                    DO N = 1,NEXP
                       TEMP1 (N) = SCALE (N) * INT1DX (N,XE)
                    END DO
                END IF

                XYEP = XEP - XEMAX
                DO 200 YE = 0,YEEND
                   XYE = XE + YE
                   XYEP = XYEP - 1
                   SEEND = MAX0 (SHELLA,XYE)

                   IF (MOD (YE,2).EQ.1) GOTO 200

                   IF (YE.EQ.0) THEN
                       DO N = 1,NEXP
                          TEMP2 (N) = TEMP1 (N)
                       END DO
                   ELSE
                       DO N = 1,NEXP
                          TEMP2 (N) = TEMP1 (N) * INT1DY (N,YE)
                       END DO
                   END IF

                   I = XYEP
                   NXYZE = NXYZP
                   DO 300 SE = SHELLP,SEEND,-1
                      ZE = SE - XYE

                      IF (MOD (ZE,2).EQ.0) THEN
                          IF (ZE.EQ.0) THEN
                              DO N = 1,NEXP
                                 BATCH (N,I) = TEMP2 (N)
                              END DO
                          ELSE
                              DO N = 1,NEXP
                                 BATCH (N,I) = TEMP2 (N) * INT1DZ (N,ZE)
                              END DO
                          END IF
                      END IF

                      I = I - NXYZE + XE
                      NXYZE = NXYZE - SE - 1
  300              CONTINUE

  200           CONTINUE
  100        CONTINUE

         ELSE
C
C
C             ...the nonatomic case. Outer loop is over x contributions,
C                middle loop over y contributions and inner loop over
C                z contributions for different E shells. Skip
C                multiplications for x,y,z = 0, since the corresponding
C                1D integrals are equal to 1.
C
C
             XEP = NXYZET + 3
             DO 110 XE = 0,SHELLP
                XEP = XEP + XE - 2
                XEMAX = XE * SHELLP
                YEEND = SHELLP - XE

                IF (XE.EQ.0) THEN
                    DO N = 1,NEXP
                       TEMP1 (N) = SCALE (N)
                    END DO
                ELSE
                    DO N = 1,NEXP
                       TEMP1 (N) = SCALE (N) * INT1DX (N,XE)
                    END DO
                END IF

                XYEP = XEP - XEMAX
                DO 210 YE = 0,YEEND
                   XYE = XE + YE
                   XYEP = XYEP - 1
                   SEEND = MAX0 (SHELLA,XYE)

                   IF (YE.EQ.0) THEN
                       DO N = 1,NEXP
                          TEMP2 (N) = TEMP1 (N)
                       END DO
                   ELSE
                       DO N = 1,NEXP
                          TEMP2 (N) = TEMP1 (N) * INT1DY (N,YE)
                       END DO
                   END IF

                   I = XYEP
                   NXYZE = NXYZP
                   DO 310 SE = SHELLP,SEEND,-1
                      ZE = SE - XYE

                      IF (ZE.EQ.0) THEN
                          DO N = 1,NEXP
                             BATCH (N,I) = TEMP2 (N)
                          END DO
                      ELSE
                          DO N = 1,NEXP
                             BATCH (N,I) = TEMP2 (N) * INT1DZ (N,ZE)
                          END DO
                      END IF

                      I = I - NXYZE + XE
                      NXYZE = NXYZE - SE - 1
  310              CONTINUE

  210           CONTINUE
  110        CONTINUE

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
