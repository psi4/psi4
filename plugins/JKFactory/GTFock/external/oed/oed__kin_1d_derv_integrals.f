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
         SUBROUTINE  OED__KIN_1D_DERV_INTEGRALS
     +
     +                    ( SHELLA,SHELLB,
     +                      NEXP,
     +                      EA,EB,E2AB,
     +                      OVL1DX,OVL1DY,OVL1DZ,
     +
     +                               KIN1DX,
     +                               KIN1DY,
     +                               KIN1DZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__KIN_1D_DERV_INTEGRALS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation calculates a full table of kinetic 1D
C                X,Y,Z derivative integrals using the previously
C                established overlap 1D X,Y,Z derivative integrals.
C
C                The I,J-th 1D kinetic derivative integral is given
C                by the following formula:
C
C                    KIN1D (I,J)  =    1/2 * I * J * OVL1D (I-1,J-1)
C                                   -   I * ALPHAB * OVL1D (I-1,J+1)
C                                   -   J * ALPHAA * OVL1D (I+1,J-1)
C                           +  2 * ALPHAA * ALPHAB * OVL1D (I+1,J+1)
C
C
C                The 1D kinetic derivative integrals are calculated
C                for all exponent pairs simultaneously and placed into
C                a 3-dimensional array with the exponent pair index
C                varying fastest.
C
C
C                  Input:
C
C                    SHELLx       =  the shell type for contraction
C                                    shells x = A and B
C                    NEXP         =  current # of exponent pairs
C                    Ex           =  current MIJ exponents from centers
C                                    x = A,B
C                    E2AB         =  current MIJ double exponent
C                                    products between centers A and B
C                    OVL1Dx       =  current 1D overlap derivative
C                                    integrals for each cartesian
C                                    component (x=X,Y,Z)
C
C                  Output:
C
C                    KIN1Dx       =  current 1D kinetic derivative
C                                    integrals for each cartesian
C                                    component (x=X,Y,Z)
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT  NONE

         INTEGER   CASE1D
         INTEGER   I,J,N
         INTEGER   IM1,IP1,JM1,JP1
         INTEGER   NEXP
         INTEGER   SHELLA,SHELLB

         DOUBLE PRECISION  HALF,HALFJ,HALFIJ
         DOUBLE PRECISION  XI,XJ,XIB,XJA,X2AB

         DOUBLE PRECISION  EA     (1:NEXP)
         DOUBLE PRECISION  EB     (1:NEXP)
         DOUBLE PRECISION  E2AB   (1:NEXP)

         DOUBLE PRECISION  KIN1DX (1:NEXP,0:SHELLA  ,0:SHELLB  )
         DOUBLE PRECISION  KIN1DY (1:NEXP,0:SHELLA  ,0:SHELLB  )
         DOUBLE PRECISION  KIN1DZ (1:NEXP,0:SHELLA  ,0:SHELLB  )
         DOUBLE PRECISION  OVL1DX (1:NEXP,0:SHELLA+1,0:SHELLB+1)
         DOUBLE PRECISION  OVL1DY (1:NEXP,0:SHELLA+1,0:SHELLB+1)
         DOUBLE PRECISION  OVL1DZ (1:NEXP,0:SHELLA+1,0:SHELLB+1)

         PARAMETER  (HALF = 0.5D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...the case A = s-shell and B = s-shell.
C
C
         DO N = 1,NEXP
            X2AB = E2AB (N)
            KIN1DX (N,0,0) = X2AB * OVL1DX (N,1,1)
            KIN1DY (N,0,0) = X2AB * OVL1DY (N,1,1)
            KIN1DZ (N,0,0) = X2AB * OVL1DZ (N,1,1)
         END DO
C
C
C             ...jump according to the 4 different cases that can arise:
C
C                  A-shell = s- or higher angular momentum
C                  B-shell = s- or higher angular momentum
C
C                each leading to specific simplifications.
C
C
         CASE1D = 2 * MIN (1,SHELLA) + MIN (1,SHELLB) + 1

         GOTO (1,2,3,4) CASE1D
C
C
C             ...the case A = s-shell and B = s-shell (done already,
C                immediate return).
C
C
    1    RETURN
C
C
C             ...the cases A = s-shell and B >= p-shell.
C
C
    2    DO J = 1,SHELLB
            JM1 = J - 1
            JP1 = J + 1
            XJ = DFLOAT (J)

            DO N = 1,NEXP
               XJA = XJ * EA (N)
               X2AB = E2AB (N)
               KIN1DX (N,0,J)  =  -  XJA * OVL1DX (N,1,JM1)
     +                            + X2AB * OVL1DX (N,1,JP1)
               KIN1DY (N,0,J)  =  -  XJA * OVL1DY (N,1,JM1)
     +                            + X2AB * OVL1DY (N,1,JP1)
               KIN1DZ (N,0,J)  =  -  XJA * OVL1DZ (N,1,JM1)
     +                            + X2AB * OVL1DZ (N,1,JP1)
            END DO
         END DO

         RETURN
C
C
C             ...the cases A >= p-shell and B = s-shell.
C
C
    3    DO I = 1,SHELLA
            IM1 = I - 1
            IP1 = I + 1
            XI = DFLOAT (I)
            DO N = 1,NEXP
               XIB = XI * EB (N)
               X2AB = E2AB (N)
               KIN1DX (N,I,0)  =  -  XIB * OVL1DX (N,IM1,1)
     +                            + X2AB * OVL1DX (N,IP1,1)
               KIN1DY (N,I,0)  =  -  XIB * OVL1DY (N,IM1,1)
     +                            + X2AB * OVL1DY (N,IP1,1)
               KIN1DZ (N,I,0)  =  -  XIB * OVL1DZ (N,IM1,1)
     +                            + X2AB * OVL1DZ (N,IP1,1)
            END DO
         END DO

         RETURN
C
C
C             ...the cases A >= p-shell and B >= p-shell.
C
C
    4    DO I = 1,SHELLA
            IM1 = I - 1
            IP1 = I + 1
            XI = DFLOAT (I)
            DO N = 1,NEXP
               XIB = XI * EB (N)
               X2AB = E2AB (N)
               KIN1DX (N,I,0)  =  -  XIB * OVL1DX (N,IM1,1)
     +                            + X2AB * OVL1DX (N,IP1,1)
               KIN1DY (N,I,0)  =  -  XIB * OVL1DY (N,IM1,1)
     +                            + X2AB * OVL1DY (N,IP1,1)
               KIN1DZ (N,I,0)  =  -  XIB * OVL1DZ (N,IM1,1)
     +                            + X2AB * OVL1DZ (N,IP1,1)
            END DO
         END DO

         DO J = 1,SHELLB
            JM1 = J - 1
            JP1 = J + 1
            XJ = DFLOAT (J)

            DO N = 1,NEXP
               XJA = XJ * EA (N)
               X2AB = E2AB (N)
               KIN1DX (N,0,J)  =  -  XJA * OVL1DX (N,1,JM1)
     +                            + X2AB * OVL1DX (N,1,JP1)
               KIN1DY (N,0,J)  =  -  XJA * OVL1DY (N,1,JM1)
     +                            + X2AB * OVL1DY (N,1,JP1)
               KIN1DZ (N,0,J)  =  -  XJA * OVL1DZ (N,1,JM1)
     +                            + X2AB * OVL1DZ (N,1,JP1)
            END DO

            HALFJ = HALF * XJ
            DO I = 1,SHELLA
               IM1 = I - 1
               IP1 = I + 1
               XI = DFLOAT (I)
               HALFIJ = HALFJ * XI
               DO N = 1,NEXP
                  XIB = XI * EB (N)
                  XJA = XJ * EA (N)
                  X2AB = E2AB (N)

                  KIN1DX (N,I,J)  =  HALFIJ * OVL1DX (N,IM1,JM1)
     +                               -  XIB * OVL1DX (N,IM1,JP1)
     +                               -  XJA * OVL1DX (N,IP1,JM1)
     +                               + X2AB * OVL1DX (N,IP1,JP1)
                  KIN1DY (N,I,J)  =  HALFIJ * OVL1DY (N,IM1,JM1)
     +                               -  XIB * OVL1DY (N,IM1,JP1)
     +                               -  XJA * OVL1DY (N,IP1,JM1)
     +                               + X2AB * OVL1DY (N,IP1,JP1)
                  KIN1DZ (N,I,J)  =  HALFIJ * OVL1DZ (N,IM1,JM1)
     +                               -  XIB * OVL1DZ (N,IM1,JP1)
     +                               -  XJA * OVL1DZ (N,IP1,JM1)
     +                               + X2AB * OVL1DZ (N,IP1,JP1)
               END DO
            END DO
         END DO
C
C
C             ...ready!
C
C
         RETURN
         END
