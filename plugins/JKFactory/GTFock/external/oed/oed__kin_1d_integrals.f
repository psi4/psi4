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
         SUBROUTINE  OED__KIN_1D_INTEGRALS
     +
     +                    ( SHELLA,SHELLB,SHELLP,
     +                      ATOMIC,
     +                      NEXP,
     +                      EA,EB,E2AB,
     +                      OVL1DX,OVL1DY,OVL1DZ,
     +
     +                               KIN1DX,
     +                               KIN1DY,
     +                               KIN1DZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__KIN_1D_INTEGRALS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation calculates a full table of kinetic 1D
C                X,Y,Z integrals using the previously established
C                overlap 1D X,Y,Z overlap integrals.
C
C                The I,J-th 1D kinetic integral is given by the
C                following formula:
C
C                    KIN1D (I,J)  =    1/2 * I * J * OVL1D (I-1,J-1)
C                                   -   I * ALPHAB * OVL1D (I-1,J+1)
C                                   -   J * ALPHAA * OVL1D (I+1,J-1)
C                           +  2 * ALPHAA * ALPHAB * OVL1D (I+1,J+1)
C
C
C                The 1D kinetic integrals are calculated for all
C                exponent pairs simultaneously and placed into a
C                3-dimensional array with the exponent pair index
C                varying fastest.
C
C
C                  Input:
C
C                    SHELLx       =  the shell type for contraction
C                                    shells x = A,B and P=A+B
C                    ATOMIC       =  indicates, if purely atomic
C                                    integrals will be evaluated
C                    NEXP         =  current # of exponent pairs
C                    Ex           =  current MIJ exponents from centers
C                                    x = A,B
C                    E2AB         =  current MIJ double exponent
C                                    products between centers A and B
C                    OVL1Dx       =  current overlap 1D integrals for
C                                    each cartesian component (x=X,Y,Z)
C
C                  Output:
C
C                    KIN1Dx       =  current kinetic 1D integrals for
C                                    each cartesian component (x=X,Y,Z)
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT  NONE

         LOGICAL   ATOMIC

         INTEGER   I,J,N
         INTEGER   IM1,IP1,JM1,JP1
         INTEGER   NEXP
         INTEGER   SHELLA,SHELLB,SHELLP

         DOUBLE PRECISION  HALF,HALFJ,HALFIJ
         DOUBLE PRECISION  XI,XJ,XIB,XJA,X2AB

         DOUBLE PRECISION  EA     (1:NEXP)
         DOUBLE PRECISION  EB     (1:NEXP)
         DOUBLE PRECISION  E2AB   (1:NEXP)

         DOUBLE PRECISION  KIN1DX (1:NEXP,0:SHELLA  ,0:SHELLB  )
         DOUBLE PRECISION  KIN1DY (1:NEXP,0:SHELLA  ,0:SHELLB  )
         DOUBLE PRECISION  KIN1DZ (1:NEXP,0:SHELLA  ,0:SHELLB  )
         DOUBLE PRECISION  OVL1DX (1:NEXP,0:SHELLP+2,0:SHELLB+1)
         DOUBLE PRECISION  OVL1DY (1:NEXP,0:SHELLP+2,0:SHELLB+1)
         DOUBLE PRECISION  OVL1DZ (1:NEXP,0:SHELLP+2,0:SHELLB+1)

         PARAMETER  (HALF = 0.5D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...pure atomic/nonatomic s-shell integrals first.
C
C
         DO 10 N = 1,NEXP
            X2AB = E2AB (N)
            KIN1DX (N,0,0) = X2AB * OVL1DX (N,1,1)
            KIN1DY (N,0,0) = X2AB * OVL1DY (N,1,1)
            KIN1DZ (N,0,0) = X2AB * OVL1DZ (N,1,1)
   10    CONTINUE

         IF (SHELLP.EQ.0) RETURN
C
C
C             ...check atomic or nonatomic case for extra integrals
C                (if any).
C
C
         IF (ATOMIC) THEN
C
C
C             ...atomic even-even column s-shell integrals.
C
C
             DO 100 I = 2,SHELLA,2
                IM1 = I - 1
                IP1 = I + 1
                XI = DFLOAT (I)

                DO 102 N = 1,NEXP
                   XIB = XI * EB (N)
                   X2AB = E2AB (N)
                   KIN1DX (N,I,0)  =  -  XIB * OVL1DX (N,IM1,1)
     +                                + X2AB * OVL1DX (N,IP1,1)
                   KIN1DY (N,I,0)  =  -  XIB * OVL1DY (N,IM1,1)
     +                                + X2AB * OVL1DY (N,IP1,1)
                   KIN1DZ (N,I,0)  =  -  XIB * OVL1DZ (N,IM1,1)
     +                                + X2AB * OVL1DZ (N,IP1,1)
  102           CONTINUE
  100        CONTINUE
C
C
C             ...atomic odd-odd column non s-shell integrals.
C
C
             DO 110 J = 1,SHELLB,2
                JM1 = J - 1
                JP1 = J + 1
                XJ = DFLOAT (J)
                HALFJ = HALF * XJ

                DO 112 I = 1,SHELLA,2
                   IM1 = I - 1
                   IP1 = I + 1
                   XI = DFLOAT (I)
                   HALFIJ = HALFJ * XI

                   DO 114 N = 1,NEXP
                      XIB = XI * EB (N)
                      XJA = XJ * EA (N)
                      X2AB = E2AB (N)

                      KIN1DX (N,I,J)  =  HALFIJ * OVL1DX (N,IM1,JM1)
     +                                   -  XIB * OVL1DX (N,IM1,JP1)
     +                                   -  XJA * OVL1DX (N,IP1,JM1)
     +                                   + X2AB * OVL1DX (N,IP1,JP1)
                      KIN1DY (N,I,J)  =  HALFIJ * OVL1DY (N,IM1,JM1)
     +                                   -  XIB * OVL1DY (N,IM1,JP1)
     +                                   -  XJA * OVL1DY (N,IP1,JM1)
     +                                   + X2AB * OVL1DY (N,IP1,JP1)
                      KIN1DZ (N,I,J)  =  HALFIJ * OVL1DZ (N,IM1,JM1)
     +                                   -  XIB * OVL1DZ (N,IM1,JP1)
     +                                   -  XJA * OVL1DZ (N,IP1,JM1)
     +                                   + X2AB * OVL1DZ (N,IP1,JP1)
  114              CONTINUE
  112           CONTINUE
  110        CONTINUE
C
C
C             ...atomic even-even column non s-shell integrals.
C
C
             DO 120 J = 2,SHELLB,2
                JM1 = J - 1
                JP1 = J + 1
                XJ = DFLOAT (J)

                DO 122 N = 1,NEXP
                   XJA = XJ * EA (N)
                   X2AB = E2AB (N)
                   KIN1DX (N,0,J)  =  -  XJA * OVL1DX (N,1,JM1)
     +                                + X2AB * OVL1DX (N,1,JP1)
                   KIN1DY (N,0,J)  =  -  XJA * OVL1DY (N,1,JM1)
     +                                + X2AB * OVL1DY (N,1,JP1)
                   KIN1DZ (N,0,J)  =  -  XJA * OVL1DZ (N,1,JM1)
     +                                + X2AB * OVL1DZ (N,1,JP1)
  122           CONTINUE

                HALFJ = HALF * XJ
                DO 124 I = 2,SHELLA,2
                   IM1 = I - 1
                   IP1 = I + 1
                   XI = DFLOAT (I)
                   HALFIJ = HALFJ * XI
                   DO 126 N = 1,NEXP
                      XIB = XI * EB (N)
                      XJA = XJ * EA (N)
                      X2AB = E2AB (N)

                      KIN1DX (N,I,J)  =  HALFIJ * OVL1DX (N,IM1,JM1)
     +                                   -  XIB * OVL1DX (N,IM1,JP1)
     +                                   -  XJA * OVL1DX (N,IP1,JM1)
     +                                   + X2AB * OVL1DX (N,IP1,JP1)
                      KIN1DY (N,I,J)  =  HALFIJ * OVL1DY (N,IM1,JM1)
     +                                   -  XIB * OVL1DY (N,IM1,JP1)
     +                                   -  XJA * OVL1DY (N,IP1,JM1)
     +                                   + X2AB * OVL1DY (N,IP1,JP1)
                      KIN1DZ (N,I,J)  =  HALFIJ * OVL1DZ (N,IM1,JM1)
     +                                   -  XIB * OVL1DZ (N,IM1,JP1)
     +                                   -  XJA * OVL1DZ (N,IP1,JM1)
     +                                   + X2AB * OVL1DZ (N,IP1,JP1)
  126              CONTINUE
  124           CONTINUE
  120        CONTINUE

         ELSE
C
C
C             ...nonatomic column s-shell integrals.
C
C
             DO 200 I = 1,SHELLA
                IM1 = I - 1
                IP1 = I + 1
                XI = DFLOAT (I)
                DO 202 N = 1,NEXP
                   XIB = XI * EB (N)
                   X2AB = E2AB (N)
                   KIN1DX (N,I,0)  =  -  XIB * OVL1DX (N,IM1,1)
     +                                + X2AB * OVL1DX (N,IP1,1)
                   KIN1DY (N,I,0)  =  -  XIB * OVL1DY (N,IM1,1)
     +                                + X2AB * OVL1DY (N,IP1,1)
                   KIN1DZ (N,I,0)  =  -  XIB * OVL1DZ (N,IM1,1)
     +                                + X2AB * OVL1DZ (N,IP1,1)
  202           CONTINUE
  200        CONTINUE
C
C
C             ...nonatomic column non s-shell integrals.
C
C
             DO 210 J = 1,SHELLB
                JM1 = J - 1
                JP1 = J + 1
                XJ = DFLOAT (J)

                DO 212 N = 1,NEXP
                   XJA = XJ * EA (N)
                   X2AB = E2AB (N)
                   KIN1DX (N,0,J)  =  -  XJA * OVL1DX (N,1,JM1)
     +                                + X2AB * OVL1DX (N,1,JP1)
                   KIN1DY (N,0,J)  =  -  XJA * OVL1DY (N,1,JM1)
     +                                + X2AB * OVL1DY (N,1,JP1)
                   KIN1DZ (N,0,J)  =  -  XJA * OVL1DZ (N,1,JM1)
     +                                + X2AB * OVL1DZ (N,1,JP1)
  212           CONTINUE

                HALFJ = HALF * XJ
                DO 214 I = 1,SHELLA
                   IM1 = I - 1
                   IP1 = I + 1
                   XI = DFLOAT (I)
                   HALFIJ = HALFJ * XI
                   DO 216 N = 1,NEXP
                      XIB = XI * EB (N)
                      XJA = XJ * EA (N)
                      X2AB = E2AB (N)

                      KIN1DX (N,I,J)  =  HALFIJ * OVL1DX (N,IM1,JM1)
     +                                   -  XIB * OVL1DX (N,IM1,JP1)
     +                                   -  XJA * OVL1DX (N,IP1,JM1)
     +                                   + X2AB * OVL1DX (N,IP1,JP1)
                      KIN1DY (N,I,J)  =  HALFIJ * OVL1DY (N,IM1,JM1)
     +                                   -  XIB * OVL1DY (N,IM1,JP1)
     +                                   -  XJA * OVL1DY (N,IP1,JM1)
     +                                   + X2AB * OVL1DY (N,IP1,JP1)
                      KIN1DZ (N,I,J)  =  HALFIJ * OVL1DZ (N,IM1,JM1)
     +                                   -  XIB * OVL1DZ (N,IM1,JP1)
     +                                   -  XJA * OVL1DZ (N,IP1,JM1)
     +                                   + X2AB * OVL1DZ (N,IP1,JP1)
  216              CONTINUE
  214           CONTINUE
  210        CONTINUE

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
