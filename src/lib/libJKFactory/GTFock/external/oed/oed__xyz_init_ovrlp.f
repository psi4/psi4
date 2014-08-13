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
         SUBROUTINE  OED__XYZ_INIT_OVRLP
     +
     +                    ( SHELLP,SHELLB,
     +                      ATOMIC,
     +                      NEXP,
     +                      PAD,
     +                      DP,
     +                      PINVHF,
     +                      ABM,
     +
     +                               OVRLAP )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__XYZ_INIT_OVRLP
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation calculates a full table of 1-dimensional
C                overlap integrals of x,y or z type using the Rys
C                vertical recurrence scheme RI explained below.
C
C                The vertical recurrence scheme RI is due to Rys, Dupuis
C                and King, J. Comp. Chem. 4, p.154-157 (1983).
C
C                It's purpose is to generate the first set of integrals
C                used in the recurrence relations for moment integrals
C
C
C                   OVRLAP (0,0) = 1.D0
C                   OVRLAP (1,0) = PA
C
C                   For I = 1,...,SHELLP-1
C                       OVRLAP (I+1,0) = I * (1/2P) * OVRLAP (I-1,0)
C                                              + PA * OVRLAP (I,0)
C
C
C                If SHELLB is > 0, then the horizontal transfer scheme
C                HI is also applied:
C
C                   For J = 1,...,SHELLB
C                   For I = 0,...,SHELLP-J
C                       OVRLAP (I,J) =      OVRLAP (I+1,J-1)
C                                    + AB * OVRLAP (I,J-1)
C
C                If the HI is not needed, the values in ABX,ABY and ABZ
C                have no meaning and can be set as wished.
C
C                The 1D integrals are calculated for all exponent pairs
C                simultaneously and placed into a 3-dimensional array
C                with the exponent pair index varying fastest.
C
C                              !!!IMPORTANT!!!
C
C                In case the 1D integrals are needed for an atom, all
C                even-odd and odd-even shell combinations lead to 1D
C                integrals whose value is zero.
C
C
C                  Input:
C
C                    SHELLx       =  the shell types for shells
C                                    x=P=A+B and x=B
C                    ATOMIC       =  indicates, if purely atomic 1D
C                                    integrals will be evaluated
C                    NEXP         =  # of exponent pairs
C                    PAD          =  current NEXP coordinate D=X,Y,Z
C                                    differences P-A between centers
C                                    P and A
C                    DP           =  current NEXP values of Px, Py, 
C                                    or Pz
C                    PINVHF       =  current NEXP values of 1/(2*P),
C                                    where P are the exponent sums
C                                    for contraction shells A and B
C                    ABM          =  coordinate M=X,Y,Z differences
C                                    A-B between centers A and B
C
C                  Output:
C
C                    OVRLAP       =  1D overlap integrals for one
C                                    cartesian component x = X,Y,Z
C                                    for the current exponent pairs
C
C
C  AUTHOR      : Norbert Flocke
C                  - Wrote original OED package
C
C  MODIFIED    : Thomas Watson Jr.                   p  q  r
C                  - Modified OED package to handle X, Y, Z integrals
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         LOGICAL    ATOMIC

         INTEGER    I,J,N
         INTEGER    I1,I2,J1
         INTEGER    IEND,JEND
         INTEGER    NEXP
         INTEGER    SHELLP,SHELLB

         DOUBLE PRECISION  ABM
         DOUBLE PRECISION  F
         DOUBLE PRECISION  ZERO,ONE,TWO
         DOUBLE PRECISION  P1

         DOUBLE PRECISION  DP     (1:NEXP)
         DOUBLE PRECISION  PAD    (1:NEXP)
         DOUBLE PRECISION  PINVHF (1:NEXP)

         DOUBLE PRECISION  OVRLAP (1:NEXP,0:SHELLP,0:SHELLB)

         PARAMETER  (ZERO  =  0.D0)
         PARAMETER  (ONE   =  1.D0)
         PARAMETER  (TWO   =  2.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...initialize 1D overlap integrals.
C
C
         IF (ATOMIC) THEN

             IF (SHELLP.GT.1) THEN
                 F = ONE
                 DO 200 I = 2,SHELLP,2
                    I2 = I - 2
                    DO 210 N = 1,NEXP
                       P1 = F * PINVHF (N)
                       OVRLAP (N,I,0) = P1 * OVRLAP (N,I2,0)
  210               CONTINUE
                    F = F + TWO
  200            CONTINUE
             END IF

             IF (SHELLB.GT.0) THEN

                 JEND = SHELLB - 1
                 DO 220 J = 1,JEND,2
                    J1 = J - 1
                    IEND = SHELLP - J
                    DO 230 I = 1,IEND,2
                       I1 = I + 1
                       DO 240 N = 1,NEXP
                          OVRLAP (N,I,J) = OVRLAP (N,I1,J1)
  240                  CONTINUE
  230               CONTINUE
                    J1 = J + 1
                    DO 250 I = 0,IEND,2
                       I1 = I + 1
                       DO 260 N = 1,NEXP
                          OVRLAP (N,I,J1) = OVRLAP (N,I1,J)
  260                  CONTINUE
  250               CONTINUE
  220            CONTINUE

                 IF (MOD (SHELLB,2).NE.0) THEN
                     J = SHELLB
                     J1 = J - 1
                     IEND = SHELLP - J
                     DO 270 I = 1,IEND,2
                        I1 = I + 1
                        DO 280 N = 1,NEXP
                           OVRLAP (N,I,J) = OVRLAP (N,I1,J1)
  280                   CONTINUE
  270                CONTINUE
                 END IF

             END IF

         ELSE
C
C
C             ...the general non-atomic case.
C
C
             DO 300 N = 1,NEXP
                OVRLAP (N,1,0) = PAD (N)
  300        CONTINUE

             F = ONE
             DO 310 I = 2,SHELLP
                I1 = I - 1
                I2 = I - 2
                DO 320 N = 1,NEXP
                   P1 = F * PINVHF (N)
                   OVRLAP (N,I,0) =        P1 * OVRLAP (N,I2,0)
     +                              + PAD (N) * OVRLAP (N,I1,0)
  320           CONTINUE
                F = F + ONE
  310        CONTINUE

             IF (SHELLB.GT.0) THEN

                 DO 330 J = 1,SHELLB
                    J1 = J - 1
                    IEND = SHELLP-J
                    DO 340 I = 0,IEND
                       I1 = I + 1
                       DO 350 N = 1,NEXP
                          OVRLAP (N,I,J) =  ABM * OVRLAP (N,I,J1)
     +                                          + OVRLAP (N,I1,J1)
  350                  CONTINUE
  340               CONTINUE
  330            CONTINUE

             END IF

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END

