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
         SUBROUTINE  OED__XYZ_MOM_INTEGRALS
     +
     +                    ( SHELLP,SHELLB,
     +                      MOMENT,
     +                      ATOMIC,
     +                      NEXP,
     +                      PAM,
     +                      PINVHF,
     +                      AB,DA,DP,
     +                      OLDMOM,
     +
     +                               MOMINT )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__XYZ_MOM_INTEGRALS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation calculates a full table of 1-dimensional
C                moment integrals of x,y or z type using the Rys
C                vertical recurrence scheme RI explained below.
C
C                The vertical recurrence scheme RI is due to Rys, Dupuis
C                and King, J. Comp. Chem. 4, p.154-157 (1983).
C
C
C                   For I = 1,...,SHELLP-1
C                       INT1D (I+1,0) = I * (1/2P) * INT1D  (I-1,0)
C                                             + PA * INT1D  (I,0)
C                                        MU (1/2P) * OLDMOM (I,0)
C
C
C                If SHELLB is > 0, then the horizontal transfer scheme
C                HI is also applied:
C
C                   For J = 1,...,SHELLB
C                   For I = 0,...,SHELLP-J
C                       INT1D (I,J) =       INT1D (I+1,J-1)
C                                    + AB * INT1D (I,J-1)
C
C                If the HI is not needed, the values in ABX,ABY and ABZ
C                have no meaning and can be set as wished.
C
C                The 1D integrals are calculated for all exponent pairs
C                simultaneously and placed into a 3-dimensional array
C                with the exponent pair index varying fastest.
C
C
C                  Input:
C
C                    SHELLx       =  the shell types for shells
C                                    x=P=A+B and x=B
C                    ATOMIC       =  indicates, if purely atomic 1D
C                                    integrals will be evaluated
C                    NEXP         =  # of exponent pairs
C                    PAx          =  current NEXP coordinate x=X,Y,Z
C                                    differences P-A between centers
C                                    P and A
C                    xP           =  current NEXP coordinate x=X,Y,Z
C                                    Gaussian Product centers
C                    PINVHF       =  current NEXP values of 1/(2*P),
C                                    where P are the exponent sums
C                                    for contraction shells A and B
C                    ABx          =  coordinate x=X,Y,Z differences
C                                    A-B between centers A and B
C                    MOMENT       =  power of the integral to compute
C                    OLDMOM       =  all 1d overlap integrals for
C                                    the moment of interest
C
C                  Output:
C
C                    MOMINT       =  all 1D moment integrals for the
C                                    cartesian component of interest
C                                    for the current exponent pairs
C
C
C  AUTHOR      : Norbert Flocke
C                  - Wrote original OED Package
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
         INTEGER    NEXP,MOMENT
         INTEGER    SHELLP,SHELLB

         DOUBLE PRECISION  AB,DA
         DOUBLE PRECISION  F,G
         DOUBLE PRECISION  ZERO,ONE,TWO
         DOUBLE PRECISION  P1,P2

         DOUBLE PRECISION  PAM    (1:NEXP)
         DOUBLE PRECISION  PINVHF (1:NEXP)
         DOUBLE PRECISION  DP     (1:NEXP)

         DOUBLE PRECISION  MOMINT (1:NEXP,0:SHELLP,0:SHELLB)
         DOUBLE PRECISION  OLDMOM (1:NEXP,0:SHELLP,0:SHELLB)

         PARAMETER  (ZERO  =  0.D0)
         PARAMETER  (ONE   =  1.D0)
         PARAMETER  (TWO   =  2.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...Now that we have (a|mu-1|b) for the atomic case,
C                stored in OLDMOM, we can begin computing
C                (a|mu|b) since the recursion relations are as
C                follows:
C
C                  (I+1|mu|0) =        I * (1/2P) * (I-1| mu |0)
C                                            + PA * (I  | mu |0)
C                               + MOMENT * (1/2P) * (I  |mu-1|0)
C
C
C                If SHELLB is > 0, then the horizontal transfer scheme
C                HI is also applied:
C
C                   For J = 1,...,SHELLB
C                   For I = 0,...,SHELLP-J
C                       (I|X|J) = (I+1|X|J-1) + AB * (I|X|J-1)
C
C
C-----------------------------------------------------------------------
C
C
C             ...initialize the (1|0) integrals
C
C
         DO N = 1,NEXP
            P1 = MOMENT * PINVHF (N)
            MOMINT (N,1,0) = PAM (N)*MOMINT (N,0,0) + P1*OLDMOM (N,0,0)
         END DO
C
C
C             ...begin the recursion relations.
C
C
         IF (ATOMIC) THEN
             F = ONE
             G = MOMENT
             DO 200 I = 2,SHELLP
                I1 = I - 1
                I2 = I - 2
                DO 210 N = 1,NEXP
                   P1 = F * PINVHF (N)
                   P2 = G * PINVHF (N)
                   MOMINT (N,I,0) =   P1 * MOMINT (N,I2,0)
     +                              + P2 * OLDMOM (N,I1,0)
  210           CONTINUE
                F = F + ONE
  200        CONTINUE

             IF (SHELLB.GT.0) THEN

                 DO 220 J = 1,SHELLB
                    J1 = J - 1
                    IEND = SHELLP-J
                    DO 230 I = 0,IEND
                       I1 = I + 1
                       DO 240 N = 1,NEXP
                          MOMINT (N,I,J) = MOMINT (N,I1,J1)
  240                  CONTINUE
  230               CONTINUE
  220            CONTINUE

             END IF

         ELSE

             F = ONE
             G = MOMENT
             DO 250 I = 2,SHELLP
                I1 = I - 1
                I2 = I - 2
                DO 260 N = 1,NEXP
                   P1 = F * PINVHF (N)
                   P2 = G * PINVHF (N)
                   MOMINT (N,I,0) =        P1 * MOMINT (N,I2,0)
     +                              + PAM (N) * MOMINT (N,I1,0)
     +                                   + P2 * OLDMOM (N,I1,0)
  260           CONTINUE
                F = F + ONE
  250        CONTINUE

             IF (SHELLB.GT.0) THEN

                 DO 270 J = 1,SHELLB
                    J1 = J - 1
                    IEND = SHELLP-J
                    DO 280 I = 0,IEND
                       I1 = I + 1
                       DO 290 N = 1,NEXP
                          MOMINT (N,I,J) =  AB * MOMINT (N,I,J1)
     +                                         + MOMINT (N,I1,J1)
  290                  CONTINUE
  280               CONTINUE
  270            CONTINUE

             END IF

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
