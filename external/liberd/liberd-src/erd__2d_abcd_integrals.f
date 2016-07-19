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
         SUBROUTINE  ERD__2D_ABCD_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,SHELLX,SHELLY,
     +                      SHELLA,SHELLB,SHELLC,SHELLD,
     +                      NGQEXQ,
     +                      WTS,
     +                      B00,B01,B10,C00,D00,
     +                      AB,CD,
     +                      WTAKE,
     +                      CASE2D,
     +                      INTSCR,
     +
     +                               INT2D )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__2D_ABCD_INTEGRALS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation calculates a full table of 2D ABCD
C                integrals using the Rys vertical recurrence scheme
C                VRR and the Rys horizontal transfer scheme explained
C                below.
C
C                Scheme:
C
C                    i) VRR => generate the (P,0,Q,0) integrals
C                   ii) HRR => generate the (A,B,C,D) integrals
C
C                If WTAKE is set true, the Rys weight is being
C                multiplied to the 2D ABCD integral to reduce overall
C                FLOP count. Note, that the Rys weight factor needs
C                to be introduced only three times for the starting
C                2D ABCD integrals for the VRR recurrence scheme,
C                namely to the (0,0,0,0), (1,0,0,0) and (0,0,1,0)
C                elements. The weight factor is then automatically
C                propagated through the vertical and horizontal
C                transfer equations (see below).
C
C                The recurrence schemes VRR and HRR is due to Rys,
C                Dupuis and King, J. Comp. Chem. 4, p.154-157 (1983).
C
C
C                 i) VRR Scheme:
C
C
C                 INT2D (0,0,0,0) = WEIGHT
C                 INT2D (1,0,0,0) = C00 * WEIGHT
C                 INT2D (0,0,1,0) = D00 * WEIGHT
C
C                 For I = 1,...,SHELLP-1
C                     INT2D (I+1,0,0,0) = I * B10 * INT2D (I-1,0,0,0)
C                                           + C00 * INT2D (I,0,0,0)
C                 For K = 1,...,SHELLQ-1
C                     INT2D (0,0,K+1,0) = K * B01 * INT2D (0,0,K-1,0)
C                                           + D00 * INT2D (0,0,K,0)
C                 For I = 1,...,SHELLP
C                     INT2D (I,0,1,0)   = I * B00 * INT2D (I-1,0,0,0)
C                                           + D00 * INT2D (I,0,0,0)
C                 For K = 2,...,SHELLQ
C                     INT2D (1,0,K,0)   = K * B00 * INT2D (0,0,K-1,0)
C                                           + C00 * INT2D (0,0,K,0)
C                 For K = 2,...,SHELLQ
C                 For I = 2,...,SHELLP
C                     INT2D (I,0,K,0) = (I-1) * B10 * INT2D (I-2,0,K,0)
C                                         + K * B00 * INT2D (I-1,0,K-1,0)
C                                             + C00 * INT2D (I-1,0,K,0)
C
C
C                ii) HRR Scheme:
C
C
C                   For K = 0,SHELLQ
C                   For J = 1,...,SHELLX   (SHELLX = Min (SHELLA,SHELLB)
C                   For I = SHELLP-J,...,0
C                       INT2D (I,J,K,0) =  +/- AB * INT2D (I,J-1,K,0)
C                                                 + INT2D (I+1,J-1,K,0)
C
C                   For L = 1,...,SHELLY   (SHELLY = Min (SHELLC,SHELLD)
C                   For K = SHELLQ-L,...,0
C                   For J = 0,...,SHELLB
C                   For I = 0,...,SHELLA
C                       INT2D (I,J,K,L) =  +/- CD * INT2D (I,J,K,L-1)
C                                                 + INT2D (I,J,K+1,L-1)
C
C
C                The 2D integrals are calculated for all roots (info
C                already present in transmitted VRR coefficients!) and
C                for all exponent quadruples simultaneously and placed
C                into a 5-dimensional array with the combined root and
C                exponent quadruplet index varying fastest.
C
C                An important feature of this routine is that it can
C                be called with shell magnitudes A < B and/or C < D.
C                If these two cases happen, the VRR and HRR have to
C                be processed differently. Lets examine in detail the
C                A < B case (analogous is the C < D situation). If
C                shell A < shell B, the HRR must be defined such that
C                P,0 -> B,A due to efficiency and numerical stability
C                reasons. This in turn needs a redefinition of the
C                C00 coefficients, which now must be such that they
C                are defined with the P shell accumulated on center B
C                instead of A. The C00 coefficients are (example x
C                coordinate):
C
C                  i) P accumulated on center A:
C
C                        C00x = (Px - Ax) + terms independent of Ax,Bx
C
C                 ii) P accumulated on center B:
C
C                        C00x = (Px - Bx) + terms independent of Ax,Bx
C
C                When entering the present routine however, all we have
C                are the C00 values based on case i). If case ii)
C                applies, we have to form new C00 via:
C
C                          C00 -> C00 + (A - B) = C00 + AB
C
C                using the cartesian coordinate differences AB = A - B
C                between centers A and B. The constant AB is thus added
C                to all C00 values presently transmitted. Same applies
C                for D00 -> D00 + CD if shell C < shell D. The B00,B01
C                and B10 coefficients are independent of Ax and Bx and
C                thus need not be modified.
C
C
C                  Input:
C
C                    SHELLx      =  maximum shell sum for electron
C                                   1 and 2 (x = P=A+B,Q=C+D)
C                    SHELLx      =  minimum shell type between A and
C                                   B and between C and D for electron
C                                   1 and 2 (x = X,Y)
C                    SHELLy      =  maximum shell type for electron
C                                   1 on sites y = A,B and electron
C                                   2 on sites y = C,D
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    WTS         =  all quadrature Rys weights
C                    B00,B01,B10 =  VRR expansion coefficients
C                                   (cartesian coordinate independent)
C                    C00,D00     =  cartesian coordinate dependent
C                                   VRR expansion coefficients based
C                                   on centers A and C
C                    AB,CD       =  cartesian coordinate differences
C                                   A - B and C - D between sites A
C                                   and B and between sites C and D,
C                                   respectively
C                    WTAKE       =  if true, the Rys weights will be
C                                   build into the 2D ABCD integrals;
C                                   if false, the Rys weights will not
C                                   be considered
C                    CASE2D      =  logical flag for simplifications
C                                   in 2D integral evaluation for
C                                   low quantum numbers
C                    INTSCR      =  scratch array to perform the VRR
C                                   and HRR schemes.
C
C
C                  Output:
C
C                    INT2D       =  all 2D ABCD integrals with reduced
C                                   PB -> AB and QD -> CD dimensions.
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

         LOGICAL   COPYIJ,COPYKL
         LOGICAL   HRRBA,HRRDC
         LOGICAL   WTAKE

         INTEGER   CASE2D
         INTEGER   I,J,K,L,N
         INTEGER   I1,I2,IP1,JM1,K1,K2,KP1,LM1
         INTEGER   NGQEXQ
         INTEGER   SHELLA,SHELLB,SHELLC,SHELLD
         INTEGER   SHELLI,SHELLJ
         INTEGER   SHELLP,SHELLQ,SHELLX,SHELLY

         DOUBLE PRECISION  AB,CD
         DOUBLE PRECISION  F,F1,F2
         DOUBLE PRECISION  WEIGHT
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  B00  (1:NGQEXQ)
         DOUBLE PRECISION  B01  (1:NGQEXQ)
         DOUBLE PRECISION  B10  (1:NGQEXQ)
         DOUBLE PRECISION  C00  (1:NGQEXQ)
         DOUBLE PRECISION  D00  (1:NGQEXQ)
         DOUBLE PRECISION  WTS  (1:NGQEXQ)

         DOUBLE PRECISION  INT2D  (1:NGQEXQ,0:SHELLA,
     +                                      0:SHELLB,
     +                                      0:SHELLC,
     +                                      0:SHELLD)
         DOUBLE PRECISION  INTSCR (1:NGQEXQ,0:SHELLP,
     +                                      0:SHELLX,
     +                                      0:SHELLQ,
     +                                      0:SHELLY)

         PARAMETER  (ONE  = 1.D0)
         PARAMETER  (ZERO = 0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...jump according to the 4 different cases that can arise:
C
C                  P-shell = s- or higher angular momentum
C                  Q-shell = s- or higher angular momentum
C
C                each leading to simplifications in the VRR and HRR
C                formulas. Note, that the HRR is only required if the
C                P-shell and/or Q-shell is higher than s-shell. The
C                case present has been evaluated outside this routine
C                and is transmitted via argument.
C
C
         GOTO (1,2,3,4) CASE2D
C
C
C             ...the case P = s-shell and Q = s-shell (no HRR!).
C                i) VRR => Evaluate: I = 0
C                                    J = 0
C                                    K = 0
C                                    L = 0
C
C
    1    IF (WTAKE) THEN
             DO 100 N = 1,NGQEXQ
                INT2D (N,0,0,0,0) = WTS (N)
  100        CONTINUE
         ELSE
             DO 110 N = 1,NGQEXQ
                INT2D (N,0,0,0,0) = ONE
  110        CONTINUE
         END IF

         RETURN
C
C
C             ...the cases P = s-shell and Q >= p-shell.
C                i) VRR => Evaluate: I = 0
C                                    J = 0
C                                    K = 0,SHELLQ
C                                    L = 0
C
C
    2    HRRDC = SHELLC .LT. SHELLD

         IF (HRRDC .AND. CD.NE.ZERO) THEN
             DO 200 N = 1,NGQEXQ
                D00 (N) = D00 (N) + CD
  200        CONTINUE
         END IF

         IF (WTAKE) THEN
             DO 210 N = 1,NGQEXQ
                WEIGHT = WTS (N)
                INTSCR (N,0,0,0,0) = WEIGHT
                INTSCR (N,0,0,1,0) = D00 (N) * WEIGHT
  210        CONTINUE
         ELSE
             DO 212 N = 1,NGQEXQ
                INTSCR (N,0,0,0,0) = ONE
                INTSCR (N,0,0,1,0) = D00 (N)
  212        CONTINUE
         END IF

         F = ONE
         DO 220 K = 2,SHELLQ
            K1 = K - 1
            K2 = K - 2
            DO 222 N = 1,NGQEXQ
               INTSCR (N,0,0,K,0) = F * B01 (N) * INTSCR (N,0,0,K2,0)
     +                                + D00 (N) * INTSCR (N,0,0,K1,0)
  222       CONTINUE
            F = F + ONE
  220    CONTINUE
C
C
C             ...ii) HRR => Evaluate: I = 0
C                                     J = 0
C                                     K = 0,SHELLQ-SHELLC/SHELLD
C                                     L = 0,SHELLD/SHELLC
C                           and
C
C               iii) Copy INTSCR to INT2D: I = 0
C                                          J = 0
C                                          K = 0,SHELLC
C                                          L = 0,SHELLD
C
C
         IF (.NOT.HRRDC) THEN

             IF (CD.EQ.ZERO) THEN
                 DO 230 L = 1,SHELLD
                    LM1 = L - 1
                    DO 231 K = 0,SHELLQ-L
                       KP1 = K + 1
                       DO 232 N = 1,NGQEXQ
                          INTSCR (N,0,0,K,L) = INTSCR (N,0,0,KP1,LM1)
  232                  CONTINUE
  231               CONTINUE
  230            CONTINUE
             ELSE
                 DO 233 L = 1,SHELLD
                    LM1 = L - 1
                    DO 234 K = 0,SHELLQ-L
                       KP1 = K + 1
                       DO 235 N = 1,NGQEXQ
                          INTSCR (N,0,0,K,L) = INTSCR (N,0,0,KP1,LM1)
     +                                  + CD * INTSCR (N,0,0,K,LM1)
  235                  CONTINUE
  234               CONTINUE
  233            CONTINUE
             END IF

             DO 236 L = 0,SHELLD
             DO 236 K = 0,SHELLC
             DO 236 N = 1,NGQEXQ
                INT2D (N,0,0,K,L) = INTSCR (N,0,0,K,L)
  236        CONTINUE

         ELSE

             IF (CD.EQ.ZERO) THEN
                 DO 240 L = 1,SHELLC
                    LM1 = L - 1
                    DO 241 K = 0,SHELLQ-L
                       KP1 = K + 1
                       DO 242 N = 1,NGQEXQ
                          INTSCR (N,0,0,K,L) = INTSCR (N,0,0,KP1,LM1)
  242                  CONTINUE
  241               CONTINUE
  240            CONTINUE
             ELSE
                 DO 243 L = 1,SHELLC
                    LM1 = L - 1
                    DO 244 K = 0,SHELLQ-L
                       KP1 = K + 1
                       DO 245 N = 1,NGQEXQ
                          INTSCR (N,0,0,K,L) = INTSCR (N,0,0,KP1,LM1)
     +                                  - CD * INTSCR (N,0,0,K,LM1)
  245                  CONTINUE
  244               CONTINUE
  243            CONTINUE
             END IF

             DO 246 L = 0,SHELLD
             DO 246 K = 0,SHELLC
             DO 246 N = 1,NGQEXQ
                INT2D (N,0,0,K,L) = INTSCR (N,0,0,L,K)
  246        CONTINUE

         END IF

         RETURN
C
C
C             ...the cases P >= p-shell and Q = s-shell.
C                i) VRR => Evaluate: I = 0,SHELLP
C                                    J = 0
C                                    K = 0
C                                    L = 0
C
C
    3    HRRBA = SHELLA .LT. SHELLB

         IF (HRRBA .AND. AB.NE.ZERO) THEN
             DO 300 N = 1,NGQEXQ
                C00 (N) = C00 (N) + AB
  300        CONTINUE
         END IF

         IF (WTAKE) THEN
             DO 310 N = 1,NGQEXQ
                WEIGHT = WTS (N)
                INTSCR (N,0,0,0,0) = WEIGHT
                INTSCR (N,1,0,0,0) = C00 (N) * WEIGHT
  310        CONTINUE
         ELSE
             DO 312 N = 1,NGQEXQ
                INTSCR (N,0,0,0,0) = ONE
                INTSCR (N,1,0,0,0) = C00 (N)
  312        CONTINUE
         END IF

         F = ONE
         DO 320 I = 2,SHELLP
            I1 = I - 1
            I2 = I - 2
            DO 322 N = 1,NGQEXQ
               INTSCR (N,I,0,0,0) = F * B10 (N) * INTSCR (N,I2,0,0,0)
     +                                + C00 (N) * INTSCR (N,I1,0,0,0)
  322       CONTINUE
            F = F + ONE
  320    CONTINUE
C
C
C             ...ii) HRR => Evaluate: I = 0,SHELLP-SHELLA/SHELLB
C                                     J = 0,SHELLB/SHELLA
C                                     K = 0
C                                     L = 0
C                           and
C
C               iii) Copy INTSCR to INT2D: I = 0,SHELLA
C                                          J = 0,SHELLB
C                                          K = 0
C                                          L = 0
C
C
         IF (.NOT.HRRBA) THEN

             IF (AB.EQ.ZERO) THEN
                 DO 330 J = 1,SHELLB
                    JM1 = J - 1
                    DO 331 I = 0,SHELLP-J
                       IP1 = I + 1
                       DO 332 N = 1,NGQEXQ
                          INTSCR (N,I,J,0,0) = INTSCR (N,IP1,JM1,0,0)
  332                  CONTINUE
  331               CONTINUE
  330            CONTINUE
             ELSE
                 DO 333 J = 1,SHELLB
                    JM1 = J - 1
                    DO 334 I = 0,SHELLP-J
                       IP1 = I + 1
                       DO 335 N = 1,NGQEXQ
                          INTSCR (N,I,J,0,0) = INTSCR (N,IP1,JM1,0,0)
     +                                  + AB * INTSCR (N,I,JM1,0,0)
  335                  CONTINUE
  334               CONTINUE
  333            CONTINUE
             END IF

             DO 336 J = 0,SHELLB
             DO 336 I = 0,SHELLA
             DO 336 N = 1,NGQEXQ
                INT2D (N,I,J,0,0) = INTSCR (N,I,J,0,0)
  336        CONTINUE

         ELSE

             IF (AB.EQ.ZERO) THEN
                 DO 340 J = 1,SHELLA
                    JM1 = J - 1
                    DO 341 I = 0,SHELLP-J
                       IP1 = I + 1
                       DO 342 N = 1,NGQEXQ
                          INTSCR (N,I,J,0,0) = INTSCR (N,IP1,JM1,0,0)
  342                  CONTINUE
  341               CONTINUE
  340            CONTINUE
             ELSE
                 DO 343 J = 1,SHELLA
                    JM1 = J - 1
                    DO 344 I = 0,SHELLP-J
                       IP1 = I + 1
                       DO 345 N = 1,NGQEXQ
                          INTSCR (N,I,J,0,0) = INTSCR (N,IP1,JM1,0,0)
     +                                  - AB * INTSCR (N,I,JM1,0,0)
  345                  CONTINUE
  344               CONTINUE
  343            CONTINUE
             END IF

             DO 346 J = 0,SHELLB
             DO 346 I = 0,SHELLA
             DO 346 N = 1,NGQEXQ
                INT2D (N,I,J,0,0) = INTSCR (N,J,I,0,0)
  346        CONTINUE

         END IF

         RETURN
C
C
C             ...the cases P >= p-shell and Q >= p-shell.
C                i) VRR => Evaluate: I = 0,SHELLP       I = 0
C                                    J = 0              J = 0
C                                    K = 0         and  K = 0,SHELLQ
C                                    L = 0              L = 0
C
C
    4    HRRBA = SHELLA .LT. SHELLB
         HRRDC = SHELLC .LT. SHELLD

         IF (HRRBA .AND. AB.NE.ZERO) THEN
             DO 400 N = 1,NGQEXQ
                C00 (N) = C00 (N) + AB
  400        CONTINUE
         END IF

         IF (HRRDC .AND. CD.NE.ZERO) THEN
             DO 402 N = 1,NGQEXQ
                D00 (N) = D00 (N) + CD
  402        CONTINUE
         END IF

         IF (WTAKE) THEN
             DO 410 N = 1,NGQEXQ
                WEIGHT = WTS (N)
                INTSCR (N,0,0,0,0) = WEIGHT
                INTSCR (N,1,0,0,0) = C00 (N) * WEIGHT
                INTSCR (N,0,0,1,0) = D00 (N) * WEIGHT
  410        CONTINUE
         ELSE
             DO 411 N = 1,NGQEXQ
                INTSCR (N,0,0,0,0) = ONE
                INTSCR (N,1,0,0,0) = C00 (N)
                INTSCR (N,0,0,1,0) = D00 (N)
  411        CONTINUE
         END IF

         F = ONE
         DO 412 I = 2,SHELLP
            I1 = I - 1
            I2 = I - 2
            DO 413 N = 1,NGQEXQ
               INTSCR (N,I,0,0,0) = F * B10 (N) * INTSCR (N,I2,0,0,0)
     +                                + C00 (N) * INTSCR (N,I1,0,0,0)
  413       CONTINUE
            F = F + ONE
  412    CONTINUE

         F = ONE
         DO 414 K = 2,SHELLQ
            K1 = K - 1
            K2 = K - 2
            DO 415 N = 1,NGQEXQ
               INTSCR (N,0,0,K,0) = F * B01 (N) * INTSCR (N,0,0,K2,0)
     +                                + D00 (N) * INTSCR (N,0,0,K1,0)
  415       CONTINUE
            F = F + ONE
  414    CONTINUE
C
C
C             ...i) VRR => Evaluate (if any): I = 1,SHELLP
C                                             J = 0
C                                             K = 1,SHELLQ
C                                             L = 0
C
C
         IF (SHELLQ.LE.SHELLP) THEN

             F1 = ONE
             DO 420 K = 1,SHELLQ
                K1 = K - 1
                DO 421 N = 1,NGQEXQ
                   INTSCR (N,1,0,K,0) = 
     +                              F1 * B00 (N) * INTSCR (N,0,0,K1,0)
     +                                 + C00 (N) * INTSCR (N,0,0,K,0)
  421           CONTINUE
                F2 = ONE
                DO 422 I = 2,SHELLP
                   I1 = I - 1
                   I2 = I - 2
                   DO 423 N = 1,NGQEXQ
                      INTSCR (N,I,0,K,0) =
     +                              F1 * B00 (N) * INTSCR (N,I1,0,K1,0)
     +                            + F2 * B10 (N) * INTSCR (N,I2,0,K,0)
     +                                 + C00 (N) * INTSCR (N,I1,0,K,0)
  423              CONTINUE
                   F2 = F2 + ONE
  422           CONTINUE
                F1 = F1 + ONE
  420        CONTINUE

         ELSE

             F1 = ONE
             DO 430 I = 1,SHELLP
                I1 = I - 1
                DO 431 N = 1,NGQEXQ
                   INTSCR (N,I,0,1,0) = 
     +                               F1 * B00 (N) * INTSCR (N,I1,0,0,0)
     +                                  + D00 (N) * INTSCR (N,I,0,0,0)
  431           CONTINUE
                F2 = ONE
                DO 432 K = 2,SHELLQ
                   K1 = K - 1
                   K2 = K - 2
                   DO 433 N = 1,NGQEXQ
                      INTSCR (N,I,0,K,0) =
     +                               F1 * B00 (N) * INTSCR (N,I1,0,K1,0)
     +                             + F2 * B01 (N) * INTSCR (N,I,0,K2,0)
     +                                  + D00 (N) * INTSCR (N,I,0,K1,0)
  433              CONTINUE
                   F2 = F2 + ONE
  432           CONTINUE
                F1 = F1 + ONE
  430        CONTINUE

         END IF
C
C
C             ...ii) HRR => Evaluate: I = 0,SHELLP-SHELLA/SHELLB
C                                     J = 0,SHELLB/SHELLA
C                                     K = 0,SHELLQ
C                                     L = 0
C
C
         IF (.NOT.HRRBA) THEN

             IF (AB.EQ.ZERO) THEN
                 DO 440 K = 0,SHELLQ
                 DO 440 J = 1,SHELLB
                    JM1 = J - 1
                    DO 441 I = 0,SHELLP-J
                       IP1 = I + 1
                       DO 442 N = 1,NGQEXQ
                          INTSCR (N,I,J,K,0) = INTSCR (N,IP1,JM1,K,0)
  442                  CONTINUE
  441               CONTINUE
  440            CONTINUE
             ELSE
                 DO 443 K = 0,SHELLQ
                 DO 443 J = 1,SHELLB
                    JM1 = J - 1
                    DO 444 I = 0,SHELLP-J
                       IP1 = I + 1
                       DO 445 N = 1,NGQEXQ
                          INTSCR (N,I,J,K,0) = INTSCR (N,IP1,JM1,K,0)
     +                                  + AB * INTSCR (N,I,JM1,K,0)
  445                  CONTINUE
  444               CONTINUE
  443            CONTINUE
             END IF

             COPYIJ = .TRUE.
             SHELLI = SHELLA
             SHELLJ = SHELLB

         ELSE

             IF (AB.EQ.ZERO) THEN
                 DO 450 K = 0,SHELLQ
                 DO 450 J = 1,SHELLA
                    JM1 = J - 1
                    DO 451 I = 0,SHELLP-J
                       IP1 = I + 1
                       DO 452 N = 1,NGQEXQ
                          INTSCR (N,I,J,K,0) = INTSCR (N,IP1,JM1,K,0)
  452                  CONTINUE
  451               CONTINUE
  450            CONTINUE
             ELSE
                 DO 453 K = 0,SHELLQ
                 DO 453 J = 1,SHELLA
                    JM1 = J - 1
                    DO 454 I = 0,SHELLP-J
                       IP1 = I + 1
                       DO 455 N = 1,NGQEXQ
                          INTSCR (N,I,J,K,0) = INTSCR (N,IP1,JM1,K,0)
     +                                  - AB * INTSCR (N,I,JM1,K,0)
  455                  CONTINUE
  454               CONTINUE
  453            CONTINUE
             END IF

             COPYIJ = .FALSE.
             SHELLI = SHELLB
             SHELLJ = SHELLA

         END IF
C
C
C             ...ii) HRR => Evaluate: I = 0,SHELLI
C                                     J = 0,SHELLJ
C                                     K = 0,SHELLQ-SHELLC/SHELLD
C                                     L = 0,SHELLD/SHELLC
C
C
         IF (.NOT.HRRDC) THEN

             IF (CD.EQ.ZERO) THEN
                 DO 460 L = 1,SHELLD
                    LM1 = L - 1
                    DO 461 K = 0,SHELLQ-L
                       KP1 = K + 1
                       DO 462 J = 0,SHELLJ
                       DO 462 I = 0,SHELLI
                       DO 462 N = 1,NGQEXQ
                          INTSCR (N,I,J,K,L) = INTSCR (N,I,J,KP1,LM1)
  462                  CONTINUE
  461               CONTINUE
  460            CONTINUE
             ELSE
                 DO 463 L = 1,SHELLD
                    LM1 = L - 1
                    DO 464 K = 0,SHELLQ-L
                       KP1 = K + 1
                       DO 465 J = 0,SHELLJ
                       DO 465 I = 0,SHELLI
                       DO 465 N = 1,NGQEXQ
                          INTSCR (N,I,J,K,L) = INTSCR (N,I,J,KP1,LM1)
     +                                  + CD * INTSCR (N,I,J,K,LM1)
  465                  CONTINUE
  464               CONTINUE
  463            CONTINUE
             END IF

             COPYKL = .TRUE.

         ELSE

             IF (CD.EQ.ZERO) THEN
                 DO 470 L = 1,SHELLC
                    LM1 = L - 1
                    DO 471 K = 0,SHELLQ-L
                       KP1 = K + 1
                       DO 472 J = 0,SHELLJ
                       DO 472 I = 0,SHELLI
                       DO 472 N = 1,NGQEXQ
                          INTSCR (N,I,J,K,L) = INTSCR (N,I,J,KP1,LM1)
  472                  CONTINUE
  471               CONTINUE
  470            CONTINUE
             ELSE

                 DO 473 L = 1,SHELLC
                    LM1 = L - 1
                    DO 474 K = 0,SHELLQ-L
                       KP1 = K + 1
                       DO 475 J = 0,SHELLJ
                       DO 475 I = 0,SHELLI
                       DO 475 N = 1,NGQEXQ
                          INTSCR (N,I,J,K,L) = INTSCR (N,I,J,KP1,LM1)
     +                                  - CD * INTSCR (N,I,J,K,LM1)
  475                  CONTINUE
  474               CONTINUE
  473            CONTINUE
             END IF

             COPYKL = .FALSE.

         END IF
C
C
C             ...iii) Copy INTSCR to INT2D: I = 0,SHELLA
C                                           J = 0,SHELLB
C                                           K = 0,SHELLC
C                                           L = 0,SHELLD
C
C
         IF (COPYIJ .AND. COPYKL) THEN

             DO 480 L = 0,SHELLD
             DO 480 K = 0,SHELLC
             DO 480 J = 0,SHELLB
             DO 480 I = 0,SHELLA
             DO 480 N = 1,NGQEXQ
                INT2D (N,I,J,K,L) = INTSCR (N,I,J,K,L)
  480        CONTINUE

         ELSE IF (COPYIJ) THEN

             DO 482 L = 0,SHELLD
             DO 482 K = 0,SHELLC
             DO 482 J = 0,SHELLB
             DO 482 I = 0,SHELLA
             DO 482 N = 1,NGQEXQ
                INT2D (N,I,J,K,L) = INTSCR (N,I,J,L,K)
  482        CONTINUE

         ELSE IF (COPYKL) THEN

             DO 484 L = 0,SHELLD
             DO 484 K = 0,SHELLC
             DO 484 J = 0,SHELLB
             DO 484 I = 0,SHELLA
             DO 484 N = 1,NGQEXQ
                INT2D (N,I,J,K,L) = INTSCR (N,J,I,K,L)
  484        CONTINUE

         ELSE

             DO 486 L = 0,SHELLD
             DO 486 K = 0,SHELLC
             DO 486 J = 0,SHELLB
             DO 486 I = 0,SHELLA
             DO 486 N = 1,NGQEXQ
                INT2D (N,I,J,K,L) = INTSCR (N,J,I,L,K)
  486        CONTINUE

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
