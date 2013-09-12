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
         SUBROUTINE  ERD__2D_PCD_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,SHELLY,
     +                      SHELLC,SHELLD,
     +                      NGQEXQ,
     +                      WTS,
     +                      B00,B01,B10,C00,D00,
     +                      CD,
     +                      WTAKE,
     +                      CASE2D,
     +                      INTSCR,
     +
     +                               INT2D )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__2D_PCD_INTEGRALS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation calculates a full table of 2D PCD
C                integrals using the Rys vertical recurrence scheme
C                VRR and the Rys horizontal transfer scheme explained
C                below.
C
C                Scheme:
C
C                    i) VRR => generate the (P,Q,0) integrals
C                   ii) HRR => generate the (P,C,D) integrals
C
C                If WTAKE is set true, the Rys weight is being
C                multiplied to the 2D PCD integral to reduce overall
C                FLOP count. Note, that the Rys weight factor needs
C                to be introduced only three times for the starting
C                2D PCD integrals for the VRR recurrence scheme,
C                namely to the (0,0,0), (1,0,0) and (0,1,0) elements.
C                The weight factor is then automatically propagated
C                through the vertical and horizontal transfer
C                equations (see below).
C
C                The recurrence schemes VRR and HRR is due to Rys,
C                Dupuis and King, J. Comp. Chem. 4, p.154-157 (1983).
C
C
C                 i) VRR Scheme:
C
C
C                 INT2D (0,0,0) = WEIGHT
C                 INT2D (1,0,0) = C00 * WEIGHT
C                 INT2D (0,1,0) = D00 * WEIGHT
C
C                 For I = 1,...,SHELLP-1
C                     INT2D (I+1,0,0) = I * B10 * INT2D (I-1,0,0)
C                                         + C00 * INT2D (I,0,0)
C                 For K = 1,...,SHELLQ-1
C                     INT2D (0,K+1,0) = K * B01 * INT2D (0,K-1,0)
C                                         + D00 * INT2D (0,K,0)
C                 For I = 1,...,SHELLP
C                     INT2D (I,1,0)   = I * B00 * INT2D (I-1,0,0)
C                                         + D00 * INT2D (I,0,0)
C                 For K = 2,...,SHELLQ
C                     INT2D (1,K,0)   = K * B00 * INT2D (0,K-1,0)
C                                         + C00 * INT2D (0,K,0)
C                 For K = 2,...,SHELLQ
C                 For I = 2,...,SHELLP
C                     INT2D (I,K,0) = (I-1) * B10 * INT2D (I-2,K,0)
C                                       + K * B00 * INT2D (I-1,K-1,0)
C                                           + C00 * INT2D (I-1,K,0)
C
C
C                ii) HRR Scheme:
C
C
C                   For L = 1,...,SHELLY   (SHELLY = Min (SHELLC,SHELLD)
C                   For K = SHELLQ-L,...,0
C                   For I = 0,...,SHELLP
C                       INT2D (I,K,L) =  +/- CD * INT2D (I,K,L-1)
C                                               + INT2D (I,K+1,L-1)
C
C
C                The 2D integrals are calculated for all roots (info
C                already present in transmitted VRR coefficients!) and
C                for all exponent quadruples simultaneously and placed
C                into a 3-dimensional array with the combined root and
C                exponent quadruplet index varying fastest.
C
C                An important feature of this routine is that it can
C                be called with shell magnitudes C < D. If this is the
C                case, the VRR and HRR have to be processed differently.
C                Lets examine in detail what is different. If C < D,
C                the HRR must be defined such that Q,0 -> D,C due to
C                efficiency and numerical stability reasons. This in
C                turn needs a redefinition of the D00 coefficients,
C                which now must be such that they are defined with the
C                Q shell accumulated on center D instead of C. The D00
C                coefficients are (example x coordinate):
C
C                  i) Q accumulated on center C:
C
C                        D00x = (Qx - Cx) + terms independent of Cx,Dx
C
C                 ii) Q accumulated on center D:
C
C                        D00x = (Qx - Dx) + terms independent of Cx,Dx
C
C                When entering the present routine however, all we have
C                are the D00 values based on case i). If case ii)
C                applies, we have to form new D00 via:
C
C                          D00 -> D00 + (C - D) = D00 + CD
C
C                using the cartesian coordinate differences CD = C - D
C                between centers C and D. The constant CD is thus added
C                to all D00 values presently transmitted. The B00,B01
C                and B10 coefficients are independent of Cx and Dx and
C                thus need not be modified.
C
C
C                  Input:
C
C                    SHELLx      =  maximum shell type for electron
C                                   1 and 2 (x = P=A+B,Q=C+D)
C                    SHELLY      =  minimum shell type between C and D
C                                   for electron 2
C                    SHELLy      =  maximum shell type for electron
C                                   2 on sites y = C,D
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    WTS         =  all quadrature Rys weights
C                    B00,B01,B10 =  VRR expansion coefficients
C                                   (cartesian coordinate independent)
C                    C00,D00     =  cartesian coordinate dependent
C                                   VRR expansion coefficients based
C                                   on centers A and C
C                    CD          =  cartesian coordinate difference
C                                   C - D between sites C and D
C                    WTAKE       =  if true, the Rys weights will be
C                                   build into the 2D PCD integrals;
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
C                    INT2D       =  all 2D PCD integrals with reduced
C                                   QD -> CD dimensions.
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

         LOGICAL   HRRDC
         LOGICAL   WTAKE

         INTEGER   CASE2D
         INTEGER   I,K,L,N
         INTEGER   I1,I2,K1,K2,KP1,LM1
         INTEGER   NGQEXQ
         INTEGER   SHELLP,SHELLQ,SHELLY
         INTEGER   SHELLC,SHELLD

         DOUBLE PRECISION  CD
         DOUBLE PRECISION  F,F1,F2
         DOUBLE PRECISION  WEIGHT
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  B00  (1:NGQEXQ)
         DOUBLE PRECISION  B01  (1:NGQEXQ)
         DOUBLE PRECISION  B10  (1:NGQEXQ)
         DOUBLE PRECISION  C00  (1:NGQEXQ)
         DOUBLE PRECISION  D00  (1:NGQEXQ)
         DOUBLE PRECISION  WTS  (1:NGQEXQ)

         DOUBLE PRECISION  INT2D  (1:NGQEXQ,0:SHELLP,0:SHELLC,0:SHELLD)
         DOUBLE PRECISION  INTSCR (1:NGQEXQ,0:SHELLP,0:SHELLQ,0:SHELLY)

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
C                Q-shell is higher than s-shell. The case present
C                has been evaluated outside this routine and is
C                transmitted via argument.
C
C
         GOTO (1,2,3,4) CASE2D
C
C
C             ...the case P = s-shell and Q = s-shell (no HRR!).
C                i) VRR => Evaluate: I = 0
C                                    K = 0
C                                    L = 0
C
C
    1    IF (WTAKE) THEN
             DO 100 N = 1,NGQEXQ
                INT2D (N,0,0,0) = WTS (N)
  100        CONTINUE
         ELSE
             DO 110 N = 1,NGQEXQ
                INT2D (N,0,0,0) = ONE
  110        CONTINUE
         END IF

         RETURN
C
C
C             ...the cases P = s-shell and Q >= p-shell.
C                i) VRR => Evaluate: I = 0
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
                INTSCR (N,0,0,0) = WEIGHT
                INTSCR (N,0,1,0) = D00 (N) * WEIGHT
  210        CONTINUE
         ELSE
             DO 212 N = 1,NGQEXQ
                INTSCR (N,0,0,0) = ONE
                INTSCR (N,0,1,0) = D00 (N)
  212        CONTINUE
         END IF

         F = ONE
         DO 220 K = 2,SHELLQ
            K1 = K - 1
            K2 = K - 2
            DO 222 N = 1,NGQEXQ
               INTSCR (N,0,K,0) = F * B01 (N) * INTSCR (N,0,K2,0)
     +                              + D00 (N) * INTSCR (N,0,K1,0)
  222       CONTINUE
            F = F + ONE
  220    CONTINUE
C
C
C             ...ii) HRR => Evaluate: I = 0
C                                     K = 0,SHELLQ-SHELLC/SHELLD
C                                     L = 0,SHELLD/SHELLC
C                           and
C
C               iii) Copy INTSCR to INT2D: I = 0
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
                          INTSCR (N,0,K,L) = INTSCR (N,0,KP1,LM1)
  232                  CONTINUE
  231               CONTINUE
  230            CONTINUE
             ELSE
                 DO 233 L = 1,SHELLD
                    LM1 = L - 1
                    DO 234 K = 0,SHELLQ-L
                       KP1 = K + 1
                       DO 235 N = 1,NGQEXQ
                          INTSCR (N,0,K,L) = INTSCR (N,0,KP1,LM1)
     +                                + CD * INTSCR (N,0,K,LM1)
  235                  CONTINUE
  234               CONTINUE
  233            CONTINUE
             END IF

             DO 236 L = 0,SHELLD
             DO 236 K = 0,SHELLC
             DO 236 N = 1,NGQEXQ
                INT2D (N,0,K,L) = INTSCR (N,0,K,L)
  236        CONTINUE

         ELSE

             IF (CD.EQ.ZERO) THEN
                 DO 240 L = 1,SHELLC
                    LM1 = L - 1
                    DO 241 K = 0,SHELLQ-L
                       KP1 = K + 1
                       DO 242 N = 1,NGQEXQ
                          INTSCR (N,0,K,L) = INTSCR (N,0,KP1,LM1)
  242                  CONTINUE
  241               CONTINUE
  240            CONTINUE
             ELSE
                 DO 243 L = 1,SHELLC
                    LM1 = L - 1
                    DO 244 K = 0,SHELLQ-L
                       KP1 = K + 1
                       DO 245 N = 1,NGQEXQ
                          INTSCR (N,0,K,L) = INTSCR (N,0,KP1,LM1)
     +                                - CD * INTSCR (N,0,K,LM1)
  245                  CONTINUE
  244               CONTINUE
  243            CONTINUE
             END IF

             DO 246 L = 0,SHELLD
             DO 246 K = 0,SHELLC
             DO 246 N = 1,NGQEXQ
                INT2D (N,0,K,L) = INTSCR (N,0,L,K)
  246        CONTINUE

         END IF

         RETURN
C
C
C             ...the cases P >= p-shell and Q = s-shell (no HRR!).
C                i) VRR => Evaluate: I = 0,SHELLP
C                                    K = 0
C                                    L = 0
C
C
    3    IF (WTAKE) THEN
             DO 300 N = 1,NGQEXQ
                WEIGHT = WTS (N)
                INT2D (N,0,0,0) = WEIGHT
                INT2D (N,1,0,0) = C00 (N) * WEIGHT
  300        CONTINUE
         ELSE
             DO 310 N = 1,NGQEXQ
                INT2D (N,0,0,0) = ONE
                INT2D (N,1,0,0) = C00 (N)
  310        CONTINUE
         END IF

         F = ONE
         DO 320 I = 2,SHELLP
            I1 = I - 1
            I2 = I - 2
            DO 322 N = 1,NGQEXQ
               INT2D (N,I,0,0) = F * B10 (N) * INT2D (N,I2,0,0)
     +                             + C00 (N) * INT2D (N,I1,0,0)
  322       CONTINUE
            F = F + ONE
  320    CONTINUE

         RETURN
C
C
C             ...the cases P >= p-shell and Q >= p-shell.
C                i) VRR => Evaluate: I = 0,SHELLP       I = 0
C                                    K = 0         and  K = 0,SHELLQ
C                                    L = 0              L = 0
C
C
    4    HRRDC = SHELLC .LT. SHELLD

         IF (HRRDC .AND. CD.NE.ZERO) THEN
             DO 400 N = 1,NGQEXQ
                D00 (N) = D00 (N) + CD
  400        CONTINUE
         END IF

         IF (WTAKE) THEN
             DO 410 N = 1,NGQEXQ
                WEIGHT = WTS (N)
                INTSCR (N,0,0,0) = WEIGHT
                INTSCR (N,1,0,0) = C00 (N) * WEIGHT
                INTSCR (N,0,1,0) = D00 (N) * WEIGHT
  410        CONTINUE
         ELSE
             DO 411 N = 1,NGQEXQ
                INTSCR (N,0,0,0) = ONE
                INTSCR (N,1,0,0) = C00 (N)
                INTSCR (N,0,1,0) = D00 (N)
  411        CONTINUE
         END IF

         F = ONE
         DO 412 I = 2,SHELLP
            I1 = I - 1
            I2 = I - 2
            DO 413 N = 1,NGQEXQ
               INTSCR (N,I,0,0) = F * B10 (N) * INTSCR (N,I2,0,0)
     +                              + C00 (N) * INTSCR (N,I1,0,0)
  413       CONTINUE
            F = F + ONE
  412    CONTINUE

         F = ONE
         DO 414 K = 2,SHELLQ
            K1 = K - 1
            K2 = K - 2
            DO 415 N = 1,NGQEXQ
               INTSCR (N,0,K,0) = F * B01 (N) * INTSCR (N,0,K2,0)
     +                              + D00 (N) * INTSCR (N,0,K1,0)
  415       CONTINUE
            F = F + ONE
  414    CONTINUE
C
C
C             ...i) VRR => Evaluate (if any): I = 1,SHELLP
C                                             K = 1,SHELLQ
C                                             L = 0
C
C
         IF (SHELLQ.LE.SHELLP) THEN

             F1 = ONE
             DO 420 K = 1,SHELLQ
                K1 = K - 1
                DO 421 N = 1,NGQEXQ
                   INTSCR (N,1,K,0) = 
     +                              F1 * B00 (N) * INTSCR (N,0,K1,0)
     +                                 + C00 (N) * INTSCR (N,0,K,0)
  421           CONTINUE
                F2 = ONE
                DO 422 I = 2,SHELLP
                   I1 = I - 1
                   I2 = I - 2
                   DO 423 N = 1,NGQEXQ
                      INTSCR (N,I,K,0) =
     +                              F1 * B00 (N) * INTSCR (N,I1,K1,0)
     +                            + F2 * B10 (N) * INTSCR (N,I2,K,0)
     +                                 + C00 (N) * INTSCR (N,I1,K,0)
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
                   INTSCR (N,I,1,0) = 
     +                               F1 * B00 (N) * INTSCR (N,I1,0,0)
     +                                  + D00 (N) * INTSCR (N,I,0,0)
  431           CONTINUE
                F2 = ONE
                DO 432 K = 2,SHELLQ
                   K1 = K - 1
                   K2 = K - 2
                   DO 433 N = 1,NGQEXQ
                      INTSCR (N,I,K,0) =
     +                               F1 * B00 (N) * INTSCR (N,I1,K1,0)
     +                             + F2 * B01 (N) * INTSCR (N,I,K2,0)
     +                                  + D00 (N) * INTSCR (N,I,K1,0)
  433              CONTINUE
                   F2 = F2 + ONE
  432           CONTINUE
                F1 = F1 + ONE
  430        CONTINUE

         END IF
C
C
C             ...ii) HRR => Evaluate: I = 0,SHELLP
C                                     K = 0,SHELLQ-SHELLC/SHELLD
C                                     L = 0,SHELLD/SHELLC
C                           and
C
C               iii) Copy INTSCR to INT2D: I = 0,SHELLP
C                                          K = 0,SHELLC
C                                          L = 0,SHELLD
C
C
         IF (.NOT.HRRDC) THEN

             IF (CD.EQ.ZERO) THEN
                 DO 440 L = 1,SHELLD
                    LM1 = L - 1
                    DO 441 K = 0,SHELLQ-L
                       KP1 = K + 1
                       DO 442 I = 0,SHELLP
                       DO 442 N = 1,NGQEXQ
                          INTSCR (N,I,K,L) = INTSCR (N,I,KP1,LM1)
  442                  CONTINUE
  441               CONTINUE
  440            CONTINUE
             ELSE
                 DO 443 L = 1,SHELLD
                    LM1 = L - 1
                    DO 444 K = 0,SHELLQ-L
                       KP1 = K + 1
                       DO 445 I = 0,SHELLP
                       DO 445 N = 1,NGQEXQ
                          INTSCR (N,I,K,L) = INTSCR (N,I,KP1,LM1)
     +                                + CD * INTSCR (N,I,K,LM1)
  445                  CONTINUE
  444               CONTINUE
  443            CONTINUE
             END IF

             DO 446 L = 0,SHELLD
             DO 446 K = 0,SHELLC
             DO 446 I = 0,SHELLP
             DO 446 N = 1,NGQEXQ
                INT2D (N,I,K,L) = INTSCR (N,I,K,L)
  446        CONTINUE

         ELSE

             IF (CD.EQ.ZERO) THEN
                 DO 450 L = 1,SHELLC
                    LM1 = L - 1
                    DO 451 K = 0,SHELLQ-L
                       KP1 = K + 1
                       DO 452 I = 0,SHELLP
                       DO 452 N = 1,NGQEXQ
                          INTSCR (N,I,K,L) = INTSCR (N,I,KP1,LM1)
  452                  CONTINUE
  451               CONTINUE
  450            CONTINUE
             ELSE
                 DO 453 L = 1,SHELLC
                    LM1 = L - 1
                    DO 454 K = 0,SHELLQ-L
                       KP1 = K + 1
                       DO 455 I = 0,SHELLP
                       DO 455 N = 1,NGQEXQ
                          INTSCR (N,I,K,L) = INTSCR (N,I,KP1,LM1)
     +                                - CD * INTSCR (N,I,K,LM1)
  455                  CONTINUE
  454               CONTINUE
  453            CONTINUE
             END IF

             DO 456 L = 0,SHELLD
             DO 456 K = 0,SHELLC
             DO 456 I = 0,SHELLP
             DO 456 N = 1,NGQEXQ
                INT2D (N,I,K,L) = INTSCR (N,I,L,K)
  456        CONTINUE

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
