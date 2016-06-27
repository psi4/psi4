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
         SUBROUTINE  ERD__2D_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      NGQEXQ,
     +                      WTS,
     +                      B00,B01,B10,
     +                      C00X,C00Y,C00Z,
     +                      D00X,D00Y,D00Z,
     +                      CASE2D,
     +
     +                               INT2DX,
     +                               INT2DY,
     +                               INT2DZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__2D_PQ_INTEGRALS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation calculates a full table of 2D PQ X,Y,Z
C                integrals using the Rys vertical recurrence scheme
C                VRR explained below.
C
C                The Rys weight is multiplied to the 2DX PQ integral
C                to reduce overall FLOP count. Note, that the Rys weight
C                factor needs to be introduced only three times for the
C                starting 2DX PQ integrals for the recurrence scheme,
C                namely to the (0,0), (1,0) and (0,1) elements. The
C                weight factor is then automatically propagated
C                through the vertical transfer equations (see below).
C                The recurrence scheme VRR is due to Rys, Dupuis and
C                King, J. Comp. Chem. 4, p.154-157 (1983).
C
C
C                   INT2D (0,0) = 1.D0    (* WEIGHT for the 2DX case)
C                   INT2D (1,0) = C00     (* WEIGHT for the 2DX case)
C                   INT2D (0,1) = D00     (* WEIGHT for the 2DX case)
C
C                   For I = 1,...,SHELLP-1
C                       INT2D (I+1,0) = I * B10 * INT2D (I-1,0)
C                                         + C00 * INT2D (I,0)
C                   For K = 1,...,SHELLQ-1
C                       INT2D (0,K+1) = K * B01 * INT2D (0,K-1)
C                                         + D00 * INT2D (0,K)
C                   For I = 1,...,SHELLP
C                       INT2D (I,1)   = I * B00 * INT2D (I-1,0)
C                                         + D00 * INT2D (I,0)
C                   For K = 2,...,SHELLQ
C                       INT2D (1,K)   = K * B00 * INT2D (0,K-1)
C                                         + C00 * INT2D (0,K)
C                   For K = 2,...,SHELLQ
C                   For I = 2,...,SHELLP
C                       INT2D (I,K)   = (I-1) * B10 * INT2D (I-2,K)
C                                         + K * B00 * INT2D (I-1,K-1)
C                                             + C00 * INT2D (I-1,K)
C
C
C                The 2D PQ integrals are calculated for all roots (info
C                already present in transmitted VRR coefficients!) and
C                for all exponent quadruples simultaneously and placed
C                into a 3-dimensional array.
C
C
C                  Input:
C
C                    SHELLx      =  maximum shell type for electrons
C                                   1 and 2 (x = P,Q)
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    WTS         =  all quadrature weights
C                    B00,B01,B10 =  VRR expansion coefficients
C                                   (cartesian coordinate independent)
C                    C00x,D00x   =  cartesian coordinate dependent
C                                   VRR expansion coefficients
C                                   (x = X,Y,Z)
C                    CASE2D      =  logical flag for simplifications
C                                   in 2D integral evaluation for
C                                   low quantum numbers
C
C
C                  Output:
C
C                    INT2Dx      =  all 2D PQ integrals for each
C                                   cartesian component (x = X,Y,Z)
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT  NONE

         INTEGER   CASE2D
         INTEGER   I,K,N
         INTEGER   I1,I2,K1,K2
         INTEGER   NGQEXQ
         INTEGER   SHELLP,SHELLQ

         DOUBLE PRECISION  B0,B1
         DOUBLE PRECISION  F,F1,F2
         DOUBLE PRECISION  ONE
         DOUBLE PRECISION  WEIGHT

         DOUBLE PRECISION  B00  (1:NGQEXQ)
         DOUBLE PRECISION  B01  (1:NGQEXQ)
         DOUBLE PRECISION  B10  (1:NGQEXQ)
         DOUBLE PRECISION  C00X (1:NGQEXQ)
         DOUBLE PRECISION  C00Y (1:NGQEXQ)
         DOUBLE PRECISION  C00Z (1:NGQEXQ)
         DOUBLE PRECISION  D00X (1:NGQEXQ)
         DOUBLE PRECISION  D00Y (1:NGQEXQ)
         DOUBLE PRECISION  D00Z (1:NGQEXQ)
         DOUBLE PRECISION  WTS  (1:NGQEXQ)

         DOUBLE PRECISION  INT2DX (1:NGQEXQ,0:SHELLP,0:SHELLQ)
         DOUBLE PRECISION  INT2DY (1:NGQEXQ,0:SHELLP,0:SHELLQ)
         DOUBLE PRECISION  INT2DZ (1:NGQEXQ,0:SHELLP,0:SHELLQ)

         PARAMETER  (ONE = 1.D0)
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
C                each leading to simplifications in the VRR formulas.
C                The case present has been evaluated outside this
C                routine and is transmitted via argument.
C
C
         GOTO (1,3,3,2,4,4,2,4,4) CASE2D
C
C
C             ...the case P = s-shell and Q = s-shell.
C
C
    1    DO 100 N = 1,NGQEXQ
            INT2DX (N,0,0) = WTS (N)
            INT2DY (N,0,0) = ONE
            INT2DZ (N,0,0) = ONE
  100    CONTINUE

         RETURN
C
C
C             ...the cases P = s-shell and Q >= p-shell.
C                Evaluate I=0 and K=0,1.
C
C
    2    DO 200 N = 1,NGQEXQ
            WEIGHT = WTS (N)
            INT2DX (N,0,0) = WEIGHT
            INT2DX (N,0,1) = D00X (N) * WEIGHT
            INT2DY (N,0,0) = ONE
            INT2DY (N,0,1) = D00Y (N)
            INT2DZ (N,0,0) = ONE
            INT2DZ (N,0,1) = D00Z (N)
  200    CONTINUE
C
C
C             ...evaluate I=0 and K=2,SHELLQ (if any).
C
C
         F = ONE
         DO 210 K = 2,SHELLQ
            K1 = K - 1
            K2 = K - 2
            DO 212 N = 1,NGQEXQ
               B1 = F * B01 (N)
               INT2DX (N,0,K) =        B1 * INT2DX (N,0,K2)
     +                         + D00X (N) * INT2DX (N,0,K1)
               INT2DY (N,0,K) =        B1 * INT2DY (N,0,K2)
     +                         + D00Y (N) * INT2DY (N,0,K1)
               INT2DZ (N,0,K) =        B1 * INT2DZ (N,0,K2)
     +                         + D00Z (N) * INT2DZ (N,0,K1)
  212       CONTINUE
            F = F + ONE
  210    CONTINUE

         RETURN
C
C
C             ...the cases P >= p-shell and Q = s-shell.
C                Evaluate I=0,1 and K=0.
C
C
    3    DO 300 N = 1,NGQEXQ
            WEIGHT = WTS (N)
            INT2DX (N,0,0) = WEIGHT
            INT2DX (N,1,0) = C00X (N) * WEIGHT
            INT2DY (N,0,0) = ONE
            INT2DY (N,1,0) = C00Y (N)
            INT2DZ (N,0,0) = ONE
            INT2DZ (N,1,0) = C00Z (N)
  300    CONTINUE
C
C
C             ...evaluate I=2,SHELLP (if any) and K=0.
C
C
         F = ONE
         DO 310 I = 2,SHELLP
            I1 = I - 1
            I2 = I - 2
            DO 312 N = 1,NGQEXQ
               B1 = F * B10 (N)
               INT2DX (N,I,0) =        B1 * INT2DX (N,I2,0)
     +                         + C00X (N) * INT2DX (N,I1,0)
               INT2DY (N,I,0) =        B1 * INT2DY (N,I2,0)
     +                         + C00Y (N) * INT2DY (N,I1,0)
               INT2DZ (N,I,0) =        B1 * INT2DZ (N,I2,0)
     +                         + C00Z (N) * INT2DZ (N,I1,0)

  312       CONTINUE
            F = F + ONE
  310    CONTINUE

         RETURN
C
C
C             ...the cases P >= p-shell and Q >= p-shell.
C                Evaluate I=0,SHELLP       I=0
C                         K=0        and   K=0,SHELLQ
C
C
    4    DO 400 N = 1,NGQEXQ
            WEIGHT = WTS (N)
            INT2DX (N,0,0) = WEIGHT
            INT2DX (N,1,0) = C00X (N) * WEIGHT
            INT2DX (N,0,1) = D00X (N) * WEIGHT
            INT2DY (N,0,0) = ONE
            INT2DY (N,1,0) = C00Y (N)
            INT2DY (N,0,1) = D00Y (N)
            INT2DZ (N,0,0) = ONE
            INT2DZ (N,1,0) = C00Z (N)
            INT2DZ (N,0,1) = D00Z (N)
  400    CONTINUE

         F = ONE
         DO 410 I = 2,SHELLP
            I1 = I - 1
            I2 = I - 2
            DO 412 N = 1,NGQEXQ
               B1 = F * B10 (N)
               INT2DX (N,I,0) =        B1 * INT2DX (N,I2,0)
     +                         + C00X (N) * INT2DX (N,I1,0)
               INT2DY (N,I,0) =        B1 * INT2DY (N,I2,0)
     +                         + C00Y (N) * INT2DY (N,I1,0)
               INT2DZ (N,I,0) =        B1 * INT2DZ (N,I2,0)
     +                         + C00Z (N) * INT2DZ (N,I1,0)

  412       CONTINUE
            F = F + ONE
  410    CONTINUE

         F = ONE
         DO 414 K = 2,SHELLQ
            K1 = K - 1
            K2 = K - 2
            DO 416 N = 1,NGQEXQ
               B1 = F * B01 (N)
               INT2DX (N,0,K) =        B1 * INT2DX (N,0,K2)
     +                         + D00X (N) * INT2DX (N,0,K1)
               INT2DY (N,0,K) =        B1 * INT2DY (N,0,K2)
     +                         + D00Y (N) * INT2DY (N,0,K1)
               INT2DZ (N,0,K) =        B1 * INT2DZ (N,0,K2)
     +                         + D00Z (N) * INT2DZ (N,0,K1)
  416       CONTINUE
            F = F + ONE
  414    CONTINUE
C
C
C             ...evaluate I=1,SHELLP and K=1,SHELLQ (if any)
C                in most economical way.
C
C
         IF (SHELLQ.LE.SHELLP) THEN

             F1 = ONE
             DO 420 K = 1,SHELLQ
                K1 = K - 1
                DO 421 N = 1,NGQEXQ
                   B0 = F1 * B00 (N)
                   INT2DX (N,1,K) =        B0 * INT2DX (N,0,K1)
     +                             + C00X (N) * INT2DX (N,0,K)
                   INT2DY (N,1,K) =        B0 * INT2DY (N,0,K1)
     +                             + C00Y (N) * INT2DY (N,0,K)
                   INT2DZ (N,1,K) =        B0 * INT2DZ (N,0,K1)
     +                             + C00Z (N) * INT2DZ (N,0,K)
  421           CONTINUE
                F2 = ONE
                DO 422 I = 2,SHELLP
                   I1 = I - 1
                   I2 = I - 2
                   DO 423 N = 1,NGQEXQ
                      B0 = F1 * B00 (N)
                      B1 = F2 * B10 (N)
                      INT2DX (N,I,K) =        B0 * INT2DX (N,I1,K1)
     +                                      + B1 * INT2DX (N,I2,K)
     +                                + C00X (N) * INT2DX (N,I1,K)
                      INT2DY (N,I,K) =        B0 * INT2DY (N,I1,K1)
     +                                      + B1 * INT2DY (N,I2,K)
     +                                + C00Y (N) * INT2DY (N,I1,K)
                      INT2DZ (N,I,K) =        B0 * INT2DZ (N,I1,K1)
     +                                      + B1 * INT2DZ (N,I2,K)
     +                                + C00Z (N) * INT2DZ (N,I1,K)
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
                   B0 = F1 * B00 (N)
                   INT2DX (N,I,1) =        B0 * INT2DX (N,I1,0)
     +                             + D00X (N) * INT2DX (N,I,0)
                   INT2DY (N,I,1) =        B0 * INT2DY (N,I1,0)
     +                             + D00Y (N) * INT2DY (N,I,0)
                   INT2DZ (N,I,1) =        B0 * INT2DZ (N,I1,0)
     +                             + D00Z (N) * INT2DZ (N,I,0)
  431           CONTINUE
                F2 = ONE
                DO 432 K = 2,SHELLQ
                   K1 = K - 1
                   K2 = K - 2
                   DO 433 N = 1,NGQEXQ
                      B0 = F1 * B00 (N)
                      B1 = F2 * B01 (N)
                      INT2DX (N,I,K) =        B0 * INT2DX (N,I1,K1)
     +                                      + B1 * INT2DX (N,I,K2)
     +                                + D00X (N) * INT2DX (N,I,K1)
                      INT2DY (N,I,K) =        B0 * INT2DY (N,I1,K1)
     +                                      + B1 * INT2DY (N,I,K2)
     +                                + D00Y (N) * INT2DY (N,I,K1)
                      INT2DZ (N,I,K) =        B0 * INT2DZ (N,I1,K1)
     +                                      + B1 * INT2DZ (N,I,K2)
     +                                + D00Z (N) * INT2DZ (N,I,K1)
  433              CONTINUE
                   F2 = F2 + ONE
  432           CONTINUE
                F1 = F1 + ONE
  430        CONTINUE

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
