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
         SUBROUTINE  ERD__2D_ATOM_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      NGQEXQ,
     +                      WTS,
     +                      B00,B01,B10,
     +                      CASE2D,
     +
     +                               INT2DX,
     +                               INT2DY,
     +                               INT2DZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__2D_ATOM_PQ_INTEGRALS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation calculates a full table of 2D PQ X,Y,Z
C                atomic integrals using the Rys vertical recurrence
C                scheme VRR explained below. The present atomic case
C                differs from the general nonatomic case in the sense
C                that only B00,B01 and B10 coefficients are required
C                (hence no C00x and D00x with x=X,Y,Z coefficients
C                need to be ever calculated) and only nonzero 2D PQ
C                integrals arise in case the sum of both indices
C                corresponding to P and Q is even.
C
C                The Rys weight is multiplied to the 2DX PQ integral
C                to reduce overall FLOP count. Note, that the Rys weight
C                factor needs to be introduced only two times for
C                the starting atomic 2DX PQ integrals for the recurrence
C                scheme, namely to the (0,0) and (1,1) elements. The
C                weight factor is then automatically propagated
C                through the vertical transfer equations.
C
C                The recurrence scheme VRR is due to Rys, Dupuis and
C                King, J. Comp. Chem. 4, p.154-157 (1983) and details
C                about the scheme can be found in the general nonatomic
C                2D PQ X,Y,Z integral generation routine.
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
C                    B00,B01,B10 =  cartesian coordinate independent
C                                   VRR expansion coefficients
C                    CASE2D      =  logical flag for simplifications
C                                   in 2D integral evaluation for
C                                   low quantum numbers
C
C
C                  Output:
C
C                    INT2Dx      =  all 2D PQ atomic integrals for each
C                                   cartesian component (x = X,Y,Z)
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

         INTEGER   CASE2D
         INTEGER   I,K,N
         INTEGER   I1,I2,K1,K2,KP1
         INTEGER   NGQEXQ
         INTEGER   SHELLP,SHELLQ

         DOUBLE PRECISION  B0,B1
         DOUBLE PRECISION  ONE
         DOUBLE PRECISION  WEIGHT
         DOUBLE PRECISION  XI,XK,XI1,XK1,XKP1

         DOUBLE PRECISION  B00 (1:NGQEXQ)
         DOUBLE PRECISION  B01 (1:NGQEXQ)
         DOUBLE PRECISION  B10 (1:NGQEXQ)
         DOUBLE PRECISION  WTS (1:NGQEXQ)

         DOUBLE PRECISION  INT2DX (1:NGQEXQ,0:SHELLP,0:SHELLQ)
         DOUBLE PRECISION  INT2DY (1:NGQEXQ,0:SHELLP,0:SHELLQ)
         DOUBLE PRECISION  INT2DZ (1:NGQEXQ,0:SHELLP,0:SHELLQ)

         PARAMETER  (ONE   =  1.D0)
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
C                Evaluate I=0 and K=0.
C
C
    2    DO 200 N = 1,NGQEXQ
            INT2DX (N,0,0) = WTS (N)
            INT2DY (N,0,0) = ONE
            INT2DZ (N,0,0) = ONE
  200    CONTINUE
C
C
C             ...evaluate I=0 and even K's in range K=2,SHELLQ
C                (if any).
C
C
         DO 210 K = 2,SHELLQ,2
            K2 = K - 2
            XK1 = DFLOAT (K-1)

            DO 212 N = 1,NGQEXQ
               B1 = B01 (N) * XK1
               INT2DX (N,0,K) = B1 * INT2DX (N,0,K2)
               INT2DY (N,0,K) = B1 * INT2DY (N,0,K2)
               INT2DZ (N,0,K) = B1 * INT2DZ (N,0,K2)
  212       CONTINUE
  210    CONTINUE

         RETURN
C
C
C             ...the cases P >= p-shell and Q = s-shell.
C                Evaluate I=0 and K=0.
C
C
    3    DO 300 N = 1,NGQEXQ
            INT2DX (N,0,0) = WTS (N)
            INT2DY (N,0,0) = ONE
            INT2DZ (N,0,0) = ONE
  300    CONTINUE
C
C
C             ...evaluate even I's in range I=2,SHELLP (if any)
C                and K=0.
C
C
         DO 310 I = 2,SHELLP,2
            I2 = I - 2
            XI1 = DFLOAT (I-1)

            DO 312 N = 1,NGQEXQ
               B1 = B10 (N) * XI1
               INT2DX (N,I,0) = B1 * INT2DX (N,I2,0)
               INT2DY (N,I,0) = B1 * INT2DY (N,I2,0)
               INT2DZ (N,I,0) = B1 * INT2DZ (N,I2,0)
  312       CONTINUE
  310    CONTINUE

         RETURN
C
C
C             ...the cases P >= p-shell and Q >= p-shell.
C                Evaluate I=0,K=0 and I=1,K=1 pairs.
C
C
    4    DO 400 N = 1,NGQEXQ
            B0 = B00  (N)
            WEIGHT = WTS (N)
            INT2DX (N,0,0) = WEIGHT
            INT2DX (N,1,1) = WEIGHT * B0
            INT2DY (N,0,0) = ONE
            INT2DY (N,1,1) = B0
            INT2DZ (N,0,0) = ONE
            INT2DZ (N,1,1) = B0
  400    CONTINUE
C
C
C             ...evaluate: i) even I's in range I=2,SHELLP (if any)
C                             and K=0.
C                         ii) odd I's in range I=3,SHELLP (if any)
C                             and K=1.
C
C                Note, that step i) must be evaluated before step ii)
C                since the even I integrals are needed to evaluate
C                the odd I ones.
C
C
         DO 410 I = 2,SHELLP,2
            I2 = I - 2
            XI1 = DFLOAT (I-1)
            DO 412 N = 1,NGQEXQ
               B1 = B10 (N) * XI1
               INT2DX (N,I,0) = B1 * INT2DX (N,I2,0)
               INT2DY (N,I,0) = B1 * INT2DY (N,I2,0)
               INT2DZ (N,I,0) = B1 * INT2DZ (N,I2,0)
  412       CONTINUE
  410    CONTINUE

         DO 414 I = 3,SHELLP,2
            I1 = I - 1
            XI = DFLOAT (I)
            DO 416 N = 1,NGQEXQ
               B0 = B00 (N) * XI
               INT2DX (N,I,1) = B0 * INT2DX (N,I1,0)
               INT2DY (N,I,1) = B0 * INT2DY (N,I1,0)
               INT2DZ (N,I,1) = B0 * INT2DZ (N,I1,0)
  416       CONTINUE
  414    CONTINUE
C
C
C             ...evaluate I=0,SHELLP and K=2,SHELLQ (if any),
C                matching even I's with even K's and odd I's
C                with odd K's. Note that the even K integrals
C                are needed to evaluate the odd K ones, hence
C                the K loop needs to consider both even K and
C                odd K+1 cases consecutively.
C
C
         DO 420 K = 2,SHELLQ,2
            K1 = K - 1
            K2 = K - 2
            XK = DFLOAT (K)
            XK1 = DFLOAT (K1)

            DO 422 N = 1,NGQEXQ
               B1 = B01 (N) * XK1
               INT2DX (N,0,K) = B1 * INT2DX (N,0,K2)
               INT2DY (N,0,K) = B1 * INT2DY (N,0,K2)
               INT2DZ (N,0,K) = B1 * INT2DZ (N,0,K2)
  422       CONTINUE

            DO 424 I = 2,SHELLP,2
               I1 = I - 1
               I2 = I - 2
               XI1 = DFLOAT (I1)

               DO 426 N = 1,NGQEXQ
                  B0 = B00 (N) * XK
                  B1 = B10 (N) * XI1
                  INT2DX (N,I,K) =   B1 * INT2DX (N,I2,K)
     +                             + B0 * INT2DX (N,I1,K1)
                  INT2DY (N,I,K) =   B1 * INT2DY (N,I2,K)
     +                             + B0 * INT2DY (N,I1,K1)
                  INT2DZ (N,I,K) =   B1 * INT2DZ (N,I2,K)
     +                             + B0 * INT2DZ (N,I1,K1)
  426          CONTINUE
  424       CONTINUE

            IF (K.EQ.SHELLQ) GOTO 420

            KP1 = K + 1
            XKP1 = DFLOAT (KP1)

            DO 421 N = 1,NGQEXQ
               B0 = B00 (N) * XKP1
               INT2DX (N,1,KP1) = B0 * INT2DX (N,0,K)
               INT2DY (N,1,KP1) = B0 * INT2DY (N,0,K)
               INT2DZ (N,1,KP1) = B0 * INT2DZ (N,0,K)
  421       CONTINUE

            DO 423 I = 3,SHELLP,2
               I1 = I - 1
               I2 = I - 2
               XI1 = DFLOAT (I1)

               DO 425 N = 1,NGQEXQ
                  B0 = B00 (N) * XKP1
                  B1 = B10 (N) * XI1
                  INT2DX (N,I,KP1) =   B1 * INT2DX (N,I2,KP1)
     +                               + B0 * INT2DX (N,I1,K)
                  INT2DY (N,I,KP1) =   B1 * INT2DY (N,I2,KP1)
     +                               + B0 * INT2DY (N,I1,K)
                  INT2DZ (N,I,KP1) =   B1 * INT2DZ (N,I2,KP1)
     +                               + B0 * INT2DZ (N,I1,K)
  425          CONTINUE
  423       CONTINUE

  420    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
