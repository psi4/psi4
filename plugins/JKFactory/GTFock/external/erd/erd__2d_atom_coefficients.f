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
         SUBROUTINE  ERD__2D_ATOM_COEFFICIENTS
     +
     +                    ( MIJ,MKL,MIJKL,
     +                      NGQP,MGQIJKL,
     +                      P,
     +                      PINVHF,QINVHF,PQPINV,
     +                      RTS,
     +                      CASE2D,
     +
     +                               B00,B01,B10 )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__2D_ATOM_COEFFICIENTS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation evaluates the atomic VRR coefficients
C                for the 2D integrals for the present set of NGQP roots
C                corresponding to all i,j,k,l exponent quadruplets.
C
C
C                  Input:
C
C                    MIJ(KL)      =  current # of ij (kl) primitive
C                                    index pairs corresponding to
C                                    the csh pairs A,B (C,D)
C                    MIJKL        =  current # of ijkl primitive
C                                    index quadruplets (= MIJ*MKL)
C                    NGQP         =  # of gaussian quadrature points
C                                    (roots)
C                    MGQIJKL      =  # of roots times # of ijkl
C                                    quadruplets (= NGQP*MIJKL)
C                    P            =  current MIJ exponent sums for
C                                    csh A and B
C                    P(Q)INVHF    =  current MIJ (MKL) values of
C                                    1/(2*P(Q)), where P and Q are
C                                    the corresponding exponent sums
C                                    for csh A and B (C and D)
C                    PQPINV       =  current MIJKL values of 1/(P+Q)
C                    RTS          =  current MGQIJKL values of all
C                                    quadrature roots
C                    CASE2D       =  integer value within the range
C                                    from 1 to 9, indicating which
C                                    2D coefficient evaluation case
C                                    is present to trigger specific
C                                    simplified sections of the code
C
C                  Output:
C
C                    Bxx          =  the coordinate independent atomic
C                                    B-coefficients (xx=00,01,10)
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
         INTEGER   IJ,KL
         INTEGER   M,N
         INTEGER   MIJ,MKL,MIJKL
         INTEGER   NG,NGQP,MGQIJKL

         DOUBLE PRECISION  PIJ
         DOUBLE PRECISION  PSCALE,QSCALE
         DOUBLE PRECISION  ROOT
         DOUBLE PRECISION  TWOP,TWOQ,TWOPQ
         DOUBLE PRECISION  HALF,ONE

         DOUBLE PRECISION  P      (1:MIJ)
         DOUBLE PRECISION  PINVHF (1:MIJ)
         DOUBLE PRECISION  QINVHF (1:MKL)

         DOUBLE PRECISION  PQPINV (1:MIJKL)

         DOUBLE PRECISION  B00  (1:MGQIJKL)
         DOUBLE PRECISION  B01  (1:MGQIJKL)
         DOUBLE PRECISION  B10  (1:MGQIJKL)
         DOUBLE PRECISION  RTS  (1:MGQIJKL)

         PARAMETER  (HALF  =  0.5D0)
         PARAMETER  (ONE   =  1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...jump according to the 9 different cases that can arise:
C
C                  P-shell = s- ,p- or higher angular momentum
C                  Q-shell = s- ,p- or higher angular momentum
C
C                each leading to simplifications in the VRR formulas.
C                The case present has been evaluated outside this
C                routine and is transmitted via argument.
C
C
         GOTO (1,1,5,1,4,7,6,8,9) CASE2D
C
C
C             ...the cases: P = s-shell and Q = s-shell
C                           P = p-shell and Q = s-shell
C                           P = s-shell and Q = p-shell
C                (no coefficients here)
C
C
    1    RETURN
C
C
C             ...the case P = p-shell and Q = p-shell.
C                (here we know that NGQP = 2)
C
C
    4    N = 1
         DO 400 M = 1,MIJKL
            TWOPQ = HALF * PQPINV (M)
            B00 (N)   = RTS (N)   * TWOPQ
            B00 (N+1) = RTS (N+1) * TWOPQ
            N = N + 2
  400    CONTINUE

         RETURN
C
C
C             ...the case P > p-shell and Q = s-shell.
C
C
    5    M = 0
         N = 0
         DO 500 IJ = 1,MIJ
            PIJ = P (IJ)
            TWOP = PINVHF (IJ)
            DO 510 KL = 1,MKL
               M = M + 1
               QSCALE = ONE - PIJ * PQPINV (M)

               DO 520 NG = 1,NGQP
                  N = N + 1
                  B10 (N)= (ONE - QSCALE * RTS (N)) * TWOP
  520          CONTINUE

  510       CONTINUE
  500    CONTINUE

         RETURN
C
C
C             ...the case P = s-shell and Q > p-shell.
C
C
    6    M = 0
         N = 0
         DO 600 IJ = 1,MIJ
            PIJ = P (IJ)
            DO 610 KL = 1,MKL
               TWOQ = QINVHF (KL)
               M = M + 1
               PSCALE = PIJ * PQPINV (M)

               DO 620 NG = 1,NGQP
                  N = N + 1
                  B01 (N) = (ONE - PSCALE * RTS (N)) * TWOQ
  620          CONTINUE

  610       CONTINUE
  600    CONTINUE

         RETURN
C
C
C             ...the case P > p-shell and Q = p-shell.
C
C
    7    M = 0
         N = 0
         DO 700 IJ = 1,MIJ
            PIJ = P (IJ)
            TWOP = PINVHF (IJ)
            DO 710 KL = 1,MKL
               M = M + 1
               TWOPQ = HALF * PQPINV (M)
               QSCALE = ONE - PIJ * PQPINV (M)

               DO 720 NG = 1,NGQP
                  N = N + 1
                  ROOT = RTS (N)
                  B00 (N) = ROOT * TWOPQ
                  B10 (N) = (ONE - QSCALE * ROOT) * TWOP
  720          CONTINUE

  710       CONTINUE
  700    CONTINUE

         RETURN
C
C
C             ...the case P = p-shell and Q > p-shell.
C
C
    8    M = 0
         N = 0
         DO 800 IJ = 1,MIJ
            PIJ = P (IJ)
            DO 810 KL = 1,MKL
               TWOQ  = QINVHF (KL)
               M = M + 1
               TWOPQ = HALF * PQPINV (M)
               PSCALE = PIJ * PQPINV (M)

               DO 820 NG = 1,NGQP
                  N = N + 1
                  ROOT = RTS (N)
                  B00 (N) = ROOT * TWOPQ
                  B01 (N) = (ONE - PSCALE * ROOT) * TWOQ
  820          CONTINUE

  810       CONTINUE
  800    CONTINUE

         RETURN
C
C
C             ...the case P > p-shell and Q > p-shell.
C
C
    9    M = 0
         N = 0
         DO 900 IJ = 1,MIJ
            PIJ = P (IJ)
            TWOP = PINVHF (IJ)
            DO 910 KL = 1,MKL
               TWOQ  = QINVHF (KL)
               M = M + 1
               TWOPQ = HALF * PQPINV (M)
               PSCALE = PIJ * PQPINV (M)
               QSCALE = ONE - PSCALE

               DO 920 NG = 1,NGQP
                  N = N + 1
                  ROOT = RTS (N)
                  B00 (N) = ROOT * TWOPQ
                  B01 (N) = (ONE - PSCALE * ROOT) * TWOQ
                  B10 (N) = (ONE - QSCALE * ROOT) * TWOP
  920          CONTINUE

  910       CONTINUE
  900    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
