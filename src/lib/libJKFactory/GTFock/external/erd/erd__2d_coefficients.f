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
         SUBROUTINE  ERD__2D_COEFFICIENTS
     +
     +                    ( MIJ,MKL,MIJKL,
     +                      NGQP,MGQIJKL,
     +                      ATOMAB,ATOMCD,
     +                      P,Q,
     +                      PX,PY,PZ,QX,QY,QZ,
     +                      PAX,PAY,PAZ,QCX,QCY,QCZ,
     +                      PINVHF,QINVHF,PQPINV,
     +                      RTS,
     +                      CASE2D,
     +
     +                               B00,B01,B10,
     +                               C00X,C00Y,C00Z,
     +                               D00X,D00Y,D00Z )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__2D_COEFFICIENTS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation evaluates the VRR coefficients for
C                the 2D integrals for the present set of NGQP roots
C                corresponding to all i,j,k,l exponent quadruplets.
C
C                Not all coefficients are needed in case there are
C                s- or p-shells present on either P or Q side. Also
C                observe that the cartesian distance components PAX,
C                PAY and PAZ are all equal to zero in case the two
C                atomic centers A and B coincide, in which case they
C                need not to be addressed inside the algorithm.
C                Likewise for QCX,QCY and QCZ, if centers C and D
C                coincide.
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
C                    ATOMAB(CD)   =  is true, if centers A and B
C                                    (C and D) coincide
C                    P(Q)         =  current MIJ (MKL) exponent sums
C                                    for csh A and B (C and D)
C                    Px(Qx)       =  current MIJ (MKL) coordinates
C                                    x=X,Y,Z for gaussian product
C                                    centers P=A+B (Q=C+D)
C                    PAx(QCx)     =  current MIJ (MKL) coordinate
C                                    x=X,Y,Z differences P-A (Q-C)
C                                    between centers P and A (Q and C)
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
C                    Bxx          =  the coordinate independent
C                                    B-coefficients (xx=00,01,10)
C                    C00x         =  the C-coefficients (individual
C                                    cartesian components x=X,Y,Z) for
C                                    shell expansion on center P
C                    D00x         =  the D-coefficients (individual
C                                    cartesian components x=X,Y,Z) for
C                                    shell expansion on center Q
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

         LOGICAL   ATOMAB,ATOMCD

         INTEGER   CASE2D
         INTEGER   IJ,KL
         INTEGER   M,N
         INTEGER   MIJ,MKL,MIJKL
         INTEGER   NG,NGQP,MGQIJKL

         DOUBLE PRECISION  PIJ,PXIJ,PYIJ,PZIJ,PAXIJ,PAYIJ,PAZIJ
         DOUBLE PRECISION  PQX,PQY,PQZ
         DOUBLE PRECISION  PSCALE,QSCALE
         DOUBLE PRECISION  QCXKL,QCYKL,QCZKL
         DOUBLE PRECISION  ROOT,PROOT,QROOT
         DOUBLE PRECISION  TWOP,TWOQ,TWOPQ
         DOUBLE PRECISION  HALF,ONE

         DOUBLE PRECISION  P      (1:MIJ)
         DOUBLE PRECISION  PX     (1:MIJ)
         DOUBLE PRECISION  PY     (1:MIJ)
         DOUBLE PRECISION  PZ     (1:MIJ)
         DOUBLE PRECISION  PAX    (1:MIJ)
         DOUBLE PRECISION  PAY    (1:MIJ)
         DOUBLE PRECISION  PAZ    (1:MIJ)
         DOUBLE PRECISION  PINVHF (1:MIJ)

         DOUBLE PRECISION  Q      (1:MKL)
         DOUBLE PRECISION  QX     (1:MKL)
         DOUBLE PRECISION  QY     (1:MKL)
         DOUBLE PRECISION  QZ     (1:MKL)
         DOUBLE PRECISION  QCX    (1:MKL)
         DOUBLE PRECISION  QCY    (1:MKL)
         DOUBLE PRECISION  QCZ    (1:MKL)
         DOUBLE PRECISION  QINVHF (1:MKL)

         DOUBLE PRECISION  PQPINV (1:MIJKL)

         DOUBLE PRECISION  B00    (1:MGQIJKL)
         DOUBLE PRECISION  B01    (1:MGQIJKL)
         DOUBLE PRECISION  B10    (1:MGQIJKL)
         DOUBLE PRECISION  C00X   (1:MGQIJKL)
         DOUBLE PRECISION  C00Y   (1:MGQIJKL)
         DOUBLE PRECISION  C00Z   (1:MGQIJKL)
         DOUBLE PRECISION  D00X   (1:MGQIJKL)
         DOUBLE PRECISION  D00Y   (1:MGQIJKL)
         DOUBLE PRECISION  D00Z   (1:MGQIJKL)
         DOUBLE PRECISION  RTS    (1:MGQIJKL)

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
         GOTO (1,2,5,3,4,7,6,8,9) CASE2D
C
C
C             ...the case P = s-shell and Q = s-shell.
C                (no coefficients here)
C
C
    1    RETURN
C
C
C             ...the case P = p-shell and Q = s-shell.
C                (no B00,B01,B10,D00)
C
C
    2    IF (ATOMAB) THEN
             M = 0
             N = 0
             DO 200 IJ = 1,MIJ
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                DO 210 KL = 1,MKL
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   QSCALE = - Q (KL) * PQPINV (M)

                   DO 215 NG = 1,NGQP
                      N = N + 1
                      QROOT = QSCALE * RTS (N)
                      C00X (N) = QROOT * PQX
                      C00Y (N) = QROOT * PQY
                      C00Z (N) = QROOT * PQZ
  215              CONTINUE

  210           CONTINUE
  200        CONTINUE
         ELSE
             M = 0
             N = 0
             DO 220 IJ = 1,MIJ
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                PAXIJ = PAX (IJ)
                PAYIJ = PAY (IJ)
                PAZIJ = PAZ (IJ)
                DO 230 KL = 1,MKL
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   QSCALE = Q (KL) * PQPINV (M)

                   DO 235 NG = 1,NGQP
                      N = N + 1
                      QROOT = QSCALE * RTS (N)
                      C00X (N) = PAXIJ - QROOT * PQX
                      C00Y (N) = PAYIJ - QROOT * PQY
                      C00Z (N) = PAZIJ - QROOT * PQZ
  235              CONTINUE

  230           CONTINUE
  220        CONTINUE
         END IF

         RETURN
C
C
C             ...the case P = s-shell and Q = p-shell.
C                (no B00,B01,B10,C00)
C
C
    3    IF (ATOMCD) THEN
             M = 0
             N = 0
             DO 300 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                DO 310 KL = 1,MKL
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   PSCALE = PIJ * PQPINV (M)

                   DO 315 NG = 1,NGQP
                      N = N + 1
                      PROOT  = PSCALE * RTS (N)
                      D00X (N) = PROOT * PQX
                      D00Y (N) = PROOT * PQY
                      D00Z (N) = PROOT * PQZ
  315              CONTINUE

  310           CONTINUE
  300        CONTINUE
         ELSE
             M = 0
             N = 0
             DO 320 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                DO 330 KL = 1,MKL
                   QCXKL = QCX (KL)
                   QCYKL = QCY (KL)
                   QCZKL = QCZ (KL)
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   PSCALE = PIJ * PQPINV (M)

                   DO 335 NG = 1,NGQP
                      N = N + 1
                      PROOT  = PSCALE * RTS (N)
                      D00X (N) = QCXKL + PROOT * PQX
                      D00Y (N) = QCYKL + PROOT * PQY
                      D00Z (N) = QCZKL + PROOT * PQZ
  335              CONTINUE

  330           CONTINUE
  320        CONTINUE
         END IF

         RETURN
C
C
C             ...the case P = p-shell and Q = p-shell.
C                (no B01,B10)
C
C
    4    IF (ATOMAB.AND.ATOMCD) THEN

             M = 0
             N = 0
             DO 400 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                DO 410 KL = 1,MKL
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   TWOPQ = HALF * PQPINV (M)
                   PSCALE = PIJ * PQPINV (M)
                   QSCALE = ONE - PSCALE

                   DO 415 NG = 1,NGQP
                      N = N + 1
                      ROOT = RTS (N)
                      PROOT = PSCALE * ROOT
                      QROOT = - QSCALE * ROOT
                      B00  (N) = ROOT * TWOPQ
                      C00X (N) = QROOT * PQX
                      C00Y (N) = QROOT * PQY
                      C00Z (N) = QROOT * PQZ
                      D00X (N) = PROOT * PQX
                      D00Y (N) = PROOT * PQY
                      D00Z (N) = PROOT * PQZ
  415              CONTINUE

  410           CONTINUE
  400        CONTINUE

         ELSE IF (ATOMAB) THEN

             M = 0
             N = 0
             DO 420 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                DO 430 KL = 1,MKL
                   QCXKL = QCX (KL)
                   QCYKL = QCY (KL)
                   QCZKL = QCZ (KL)
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   TWOPQ = HALF * PQPINV (M)
                   PSCALE = PIJ * PQPINV (M)
                   QSCALE = ONE - PSCALE

                   DO 435 NG = 1,NGQP
                      N = N + 1
                      ROOT = RTS (N)
                      PROOT = PSCALE * ROOT
                      QROOT = - QSCALE * ROOT
                      B00  (N) = ROOT * TWOPQ
                      C00X (N) = QROOT * PQX
                      C00Y (N) = QROOT * PQY
                      C00Z (N) = QROOT * PQZ
                      D00X (N) = QCXKL + PROOT * PQX
                      D00Y (N) = QCYKL + PROOT * PQY
                      D00Z (N) = QCZKL + PROOT * PQZ
  435              CONTINUE

  430           CONTINUE
  420        CONTINUE

         ELSE IF (ATOMCD) THEN

             M = 0
             N = 0
             DO 440 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                PAXIJ = PAX (IJ)
                PAYIJ = PAY (IJ)
                PAZIJ = PAZ (IJ)
                DO 450 KL = 1,MKL
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   TWOPQ = HALF * PQPINV (M)
                   PSCALE = PIJ * PQPINV (M)
                   QSCALE = ONE - PSCALE

                   DO 455 NG = 1,NGQP
                      N = N + 1
                      ROOT = RTS (N)
                      PROOT = PSCALE * ROOT
                      QROOT = QSCALE * ROOT
                      B00  (N) = ROOT * TWOPQ
                      C00X (N) = PAXIJ - QROOT * PQX
                      C00Y (N) = PAYIJ - QROOT * PQY
                      C00Z (N) = PAZIJ - QROOT * PQZ
                      D00X (N) = PROOT * PQX
                      D00Y (N) = PROOT * PQY
                      D00Z (N) = PROOT * PQZ
  455              CONTINUE

  450           CONTINUE
  440        CONTINUE
         ELSE
             M = 0
             N = 0
             DO 460 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                PAXIJ = PAX (IJ)
                PAYIJ = PAY (IJ)
                PAZIJ = PAZ (IJ)
                DO 470 KL = 1,MKL
                   QCXKL = QCX (KL)
                   QCYKL = QCY (KL)
                   QCZKL = QCZ (KL)
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   TWOPQ = HALF * PQPINV (M)
                   PSCALE = PIJ * PQPINV (M)
                   QSCALE = ONE - PSCALE

                   DO 475 NG = 1,NGQP
                      N = N + 1
                      ROOT = RTS (N)
                      PROOT = PSCALE * ROOT
                      QROOT = QSCALE * ROOT
                      B00  (N) = ROOT * TWOPQ
                      C00X (N) = PAXIJ - QROOT * PQX
                      C00Y (N) = PAYIJ - QROOT * PQY
                      C00Z (N) = PAZIJ - QROOT * PQZ
                      D00X (N) = QCXKL + PROOT * PQX
                      D00Y (N) = QCYKL + PROOT * PQY
                      D00Z (N) = QCZKL + PROOT * PQZ
  475              CONTINUE

  470           CONTINUE
  460        CONTINUE
         END IF

         RETURN
C
C
C             ...the case P > p-shell and Q = s-shell.
C                (no B00,B01,D00)
C
C
    5    IF (ATOMAB) THEN
             M = 0
             N = 0
             DO 500 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                TWOP = PINVHF (IJ)
                DO 510 KL = 1,MKL
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   QSCALE = ONE - PIJ * PQPINV (M)

                   DO 515 NG = 1,NGQP
                      N = N + 1
                      QROOT  = - QSCALE * RTS (N)
                      B10  (N)= (ONE + QROOT) * TWOP
                      C00X (N) = QROOT * PQX
                      C00Y (N) = QROOT * PQY
                      C00Z (N) = QROOT * PQZ
  515              CONTINUE

  510           CONTINUE
  500        CONTINUE
         ELSE
             M = 0
             N = 0
             DO 520 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                PAXIJ = PAX (IJ)
                PAYIJ = PAY (IJ)
                PAZIJ = PAZ (IJ)
                TWOP = PINVHF (IJ)
                DO 530 KL = 1,MKL
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   QSCALE = ONE - PIJ * PQPINV (M)

                   DO 535 NG = 1,NGQP
                      N = N + 1
                      QROOT  = QSCALE * RTS (N)
                      B10  (N)= (ONE - QROOT) * TWOP
                      C00X (N) = PAXIJ - QROOT * PQX
                      C00Y (N) = PAYIJ - QROOT * PQY
                      C00Z (N) = PAZIJ - QROOT * PQZ
  535              CONTINUE

  530           CONTINUE
  520        CONTINUE
         END IF

         RETURN
C
C
C             ...the case P = s-shell and Q > p-shell.
C                (no B00,B10,C00)
C
C
    6    IF (ATOMCD) THEN
             M = 0
             N = 0
             DO 600 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                DO 610 KL = 1,MKL
                   TWOQ = QINVHF (KL)
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   PSCALE = PIJ * PQPINV (M)

                   DO 615 NG = 1,NGQP
                      N = N + 1
                      PROOT  = PSCALE * RTS (N)
                      B01  (N) = (ONE - PROOT) * TWOQ
                      D00X (N) = PROOT * PQX
                      D00Y (N) = PROOT * PQY
                      D00Z (N) = PROOT * PQZ
  615              CONTINUE

  610           CONTINUE
  600        CONTINUE
         ELSE
             M = 0
             N = 0
             DO 620 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                DO 630 KL = 1,MKL
                   QCXKL = QCX (KL)
                   QCYKL = QCY (KL)
                   QCZKL = QCZ (KL)
                   TWOQ = QINVHF (KL)
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   PSCALE = PIJ * PQPINV (M)

                   DO 635 NG = 1,NGQP
                      N = N + 1
                      PROOT  = PSCALE * RTS (N)
                      B01  (N) = (ONE - PROOT) * TWOQ
                      D00X (N) = QCXKL + PROOT * PQX
                      D00Y (N) = QCYKL + PROOT * PQY
                      D00Z (N) = QCZKL + PROOT * PQZ
  635              CONTINUE

  630           CONTINUE
  620        CONTINUE
         END IF

         RETURN
C
C
C             ...the case P > p-shell and Q = p-shell.
C                (no B01)
C
C
    7    IF (ATOMAB.AND.ATOMCD) THEN

             M = 0
             N = 0
             DO 700 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                TWOP = PINVHF (IJ)
                DO 710 KL = 1,MKL
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   TWOPQ = HALF * PQPINV (M)
                   PSCALE = PIJ * PQPINV (M)
                   QSCALE = ONE - PSCALE

                   DO 715 NG = 1,NGQP
                      N = N + 1
                      ROOT = RTS (N)
                      PROOT = PSCALE * ROOT
                      QROOT = - QSCALE * ROOT
                      B00  (N) = ROOT * TWOPQ
                      B10  (N) = (ONE + QROOT) * TWOP
                      C00X (N) = QROOT * PQX
                      C00Y (N) = QROOT * PQY
                      C00Z (N) = QROOT * PQZ
                      D00X (N) = PROOT * PQX
                      D00Y (N) = PROOT * PQY
                      D00Z (N) = PROOT * PQZ
  715              CONTINUE

  710           CONTINUE
  700        CONTINUE

         ELSE IF (ATOMAB) THEN

             M = 0
             N = 0
             DO 720 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                TWOP = PINVHF (IJ)
                DO 730 KL = 1,MKL
                   QCXKL = QCX (KL)
                   QCYKL = QCY (KL)
                   QCZKL = QCZ (KL)
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   TWOPQ = HALF * PQPINV (M)
                   PSCALE = PIJ * PQPINV (M)
                   QSCALE = ONE - PSCALE

                   DO 735 NG = 1,NGQP
                      N = N + 1
                      ROOT = RTS (N)
                      PROOT = PSCALE * ROOT
                      QROOT = - QSCALE * ROOT
                      B00  (N) = ROOT * TWOPQ
                      B10  (N) = (ONE + QROOT) * TWOP
                      C00X (N) = QROOT * PQX
                      C00Y (N) = QROOT * PQY
                      C00Z (N) = QROOT * PQZ
                      D00X (N) = QCXKL + PROOT * PQX
                      D00Y (N) = QCYKL + PROOT * PQY
                      D00Z (N) = QCZKL + PROOT * PQZ
  735              CONTINUE

  730           CONTINUE
  720        CONTINUE

         ELSE IF (ATOMCD) THEN

             M = 0
             N = 0
             DO 740 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                PAXIJ = PAX (IJ)
                PAYIJ = PAY (IJ)
                PAZIJ = PAZ (IJ)
                TWOP = PINVHF (IJ)
                DO 750 KL = 1,MKL
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   TWOPQ = HALF * PQPINV (M)
                   PSCALE = PIJ * PQPINV (M)
                   QSCALE = ONE - PSCALE

                   DO 755 NG = 1,NGQP
                      N = N + 1
                      ROOT = RTS (N)
                      PROOT = PSCALE * ROOT
                      QROOT = QSCALE * ROOT
                      B00  (N) = ROOT * TWOPQ
                      B10  (N) = (ONE - QROOT) * TWOP
                      C00X (N) = PAXIJ - QROOT * PQX
                      C00Y (N) = PAYIJ - QROOT * PQY
                      C00Z (N) = PAZIJ - QROOT * PQZ
                      D00X (N) = PROOT * PQX
                      D00Y (N) = PROOT * PQY
                      D00Z (N) = PROOT * PQZ
  755              CONTINUE

  750           CONTINUE
  740        CONTINUE
         ELSE
             M = 0
             N = 0
             DO 760 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                PAXIJ = PAX (IJ)
                PAYIJ = PAY (IJ)
                PAZIJ = PAZ (IJ)
                TWOP = PINVHF (IJ)
                DO 770 KL = 1,MKL
                   QCXKL = QCX (KL)
                   QCYKL = QCY (KL)
                   QCZKL = QCZ (KL)
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   TWOPQ = HALF * PQPINV (M)
                   PSCALE = PIJ * PQPINV (M)
                   QSCALE = ONE - PSCALE

                   DO 775 NG = 1,NGQP
                      N = N + 1
                      ROOT = RTS (N)
                      PROOT = PSCALE * ROOT
                      QROOT = QSCALE * ROOT
                      B00  (N) = ROOT * TWOPQ
                      B10  (N) = (ONE - QROOT) * TWOP
                      C00X (N) = PAXIJ - QROOT * PQX
                      C00Y (N) = PAYIJ - QROOT * PQY
                      C00Z (N) = PAZIJ - QROOT * PQZ
                      D00X (N) = QCXKL + PROOT * PQX
                      D00Y (N) = QCYKL + PROOT * PQY
                      D00Z (N) = QCZKL + PROOT * PQZ
  775              CONTINUE

  770           CONTINUE
  760        CONTINUE
         END IF

         RETURN
C
C
C             ...the case P = p-shell and Q > p-shell.
C                (no B10)
C
C
    8    IF (ATOMAB.AND.ATOMCD) THEN

             M = 0
             N = 0
             DO 800 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                DO 810 KL = 1,MKL
                   TWOQ  = QINVHF (KL)
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   TWOPQ = HALF * PQPINV (M)
                   PSCALE = PIJ * PQPINV (M)
                   QSCALE = ONE - PSCALE

                   DO 815 NG = 1,NGQP
                      N = N + 1
                      ROOT = RTS (N)
                      PROOT = PSCALE * ROOT
                      QROOT = - QSCALE * ROOT
                      B00  (N) = ROOT * TWOPQ
                      B01  (N) = (ONE - PROOT) * TWOQ
                      C00X (N) = QROOT * PQX
                      C00Y (N) = QROOT * PQY
                      C00Z (N) = QROOT * PQZ
                      D00X (N) = PROOT * PQX
                      D00Y (N) = PROOT * PQY
                      D00Z (N) = PROOT * PQZ
  815              CONTINUE

  810           CONTINUE
  800        CONTINUE

         ELSE IF (ATOMAB) THEN

             M = 0
             N = 0
             DO 820 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                DO 830 KL = 1,MKL
                   QCXKL = QCX (KL)
                   QCYKL = QCY (KL)
                   QCZKL = QCZ (KL)
                   TWOQ  = QINVHF (KL)
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   TWOPQ = HALF * PQPINV (M)
                   PSCALE = PIJ * PQPINV (M)
                   QSCALE = ONE - PSCALE

                   DO 835 NG = 1,NGQP
                      N = N + 1
                      ROOT = RTS (N)
                      PROOT = PSCALE * ROOT
                      QROOT = - QSCALE * ROOT
                      B00  (N) = ROOT * TWOPQ
                      B01  (N) = (ONE - PROOT) * TWOQ
                      C00X (N) = QROOT * PQX
                      C00Y (N) = QROOT * PQY
                      C00Z (N) = QROOT * PQZ
                      D00X (N) = QCXKL + PROOT * PQX
                      D00Y (N) = QCYKL + PROOT * PQY
                      D00Z (N) = QCZKL + PROOT * PQZ
  835              CONTINUE

  830           CONTINUE
  820        CONTINUE

         ELSE IF (ATOMCD) THEN

             M = 0
             N = 0
             DO 840 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                PAXIJ = PAX (IJ)
                PAYIJ = PAY (IJ)
                PAZIJ = PAZ (IJ)
                DO 850 KL = 1,MKL
                   TWOQ  = QINVHF (KL)
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   TWOPQ = HALF * PQPINV (M)
                   PSCALE = PIJ * PQPINV (M)
                   QSCALE = ONE - PSCALE

                   DO 855 NG = 1,NGQP
                      N = N + 1
                      ROOT = RTS (N)
                      PROOT = PSCALE * ROOT
                      QROOT = QSCALE * ROOT
                      B00  (N) = ROOT * TWOPQ
                      B01  (N) = (ONE - PROOT) * TWOQ
                      C00X (N) = PAXIJ - QROOT * PQX
                      C00Y (N) = PAYIJ - QROOT * PQY
                      C00Z (N) = PAZIJ - QROOT * PQZ
                      D00X (N) = PROOT * PQX
                      D00Y (N) = PROOT * PQY
                      D00Z (N) = PROOT * PQZ
  855              CONTINUE

  850           CONTINUE
  840        CONTINUE
         ELSE
             M = 0
             N = 0
             DO 860 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                PAXIJ = PAX (IJ)
                PAYIJ = PAY (IJ)
                PAZIJ = PAZ (IJ)
                DO 870 KL = 1,MKL
                   QCXKL = QCX (KL)
                   QCYKL = QCY (KL)
                   QCZKL = QCZ (KL)
                   TWOQ  = QINVHF (KL)
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   TWOPQ = HALF * PQPINV (M)
                   PSCALE = PIJ * PQPINV (M)
                   QSCALE = ONE - PSCALE

                   DO 875 NG = 1,NGQP
                      N = N + 1
                      ROOT = RTS (N)
                      PROOT = PSCALE * ROOT
                      QROOT = QSCALE * ROOT
                      B00  (N) = ROOT * TWOPQ
                      B01  (N) = (ONE - PROOT) * TWOQ
                      C00X (N) = PAXIJ - QROOT * PQX
                      C00Y (N) = PAYIJ - QROOT * PQY
                      C00Z (N) = PAZIJ - QROOT * PQZ
                      D00X (N) = QCXKL + PROOT * PQX
                      D00Y (N) = QCYKL + PROOT * PQY
                      D00Z (N) = QCZKL + PROOT * PQZ
  875              CONTINUE

  870           CONTINUE
  860        CONTINUE
         END IF

         RETURN
C
C
C             ...the case P > p-shell and Q > p-shell.
C
C
    9    IF (ATOMAB.AND.ATOMCD) THEN

             M = 0
             N = 0
             DO 900 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                TWOP = PINVHF (IJ)
                DO 910 KL = 1,MKL
                   TWOQ  = QINVHF (KL)
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   TWOPQ = HALF * PQPINV (M)
                   PSCALE = PIJ * PQPINV (M)
                   QSCALE = ONE - PSCALE

                   DO 915 NG = 1,NGQP
                      N = N + 1
                      ROOT = RTS (N)
                      PROOT = PSCALE * ROOT
                      QROOT = - QSCALE * ROOT
                      B00  (N) = ROOT * TWOPQ
                      B01  (N) = (ONE - PROOT) * TWOQ
                      B10  (N) = (ONE + QROOT) * TWOP
                      C00X (N) = QROOT * PQX
                      C00Y (N) = QROOT * PQY
                      C00Z (N) = QROOT * PQZ
                      D00X (N) = PROOT * PQX
                      D00Y (N) = PROOT * PQY
                      D00Z (N) = PROOT * PQZ
  915              CONTINUE

  910           CONTINUE
  900        CONTINUE

         ELSE IF (ATOMAB) THEN

             M = 0
             N = 0
             DO 920 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                TWOP = PINVHF (IJ)
                DO 930 KL = 1,MKL
                   QCXKL = QCX (KL)
                   QCYKL = QCY (KL)
                   QCZKL = QCZ (KL)
                   TWOQ  = QINVHF (KL)
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   TWOPQ = HALF * PQPINV (M)
                   PSCALE = PIJ * PQPINV (M)
                   QSCALE = ONE - PSCALE

                   DO 935 NG = 1,NGQP
                      N = N + 1
                      ROOT = RTS (N)
                      PROOT = PSCALE * ROOT
                      QROOT = - QSCALE * ROOT
                      B00  (N) = ROOT * TWOPQ
                      B01  (N) = (ONE - PROOT) * TWOQ
                      B10  (N) = (ONE + QROOT) * TWOP
                      C00X (N) = QROOT * PQX
                      C00Y (N) = QROOT * PQY
                      C00Z (N) = QROOT * PQZ
                      D00X (N) = QCXKL + PROOT * PQX
                      D00Y (N) = QCYKL + PROOT * PQY
                      D00Z (N) = QCZKL + PROOT * PQZ
  935              CONTINUE

  930           CONTINUE
  920        CONTINUE

         ELSE IF (ATOMCD) THEN

             M = 0
             N = 0
             DO 940 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                PAXIJ = PAX (IJ)
                PAYIJ = PAY (IJ)
                PAZIJ = PAZ (IJ)
                TWOP = PINVHF (IJ)
                DO 950 KL = 1,MKL
                   TWOQ  = QINVHF (KL)
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   TWOPQ = HALF * PQPINV (M)
                   PSCALE = PIJ * PQPINV (M)
                   QSCALE = ONE - PSCALE

                   DO 955 NG = 1,NGQP
                      N = N + 1
                      ROOT = RTS (N)
                      PROOT = PSCALE * ROOT
                      QROOT = QSCALE * ROOT
                      B00  (N) = ROOT * TWOPQ
                      B01  (N) = (ONE - PROOT) * TWOQ
                      B10  (N) = (ONE - QROOT) * TWOP
                      C00X (N) = PAXIJ - QROOT * PQX
                      C00Y (N) = PAYIJ - QROOT * PQY
                      C00Z (N) = PAZIJ - QROOT * PQZ
                      D00X (N) = PROOT * PQX
                      D00Y (N) = PROOT * PQY
                      D00Z (N) = PROOT * PQZ
  955              CONTINUE

  950           CONTINUE
  940        CONTINUE
         ELSE
             M = 0
             N = 0
             DO 960 IJ = 1,MIJ
                PIJ = P (IJ)
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                PAXIJ = PAX (IJ)
                PAYIJ = PAY (IJ)
                PAZIJ = PAZ (IJ)
                TWOP = PINVHF (IJ)
                DO 970 KL = 1,MKL
                   QCXKL = QCX (KL)
                   QCYKL = QCY (KL)
                   QCZKL = QCZ (KL)
                   TWOQ  = QINVHF (KL)
                   PQX = PXIJ - QX (KL)
                   PQY = PYIJ - QY (KL)
                   PQZ = PZIJ - QZ (KL)

                   M = M + 1
                   TWOPQ = HALF * PQPINV (M)
                   PSCALE = PIJ * PQPINV (M)
                   QSCALE = ONE - PSCALE

                   DO 975 NG = 1,NGQP
                      N = N + 1
                      ROOT = RTS (N)
                      PROOT = PSCALE * ROOT
                      QROOT = QSCALE * ROOT
                      B00  (N) = ROOT * TWOPQ
                      B01  (N) = (ONE - PROOT) * TWOQ
                      B10  (N) = (ONE - QROOT) * TWOP
                      C00X (N) = PAXIJ - QROOT * PQX
                      C00Y (N) = PAYIJ - QROOT * PQY
                      C00Z (N) = PAZIJ - QROOT * PQZ
                      D00X (N) = QCXKL + PROOT * PQX
                      D00Y (N) = QCYKL + PROOT * PQY
                      D00Z (N) = QCZKL + PROOT * PQZ
  975              CONTINUE

  970           CONTINUE
  960        CONTINUE
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
