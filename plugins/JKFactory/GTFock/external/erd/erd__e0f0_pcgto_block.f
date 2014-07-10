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
         SUBROUTINE  ERD__E0F0_PCGTO_BLOCK
     +
     +                    ( NBATCH,NINT2D,
     +                      ATOMIC,ATOMAB,ATOMCD,
     +                      MIJ,MKL,MIJKL,
     +                      NIJ,NIJBEG,NIJEND,NKL,NKLBEG,NKLEND,
     +                      NGQP,NMOM,NGQSCR,MGQIJKL,
     +                      NPGTOA,NPGTOB,NPGTOC,NPGTOD,
     +                      NXYZET,NXYZFT,NXYZP,NXYZQ,
     +                      SHELLA,SHELLP,SHELLC,SHELLQ,
     +                      XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     +                      ABX,ABY,ABZ,CDX,CDY,CDZ,
     +                      ALPHAA,ALPHAB,ALPHAC,ALPHAD,
     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                      PRIMA,PRIMB,PRIMC,PRIMD,
     +                      NORMA,NORMB,NORMC,NORMD,
     +                      RHOAB,RHOCD,
     +                      P,PX,PY,PZ,PAX,PAY,PAZ,PINVHF,SCALEP,
     +                      Q,QX,QY,QZ,QCX,QCY,QCZ,QINVHF,SCALEQ,
     +                      RTS,WTS,GQSCR,TVAL,PQPINV,SCALEPQ,
     +                      B00,B01,B10,C00X,C00Y,C00Z,D00X,D00Y,D00Z,
     +                      INT2DX,INT2DY,INT2DZ,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__E0F0_PCGTO_BLOCK
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : ERD__RYS_ROOTS_WEIGHTS
C                ERD__2D_COEFFICIENTS
C                ERD__2D_PQ_INTEGRALS
C                ERD__INT2D_TO_E000
C                ERD__INT2D_TO_E0F0
C                ERD__2D_ATOM_COEFFICIENTS
C                ERD__2D_ATOM_PQ_INTEGRALS
C                ERD__ATOM_INT2D_TO_E000
C                ERD__ATOM_INT2D_TO_E0F0
C  DESCRIPTION : This operation calculates a batch of unnormed electron
C                repulsion integrals between primitive cartesian
C                gaussians for the shell quadruplet range:
C
C                    [E0|F0]     , E = A to P, F = C to Q
C                           ijkl
C
C                and the block of ij and kl exponent pairs. The total
C                number of eris generated here is thus given by the
C                total number of cartesian monomials NXYZET*NXYZFT
C                times the total number of exponent pairs MIJKL in the
C                present block.
C
C                On exit, the batch elements will be stored as:
C
C                             batch (kl,ij,nxyzt)
C
C
C                  Input:
C
C                    NBATCH       =  size of the primitive cartesian
C                                    integral batch
C                    NINT2D       =  space needed for each of the 2D
C                                    X,Y,Z integral arrays
C                    ATOMIC       =  indicates, if purely atomic
C                                    integrals will be evaluated
C                    ATOMAB(CD)   =  indicates, if centers A and B
C                                    (C and D) coincide
C                    MIJ(KL)      =  current # of ij (kl) primitive
C                                    index pairs corresponding to
C                                    the contracted shell pairs A,B
C                                    (C,D)
C                    MIJKL        =  current # of ijkl primitive
C                                    index quadruplets (= MIJ*MKL)
C                    NIJ          =  total # of ij primitive index
C                                    pairs for the contracted shell
C                                    pair A,B
C                    NIJBEG(END)  =  first(last) ij primitive index
C                                    defining the ij block
C                    NKL          =  total # of kl primitive index
C                                    pairs for the contracted shell
C                                    pair C,D
C                    NKLBEG(END)  =  first(last) kl primitive index
C                                    defining the kl block
C                    NGQP         =  # of gaussian quadrature points
C                                    (roots)
C                    NMOM         =  # of necessary moment integrals
C                                    to calculate the quadrature roots
C                    NGQSCR       =  size of gaussian quadrature
C                                    scratch space needed to calculate
C                                    all the quadrature roots
C                    MGQIJKL      =  # of roots times # of ijkl
C                                    quadruplets (= NGQP*MIJKL)
C                    NPGTOx       =  # of primitives per contraction
C                                    for contraction shells x = A,B,C,D
C                    NXYZE(F)T    =  sum of # of cartesian monomials
C                                    for all shells in the range
C                                    E = A,...,P=A+B and in the range
C                                    F = C,...,Q=C+D
C                    NXYZP(Q)     =  # of cartesian monomials for
C                                    the P=A+B and Q=C+D shells
C                    SHELLx       =  the shell type for contraction
C                                    shells x = A,P=A+B,C,Q=C+D
C                    Xx,Yx,Zx     =  the x,y,z-coordinates for centers
C                                    x = A,B,C,D
C                    ABm(CDm)     =  the m=x,y,z-coordinate differences
C                                    between centers A and B (C and D)
C                    ALPHAx       =  the primitive exponents for
C                                    contraction shells x = A,B,C,D
C                    FTABLE       =  Fm (T) table for interpolation
C                                    in low T region
C                    MGRID        =  maximum m in Fm (T) table
C                    NGRID        =  # of T's for which Fm (T) table
C                                    was set up
C                    TMAX         =  maximum T in Fm (T) table
C                    TSTEP        =  difference between two consecutive
C                                    T's in Fm (T) table
C                    TVSTEP       =  Inverse of TSTEP
C                    PRIMx        =  i,j,k,l labels of primitives for
C                                    the respective contraction shells
C                                    x = A,B,C,D
C                    NORMx        =  the normalization factors due to
C                                    the primitive exponents for the
C                                    contraction shells x = A,B,C,D
C                    RHOAB(CD)    =  the complete set of NIJ (NKL)
C                                    exponential prefactors between
C                                    contraction shells A and B
C                                    (C and D)
C                    P            =  will hold current MIJ exponent
C                                    sums for contraction shells A
C                                    and B
C                    Px           =  will hold current MIJ coordinates
C                                    x=X,Y,Z for the gaussian product
C                                    centers P=A+B
C                    PAx          =  will hold current MIJ coordinate
C                                    x=X,Y,Z differences P-A between
C                                    centers P and A
C                    PINVHF       =  will hold current MIJ values of
C                                    1/(2*P), where P are the exponent
C                                    sums for contraction shells A
C                                    and B
C                    SCALEP       =  will hold current MIJ values of
C                                    scaling factors related to point P
C                    Q            =  will hold current MKL exponent
C                                    sums for contraction shells C
C                                    and D
C                    Qx           =  will hold current MKL coordinates
C                                    x=X,Y,Z for the gaussian product
C                                    centers Q=C+D
C                    QCx          =  will hold current MKL coordinate
C                                    x=X,Y,Z differences Q-C between
C                                    centers Q and C
C                    QINVHF       =  will hold current MKL values of
C                                    1/(2*Q), where Q are the exponent
C                                    sums for contraction shells C
C                                    and D
C                    SCALEQ       =  will hold current MKL values of
C                                    scaling factors related to point Q
C                    RTS          =  will hold all current MGQIJKL
C                                    quadrature roots
C                    WTS          =  will hold all current MGQIJKL
C                                    quadrature weights
C                    GQSCR        =  will be used as scratch space
C                                    for determining the quadrature
C                                    roots and weights
C                    TVAL         =  will hold current MIJKL values
C                                    of T-exponents defining the Rys
C                                    weight functions
C                    PQPINV       =  will hold current MIJKL values
C                                    of 1/(P+Q), i.e. the inverses
C                                    of all total exponent sums
C                    SCALEPQ      =  will hold current distinct MIJKL
C                                    (expanded to MGQIJKL) values of
C                                    the overal scaling factors for
C                                    the integrals
C                    Bxx          =  will hold the current MGQIJKL
C                                    coordinate independent VRR
C                                    B-coefficients (xx=00,01,10)
C                    C00x         =  will hold the current MGQIJKL
C                                    VRR C-coefficients (individual
C                                    cartesian components x=X,Y,Z) for
C                                    shell expansion on center P
C                    D00x         =  will hold the current MGQIJKL
C                                    VRR D-coefficients (individual
C                                    cartesian components x=X,Y,Z) for
C                                    shell expansion on center Q
C                    INT2Dx       =  will hold all current 2D integrals
C                                    for each cartesian component
C                                    (x = X,Y,Z)
C
C                  Output:
C
C                    BATCH        =  current batch of primitive
C                                    cartesian [E0|F0] integrals
C
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     ATOMIC,ATOMAB,ATOMCD

         INTEGER     CASE2D,CASEAT
         INTEGER     G000,G010,G020,G030,G040,G050,G060
         INTEGER     I,J,K,L,M,N
         INTEGER     IJ,KL
         INTEGER     MGRID,NGRID
         INTEGER     MIJ,MKL,MIJKL,MGQIJKL
         INTEGER     NBATCH,NINT2D
         INTEGER     NGQP,NMOM,NGQSCR
         INTEGER     NIJ,NKL
         INTEGER     NIJBEG,NIJEND,NKLBEG,NKLEND
         INTEGER     NPGTOA,NPGTOB,NPGTOC,NPGTOD
         INTEGER     NXYZET,NXYZFT,NXYZP,NXYZQ
         INTEGER     SHELLA,SHELLP,SHELLC,SHELLQ

         INTEGER     PRIMA (1:MIJ)
         INTEGER     PRIMB (1:MIJ)
         INTEGER     PRIMC (1:MKL)
         INTEGER     PRIMD (1:MKL)

         DOUBLE PRECISION  ABX,ABY,ABZ,CDX,CDY,CDZ
         DOUBLE PRECISION  EXPA,EXPB,EXPC,EXPD
         DOUBLE PRECISION  INVERS
         DOUBLE PRECISION  PINV,QINV,PVAL,QVAL
         DOUBLE PRECISION  PQPLUS,PQMULT
         DOUBLE PRECISION  PQX,PQY,PQZ
         DOUBLE PRECISION  PSCALE
         DOUBLE PRECISION  PXVAL,PYVAL,PZVAL
         DOUBLE PRECISION  RNPQSQ
         DOUBLE PRECISION  TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
         DOUBLE PRECISION  ZERO,HALF,ONE

         DOUBLE PRECISION  ALPHAA  (1:NPGTOA)
         DOUBLE PRECISION  ALPHAB  (1:NPGTOB)
         DOUBLE PRECISION  ALPHAC  (1:NPGTOC)
         DOUBLE PRECISION  ALPHAD  (1:NPGTOD)

         DOUBLE PRECISION  NORMA   (1:NPGTOA)
         DOUBLE PRECISION  NORMB   (1:NPGTOB)
         DOUBLE PRECISION  NORMC   (1:NPGTOC)
         DOUBLE PRECISION  NORMD   (1:NPGTOD)

         DOUBLE PRECISION  RHOAB   (1:NIJ)
         DOUBLE PRECISION  RHOCD   (1:NKL)

         DOUBLE PRECISION  BATCH   (1:NBATCH)

         DOUBLE PRECISION  P       (1:MIJ)
         DOUBLE PRECISION  PX      (1:MIJ)
         DOUBLE PRECISION  PY      (1:MIJ)
         DOUBLE PRECISION  PZ      (1:MIJ)
         DOUBLE PRECISION  PAX     (1:MIJ)
         DOUBLE PRECISION  PAY     (1:MIJ)
         DOUBLE PRECISION  PAZ     (1:MIJ)
         DOUBLE PRECISION  PINVHF  (1:MIJ)
         DOUBLE PRECISION  SCALEP  (1:MIJ)

         DOUBLE PRECISION  Q       (1:MKL)
         DOUBLE PRECISION  QX      (1:MKL)
         DOUBLE PRECISION  QY      (1:MKL)
         DOUBLE PRECISION  QZ      (1:MKL)
         DOUBLE PRECISION  QCX     (1:MKL)
         DOUBLE PRECISION  QCY     (1:MKL)
         DOUBLE PRECISION  QCZ     (1:MKL)
         DOUBLE PRECISION  QINVHF  (1:MKL)
         DOUBLE PRECISION  SCALEQ  (1:MKL)

         DOUBLE PRECISION  B00     (1:MGQIJKL)
         DOUBLE PRECISION  B01     (1:MGQIJKL)
         DOUBLE PRECISION  B10     (1:MGQIJKL)
         DOUBLE PRECISION  C00X    (1:MGQIJKL)
         DOUBLE PRECISION  C00Y    (1:MGQIJKL)
         DOUBLE PRECISION  C00Z    (1:MGQIJKL)
         DOUBLE PRECISION  D00X    (1:MGQIJKL)
         DOUBLE PRECISION  D00Y    (1:MGQIJKL)
         DOUBLE PRECISION  D00Z    (1:MGQIJKL)
         DOUBLE PRECISION  RTS     (1:MGQIJKL)
         DOUBLE PRECISION  SCALEPQ (1:MGQIJKL)
         DOUBLE PRECISION  WTS     (1:MGQIJKL)
         DOUBLE PRECISION  GQSCR   (1:NGQSCR)
         DOUBLE PRECISION  TVAL    (1:MIJKL)
         DOUBLE PRECISION  PQPINV  (1:MIJKL)

         DOUBLE PRECISION  INT2DX  (1:NINT2D)
         DOUBLE PRECISION  INT2DY  (1:NINT2D)
         DOUBLE PRECISION  INT2DZ  (1:NINT2D)

         DOUBLE PRECISION  FTABLE  (0:MGRID,0:NGRID)

         PARAMETER  (ZERO  =  0.D0)
         PARAMETER  (HALF  =  0.5D0)
         PARAMETER  (ONE   =  1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...predetermine 2D integral case. This is done in
C                order to distinguish the P- and Q-shell combinations
C                for efficient evaluation of the 2D integrals and
C                their VRR coefficients. The cases distinguished
C                are summarized in the following table, indicating
C                the value of CASE2D:
C
C                                Q-shell
C                              s    p   >p
C                            ---------------
C                           |
C                         s |  1    4    7
C                           |
C                P-shell  p |  2    5    8
C                           |
C                        >p |  3    6    9
C                           |
C
C
C
         CASE2D = 3 * MIN0 (2,SHELLQ) + MIN0 (2,SHELLP) + 1
C
C
C             ...predetermine in 'K2' loops the quantities associated
C                with the A,B-part and C,D-part. Set the atom equality
C                case CASEAT here to exploit simplifications due to
C                center equalities later on:
C
C                     CASEAT = 1  -->    atomic (AA|AA) integrals
C                            = 2  -->  2-center (AA|CC) integrals
C                            = 3  -->  3-center (AB|CC) integrals
C                            = 4  -->  3-center (AA|CD) integrals
C                            = 5  -->  4-center (AB|CD) integrals
C
C
         CASEAT = 5

         IF (ATOMAB) THEN
             M = 0
             DO 100 IJ = NIJBEG,NIJEND
                M = M + 1
                I = PRIMA (M)
                J = PRIMB (M)
                PVAL = ALPHAA (I) + ALPHAB (J)
                P (M) = PVAL
                PX (M) = XA
                PY (M) = YA
                PZ (M) = ZA
                PINVHF (M) = HALF / PVAL
                SCALEP (M) = NORMA (I) * NORMB (J)
  100        CONTINUE
             CASEAT = CASEAT - 1
         ELSE
             M = 0
             DO 110 IJ = NIJBEG,NIJEND
                M = M + 1
                I = PRIMA (M)
                J = PRIMB (M)
                EXPA = ALPHAA (I)
                EXPB = ALPHAB (J)
                PVAL = EXPA + EXPB
                PINV = ONE / PVAL
                P (M) = PVAL
                PVAL = - EXPB * PINV
                PAX (M) = PVAL * ABX
                PAY (M) = PVAL * ABY
                PAZ (M) = PVAL * ABZ
                PX (M) = PAX (M) + XA
                PY (M) = PAY (M) + YA
                PZ (M) = PAZ (M) + ZA
                PINVHF (M) = HALF * PINV
                SCALEP (M) = NORMA (I) * NORMB (J) * RHOAB (IJ)
  110        CONTINUE
         END IF

         IF (ATOMCD) THEN
             M = 0
             DO 200 KL = NKLBEG,NKLEND
                M = M + 1
                K = PRIMC (M)
                L = PRIMD (M)
                QVAL = ALPHAC (K) + ALPHAD (L)
                Q (M) = QVAL
                QX (M) = XC
                QY (M) = YC
                QZ (M) = ZC
                QINVHF (M) = HALF / QVAL
                SCALEQ (M) = NORMC (K) * NORMD (L)
  200        CONTINUE
             CASEAT = CASEAT - 2
         ELSE
             M = 0
             DO 220 KL = NKLBEG,NKLEND
                M = M + 1
                K = PRIMC (M)
                L = PRIMD (M)
                EXPC = ALPHAC (K)
                EXPD = ALPHAD (L)
                QVAL = EXPC + EXPD
                QINV = ONE / QVAL
                Q (M) = QVAL
                QVAL = - EXPD * QINV
                QCX (M) = QVAL * CDX
                QCY (M) = QVAL * CDY
                QCZ (M) = QVAL * CDZ
                QX (M) = QCX (M) + XC
                QY (M) = QCY (M) + YC
                QZ (M) = QCZ (M) + ZC
                QINVHF (M) = HALF * QINV
                SCALEQ (M) = NORMC (K) * NORMD (L) * RHOCD (KL)
  220        CONTINUE
         END IF

         IF (ATOMIC) THEN
             CASEAT = CASEAT - 1
         END IF
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block to calculate all T's
C                and scaling factors for the cases:
C
C                     CASEAT = 1  -->    atomic (AA|AA) integrals
C                            = 2  -->  2-center (AA|CC) integrals
C                            = 3  -->  3-center (AB|CC) integrals
C                            = 4  -->  3-center (AA|CD) integrals
C                            = 5  -->  4-center (AB|CD) integrals
C
C                4-center (AB|CD) integrals are checked first
C                (most common occurence in large systems).
C
C
         IF (CASEAT.EQ.5) THEN

             M = 1
             DO 5000 IJ = 1,MIJ
                PVAL = P (IJ)
                PXVAL = PX (IJ)
                PYVAL = PY (IJ)
                PZVAL = PZ (IJ)
                PSCALE = SCALEP (IJ)
                DO 5500 KL = 1,MKL
                   QVAL = Q (KL)
                   PQMULT = PVAL * QVAL
                   PQPLUS = PVAL + QVAL
                   INVERS = ONE / PQPLUS
                   PQX = PXVAL - QX (KL)
                   PQY = PYVAL - QY (KL)
                   PQZ = PZVAL - QZ (KL)

                   TVAL (M) = (PQX * PQX + PQY * PQY + PQZ * PQZ)
     +                                  * PQMULT * INVERS
                   PQPINV (M) = INVERS
                   SCALEPQ (M) = PSCALE * SCALEQ (KL)
     +                                  / (PQMULT * DSQRT (PQPLUS))
                   M = M + 1
 5500           CONTINUE
 5000        CONTINUE

         ELSE IF (CASEAT.EQ.4) THEN

             M = 1
             DO 4000 IJ = 1,MIJ
                PVAL = P (IJ)
                PSCALE = SCALEP (IJ)
                DO 4400 KL = 1,MKL
                   QVAL = Q (KL)
                   PQMULT = PVAL * QVAL
                   PQPLUS = PVAL + QVAL
                   INVERS = ONE / PQPLUS
                   PQX = XA - QX (KL)
                   PQY = YA - QY (KL)
                   PQZ = ZA - QZ (KL)

                   TVAL (M) = (PQX * PQX + PQY * PQY + PQZ * PQZ)
     +                                  * PQMULT * INVERS
                   PQPINV (M) = INVERS
                   SCALEPQ (M) = PSCALE * SCALEQ (KL)
     +                                  / (PQMULT * DSQRT (PQPLUS))
                   M = M + 1
 4400           CONTINUE
 4000        CONTINUE

         ELSE IF (CASEAT.EQ.3) THEN

             M = 1
             DO 3000 IJ = 1,MIJ
                PVAL = P (IJ)
                PXVAL = PX (IJ)
                PYVAL = PY (IJ)
                PZVAL = PZ (IJ)
                PQX = PXVAL - XC
                PQY = PYVAL - YC
                PQZ = PZVAL - ZC
                RNPQSQ = PQX * PQX + PQY * PQY + PQZ * PQZ
                PSCALE = SCALEP (IJ)
                DO 3300 KL = 1,MKL
                   QVAL = Q (KL)
                   PQMULT = PVAL * QVAL
                   PQPLUS = PVAL + QVAL
                   INVERS = ONE / PQPLUS

                   TVAL (M) = RNPQSQ * PQMULT * INVERS
                   PQPINV (M) = INVERS
                   SCALEPQ (M) = PSCALE * SCALEQ (KL)
     +                                  / (PQMULT * DSQRT (PQPLUS))
                   M = M + 1
 3300           CONTINUE
 3000        CONTINUE

         ELSE IF (CASEAT.EQ.2) THEN

             PQX = XA - XC
             PQY = YA - YC
             PQZ = ZA - ZC
             RNPQSQ = PQX * PQX + PQY * PQY + PQZ * PQZ

             M = 1
             DO 2000 IJ = 1,MIJ
                PVAL = P (IJ)
                PSCALE = SCALEP (IJ)
                DO 2200 KL = 1,MKL
                   QVAL = Q (KL)
                   PQMULT = PVAL * QVAL
                   PQPLUS = PVAL + QVAL
                   INVERS = ONE / PQPLUS

                   TVAL (M) = RNPQSQ * PQMULT * INVERS
                   PQPINV (M) = INVERS
                   SCALEPQ (M) = PSCALE * SCALEQ (KL)
     +                                  / (PQMULT * DSQRT (PQPLUS))
                   M = M + 1
 2200           CONTINUE
 2000        CONTINUE

         ELSE

             M = 1
             DO 1000 IJ = 1,MIJ
                PVAL = P (IJ)
                PSCALE = SCALEP (IJ)
                DO 1100 KL = 1,MKL
                   QVAL = Q (KL)
                   PQPLUS = PVAL + QVAL
                   TVAL (M) = ZERO
                   PQPINV (M) = ONE / PQPLUS
                   SCALEPQ (M) = PSCALE * SCALEQ (KL)
     +                           / (PVAL * QVAL * DSQRT (PVAL + QVAL))
                   M = M + 1
 1100           CONTINUE
 1000        CONTINUE

         END IF
C
C
C             ...if necessary, expand the scaling array size from
C                MIJKL to MGQIJKL starting from the last elements.
C
C
         IF (NGQP.GT.1) THEN
             N = MGQIJKL + 1
             DO 10 M = MIJKL,1,-1
                DO 20 I = 1,NGQP
                   SCALEPQ (N-I) = SCALEPQ (M)
   20           CONTINUE
                N = N - NGQP
   10        CONTINUE
         END IF
C
C
C             ...determine memory allocation offsets for the scratch
C                arrays used to calculate the quadrature roots +
C                weights:
C
C                   G000 = offset for A coefficients (Jacobi/Laguerre)
C                   G010 = offset for B coefficients (Jacobi/Laguerre)
C                   G020 = offset for moments (Jacobi/Laguerre)
C                   G030 = offset for diagonals of symmetric termat
C                   G040 = offset for offdiagonals of symmetric termat
C                   G050 = offset for first row intermediates during
C                          evaluation of symmetric termat
C                   G060 = offset for second row intermediates during
C                          evaluation of symmetric termat
C
C
         G000 = 1
         G010 = G000 + NMOM
         G020 = G010 + NMOM - 1
         G030 = G020 + NMOM
         G040 = G030 + NGQP
         G050 = G040 + NGQP
         G060 = G050 + NMOM
C
C
C             ...calculate all roots and weights. Array B00 is passed
C                as a scratch array.
C
C
         CALL    ERD__RYS_ROOTS_WEIGHTS
     +
     +                ( MIJKL,MGQIJKL,
     +                  NGQP,NMOM,
     +                  TVAL,B00,
     +                  FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                  GQSCR(G000),GQSCR(G010),
     +                  GQSCR(G020),
     +                  GQSCR(G030),GQSCR(G040),
     +                  GQSCR(G050),GQSCR(G060),
     +
     +                           RTS,
     +                           WTS )
     +
     +
C
C
C             ...perform the following steps:
C
C                1) generate all VRR coefficients.
C
C                2) construct all 2D PQ x,y,z integrals using all the
C                   weights and all the generated VRR coefficients for
C                   all exponent quadruples.
C
C                3) assemble the complete [E0|F0] batch for all ij and
C                   kl pairs using the 2D integrals. Arrays B00 and B01
C                   are passed as scratch arrays.
C
C                The last step 3) is the most compute intensive and
C                separate routines are provided depending on presence
C                of s-shells. The gainings are in the innermost loops
C                of these routines, which are considerably simplified
C                for the special s-shell cases. Note, that the case
C                in which both P- and Q-shells are s-shells cannot
C                arise, as this case is dealt with in separate routines.
C
C
         IF (.NOT.ATOMIC) THEN

             CALL    ERD__2D_COEFFICIENTS
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
     +                                B00,B01,B10,
     +                                C00X,C00Y,C00Z,
     +                                D00X,D00Y,D00Z )
     +
     +
             CALL    ERD__2D_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MGQIJKL,
     +                      WTS,
     +                      B00,B01,B10,
     +                      C00X,C00Y,C00Z,
     +                      D00X,D00Y,D00Z,
     +                      CASE2D,
     +
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ )
     +
     +
             IF (SHELLQ.EQ.0) THEN

                 CALL    ERD__INT2D_TO_E000
     +
     +                        ( SHELLA,SHELLP,
     +                          NGQP,MIJKL,MGQIJKL,
     +                          NXYZET,NXYZP,
     +                          INT2DX,INT2DY,INT2DZ,
     +                          B00,B01,
     +                          SCALEPQ,
     +
     +                                    BATCH )
     +
     +
             ELSE IF (SHELLP.EQ.0) THEN

                 CALL    ERD__INT2D_TO_E000
     +
     +                        ( SHELLC,SHELLQ,
     +                          NGQP,MIJKL,MGQIJKL,
     +                          NXYZFT,NXYZQ,
     +                          INT2DX,INT2DY,INT2DZ,
     +                          B00,B01,
     +                          SCALEPQ,
     +
     +                                    BATCH )
     +
     +
             ELSE

                 CALL    ERD__INT2D_TO_E0F0
     +
     +                        ( SHELLA,SHELLP,SHELLC,SHELLQ,
     +                          NGQP,MIJKL,MGQIJKL,
     +                          NXYZET,NXYZFT,NXYZP,NXYZQ,
     +                          INT2DX,INT2DY,INT2DZ,
     +                          B00,B01,
     +                          SCALEPQ,
     +
     +                                    BATCH )
     +
     +
             END IF

         ELSE

             CALL    ERD__2D_ATOM_COEFFICIENTS
     +
     +                    ( MIJ,MKL,MIJKL,
     +                      NGQP,MGQIJKL,
     +                      P,
     +                      PINVHF,QINVHF,PQPINV,
     +                      RTS,
     +                      CASE2D,
     +
     +                                B00,B01,B10 )
     +
     +
             CALL    ERD__2D_ATOM_PQ_INTEGRALS
     +
     +                    ( SHELLP,SHELLQ,
     +                      MGQIJKL,
     +                      WTS,
     +                      B00,B01,B10,
     +                      CASE2D,
     +
     +                                INT2DX,
     +                                INT2DY,
     +                                INT2DZ )
     +
     +
             IF (SHELLQ.EQ.0) THEN

                 CALL    ERD__ATOM_INT2D_TO_E000
     +
     +                        ( SHELLA,SHELLP,
     +                          NGQP,MIJKL,MGQIJKL,
     +                          NXYZET,NXYZP,
     +                          INT2DX,INT2DY,INT2DZ,
     +                          B00,B01,
     +                          SCALEPQ,
     +
     +                                     BATCH )
     +
     +
             ELSE IF (SHELLP.EQ.0) THEN

                 CALL    ERD__ATOM_INT2D_TO_E000
     +
     +                        ( SHELLC,SHELLQ,
     +                          NGQP,MIJKL,MGQIJKL,
     +                          NXYZFT,NXYZQ,
     +                          INT2DX,INT2DY,INT2DZ,
     +                          B00,B01,
     +                          SCALEPQ,
     +
     +                                     BATCH )
     +
     +
             ELSE

                 CALL    ERD__ATOM_INT2D_TO_E0F0
     +
     +                        ( SHELLA,SHELLP,SHELLC,SHELLQ,
     +                          NGQP,MIJKL,MGQIJKL,
     +                          NXYZET,NXYZFT,NXYZP,NXYZQ,
     +                          INT2DX,INT2DY,INT2DZ,
     +                          B00,B01,
     +                          SCALEPQ,
     +
     +                                     BATCH )
     +
     +
             END IF

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
