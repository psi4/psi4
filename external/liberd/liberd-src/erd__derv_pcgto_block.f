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
         SUBROUTINE  ERD__DERV_PCGTO_BLOCK
     +
     +                    ( NBATCH,
     +                      NINT2DX,NINT2DY,NINT2DZ,
     +                      ATOMAB,ATOMCD,
     +                      MIJ,MKL,MIJKL,
     +                      NIJ,NIJBEG,NIJEND,NKL,NKLBEG,NKLEND,
     +                      NGQP,NMOM,NGQSCR,MGQIJKL,
     +                      NPGTOA,NPGTOB,NPGTOC,NPGTOD,
     +                      NZSHELL,NZNXYZ,
     +                      NXYZET,NXYZP,
     +                      NXYZA,NXYZB,NXYZC,NXYZD,
     +                      SHELLA,SHELLB,SHELLC,SHELLD,
     +                      SHELLP,SHELLQ,SHELLT,
     +                      XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     +                      ABX,ABY,ABZ,CDX,CDY,CDZ,
     +                      NDERX,NDERY,NDERZ,
     +                      DERAX,DERAY,DERAZ,
     +                      DERBX,DERBY,DERBZ,
     +                      DERCX,DERCY,DERCZ,
     +                      DERDX,DERDY,DERDZ,
     +                      DERPX,DERPY,DERPZ,
     +                      DERQX,DERQY,DERQZ,
     +                      DIFFX,DIFFY,DIFFZ,
     +                      DIFFA,DIFFB,DIFFC,DIFFD,
     +                      PRIMTYP,ANGMTYP,
     +                      ALPHAA,ALPHAB,ALPHAC,ALPHAD,
     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                      CENEQS,
     +                      CENSQX,CENSQY,CENSQZ,
     +                      PRIMA,PRIMB,PRIMC,PRIMD,
     +                      NORMA,NORMB,NORMC,NORMD,
     +                      RHOAB,RHOCD,
     +                      P,PX,PY,PZ,PAX,PAY,PAZ,PINVHF,SCALEP,
     +                      Q,QX,QY,QZ,QCX,QCY,QCZ,QINVHF,SCALEQ,
     +                      RTS,WTS,GQSCR,TVAL,PQPINV,SCALEPQ,
     +                      B00,B01,B10,C00X,C00Y,C00Z,D00X,D00Y,D00Z,
     +                      EXP2A,EXP2B,EXP2C,EXP2D,
     +                      INT2DX,INT2DY,INT2DZ,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__DERV_PCGTO_BLOCK
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : ERD__RYS_ROOTS_WEIGHTS
C                ERD__2D_COEFFICIENTS
C                ERD__2D_PCD_INTEGRALS
C                ERD__2D_DERV_PCD_INTEGRALS
C                ERD__2D_ABCD_INTEGRALS
C                ERD__2D_DERV_ABCD_INTEGRALS
C                ERD__DERV_INT2D_TO_0000
C                ERD__DERV_INT2D_TO_E0C0
C                ERD__DERV_INT2D_TO_E0CD
C                ERD__DERV_INT2D_TO_A000
C                ERD__DERV_INT2D_TO_AB00
C                ERD__DERV_INT2D_TO_ABC0
C                ERD__DERV_INT2D_TO_ABCD
C  DESCRIPTION : This operation calculates a batch of derivated
C                unnormed electron repulsion integrals between
C                primitive cartesian gaussians.
C
C                The shell quadruplet ranges depend on where the
C                derivation is being performed:
C
C                   i) Derivation(s) on centers C and/or D:
C
C                        [E0|CD]     , E = A to P
C                               ijkl
C
C                  ii) Derivation(s) on centers (A and/or B)
C                      and on centers (C and/or D):
C
C                        [AB|CD]
C                               ijkl
C
C                The index ijkl stands for the block of ij and kl
C                exponent pairs used. The total number of eris
C                generated here is thus also determined by the
C                center(s) of derivation(s):
C
C                   Case i) # of eris: NXYZT = NXYZET*NXYZC*NXYYD
C                  Case ii) # of eris: NXYZT = NXYZA*NXYZB*NXYZC*NXYYD
C
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
C                    NBATCH       =  size of the array that will hold
C                                    the final primitive cartesian
C                                    derivative integral batch as well
C                                    as intermediate differentiated
C                                    2D integrals
C                    NINT2Dx      =  space needed for each of the 2D
C                                    x = X,Y,Z integral arrays (they
C                                    might be different due to different
C                                    orders of differentiation for
C                                    each cartesian component)
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
C                    NXYZET       =  sum of # of cartesian monomials
C                                    for all shells in the range
C                                    E = A,...,P=A+B
C                    NXYZP        =  # of cartesian monomials for
C                                    the P=A+B shell
C                    NXYZx        =  # of cartesian monomials for
C                                    each contraction shell x = A,B,C,D
C                    SHELLx       =  the shell type for each contraction
C                                    shell x = A,B,C,D
C                    SHELLx       =  the shell type summations for
C                                    shells x = P=A+B,Q=C+D,T=A+B+C+D
C                    Xx,Yx,Zx     =  the x,y,z-coordinates for centers
C                                    x = A,B,C,D
C                    ABm(CDm)     =  the m=x,y,z-coordinate differences
C                                    between centers A and B (C and D)
C                    NDERp        =  the total order of differentiation
C                                    with respect to the p = x,y,z
C                                    coordinates
C                    DERxp        =  the order of differentiation on
C                                    centers x = A,B,C,D,P=A+B,Q=C+D
C                                    with respect to the p = x,y,z
C                                    coordinates
C                    DIFFp        =  is true, if differentiation will be
C                                    performed on centers p = A,B,C,D
C                                    involving the p = x,y,z coordinates
C                    PRIMTYP      =  character variable, indicating
C                                    which type of primitives will
C                                    be generated and contracted. Can
C                                    be only 'E0CD' or 'ABCD'
C                    ANGMTYP      =  character variable, indicating
C                                    the overall angular momentum type
C                                    combination without observing
C                                    the order. Can be only 'SSSS',
C                                    'SSSX','SSXX','SXXX' or 'XXXX',
C                                    where the S-symbol indicates the
C                                    presence of an s-shell and the
C                                    X-symbol the presence of a shell
C                                    >= p-shell
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
C                    CENEQS (I,J) =  center equality indicator of size
C                                    4 x 4. The matrix is defined as
C                                    follows: if the centers indexed by
C                                    I and J are equal => value = 1,
C                                    if not => value = 0. The indices
C                                    correspond to the A,B,C,D ordering,
C                                    i.e. 1st index -> A, 2nd -> B, etc
C                    CENSQp       =  center differentiation sequence
C                                    array for the p = x,y,z coordinates
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
C                                    the derivative integrals
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
C                    EXP2x        =  will hold current double primitive
C                                    exponent values in MIJKL order
C                                    (expanded to MGQIJKL) for each
C                                    contracted shell x = A,B,C,D
C                    INT2Dx       =  will hold all current derivated
C                                    2D integrals for each cartesian
C                                    component (x = X,Y,Z) during all
C                                    stages of differentiation
C
C                  Output:
C
C                    BATCH        =  current batch of primitive
C                                    derivative cartesian [E0|CD] or
C                                    [AB|CD] integrals
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT      NONE

         LOGICAL       ATOMAB,ATOMCD
         LOGICAL       DIFFA,DIFFB,DIFFC,DIFFD
         LOGICAL       DIFFX,DIFFY,DIFFZ

         CHARACTER*4   ANGMTYP
         CHARACTER*4   PRIMTYP

         INTEGER       CASE2D,CASEAT
         INTEGER       CENTER
         INTEGER       DER1,DER2,DER3,DER4
         INTEGER       DERAX,DERBX,DERCX,DERDX
         INTEGER       DERAY,DERBY,DERCY,DERDY
         INTEGER       DERAZ,DERBZ,DERCZ,DERDZ
         INTEGER       DEROP1,DEROP2,DEROP3,DEROP4
         INTEGER       DERPX,DERPY,DERPZ
         INTEGER       DERQX,DERQY,DERQZ
         INTEGER       G000,G010,G020,G030,G040,G050,G060
         INTEGER       I,J,K,L,M,N
         INTEGER       IJ,KL
         INTEGER       MGRID,NGRID
         INTEGER       MIJ,MKL,MIJKL,MGQKL,MGQIJKL
         INTEGER       NBATCH
         INTEGER       NDERX,NDERY,NDERZ
         INTEGER       NGQP,NMOM,NGQSCR
         INTEGER       NIJ,NKL
         INTEGER       NIJBEG,NIJEND,NKLBEG,NKLEND
         INTEGER       NINT,NINT2DX,NINT2DY,NINT2DZ
         INTEGER       NPGTOA,NPGTOB,NPGTOC,NPGTOD
         INTEGER       NXYZET,NXYZP
         INTEGER       NXYZA,NXYZB,NXYZC,NXYZD
         INTEGER       SHELLA,SHELLB,SHELLC,SHELLD
         INTEGER       SHELLP,SHELLQ,SHELLT

         INTEGER       CENSQX (0:NDERX)
         INTEGER       CENSQY (0:NDERY)
         INTEGER       CENSQZ (0:NDERZ)

         INTEGER       NZNXYZ  (1:4)
         INTEGER       NZSHELL (1:4)

         INTEGER       PRIMA  (1:MIJ)
         INTEGER       PRIMB  (1:MIJ)
         INTEGER       PRIMC  (1:MKL)
         INTEGER       PRIMD  (1:MKL)

         INTEGER       CENEQS (1:4,1:4)

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
         DOUBLE PRECISION  EXP2A   (1:MGQIJKL)
         DOUBLE PRECISION  EXP2B   (1:MGQIJKL)
         DOUBLE PRECISION  EXP2C   (1:MGQIJKL)
         DOUBLE PRECISION  EXP2D   (1:MGQIJKL)
         DOUBLE PRECISION  RTS     (1:MGQIJKL)
         DOUBLE PRECISION  SCALEPQ (1:MGQIJKL)
         DOUBLE PRECISION  WTS     (1:MGQIJKL)
         DOUBLE PRECISION  GQSCR   (1:NGQSCR)
         DOUBLE PRECISION  TVAL    (1:MIJKL)
         DOUBLE PRECISION  PQPINV  (1:MIJKL)

         DOUBLE PRECISION  INT2DX  (1:NINT2DX)
         DOUBLE PRECISION  INT2DY  (1:NINT2DY)
         DOUBLE PRECISION  INT2DZ  (1:NINT2DZ)

         DOUBLE PRECISION  FTABLE  (0:MGRID,0:NGRID)

         PARAMETER  (ZERO  =  0.D0)
         PARAMETER  (HALF  =  0.5D0)
         PARAMETER  (ONE   =  1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...predetermine in 'K2' loops the quantities associated
C                with the A,B-part and C,D-part. Set the atom equality
C                case CASEAT here to exploit simplifications due to
C                center equalities later on:
C
C                CASEAT = 1  -->  2-center (AA|CC) derivative integrals
C                       = 2  -->  3-center (AB|CC) derivative integrals
C                       = 3  -->  3-center (AA|CD) derivative integrals
C                       = 4  -->  4-center (AB|CD) derivative integrals
C
C
         CASEAT = 4

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
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block to calculate all T's
C                and scaling factors for the cases:
C
C                CASEAT = 1  -->  2-center (AA|CC) derivative integrals
C                       = 2  -->  3-center (AB|CC) derivative integrals
C                       = 3  -->  3-center (AA|CD) derivative integrals
C                       = 4  -->  4-center (AB|CD) derivative integrals
C
C                4-center (AB|CD) integrals are checked first
C                (most common occurence in large systems).
C
C
         IF (CASEAT.EQ.4) THEN

             M = 1
             DO 4000 IJ = 1,MIJ
                PVAL = P (IJ)
                PXVAL = PX (IJ)
                PYVAL = PY (IJ)
                PZVAL = PZ (IJ)
                PSCALE = SCALEP (IJ)
                DO 4400 KL = 1,MKL
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
 4400           CONTINUE
 4000        CONTINUE

         ELSE IF (CASEAT.EQ.3) THEN

             M = 1
             DO 3000 IJ = 1,MIJ
                PVAL = P (IJ)
                PSCALE = SCALEP (IJ)
                DO 3300 KL = 1,MKL
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
 3300           CONTINUE
 3000        CONTINUE

         ELSE IF (CASEAT.EQ.2) THEN

             M = 1
             DO 2000 IJ = 1,MIJ
                PVAL = P (IJ)
                PXVAL = PX (IJ)
                PYVAL = PY (IJ)
                PZVAL = PZ (IJ)
                PQX = PXVAL - XC
                PQY = PYVAL - YC
                PQZ = PZVAL - ZC
                RNPQSQ = PQX * PQX + PQY * PQY + PQZ * PQZ
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

             PQX = XA - XC
             PQY = YA - YC
             PQZ = ZA - ZC
             RNPQSQ = PQX * PQX + PQY * PQY + PQZ * PQZ

             M = 1
             DO 1000 IJ = 1,MIJ
                PVAL = P (IJ)
                PSCALE = SCALEP (IJ)
                DO 1100 KL = 1,MKL
                   QVAL = Q (KL)
                   PQMULT = PVAL * QVAL
                   PQPLUS = PVAL + QVAL
                   INVERS = ONE / PQPLUS

                   TVAL (M) = RNPQSQ * PQMULT * INVERS
                   PQPINV (M) = INVERS
                   SCALEPQ (M) = PSCALE * SCALEQ (KL)
     +                                  / (PQMULT * DSQRT (PQPLUS))
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
C             ...calculate all needed double exponent values for
C                differentiation. First the double exponents on
C                centers A and/or B (if needed).
C
C
         MGQKL = NGQP * MKL

         IF (DIFFA) THEN
             M = 0
             DO 300 IJ = 1,MIJ
                I = PRIMA (IJ)
                EXPA = ALPHAA (I) + ALPHAA (I)
                DO 310 N = 1,MGQKL
                   EXP2A (M+N) = EXPA
  310           CONTINUE
                M = M + MGQKL
  300        CONTINUE
         END IF

         IF (DIFFB) THEN
             M = 0
             DO 320 IJ = 1,MIJ
                J = PRIMB (IJ)
                EXPB = ALPHAB (J) + ALPHAB (J)
                DO 330 N = 1,MGQKL
                   EXP2B (M+N) = EXPB
  330           CONTINUE
                M = M + MGQKL
  320        CONTINUE
         END IF
C
C
C             ...next the double exponents on centers C and/or D
C                (if needed).
C
C
         IF (DIFFC) THEN
             M = 0
             DO 400 KL = 1,MKL
                K = PRIMC (KL)
                EXPC = ALPHAC (K) + ALPHAC (K)
                DO 410 N = 1,NGQP
                   EXP2C (M+N) = EXPC
  410           CONTINUE
                M = M + NGQP
  400        CONTINUE
             M = 0
             DO 420 IJ = 2,MIJ
                DO 430 N = M+1,M+MGQKL
                   EXP2C (MGQKL+N) = EXP2C (N)
  430           CONTINUE
                M = M + MGQKL
  420        CONTINUE
         END IF

         IF (DIFFD) THEN
             M = 0
             DO 440 KL = 1,MKL
                L = PRIMD (KL)
                EXPD = ALPHAD (L) + ALPHAD (L)
                DO 450 N = 1,NGQP
                   EXP2D (M+N) = EXPD
  450           CONTINUE
                M = M + NGQP
  440        CONTINUE
             M = 0
             DO 460 IJ = 2,MIJ
                DO 470 N = M+1,M+MGQKL
                   EXP2D (MGQKL+N) = EXP2D (N)
  470           CONTINUE
                M = M + MGQKL
  460        CONTINUE
         END IF
C
C
C             ...determine the 2D integral VRR coefficent evaluation
C                case and generate all VRR coefficients. The case is
C                determined as follows:
C
C                    This is done in order to distinguish the maximum
C                    P- and maximum Q-shell combinations for efficient
C                    evaluation of the needed initial 2D integrals
C                    before differentiation and their VRR coefficients.
C                    The cases distinguished are summarized in the
C                    following table, indicating the value of CASE2D:
C
C                                 maximum Q-shell
C                                   s    p   >p
C                                 ---------------
C                                |
C                              s |  1    4    7
C                                |
C             maximum P-shell  p |  2    5    8
C                                |
C                             >p |  3    6    9
C                                |
C
C                    Note, that maximum P- and Q-shells refer to the
C                    overall maximum value of these shells between
C                    the different coordinate components.
C
C
         CASE2D = 3 * MIN0 (2,SHELLQ + MAX0 (DERQX,DERQY,DERQZ))
     +              + MIN0 (2,SHELLP + MAX0 (DERPX,DERPY,DERPZ)) + 1

         CALL    ERD__2D_COEFFICIENTS
     +
     +                ( MIJ,MKL,MIJKL,
     +                  NGQP,MGQIJKL,
     +                  ATOMAB,ATOMCD,
     +                  P,Q,
     +                  PX,PY,PZ,QX,QY,QZ,
     +                  PAX,PAY,PAZ,QCX,QCY,QCZ,
     +                  PINVHF,QINVHF,PQPINV,
     +                  RTS,
     +                  CASE2D,
     +
     +                            B00,B01,B10,
     +                            C00X,C00Y,C00Z,
     +                            D00X,D00Y,D00Z )
     +
     +
C
C
C             ...VRR coefficients are ready. We now assemble the
C                initial 2D integrals and differentiate them for
C                each x,y,z coordinate separately. For each step
C                we determine the 2D integral evaluation cases
C                CASE2D individually. Note, that since the VRR
C                coefficients were evaluated using the overall
C                maximum P- and Q-shell values that can possibly
C                occur, there is no danger in having some VRR
C                coefficients left undefined.
C
C
         IF (PRIMTYP .EQ. 'E0CD') THEN
C
C
C             ...the case where differentiation is to be performed
C                on centers C and/or D only. Generate initial 2DX PCD
C                integrals. If differentiation on x-component is to
C                be done, perform sequence of differentiation to
C                get final differentiated 2DX PCD integrals. Use the
C                BATCH array as scratch.
C
C
             CASE2D = 2 * MIN0 (1,SHELLP) + MIN0 (1,SHELLQ+DERQX) + 1

             CALL    ERD__2D_PCD_INTEGRALS
     +
     +                    ( SHELLP,
     +                      SHELLQ+DERQX,
     +                      MIN0 (SHELLC+DERCX,SHELLD+DERDX),
     +                      SHELLC+DERCX,
     +                      SHELLD+DERDX,
     +                      MGQIJKL,
     +                      WTS,
     +                      B00,B01,B10,C00X,D00X,
     +                      CDX,
     +                      .TRUE.,
     +                      CASE2D,
     +                      BATCH,
     +
     +                                INT2DX )
     +
     +
             IF (DIFFX) THEN

                 DER3 = DERCX
                 DER4 = DERDX

                 DO 500 N = 1,NDERX

                    NINT = MGQIJKL * (SHELLP+1)
     +                             * (SHELLC+DER3+1)
     +                             * (SHELLD+DER4+1)

                    DO 505 M = 1,NINT
                       BATCH (M) = INT2DX (M)
  505               CONTINUE

                    CENTER = CENSQX (N)

                    DEROP3 = CENEQS (3,CENTER)
                    DEROP4 = CENEQS (4,CENTER)

                    DER3 = DER3 - DEROP3
                    DER4 = DER4 - DEROP4

                    CALL    ERD__2D_DERV_PCD_INTEGRALS
     +
     +                           ( SHELLP,
     +                             SHELLC+DER3,
     +                             SHELLD+DER4,
     +                             DEROP3,DEROP4,
     +                             MGQIJKL,
     +                             EXP2C,EXP2D,
     +                             BATCH,
     +
     +                                      INT2DX )
     +
     +
  500            CONTINUE

             END IF
C
C
C             ...generate initial 2DY PCD integrals (if needed).
C                If differentiation on y-component is to be done,
C                perform sequence of differentiation to get final
C                differentiated 2DY PCD integrals. Use the BATCH
C                array as scratch.
C
C
             IF (SHELLT.GT.0 .OR. DIFFY) THEN

                 CASE2D =   2 * MIN0 (1,SHELLP)
     +                        + MIN0 (1,SHELLQ+DERQY) + 1

                 CALL    ERD__2D_PCD_INTEGRALS
     +
     +                        ( SHELLP,
     +                          SHELLQ+DERQY,
     +                          MIN0 (SHELLC+DERCY,SHELLD+DERDY),
     +                          SHELLC+DERCY,
     +                          SHELLD+DERDY,
     +                          MGQIJKL,
     +                          WTS,
     +                          B00,B01,B10,C00Y,D00Y,
     +                          CDY,
     +                          .FALSE.,
     +                          CASE2D,
     +                          BATCH,
     +
     +                                    INT2DY )
     +
     +
             END IF

             IF (DIFFY) THEN

                 DER3 = DERCY
                 DER4 = DERDY

                 DO 510 N = 1,NDERY

                    NINT = MGQIJKL * (SHELLP+1)
     +                             * (SHELLC+DER3+1)
     +                             * (SHELLD+DER4+1)

                    DO 515 M = 1,NINT
                       BATCH (M) = INT2DY (M)
  515               CONTINUE

                    CENTER = CENSQY (N)

                    DEROP3 = CENEQS (3,CENTER)
                    DEROP4 = CENEQS (4,CENTER)

                    DER3 = DER3 - DEROP3
                    DER4 = DER4 - DEROP4

                    CALL    ERD__2D_DERV_PCD_INTEGRALS
     +
     +                           ( SHELLP,
     +                             SHELLC+DER3,
     +                             SHELLD+DER4,
     +                             DEROP3,DEROP4,
     +                             MGQIJKL,
     +                             EXP2C,EXP2D,
     +                             BATCH,
     +
     +                                      INT2DY )
     +
     +
  510            CONTINUE

             END IF
C
C
C             ...generate initial 2DZ PCD integrals. If differentiation
C                on z-component is to be done, perform sequence of
C                differentiation to get final differentiated 2DZ PCD
C                integrals. Use the BATCH array as scratch.
C
C
             IF (SHELLT.GT.0 .OR. DIFFZ) THEN

                 CASE2D =   2 * MIN0 (1,SHELLP)
     +                        + MIN0 (1,SHELLQ+DERQZ) + 1

                 CALL    ERD__2D_PCD_INTEGRALS
     +
     +                        ( SHELLP,
     +                          SHELLQ+DERQZ,
     +                          MIN0 (SHELLC+DERCZ,SHELLD+DERDZ),
     +                          SHELLC+DERCZ,
     +                          SHELLD+DERDZ,
     +                          MGQIJKL,
     +                          WTS,
     +                          B00,B01,B10,C00Z,D00Z,
     +                          CDZ,
     +                          .FALSE.,
     +                          CASE2D,
     +                          BATCH,
     +
     +                                    INT2DZ )
     +
     +
             END IF

             IF (DIFFZ) THEN

                 DER3 = DERCZ
                 DER4 = DERDZ

                 DO 520 N = 1,NDERZ

                    NINT = MGQIJKL * (SHELLP+1)
     +                             * (SHELLC+DER3+1)
     +                             * (SHELLD+DER4+1)

                    DO 525 M = 1,NINT
                       BATCH (M) = INT2DZ (M)
  525               CONTINUE

                    CENTER = CENSQZ (N)

                    DEROP3 = CENEQS (3,CENTER)
                    DEROP4 = CENEQS (4,CENTER)

                    DER3 = DER3 - DEROP3
                    DER4 = DER4 - DEROP4

                    CALL    ERD__2D_DERV_PCD_INTEGRALS
     +
     +                           ( SHELLP,
     +                             SHELLC+DER3,
     +                             SHELLD+DER4,
     +                             DEROP3,DEROP4,
     +                             MGQIJKL,
     +                             EXP2C,EXP2D,
     +                             BATCH,
     +
     +                                      INT2DZ )
     +
     +
  520            CONTINUE

             END IF
C
C
C             ...assemble the 2D PCD integrals to the [E0|CD] batch
C                according to the angular shell case present.
C
C
             IF (SHELLT.EQ.0 .AND. NGQP.EQ.1) THEN

                 DO 530 M = 1,MIJKL
                    BATCH (M) = INT2DX (M) * SCALEPQ (M)
  530            CONTINUE

                 IF (DIFFY) THEN
                     DO 532 M = 1,MIJKL
                        BATCH (M) = BATCH (M) * INT2DY (M)
  532                CONTINUE
                 END IF

                 IF (DIFFZ) THEN
                     DO 534 M = 1,MIJKL
                        BATCH (M) = BATCH (M) * INT2DZ (M)
  534                CONTINUE
                 END IF

             ELSE IF (SHELLT.EQ.0) THEN

                 CALL    ERD__DERV_INT2D_TO_0000
     +
     +                        ( NGQP,MIJKL,MGQIJKL,
     +                          INT2DX,INT2DY,INT2DZ,
     +                          DIFFY,DIFFZ,
     +                          B00,
     +                          SCALEPQ,
     +
     +                                     BATCH )
     +
     +
             ELSE IF (ANGMTYP.EQ.'SSXX') THEN

                 CALL    ERD__DERV_INT2D_TO_E000
     +
     +                        ( SHELLA,SHELLP,
     +                          NGQP,MIJKL,MGQIJKL,
     +                          NXYZET,NXYZP,
     +                          INT2DX,INT2DY,INT2DZ,
     +                          DIFFY,DIFFZ,
     +                          B00,B01,
     +                          SCALEPQ,
     +
     +                                    BATCH )
     +
     +
             ELSE IF (ANGMTYP.EQ.'SXXX') THEN

                 CALL    ERD__DERV_INT2D_TO_E0C0
     +
     +                        ( SHELLA,SHELLP,NZSHELL (1),
     +                          NGQP,MIJKL,MGQIJKL,
     +                          NXYZET,NXYZP,NZNXYZ (1),
     +                          INT2DX,INT2DY,INT2DZ,
     +                          DIFFY,DIFFZ,
     +                          B00,B01,
     +                          SCALEPQ,
     +
     +                                    BATCH )
     +
     +
             ELSE

                 CALL    ERD__DERV_INT2D_TO_E0CD
     +
     +                        ( SHELLA,SHELLP,SHELLC,SHELLD,
     +                          NGQP,MIJKL,MGQIJKL,
     +                          NXYZET,NXYZP,NXYZC,NXYZD,
     +                          INT2DX,INT2DY,INT2DZ,
     +                          DIFFY,DIFFZ,
     +                          B00,B01,
     +                          SCALEPQ,
     +
     +                                    BATCH )
     +
     +
             END IF

         ELSE IF (PRIMTYP .EQ. 'ABCD') THEN
C
C
C             ...the case where differentiation is to be performed on
C                centers A and/or B and on centers C and/or D. Generate
C                initial 2DX ABCD integrals. If differentiation on
C                x-component is to be done, perform sequence of
C                differentiation to get final differentiated 2DX ABCD
C                integrals. Use the BATCH array as scratch.
C
C
             CASE2D = 2*MIN0(1,SHELLP+DERPX) + MIN0(1,SHELLQ+DERQX) + 1

             CALL    ERD__2D_ABCD_INTEGRALS
     +
     +                    ( SHELLP+DERPX,
     +                      SHELLQ+DERQX,
     +                      MIN0 (SHELLA+DERAX,SHELLB+DERBX),
     +                      MIN0 (SHELLC+DERCX,SHELLD+DERDX),
     +                      SHELLA+DERAX,
     +                      SHELLB+DERBX,
     +                      SHELLC+DERCX,
     +                      SHELLD+DERDX,
     +                      MGQIJKL,
     +                      WTS,
     +                      B00,B01,B10,C00X,D00X,
     +                      ABX,CDX,
     +                      .TRUE.,
     +                      CASE2D,
     +                      BATCH,
     +
     +                                INT2DX )
     +
     +
             IF (DIFFX) THEN

                 DER1 = DERAX
                 DER2 = DERBX
                 DER3 = DERCX
                 DER4 = DERDX

                 DO 600 N = 1,NDERX

                    NINT = MGQIJKL * (SHELLA+DER1+1)
     +                             * (SHELLB+DER2+1)
     +                             * (SHELLC+DER3+1)
     +                             * (SHELLD+DER4+1)

                    DO 605 M = 1,NINT
                       BATCH (M) = INT2DX (M)
  605               CONTINUE

                    CENTER = CENSQX (N)

                    DEROP1 = CENEQS (1,CENTER)
                    DEROP2 = CENEQS (2,CENTER)
                    DEROP3 = CENEQS (3,CENTER)
                    DEROP4 = CENEQS (4,CENTER)

                    DER1 = DER1 - DEROP1
                    DER2 = DER2 - DEROP2
                    DER3 = DER3 - DEROP3
                    DER4 = DER4 - DEROP4

                    CALL    ERD__2D_DERV_ABCD_INTEGRALS
     +
     +                           ( SHELLA+DER1,
     +                             SHELLB+DER2,
     +                             SHELLC+DER3,
     +                             SHELLD+DER4,
     +                             DEROP1,DEROP2,DEROP3,DEROP4,
     +                             MGQIJKL,
     +                             EXP2A,EXP2B,EXP2C,EXP2D,
     +                             BATCH,
     +
     +                                      INT2DX )
     +
     +
  600            CONTINUE

             END IF
C
C
C             ...generate initial 2DY ABCD integrals (if needed).
C                If differentiation on y-component is to be done,
C                perform sequence of differentiation to get final
C                differentiated 2DY ABCD integrals. Use the BATCH
C                array as scratch.
C
C
             IF (SHELLT.GT.0 .OR. DIFFY) THEN

                 CASE2D =   2 * MIN0 (1,SHELLP+DERPY)
     +                        + MIN0 (1,SHELLQ+DERQY) + 1

                 CALL    ERD__2D_ABCD_INTEGRALS
     +
     +                        ( SHELLP+DERPY,
     +                          SHELLQ+DERQY,
     +                          MIN0 (SHELLA+DERAY,SHELLB+DERBY),
     +                          MIN0 (SHELLC+DERCY,SHELLD+DERDY),
     +                          SHELLA+DERAY,
     +                          SHELLB+DERBY,
     +                          SHELLC+DERCY,
     +                          SHELLD+DERDY,
     +                          MGQIJKL,
     +                          WTS,
     +                          B00,B01,B10,C00Y,D00Y,
     +                          ABY,CDY,
     +                          .FALSE.,
     +                          CASE2D,
     +                          BATCH,
     +
     +                                    INT2DY )
     +
     +
             END IF

             IF (DIFFY) THEN

                 DER1 = DERAY
                 DER2 = DERBY
                 DER3 = DERCY
                 DER4 = DERDY

                 DO 610 N = 1,NDERY

                    NINT = MGQIJKL * (SHELLA+DER1+1)
     +                             * (SHELLB+DER2+1)
     +                             * (SHELLC+DER3+1)
     +                             * (SHELLD+DER4+1)

                    DO 615 M = 1,NINT
                       BATCH (M) = INT2DY (M)
  615               CONTINUE

                    CENTER = CENSQY (N)

                    DEROP1 = CENEQS (1,CENTER)
                    DEROP2 = CENEQS (2,CENTER)
                    DEROP3 = CENEQS (3,CENTER)
                    DEROP4 = CENEQS (4,CENTER)

                    DER1 = DER1 - DEROP1
                    DER2 = DER2 - DEROP2
                    DER3 = DER3 - DEROP3
                    DER4 = DER4 - DEROP4

                    CALL    ERD__2D_DERV_ABCD_INTEGRALS
     +
     +                           ( SHELLA+DER1,
     +                             SHELLB+DER2,
     +                             SHELLC+DER3,
     +                             SHELLD+DER4,
     +                             DEROP1,DEROP2,DEROP3,DEROP4,
     +                             MGQIJKL,
     +                             EXP2A,EXP2B,EXP2C,EXP2D,
     +                             BATCH,
     +
     +                                      INT2DY )
     +
     +
  610            CONTINUE

             END IF
C
C
C             ...generate initial 2DZ ABCD integrals (if needed).
C                If differentiation on z-component is to be done,
C                perform sequence of differentiation to get final
C                differentiated 2DZ ABCD integrals. Use the BATCH
C                array as scratch.
C
C
             IF (SHELLT.GT.0 .OR. DIFFZ) THEN

                 CASE2D =   2 * MIN0 (1,SHELLP+DERPZ)
     +                        + MIN0 (1,SHELLQ+DERQZ) + 1

                 CALL    ERD__2D_ABCD_INTEGRALS
     +
     +                        ( SHELLP+DERPZ,
     +                          SHELLQ+DERQZ,
     +                          MIN0 (SHELLA+DERAZ,SHELLB+DERBZ),
     +                          MIN0 (SHELLC+DERCZ,SHELLD+DERDZ),
     +                          SHELLA+DERAZ,
     +                          SHELLB+DERBZ,
     +                          SHELLC+DERCZ,
     +                          SHELLD+DERDZ,
     +                          MGQIJKL,
     +                          WTS,
     +                          B00,B01,B10,C00Z,D00Z,
     +                          ABZ,CDZ,
     +                          .FALSE.,
     +                          CASE2D,
     +                          BATCH,
     +
     +                                    INT2DZ )
     +
     +
             END IF

             IF (DIFFZ) THEN

                 DER1 = DERAZ
                 DER2 = DERBZ
                 DER3 = DERCZ
                 DER4 = DERDZ

                 DO 620 N = 1,NDERZ

                    NINT = MGQIJKL * (SHELLA+DER1+1)
     +                             * (SHELLB+DER2+1)
     +                             * (SHELLC+DER3+1)
     +                             * (SHELLD+DER4+1)

                    DO 625 M = 1,NINT
                       BATCH (M) = INT2DZ (M)
  625               CONTINUE

                    CENTER = CENSQZ (N)

                    DEROP1 = CENEQS (1,CENTER)
                    DEROP2 = CENEQS (2,CENTER)
                    DEROP3 = CENEQS (3,CENTER)
                    DEROP4 = CENEQS (4,CENTER)

                    DER1 = DER1 - DEROP1
                    DER2 = DER2 - DEROP2
                    DER3 = DER3 - DEROP3
                    DER4 = DER4 - DEROP4

                    CALL    ERD__2D_DERV_ABCD_INTEGRALS
     +
     +                           ( SHELLA+DER1,
     +                             SHELLB+DER2,
     +                             SHELLC+DER3,
     +                             SHELLD+DER4,
     +                             DEROP1,DEROP2,DEROP3,DEROP4,
     +                             MGQIJKL,
     +                             EXP2A,EXP2B,EXP2C,EXP2D,
     +                             BATCH,
     +
     +                                      INT2DZ )
     +
     +
  620            CONTINUE

             END IF
C
C
C             ...assemble the 2D ABCD integrals to the [AB|CD] batch
C                according to the angular shell case present.
C
C
             IF (ANGMTYP.EQ.'SSSS' .AND. NGQP.EQ.1) THEN

                 DO 630 M = 1,MIJKL
                    BATCH (M) = INT2DX (M) * SCALEPQ (M)
  630            CONTINUE

                 IF (DIFFY) THEN
                     DO 632 M = 1,MIJKL
                        BATCH (M) = BATCH (M) * INT2DY (M)
  632                CONTINUE
                 END IF

                 IF (DIFFZ) THEN
                     DO 634 M = 1,MIJKL
                        BATCH (M) = BATCH (M) * INT2DZ (M)
  634                CONTINUE
                 END IF

             ELSE IF (ANGMTYP.EQ.'SSSS') THEN

                 CALL    ERD__DERV_INT2D_TO_0000
     +
     +                        ( NGQP,MIJKL,MGQIJKL,
     +                          INT2DX,INT2DY,INT2DZ,
     +                          DIFFY,DIFFZ,
     +                          B00,
     +                          SCALEPQ,
     +
     +                                     BATCH )
     +
     +
             ELSE IF (ANGMTYP.EQ.'SSSX') THEN

                 CALL    ERD__DERV_INT2D_TO_A000
     +
     +                        ( NZSHELL (1),
     +                          NGQP,MIJKL,MGQIJKL,
     +                          NZNXYZ (1),
     +                          INT2DX,INT2DY,INT2DZ,
     +                          DIFFY,DIFFZ,
     +                          B00,B01,
     +                          SCALEPQ,
     +
     +                                     BATCH )
     +
     +
             ELSE IF (ANGMTYP.EQ.'SSXX') THEN

                 CALL    ERD__DERV_INT2D_TO_AB00
     +
     +                        ( NZSHELL (1),NZSHELL (2),
     +                          NGQP,MIJKL,MGQIJKL,
     +                          NZNXYZ (1),NZNXYZ (2),
     +                          INT2DX,INT2DY,INT2DZ,
     +                          DIFFY,DIFFZ,
     +                          B00,B01,
     +                          SCALEPQ,
     +
     +                                     BATCH )
     +
     +
             ELSE IF (ANGMTYP.EQ.'SXXX') THEN

                 CALL    ERD__DERV_INT2D_TO_ABC0
     +
     +                        ( NZSHELL (1),NZSHELL (2),NZSHELL (3),
     +                          NGQP,MIJKL,MGQIJKL,
     +                          NZNXYZ (1),NZNXYZ (2),NZNXYZ (3),
     +                          INT2DX,INT2DY,INT2DZ,
     +                          DIFFY,DIFFZ,
     +                          B00,B01,
     +                          SCALEPQ,
     +
     +                                     BATCH )
     +
     +
             ELSE

                 CALL    ERD__DERV_INT2D_TO_ABCD
     +
     +                        ( SHELLA,SHELLB,SHELLC,SHELLD,
     +                          NGQP,MIJKL,MGQIJKL,
     +                          NXYZA,NXYZB,NXYZC,NXYZD,
     +                          INT2DX,INT2DY,INT2DZ,
     +                          DIFFY,DIFFZ,
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
