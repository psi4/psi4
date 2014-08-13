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
         SUBROUTINE  ERD__SPPP_PCGTO_BLOCK
     +
     +                    ( NBATCH,
     +                      ATOMIC,ATOM12,ATOM34,
     +                      MIJ,MKL,
     +                      NIJ,NIJBEG,NIJEND,NKL,NKLBEG,NKLEND,
     +                      NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                      SHELL1,SHELL3,SHELLP,
     +                      X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                      X12,Y12,Z12,X34,Y34,Z34,
     +                      ALPHA1,ALPHA2,ALPHA3,ALPHA4,
     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                      PRIM1,PRIM2,PRIM3,PRIM4,
     +                      NORM1,NORM2,NORM3,NORM4,
     +                      RHO12,RHO34,
     +                      P,PX,PY,PZ,SCALEP,
     +                      Q,QX,QY,QZ,SCALEQ,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__SPPP_PCGTO_BLOCK
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation is designed to provide ultrafast block
C                evaluation of a batch of normalized electron repulsion
C                integrals between s-shell and p-shell primitive
C                spherical gaussian type orbitals.
C
C                A batch is defined here as containing all possible
C                integrals, that is its dimension is determined by
C                the total number of primitive functions (here = 27)
C                times the total number of ij and kl exponent pair
C                combinations.
C
C                The integrals are ordered in the batch the following
C                way (first index varying fastest):
C
C                    batch (nxyz1,nxyz2,nxyz3,nxyz4,kl,ij)
C
C                where ij and kl indicates alpha exponent pairs
C                defining the present block.
C
C                The present routine evaluates batches of the type:
C
C                            sppp , pspp , ppsp , ppps
C
C                The cartesian primitive integrals are evaluated using
C                the auxillary functions technique, described by
C                V.R. Saunders, "An Introduction to Molecular Integral
C                Evaluation" in "Computational Techniques in Quantum
C                Chemistry and Molecular Physics" edited by
C                GHF Diercksen, BT Sutcliff and A Veillard, D. Reidel
C                Publ. Comp. Dordrecht (1975), p. 347.
C
C                The cartesian primitive integrals are each evaluated
C                explicitely using common multipliers if possible and
C                avoiding array addresses.
C
C
C                  Input:
C
C                    NBATCH       =  size of the primitive cartesian
C                                    sppp/pspp/ppsp/ppps integral batch
C                    ATOMIC       =  indicates, if purely atomic
C                                    integrals will be evaluated
C                    ATOM12(34)   =  indicates, if centers 1 and 2
C                                    (3 and 4) coincide
C                    MIJ(KL)      =  current # of ij (kl) primitive
C                                    index pairs corresponding to
C                                    the contracted shell pairs 1,2
C                                    (3,4)
C                    NIJ          =  total # of ij primitive index
C                                    pairs for the contracted shell
C                                    pair 1,2
C                    NIJBEG(END)  =  first(last) ij primitive index
C                                    defining the ij block
C                    NKL          =  total # of kl primitive index
C                                    pairs for the contracted shell
C                                    pair 3,4
C                    NKLBEG(END)  =  first(last) kl primitive index
C                                    defining the kl block
C                    NPGTOx       =  # of primitives per contraction
C                                    for contraction shells x = 1,2,3,4
C                    SHELLx       =  the shell type for contraction
C                                    shells x = 1,3,P=1+2
C                    Xx,Yx,Zx     =  the x,y,z-coordinates for centers
C                                    x = 1,2,3,4
C                    Xxx,Yxx,Zxx  =  the x,y,z-coordinate differences
C                                    between centers xx = 12 and 34
C                    ALPHAx       =  the primitive exponents for
C                                    contraction shells x = 1,2,3,4
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
C                                    x = 1,2,3,4
C                    NORMx        =  the normalization factors due to
C                                    the primitive exponents for the
C                                    contraction shells x = 1,2,3,4
C                    RHO12(34)    =  the complete set of NIJ (NKL)
C                                    exponential prefactors between
C                                    contraction shells 1 and 2
C                                    (3 and 4)
C                    P            =  will hold current MIJ exponent
C                                    sums for contraction shells 1
C                                    and 2
C                    Px           =  will hold current MIJ coordinates
C                                    x=X,Y,Z for the gaussian product
C                                    centers P=1+2
C                    SCALEP       =  will hold current MIJ values of
C                                    scaling factors related to point P
C                    Q            =  will hold current MKL exponent
C                                    sums for contraction shells 3
C                                    and 4
C                    Qx           =  will hold current MKL coordinates
C                                    x=X,Y,Z for the gaussian product
C                                    centers Q=3+4
C                    SCALEQ       =  will hold current MKL values of
C                                    scaling factors related to point Q
C
C                  Output:
C
C                    BATCH        =  current batch of primitive
C                                    cartesian sppp/pspp/ppsp/ppps
C                                    integrals
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
         IMPLICIT    NONE

         LOGICAL     ATOMIC,ATOM12,ATOM34

         INTEGER     CASE
         INTEGER     I,J,K,L,M
         INTEGER     IJ,KL
         INTEGER     MGRID,NGRID,TGRID
         INTEGER     MIJ,MKL
         INTEGER     NBATCH
         INTEGER     NIJ,NKL
         INTEGER     NIJBEG,NIJEND,NKLBEG,NKLEND
         INTEGER     NPGTO1,NPGTO2,NPGTO3,NPGTO4
         INTEGER     SHELL1,SHELL3,SHELLP

         INTEGER     PRIM1 (1:MIJ)
         INTEGER     PRIM2 (1:MIJ)
         INTEGER     PRIM3 (1:MKL)
         INTEGER     PRIM4 (1:MKL)

         DOUBLE PRECISION  A,B,C,D,E,F,G,H,R
         DOUBLE PRECISION  AA,BB
         DOUBLE PRECISION  DELTA1,DELTA2,DELTA3,DELTA4,DELTA5,DELTA6
         DOUBLE PRECISION  EXP1,EXP2,EXP3,EXP4
         DOUBLE PRECISION  F0,F1,F2,F3
         DOUBLE PRECISION  GXXX,GXXY,GXXZ,GXYX,GXYY,GXYZ,GXZX,GXZY,GXZZ,
     +                     GYXX,GYXY,GYXZ,GYYX,GYYY,GYYZ,GYZX,GYZY,GYZZ,
     +                     GZXX,GZXY,GZXZ,GZYX,GZYY,GZYZ,GZZX,GZZY,GZZZ
         DOUBLE PRECISION  PVAL,QVAL,PINV,QINV
         DOUBLE PRECISION  PQPLUS,PQMULT,PQPINV
         DOUBLE PRECISION  PQX,PQY,PQZ
         DOUBLE PRECISION  PSCALE
         DOUBLE PRECISION  PXSUB,PYSUB,PZSUB,QXSUB,QYSUB,QZSUB
         DOUBLE PRECISION  PXVAL,PYVAL,PZVAL,QXVAL,QYVAL,QZVAL
         DOUBLE PRECISION  RNPQSQ
         DOUBLE PRECISION  SCALE
         DOUBLE PRECISION  T,TINV,T2INV,TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION  U0,U1,U2,U3,U4,U5
         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4
         DOUBLE PRECISION  X12,Y12,Z12,X34,Y34,Z34
         DOUBLE PRECISION  XSSP1,XSPS1,XPSS1,XSSP2,XSPS2,XPSS2,
     +                     XSPP1,XPSP1,XPPS1,XSPP2,XPSP2,XPPS2,
     +                     XSPP3,XPSP3,XPPS3,XPPP1,XPPP2,XPPP3,XPPP4,
     +                     YSSP1,YSPS1,YPSS1,YSSP2,YSPS2,YPSS2,
     +                     YSPP1,YPSP1,YPPS1,YSPP2,YPSP2,YPPS2,
     +                     YSPP3,YPSP3,YPPS3,YPPP1,YPPP2,YPPP3,YPPP4,
     +                     ZSSP1,ZSPS1,ZPSS1,ZSSP2,ZSPS2,ZPSS2,
     +                     ZSPP1,ZPSP1,ZPPS1,ZSPP2,ZPSP2,ZPPS2,
     +                     ZSPP3,ZPSP3,ZPPS3,ZPPP1,ZPPP2,ZPPP3,ZPPP4
         DOUBLE PRECISION  FIFTH,THIRD,HALF,ONE,THREE,PI,FIVE

         DOUBLE PRECISION  ALPHA1 (1:NPGTO1)
         DOUBLE PRECISION  ALPHA2 (1:NPGTO2)
         DOUBLE PRECISION  ALPHA3 (1:NPGTO3)
         DOUBLE PRECISION  ALPHA4 (1:NPGTO4)

         DOUBLE PRECISION  NORM1  (1:NPGTO1)
         DOUBLE PRECISION  NORM2  (1:NPGTO2)
         DOUBLE PRECISION  NORM3  (1:NPGTO3)
         DOUBLE PRECISION  NORM4  (1:NPGTO4)

         DOUBLE PRECISION  RHO12  (1:NIJ)
         DOUBLE PRECISION  RHO34  (1:NKL)

         DOUBLE PRECISION  BATCH  (1:NBATCH)

         DOUBLE PRECISION  P      (1:MIJ)
         DOUBLE PRECISION  PX     (1:MIJ)
         DOUBLE PRECISION  PY     (1:MIJ)
         DOUBLE PRECISION  PZ     (1:MIJ)
         DOUBLE PRECISION  SCALEP (1:MIJ)

         DOUBLE PRECISION  Q      (1:MKL)
         DOUBLE PRECISION  QX     (1:MKL)
         DOUBLE PRECISION  QY     (1:MKL)
         DOUBLE PRECISION  QZ     (1:MKL)
         DOUBLE PRECISION  SCALEQ (1:MKL)

         DOUBLE PRECISION  FTABLE (0:MGRID,0:NGRID)

         PARAMETER  (FIFTH   =  0.2D0)
         PARAMETER  (THIRD   =  0.333333333333333D0)
         PARAMETER  (HALF    =  0.5D0)
         PARAMETER  (ONE     =  1.D0 )
         PARAMETER  (THREE   =  3.D0 )
         PARAMETER  (PI      =  3.141592653589793D0)
         PARAMETER  (FIVE    =  5.D0 )
C
C
C------------------------------------------------------------------------
C
C
C             ...predetermine in 'K2' loops the quantities associated
C                with the 1,2-part and 3,4-part.
C
C
         CASE = 5

         IF (ATOMIC.OR.ATOM12) THEN
             M = 0
             DO 100 IJ = NIJBEG,NIJEND
                M = M + 1
                I = PRIM1 (M)
                J = PRIM2 (M)
                P (M) = ALPHA1 (I) + ALPHA2 (J)
                SCALEP (M) = NORM1 (I) * NORM2 (J)
  100        CONTINUE
             CASE = CASE - 1
         ELSE
             M = 0
             DO 110 IJ = NIJBEG,NIJEND
                M = M + 1
                I = PRIM1 (M)
                J = PRIM2 (M)
                EXP1 = ALPHA1 (I)
                EXP2 = ALPHA2 (J)
                PVAL = EXP1 + EXP2
                P (M) = PVAL
                PVAL = EXP1 / PVAL
                PX (M) = PVAL * X12 + X2
                PY (M) = PVAL * Y12 + Y2
                PZ (M) = PVAL * Z12 + Z2
                SCALEP (M) = NORM1 (I) * NORM2 (J) * RHO12 (IJ)
  110        CONTINUE
         END IF

         IF (ATOMIC.OR.ATOM34) THEN
             M = 0
             DO 200 KL = NKLBEG,NKLEND
                M = M + 1
                K = PRIM3 (M)
                L = PRIM4 (M)
                Q (M) = ALPHA3 (K) + ALPHA4 (L)
                SCALEQ (M) = NORM3 (K) * NORM4 (L)
  200        CONTINUE
             CASE = CASE - 2
         ELSE
             M = 0
             DO 220 KL = NKLBEG,NKLEND
                M = M + 1
                K = PRIM3 (M)
                L = PRIM4 (M)
                EXP3 = ALPHA3 (K)
                EXP4 = ALPHA4 (L)
                QVAL = EXP3 + EXP4
                Q (M) = QVAL
                QVAL = EXP3 / QVAL
                QX (M) = QVAL * X34 + X4
                QY (M) = QVAL * Y34 + Y4
                QZ (M) = QVAL * Z34 + Z4
                SCALEQ (M) = NORM3 (K) * NORM4 (L) * RHO34 (KL)
  220        CONTINUE
         END IF

         IF (ATOMIC) THEN
             CASE = CASE - 1
         END IF
C
C
C             ...jump according to where the s-type function is
C                located and the type of 'K4' loop:
C
C                 SHELLP CASE  |        Integral center type
C                --------------|------------------------------------
C                    1     1   |  (AA|AA)  atomic     sppp and pspp
C                    1     2   |  (AA|CC)  2-center   sppp and pspp
C                    1     3   |  (AB|CC)  3-center   sppp and pspp
C                    1     4   |  (AA|CD)  3-center   sppp and pspp
C                    1     5   |  (AB|CD)  4-center   sppp and pspp
C
C                    2     1   |  (AA|AA)  atomic     ppsp and ppps
C                    2     2   |  (AA|CC)  2-center   ppsp and ppps
C                    2     3   |  (AB|CC)  3-center   ppsp and ppps
C                    2     4   |  (AA|CD)  3-center   ppsp and ppps
C                    2     5   |  (AB|CD)  4-center   ppsp and ppps
C
C
         GOTO (11,12,13,14,15,21,22,23,24,25) 5*(SHELLP-1)+CASE
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for atomic sppp and
C                pspp integrals. All these integrals are zero by
C                symmetry.
C
C
   11    RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 2-center sppp and
C                pspp (AA|CC) type integrals.
C
C
   12    PQX = X1 - X3
         PQY = Y1 - Y3
         PQZ = Z1 - Z3
         RNPQSQ = PQX * PQX + PQY * PQY + PQZ * PQZ

         M = 0
         DO 1200 IJ = 1,MIJ
            PVAL = P (IJ)
            PSCALE = SCALEP (IJ)
            DO 1220 KL = 1,MKL
               QVAL = Q (KL)
               QINV = ONE / QVAL
               PQMULT = PVAL * QVAL
               PQPLUS = PVAL + QVAL
               PQPINV = ONE / PQPLUS
               T = RNPQSQ * PQMULT * PQPINV
               SCALE = PSCALE * SCALEQ (KL) / (PQMULT * DSQRT (PQPLUS))

               IF (T.LE.TMAX) THEN
                   TGRID = INT (T * TVSTEP + HALF)
                   DELTA1 = TGRID * TSTEP - T
                   DELTA2 = DELTA1 * HALF
                   DELTA3 = DELTA1 * THIRD
                   DELTA4 = DELTA2 * HALF
                   DELTA5 = DELTA1 * FIFTH
                   DELTA6 = DELTA3 * HALF

                   F1 = ((((( FTABLE (7,TGRID)  * DELTA6
     +                      + FTABLE (6,TGRID)) * DELTA5
     +                      + FTABLE (5,TGRID)) * DELTA4
     +                      + FTABLE (4,TGRID)) * DELTA3
     +                      + FTABLE (3,TGRID)) * DELTA2
     +                      + FTABLE (2,TGRID)) * DELTA1
     +                      + FTABLE (1,TGRID)

                   F2 = ((((( FTABLE (8,TGRID)  * DELTA6
     +                      + FTABLE (7,TGRID)) * DELTA5
     +                      + FTABLE (6,TGRID)) * DELTA4
     +                      + FTABLE (5,TGRID)) * DELTA3
     +                      + FTABLE (4,TGRID)) * DELTA2
     +                      + FTABLE (3,TGRID)) * DELTA1
     +                      + FTABLE (2,TGRID)

                   F3 = ((((( FTABLE (9,TGRID)  * DELTA6
     +                      + FTABLE (8,TGRID)) * DELTA5
     +                      + FTABLE (7,TGRID)) * DELTA4
     +                      + FTABLE (6,TGRID)) * DELTA3
     +                      + FTABLE (5,TGRID)) * DELTA2
     +                      + FTABLE (4,TGRID)) * DELTA1
     +                      + FTABLE (3,TGRID)

                   F1 = SCALE * F1
                   F2 = SCALE * F2
                   F3 = SCALE * F3

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F1 = T2INV * SCALE * HALF * DSQRT (PI*TINV)
                   F2 = THREE * T2INV * F1
                   F3 = FIVE  * T2INV * F2

               END IF

               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               U2 = HALF * PQPINV
               U3 = U2 + PQPINV
               U4 = HALF * QINV
C
C
C             ...the X-terms.
C
C
               XSSP2 = PQX * U0
               XPSS2 = PQX * U1
               XSPP3 = XSSP2 * XSSP2
               XPSP3 = XSSP2 * XPSS2
               XPPP2 = U4 * XPSS2
               XPPP3 = U3 * XSSP2
               XPPP4 = XPSS2 * XSPP3
C
C
C             ...the Y-terms.
C
C
               YSSP2 = PQY * U0
               YPSS2 = PQY * U1
               YSPP3 = YSSP2 * YSSP2
               YPSP3 = YSSP2 * YPSS2
               YPPP2 = U4 * YPSS2
               YPPP3 = U3 * YSSP2
               YPPP4 = YPSS2 * YSPP3
C
C
C             ...the Z-terms.
C
C
               ZSSP2 = PQZ * U0
               ZPSS2 = PQZ * U1
               ZSPP3 = ZSSP2 * ZSSP2
               ZPSP3 = ZSSP2 * ZPSS2
               ZPPP2 = U4 * ZPSS2
               ZPPP3 = U3 * ZSSP2
               ZPPP4 = ZPSS2 * ZSPP3
C
C
C             ...assemble the 2-center (AA|CC) type integrals.
C
C
               GXXX = XPPP2 * F1 + XPPP3 * F2 + XPPP4 * F3
               GYYY = YPPP2 * F1 + YPPP3 * F2 + YPPP4 * F3
               GZZZ = ZPPP2 * F1 + ZPPP3 * F2 + ZPPP4 * F3

               A = U4 * (F1 - U0 * F2)
               B = U2 * F2
               C = B + XPSP3 * F3
               D = A + XSPP3 * F3

               GXXY = YSSP2 * C
               GXXZ = ZSSP2 * C
               GYXX = YPSS2 * D
               GZXX = ZPSS2 * D

               C = B + YPSP3 * F3
               D = A + YSPP3 * F3

               GYXY = XSSP2 * C
               GYYZ = ZSSP2 * C
               GXYY = XPSS2 * D
               GZYY = ZPSS2 * D

               C = B + ZPSP3 * F3
               D = A + ZSPP3 * F3

               GZXZ = XSSP2 * C
               GZYZ = YSSP2 * C
               GXZZ = XPSS2 * D
               GYZZ = YPSS2 * D

               GXYZ = YSSP2 * ZSSP2 * XPSS2 * F3
               GYXZ = XSSP2 * ZSSP2 * YPSS2 * F3
               GZXY = XSSP2 * YSSP2 * ZPSS2 * F3

               BATCH (M+ 1) = GXXX
               BATCH (M+ 2) = GYXX
               BATCH (M+ 3) = GZXX
               BATCH (M+ 4) = GXXY
               BATCH (M+ 5) = GYXY
               BATCH (M+ 6) = GZXY
               BATCH (M+ 7) = GXXZ
               BATCH (M+ 8) = GYXZ
               BATCH (M+ 9) = GZXZ
               BATCH (M+10) = GXXY
               BATCH (M+11) = GYXY
               BATCH (M+12) = GZXY
               BATCH (M+13) = GXYY
               BATCH (M+14) = GYYY
               BATCH (M+15) = GZYY
               BATCH (M+16) = GXYZ
               BATCH (M+17) = GYYZ
               BATCH (M+18) = GZYZ
               BATCH (M+19) = GXXZ
               BATCH (M+20) = GYXZ
               BATCH (M+21) = GZXZ
               BATCH (M+22) = GXYZ
               BATCH (M+23) = GYYZ
               BATCH (M+24) = GZYZ
               BATCH (M+25) = GXZZ
               BATCH (M+26) = GYZZ
               BATCH (M+27) = GZZZ

               M = M + 27

 1220       CONTINUE
 1200    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 3-center sppp and
C                pspp (AB|CC) type integrals.
C
C
   13    IF (SHELL1.EQ.1) THEN
             PXSUB = X1
             PYSUB = Y1
             PZSUB = Z1
         ELSE
             PXSUB = X2
             PYSUB = Y2
             PZSUB = Z2
         END IF

         M = 0
         DO 1300 IJ = 1,MIJ
            PVAL = P (IJ)
            PXVAL = PX (IJ)
            PYVAL = PY (IJ)
            PZVAL = PZ (IJ)
            PQX = PXVAL - X3
            PQY = PYVAL - Y3
            PQZ = PZVAL - Z3
            RNPQSQ = PQX * PQX + PQY * PQY + PQZ * PQZ 
            PSCALE = SCALEP (IJ)
            XPSS1 = PXVAL - PXSUB
            YPSS1 = PYVAL - PYSUB
            ZPSS1 = PZVAL - PZSUB
            DO 1330 KL = 1,MKL
               QVAL = Q (KL)
               QINV = ONE / QVAL
               PQMULT = PVAL * QVAL
               PQPLUS = PVAL + QVAL
               PQPINV = ONE / PQPLUS
               T = RNPQSQ * PQMULT * PQPINV
               SCALE = PSCALE * SCALEQ (KL) / (PQMULT * DSQRT (PQPLUS))

               IF (T.LE.TMAX) THEN
                   TGRID = INT (T * TVSTEP + HALF)
                   DELTA1 = TGRID * TSTEP - T
                   DELTA2 = DELTA1 * HALF
                   DELTA3 = DELTA1 * THIRD
                   DELTA4 = DELTA2 * HALF
                   DELTA5 = DELTA1 * FIFTH
                   DELTA6 = DELTA3 * HALF

                   F0 = ((((( FTABLE (6,TGRID)  * DELTA6
     +                      + FTABLE (5,TGRID)) * DELTA5
     +                      + FTABLE (4,TGRID)) * DELTA4
     +                      + FTABLE (3,TGRID)) * DELTA3
     +                      + FTABLE (2,TGRID)) * DELTA2
     +                      + FTABLE (1,TGRID)) * DELTA1
     +                      + FTABLE (0,TGRID)

                   F1 = ((((( FTABLE (7,TGRID)  * DELTA6
     +                      + FTABLE (6,TGRID)) * DELTA5
     +                      + FTABLE (5,TGRID)) * DELTA4
     +                      + FTABLE (4,TGRID)) * DELTA3
     +                      + FTABLE (3,TGRID)) * DELTA2
     +                      + FTABLE (2,TGRID)) * DELTA1
     +                      + FTABLE (1,TGRID)

                   F2 = ((((( FTABLE (8,TGRID)  * DELTA6
     +                      + FTABLE (7,TGRID)) * DELTA5
     +                      + FTABLE (6,TGRID)) * DELTA4
     +                      + FTABLE (5,TGRID)) * DELTA3
     +                      + FTABLE (4,TGRID)) * DELTA2
     +                      + FTABLE (3,TGRID)) * DELTA1
     +                      + FTABLE (2,TGRID)

                   F3 = ((((( FTABLE (9,TGRID)  * DELTA6
     +                      + FTABLE (8,TGRID)) * DELTA5
     +                      + FTABLE (7,TGRID)) * DELTA4
     +                      + FTABLE (6,TGRID)) * DELTA3
     +                      + FTABLE (5,TGRID)) * DELTA2
     +                      + FTABLE (4,TGRID)) * DELTA1
     +                      + FTABLE (3,TGRID)

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2
                   F3 = SCALE * F3

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1
                   F3 = FIVE  * T2INV * F2

               END IF

               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               U2 = HALF * PQPINV
               U3 = U2 + PQPINV
               U4 = HALF * QINV
               U5 = - U0 * U4
C
C
C             ...the X-terms.
C
C
               XSSP2 = PQX * U0
               XPSS2 = PQX * U1
               XPSP2 = XPSS1 * XSSP2 + U2
               XSPP3 = XSSP2 * XSSP2
               XPSP3 = XSSP2 * XPSS2
               XPPP1 = XPSS1 * U4
               XPPP2 = XPSS1 * U5 + U4 * XPSS2
               XPPP3 = XPSS1 * XSPP3 + U3 * XSSP2
               XPPP4 = XPSS2 * XSPP3
C
C
C             ...the Y-terms.
C
C
               YSSP2 = PQY * U0
               YPSS2 = PQY * U1
               YPSP2 = YPSS1 * YSSP2 + U2
               YSPP3 = YSSP2 * YSSP2
               YPSP3 = YSSP2 * YPSS2
               YPPP1 = YPSS1 * U4
               YPPP2 = YPSS1 * U5 + U4 * YPSS2
               YPPP3 = YPSS1 * YSPP3 + U3 * YSSP2
               YPPP4 = YPSS2 * YSPP3
C
C
C             ...the Z-terms.
C
C
               ZSSP2 = PQZ * U0
               ZPSS2 = PQZ * U1
               ZPSP2 = ZPSS1 * ZSSP2 + U2
               ZSPP3 = ZSSP2 * ZSSP2
               ZPSP3 = ZSSP2 * ZPSS2
               ZPPP1 = ZPSS1 * U4
               ZPPP2 = ZPSS1 * U5 + U4 * ZPSS2
               ZPPP3 = ZPSS1 * ZSPP3 + U3 * ZSSP2
               ZPPP4 = ZPSS2 * ZSPP3
C
C
C             ...assemble the 3-center (AB|CC) type integrals.
C
C
               GXXX = XPPP1 * F0 + XPPP2 * F1 + XPPP3 * F2 + XPPP4 * F3
               GYYY = YPPP1 * F0 + YPPP2 * F1 + YPPP3 * F2 + YPPP4 * F3
               GZZZ = ZPPP1 * F0 + ZPPP2 * F1 + ZPPP3 * F2 + ZPPP4 * F3

               A = U4 * F0 + U5 * F1
               B = U4 * F1 + U5 * F2
               C = XPSP2 * F2 + XPSP3 * F3
               D = A + XSPP3 * F2
               E = B + XSPP3 * F3

               GXXY = YSSP2 * C
               GXXZ = ZSSP2 * C
               GYXX = YPSS1 * D + YPSS2 * E
               GZXX = ZPSS1 * D + ZPSS2 * E

               C = YPSP2 * F2 + YPSP3 * F3
               D = A + YSPP3 * F2
               E = B + YSPP3 * F3

               GYXY = XSSP2 * C
               GYYZ = ZSSP2 * C
               GXYY = XPSS1 * D + XPSS2 * E
               GZYY = ZPSS1 * D + ZPSS2 * E

               C = ZPSP2 * F2 + ZPSP3 * F3
               D = A + ZSPP3 * F2
               E = B + ZSPP3 * F3

               GZXZ = XSSP2 * C
               GZYZ = YSSP2 * C
               GXZZ = XPSS1 * D + XPSS2 * E
               GYZZ = YPSS1 * D + YPSS2 * E

               GXYZ = YSSP2 * ZSSP2 * (XPSS1 * F2 + XPSS2 * F3)
               GYXZ = XSSP2 * ZSSP2 * (YPSS1 * F2 + YPSS2 * F3)
               GZXY = XSSP2 * YSSP2 * (ZPSS1 * F2 + ZPSS2 * F3)

               BATCH (M+ 1) = GXXX
               BATCH (M+ 2) = GYXX
               BATCH (M+ 3) = GZXX
               BATCH (M+ 4) = GXXY
               BATCH (M+ 5) = GYXY
               BATCH (M+ 6) = GZXY
               BATCH (M+ 7) = GXXZ
               BATCH (M+ 8) = GYXZ
               BATCH (M+ 9) = GZXZ
               BATCH (M+10) = GXXY
               BATCH (M+11) = GYXY
               BATCH (M+12) = GZXY
               BATCH (M+13) = GXYY
               BATCH (M+14) = GYYY
               BATCH (M+15) = GZYY
               BATCH (M+16) = GXYZ
               BATCH (M+17) = GYYZ
               BATCH (M+18) = GZYZ
               BATCH (M+19) = GXXZ
               BATCH (M+20) = GYXZ
               BATCH (M+21) = GZXZ
               BATCH (M+22) = GXYZ
               BATCH (M+23) = GYYZ
               BATCH (M+24) = GZYZ
               BATCH (M+25) = GXZZ
               BATCH (M+26) = GYZZ
               BATCH (M+27) = GZZZ

               M = M + 27

 1330       CONTINUE
 1300    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 3-center sppp and
C                pspp (AA|CD) type integrals.
C
C
   14    M = 0
         DO 1400 IJ = 1,MIJ
            PVAL = P (IJ)
            PSCALE = SCALEP (IJ)
            DO 1440 KL = 1,MKL
               QVAL = Q (KL)
               QINV = ONE / QVAL
               QXVAL = QX (KL)
               QYVAL = QY (KL)
               QZVAL = QZ (KL)
               PQMULT = PVAL * QVAL
               PQPLUS = PVAL + QVAL
               PQPINV = ONE / PQPLUS
               PQX = X1 - QXVAL
               PQY = Y1 - QYVAL
               PQZ = Z1 - QZVAL
               T = (PQX*PQX + PQY*PQY + PQZ*PQZ) * PQMULT * PQPINV
               SCALE = PSCALE * SCALEQ (KL) / (PQMULT * DSQRT (PQPLUS))

               IF (T.LE.TMAX) THEN
                   TGRID = INT (T * TVSTEP + HALF)
                   DELTA1 = TGRID * TSTEP - T
                   DELTA2 = DELTA1 * HALF
                   DELTA3 = DELTA1 * THIRD
                   DELTA4 = DELTA2 * HALF
                   DELTA5 = DELTA1 * FIFTH
                   DELTA6 = DELTA3 * HALF

                   F1 = ((((( FTABLE (7,TGRID)  * DELTA6
     +                      + FTABLE (6,TGRID)) * DELTA5
     +                      + FTABLE (5,TGRID)) * DELTA4
     +                      + FTABLE (4,TGRID)) * DELTA3
     +                      + FTABLE (3,TGRID)) * DELTA2
     +                      + FTABLE (2,TGRID)) * DELTA1
     +                      + FTABLE (1,TGRID)

                   F2 = ((((( FTABLE (8,TGRID)  * DELTA6
     +                      + FTABLE (7,TGRID)) * DELTA5
     +                      + FTABLE (6,TGRID)) * DELTA4
     +                      + FTABLE (5,TGRID)) * DELTA3
     +                      + FTABLE (4,TGRID)) * DELTA2
     +                      + FTABLE (3,TGRID)) * DELTA1
     +                      + FTABLE (2,TGRID)

                   F3 = ((((( FTABLE (9,TGRID)  * DELTA6
     +                      + FTABLE (8,TGRID)) * DELTA5
     +                      + FTABLE (7,TGRID)) * DELTA4
     +                      + FTABLE (6,TGRID)) * DELTA3
     +                      + FTABLE (5,TGRID)) * DELTA2
     +                      + FTABLE (4,TGRID)) * DELTA1
     +                      + FTABLE (3,TGRID)

                   F1 = SCALE * F1
                   F2 = SCALE * F2
                   F3 = SCALE * F3

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F1 = T2INV * SCALE * HALF * DSQRT (PI*TINV)
                   F2 = THREE * T2INV * F1
                   F3 = FIVE  * T2INV * F2

               END IF

               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               U2 = HALF * PQPINV
               U3 = U2 + PQPINV
               U4 = HALF * QINV
               U5 = U0 * U4
C
C
C             ...the X-terms.
C
C
               XSSP1 = QXVAL - X4
               XSPS1 = QXVAL - X3
               XSSP2 = PQX * U0
               XPSS2 = PQX * U1
               A     = XSPS1 + XSSP1
               XSPP1 = XSPS1 * XSSP1 + U4
               XSPP2 = A * XSSP2 - U5
               XSPP3 = XSSP2 * XSSP2
               XPSP2 = XSSP1 * XPSS2 + U2
               XPSP3 = XSSP2 * XPSS2
               XPPS2 = XSPS1 * XPSS2 + U2
               XPPP2 = XSPP1 * XPSS2 + A * U2
               XPPP3 = A * XPSP3 + U3 * XSSP2
               XPPP4 = XPSS2 * XSPP3
C
C
C             ...the Y-terms.
C
C
               YSSP1 = QYVAL - Y4
               YSPS1 = QYVAL - Y3
               YSSP2 = PQY * U0
               YPSS2 = PQY * U1
               A     = YSPS1 + YSSP1
               YSPP1 = YSPS1 * YSSP1 + U4
               YSPP2 = A * YSSP2 - U5
               YSPP3 = YSSP2 * YSSP2
               YPSP2 = YSSP1 * YPSS2 + U2
               YPSP3 = YSSP2 * YPSS2
               YPPS2 = YSPS1 * YPSS2 + U2
               YPPP2 = YSPP1 * YPSS2 + A * U2
               YPPP3 = A * YPSP3 + U3 * YSSP2
               YPPP4 = YPSS2 * YSPP3
C
C
C             ...the Z-terms.
C
C
               ZSSP1 = QZVAL - Z4
               ZSPS1 = QZVAL - Z3
               ZSSP2 = PQZ * U0
               ZPSS2 = PQZ * U1
               A     = ZSPS1 + ZSSP1
               ZSPP1 = ZSPS1 * ZSSP1 + U4
               ZSPP2 = A * ZSSP2 - U5
               ZSPP3 = ZSSP2 * ZSSP2
               ZPSP2 = ZSSP1 * ZPSS2 + U2
               ZPSP3 = ZSSP2 * ZPSS2
               ZPPS2 = ZSPS1 * ZPSS2 + U2
               ZPPP2 = ZSPP1 * ZPSS2 + A * U2
               ZPPP3 = A * ZPSP3 + U3 * ZSSP2
               ZPPP4 = ZPSS2 * ZSPP3
C
C
C             ...assemble the 3-center (AA|CD) type integrals.
C
C
               GXXX = XPPP2 * F1 + XPPP3 * F2 + XPPP4 * F3
               GYYY = YPPP2 * F1 + YPPP3 * F2 + YPPP4 * F3
               GZZZ = ZPPP2 * F1 + ZPPP3 * F2 + ZPPP4 * F3

               AA = XPSP3 * F2
               BB = XPSP3 * F3

               A = XPPS2 * F1 + AA
               B = XPPS2 * F2 + BB
               C = XPSP2 * F1 + AA
               D = XPSP2 * F2 + BB
               F = XSPP1 * F1 + XSPP2 * F2 + XSPP3 * F3

               GXXY = YSSP1 * A + YSSP2 * B
               GXXZ = ZSSP1 * A + ZSSP2 * B
               GXYX = YSPS1 * C + YSSP2 * D
               GXZX = ZSPS1 * C + ZSSP2 * D
               GYXX = YPSS2 * F
               GZXX = ZPSS2 * F

               AA = YPSP3 * F2
               BB = YPSP3 * F3

               A = YPPS2 * F1 + AA
               B = YPPS2 * F2 + BB
               C = YPSP2 * F1 + AA
               D = YPSP2 * F2 + BB
               F = YSPP1 * F1 + YSPP2 * F2 + YSPP3 * F3

               GYYX = XSSP1 * A + XSSP2 * B
               GYYZ = ZSSP1 * A + ZSSP2 * B
               GYXY = XSPS1 * C + XSSP2 * D
               GYZY = ZSPS1 * C + ZSSP2 * D
               GXYY = XPSS2 * F
               GZYY = ZPSS2 * F

               AA = ZPSP3 * F2
               BB = ZPSP3 * F3

               A = ZPPS2 * F1 + AA
               B = ZPPS2 * F2 + BB
               C = ZPSP2 * F1 + AA
               D = ZPSP2 * F2 + BB
               F = ZSPP1 * F1 + ZSPP2 * F2 + ZSPP3 * F3

               GZZX = XSSP1 * A + XSSP2 * B
               GZZY = YSSP1 * A + YSSP2 * B
               GZXZ = XSPS1 * C + XSSP2 * D
               GZYZ = YSPS1 * C + YSSP2 * D
               GXZZ = XPSS2 * F
               GYZZ = YPSS2 * F

               A = XSSP1 * F1 + XSSP2 * F2
               B = XSSP1 * F2 + XSSP2 * F3
               C = YSSP1 * F1 + YSSP2 * F2
               D = YSSP1 * F2 + YSSP2 * F3
               E = ZSSP1 * F1 + ZSSP2 * F2
               F = ZSSP1 * F2 + ZSSP2 * F3

               GXYZ = (YSPS1 * E + YSSP2 * F) * XPSS2
               GXZY = (ZSPS1 * C + ZSSP2 * D) * XPSS2
               GYXZ = (XSPS1 * E + XSSP2 * F) * YPSS2
               GYZX = (ZSPS1 * A + ZSSP2 * B) * YPSS2
               GZXY = (XSPS1 * C + XSSP2 * D) * ZPSS2
               GZYX = (YSPS1 * A + YSSP2 * B) * ZPSS2

               BATCH (M+ 1) = GXXX
               BATCH (M+ 2) = GYXX
               BATCH (M+ 3) = GZXX
               BATCH (M+ 4) = GXYX
               BATCH (M+ 5) = GYYX
               BATCH (M+ 6) = GZYX
               BATCH (M+ 7) = GXZX
               BATCH (M+ 8) = GYZX
               BATCH (M+ 9) = GZZX
               BATCH (M+10) = GXXY
               BATCH (M+11) = GYXY
               BATCH (M+12) = GZXY
               BATCH (M+13) = GXYY
               BATCH (M+14) = GYYY
               BATCH (M+15) = GZYY
               BATCH (M+16) = GXZY
               BATCH (M+17) = GYZY
               BATCH (M+18) = GZZY
               BATCH (M+19) = GXXZ
               BATCH (M+20) = GYXZ
               BATCH (M+21) = GZXZ
               BATCH (M+22) = GXYZ
               BATCH (M+23) = GYYZ
               BATCH (M+24) = GZYZ
               BATCH (M+25) = GXZZ
               BATCH (M+26) = GYZZ
               BATCH (M+27) = GZZZ

               M = M + 27

 1440       CONTINUE
 1400    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 4-center sppp and
C                pspp (AB|CD) type integrals.
C
C
   15    IF (SHELL1.EQ.1) THEN
             PXSUB = X1
             PYSUB = Y1
             PZSUB = Z1
         ELSE
             PXSUB = X2
             PYSUB = Y2
             PZSUB = Z2
         END IF

         M = 0
         DO 1500 IJ = 1,MIJ
            PVAL = P (IJ)
            PXVAL = PX (IJ)
            PYVAL = PY (IJ)
            PZVAL = PZ (IJ)
            PSCALE = SCALEP (IJ)
            XPSS1 = PXVAL - PXSUB
            YPSS1 = PYVAL - PYSUB
            ZPSS1 = PZVAL - PZSUB
            DO 1550 KL = 1,MKL
               QVAL = Q (KL)
               QINV = ONE / QVAL
               QXVAL = QX (KL)
               QYVAL = QY (KL)
               QZVAL = QZ (KL)
               PQMULT = PVAL * QVAL
               PQPLUS = PVAL + QVAL
               PQPINV = ONE / PQPLUS
               PQX = PXVAL - QXVAL
               PQY = PYVAL - QYVAL
               PQZ = PZVAL - QZVAL
               T = (PQX*PQX + PQY*PQY + PQZ*PQZ) * PQMULT * PQPINV
               SCALE = PSCALE * SCALEQ (KL) / (PQMULT * DSQRT (PQPLUS))

               IF (T.LE.TMAX) THEN
                   TGRID = INT (T * TVSTEP + HALF)
                   DELTA1 = TGRID * TSTEP - T
                   DELTA2 = DELTA1 * HALF
                   DELTA3 = DELTA1 * THIRD
                   DELTA4 = DELTA2 * HALF
                   DELTA5 = DELTA1 * FIFTH
                   DELTA6 = DELTA3 * HALF

                   F0 = ((((( FTABLE (6,TGRID)  * DELTA6
     +                      + FTABLE (5,TGRID)) * DELTA5
     +                      + FTABLE (4,TGRID)) * DELTA4
     +                      + FTABLE (3,TGRID)) * DELTA3
     +                      + FTABLE (2,TGRID)) * DELTA2
     +                      + FTABLE (1,TGRID)) * DELTA1
     +                      + FTABLE (0,TGRID)

                   F1 = ((((( FTABLE (7,TGRID)  * DELTA6
     +                      + FTABLE (6,TGRID)) * DELTA5
     +                      + FTABLE (5,TGRID)) * DELTA4
     +                      + FTABLE (4,TGRID)) * DELTA3
     +                      + FTABLE (3,TGRID)) * DELTA2
     +                      + FTABLE (2,TGRID)) * DELTA1
     +                      + FTABLE (1,TGRID)

                   F2 = ((((( FTABLE (8,TGRID)  * DELTA6
     +                      + FTABLE (7,TGRID)) * DELTA5
     +                      + FTABLE (6,TGRID)) * DELTA4
     +                      + FTABLE (5,TGRID)) * DELTA3
     +                      + FTABLE (4,TGRID)) * DELTA2
     +                      + FTABLE (3,TGRID)) * DELTA1
     +                      + FTABLE (2,TGRID)

                   F3 = ((((( FTABLE (9,TGRID)  * DELTA6
     +                      + FTABLE (8,TGRID)) * DELTA5
     +                      + FTABLE (7,TGRID)) * DELTA4
     +                      + FTABLE (6,TGRID)) * DELTA3
     +                      + FTABLE (5,TGRID)) * DELTA2
     +                      + FTABLE (4,TGRID)) * DELTA1
     +                      + FTABLE (3,TGRID)

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2
                   F3 = SCALE * F3

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1
                   F3 = FIVE  * T2INV * F2

               END IF

               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               U2 = HALF * PQPINV
               U3 = U2 + PQPINV
               U4 = HALF * QINV
               U5 = U0 * U4
C
C
C             ...the X-terms.
C
C
               XSSP1 = QXVAL - X4
               XSPS1 = QXVAL - X3
               XSSP2 = PQX * U0
               XPSS2 = PQX * U1
               A     = XSPS1 + XSSP1
               B     = XPSS1 * XSSP2 + U2
               XSPP1 = XSPS1 * XSSP1 + U4
               XSPP2 = A * XSSP2 - U5
               XSPP3 = XSSP2 * XSSP2
               XPSP1 = XPSS1 * XSSP1
               XPSP2 = XSSP1 * XPSS2 + B
               XPSP3 = XSSP2 * XPSS2
               XPPS1 = XPSS1 * XSPS1
               XPPS2 = XSPS1 * XPSS2 + B
               XPPP1 = XPSS1 * XSPP1
               XPPP2 = XPSS1 * XSPP2 + XSPP1 * XPSS2 + A * U2
               XPPP3 = XPSS1 * XSPP3 + A * XPSP3 + U3 * XSSP2
               XPPP4 = XPSS2 * XSPP3
C
C
C             ...the Y-terms.
C
C
               YSSP1 = QYVAL - Y4
               YSPS1 = QYVAL - Y3
               YSSP2 = PQY * U0
               YPSS2 = PQY * U1
               A     = YSPS1 + YSSP1
               B     = YPSS1 * YSSP2 + U2
               YSPP1 = YSPS1 * YSSP1 + U4
               YSPP2 = A * YSSP2 - U5
               YSPP3 = YSSP2 * YSSP2
               YPSP1 = YPSS1 * YSSP1
               YPSP2 = YSSP1 * YPSS2 + B
               YPSP3 = YSSP2 * YPSS2
               YPPS1 = YPSS1 * YSPS1
               YPPS2 = YSPS1 * YPSS2 + B
               YPPP1 = YPSS1 * YSPP1
               YPPP2 = YPSS1 * YSPP2 + YSPP1 * YPSS2 + A * U2
               YPPP3 = YPSS1 * YSPP3 + A * YPSP3 + U3 * YSSP2
               YPPP4 = YPSS2 * YSPP3
C
C
C             ...the Z-terms.
C
C
               ZSSP1 = QZVAL - Z4
               ZSPS1 = QZVAL - Z3
               ZSSP2 = PQZ * U0
               ZPSS2 = PQZ * U1
               A     = ZSPS1 + ZSSP1
               B     = ZPSS1 * ZSSP2 + U2
               ZSPP1 = ZSPS1 * ZSSP1 + U4
               ZSPP2 = A * ZSSP2 - U5
               ZSPP3 = ZSSP2 * ZSSP2
               ZPSP1 = ZPSS1 * ZSSP1
               ZPSP2 = ZSSP1 * ZPSS2 + B
               ZPSP3 = ZSSP2 * ZPSS2
               ZPPS1 = ZPSS1 * ZSPS1
               ZPPS2 = ZSPS1 * ZPSS2 + B
               ZPPP1 = ZPSS1 * ZSPP1
               ZPPP2 = ZPSS1 * ZSPP2 + ZSPP1 * ZPSS2 + A * U2
               ZPPP3 = ZPSS1 * ZSPP3 + A * ZPSP3 + U3 * ZSSP2
               ZPPP4 = ZPSS2 * ZSPP3
C
C
C             ...assemble the 4-center (AB|CD) type integrals.
C
C
               GXXX = XPPP1 * F0 + XPPP2 * F1 + XPPP3 * F2 + XPPP4 * F3
               GYYY = YPPP1 * F0 + YPPP2 * F1 + YPPP3 * F2 + YPPP4 * F3
               GZZZ = ZPPP1 * F0 + ZPPP2 * F1 + ZPPP3 * F2 + ZPPP4 * F3

               AA = XPSP3 * F2
               BB = XPSP3 * F3

               A = XPPS1 * F0 + XPPS2 * F1 + AA
               B = XPPS1 * F1 + XPPS2 * F2 + BB
               C = XPSP1 * F0 + XPSP2 * F1 + AA
               D = XPSP1 * F1 + XPSP2 * F2 + BB
               E = XSPP1 * F0 + XSPP2 * F1 + XSPP3 * F2
               F = XSPP1 * F1 + XSPP2 * F2 + XSPP3 * F3

               GXXY = YSSP1 * A + YSSP2 * B
               GXXZ = ZSSP1 * A + ZSSP2 * B
               GXYX = YSPS1 * C + YSSP2 * D
               GXZX = ZSPS1 * C + ZSSP2 * D
               GYXX = YPSS1 * E + YPSS2 * F
               GZXX = ZPSS1 * E + ZPSS2 * F

               AA = YPSP3 * F2
               BB = YPSP3 * F3

               A = YPPS1 * F0 + YPPS2 * F1 + AA
               B = YPPS1 * F1 + YPPS2 * F2 + BB
               C = YPSP1 * F0 + YPSP2 * F1 + AA
               D = YPSP1 * F1 + YPSP2 * F2 + BB
               E = YSPP1 * F0 + YSPP2 * F1 + YSPP3 * F2
               F = YSPP1 * F1 + YSPP2 * F2 + YSPP3 * F3

               GYYX = XSSP1 * A + XSSP2 * B
               GYYZ = ZSSP1 * A + ZSSP2 * B
               GYXY = XSPS1 * C + XSSP2 * D
               GYZY = ZSPS1 * C + ZSSP2 * D
               GXYY = XPSS1 * E + XPSS2 * F
               GZYY = ZPSS1 * E + ZPSS2 * F

               AA = ZPSP3 * F2
               BB = ZPSP3 * F3

               A = ZPPS1 * F0 + ZPPS2 * F1 + AA
               B = ZPPS1 * F1 + ZPPS2 * F2 + BB
               C = ZPSP1 * F0 + ZPSP2 * F1 + AA
               D = ZPSP1 * F1 + ZPSP2 * F2 + BB
               E = ZSPP1 * F0 + ZSPP2 * F1 + ZSPP3 * F2
               F = ZSPP1 * F1 + ZSPP2 * F2 + ZSPP3 * F3

               GZZX = XSSP1 * A + XSSP2 * B
               GZZY = YSSP1 * A + YSSP2 * B
               GZXZ = XSPS1 * C + XSSP2 * D
               GZYZ = YSPS1 * C + YSSP2 * D
               GXZZ = XPSS1 * E + XPSS2 * F
               GYZZ = YPSS1 * E + YPSS2 * F

               A = XPSS1 * F0 + XPSS2 * F1
               B = XPSS1 * F1 + XPSS2 * F2
               C = XPSS1 * F2 + XPSS2 * F3
               D = YPSS1 * F0 + YPSS2 * F1
               E = YPSS1 * F1 + YPSS2 * F2
               F = YPSS1 * F2 + YPSS2 * F3
               G = ZPSS1 * F0 + ZPSS2 * F1
               H = ZPSS1 * F1 + ZPSS2 * F2
               R = ZPSS1 * F2 + ZPSS2 * F3

               GXYZ = YSPS1*(ZSSP1*A+ZSSP2*B) + YSSP2*(ZSSP1*B+ZSSP2*C)
               GXZY = ZSPS1*(YSSP1*A+YSSP2*B) + ZSSP2*(YSSP1*B+YSSP2*C)
               GYXZ = XSPS1*(ZSSP1*D+ZSSP2*E) + XSSP2*(ZSSP1*E+ZSSP2*F)
               GYZX = ZSPS1*(XSSP1*D+XSSP2*E) + ZSSP2*(XSSP1*E+XSSP2*F)
               GZXY = XSPS1*(YSSP1*G+YSSP2*H) + XSSP2*(YSSP1*H+YSSP2*R)
               GZYX = YSPS1*(XSSP1*G+XSSP2*H) + YSSP2*(XSSP1*H+XSSP2*R)

               BATCH (M+ 1) = GXXX
               BATCH (M+ 2) = GYXX
               BATCH (M+ 3) = GZXX
               BATCH (M+ 4) = GXYX
               BATCH (M+ 5) = GYYX
               BATCH (M+ 6) = GZYX
               BATCH (M+ 7) = GXZX
               BATCH (M+ 8) = GYZX
               BATCH (M+ 9) = GZZX
               BATCH (M+10) = GXXY
               BATCH (M+11) = GYXY
               BATCH (M+12) = GZXY
               BATCH (M+13) = GXYY
               BATCH (M+14) = GYYY
               BATCH (M+15) = GZYY
               BATCH (M+16) = GXZY
               BATCH (M+17) = GYZY
               BATCH (M+18) = GZZY
               BATCH (M+19) = GXXZ
               BATCH (M+20) = GYXZ
               BATCH (M+21) = GZXZ
               BATCH (M+22) = GXYZ
               BATCH (M+23) = GYYZ
               BATCH (M+24) = GZYZ
               BATCH (M+25) = GXZZ
               BATCH (M+26) = GYZZ
               BATCH (M+27) = GZZZ

               M = M + 27

 1550       CONTINUE
 1500    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for atomic ppsp and
C                ppps integrals. All these integrals are zero by
C                symmetry.
C
C
   21    RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 2-center ppsp and
C                ppps (AA|CC) type integrals.
C
C
   22    PQX = X1 - X3
         PQY = Y1 - Y3
         PQZ = Z1 - Z3
         RNPQSQ = PQX * PQX + PQY * PQY + PQZ * PQZ

         M = 0
         DO 2200 IJ = 1,MIJ
            PVAL = P (IJ)
            PINV = ONE / PVAL
            PSCALE = SCALEP (IJ)
            U4 = HALF * PINV
            DO 2220 KL = 1,MKL
               QVAL = Q (KL)
               PQMULT = PVAL * QVAL
               PQPLUS = PVAL + QVAL
               PQPINV = ONE / PQPLUS
               T = RNPQSQ * PQMULT * PQPINV
               SCALE = PSCALE * SCALEQ (KL) / (PQMULT * DSQRT (PQPLUS))

               IF (T.LE.TMAX) THEN
                   TGRID = INT (T * TVSTEP + HALF)
                   DELTA1 = TGRID * TSTEP - T
                   DELTA2 = DELTA1 * HALF
                   DELTA3 = DELTA1 * THIRD
                   DELTA4 = DELTA2 * HALF
                   DELTA5 = DELTA1 * FIFTH
                   DELTA6 = DELTA3 * HALF

                   F1 = ((((( FTABLE (7,TGRID)  * DELTA6
     +                      + FTABLE (6,TGRID)) * DELTA5
     +                      + FTABLE (5,TGRID)) * DELTA4
     +                      + FTABLE (4,TGRID)) * DELTA3
     +                      + FTABLE (3,TGRID)) * DELTA2
     +                      + FTABLE (2,TGRID)) * DELTA1
     +                      + FTABLE (1,TGRID)

                   F2 = ((((( FTABLE (8,TGRID)  * DELTA6
     +                      + FTABLE (7,TGRID)) * DELTA5
     +                      + FTABLE (6,TGRID)) * DELTA4
     +                      + FTABLE (5,TGRID)) * DELTA3
     +                      + FTABLE (4,TGRID)) * DELTA2
     +                      + FTABLE (3,TGRID)) * DELTA1
     +                      + FTABLE (2,TGRID)

                   F3 = ((((( FTABLE (9,TGRID)  * DELTA6
     +                      + FTABLE (8,TGRID)) * DELTA5
     +                      + FTABLE (7,TGRID)) * DELTA4
     +                      + FTABLE (6,TGRID)) * DELTA3
     +                      + FTABLE (5,TGRID)) * DELTA2
     +                      + FTABLE (4,TGRID)) * DELTA1
     +                      + FTABLE (3,TGRID)

                   F1 = SCALE * F1
                   F2 = SCALE * F2
                   F3 = SCALE * F3

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F1 = T2INV * SCALE * HALF * DSQRT (PI*TINV)
                   F2 = THREE * T2INV * F1
                   F3 = FIVE  * T2INV * F2

               END IF

               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               U2 = HALF * PQPINV
               U3 = U2 + PQPINV
               U5 = U1 * U4
C
C
C             ...the X-terms.
C
C
               XSSP2 = PQX * U0
               XSPS2 = PQX * U1
               XSPP3 = XSSP2 * XSPS2
               XPPS3 = XSPS2 * XSPS2
               XPPP2 = U4 * XSSP2
               XPPP3 = U3 * XSPS2
               XPPP4 = XSSP2 * XPPS3
C
C
C             ...the Y-terms.
C
C
               YSSP2 = PQY * U0
               YSPS2 = PQY * U1
               YSPP3 = YSSP2 * YSPS2
               YPPS3 = YSPS2 * YSPS2
               YPPP2 = U4 * YSSP2
               YPPP3 = U3 * YSPS2
               YPPP4 = YSSP2 * YPPS3
C
C
C             ...the Z-terms.
C
C
               ZSSP2 = PQZ * U0
               ZSPS2 = PQZ * U1
               ZSPP3 = ZSSP2 * ZSPS2
               ZPPS3 = ZSPS2 * ZSPS2
               ZPPP2 = U4 * ZSSP2
               ZPPP3 = U3 * ZSPS2
               ZPPP4 = ZSSP2 * ZPPS3
C
C
C             ...assemble the 2-center (AA|CC) type integrals.
C
C
               GXXX = XPPP2 * F1 + XPPP3 * F2 + XPPP4 * F3
               GYYY = YPPP2 * F1 + YPPP3 * F2 + YPPP4 * F3
               GZZZ = ZPPP2 * F1 + ZPPP3 * F2 + ZPPP4 * F3

               A = U4 * (F1 + U1 * F2)
               B = U2 * F2
               C = A + XPPS3 * F3
               D = B + XSPP3 * F3

               GXXY = YSSP2 * C
               GXXZ = ZSSP2 * C
               GXYX = YSPS2 * D
               GXZX = ZSPS2 * D

               C = A + YPPS3 * F3
               D = B + YSPP3 * F3

               GYYX = XSSP2 * C
               GYYZ = ZSSP2 * C
               GXYY = XSPS2 * D
               GYZY = ZSPS2 * D

               C = A + ZPPS3 * F3
               D = B + ZSPP3 * F3

               GZZX = XSSP2 * C
               GZZY = YSSP2 * C
               GXZZ = XSPS2 * D
               GYZZ = YSPS2 * D

               GXYZ = XSPS2 * YSPS2 * ZSSP2 * F3
               GXZY = XSPS2 * ZSPS2 * YSSP2 * F3
               GYZX = YSPS2 * ZSPS2 * XSSP2 * F3

               BATCH (M+ 1) = GXXX
               BATCH (M+ 2) = GXYX
               BATCH (M+ 3) = GXZX
               BATCH (M+ 4) = GXYX
               BATCH (M+ 5) = GYYX
               BATCH (M+ 6) = GYZX
               BATCH (M+ 7) = GXZX
               BATCH (M+ 8) = GYZX
               BATCH (M+ 9) = GZZX
               BATCH (M+10) = GXXY
               BATCH (M+11) = GXYY
               BATCH (M+12) = GXZY
               BATCH (M+13) = GXYY
               BATCH (M+14) = GYYY
               BATCH (M+15) = GYZY
               BATCH (M+16) = GXZY
               BATCH (M+17) = GYZY
               BATCH (M+18) = GZZY
               BATCH (M+19) = GXXZ
               BATCH (M+20) = GXYZ
               BATCH (M+21) = GXZZ
               BATCH (M+22) = GXYZ
               BATCH (M+23) = GYYZ
               BATCH (M+24) = GYZZ
               BATCH (M+25) = GXZZ
               BATCH (M+26) = GYZZ
               BATCH (M+27) = GZZZ

               M = M + 27

 2220       CONTINUE
 2200    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 3-center ppsp and
C                ppps (AB|CC) type integrals.
C
C
   23    M = 0
         DO 2300 IJ = 1,MIJ
            PVAL = P (IJ)
            PINV = ONE / PVAL
            PXVAL = PX (IJ)
            PYVAL = PY (IJ)
            PZVAL = PZ (IJ)
            PQX = PXVAL - X3
            PQY = PYVAL - Y3
            PQZ = PZVAL - Z3
            RNPQSQ = PQX * PQX + PQY * PQY + PQZ * PQZ
            PSCALE = SCALEP (IJ)
            U4 = HALF * PINV
            XSPS1 = PXVAL - X2
            XPSS1 = PXVAL - X1
            YSPS1 = PYVAL - Y2
            YPSS1 = PYVAL - Y1
            ZSPS1 = PZVAL - Z2
            ZPSS1 = PZVAL - Z1
            DO 2330 KL = 1,MKL
               QVAL = Q (KL)
               PQMULT = PVAL * QVAL
               PQPLUS = PVAL + QVAL
               PQPINV = ONE / PQPLUS
               T = RNPQSQ * PQMULT * PQPINV
               SCALE = PSCALE * SCALEQ (KL) / (PQMULT * DSQRT (PQPLUS))

               IF (T.LE.TMAX) THEN
                   TGRID = INT (T * TVSTEP + HALF)
                   DELTA1 = TGRID * TSTEP - T
                   DELTA2 = DELTA1 * HALF
                   DELTA3 = DELTA1 * THIRD
                   DELTA4 = DELTA2 * HALF
                   DELTA5 = DELTA1 * FIFTH
                   DELTA6 = DELTA3 * HALF

                   F1 = ((((( FTABLE (7,TGRID)  * DELTA6
     +                      + FTABLE (6,TGRID)) * DELTA5
     +                      + FTABLE (5,TGRID)) * DELTA4
     +                      + FTABLE (4,TGRID)) * DELTA3
     +                      + FTABLE (3,TGRID)) * DELTA2
     +                      + FTABLE (2,TGRID)) * DELTA1
     +                      + FTABLE (1,TGRID)

                   F2 = ((((( FTABLE (8,TGRID)  * DELTA6
     +                      + FTABLE (7,TGRID)) * DELTA5
     +                      + FTABLE (6,TGRID)) * DELTA4
     +                      + FTABLE (5,TGRID)) * DELTA3
     +                      + FTABLE (4,TGRID)) * DELTA2
     +                      + FTABLE (3,TGRID)) * DELTA1
     +                      + FTABLE (2,TGRID)

                   F3 = ((((( FTABLE (9,TGRID)  * DELTA6
     +                      + FTABLE (8,TGRID)) * DELTA5
     +                      + FTABLE (7,TGRID)) * DELTA4
     +                      + FTABLE (6,TGRID)) * DELTA3
     +                      + FTABLE (5,TGRID)) * DELTA2
     +                      + FTABLE (4,TGRID)) * DELTA1
     +                      + FTABLE (3,TGRID)

                   F1 = SCALE * F1
                   F2 = SCALE * F2
                   F3 = SCALE * F3

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F1 = T2INV * SCALE * HALF * DSQRT (PI*TINV)
                   F2 = THREE * T2INV * F1
                   F3 = FIVE  * T2INV * F2

               END IF

               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               U2 = HALF * PQPINV
               U3 = U2 + PQPINV
               U5 = U1 * U4
C
C
C             ...the X-terms.
C
C
               XSSP2 = PQX * U0
               XSPS2 = PQX * U1
               A     = XPSS1 + XSPS1
               XSPP2 = XSPS1 * XSSP2 + U2
               XSPP3 = XSSP2 * XSPS2
               XPSP2 = XPSS1 * XSSP2 + U2
               XPPS1 = XPSS1 * XSPS1 + U4
               XPPS2 = A * XSPS2 + U5
               XPPS3 = XSPS2 * XSPS2
               XPPP2 = XPPS1 * XSSP2 + A * U2
               XPPP3 = A * XSPP3 + U3 * XSPS2
               XPPP4 = XSSP2 * XPPS3
C
C
C             ...the Y-terms.
C
C
               YSSP2 = PQY * U0
               YSPS2 = PQY * U1
               A     = YPSS1 + YSPS1
               YSPP2 = YSPS1 * YSSP2 + U2
               YSPP3 = YSSP2 * YSPS2
               YPSP2 = YPSS1 * YSSP2 + U2
               YPPS1 = YPSS1 * YSPS1 + U4
               YPPS2 = A * YSPS2 + U5
               YPPS3 = YSPS2 * YSPS2
               YPPP2 = YPPS1 * YSSP2 + A * U2
               YPPP3 = A * YSPP3 + U3 * YSPS2
               YPPP4 = YSSP2 * YPPS3
C
C
C             ...the Z-terms.
C
C
               ZSSP2 = PQZ * U0
               ZSPS2 = PQZ * U1
               A     = ZPSS1 + ZSPS1
               ZSPP2 = ZSPS1 * ZSSP2 + U2
               ZSPP3 = ZSSP2 * ZSPS2
               ZPSP2 = ZPSS1 * ZSSP2 + U2
               ZPPS1 = ZPSS1 * ZSPS1 + U4
               ZPPS2 = A * ZSPS2 + U5
               ZPPS3 = ZSPS2 * ZSPS2
               ZPPP2 = ZPPS1 * ZSSP2 + A * U2
               ZPPP3 = A * ZSPP3 + U3 * ZSPS2
               ZPPP4 = ZSSP2 * ZPPS3
C
C
C             ...assemble the 3-center (AB|CC) type integrals.
C
C
               GXXX = XPPP2 * F1 + XPPP3 * F2 + XPPP4 * F3
               GYYY = YPPP2 * F1 + YPPP3 * F2 + YPPP4 * F3
               GZZZ = ZPPP2 * F1 + ZPPP3 * F2 + ZPPP4 * F3

               AA = XSPP3 * F2
               BB = XSPP3 * F3

               A = XPPS1 * F1 + XPPS2 * F2 + XPPS3 * F3
               B = XPSP2 * F1 + AA
               C = XPSP2 * F2 + BB
               D = XSPP2 * F1 + AA
               E = XSPP2 * F2 + BB

               GXXY = YSSP2 * A
               GXXZ = ZSSP2 * A
               GXYX = YSPS1 * B + YSPS2 * C
               GXZX = ZSPS1 * B + ZSPS2 * C
               GYXX = YPSS1 * D + YSPS2 * E
               GZXX = ZPSS1 * D + ZSPS2 * E

               AA = YSPP3 * F2
               BB = YSPP3 * F3

               A = YPPS1 * F1 + YPPS2 * F2 + YPPS3 * F3
               B = YPSP2 * F1 + AA
               C = YPSP2 * F2 + BB
               D = YSPP2 * F1 + AA
               E = YSPP2 * F2 + BB

               GYYX = XSSP2 * A
               GYYZ = ZSSP2 * A
               GYXY = XSPS1 * B + XSPS2 * C
               GYZY = ZSPS1 * B + ZSPS2 * C
               GXYY = XPSS1 * D + XSPS2 * E
               GZYY = ZPSS1 * D + ZSPS2 * E

               AA = ZSPP3 * F2
               BB = ZSPP3 * F3

               A = ZPPS1 * F1 + ZPPS2 * F2 + ZPPS3 * F3
               B = ZPSP2 * F1 + AA
               C = ZPSP2 * F2 + BB
               D = ZSPP2 * F1 + AA
               E = ZSPP2 * F2 + BB

               GZZX = XSSP2 * A
               GZZY = YSSP2 * A
               GZXZ = XSPS1 * B + XSPS2 * C
               GZYZ = YSPS1 * B + YSPS2 * C
               GXZZ = XPSS1 * D + XSPS2 * E
               GYZZ = YPSS1 * D + YSPS2 * E

               A = XPSS1 * F1 + XSPS2 * F2
               B = XPSS1 * F2 + XSPS2 * F3
               C = YPSS1 * F1 + YSPS2 * F2
               D = YPSS1 * F2 + YSPS2 * F3
               E = ZPSS1 * F1 + ZSPS2 * F2
               F = ZPSS1 * F2 + ZSPS2 * F3

               GXYZ = (YSPS1 * A + YSPS2 * B) * ZSSP2
               GXZY = (ZSPS1 * A + ZSPS2 * B) * YSSP2
               GYXZ = (XSPS1 * C + XSPS2 * D) * ZSSP2
               GYZX = (ZSPS1 * C + ZSPS2 * D) * XSSP2
               GZXY = (XSPS1 * E + XSPS2 * F) * YSSP2
               GZYX = (YSPS1 * E + YSPS2 * F) * XSSP2

               BATCH (M+ 1) = GXXX
               BATCH (M+ 2) = GYXX
               BATCH (M+ 3) = GZXX
               BATCH (M+ 4) = GXYX
               BATCH (M+ 5) = GYYX
               BATCH (M+ 6) = GZYX
               BATCH (M+ 7) = GXZX
               BATCH (M+ 8) = GYZX
               BATCH (M+ 9) = GZZX
               BATCH (M+10) = GXXY
               BATCH (M+11) = GYXY
               BATCH (M+12) = GZXY
               BATCH (M+13) = GXYY
               BATCH (M+14) = GYYY
               BATCH (M+15) = GZYY
               BATCH (M+16) = GXZY
               BATCH (M+17) = GYZY
               BATCH (M+18) = GZZY
               BATCH (M+19) = GXXZ
               BATCH (M+20) = GYXZ
               BATCH (M+21) = GZXZ
               BATCH (M+22) = GXYZ
               BATCH (M+23) = GYYZ
               BATCH (M+24) = GZYZ
               BATCH (M+25) = GXZZ
               BATCH (M+26) = GYZZ
               BATCH (M+27) = GZZZ

               M = M + 27

 2330       CONTINUE
 2300    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 3-center ppsp and
C                ppps (AA|CD) type integrals.
C
C
   24    IF (SHELL3.EQ.1) THEN
             QXSUB = X3
             QYSUB = Y3
             QZSUB = Z3
         ELSE
             QXSUB = X4
             QYSUB = Y4
             QZSUB = Z4
         END IF

         M = 0
         DO 2400 IJ = 1,MIJ
            PVAL = P (IJ)
            PINV = ONE / PVAL
            PSCALE = SCALEP (IJ)
            U4 = HALF * PINV
            DO 2440 KL = 1,MKL
               QVAL = Q (KL)
               QXVAL = QX (KL)
               QYVAL = QY (KL)
               QZVAL = QZ (KL)
               PQMULT = PVAL * QVAL
               PQPLUS = PVAL + QVAL
               PQPINV = ONE / PQPLUS
               PQX = X1 - QXVAL
               PQY = Y1 - QYVAL
               PQZ = Z1 - QZVAL
               T = (PQX*PQX + PQY*PQY + PQZ*PQZ) * PQMULT * PQPINV
               SCALE = PSCALE * SCALEQ (KL) / (PQMULT * DSQRT (PQPLUS))

               IF (T.LE.TMAX) THEN
                   TGRID = INT (T * TVSTEP + HALF)
                   DELTA1 = TGRID * TSTEP - T
                   DELTA2 = DELTA1 * HALF
                   DELTA3 = DELTA1 * THIRD
                   DELTA4 = DELTA2 * HALF
                   DELTA5 = DELTA1 * FIFTH
                   DELTA6 = DELTA3 * HALF

                   F0 = ((((( FTABLE (6,TGRID)  * DELTA6
     +                      + FTABLE (5,TGRID)) * DELTA5
     +                      + FTABLE (4,TGRID)) * DELTA4
     +                      + FTABLE (3,TGRID)) * DELTA3
     +                      + FTABLE (2,TGRID)) * DELTA2
     +                      + FTABLE (1,TGRID)) * DELTA1
     +                      + FTABLE (0,TGRID)

                   F1 = ((((( FTABLE (7,TGRID)  * DELTA6
     +                      + FTABLE (6,TGRID)) * DELTA5
     +                      + FTABLE (5,TGRID)) * DELTA4
     +                      + FTABLE (4,TGRID)) * DELTA3
     +                      + FTABLE (3,TGRID)) * DELTA2
     +                      + FTABLE (2,TGRID)) * DELTA1
     +                      + FTABLE (1,TGRID)

                   F2 = ((((( FTABLE (8,TGRID)  * DELTA6
     +                      + FTABLE (7,TGRID)) * DELTA5
     +                      + FTABLE (6,TGRID)) * DELTA4
     +                      + FTABLE (5,TGRID)) * DELTA3
     +                      + FTABLE (4,TGRID)) * DELTA2
     +                      + FTABLE (3,TGRID)) * DELTA1
     +                      + FTABLE (2,TGRID)

                   F3 = ((((( FTABLE (9,TGRID)  * DELTA6
     +                      + FTABLE (8,TGRID)) * DELTA5
     +                      + FTABLE (7,TGRID)) * DELTA4
     +                      + FTABLE (6,TGRID)) * DELTA3
     +                      + FTABLE (5,TGRID)) * DELTA2
     +                      + FTABLE (4,TGRID)) * DELTA1
     +                      + FTABLE (3,TGRID)

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2
                   F3 = SCALE * F3

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1
                   F3 = FIVE  * T2INV * F2

               END IF

               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               U2 = HALF * PQPINV
               U3 = U2 + PQPINV
               U5 = U1 * U4
C
C
C             ...the X-terms.
C
C
               XSSP1 = QXVAL - QXSUB
               XSSP2 = PQX * U0
               XSPS2 = PQX * U1
               XSPP2 = XSSP1 * XSPS2 + U2
               XSPP3 = XSSP2 * XSPS2
               XPPS3 = XSPS2 * XSPS2
               XPPP1 = XSSP1 * U4
               XPPP2 = XSSP1 * U5 + U4 * XSSP2
               XPPP3 = XSSP1 * XPPS3 + U3 * XSPS2
               XPPP4 = XSSP2 * XPPS3
C
C
C             ...the Y-terms.
C
C
               YSSP1 = QYVAL - QYSUB
               YSSP2 = PQY * U0
               YSPS2 = PQY * U1
               YSPP2 = YSSP1 * YSPS2 + U2
               YSPP3 = YSSP2 * YSPS2
               YPPS3 = YSPS2 * YSPS2
               YPPP1 = YSSP1 * U4
               YPPP2 = YSSP1 * U5 + U4 * YSSP2
               YPPP3 = YSSP1 * YPPS3 + U3 * YSPS2
               YPPP4 = YSSP2 * YPPS3
C
C
C             ...the Z-terms.
C
C
               ZSSP1 = QZVAL - QZSUB
               ZSSP2 = PQZ * U0
               ZSPS2 = PQZ * U1
               ZSPP2 = ZSSP1 * ZSPS2 + U2
               ZSPP3 = ZSSP2 * ZSPS2
               ZPPS3 = ZSPS2 * ZSPS2
               ZPPP1 = ZSSP1 * U4
               ZPPP2 = ZSSP1 * U5 + U4 * ZSSP2
               ZPPP3 = ZSSP1 * ZPPS3 + U3 * ZSPS2
               ZPPP4 = ZSSP2 * ZPPS3
C
C
C             ...assemble the 3-center (AA|CD) type integrals.
C
C
               GXXX = XPPP1 * F0 + XPPP2 * F1 + XPPP3 * F2 + XPPP4 * F3
               GYYY = YPPP1 * F0 + YPPP2 * F1 + YPPP3 * F2 + YPPP4 * F3
               GZZZ = ZPPP1 * F0 + ZPPP2 * F1 + ZPPP3 * F2 + ZPPP4 * F3

               A = U4 * F0 + U5 * F1
               B = U4 * F1 + U5 * F2
               C = A + XPPS3 * F2
               D = B + XPPS3 * F3
               E = XSPP2 * F2 + XSPP3 * F3

               GXXY = YSSP1 * C + YSSP2 * D
               GXXZ = ZSSP1 * C + ZSSP2 * D
               GXYX = YSPS2 * E
               GXZX = ZSPS2 * E

               C = A + YPPS3 * F2
               D = B + YPPS3 * F3
               E = YSPP2 * F2 + YSPP3 * F3

               GYYX = XSSP1 * C + XSSP2 * D
               GYYZ = ZSSP1 * C + ZSSP2 * D
               GXYY = XSPS2 * E
               GYZY = ZSPS2 * E

               C = A + ZPPS3 * F2
               D = B + ZPPS3 * F3
               E = ZSPP2 * F2 + ZSPP3 * F3

               GZZX = XSSP1 * C + XSSP2 * D
               GZZY = YSSP1 * C + YSSP2 * D
               GXZZ = XSPS2 * E
               GYZZ = YSPS2 * E

               GXYZ = XSPS2 * YSPS2 * (ZSSP1 * F2 + ZSSP2 * F3)
               GXZY = XSPS2 * ZSPS2 * (YSSP1 * F2 + YSSP2 * F3)
               GYZX = YSPS2 * ZSPS2 * (XSSP1 * F2 + XSSP2 * F3)

               BATCH (M+ 1) = GXXX
               BATCH (M+ 2) = GXYX
               BATCH (M+ 3) = GXZX
               BATCH (M+ 4) = GXYX
               BATCH (M+ 5) = GYYX
               BATCH (M+ 6) = GYZX
               BATCH (M+ 7) = GXZX
               BATCH (M+ 8) = GYZX
               BATCH (M+ 9) = GZZX
               BATCH (M+10) = GXXY
               BATCH (M+11) = GXYY
               BATCH (M+12) = GXZY
               BATCH (M+13) = GXYY
               BATCH (M+14) = GYYY
               BATCH (M+15) = GYZY
               BATCH (M+16) = GXZY
               BATCH (M+17) = GYZY
               BATCH (M+18) = GZZY
               BATCH (M+19) = GXXZ
               BATCH (M+20) = GXYZ
               BATCH (M+21) = GXZZ
               BATCH (M+22) = GXYZ
               BATCH (M+23) = GYYZ
               BATCH (M+24) = GYZZ
               BATCH (M+25) = GXZZ
               BATCH (M+26) = GYZZ
               BATCH (M+27) = GZZZ

               M = M + 27

 2440       CONTINUE
 2400    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 4-center ppsp and
C                ppps (AB|CD) type integrals.
C
C
   25    IF (SHELL3.EQ.1) THEN
             QXSUB = X3
             QYSUB = Y3
             QZSUB = Z3
         ELSE
             QXSUB = X4
             QYSUB = Y4
             QZSUB = Z4
         END IF

         M = 0
         DO 2500 IJ = 1,MIJ
            PVAL = P (IJ)
            PINV = ONE / PVAL
            PXVAL = PX (IJ)
            PYVAL = PY (IJ)
            PZVAL = PZ (IJ)
            PSCALE = SCALEP (IJ)
            U4 = HALF * PINV
            XSPS1 = PXVAL - X2
            XPSS1 = PXVAL - X1
            YSPS1 = PYVAL - Y2
            YPSS1 = PYVAL - Y1
            ZSPS1 = PZVAL - Z2
            ZPSS1 = PZVAL - Z1
            DO 2550 KL = 1,MKL
               QVAL = Q (KL)
               QXVAL = QX (KL)
               QYVAL = QY (KL)
               QZVAL = QZ (KL)
               PQMULT = PVAL * QVAL
               PQPLUS = PVAL + QVAL
               PQPINV = ONE / PQPLUS
               PQX = PXVAL - QXVAL
               PQY = PYVAL - QYVAL
               PQZ = PZVAL - QZVAL
               T = (PQX*PQX + PQY*PQY + PQZ*PQZ) * PQMULT * PQPINV
               SCALE = PSCALE * SCALEQ (KL) / (PQMULT * DSQRT (PQPLUS))

               IF (T.LE.TMAX) THEN
                   TGRID = INT (T * TVSTEP + HALF)
                   DELTA1 = TGRID * TSTEP - T
                   DELTA2 = DELTA1 * HALF
                   DELTA3 = DELTA1 * THIRD
                   DELTA4 = DELTA2 * HALF
                   DELTA5 = DELTA1 * FIFTH
                   DELTA6 = DELTA3 * HALF

                   F0 = ((((( FTABLE (6,TGRID)  * DELTA6
     +                      + FTABLE (5,TGRID)) * DELTA5
     +                      + FTABLE (4,TGRID)) * DELTA4
     +                      + FTABLE (3,TGRID)) * DELTA3
     +                      + FTABLE (2,TGRID)) * DELTA2
     +                      + FTABLE (1,TGRID)) * DELTA1
     +                      + FTABLE (0,TGRID)

                   F1 = ((((( FTABLE (7,TGRID)  * DELTA6
     +                      + FTABLE (6,TGRID)) * DELTA5
     +                      + FTABLE (5,TGRID)) * DELTA4
     +                      + FTABLE (4,TGRID)) * DELTA3
     +                      + FTABLE (3,TGRID)) * DELTA2
     +                      + FTABLE (2,TGRID)) * DELTA1
     +                      + FTABLE (1,TGRID)

                   F2 = ((((( FTABLE (8,TGRID)  * DELTA6
     +                      + FTABLE (7,TGRID)) * DELTA5
     +                      + FTABLE (6,TGRID)) * DELTA4
     +                      + FTABLE (5,TGRID)) * DELTA3
     +                      + FTABLE (4,TGRID)) * DELTA2
     +                      + FTABLE (3,TGRID)) * DELTA1
     +                      + FTABLE (2,TGRID)

                   F3 = ((((( FTABLE (9,TGRID)  * DELTA6
     +                      + FTABLE (8,TGRID)) * DELTA5
     +                      + FTABLE (7,TGRID)) * DELTA4
     +                      + FTABLE (6,TGRID)) * DELTA3
     +                      + FTABLE (5,TGRID)) * DELTA2
     +                      + FTABLE (4,TGRID)) * DELTA1
     +                      + FTABLE (3,TGRID)

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2
                   F3 = SCALE * F3

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1
                   F3 = FIVE  * T2INV * F2

               END IF

               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               U2 = HALF * PQPINV
               U3 = U2 + PQPINV
               U5 = U1 * U4
C
C
C             ...the X-terms.
C
C
               XSSP1 = QXVAL - QXSUB
               XSSP2 = PQX * U0
               XSPS2 = PQX * U1
               A     = XPSS1 + XSPS1
               B     = XSSP1 * XSPS2 + U2
               XSPP1 = XSPS1 * XSSP1
               XSPP2 = XSPS1 * XSSP2 + B
               XSPP3 = XSSP2 * XSPS2
               XPSP1 = XPSS1 * XSSP1
               XPSP2 = XPSS1 * XSSP2 + B
               XPPS1 = XPSS1 * XSPS1 + U4
               XPPS2 = A * XSPS2 + U5
               XPPS3 = XSPS2 * XSPS2
               XPPP1 = XSSP1 * XPPS1
               XPPP2 = XSSP1 * XPPS2 + XPPS1 * XSSP2 + A * U2
               XPPP3 = XSSP1 * XPPS3 + A * XSPP3 + U3 * XSPS2
               XPPP4 = XSSP2 * XPPS3
C
C
C             ...the Y-terms.
C
C
               YSSP1 = QYVAL - QYSUB
               YSSP2 = PQY * U0
               YSPS2 = PQY * U1
               A     = YPSS1 + YSPS1
               B     = YSSP1 * YSPS2 + U2
               YSPP1 = YSPS1 * YSSP1
               YSPP2 = YSPS1 * YSSP2 + B
               YSPP3 = YSSP2 * YSPS2
               YPSP1 = YPSS1 * YSSP1
               YPSP2 = YPSS1 * YSSP2 + B
               YPPS1 = YPSS1 * YSPS1 + U4
               YPPS2 = A * YSPS2 + U5
               YPPS3 = YSPS2 * YSPS2
               YPPP1 = YSSP1 * YPPS1
               YPPP2 = YSSP1 * YPPS2 + YPPS1 * YSSP2 + A * U2
               YPPP3 = YSSP1 * YPPS3 + A * YSPP3 + U3 * YSPS2
               YPPP4 = YSSP2 * YPPS3
C
C
C             ...the Z-terms.
C
C
               ZSSP1 = QZVAL - QZSUB
               ZSSP2 = PQZ * U0
               ZSPS2 = PQZ * U1
               A     = ZPSS1 + ZSPS1
               B     = ZSSP1 * ZSPS2 + U2
               ZSPP1 = ZSPS1 * ZSSP1
               ZSPP2 = ZSPS1 * ZSSP2 + B
               ZSPP3 = ZSSP2 * ZSPS2
               ZPSP1 = ZPSS1 * ZSSP1
               ZPSP2 = ZPSS1 * ZSSP2 + B
               ZPPS1 = ZPSS1 * ZSPS1 + U4
               ZPPS2 = A * ZSPS2 + U5
               ZPPS3 = ZSPS2 * ZSPS2
               ZPPP1 = ZSSP1 * ZPPS1
               ZPPP2 = ZSSP1 * ZPPS2 + ZPPS1 * ZSSP2 + A * U2
               ZPPP3 = ZSSP1 * ZPPS3 + A * ZSPP3 + U3 * ZSPS2
               ZPPP4 = ZSSP2 * ZPPS3
C
C
C             ...assemble the 4-center (AB|CD) type integrals.
C
C
               GXXX = XPPP1 * F0 + XPPP2 * F1 + XPPP3 * F2 + XPPP4 * F3
               GYYY = YPPP1 * F0 + YPPP2 * F1 + YPPP3 * F2 + YPPP4 * F3
               GZZZ = ZPPP1 * F0 + ZPPP2 * F1 + ZPPP3 * F2 + ZPPP4 * F3

               AA = XSPP3 * F2
               BB = XSPP3 * F3

               A = XPPS1 * F0 + XPPS2 * F1 + XPPS3 * F2
               B = XPPS1 * F1 + XPPS2 * F2 + XPPS3 * F3
               C = XPSP1 * F0 + XPSP2 * F1 + AA
               D = XPSP1 * F1 + XPSP2 * F2 + BB
               E = XSPP1 * F0 + XSPP2 * F1 + AA
               F = XSPP1 * F1 + XSPP2 * F2 + BB

               GXXY = YSSP1 * A + YSSP2 * B
               GXXZ = ZSSP1 * A + ZSSP2 * B
               GXYX = YSPS1 * C + YSPS2 * D
               GXZX = ZSPS1 * C + ZSPS2 * D
               GYXX = YPSS1 * E + YSPS2 * F
               GZXX = ZPSS1 * E + ZSPS2 * F

               AA = YSPP3 * F2
               BB = YSPP3 * F3

               A = YPPS1 * F0 + YPPS2 * F1 + YPPS3 * F2
               B = YPPS1 * F1 + YPPS2 * F2 + YPPS3 * F3
               C = YPSP1 * F0 + YPSP2 * F1 + AA
               D = YPSP1 * F1 + YPSP2 * F2 + BB
               E = YSPP1 * F0 + YSPP2 * F1 + AA
               F = YSPP1 * F1 + YSPP2 * F2 + BB

               GYYX = XSSP1 * A + XSSP2 * B
               GYYZ = ZSSP1 * A + ZSSP2 * B
               GYXY = XSPS1 * C + XSPS2 * D
               GYZY = ZSPS1 * C + ZSPS2 * D
               GXYY = XPSS1 * E + XSPS2 * F
               GZYY = ZPSS1 * E + ZSPS2 * F

               AA = ZSPP3 * F2
               BB = ZSPP3 * F3

               A = ZPPS1 * F0 + ZPPS2 * F1 + ZPPS3 * F2
               B = ZPPS1 * F1 + ZPPS2 * F2 + ZPPS3 * F3
               C = ZPSP1 * F0 + ZPSP2 * F1 + AA
               D = ZPSP1 * F1 + ZPSP2 * F2 + BB
               E = ZSPP1 * F0 + ZSPP2 * F1 + AA
               F = ZSPP1 * F1 + ZSPP2 * F2 + BB

               GZZX = XSSP1 * A + XSSP2 * B
               GZZY = YSSP1 * A + YSSP2 * B
               GZXZ = XSPS1 * C + XSPS2 * D
               GZYZ = YSPS1 * C + YSPS2 * D
               GXZZ = XPSS1 * E + XSPS2 * F
               GYZZ = YPSS1 * E + YSPS2 * F

               A = XPSS1 * F0 + XSPS2 * F1
               B = XPSS1 * F1 + XSPS2 * F2
               C = XPSS1 * F2 + XSPS2 * F3
               D = YPSS1 * F0 + YSPS2 * F1
               E = YPSS1 * F1 + YSPS2 * F2
               F = YPSS1 * F2 + YSPS2 * F3
               G = ZPSS1 * F0 + ZSPS2 * F1
               H = ZPSS1 * F1 + ZSPS2 * F2
               R = ZPSS1 * F2 + ZSPS2 * F3

               GXYZ = YSPS1*(ZSSP1*A+ZSSP2*B) + YSPS2*(ZSSP1*B+ZSSP2*C)
               GXZY = ZSPS1*(YSSP1*A+YSSP2*B) + ZSPS2*(YSSP1*B+YSSP2*C)
               GYXZ = XSPS1*(ZSSP1*D+ZSSP2*E) + XSPS2*(ZSSP1*E+ZSSP2*F)
               GYZX = ZSPS1*(XSSP1*D+XSSP2*E) + ZSPS2*(XSSP1*E+XSSP2*F)
               GZXY = XSPS1*(YSSP1*G+YSSP2*H) + XSPS2*(YSSP1*H+YSSP2*R)
               GZYX = YSPS1*(XSSP1*G+XSSP2*H) + YSPS2*(XSSP1*H+XSSP2*R)

               BATCH (M+ 1) = GXXX
               BATCH (M+ 2) = GYXX
               BATCH (M+ 3) = GZXX
               BATCH (M+ 4) = GXYX
               BATCH (M+ 5) = GYYX
               BATCH (M+ 6) = GZYX
               BATCH (M+ 7) = GXZX
               BATCH (M+ 8) = GYZX
               BATCH (M+ 9) = GZZX
               BATCH (M+10) = GXXY
               BATCH (M+11) = GYXY
               BATCH (M+12) = GZXY
               BATCH (M+13) = GXYY
               BATCH (M+14) = GYYY
               BATCH (M+15) = GZYY
               BATCH (M+16) = GXZY
               BATCH (M+17) = GYZY
               BATCH (M+18) = GZZY
               BATCH (M+19) = GXXZ
               BATCH (M+20) = GYXZ
               BATCH (M+21) = GZXZ
               BATCH (M+22) = GXYZ
               BATCH (M+23) = GYYZ
               BATCH (M+24) = GZYZ
               BATCH (M+25) = GXZZ
               BATCH (M+26) = GYZZ
               BATCH (M+27) = GZZZ

               M = M + 27

 2550       CONTINUE
 2500    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
