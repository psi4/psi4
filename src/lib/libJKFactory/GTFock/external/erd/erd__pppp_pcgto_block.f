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
         SUBROUTINE  ERD__PPPP_PCGTO_BLOCK
     +
     +                    ( NBATCH,
     +                      ATOMIC,ATOM12,ATOM34,
     +                      MIJ,MKL,
     +                      NIJ,NIJBEG,NIJEND,NKL,NKLBEG,NKLEND,
     +                      NPGTO1,NPGTO2,NPGTO3,NPGTO4,
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
C  OPERATION   : ERD__PPPP_PCGTO_BLOCK
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation is designed to provide ultrafast block
C                evaluation of a batch of normalized electron repulsion
C                integrals between p-shell primitive spherical gaussian
C                type orbitals.
C
C                A batch is defined here as containing all possible
C                integrals, that is its dimension is determined by
C                the total number of primitive functions (here = 81)
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
C                                         pppp
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
C                                    pppp integral batch
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
C                                    cartesian pppp integrals
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

         INTEGER     PRIM1 (1:MIJ)
         INTEGER     PRIM2 (1:MIJ)
         INTEGER     PRIM3 (1:MKL)
         INTEGER     PRIM4 (1:MKL)

         DOUBLE PRECISION  A,B,C,D,E,F,G,H,R
         DOUBLE PRECISION  AA,BB,CC,DD,EE,FF,GG,HH,RR
         DOUBLE PRECISION  DELTA1,DELTA2,DELTA3,DELTA4,DELTA5,DELTA6
         DOUBLE PRECISION  EXP1,EXP2,EXP3,EXP4
         DOUBLE PRECISION  F0,F1,F2,F3,F4
         DOUBLE PRECISION  GXXXX,GXXXY,GXXXZ,GXXYX,GXXYY,GXXYZ,
     +                     GXXZX,GXXZY,GXXZZ,GXYXX,GXYXY,GXYXZ,
     +                     GXYYX,GXYYY,GXYYZ,GXYZX,GXYZY,GXYZZ,
     +                     GXZXX,GXZXY,GXZXZ,GXZYX,GXZYY,GXZYZ,
     +                     GXZZX,GXZZY,GXZZZ,GYXXX,GYXXY,GYXXZ,
     +                     GYXYX,GYXYY,GYXYZ,GYXZX,GYXZY,GYXZZ,
     +                     GYYXX,GYYXY,GYYXZ,GYYYX,GYYYY,GYYYZ,
     +                     GYYZX,GYYZY,GYYZZ,GYZXX,GYZXY,GYZXZ,
     +                     GYZYX,GYZYY,GYZYZ,GYZZX,GYZZY,GYZZZ,
     +                     GZXXX,GZXXY,GZXXZ,GZXYX,GZXYY,GZXYZ,
     +                     GZXZX,GZXZY,GZXZZ,GZYXX,GZYXY,GZYXZ,
     +                     GZYYX,GZYYY,GZYYZ,GZYZX,GZYZY,GZYZZ,
     +                     GZZXX,GZZXY,GZZXZ,GZZYX,GZZYY,GZZYZ,
     +                     GZZZX,GZZZY,GZZZZ
         DOUBLE PRECISION  PVAL,QVAL
         DOUBLE PRECISION  PQPLUS,PQMULT,PQPINV
         DOUBLE PRECISION  PQX,PQY,PQZ
         DOUBLE PRECISION  PSCALE
         DOUBLE PRECISION  PXVAL,PYVAL,PZVAL,QXVAL,QYVAL,QZVAL
         DOUBLE PRECISION  RNPQSQ
         DOUBLE PRECISION  SCALE
         DOUBLE PRECISION  T,TINV,T2INV,TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION  U0,U1,U2,U3,U4,U5,U6,U7,U8,U9,U10
         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4
         DOUBLE PRECISION  X12,Y12,Z12,X34,Y34,Z34
         DOUBLE PRECISION  XSSSP1,XSSPS1,XSPSS1,XPSSS1,XSSPP1,XSPSP1,
     +                     XPSSP1,XSPPS1,XPSPS1,XPPSS1,XSPPP1,XPSPP1,
     +                     XPPSP1,XPPPS1,XPPPP1,XSSSP2,XSPSS2,XSSPP2,
     +                     XSPSP2,XPSSP2,XSPPS2,XPSPS2,XPPSS2,XSPPP2,
     +                     XPSPP2,XPPSP2,XPPPS2,XPPPP2,XSSPP3,XSPSP3,
     +                     XPPSS3,XSPPP3,XPSPP3,XPPSP3,XPPPS3,XPPPP3,
     +                     XSPPP4,XPPSP4,XPPPP4,XPPPP5,
     +                     YSSSP1,YSSPS1,YSPSS1,YPSSS1,YSSPP1,YSPSP1,
     +                     YPSSP1,YSPPS1,YPSPS1,YPPSS1,YSPPP1,YPSPP1,
     +                     YPPSP1,YPPPS1,YPPPP1,YSSSP2,YSPSS2,YSSPP2,
     +                     YSPSP2,YPSSP2,YSPPS2,YPSPS2,YPPSS2,YSPPP2,
     +                     YPSPP2,YPPSP2,YPPPS2,YPPPP2,YSSPP3,YSPSP3,
     +                     YPPSS3,YSPPP3,YPSPP3,YPPSP3,YPPPS3,YPPPP3,
     +                     YSPPP4,YPPSP4,YPPPP4,YPPPP5,
     +                     ZSSSP1,ZSSPS1,ZSPSS1,ZPSSS1,ZSSPP1,ZSPSP1,
     +                     ZPSSP1,ZSPPS1,ZPSPS1,ZPPSS1,ZSPPP1,ZPSPP1,
     +                     ZPPSP1,ZPPPS1,ZPPPP1,ZSSSP2,ZSPSS2,ZSSPP2,
     +                     ZSPSP2,ZPSSP2,ZSPPS2,ZPSPS2,ZPPSS2,ZSPPP2,
     +                     ZPSPP2,ZPPSP2,ZPPPS2,ZPPPP2,ZSSPP3,ZSPSP3,
     +                     ZPPSS3,ZSPPP3,ZPSPP3,ZPPSP3,ZPPPS3,ZPPPP3,
     +                     ZSPPP4,ZPPSP4,ZPPPP4,ZPPPP5
         DOUBLE PRECISION  ZERO,FIFTH,THIRD,HALF,ONE,THREE,PI,FIVE,
     +                     SIX,SEVEN

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

         PARAMETER  (ZERO    =  0.D0 )
         PARAMETER  (FIFTH   =  0.2D0)
         PARAMETER  (THIRD   =  0.333333333333333D0)
         PARAMETER  (HALF    =  0.5D0)
         PARAMETER  (ONE     =  1.D0 )
         PARAMETER  (THREE   =  3.D0 )
         PARAMETER  (PI      =  3.141592653589793D0)
         PARAMETER  (FIVE    =  5.D0 )
         PARAMETER  (SIX     =  6.D0 )
         PARAMETER  (SEVEN   =  7.D0 )
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
C             ...jump according to type of 'K4' loop:
C
C                       CASE  |  Integral center type
C                     --------|-----------------------
C                         1   |     (AA|AA)  atomic
C                         2   |     (AA|CC)  2-center
C                         3   |     (AB|CC)  3-center
C                         4   |     (AA|CD)  3-center
C                         5   |     (AB|CD)  4-center
C
C
         GOTO (10,20,30,40,50) CASE
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for atomic (AA|AA)
C                type integrals.
C
C
   10    M = 0
         DO 1000 IJ = 1,MIJ
            PVAL = P (IJ)
            PSCALE = SCALEP (IJ)
            DO 1100 KL = 1,MKL
               QVAL = Q (KL)
               PQMULT = PVAL * QVAL
               PQPLUS = PVAL + QVAL
               PQPINV = ONE / PQPLUS
               F0 = PSCALE * SCALEQ (KL) / (PQMULT * DSQRT (PQPLUS))
               U0 = HALF * PQPINV
C
C
C             ...assemble the non-zero atomic (AA|AA) type integrals.
C
C
               A = U0 * F0 * FIFTH
               B = U0 * A

               GXYXY = B
               GXXYY = GXYXY + F0 / (SIX * PQMULT)
               GXXXX = GXXYY + A * PQPINV

               BATCH (M+ 1) =  GXXXX
               BATCH (M+ 2) =  ZERO
               BATCH (M+ 3) =  ZERO
               BATCH (M+ 4) =  ZERO
               BATCH (M+ 5) =  GXXYY
               BATCH (M+ 6) =  ZERO
               BATCH (M+ 7) =  ZERO
               BATCH (M+ 8) =  ZERO
               BATCH (M+ 9) =  GXXYY
               BATCH (M+10) =  ZERO
               BATCH (M+11) =  GXYXY
               BATCH (M+12) =  ZERO
               BATCH (M+13) =  GXYXY
               BATCH (M+14) =  ZERO
               BATCH (M+15) =  ZERO
               BATCH (M+16) =  ZERO
               BATCH (M+17) =  ZERO
               BATCH (M+18) =  ZERO
               BATCH (M+19) =  ZERO
               BATCH (M+20) =  ZERO
               BATCH (M+21) =  GXYXY
               BATCH (M+22) =  ZERO
               BATCH (M+23) =  ZERO
               BATCH (M+24) =  ZERO
               BATCH (M+25) =  GXYXY
               BATCH (M+26) =  ZERO
               BATCH (M+27) =  ZERO
               BATCH (M+28) =  ZERO
               BATCH (M+29) =  GXYXY
               BATCH (M+30) =  ZERO
               BATCH (M+31) =  GXYXY
               BATCH (M+32) =  ZERO
               BATCH (M+33) =  ZERO
               BATCH (M+34) =  ZERO
               BATCH (M+35) =  ZERO
               BATCH (M+36) =  ZERO
               BATCH (M+37) =  GXXYY
               BATCH (M+38) =  ZERO
               BATCH (M+39) =  ZERO
               BATCH (M+40) =  ZERO
               BATCH (M+41) =  GXXXX
               BATCH (M+42) =  ZERO
               BATCH (M+43) =  ZERO
               BATCH (M+44) =  ZERO
               BATCH (M+45) =  GXXYY
               BATCH (M+46) =  ZERO
               BATCH (M+47) =  ZERO
               BATCH (M+48) =  ZERO
               BATCH (M+49) =  ZERO
               BATCH (M+50) =  ZERO
               BATCH (M+51) =  GXYXY
               BATCH (M+52) =  ZERO
               BATCH (M+53) =  GXYXY
               BATCH (M+54) =  ZERO
               BATCH (M+55) =  ZERO
               BATCH (M+56) =  ZERO
               BATCH (M+57) =  GXYXY
               BATCH (M+58) =  ZERO
               BATCH (M+59) =  ZERO
               BATCH (M+60) =  ZERO
               BATCH (M+61) =  GXYXY
               BATCH (M+62) =  ZERO
               BATCH (M+63) =  ZERO
               BATCH (M+64) =  ZERO
               BATCH (M+65) =  ZERO
               BATCH (M+66) =  ZERO
               BATCH (M+67) =  ZERO
               BATCH (M+68) =  ZERO
               BATCH (M+69) =  GXYXY
               BATCH (M+70) =  ZERO
               BATCH (M+71) =  GXYXY
               BATCH (M+72) =  ZERO
               BATCH (M+73) =  GXXYY
               BATCH (M+74) =  ZERO
               BATCH (M+75) =  ZERO
               BATCH (M+76) =  ZERO
               BATCH (M+77) =  GXXYY
               BATCH (M+78) =  ZERO
               BATCH (M+79) =  ZERO
               BATCH (M+80) =  ZERO
               BATCH (M+81) =  GXXXX

               M = M + 81

 1100       CONTINUE
 1000    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 2-center (AA|CC)
C                type integrals.
C
C
   20    PQX = X1 - X3
         PQY = Y1 - Y3
         PQZ = Z1 - Z3
         RNPQSQ = PQX*PQX + PQY*PQY + PQZ*PQZ

         M = 0
         DO 2000 IJ = 1,MIJ
            PVAL = P (IJ)
            PSCALE = SCALEP (IJ)
            U3 = HALF / PVAL
            DO 2200 KL = 1,MKL
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

                   F4 = ((((( FTABLE (10,TGRID) * DELTA6
     +                      + FTABLE (9,TGRID)) * DELTA5
     +                      + FTABLE (8,TGRID)) * DELTA4
     +                      + FTABLE (7,TGRID)) * DELTA3
     +                      + FTABLE (6,TGRID)) * DELTA2
     +                      + FTABLE (5,TGRID)) * DELTA1
     +                      + FTABLE (4,TGRID)

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2
                   F3 = SCALE * F3
                   F4 = SCALE * F4

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1
                   F3 = FIVE  * T2INV * F2
                   F4 = SEVEN * T2INV * F3

               END IF

               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               U2 = HALF * PQPINV
               U4 = HALF / QVAL
               U5 = U2 + PQPINV
               U6 = U5 + U5
               U7 = - U0 * U4
               U8 = U1 * U3
               U9 = U3 * U4
               U10 = U2 * U5
C
C
C             ...the X-terms.
C
C
               XSSSP2 = PQX * U0
               XSPSS2 = PQX * U1
               XSPSP3 = XSSSP2 * XSPSS2
               XSSPP3 = XSSSP2 * XSSSP2
               XPPSS3 = XSPSS2 * XSPSS2
               XSPPP2 = XSPSS2 * U4
               XPPSP2 = XSSSP2 * U3
               XSPPP3 = XSSSP2 * U5
               XPPSP3 = XSPSS2 * U5
               XSPPP4 = XSSPP3 * XSPSS2
               XPPSP4 = XPPSS3 * XSSSP2
               XPPPP3 = U3 * XSSPP3 + U4 * XPPSS3 + U10
               XPPPP4 = XSPSP3 * U6
               XPPPP5 = XSSPP3 * XPPSS3
C
C
C             ...the Y-terms.
C
C
               YSSSP2 = PQY * U0
               YSPSS2 = PQY * U1
               YSPSP3 = YSSSP2 * YSPSS2
               YSSPP3 = YSSSP2 * YSSSP2
               YPPSS3 = YSPSS2 * YSPSS2
               YSPPP2 = YSPSS2 * U4
               YPPSP2 = YSSSP2 * U3
               YSPPP3 = YSSSP2 * U5
               YPPSP3 = YSPSS2 * U5
               YSPPP4 = YSSPP3 * YSPSS2
               YPPSP4 = YPPSS3 * YSSSP2
               YPPPP3 = U3 * YSSPP3 + U4 * YPPSS3 + U10
               YPPPP4 = YSPSP3 * U6
               YPPPP5 = YSSPP3 * YPPSS3
C
C
C             ...the Z-terms.
C
C
               ZSSSP2 = PQZ * U0
               ZSPSS2 = PQZ * U1
               ZSPSP3 = ZSSSP2 * ZSPSS2
               ZSSPP3 = ZSSSP2 * ZSSSP2
               ZPPSS3 = ZSPSS2 * ZSPSS2
               ZSPPP2 = ZSPSS2 * U4
               ZPPSP2 = ZSSSP2 * U3
               ZSPPP3 = ZSSSP2 * U5
               ZPPSP3 = ZSPSS2 * U5
               ZSPPP4 = ZSSPP3 * ZSPSS2
               ZPPSP4 = ZPPSS3 * ZSSSP2
               ZPPPP3 = U3 * ZSSPP3 + U4 * ZPPSS3 + U10
               ZPPPP4 = ZSPSP3 * U6
               ZPPPP5 = ZSSPP3 * ZPPSS3
C
C
C             ...assemble the 2-center (AA|CC) type integrals.
C
C
               A = U9 * (F0 + (U1 - U0) * F1)

               GXXXX = A + XPPPP3 * F2 + XPPPP4 * F3 + XPPPP5 * F4
               GYYYY = A + YPPPP3 * F2 + YPPPP4 * F3 + YPPPP5 * F4
               GZZZZ = A + ZPPPP3 * F2 + ZPPPP4 * F3 + ZPPPP5 * F4

               A = XPPSP2 * F2 + XPPSP3 * F3 + XPPSP4 * F4
               B = XSPPP2 * F2 + XSPPP3 * F3 + XSPPP4 * F4
               C = YPPSP2 * F2 + YPPSP3 * F3 + YPPSP4 * F4
               D = YSPPP2 * F2 + YSPPP3 * F3 + YSPPP4 * F4
               E = ZPPSP2 * F2 + ZPPSP3 * F3 + ZPPSP4 * F4
               F = ZSPPP2 * F2 + ZSPPP3 * F3 + ZSPPP4 * F4

               GXXXY = A * YSSSP2
               GXYXX = B * YSPSS2
               GXXXZ = A * ZSSSP2
               GXZXX = B * ZSPSS2
               GYYXY = C * XSSSP2
               GXYYY = D * XSPSS2
               GYYYZ = C * ZSSSP2
               GYZYY = D * ZSPSS2
               GZZXZ = E * XSSSP2
               GXZZZ = F * XSPSS2
               GZZYZ = E * YSSSP2
               GYZZZ = F * YSPSS2

               AA = U3 * F2
               BB = U8 * F3
               CC = (U4 * F0 + U7 * F1) * U3
               DD = (U4 * F1 + U7 * F2) * U8
               EE = (U4 * F2 + U7 * F3)
               FF = AA + BB

               A = CC + XSSPP3 * AA
               B = DD + XSSPP3 * BB
               C = EE + XSSPP3 * F4
               D = CC + YSSPP3 * AA
               E = DD + YSSPP3 * BB
               F = EE + YSSPP3 * F4
               G = CC + ZSSPP3 * AA
               H = DD + ZSSPP3 * BB
               R = EE + ZSSPP3 * F4

               GXXYY = D + E + F * XPPSS3
               GXXZZ = G + H + R * XPPSS3
               GYYXX = A + B + C * YPPSS3
               GYYZZ = G + H + R * YPPSS3
               GZZXX = A + B + C * ZPPSS3
               GZZYY = D + E + F * ZPPSS3
               GYZXX = C * YSPSS2 * ZSPSS2
               GXZYY = F * XSPSS2 * ZSPSS2
               GXYZZ = R * XSPSS2 * YSPSS2
               GXXYZ = (FF + XPPSS3 * F4) * YSSSP2 * ZSSSP2
               GYYXZ = (FF + YPPSS3 * F4) * XSSSP2 * ZSSSP2
               GZZXY = (FF + ZPPSS3 * F4) * XSSSP2 * YSSSP2

               AA = U2 * F2
               BB = U2 * F3
               CC = AA * U2

               A = CC + XSPSP3 * BB
               B = BB + XSPSP3 * F4
               C = CC + YSPSP3 * BB
               D = BB + YSPSP3 * F4
               E = AA + ZSPSP3 * F3
               F = BB + ZSPSP3 * F4

               GXYXY = A + B * YSPSP3
               GXZXZ = A + B * ZSPSP3
               GYZYZ = C + D * ZSPSP3
               GXYXZ = B * YSPSS2 * ZSSSP2
               GXYYZ = D * XSPSS2 * ZSSSP2
               GXZXY = B * ZSPSS2 * YSSSP2
               GXZYZ = F * XSPSS2 * YSSSP2
               GYZXY = D * ZSPSS2 * XSSSP2
               GYZXZ = F * YSPSS2 * XSSSP2

               BATCH (M+ 1) =  GXXXX
               BATCH (M+ 2) =  GXYXX
               BATCH (M+ 3) =  GXZXX
               BATCH (M+ 4) =  GXYXX
               BATCH (M+ 5) =  GYYXX
               BATCH (M+ 6) =  GYZXX
               BATCH (M+ 7) =  GXZXX
               BATCH (M+ 8) =  GYZXX
               BATCH (M+ 9) =  GZZXX
               BATCH (M+10) =  GXXXY
               BATCH (M+11) =  GXYXY
               BATCH (M+12) =  GXZXY
               BATCH (M+13) =  GXYXY
               BATCH (M+14) =  GYYXY
               BATCH (M+15) =  GYZXY
               BATCH (M+16) =  GXZXY
               BATCH (M+17) =  GYZXY
               BATCH (M+18) =  GZZXY
               BATCH (M+19) =  GXXXZ
               BATCH (M+20) =  GXYXZ
               BATCH (M+21) =  GXZXZ
               BATCH (M+22) =  GXYXZ
               BATCH (M+23) =  GYYXZ
               BATCH (M+24) =  GYZXZ
               BATCH (M+25) =  GXZXZ
               BATCH (M+26) =  GYZXZ
               BATCH (M+27) =  GZZXZ
               BATCH (M+28) =  GXXXY
               BATCH (M+29) =  GXYXY
               BATCH (M+30) =  GXZXY
               BATCH (M+31) =  GXYXY
               BATCH (M+32) =  GYYXY
               BATCH (M+33) =  GYZXY
               BATCH (M+34) =  GXZXY
               BATCH (M+35) =  GYZXY
               BATCH (M+36) =  GZZXY
               BATCH (M+37) =  GXXYY
               BATCH (M+38) =  GXYYY
               BATCH (M+39) =  GXZYY
               BATCH (M+40) =  GXYYY
               BATCH (M+41) =  GYYYY
               BATCH (M+42) =  GYZYY
               BATCH (M+43) =  GXZYY
               BATCH (M+44) =  GYZYY
               BATCH (M+45) =  GZZYY
               BATCH (M+46) =  GXXYZ
               BATCH (M+47) =  GXYYZ
               BATCH (M+48) =  GXZYZ
               BATCH (M+49) =  GXYYZ
               BATCH (M+50) =  GYYYZ
               BATCH (M+51) =  GYZYZ
               BATCH (M+52) =  GXZYZ
               BATCH (M+53) =  GYZYZ
               BATCH (M+54) =  GZZYZ
               BATCH (M+55) =  GXXXZ
               BATCH (M+56) =  GXYXZ
               BATCH (M+57) =  GXZXZ
               BATCH (M+58) =  GXYXZ
               BATCH (M+59) =  GYYXZ
               BATCH (M+60) =  GYZXZ
               BATCH (M+61) =  GXZXZ
               BATCH (M+62) =  GYZXZ
               BATCH (M+63) =  GZZXZ
               BATCH (M+64) =  GXXYZ
               BATCH (M+65) =  GXYYZ
               BATCH (M+66) =  GXZYZ
               BATCH (M+67) =  GXYYZ
               BATCH (M+68) =  GYYYZ
               BATCH (M+69) =  GYZYZ
               BATCH (M+70) =  GXZYZ
               BATCH (M+71) =  GYZYZ
               BATCH (M+72) =  GZZYZ
               BATCH (M+73) =  GXXZZ
               BATCH (M+74) =  GXYZZ
               BATCH (M+75) =  GXZZZ
               BATCH (M+76) =  GXYZZ
               BATCH (M+77) =  GYYZZ
               BATCH (M+78) =  GYZZZ
               BATCH (M+79) =  GXZZZ
               BATCH (M+80) =  GYZZZ
               BATCH (M+81) =  GZZZZ

               M = M + 81

 2200       CONTINUE
 2000    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 3-center (AB|CC)
C                type integrals.
C
C
   30    M = 0
         DO 3000 IJ = 1,MIJ
            PVAL = P (IJ)
            PXVAL = PX (IJ)
            PYVAL = PY (IJ)
            PZVAL = PZ (IJ)
            PQX = PXVAL - X3
            PQY = PYVAL - Y3
            PQZ = PZVAL - Z3
            RNPQSQ = PQX*PQX + PQY*PQY + PQZ*PQZ
            PSCALE = SCALEP (IJ)
            U3 = HALF / PVAL
            XSPSS1 = PXVAL - X2
            XPSSS1 = PXVAL - X1
            YSPSS1 = PYVAL - Y2
            YPSSS1 = PYVAL - Y1
            ZSPSS1 = PZVAL - Z2
            ZPSSS1 = PZVAL - Z1
            DO 3300 KL = 1,MKL
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

                   F4 = ((((( FTABLE (10,TGRID) * DELTA6
     +                      + FTABLE (9,TGRID)) * DELTA5
     +                      + FTABLE (8,TGRID)) * DELTA4
     +                      + FTABLE (7,TGRID)) * DELTA3
     +                      + FTABLE (6,TGRID)) * DELTA2
     +                      + FTABLE (5,TGRID)) * DELTA1
     +                      + FTABLE (4,TGRID)

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2
                   F3 = SCALE * F3
                   F4 = SCALE * F4

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1
                   F3 = FIVE  * T2INV * F2
                   F4 = SEVEN * T2INV * F3

               END IF

               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               U2 = HALF * PQPINV
               U4 = HALF / QVAL
               U5 = U2 + PQPINV
               U6 = U5 + PQPINV
               U7 = - U0 * U4
               U8 = U1 * U3
C
C
C             ...the X-terms (with exception of XSPSS1 and XPSSS1,
C                which can be evaluated in the P-loop).
C
C
               XSSSP2 = PQX * U0
               XSPSS2 = PQX * U1
               B      = XSPSS1 + XPSSS1
               XSSPP1 = U4
               XPPSS1 = XSPSS1 * XPSSS1 + U3
               XSPSP2 = XSPSS1 * XSSSP2 + U2
               XPSSP2 = XPSSS1 * XSSSP2 + U2
               XSSPP2 = U7
               XPPSS2 = B * XSPSS2 + U8
               XSPSP3 = XSSSP2 * XSPSS2
               XSSPP3 = XSSSP2 * XSSSP2
               XPPSS3 = XSPSS2 * XSPSS2
               XSPPP1 = XSSPP1 * XSPSS1
               XPSPP1 = XSSPP1 * XPSSS1
               D      = XSPSS2 * XSSPP1
               XPPSP2 = XSSSP2 * XPPSS1 + B * U2
               XSPPP2 = XSPSS1 * XSSPP2 + D
               XPSPP2 = XPSSS1 * XSSPP2 + D
               D      = XSSSP2 * U5
               XPPSP3 = B * XSPSP3 + XSPSS2 * U5
               XSPPP3 = XSPSS1 * XSSPP3 + D
               XPSPP3 = XPSSS1 * XSSPP3 + D
               XSPPP4 = XSSPP3 * XSPSS2
               XPPSP4 = XPPSS3 * XSSSP2
               D      = B * XSSSP2 + U2
               XPPPP1 = XSSPP1 * XPPSS1
               XPPPP2 = XSSPP2 * XPPSS1 + XPPSS2 * XSSPP1
               XPPPP3 = XPPSS1 * XSSPP3 + XSSPP1 * XPPSS3 + U5 * D
               XPPPP4 = XSPSP3 * (D + U6)
               XPPPP5 = XSSPP3 * XPPSS3
C
C
C             ...the Y-terms (with exception of YSPSS1 and YPSSS1,
C                which can be evaluated in the P-loop).
C
C
               YSSSP2 = PQY * U0
               YSPSS2 = PQY * U1
               B      = YSPSS1 + YPSSS1
               YSSPP1 = U4
               YPPSS1 = YSPSS1 * YPSSS1 + U3
               YSPSP2 = YSPSS1 * YSSSP2 + U2
               YPSSP2 = YPSSS1 * YSSSP2 + U2
               YSSPP2 = U7
               YPPSS2 = B * YSPSS2 + U8
               YSPSP3 = YSSSP2 * YSPSS2
               YSSPP3 = YSSSP2 * YSSSP2
               YPPSS3 = YSPSS2 * YSPSS2
               YSPPP1 = YSSPP1 * YSPSS1
               YPSPP1 = YSSPP1 * YPSSS1
               D      = YSPSS2 * YSSPP1
               YPPSP2 = YSSSP2 * YPPSS1 + B * U2
               YSPPP2 = YSPSS1 * YSSPP2 + D
               YPSPP2 = YPSSS1 * YSSPP2 + D
               D      = YSSSP2 * U5
               YPPSP3 = B * YSPSP3 + YSPSS2 * U5
               YSPPP3 = YSPSS1 * YSSPP3 + D
               YPSPP3 = YPSSS1 * YSSPP3 + D
               YSPPP4 = YSSPP3 * YSPSS2
               YPPSP4 = YPPSS3 * YSSSP2
               D      = B * YSSSP2 + U2
               YPPPP1 = YSSPP1 * YPPSS1
               YPPPP2 = YSSPP2 * YPPSS1 + YPPSS2 * YSSPP1
               YPPPP3 = YPPSS1 * YSSPP3 + YSSPP1 * YPPSS3 + U5 * D
               YPPPP4 = YSPSP3 * (D + U6)
               YPPPP5 = YSSPP3 * YPPSS3
C
C
C             ...the Z-terms (with exception of ZSPSS1 and ZPSSS1,
C                which can be evaluated in the P-loop).
C
C
               ZSSSP2 = PQZ * U0
               ZSPSS2 = PQZ * U1
               B      = ZSPSS1 + ZPSSS1
               ZSSPP1 = U4
               ZPPSS1 = ZSPSS1 * ZPSSS1 + U3
               ZSPSP2 = ZSPSS1 * ZSSSP2 + U2
               ZPSSP2 = ZPSSS1 * ZSSSP2 + U2
               ZSSPP2 = U7
               ZPPSS2 = B * ZSPSS2 + U8
               ZSPSP3 = ZSSSP2 * ZSPSS2
               ZSSPP3 = ZSSSP2 * ZSSSP2
               ZPPSS3 = ZSPSS2 * ZSPSS2
               ZSPPP1 = ZSSPP1 * ZSPSS1
               ZPSPP1 = ZSSPP1 * ZPSSS1
               D      = ZSPSS2 * ZSSPP1
               ZPPSP2 = ZSSSP2 * ZPPSS1 + B * U2
               ZSPPP2 = ZSPSS1 * ZSSPP2 + D
               ZPSPP2 = ZPSSS1 * ZSSPP2 + D
               D      = ZSSSP2 * U5
               ZPPSP3 = B * ZSPSP3 + ZSPSS2 * U5
               ZSPPP3 = ZSPSS1 * ZSSPP3 + D
               ZPSPP3 = ZPSSS1 * ZSSPP3 + D
               ZSPPP4 = ZSSPP3 * ZSPSS2
               ZPPSP4 = ZPPSS3 * ZSSSP2
               D      = B * ZSSSP2 + U2
               ZPPPP1 = ZSSPP1 * ZPPSS1
               ZPPPP2 = ZSSPP2 * ZPPSS1 + ZPPSS2 * ZSSPP1
               ZPPPP3 = ZPPSS1 * ZSSPP3 + ZSSPP1 * ZPPSS3 + U5 * D
               ZPPPP4 = ZSPSP3 * (D + U6)
               ZPPPP5 = ZSSPP3 * ZPPSS3
C
C
C             ...assemble the 3-center (AB|CC) type integrals.
C
C
               GXXXX = XPPPP1*F0+XPPPP2*F1+XPPPP3*F2+XPPPP4*F3+XPPPP5*F4
               GYYYY = YPPPP1*F0+YPPPP2*F1+YPPPP3*F2+YPPPP4*F3+YPPPP5*F4
               GZZZZ = ZPPPP1*F0+ZPPPP2*F1+ZPPPP3*F2+ZPPPP4*F3+ZPPPP5*F4

               A = XPPSP2 * F2 + XPPSP3 * F3 + XPPSP4 * F4
               B = XPSPP1 * F0 + XPSPP2 * F1 + XPSPP3 * F2 + XSPPP4 * F3
               C = XPSPP1 * F1 + XPSPP2 * F2 + XPSPP3 * F3 + XSPPP4 * F4
               D = XSPPP1 * F0 + XSPPP2 * F1 + XSPPP3 * F2 + XSPPP4 * F3
               E = XSPPP1 * F1 + XSPPP2 * F2 + XSPPP3 * F3 + XSPPP4 * F4

               GXXXY = A * YSSSP2
               GXYXX = B * YSPSS1 + C * YSPSS2
               GYXXX = D * YPSSS1 + E * YSPSS2
               GXXXZ = A * ZSSSP2
               GXZXX = B * ZSPSS1 + C * ZSPSS2
               GZXXX = D * ZPSSS1 + E * ZSPSS2

               A = YPPSP2 * F2 + YPPSP3 * F3 + YPPSP4 * F4
               B = YPSPP1 * F0 + YPSPP2 * F1 + YPSPP3 * F2 + YSPPP4 * F3
               C = YPSPP1 * F1 + YPSPP2 * F2 + YPSPP3 * F3 + YSPPP4 * F4
               D = YSPPP1 * F0 + YSPPP2 * F1 + YSPPP3 * F2 + YSPPP4 * F3
               E = YSPPP1 * F1 + YSPPP2 * F2 + YSPPP3 * F3 + YSPPP4 * F4

               GYYXY = A * XSSSP2
               GYXYY = B * XSPSS1 + C * XSPSS2
               GXYYY = D * XPSSS1 + E * XSPSS2
               GYYYZ = A * ZSSSP2
               GYZYY = B * ZSPSS1 + C * ZSPSS2
               GZYYY = D * ZPSSS1 + E * ZSPSS2

               A = ZPPSP2 * F2 + ZPPSP3 * F3 + ZPPSP4 * F4
               B = ZPSPP1 * F0 + ZPSPP2 * F1 + ZPSPP3 * F2 + ZSPPP4 * F3
               C = ZPSPP1 * F1 + ZPSPP2 * F2 + ZPSPP3 * F3 + ZSPPP4 * F4
               D = ZSPPP1 * F0 + ZSPPP2 * F1 + ZSPPP3 * F2 + ZSPPP4 * F3
               E = ZSPPP1 * F1 + ZSPPP2 * F2 + ZSPPP3 * F3 + ZSPPP4 * F4

               GZZXZ = A * XSSSP2
               GZXZZ = B * XSPSS1 + C * XSPSS2
               GXZZZ = D * XPSSS1 + E * XSPSS2
               GZZYZ = A * YSSSP2
               GZYZZ = B * YSPSS1 + C * YSPSS2
               GYZZZ = D * YPSSS1 + E * YSPSS2

               A = XSSPP1 * F0 + XSSPP2 * F1 + XSSPP3 * F2
               B = XSSPP1 * F1 + XSSPP2 * F2 + XSSPP3 * F3
               C = XSSPP1 * F2 + XSSPP2 * F3 + XSSPP3 * F4
               D = YSSPP1 * F0 + YSSPP2 * F1 + YSSPP3 * F2
               E = YSSPP1 * F1 + YSSPP2 * F2 + YSSPP3 * F3
               F = YSSPP1 * F2 + YSSPP2 * F3 + YSSPP3 * F4
               G = ZSSPP1 * F0 + ZSSPP2 * F1 + ZSSPP3 * F2
               H = ZSSPP1 * F1 + ZSSPP2 * F2 + ZSSPP3 * F3
               R = ZSSPP1 * F2 + ZSSPP2 * F3 + ZSSPP3 * F4

               GXXYY = D * XPPSS1 + E * XPPSS2 + F * XPPSS3
               GXXZZ = G * XPPSS1 + H * XPPSS2 + R * XPPSS3
               GYYXX = A * YPPSS1 + B * YPPSS2 + C * YPPSS3
               GYYZZ = G * YPPSS1 + H * YPPSS2 + R * YPPSS3
               GZZXX = A * ZPPSS1 + B * ZPPSS2 + C * ZPPSS3
               GZZYY = D * ZPPSS1 + E * ZPPSS2 + F * ZPPSS3

               GYZXX = (A * YPSSS1 + B * YSPSS2) * ZSPSS1 +
     +                 (B * YPSSS1 + C * YSPSS2) * ZSPSS2
               GZYXX = (A * ZPSSS1 + B * ZSPSS2) * YSPSS1 +
     +                 (B * ZPSSS1 + C * ZSPSS2) * YSPSS2
               GXZYY = (D * XPSSS1 + E * XSPSS2) * ZSPSS1 +
     +                 (E * XPSSS1 + F * XSPSS2) * ZSPSS2
               GZXYY = (D * ZPSSS1 + E * ZSPSS2) * XSPSS1 +
     +                 (E * ZPSSS1 + F * ZSPSS2) * XSPSS2
               GXYZZ = (G * XPSSS1 + H * XSPSS2) * YSPSS1 +
     +                 (H * XPSSS1 + R * XSPSS2) * YSPSS2
               GYXZZ = (G * YPSSS1 + H * YSPSS2) * XSPSS1 +
     +                 (H * YPSSS1 + R * YSPSS2) * XSPSS2

               GXXYZ = (XPPSS1*F2+XPPSS2*F3+XPPSS3*F4) * YSSSP2 * ZSSSP2
               GYYXZ = (YPPSS1*F2+YPPSS2*F3+YPPSS3*F4) * XSSSP2 * ZSSSP2
               GZZXY = (ZPPSS1*F2+ZPPSS2*F3+ZPPSS3*F4) * XSSSP2 * YSSSP2

               AA = XSPSP3 * F3
               BB = XSPSP3 * F4
               CC = YSPSP3 * F3
               DD = YSPSP3 * F4
               EE = ZSPSP3 * F3
               FF = ZSPSP3 * F4

               A = XPSSP2 * F2 + AA
               B = XPSSP2 * F3 + BB
               C = YPSSP2 * F2 + CC
               D = YPSSP2 * F3 + DD
               E = ZPSSP2 * F2 + EE
               F = ZPSSP2 * F3 + FF

               GXYXY = A * YSPSP2 + B * YSPSP3
               GXZXZ = A * ZSPSP2 + B * ZSPSP3
               GYXXY = C * XSPSP2 + D * XSPSP3
               GYZYZ = C * ZSPSP2 + D * ZSPSP3
               GZXXZ = E * XSPSP2 + F * XSPSP3
               GZYYZ = E * YSPSP2 + F * YSPSP3

               GXYXZ = (A * YSPSS1 + B * YSPSS2) * ZSSSP2
               GXZXY = (A * ZSPSS1 + B * ZSPSS2) * YSSSP2
               GYXYZ = (C * XSPSS1 + D * XSPSS2) * ZSSSP2
               GYZXY = (C * ZSPSS1 + D * ZSPSS2) * XSSSP2
               GZXYZ = (E * XSPSS1 + F * XSPSS2) * YSSSP2
               GZYXZ = (E * YSPSS1 + F * YSPSS2) * XSSSP2

               A = XSPSP2 * F2 + AA
               B = XSPSP2 * F3 + BB
               C = YSPSP2 * F2 + CC
               D = YSPSP2 * F3 + DD
               E = ZSPSP2 * F2 + EE
               F = ZSPSP2 * F3 + FF

               GYXXZ = (A * YPSSS1 + B * YSPSS2) * ZSSSP2
               GZXXY = (A * ZPSSS1 + B * ZSPSS2) * YSSSP2
               GXYYZ = (C * XPSSS1 + D * XSPSS2) * ZSSSP2
               GZYXY = (C * ZPSSS1 + D * ZSPSS2) * XSSSP2
               GXZYZ = (E * XPSSS1 + F * XSPSS2) * YSSSP2
               GYZXZ = (E * YPSSS1 + F * YSPSS2) * XSSSP2

               BATCH (M+ 1) =  GXXXX
               BATCH (M+ 2) =  GYXXX
               BATCH (M+ 3) =  GZXXX
               BATCH (M+ 4) =  GXYXX
               BATCH (M+ 5) =  GYYXX
               BATCH (M+ 6) =  GZYXX
               BATCH (M+ 7) =  GXZXX
               BATCH (M+ 8) =  GYZXX
               BATCH (M+ 9) =  GZZXX
               BATCH (M+10) =  GXXXY
               BATCH (M+11) =  GYXXY
               BATCH (M+12) =  GZXXY
               BATCH (M+13) =  GXYXY
               BATCH (M+14) =  GYYXY
               BATCH (M+15) =  GZYXY
               BATCH (M+16) =  GXZXY
               BATCH (M+17) =  GYZXY
               BATCH (M+18) =  GZZXY
               BATCH (M+19) =  GXXXZ
               BATCH (M+20) =  GYXXZ
               BATCH (M+21) =  GZXXZ
               BATCH (M+22) =  GXYXZ
               BATCH (M+23) =  GYYXZ
               BATCH (M+24) =  GZYXZ
               BATCH (M+25) =  GXZXZ
               BATCH (M+26) =  GYZXZ
               BATCH (M+27) =  GZZXZ
               BATCH (M+28) =  GXXXY
               BATCH (M+29) =  GYXXY
               BATCH (M+30) =  GZXXY
               BATCH (M+31) =  GXYXY
               BATCH (M+32) =  GYYXY
               BATCH (M+33) =  GZYXY
               BATCH (M+34) =  GXZXY
               BATCH (M+35) =  GYZXY
               BATCH (M+36) =  GZZXY
               BATCH (M+37) =  GXXYY
               BATCH (M+38) =  GYXYY
               BATCH (M+39) =  GZXYY
               BATCH (M+40) =  GXYYY
               BATCH (M+41) =  GYYYY
               BATCH (M+42) =  GZYYY
               BATCH (M+43) =  GXZYY
               BATCH (M+44) =  GYZYY
               BATCH (M+45) =  GZZYY
               BATCH (M+46) =  GXXYZ
               BATCH (M+47) =  GYXYZ
               BATCH (M+48) =  GZXYZ
               BATCH (M+49) =  GXYYZ
               BATCH (M+50) =  GYYYZ
               BATCH (M+51) =  GZYYZ
               BATCH (M+52) =  GXZYZ
               BATCH (M+53) =  GYZYZ
               BATCH (M+54) =  GZZYZ
               BATCH (M+55) =  GXXXZ
               BATCH (M+56) =  GYXXZ
               BATCH (M+57) =  GZXXZ
               BATCH (M+58) =  GXYXZ
               BATCH (M+59) =  GYYXZ
               BATCH (M+60) =  GZYXZ
               BATCH (M+61) =  GXZXZ
               BATCH (M+62) =  GYZXZ
               BATCH (M+63) =  GZZXZ
               BATCH (M+64) =  GXXYZ
               BATCH (M+65) =  GYXYZ
               BATCH (M+66) =  GZXYZ
               BATCH (M+67) =  GXYYZ
               BATCH (M+68) =  GYYYZ
               BATCH (M+69) =  GZYYZ
               BATCH (M+70) =  GXZYZ
               BATCH (M+71) =  GYZYZ
               BATCH (M+72) =  GZZYZ
               BATCH (M+73) =  GXXZZ
               BATCH (M+74) =  GYXZZ
               BATCH (M+75) =  GZXZZ
               BATCH (M+76) =  GXYZZ
               BATCH (M+77) =  GYYZZ
               BATCH (M+78) =  GZYZZ
               BATCH (M+79) =  GXZZZ
               BATCH (M+80) =  GYZZZ
               BATCH (M+81) =  GZZZZ

               M = M + 81

 3300       CONTINUE
 3000    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 3-center (AA|CD)
C                type integrals.
C
C
   40    M = 0
         DO 4000 IJ = 1,MIJ
            PVAL = P (IJ)
            PSCALE = SCALEP (IJ)
            U3 = HALF / PVAL
            DO 4400 KL = 1,MKL
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

                   F4 = ((((( FTABLE (10,TGRID) * DELTA6
     +                      + FTABLE (9,TGRID)) * DELTA5
     +                      + FTABLE (8,TGRID)) * DELTA4
     +                      + FTABLE (7,TGRID)) * DELTA3
     +                      + FTABLE (6,TGRID)) * DELTA2
     +                      + FTABLE (5,TGRID)) * DELTA1
     +                      + FTABLE (4,TGRID)

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2
                   F3 = SCALE * F3
                   F4 = SCALE * F4

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1
                   F3 = FIVE  * T2INV * F2
                   F4 = SEVEN * T2INV * F3

               END IF

               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               U2 = HALF * PQPINV
               U4 = HALF / QVAL
               U5 = U2 + PQPINV
               U6 = U5 + PQPINV
               U7 = - U0 * U4
               U8 = U1 * U3
C
C
C             ...the X-terms.
C
C
               XSSSP1 = QXVAL - X4
               XSSPS1 = QXVAL - X3
               XSSSP2 = PQX * U0
               XSPSS2 = PQX * U1
               A      = XSSSP1 + XSSPS1
               XSSPP1 = XSSSP1 * XSSPS1 + U4
               XPPSS1 = U3
               XSPSP2 = XSSSP1 * XSPSS2 + U2
               XSPPS2 = XSSPS1 * XSPSS2 + U2
               XSSPP2 = A * XSSSP2 + U7
               XPPSS2 = U8
               XSPSP3 = XSSSP2 * XSPSS2
               XSSPP3 = XSSSP2 * XSSSP2
               XPPSS3 = XSPSS2 * XSPSS2
               XPPSP1 = XPPSS1 * XSSSP1
               XPPPS1 = XPPSS1 * XSSPS1
               XSPPP2 = XSPSS2 * XSSPP1 + A * U2
               E      = XSSSP2 * XPPSS1
               XPPSP2 = XSSSP1 * XPPSS2 + E
               XPPPS2 = XSSPS1 * XPPSS2 + E
               XSPPP3 = A * XSPSP3 + XSSSP2 * U5
               E      = XSPSS2 * U5
               XPPSP3 = XSSSP1 * XPPSS3 + E
               XPPPS3 = XSSPS1 * XPPSS3 + E
               XSPPP4 = XSSPP3 * XSPSS2
               XPPSP4 = XPPSS3 * XSSSP2
               D      = A * XSPSS2 + U2
               XPPPP1 = XSSPP1 * XPPSS1
               XPPPP2 = XSSPP2 * XPPSS1 + XPPSS2 * XSSPP1
               XPPPP3 = XPPSS1 * XSSPP3 + XSSPP1 * XPPSS3 + U5 * D
               XPPPP4 = XSPSP3 * (D + U6)
               XPPPP5 = XSSPP3 * XPPSS3
C
C
C             ...the Y-terms.
C
C
               YSSSP1 = QYVAL - Y4
               YSSPS1 = QYVAL - Y3
               YSSSP2 = PQY * U0
               YSPSS2 = PQY * U1
               A      = YSSSP1 + YSSPS1
               YSSPP1 = YSSSP1 * YSSPS1 + U4
               YPPSS1 = U3
               YSPSP2 = YSSSP1 * YSPSS2 + U2
               YSPPS2 = YSSPS1 * YSPSS2 + U2
               YSSPP2 = A * YSSSP2 + U7
               YPPSS2 = U8
               YSPSP3 = YSSSP2 * YSPSS2
               YSSPP3 = YSSSP2 * YSSSP2
               YPPSS3 = YSPSS2 * YSPSS2
               YPPSP1 = YPPSS1 * YSSSP1
               YPPPS1 = YPPSS1 * YSSPS1
               YSPPP2 = YSPSS2 * YSSPP1 + A * U2
               E      = YSSSP2 * YPPSS1
               YPPSP2 = YSSSP1 * YPPSS2 + E
               YPPPS2 = YSSPS1 * YPPSS2 + E
               YSPPP3 = A * YSPSP3 + YSSSP2 * U5
               E      = YSPSS2 * U5
               YPPSP3 = YSSSP1 * YPPSS3 + E
               YPPPS3 = YSSPS1 * YPPSS3 + E
               YSPPP4 = YSSPP3 * YSPSS2
               YPPSP4 = YPPSS3 * YSSSP2
               D      = A * YSPSS2 + U2
               YPPPP1 = YSSPP1 * YPPSS1
               YPPPP2 = YSSPP2 * YPPSS1 + YPPSS2 * YSSPP1
               YPPPP3 = YPPSS1 * YSSPP3 + YSSPP1 * YPPSS3 + U5 * D
               YPPPP4 = YSPSP3 * (D + U6)
               YPPPP5 = YSSPP3 * YPPSS3
C
C
C             ...the Z-terms.
C
C
               ZSSSP1 = QZVAL - Z4
               ZSSPS1 = QZVAL - Z3
               ZSSSP2 = PQZ * U0
               ZSPSS2 = PQZ * U1
               A      = ZSSSP1 + ZSSPS1
               ZSSPP1 = ZSSSP1 * ZSSPS1 + U4
               ZPPSS1 = U3
               ZSPSP2 = ZSSSP1 * ZSPSS2 + U2
               ZSPPS2 = ZSSPS1 * ZSPSS2 + U2
               ZSSPP2 = A * ZSSSP2 + U7
               ZPPSS2 = U8
               ZSPSP3 = ZSSSP2 * ZSPSS2
               ZSSPP3 = ZSSSP2 * ZSSSP2
               ZPPSS3 = ZSPSS2 * ZSPSS2
               ZPPSP1 = ZPPSS1 * ZSSSP1
               ZPPPS1 = ZPPSS1 * ZSSPS1
               ZSPPP2 = ZSPSS2 * ZSSPP1 + A * U2
               E      = ZSSSP2 * ZPPSS1
               ZPPSP2 = ZSSSP1 * ZPPSS2 + E
               ZPPPS2 = ZSSPS1 * ZPPSS2 + E
               ZSPPP3 = A * ZSPSP3 + ZSSSP2 * U5
               E      = ZSPSS2 * U5
               ZPPSP3 = ZSSSP1 * ZPPSS3 + E
               ZPPPS3 = ZSSPS1 * ZPPSS3 + E
               ZSPPP4 = ZSSPP3 * ZSPSS2
               ZPPSP4 = ZPPSS3 * ZSSSP2
               D      = A * ZSPSS2 + U2
               ZPPPP1 = ZSSPP1 * ZPPSS1
               ZPPPP2 = ZSSPP2 * ZPPSS1 + ZPPSS2 * ZSSPP1
               ZPPPP3 = ZPPSS1 * ZSSPP3 + ZSSPP1 * ZPPSS3 + U5 * D
               ZPPPP4 = ZSPSP3 * (D + U6)
               ZPPPP5 = ZSSPP3 * ZPPSS3
C
C
C             ...assemble the 3-center (AA|CD) type integrals.
C
C
               GXXXX = XPPPP1*F0+XPPPP2*F1+XPPPP3*F2+XPPPP4*F3+XPPPP5*F4
               GYYYY = YPPPP1*F0+YPPPP2*F1+YPPPP3*F2+YPPPP4*F3+YPPPP5*F4
               GZZZZ = ZPPPP1*F0+ZPPPP2*F1+ZPPPP3*F2+ZPPPP4*F3+ZPPPP5*F4

               A = XPPPS1 * F0 + XPPPS2 * F1 + XPPPS3 * F2 + XPPSP4 * F3
               B = XPPPS1 * F1 + XPPPS2 * F2 + XPPPS3 * F3 + XPPSP4 * F4
               C = XPPSP1 * F0 + XPPSP2 * F1 + XPPSP3 * F2 + XPPSP4 * F3
               D = XPPSP1 * F1 + XPPSP2 * F2 + XPPSP3 * F3 + XPPSP4 * F4
               E = XSPPP2 * F2 + XSPPP3 * F3 + XSPPP4 * F4

               GXXXY = A * YSSSP1 + B * YSSSP2
               GXXYX = C * YSSPS1 + D * YSSSP2
               GXYXX = E * YSPSS2
               GXXXZ = A * ZSSSP1 + B * ZSSSP2
               GXXZX = C * ZSSPS1 + D * ZSSSP2
               GXZXX = E * ZSPSS2

               A = YPPPS1 * F0 + YPPPS2 * F1 + YPPPS3 * F2 + YPPSP4 * F3
               B = YPPPS1 * F1 + YPPPS2 * F2 + YPPPS3 * F3 + YPPSP4 * F4
               C = YPPSP1 * F0 + YPPSP2 * F1 + YPPSP3 * F2 + YPPSP4 * F3
               D = YPPSP1 * F1 + YPPSP2 * F2 + YPPSP3 * F3 + YPPSP4 * F4
               E = YSPPP2 * F2 + YSPPP3 * F3 + YSPPP4 * F4

               GYYYX = A * XSSSP1 + B * XSSSP2
               GYYXY = C * XSSPS1 + D * XSSSP2
               GXYYY = E * XSPSS2
               GYYYZ = A * ZSSSP1 + B * ZSSSP2
               GYYZY = C * ZSSPS1 + D * ZSSSP2
               GYZYY = E * ZSPSS2

               A = ZPPPS1 * F0 + ZPPPS2 * F1 + ZPPPS3 * F2 + ZPPSP4 * F3
               B = ZPPPS1 * F1 + ZPPPS2 * F2 + ZPPPS3 * F3 + ZPPSP4 * F4
               C = ZPPSP1 * F0 + ZPPSP2 * F1 + ZPPSP3 * F2 + ZPPSP4 * F3
               D = ZPPSP1 * F1 + ZPPSP2 * F2 + ZPPSP3 * F3 + ZPPSP4 * F4
               E = ZSPPP2 * F2 + ZSPPP3 * F3 + ZSPPP4 * F4

               GZZZX = A * XSSSP1 + B * XSSSP2
               GZZXZ = C * XSSPS1 + D * XSSSP2
               GXZZZ = E * XSPSS2
               GZZZY = A * YSSSP1 + B * YSSSP2
               GZZYZ = C * YSSPS1 + D * YSSSP2
               GYZZZ = E * YSPSS2

               A = XPPSS1 * F0 + XPPSS2 * F1 + XPPSS3 * F2
               B = XPPSS1 * F1 + XPPSS2 * F2 + XPPSS3 * F3
               C = XPPSS1 * F2 + XPPSS2 * F3 + XPPSS3 * F4
               D = YPPSS1 * F0 + YPPSS2 * F1 + YPPSS3 * F2
               E = YPPSS1 * F1 + YPPSS2 * F2 + YPPSS3 * F3
               F = YPPSS1 * F2 + YPPSS2 * F3 + YPPSS3 * F4
               G = ZPPSS1 * F0 + ZPPSS2 * F1 + ZPPSS3 * F2
               H = ZPPSS1 * F1 + ZPPSS2 * F2 + ZPPSS3 * F3
               R = ZPPSS1 * F2 + ZPPSS2 * F3 + ZPPSS3 * F4

               GXXYY = A * YSSPP1 + B * YSSPP2 + C * YSSPP3
               GXXZZ = A * ZSSPP1 + B * ZSSPP2 + C * ZSSPP3
               GYYXX = D * XSSPP1 + E * XSSPP2 + F * XSSPP3
               GYYZZ = D * ZSSPP1 + E * ZSSPP2 + F * ZSSPP3
               GZZXX = G * XSSPP1 + H * XSSPP2 + R * XSSPP3
               GZZYY = G * YSSPP1 + H * YSSPP2 + R * YSSPP3

               GXXYZ = (A * YSSPS1 + B * YSSSP2) * ZSSSP1 +
     +                 (B * YSSPS1 + C * YSSSP2) * ZSSSP2
               GXXZY = (A * ZSSPS1 + B * ZSSSP2) * YSSSP1 +
     +                 (B * ZSSPS1 + C * ZSSSP2) * YSSSP2
               GYYXZ = (D * XSSPS1 + E * XSSSP2) * ZSSSP1 +
     +                 (E * XSSPS1 + F * XSSSP2) * ZSSSP2
               GYYZX = (D * ZSSPS1 + E * ZSSSP2) * XSSSP1 +
     +                 (E * ZSSPS1 + F * ZSSSP2) * XSSSP2
               GZZXY = (G * XSSPS1 + H * XSSSP2) * YSSSP1 +
     +                 (H * XSSPS1 + R * XSSSP2) * YSSSP2
               GZZYX = (G * YSSPS1 + H * YSSSP2) * XSSSP1 +
     +                 (H * YSSPS1 + R * YSSSP2) * XSSSP2

               GYZXX = (XSSPP1*F2+XSSPP2*F3+XSSPP3*F4) * YSPSS2 * ZSPSS2
               GXZYY = (YSSPP1*F2+YSSPP2*F3+YSSPP3*F4) * XSPSS2 * ZSPSS2
               GXYZZ = (ZSSPP1*F2+ZSSPP2*F3+ZSSPP3*F4) * XSPSS2 * YSPSS2

               AA = XSPSP3 * F3
               BB = XSPSP3 * F4
               CC = YSPSP3 * F3
               DD = YSPSP3 * F4
               EE = ZSPSP3 * F3
               FF = ZSPSP3 * F4

               A = XSPPS2 * F2 + AA
               B = XSPPS2 * F3 + BB
               C = YSPPS2 * F2 + CC
               D = YSPPS2 * F3 + DD
               E = ZSPPS2 * F2 + EE
               F = ZSPPS2 * F3 + FF

               GXYXY = A * YSPSP2 + B * YSPSP3
               GXZXZ = A * ZSPSP2 + B * ZSPSP3
               GXYYX = C * XSPSP2 + D * XSPSP3
               GYZYZ = C * ZSPSP2 + D * ZSPSP3
               GXZZX = E * XSPSP2 + F * XSPSP3
               GYZZY = E * YSPSP2 + F * YSPSP3

               GXYXZ = YSPSS2 * (A * ZSSSP1 + B * ZSSSP2)
               GXZXY = ZSPSS2 * (A * YSSSP1 + B * YSSSP2)
               GXYYZ = XSPSS2 * (C * ZSSSP1 + D * ZSSSP2)
               GYZYX = ZSPSS2 * (C * XSSSP1 + D * XSSSP2)
               GXZZY = XSPSS2 * (E * YSSSP1 + F * YSSSP2)
               GYZZX = YSPSS2 * (E * XSSSP1 + F * XSSSP2)

               A = XSPSP2 * F2 + AA
               B = XSPSP2 * F3 + BB
               C = YSPSP2 * F2 + CC
               D = YSPSP2 * F3 + DD
               E = ZSPSP2 * F2 + EE
               F = ZSPSP2 * F3 + FF

               GXYZX = YSPSS2 * (A * ZSSPS1 + B * ZSSSP2)
               GXZYX = ZSPSS2 * (A * YSSPS1 + B * YSSSP2)
               GXYZY = XSPSS2 * (C * ZSSPS1 + D * ZSSSP2)
               GYZXY = ZSPSS2 * (C * XSSPS1 + D * XSSSP2)
               GXZYZ = XSPSS2 * (E * YSSPS1 + F * YSSSP2)
               GYZXZ = YSPSS2 * (E * XSSPS1 + F * XSSSP2)

               BATCH (M+ 1) =  GXXXX
               BATCH (M+ 2) =  GXYXX
               BATCH (M+ 3) =  GXZXX
               BATCH (M+ 4) =  GXYXX
               BATCH (M+ 5) =  GYYXX
               BATCH (M+ 6) =  GYZXX
               BATCH (M+ 7) =  GXZXX
               BATCH (M+ 8) =  GYZXX
               BATCH (M+ 9) =  GZZXX
               BATCH (M+10) =  GXXYX
               BATCH (M+11) =  GXYYX
               BATCH (M+12) =  GXZYX
               BATCH (M+13) =  GXYYX
               BATCH (M+14) =  GYYYX
               BATCH (M+15) =  GYZYX
               BATCH (M+16) =  GXZYX
               BATCH (M+17) =  GYZYX
               BATCH (M+18) =  GZZYX
               BATCH (M+19) =  GXXZX
               BATCH (M+20) =  GXYZX
               BATCH (M+21) =  GXZZX
               BATCH (M+22) =  GXYZX
               BATCH (M+23) =  GYYZX
               BATCH (M+24) =  GYZZX
               BATCH (M+25) =  GXZZX
               BATCH (M+26) =  GYZZX
               BATCH (M+27) =  GZZZX
               BATCH (M+28) =  GXXXY
               BATCH (M+29) =  GXYXY
               BATCH (M+30) =  GXZXY
               BATCH (M+31) =  GXYXY
               BATCH (M+32) =  GYYXY
               BATCH (M+33) =  GYZXY
               BATCH (M+34) =  GXZXY
               BATCH (M+35) =  GYZXY
               BATCH (M+36) =  GZZXY
               BATCH (M+37) =  GXXYY
               BATCH (M+38) =  GXYYY
               BATCH (M+39) =  GXZYY
               BATCH (M+40) =  GXYYY
               BATCH (M+41) =  GYYYY
               BATCH (M+42) =  GYZYY
               BATCH (M+43) =  GXZYY
               BATCH (M+44) =  GYZYY
               BATCH (M+45) =  GZZYY
               BATCH (M+46) =  GXXZY
               BATCH (M+47) =  GXYZY
               BATCH (M+48) =  GXZZY
               BATCH (M+49) =  GXYZY
               BATCH (M+50) =  GYYZY
               BATCH (M+51) =  GYZZY
               BATCH (M+52) =  GXZZY
               BATCH (M+53) =  GYZZY
               BATCH (M+54) =  GZZZY
               BATCH (M+55) =  GXXXZ
               BATCH (M+56) =  GXYXZ
               BATCH (M+57) =  GXZXZ
               BATCH (M+58) =  GXYXZ
               BATCH (M+59) =  GYYXZ
               BATCH (M+60) =  GYZXZ
               BATCH (M+61) =  GXZXZ
               BATCH (M+62) =  GYZXZ
               BATCH (M+63) =  GZZXZ
               BATCH (M+64) =  GXXYZ
               BATCH (M+65) =  GXYYZ
               BATCH (M+66) =  GXZYZ
               BATCH (M+67) =  GXYYZ
               BATCH (M+68) =  GYYYZ
               BATCH (M+69) =  GYZYZ
               BATCH (M+70) =  GXZYZ
               BATCH (M+71) =  GYZYZ
               BATCH (M+72) =  GZZYZ
               BATCH (M+73) =  GXXZZ
               BATCH (M+74) =  GXYZZ
               BATCH (M+75) =  GXZZZ
               BATCH (M+76) =  GXYZZ
               BATCH (M+77) =  GYYZZ
               BATCH (M+78) =  GYZZZ
               BATCH (M+79) =  GXZZZ
               BATCH (M+80) =  GYZZZ
               BATCH (M+81) =  GZZZZ

               M = M + 81

 4400       CONTINUE
 4000    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 4-center (AB|CD)
C                type integrals.
C
C
   50    M = 0
         DO 5000 IJ = 1,MIJ
            PVAL = P (IJ)
            PXVAL = PX (IJ)
            PYVAL = PY (IJ)
            PZVAL = PZ (IJ)
            PSCALE = SCALEP (IJ)
            U3 = HALF / PVAL
            XSPSS1 = PXVAL - X2
            XPSSS1 = PXVAL - X1
            YSPSS1 = PYVAL - Y2
            YPSSS1 = PYVAL - Y1
            ZSPSS1 = PZVAL - Z2
            ZPSSS1 = PZVAL - Z1
            DO 5500 KL = 1,MKL
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

                   F4 = ((((( FTABLE (10,TGRID) * DELTA6
     +                      + FTABLE (9,TGRID)) * DELTA5
     +                      + FTABLE (8,TGRID)) * DELTA4
     +                      + FTABLE (7,TGRID)) * DELTA3
     +                      + FTABLE (6,TGRID)) * DELTA2
     +                      + FTABLE (5,TGRID)) * DELTA1
     +                      + FTABLE (4,TGRID)

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2
                   F3 = SCALE * F3
                   F4 = SCALE * F4

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1
                   F3 = FIVE  * T2INV * F2
                   F4 = SEVEN * T2INV * F3

               END IF

               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               U2 = HALF * PQPINV
               U4 = HALF / QVAL
               U5 = U2 + PQPINV
               U6 = U5 + PQPINV
               U7 = - U0 * U4
               U8 = U1 * U3
C
C
C             ...the X-terms (with exception of XSPSS1 and XPSSS1,
C                which can be evaluated in the P-loop).
C
C
               XSSSP1 = QXVAL - X4
               XSSPS1 = QXVAL - X3
               XSSSP2 = PQX * U0
               XSPSS2 = PQX * U1
               A      = XSSSP1 + XSSPS1
               B      = XSPSS1 + XPSSS1
               C      = A * B
               XSPSP1 = XSPSS1 * XSSSP1
               XPSSP1 = XPSSS1 * XSSSP1
               XSPPS1 = XSPSS1 * XSSPS1
               XPSPS1 = XPSSS1 * XSSPS1
               XSSPP1 = XSSSP1 * XSSPS1 + U4
               XPPSS1 = XSPSS1 * XPSSS1 + U3
               D      = XSSSP1 * XSPSS2 + U2
               E      = XSSPS1 * XSPSS2 + U2
               XSPSP2 = XSPSS1 * XSSSP2 + D
               XPSSP2 = XPSSS1 * XSSSP2 + D
               XSPPS2 = XSPSS1 * XSSSP2 + E
               XPSPS2 = XPSSS1 * XSSSP2 + E
               XSSPP2 = A * XSSSP2 + U7
               XPPSS2 = B * XSPSS2 + U8
               XSPSP3 = XSSSP2 * XSPSS2
               XSSPP3 = XSSSP2 * XSSSP2
               XPPSS3 = XSPSS2 * XSPSS2
               XSPPP1 = XSSPP1 * XSPSS1
               XPSPP1 = XSSPP1 * XPSSS1
               XPPSP1 = XPPSS1 * XSSSP1
               XPPPS1 = XPPSS1 * XSSPS1
               D      = XSPSS2 * XSSPP1 + A * U2
               E      = XSSSP2 * XPPSS1 + B * U2
               XSPPP2 = XSPSS1 * XSSPP2 + D
               XPSPP2 = XPSSS1 * XSSPP2 + D
               XPPSP2 = XSSSP1 * XPPSS2 + E
               XPPPS2 = XSSPS1 * XPPSS2 + E
               D      = A * XSPSP3 + XSSSP2 * U5
               E      = B * XSPSP3 + XSPSS2 * U5
               XSPPP3 = XSPSS1 * XSSPP3 + D
               XPSPP3 = XPSSS1 * XSSPP3 + D
               XPPSP3 = XSSSP1 * XPPSS3 + E
               XPPPS3 = XSSPS1 * XPPSS3 + E
               XSPPP4 = XSSPP3 * XSPSS2
               XPPSP4 = XPPSS3 * XSSSP2
               D      = B * XSSSP2 + A * XSPSS2 + U2
               XPPPP1 = XSSPP1 * XPPSS1
               XPPPP2 = XSSPP2 * XPPSS1 + XPPSS2 * XSSPP1 + C*U2
               XPPPP3 = XPPSS1 * XSSPP3 + XSSPP1 * XPPSS3 + XSPSP3 * C
     +                + U5 * D
               XPPPP4 = XSPSP3 * (D + U6)
               XPPPP5 = XSSPP3 * XPPSS3
C
C
C             ...the Y-terms (with exception of YSPSS1 and YPSSS1,
C                which can be evaluated in the P-loop).
C
C
               YSSSP1 = QYVAL - Y4
               YSSPS1 = QYVAL - Y3
               YSSSP2 = PQY * U0
               YSPSS2 = PQY * U1
               A      = YSSSP1 + YSSPS1
               B      = YSPSS1 + YPSSS1
               C      = A * B
               YSPSP1 = YSPSS1 * YSSSP1
               YPSSP1 = YPSSS1 * YSSSP1
               YSPPS1 = YSPSS1 * YSSPS1
               YPSPS1 = YPSSS1 * YSSPS1
               YSSPP1 = YSSSP1 * YSSPS1 + U4
               YPPSS1 = YSPSS1 * YPSSS1 + U3
               D      = YSSSP1 * YSPSS2 + U2
               E      = YSSPS1 * YSPSS2 + U2
               YSPSP2 = YSPSS1 * YSSSP2 + D
               YPSSP2 = YPSSS1 * YSSSP2 + D
               YSPPS2 = YSPSS1 * YSSSP2 + E
               YPSPS2 = YPSSS1 * YSSSP2 + E
               YSSPP2 = A * YSSSP2 + U7
               YPPSS2 = B * YSPSS2 + U8
               YSPSP3 = YSSSP2 * YSPSS2
               YSSPP3 = YSSSP2 * YSSSP2
               YPPSS3 = YSPSS2 * YSPSS2
               YSPPP1 = YSSPP1 * YSPSS1
               YPSPP1 = YSSPP1 * YPSSS1
               YPPSP1 = YPPSS1 * YSSSP1
               YPPPS1 = YPPSS1 * YSSPS1
               D      = YSPSS2 * YSSPP1 + A * U2
               E      = YSSSP2 * YPPSS1 + B * U2
               YSPPP2 = YSPSS1 * YSSPP2 + D
               YPSPP2 = YPSSS1 * YSSPP2 + D
               YPPSP2 = YSSSP1 * YPPSS2 + E
               YPPPS2 = YSSPS1 * YPPSS2 + E
               D      = A * YSPSP3 + YSSSP2 * U5
               E      = B * YSPSP3 + YSPSS2 * U5
               YSPPP3 = YSPSS1 * YSSPP3 + D
               YPSPP3 = YPSSS1 * YSSPP3 + D
               YPPSP3 = YSSSP1 * YPPSS3 + E
               YPPPS3 = YSSPS1 * YPPSS3 + E
               YSPPP4 = YSSPP3 * YSPSS2
               YPPSP4 = YPPSS3 * YSSSP2
               D      = B * YSSSP2 + A * YSPSS2 + U2
               YPPPP1 = YSSPP1 * YPPSS1
               YPPPP2 = YSSPP2 * YPPSS1 + YPPSS2 * YSSPP1 + C*U2
               YPPPP3 = YPPSS1 * YSSPP3 + YSSPP1 * YPPSS3 + YSPSP3 * C
     +                + U5 * D
               YPPPP4 = YSPSP3 * (D + U6)
               YPPPP5 = YSSPP3 * YPPSS3
C
C
C             ...the Z-terms (with exception of ZSPSS1 and ZPSSS1,
C                which can be evaluated in the P-loop).
C
C
               ZSSSP1 = QZVAL - Z4
               ZSSPS1 = QZVAL - Z3
               ZSSSP2 = PQZ * U0
               ZSPSS2 = PQZ * U1
               A      = ZSSSP1 + ZSSPS1
               B      = ZSPSS1 + ZPSSS1
               C      = A * B
               ZSPSP1 = ZSPSS1 * ZSSSP1
               ZPSSP1 = ZPSSS1 * ZSSSP1
               ZSPPS1 = ZSPSS1 * ZSSPS1
               ZPSPS1 = ZPSSS1 * ZSSPS1
               ZSSPP1 = ZSSSP1 * ZSSPS1 + U4
               ZPPSS1 = ZSPSS1 * ZPSSS1 + U3
               D      = ZSSSP1 * ZSPSS2 + U2
               E      = ZSSPS1 * ZSPSS2 + U2
               ZSPSP2 = ZSPSS1 * ZSSSP2 + D
               ZPSSP2 = ZPSSS1 * ZSSSP2 + D
               ZSPPS2 = ZSPSS1 * ZSSSP2 + E
               ZPSPS2 = ZPSSS1 * ZSSSP2 + E
               ZSSPP2 = A * ZSSSP2 + U7
               ZPPSS2 = B * ZSPSS2 + U8
               ZSPSP3 = ZSSSP2 * ZSPSS2
               ZSSPP3 = ZSSSP2 * ZSSSP2
               ZPPSS3 = ZSPSS2 * ZSPSS2
               ZSPPP1 = ZSSPP1 * ZSPSS1
               ZPSPP1 = ZSSPP1 * ZPSSS1
               ZPPSP1 = ZPPSS1 * ZSSSP1
               ZPPPS1 = ZPPSS1 * ZSSPS1
               D      = ZSPSS2 * ZSSPP1 + A * U2
               E      = ZSSSP2 * ZPPSS1 + B * U2
               ZSPPP2 = ZSPSS1 * ZSSPP2 + D
               ZPSPP2 = ZPSSS1 * ZSSPP2 + D
               ZPPSP2 = ZSSSP1 * ZPPSS2 + E
               ZPPPS2 = ZSSPS1 * ZPPSS2 + E
               D      = A * ZSPSP3 + ZSSSP2 * U5
               E      = B * ZSPSP3 + ZSPSS2 * U5
               ZSPPP3 = ZSPSS1 * ZSSPP3 + D
               ZPSPP3 = ZPSSS1 * ZSSPP3 + D
               ZPPSP3 = ZSSSP1 * ZPPSS3 + E
               ZPPPS3 = ZSSPS1 * ZPPSS3 + E
               ZSPPP4 = ZSSPP3 * ZSPSS2
               ZPPSP4 = ZPPSS3 * ZSSSP2
               D      = B * ZSSSP2 + A * ZSPSS2 + U2
               ZPPPP1 = ZSSPP1 * ZPPSS1
               ZPPPP2 = ZSSPP2 * ZPPSS1 + ZPPSS2 * ZSSPP1 + C*U2
               ZPPPP3 = ZPPSS1 * ZSSPP3 + ZSSPP1 * ZPPSS3 + ZSPSP3 * C
     +                + U5 * D
               ZPPPP4 = ZSPSP3 * (D + U6)
               ZPPPP5 = ZSSPP3 * ZPPSS3
C
C
C             ...assemble the 4-center (AB|CD) type integrals.
C
C
               GXXXX = XPPPP1*F0+XPPPP2*F1+XPPPP3*F2+XPPPP4*F3+XPPPP5*F4
               GYYYY = YPPPP1*F0+YPPPP2*F1+YPPPP3*F2+YPPPP4*F3+YPPPP5*F4
               GZZZZ = ZPPPP1*F0+ZPPPP2*F1+ZPPPP3*F2+ZPPPP4*F3+ZPPPP5*F4

               A = XPPPS1 * F0 + XPPPS2 * F1 + XPPPS3 * F2 + XPPSP4 * F3
               B = XPPPS1 * F1 + XPPPS2 * F2 + XPPPS3 * F3 + XPPSP4 * F4
               C = XPPSP1 * F0 + XPPSP2 * F1 + XPPSP3 * F2 + XPPSP4 * F3
               D = XPPSP1 * F1 + XPPSP2 * F2 + XPPSP3 * F3 + XPPSP4 * F4
               E = XPSPP1 * F0 + XPSPP2 * F1 + XPSPP3 * F2 + XSPPP4 * F3
               F = XPSPP1 * F1 + XPSPP2 * F2 + XPSPP3 * F3 + XSPPP4 * F4
               G = XSPPP1 * F0 + XSPPP2 * F1 + XSPPP3 * F2 + XSPPP4 * F3
               H = XSPPP1 * F1 + XSPPP2 * F2 + XSPPP3 * F3 + XSPPP4 * F4

               GXXXY = A * YSSSP1 + B * YSSSP2
               GXXYX = C * YSSPS1 + D * YSSSP2
               GXYXX = E * YSPSS1 + F * YSPSS2
               GYXXX = G * YPSSS1 + H * YSPSS2
               GXXXZ = A * ZSSSP1 + B * ZSSSP2
               GXXZX = C * ZSSPS1 + D * ZSSSP2
               GXZXX = E * ZSPSS1 + F * ZSPSS2
               GZXXX = G * ZPSSS1 + H * ZSPSS2

               A = YPPPS1 * F0 + YPPPS2 * F1 + YPPPS3 * F2 + YPPSP4 * F3
               B = YPPPS1 * F1 + YPPPS2 * F2 + YPPPS3 * F3 + YPPSP4 * F4
               C = YPPSP1 * F0 + YPPSP2 * F1 + YPPSP3 * F2 + YPPSP4 * F3
               D = YPPSP1 * F1 + YPPSP2 * F2 + YPPSP3 * F3 + YPPSP4 * F4
               E = YPSPP1 * F0 + YPSPP2 * F1 + YPSPP3 * F2 + YSPPP4 * F3
               F = YPSPP1 * F1 + YPSPP2 * F2 + YPSPP3 * F3 + YSPPP4 * F4
               G = YSPPP1 * F0 + YSPPP2 * F1 + YSPPP3 * F2 + YSPPP4 * F3
               H = YSPPP1 * F1 + YSPPP2 * F2 + YSPPP3 * F3 + YSPPP4 * F4

               GYYYX = A * XSSSP1 + B * XSSSP2
               GYYXY = C * XSSPS1 + D * XSSSP2
               GYXYY = E * XSPSS1 + F * XSPSS2
               GXYYY = G * XPSSS1 + H * XSPSS2
               GYYYZ = A * ZSSSP1 + B * ZSSSP2
               GYYZY = C * ZSSPS1 + D * ZSSSP2
               GYZYY = E * ZSPSS1 + F * ZSPSS2
               GZYYY = G * ZPSSS1 + H * ZSPSS2

               A = ZPPPS1 * F0 + ZPPPS2 * F1 + ZPPPS3 * F2 + ZPPSP4 * F3
               B = ZPPPS1 * F1 + ZPPPS2 * F2 + ZPPPS3 * F3 + ZPPSP4 * F4
               C = ZPPSP1 * F0 + ZPPSP2 * F1 + ZPPSP3 * F2 + ZPPSP4 * F3
               D = ZPPSP1 * F1 + ZPPSP2 * F2 + ZPPSP3 * F3 + ZPPSP4 * F4
               E = ZPSPP1 * F0 + ZPSPP2 * F1 + ZPSPP3 * F2 + ZSPPP4 * F3
               F = ZPSPP1 * F1 + ZPSPP2 * F2 + ZPSPP3 * F3 + ZSPPP4 * F4
               G = ZSPPP1 * F0 + ZSPPP2 * F1 + ZSPPP3 * F2 + ZSPPP4 * F3
               H = ZSPPP1 * F1 + ZSPPP2 * F2 + ZSPPP3 * F3 + ZSPPP4 * F4

               GZZZX = A * XSSSP1 + B * XSSSP2
               GZZXZ = C * XSSPS1 + D * XSSSP2
               GZXZZ = E * XSPSS1 + F * XSPSS2
               GXZZZ = G * XPSSS1 + H * XSPSS2
               GZZZY = A * YSSSP1 + B * YSSSP2
               GZZYZ = C * YSSPS1 + D * YSSSP2
               GZYZZ = E * YSPSS1 + F * YSPSS2
               GYZZZ = G * YPSSS1 + H * YSPSS2

               A = XPPSS1 * F0 + XPPSS2 * F1 + XPPSS3 * F2
               B = XPPSS1 * F1 + XPPSS2 * F2 + XPPSS3 * F3
               C = XPPSS1 * F2 + XPPSS2 * F3 + XPPSS3 * F4
               D = YPPSS1 * F0 + YPPSS2 * F1 + YPPSS3 * F2
               E = YPPSS1 * F1 + YPPSS2 * F2 + YPPSS3 * F3
               F = YPPSS1 * F2 + YPPSS2 * F3 + YPPSS3 * F4
               G = ZPPSS1 * F0 + ZPPSS2 * F1 + ZPPSS3 * F2
               H = ZPPSS1 * F1 + ZPPSS2 * F2 + ZPPSS3 * F3
               R = ZPPSS1 * F2 + ZPPSS2 * F3 + ZPPSS3 * F4

               GXXYY = A * YSSPP1 + B * YSSPP2 + C * YSSPP3
               GXXZZ = A * ZSSPP1 + B * ZSSPP2 + C * ZSSPP3
               GYYXX = D * XSSPP1 + E * XSSPP2 + F * XSSPP3
               GYYZZ = D * ZSSPP1 + E * ZSSPP2 + F * ZSSPP3
               GZZXX = G * XSSPP1 + H * XSSPP2 + R * XSSPP3
               GZZYY = G * YSSPP1 + H * YSSPP2 + R * YSSPP3

               GXXYZ = (A * YSSPS1 + B * YSSSP2) * ZSSSP1 +
     +                 (B * YSSPS1 + C * YSSSP2) * ZSSSP2
               GXXZY = (A * ZSSPS1 + B * ZSSSP2) * YSSSP1 +
     +                 (B * ZSSPS1 + C * ZSSSP2) * YSSSP2
               GYYXZ = (D * XSSPS1 + E * XSSSP2) * ZSSSP1 +
     +                 (E * XSSPS1 + F * XSSSP2) * ZSSSP2
               GYYZX = (D * ZSSPS1 + E * ZSSSP2) * XSSSP1 +
     +                 (E * ZSSPS1 + F * ZSSSP2) * XSSSP2
               GZZXY = (G * XSSPS1 + H * XSSSP2) * YSSSP1 +
     +                 (H * XSSPS1 + R * XSSSP2) * YSSSP2
               GZZYX = (G * YSSPS1 + H * YSSSP2) * XSSSP1 +
     +                 (H * YSSPS1 + R * YSSSP2) * XSSSP2

               AA = XSPSP3 * F2
               BB = XSPSP3 * F3
               CC = XSPSP3 * F4
               DD = YSPSP3 * F2
               EE = YSPSP3 * F3
               FF = YSPSP3 * F4
               GG = ZSPSP3 * F2
               HH = ZSPSP3 * F3
               RR = ZSPSP3 * F4

               A = XPSPS1 * F0 + XPSPS2 * F1 + AA
               B = XPSPS1 * F1 + XPSPS2 * F2 + BB
               C = XPSPS1 * F2 + XPSPS2 * F3 + CC
               D = YPSPS1 * F0 + YPSPS2 * F1 + DD
               E = YPSPS1 * F1 + YPSPS2 * F2 + EE
               F = YPSPS1 * F2 + YPSPS2 * F3 + FF
               G = ZPSPS1 * F0 + ZPSPS2 * F1 + GG
               H = ZPSPS1 * F1 + ZPSPS2 * F2 + HH
               R = ZPSPS1 * F2 + ZPSPS2 * F3 + RR

               GXYXY = A * YSPSP1 + B * YSPSP2 + C * YSPSP3
               GXZXZ = A * ZSPSP1 + B * ZSPSP2 + C * ZSPSP3
               GYXYX = D * XSPSP1 + E * XSPSP2 + F * XSPSP3
               GYZYZ = D * ZSPSP1 + E * ZSPSP2 + F * ZSPSP3
               GZXZX = G * XSPSP1 + H * XSPSP2 + R * XSPSP3
               GZYZY = G * YSPSP1 + H * YSPSP2 + R * YSPSP3

               GXYXZ = (A * YSPSS1 + B * YSPSS2) * ZSSSP1 +
     +                 (B * YSPSS1 + C * YSPSS2) * ZSSSP2
               GXZXY = (A * ZSPSS1 + B * ZSPSS2) * YSSSP1 +
     +                 (B * ZSPSS1 + C * ZSPSS2) * YSSSP2
               GYXYZ = (D * XSPSS1 + E * XSPSS2) * ZSSSP1 +
     +                 (E * XSPSS1 + F * XSPSS2) * ZSSSP2
               GYZYX = (D * ZSPSS1 + E * ZSPSS2) * XSSSP1 +
     +                 (E * ZSPSS1 + F * ZSPSS2) * XSSSP2
               GZXZY = (G * XSPSS1 + H * XSPSS2) * YSSSP1 +
     +                 (H * XSPSS1 + R * XSPSS2) * YSSSP2
               GZYZX = (G * YSPSS1 + H * YSPSS2) * XSSSP1 +
     +                 (H * YSPSS1 + R * YSPSS2) * XSSSP2

               A = XPSSP1 * F0 + XPSSP2 * F1 + AA
               B = XPSSP1 * F1 + XPSSP2 * F2 + BB
               C = XPSSP1 * F2 + XPSSP2 * F3 + CC
               D = YPSSP1 * F0 + YPSSP2 * F1 + DD
               E = YPSSP1 * F1 + YPSSP2 * F2 + EE
               F = YPSSP1 * F2 + YPSSP2 * F3 + FF
               G = ZPSSP1 * F0 + ZPSSP2 * F1 + GG
               H = ZPSSP1 * F1 + ZPSSP2 * F2 + HH
               R = ZPSSP1 * F2 + ZPSSP2 * F3 + RR

               GXYYX = A * YSPPS1 + B * YSPPS2 + C * YSPSP3
               GXZZX = A * ZSPPS1 + B * ZSPPS2 + C * ZSPSP3
               GYXXY = D * XSPPS1 + E * XSPPS2 + F * XSPSP3
               GYZZY = D * ZSPPS1 + E * ZSPPS2 + F * ZSPSP3
               GZXXZ = G * XSPPS1 + H * XSPPS2 + R * XSPSP3
               GZYYZ = G * YSPPS1 + H * YSPPS2 + R * YSPSP3

               GXYZX = (A * YSPSS1 + B * YSPSS2) * ZSSPS1 +
     +                 (B * YSPSS1 + C * YSPSS2) * ZSSSP2
               GXZYX = (A * ZSPSS1 + B * ZSPSS2) * YSSPS1 +
     +                 (B * ZSPSS1 + C * ZSPSS2) * YSSSP2
               GYXZY = (D * XSPSS1 + E * XSPSS2) * ZSSPS1 +
     +                 (E * XSPSS1 + F * XSPSS2) * ZSSSP2
               GYZXY = (D * ZSPSS1 + E * ZSPSS2) * XSSPS1 +
     +                 (E * ZSPSS1 + F * ZSPSS2) * XSSSP2
               GZXYZ = (G * XSPSS1 + H * XSPSS2) * YSSPS1 +
     +                 (H * XSPSS1 + R * XSPSS2) * YSSSP2
               GZYXZ = (G * YSPSS1 + H * YSPSS2) * XSSPS1 +
     +                 (H * YSPSS1 + R * YSPSS2) * XSSSP2

               A = XSPPS1 * F0 + XSPPS2 * F1 + AA
               B = XSPPS1 * F1 + XSPPS2 * F2 + BB
               C = XSPPS1 * F2 + XSPPS2 * F3 + CC
               D = YSPPS1 * F0 + YSPPS2 * F1 + DD
               E = YSPPS1 * F1 + YSPPS2 * F2 + EE
               F = YSPPS1 * F2 + YSPPS2 * F3 + FF
               G = ZSPPS1 * F0 + ZSPPS2 * F1 + GG
               H = ZSPPS1 * F1 + ZSPPS2 * F2 + HH
               R = ZSPPS1 * F2 + ZSPPS2 * F3 + RR

               GYXXZ = (A * YPSSS1 + B * YSPSS2) * ZSSSP1 +
     +                 (B * YPSSS1 + C * YSPSS2) * ZSSSP2
               GZXXY = (A * ZPSSS1 + B * ZSPSS2) * YSSSP1 +
     +                 (B * ZPSSS1 + C * ZSPSS2) * YSSSP2
               GXYYZ = (D * XPSSS1 + E * XSPSS2) * ZSSSP1 +
     +                 (E * XPSSS1 + F * XSPSS2) * ZSSSP2
               GZYYX = (D * ZPSSS1 + E * ZSPSS2) * XSSSP1 +
     +                 (E * ZPSSS1 + F * ZSPSS2) * XSSSP2
               GXZZY = (G * XPSSS1 + H * XSPSS2) * YSSSP1 +
     +                 (H * XPSSS1 + R * XSPSS2) * YSSSP2
               GYZZX = (G * YPSSS1 + H * YSPSS2) * XSSSP1 +
     +                 (H * YPSSS1 + R * YSPSS2) * XSSSP2

               A = XSPSP1 * F0 + XSPSP2 * F1 + AA
               B = XSPSP1 * F1 + XSPSP2 * F2 + BB
               C = XSPSP1 * F2 + XSPSP2 * F3 + CC
               D = YSPSP1 * F0 + YSPSP2 * F1 + DD
               E = YSPSP1 * F1 + YSPSP2 * F2 + EE
               F = YSPSP1 * F2 + YSPSP2 * F3 + FF
               G = ZSPSP1 * F0 + ZSPSP2 * F1 + GG
               H = ZSPSP1 * F1 + ZSPSP2 * F2 + HH
               R = ZSPSP1 * F2 + ZSPSP2 * F3 + RR

               GYXZX = (A * YPSSS1 + B * YSPSS2) * ZSSPS1 +
     +                 (B * YPSSS1 + C * YSPSS2) * ZSSSP2
               GZXYX = (A * ZPSSS1 + B * ZSPSS2) * YSSPS1 +
     +                 (B * ZPSSS1 + C * ZSPSS2) * YSSSP2
               GXYZY = (D * XPSSS1 + E * XSPSS2) * ZSSPS1 +
     +                 (E * XPSSS1 + F * XSPSS2) * ZSSSP2
               GZYXY = (D * ZPSSS1 + E * ZSPSS2) * XSSPS1 +
     +                 (E * ZPSSS1 + F * ZSPSS2) * XSSSP2
               GXZYZ = (G * XPSSS1 + H * XSPSS2) * YSSPS1 +
     +                 (H * XPSSS1 + R * XSPSS2) * YSSSP2
               GYZXZ = (G * YPSSS1 + H * YSPSS2) * XSSPS1 +
     +                 (H * YPSSS1 + R * YSPSS2) * XSSSP2

               A = XSSPP1 * F0 + XSSPP2 * F1 + XSSPP3 * F2
               B = XSSPP1 * F1 + XSSPP2 * F2 + XSSPP3 * F3
               C = XSSPP1 * F2 + XSSPP2 * F3 + XSSPP3 * F4
               D = YSSPP1 * F0 + YSSPP2 * F1 + YSSPP3 * F2
               E = YSSPP1 * F1 + YSSPP2 * F2 + YSSPP3 * F3
               F = YSSPP1 * F2 + YSSPP2 * F3 + YSSPP3 * F4
               G = ZSSPP1 * F0 + ZSSPP2 * F1 + ZSSPP3 * F2
               H = ZSSPP1 * F1 + ZSSPP2 * F2 + ZSSPP3 * F3
               R = ZSSPP1 * F2 + ZSSPP2 * F3 + ZSSPP3 * F4

               GYZXX = (A * YPSSS1 + B * YSPSS2) * ZSPSS1 +
     +                 (B * YPSSS1 + C * YSPSS2) * ZSPSS2
               GZYXX = (A * ZPSSS1 + B * ZSPSS2) * YSPSS1 +
     +                 (B * ZPSSS1 + C * ZSPSS2) * YSPSS2
               GXZYY = (D * XPSSS1 + E * XSPSS2) * ZSPSS1 +
     +                 (E * XPSSS1 + F * XSPSS2) * ZSPSS2
               GZXYY = (D * ZPSSS1 + E * ZSPSS2) * XSPSS1 +
     +                 (E * ZPSSS1 + F * ZSPSS2) * XSPSS2
               GXYZZ = (G * XPSSS1 + H * XSPSS2) * YSPSS1 +
     +                 (H * XPSSS1 + R * XSPSS2) * YSPSS2
               GYXZZ = (G * YPSSS1 + H * YSPSS2) * XSPSS1 +
     +                 (H * YPSSS1 + R * YSPSS2) * XSPSS2

               BATCH (M+ 1) =  GXXXX
               BATCH (M+ 2) =  GYXXX
               BATCH (M+ 3) =  GZXXX
               BATCH (M+ 4) =  GXYXX
               BATCH (M+ 5) =  GYYXX
               BATCH (M+ 6) =  GZYXX
               BATCH (M+ 7) =  GXZXX
               BATCH (M+ 8) =  GYZXX
               BATCH (M+ 9) =  GZZXX
               BATCH (M+10) =  GXXYX
               BATCH (M+11) =  GYXYX
               BATCH (M+12) =  GZXYX
               BATCH (M+13) =  GXYYX
               BATCH (M+14) =  GYYYX
               BATCH (M+15) =  GZYYX
               BATCH (M+16) =  GXZYX
               BATCH (M+17) =  GYZYX
               BATCH (M+18) =  GZZYX
               BATCH (M+19) =  GXXZX
               BATCH (M+20) =  GYXZX
               BATCH (M+21) =  GZXZX
               BATCH (M+22) =  GXYZX
               BATCH (M+23) =  GYYZX
               BATCH (M+24) =  GZYZX
               BATCH (M+25) =  GXZZX
               BATCH (M+26) =  GYZZX
               BATCH (M+27) =  GZZZX
               BATCH (M+28) =  GXXXY
               BATCH (M+29) =  GYXXY
               BATCH (M+30) =  GZXXY
               BATCH (M+31) =  GXYXY
               BATCH (M+32) =  GYYXY
               BATCH (M+33) =  GZYXY
               BATCH (M+34) =  GXZXY
               BATCH (M+35) =  GYZXY
               BATCH (M+36) =  GZZXY
               BATCH (M+37) =  GXXYY
               BATCH (M+38) =  GYXYY
               BATCH (M+39) =  GZXYY
               BATCH (M+40) =  GXYYY
               BATCH (M+41) =  GYYYY
               BATCH (M+42) =  GZYYY
               BATCH (M+43) =  GXZYY
               BATCH (M+44) =  GYZYY
               BATCH (M+45) =  GZZYY
               BATCH (M+46) =  GXXZY
               BATCH (M+47) =  GYXZY
               BATCH (M+48) =  GZXZY
               BATCH (M+49) =  GXYZY
               BATCH (M+50) =  GYYZY
               BATCH (M+51) =  GZYZY
               BATCH (M+52) =  GXZZY
               BATCH (M+53) =  GYZZY
               BATCH (M+54) =  GZZZY
               BATCH (M+55) =  GXXXZ
               BATCH (M+56) =  GYXXZ
               BATCH (M+57) =  GZXXZ
               BATCH (M+58) =  GXYXZ
               BATCH (M+59) =  GYYXZ
               BATCH (M+60) =  GZYXZ
               BATCH (M+61) =  GXZXZ
               BATCH (M+62) =  GYZXZ
               BATCH (M+63) =  GZZXZ
               BATCH (M+64) =  GXXYZ
               BATCH (M+65) =  GYXYZ
               BATCH (M+66) =  GZXYZ
               BATCH (M+67) =  GXYYZ
               BATCH (M+68) =  GYYYZ
               BATCH (M+69) =  GZYYZ
               BATCH (M+70) =  GXZYZ
               BATCH (M+71) =  GYZYZ
               BATCH (M+72) =  GZZYZ
               BATCH (M+73) =  GXXZZ
               BATCH (M+74) =  GYXZZ
               BATCH (M+75) =  GZXZZ
               BATCH (M+76) =  GXYZZ
               BATCH (M+77) =  GYYZZ
               BATCH (M+78) =  GZYZZ
               BATCH (M+79) =  GXZZZ
               BATCH (M+80) =  GYZZZ
               BATCH (M+81) =  GZZZZ

               M = M + 81

 5500       CONTINUE
 5000    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
