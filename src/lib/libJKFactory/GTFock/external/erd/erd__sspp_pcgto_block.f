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
         SUBROUTINE  ERD__SSPP_PCGTO_BLOCK
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
C  OPERATION   : ERD__SSPP_PCGTO_BLOCK
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
C                the total number of primitive functions (here = 9)
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
C                      sspp , spsp , pssp , spps , psps , ppss
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
C                                    sspp/spsp/pssp/spps/psps/ppss
C                                    integral batch
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
C                                    cartesian sspp/spsp/pssp/spps/
C                                    psps/ppss integrals
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
         INTEGER     SHELL1,SHELL3,SHELLP

         INTEGER     PRIM1 (1:MIJ)
         INTEGER     PRIM2 (1:MIJ)
         INTEGER     PRIM3 (1:MKL)
         INTEGER     PRIM4 (1:MKL)

         DOUBLE PRECISION  A,B,C,D,E,F,G
         DOUBLE PRECISION  DELTA1,DELTA2,DELTA3,DELTA4,DELTA5,DELTA6
         DOUBLE PRECISION  EXP1,EXP2,EXP3,EXP4
         DOUBLE PRECISION  F0,F1,F2
         DOUBLE PRECISION  PVAL,QVAL
         DOUBLE PRECISION  PQPLUS,PQMULT,PQPINV
         DOUBLE PRECISION  PQX,PQY,PQZ
         DOUBLE PRECISION  PSCALE
         DOUBLE PRECISION  PXSUB,PYSUB,PZSUB,QXSUB,QYSUB,QZSUB
         DOUBLE PRECISION  PXVAL,PYVAL,PZVAL,QXVAL,QYVAL,QZVAL
         DOUBLE PRECISION  RNPQSQ
         DOUBLE PRECISION  SCALE
         DOUBLE PRECISION  T,TINV,T2INV,TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION  U0,U1
         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4
         DOUBLE PRECISION  X12,Y12,Z12,X34,Y34,Z34
         DOUBLE PRECISION  XSP1,XPS1,XSP2,XPS2,
     +                     YSP1,YPS1,YSP2,YPS2,
     +                     ZSP1,ZPS1,ZSP2,ZPS2
         DOUBLE PRECISION  ZERO,SIXTH,FIFTH,THIRD,HALF,ONE,THREE,PI

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
         PARAMETER  (SIXTH   =  0.166666666666667D0)
         PARAMETER  (FIFTH   =  0.2D0)
         PARAMETER  (THIRD   =  0.333333333333333D0)
         PARAMETER  (HALF    =  0.5D0)
         PARAMETER  (ONE     =  1.D0 )
         PARAMETER  (THREE   =  3.D0 )
         PARAMETER  (PI      =  3.141592653589793D0)
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
C             ...jump according to where the s-type functions are
C                located and the type of 'K4' loop:
C
C                SHELLP CASE  |        Integral center type
C               --------------|------------------------------------
C                   0     1   |  (AA|AA)  atomic     sspp
C                   0     2   |  (AA|CC)  2-center   sspp
C                   0     3   |  (AB|CC)  3-center   sspp
C                   0     4   |  (AA|CD)  3-center   sspp
C                   0     5   |  (AB|CD)  4-center   sspp
C
C                   1     1   |  (AA|AA)  atomic     spsp,spps,pssp,psps
C                   1     2   |  (AA|CC)  2-center   spsp,spps,pssp,psps
C                   1     3   |  (AB|CC)  3-center   spsp,spps,pssp,psps
C                   1     4   |  (AA|CD)  3-center   spsp,spps,pssp,psps
C                   1     5   |  (AB|CD)  4-center   spsp,spps,pssp,psps
C
C                   2     1   |  (AA|AA)  atomic     ppss
C                   2     2   |  (AA|CC)  2-center   ppss
C                   2     3   |  (AB|CC)  3-center   ppss
C                   2     4   |  (AA|CD)  3-center   ppss
C                   2     5   |  (AB|CD)  4-center   ppss
C
C
         GOTO (11,12,13,14,15,
     +         21,22,23,24,25,
     +         31,32,33,34,35) 5*SHELLP + CASE
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for atomic sspp (AA|AA)
C                type integrals.
C
C
   11    M = 0
         DO 1100 IJ = 1,MIJ
            PVAL = P (IJ)
            PSCALE = SCALEP (IJ)
            DO 1110 KL = 1,MKL
               QVAL = Q (KL)
               PQPLUS = PVAL + QVAL
               F0 = PSCALE * SCALEQ (KL) / (PVAL*QVAL*DSQRT (PQPLUS))
               U1 = (F0 - (PVAL / PQPLUS) * F0 * THIRD) / (QVAL + QVAL)
               BATCH (M+1) = U1
               BATCH (M+2) = ZERO
               BATCH (M+3) = ZERO
               BATCH (M+4) = ZERO
               BATCH (M+5) = U1
               BATCH (M+6) = ZERO
               BATCH (M+7) = ZERO
               BATCH (M+8) = ZERO
               BATCH (M+9) = U1
               M = M + 9
 1110       CONTINUE
 1100    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 2-center sspp (AA|CC)
C                type integrals.
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

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1

               END IF

               U0 = PVAL * PQPINV
               U1 = (F0 - U0 * F1) / (QVAL + QVAL)
               XSP2 = PQX * U0
               YSP2 = PQY * U0
               ZSP2 = PQZ * U0
               A = YSP2 * F2
               B = ZSP2 * F2
               BATCH (M+1) = XSP2 * XSP2 * F2 + U1
               BATCH (M+2) = XSP2 * A
               BATCH (M+3) = XSP2 * B
               BATCH (M+4) = BATCH (M+2)
               BATCH (M+5) = YSP2 * A + U1
               BATCH (M+6) = YSP2 * B
               BATCH (M+7) = BATCH (M+3)
               BATCH (M+8) = BATCH (M+6)
               BATCH (M+9) = ZSP2 * B + U1
               M = M + 9
 1220       CONTINUE
 1200    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 3-center sspp (AB|CC)
C                type integrals.
C
C
   13    M = 0
         DO 1300 IJ = 1,MIJ
            PVAL = P (IJ)
            PQX = PX (IJ) - X3
            PQY = PY (IJ) - Y3
            PQZ = PZ (IJ) - Z3
            RNPQSQ = PQX * PQX + PQY * PQY + PQZ * PQZ
            PSCALE = SCALEP (IJ)
            DO 1330 KL = 1,MKL
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

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1

               END IF

               U0 = PVAL * PQPINV
               U1 = (F0 - U0 * F1) / (QVAL + QVAL)
               XSP2 = PQX * U0
               YSP2 = PQY * U0
               ZSP2 = PQZ * U0
               A = YSP2 * F2
               B = ZSP2 * F2
               BATCH (M+1) = XSP2 * XSP2 * F2 + U1
               BATCH (M+2) = XSP2 * A
               BATCH (M+3) = XSP2 * B
               BATCH (M+4) = BATCH (M+2)
               BATCH (M+5) = YSP2 * A + U1
               BATCH (M+6) = YSP2 * B
               BATCH (M+7) = BATCH (M+3)
               BATCH (M+8) = BATCH (M+6)
               BATCH (M+9) = ZSP2 * B + U1
               M = M + 9
 1330       CONTINUE
 1300    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 3-center sspp (AA|CD)
C                type integrals.
C
C
   14    M = 0
         DO 1400 IJ = 1,MIJ
            PVAL = P (IJ)
            PSCALE = SCALEP (IJ)
            DO 1440 KL = 1,MKL
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

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1

               END IF

               U0 = PVAL * PQPINV
               U1 = (F0 - U0 * F1) / (QVAL + QVAL)
               XPS1 = QXVAL - X3
               YPS1 = QYVAL - Y3
               ZPS1 = QZVAL - Z3
               XSP1 = QXVAL - X4
               YSP1 = QYVAL - Y4
               ZSP1 = QZVAL - Z4
               XSP2 = PQX * U0
               YSP2 = PQY * U0
               ZSP2 = PQZ * U0
               A = XPS1 * F0 +  XSP2 * F1
               B = XPS1 * F1 +  XSP2 * F2
               C = YPS1 * F0 +  YSP2 * F1
               D = YPS1 * F1 +  YSP2 * F2
               E = ZPS1 * F0 +  ZSP2 * F1
               F = ZPS1 * F1 +  ZSP2 * F2
               BATCH (M+1) = XSP1 * A + XSP2 * B + U1
               BATCH (M+2) = XSP1 * C + XSP2 * D
               BATCH (M+3) = XSP1 * E + XSP2 * F
               BATCH (M+4) = YSP1 * A + YSP2 * B
               BATCH (M+5) = YSP1 * C + YSP2 * D + U1
               BATCH (M+6) = YSP1 * E + YSP2 * F
               BATCH (M+7) = ZSP1 * A + ZSP2 * B
               BATCH (M+8) = ZSP1 * C + ZSP2 * D
               BATCH (M+9) = ZSP1 * E + ZSP2 * F + U1
               M = M + 9
 1440       CONTINUE
 1400    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 4-center sspp (AB|CD)
C                type integrals.
C
C
   15    M = 0
         DO 1500 IJ = 1,MIJ
            PVAL = P (IJ)
            PXVAL = PX (IJ)
            PYVAL = PY (IJ)
            PZVAL = PZ (IJ)
            PSCALE = SCALEP (IJ)
            DO 1550 KL = 1,MKL
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

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1

               END IF

               U0 = PVAL * PQPINV
               U1 = (F0 - U0 * F1) / (QVAL + QVAL)
               XPS1 = QXVAL - X3
               YPS1 = QYVAL - Y3
               ZPS1 = QZVAL - Z3
               XSP1 = QXVAL - X4
               YSP1 = QYVAL - Y4
               ZSP1 = QZVAL - Z4
               XSP2 = PQX * U0
               YSP2 = PQY * U0
               ZSP2 = PQZ * U0
               A = XPS1 * F0 +  XSP2 * F1
               B = XPS1 * F1 +  XSP2 * F2
               C = YPS1 * F0 +  YSP2 * F1
               D = YPS1 * F1 +  YSP2 * F2
               E = ZPS1 * F0 +  ZSP2 * F1
               F = ZPS1 * F1 +  ZSP2 * F2
               BATCH (M+1) = XSP1 * A + XSP2 * B + U1
               BATCH (M+2) = XSP1 * C + XSP2 * D
               BATCH (M+3) = XSP1 * E + XSP2 * F
               BATCH (M+4) = YSP1 * A + YSP2 * B
               BATCH (M+5) = YSP1 * C + YSP2 * D + U1
               BATCH (M+6) = YSP1 * E + YSP2 * F
               BATCH (M+7) = ZSP1 * A + ZSP2 * B
               BATCH (M+8) = ZSP1 * C + ZSP2 * D
               BATCH (M+9) = ZSP1 * E + ZSP2 * F + U1
               M = M + 9
 1550       CONTINUE
 1500    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for atomic spsp,spps,
C                pssp and psps (AA|AA) type integrals.
C
C
   21    M = 0
         DO 2100 IJ = 1,MIJ
            PVAL = P (IJ)
            PSCALE = SCALEP (IJ)
            DO 2110 KL = 1,MKL
               QVAL = Q (KL)
               PQPLUS = PVAL + QVAL
               PQPINV = ONE / PQPLUS
               F0 = PSCALE * SCALEQ (KL) / (PVAL*QVAL*DSQRT (PQPLUS))
               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               D = SIXTH * F0 * PQPINV
               BATCH (M+1) = D
               BATCH (M+2) = ZERO
               BATCH (M+3) = ZERO
               BATCH (M+4) = ZERO
               BATCH (M+5) = D
               BATCH (M+6) = ZERO
               BATCH (M+7) = ZERO
               BATCH (M+8) = ZERO
               BATCH (M+9) = D
               M = M + 9
 2110       CONTINUE
 2100    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 2-center spsp,spps,
C                pssp and psps (AA|CC) type integrals.
C
C
   22    PQX = X1 - X3
         PQY = Y1 - Y3
         PQZ = Z1 - Z3
         RNPQSQ = PQX * PQX + PQY * PQY + PQZ * PQZ

         M = 0
         DO 2200 IJ = 1,MIJ
            PVAL = P (IJ)
            PSCALE = SCALEP (IJ)
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

                   F1 = SCALE * F1
                   F2 = SCALE * F2

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F1 = T2INV * SCALE * HALF * DSQRT (PI*TINV)
                   F2 = THREE * T2INV * F1

               END IF

               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               XSP2 = PQX * U0
               YSP2 = PQY * U0
               ZSP2 = PQZ * U0
               XPS2 = PQX * U1
               YPS2 = PQY * U1
               ZPS2 = PQZ * U1
               A = XPS2 * F2
               B = YPS2 * F2
               C = ZPS2 * F2
               D = HALF * F1 * PQPINV
               BATCH (M+1) = XSP2 * A + D
               BATCH (M+2) = XSP2 * B
               BATCH (M+3) = XSP2 * C
               BATCH (M+4) = YSP2 * A
               BATCH (M+5) = YSP2 * B + D
               BATCH (M+6) = YSP2 * C
               BATCH (M+7) = ZSP2 * A
               BATCH (M+8) = ZSP2 * B
               BATCH (M+9) = ZSP2 * C + D
               M = M + 9
 2220       CONTINUE
 2200    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 3-center spsp,spps,
C                pssp and psps (AB|CC) type integrals.
C
C
   23    IF (SHELL1.EQ.1) THEN
             PXSUB = X1
             PYSUB = Y1
             PZSUB = Z1
         ELSE
             PXSUB = X2
             PYSUB = Y2
             PZSUB = Z2
         END IF

         M = 0
         DO 2300 IJ = 1,MIJ
            PVAL = P (IJ)
            PXVAL = PX (IJ)
            PYVAL = PY (IJ)
            PZVAL = PZ (IJ)
            PQX = PXVAL - X3
            PQY = PYVAL - Y3
            PQZ = PZVAL - Z3
            RNPQSQ = PQX * PQX + PQY * PQY + PQZ * PQZ
            PSCALE = SCALEP (IJ)
            XPS1 = PXVAL - PXSUB
            YPS1 = PYVAL - PYSUB
            ZPS1 = PZVAL - PZSUB
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

                   F1 = SCALE * F1
                   F2 = SCALE * F2

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F1 = T2INV * SCALE * HALF * DSQRT (PI*TINV)
                   F2 = THREE * T2INV * F1

               END IF

               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               XSP2 = PQX * U0
               YSP2 = PQY * U0
               ZSP2 = PQZ * U0
               XPS2 = PQX * U1
               YPS2 = PQY * U1
               ZPS2 = PQZ * U1
               A = XPS1 * F1 + XPS2 * F2
               B = YPS1 * F1 + YPS2 * F2
               C = ZPS1 * F1 + ZPS2 * F2
               D = HALF * F1 * PQPINV
               BATCH (M+1) = XSP2 * A + D
               BATCH (M+2) = XSP2 * B
               BATCH (M+3) = XSP2 * C
               BATCH (M+4) = YSP2 * A
               BATCH (M+5) = YSP2 * B + D
               BATCH (M+6) = YSP2 * C
               BATCH (M+7) = ZSP2 * A
               BATCH (M+8) = ZSP2 * B
               BATCH (M+9) = ZSP2 * C + D
               M = M + 9
 2330       CONTINUE
 2300    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 3-center spsp,spps,
C                pssp and psps (AA|CD) type integrals.
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
            PSCALE = SCALEP (IJ)
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

                   F1 = SCALE * F1
                   F2 = SCALE * F2

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F1 = T2INV * SCALE * HALF * DSQRT (PI*TINV)
                   F2 = THREE * T2INV * F1

               END IF

               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               XSP1 = QXVAL - QXSUB
               YSP1 = QYVAL - QYSUB
               ZSP1 = QZVAL - QZSUB
               XSP2 = PQX * U0
               YSP2 = PQY * U0
               ZSP2 = PQZ * U0
               XPS2 = PQX * U1
               YPS2 = PQY * U1
               ZPS2 = PQZ * U1
               A = XSP1 * F1 + XSP2 * F2
               B = YSP1 * F1 + YSP2 * F2
               C = ZSP1 * F1 + ZSP2 * F2
               D = HALF * F1 * PQPINV
               BATCH (M+1) = XPS2 * A + D
               BATCH (M+2) = YPS2 * A
               BATCH (M+3) = ZPS2 * A
               BATCH (M+4) = XPS2 * B
               BATCH (M+5) = YPS2 * B + D
               BATCH (M+6) = ZPS2 * B
               BATCH (M+7) = XPS2 * C
               BATCH (M+8) = YPS2 * C
               BATCH (M+9) = ZPS2 * C + D
               M = M + 9
 2440       CONTINUE
 2400    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 4-center spsp,spps,
C                pssp and psps (AB|CD) type integrals.
C
C
   25    IF (SHELL1.EQ.1) THEN
             PXSUB = X1
             PYSUB = Y1
             PZSUB = Z1
         ELSE
             PXSUB = X2
             PYSUB = Y2
             PZSUB = Z2
         END IF

         IF (SHELL3.EQ.1) THEN
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
            PXVAL = PX (IJ)
            PYVAL = PY (IJ)
            PZVAL = PZ (IJ)
            PSCALE = SCALEP (IJ)
            XPS1 = PXVAL - PXSUB
            YPS1 = PYVAL - PYSUB
            ZPS1 = PZVAL - PZSUB
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

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1

               END IF

               U0 = PVAL * PQPINV
               U1 = - QVAL * PQPINV
               XSP1 = QXVAL - QXSUB
               YSP1 = QYVAL - QYSUB
               ZSP1 = QZVAL - QZSUB
               XSP2 = PQX * U0
               YSP2 = PQY * U0
               ZSP2 = PQZ * U0
               XPS2 = PQX * U1
               YPS2 = PQY * U1
               ZPS2 = PQZ * U1
               A = XPS1 * F0 + XPS2 * F1
               B = XPS1 * F1 + XPS2 * F2
               C = YPS1 * F0 + YPS2 * F1
               D = YPS1 * F1 + YPS2 * F2
               E = ZPS1 * F0 + ZPS2 * F1
               F = ZPS1 * F1 + ZPS2 * F2
               G = HALF * F1 * PQPINV
               BATCH (M+1) = XSP1 * A + XSP2 * B + G
               BATCH (M+2) = XSP1 * C + XSP2 * D
               BATCH (M+3) = XSP1 * E + XSP2 * F
               BATCH (M+4) = YSP1 * A + YSP2 * B
               BATCH (M+5) = YSP1 * C + YSP2 * D + G
               BATCH (M+6) = YSP1 * E + YSP2 * F
               BATCH (M+7) = ZSP1 * A + ZSP2 * B
               BATCH (M+8) = ZSP1 * C + ZSP2 * D
               BATCH (M+9) = ZSP1 * E + ZSP2 * F + G
               M = M + 9
 2550       CONTINUE
 2500    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for atomic ppss (AA|AA)
C                type integrals.
C
C
   31    M = 0
         DO 3100 IJ = 1,MIJ
            PVAL = P (IJ)
            PSCALE = SCALEP (IJ)
            U1 = ONE / (PVAL + PVAL)
            DO 3110 KL = 1,MKL
               QVAL = Q (KL)
               PQPLUS = PVAL + QVAL
               F0 = PSCALE * SCALEQ (KL) / (PVAL*QVAL*DSQRT (PQPLUS))
               C = U1 * (F0 - (QVAL / PQPLUS) * F0 * THIRD)
               BATCH (M+1) = C
               BATCH (M+2) = ZERO
               BATCH (M+3) = ZERO
               BATCH (M+4) = ZERO
               BATCH (M+5) = C
               BATCH (M+6) = ZERO
               BATCH (M+7) = ZERO
               BATCH (M+8) = ZERO
               BATCH (M+9) = C
               M = M + 9
 3110       CONTINUE
 3100    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 2-center ppss (AA|CC)
C                type integrals.
C
C
   32    PQX = X1 - X3
         PQY = Y1 - Y3
         PQZ = Z1 - Z3
         RNPQSQ = PQX * PQX + PQY * PQY + PQZ * PQZ

         M = 0
         DO 3200 IJ = 1,MIJ
            PVAL = P (IJ)
            PSCALE = SCALEP (IJ)
            U1 = ONE / (PVAL + PVAL)
            DO 3220 KL = 1,MKL
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

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1

               END IF

               U0 = - QVAL * PQPINV
               XSP2 = PQX * U0
               YSP2 = PQY * U0
               ZSP2 = PQZ * U0
               A = YSP2 * F2
               B = ZSP2 * F2
               C = U1 * (F0 + U0 * F1)
               BATCH (M+1) = XSP2 * XSP2 * F2 + C
               BATCH (M+2) = XSP2 * A
               BATCH (M+3) = XSP2 * B
               BATCH (M+4) = BATCH (M+2)
               BATCH (M+5) = YSP2 * A + C
               BATCH (M+6) = YSP2 * B
               BATCH (M+7) = BATCH (M+3)
               BATCH (M+8) = BATCH (M+6)
               BATCH (M+9) = ZSP2 * B + C
               M = M + 9
 3220       CONTINUE
 3200    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 3-center ppss (AB|CC)
C                type integrals.
C
C
   33    M = 0
         DO 3300 IJ = 1,MIJ
            PVAL = P (IJ)
            PXVAL = PX (IJ)
            PYVAL = PY (IJ)
            PZVAL = PZ (IJ)
            PQX = PXVAL - X3
            PQY = PYVAL - Y3
            PQZ = PZVAL - Z3
            RNPQSQ = PQX * PQX + PQY * PQY + PQZ * PQZ
            PSCALE = SCALEP (IJ)
            XPS1 = PXVAL - X1
            YPS1 = PYVAL - Y1
            ZPS1 = PZVAL - Z1
            XSP1 = PXVAL - X2
            YSP1 = PYVAL - Y2
            ZSP1 = PZVAL - Z2
            U1 = ONE / (PVAL + PVAL)
            DO 3330 KL = 1,MKL
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

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1

               END IF

               U0 = - QVAL * PQPINV
               XSP2 = PQX * U0
               YSP2 = PQY * U0
               ZSP2 = PQZ * U0
               A = XPS1 * F0 +  XSP2 * F1
               B = XPS1 * F1 +  XSP2 * F2
               C = YPS1 * F0 +  YSP2 * F1
               D = YPS1 * F1 +  YSP2 * F2
               E = ZPS1 * F0 +  ZSP2 * F1
               F = ZPS1 * F1 +  ZSP2 * F2
               G = U1 * (F0 + U0 * F1)
               BATCH (M+1) = XSP1 * A + XSP2 * B + G
               BATCH (M+2) = XSP1 * C + XSP2 * D
               BATCH (M+3) = XSP1 * E + XSP2 * F
               BATCH (M+4) = YSP1 * A + YSP2 * B
               BATCH (M+5) = YSP1 * C + YSP2 * D + G
               BATCH (M+6) = YSP1 * E + YSP2 * F
               BATCH (M+7) = ZSP1 * A + ZSP2 * B
               BATCH (M+8) = ZSP1 * C + ZSP2 * D
               BATCH (M+9) = ZSP1 * E + ZSP2 * F + G
               M = M + 9
 3330       CONTINUE
 3300    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 3-center ppss (AA|CD)
C                type integrals.
C
C
   34    M = 0
         DO 3400 IJ = 1,MIJ
            PVAL = P (IJ)
            PSCALE = SCALEP (IJ)
            U1 = ONE / (PVAL + PVAL)
            DO 3440 KL = 1,MKL
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

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1

               END IF

               U0 = - QVAL * PQPINV
               XSP2 = PQX * U0
               YSP2 = PQY * U0
               ZSP2 = PQZ * U0
               A = YSP2 * F2
               B = ZSP2 * F2
               C = U1 * (F0 + U0 * F1)
               BATCH (M+1) = XSP2 * XSP2 * F2 + C
               BATCH (M+2) = XSP2 * A
               BATCH (M+3) = XSP2 * B
               BATCH (M+4) = BATCH (M+2)
               BATCH (M+5) = YSP2 * A + C
               BATCH (M+6) = YSP2 * B
               BATCH (M+7) = BATCH (M+3)
               BATCH (M+8) = BATCH (M+6)
               BATCH (M+9) = ZSP2 * B + C
               M = M + 9
 3440       CONTINUE
 3400    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 4-center ppss (AB|CD)
C                type integrals.
C
C
   35    M = 0
         DO 3500 IJ = 1,MIJ
            PVAL = P (IJ)
            PXVAL = PX (IJ)
            PYVAL = PY (IJ)
            PZVAL = PZ (IJ)
            PSCALE = SCALEP (IJ)
            XPS1 = PXVAL - X1
            YPS1 = PYVAL - Y1
            ZPS1 = PZVAL - Z1
            XSP1 = PXVAL - X2
            YSP1 = PYVAL - Y2
            ZSP1 = PZVAL - Z2
            U1 = ONE / (PVAL + PVAL)
            DO 3550 KL = 1,MKL
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

                   F0 = SCALE * F0
                   F1 = SCALE * F1
                   F2 = SCALE * F2

               ELSE

                   TINV = ONE / T
                   T2INV = HALF * TINV
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = T2INV * F0
                   F2 = THREE * T2INV * F1

               END IF

               U0 = - QVAL * PQPINV
               XSP2 = PQX * U0
               YSP2 = PQY * U0
               ZSP2 = PQZ * U0
               A = XPS1 * F0 +  XSP2 * F1
               B = XPS1 * F1 +  XSP2 * F2
               C = YPS1 * F0 +  YSP2 * F1
               D = YPS1 * F1 +  YSP2 * F2
               E = ZPS1 * F0 +  ZSP2 * F1
               F = ZPS1 * F1 +  ZSP2 * F2
               G = U1 * (F0 + U0 * F1)
               BATCH (M+1) = XSP1 * A + XSP2 * B + G
               BATCH (M+2) = XSP1 * C + XSP2 * D
               BATCH (M+3) = XSP1 * E + XSP2 * F
               BATCH (M+4) = YSP1 * A + YSP2 * B
               BATCH (M+5) = YSP1 * C + YSP2 * D + G
               BATCH (M+6) = YSP1 * E + YSP2 * F
               BATCH (M+7) = ZSP1 * A + ZSP2 * B
               BATCH (M+8) = ZSP1 * C + ZSP2 * D
               BATCH (M+9) = ZSP1 * E + ZSP2 * F + G
               M = M + 9
 3550       CONTINUE
 3500    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
