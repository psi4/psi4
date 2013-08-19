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
         SUBROUTINE  ERD__SSSP_PCGTO_BLOCK
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
C  OPERATION   : ERD__SSSP_PCGTO_BLOCK
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
C                the total number of primitive functions (here = 3)
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
C                        sssp , ssps , spss , psss
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
C                                    sssp/ssps/spss/psss integral batch
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
C                                    cartesian sssp/ssps/spss/psss
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

         DOUBLE PRECISION  DELTA1,DELTA2,DELTA3,DELTA4,DELTA5,DELTA6
         DOUBLE PRECISION  EXP1,EXP2,EXP3,EXP4
         DOUBLE PRECISION  F0,F1
         DOUBLE PRECISION  PVAL,QVAL
         DOUBLE PRECISION  PQPLUS,PQMULT,PQPINV
         DOUBLE PRECISION  PQX,PQY,PQZ
         DOUBLE PRECISION  PSCALE
         DOUBLE PRECISION  PXSUB,PYSUB,PZSUB,QXSUB,QYSUB,QZSUB
         DOUBLE PRECISION  PXVAL,PYVAL,PZVAL,QXVAL,QYVAL,QZVAL
         DOUBLE PRECISION  RNPQSQ
         DOUBLE PRECISION  SCALE
         DOUBLE PRECISION  T,TINV,TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4
         DOUBLE PRECISION  X12,Y12,Z12,X34,Y34,Z34
         DOUBLE PRECISION  XP,YP,ZP
         DOUBLE PRECISION  FIFTH,FOURTH,THIRD,HALF,ONE,PI

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
         PARAMETER  (FOURTH  =  0.25D0)
         PARAMETER  (THIRD   =  0.333333333333333D0)
         PARAMETER  (HALF    =  0.5D0)
         PARAMETER  (ONE     =  1.D0 )
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
C             ...jump according to where the p-type function is
C                located and the type of 'K4' loop:
C
C                 SHELLP CASE  |        Integral center type
C                --------------|------------------------------------
C                    0     1   |  (AA|AA)  atomic     sssp and ssps
C                    0     2   |  (AA|CC)  2-center   sssp and ssps
C                    0     3   |  (AB|CC)  3-center   sssp and ssps
C                    0     4   |  (AA|CD)  3-center   sssp and ssps
C                    0     5   |  (AB|CD)  4-center   sssp and ssps
C
C                    1     1   |  (AA|AA)  atomic     spss and psss
C                    1     2   |  (AA|CC)  2-center   spss and psss
C                    1     3   |  (AB|CC)  3-center   spss and psss
C                    1     4   |  (AA|CD)  3-center   spss and psss
C                    1     5   |  (AB|CD)  4-center   spss and psss
C
C
         GOTO (11,12,13,14,15,21,22,23,24,25) 5 * SHELLP + CASE
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for atomic sssp and
C                ssps integrals. All these integrals are zero by
C                symmetry.
C
C
   11    RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 2-center sssp and
C                ssps (AA|CC) type integrals.
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

                   F1 = ((((( FTABLE (7,TGRID)  * DELTA6
     +                      + FTABLE (6,TGRID)) * DELTA5
     +                      + FTABLE (5,TGRID)) * DELTA4
     +                      + FTABLE (4,TGRID)) * DELTA3
     +                      + FTABLE (3,TGRID)) * DELTA2
     +                      + FTABLE (2,TGRID)) * DELTA1
     +                      + FTABLE (1,TGRID)
               ELSE
                   TINV = ONE / T
                   F1 = TINV * FOURTH * DSQRT (PI*TINV)
               END IF

               F1 = F1 * PVAL * SCALE * PQPINV

               BATCH (M+1) = PQX * F1
               BATCH (M+2) = PQY * F1
               BATCH (M+3) = PQZ * F1
               M = M + 3
 1220       CONTINUE
 1200    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 3-center sssp and
C                ssps (AB|CC) type integrals.
C
C
   13    M = 0
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

                   F1 = ((((( FTABLE (7,TGRID)  * DELTA6
     +                      + FTABLE (6,TGRID)) * DELTA5
     +                      + FTABLE (5,TGRID)) * DELTA4
     +                      + FTABLE (4,TGRID)) * DELTA3
     +                      + FTABLE (3,TGRID)) * DELTA2
     +                      + FTABLE (2,TGRID)) * DELTA1
     +                      + FTABLE (1,TGRID)
               ELSE
                   TINV = ONE / T
                   F1 = TINV * FOURTH * DSQRT (PI*TINV)
               END IF

               F1 = F1 * PVAL * SCALE * PQPINV

               BATCH (M+1) = PQX * F1
               BATCH (M+2) = PQY * F1
               BATCH (M+3) = PQZ * F1
               M = M + 3
 1330       CONTINUE
 1300    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 3-center sssp and
C                ssps (AA|CD) type integrals.
C
C
   14    IF (SHELL3.EQ.1) THEN
             QXSUB = X3
             QYSUB = Y3
             QZSUB = Z3
         ELSE
             QXSUB = X4
             QYSUB = Y4
             QZSUB = Z4
         END IF

         M = 0
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

                   F0 = SCALE * F0
                   F1 = SCALE * F1
               ELSE
                   TINV = ONE / T
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = HALF * TINV * F0
               END IF

               F1 = F1 * PVAL * PQPINV

               BATCH (M+1) = (QXVAL-QXSUB) * F0 + PQX * F1
               BATCH (M+2) = (QYVAL-QYSUB) * F0 + PQY * F1
               BATCH (M+3) = (QZVAL-QZSUB) * F0 + PQZ * F1
               M = M + 3
 1440       CONTINUE
 1400    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 4-center sssp and
C                ssps (AB|CD) type integrals.
C
C
   15    IF (SHELL3.EQ.1) THEN
             QXSUB = X3
             QYSUB = Y3
             QZSUB = Z3
         ELSE
             QXSUB = X4
             QYSUB = Y4
             QZSUB = Z4
         END IF

         M = 0
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

                   F0 = SCALE * F0
                   F1 = SCALE * F1
               ELSE
                   TINV = ONE / T
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = HALF * TINV * F0
               END IF

               F1 = F1 * PVAL * PQPINV

               BATCH (M+1) = (QXVAL-QXSUB) * F0 + PQX * F1
               BATCH (M+2) = (QYVAL-QYSUB) * F0 + PQY * F1
               BATCH (M+3) = (QZVAL-QZSUB) * F0 + PQZ * F1
               M = M + 3
 1550       CONTINUE
 1500    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for atomic spss and
C                psss integrals. All these integrals are zero by
C                symmetry.
C
C
   21    RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 2-center spss and
C                psss (AA|CC) type integrals.
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
               ELSE
                   TINV = ONE / T
                   F1 = TINV * FOURTH * DSQRT (PI*TINV)
               END IF

               F1 = - F1 * QVAL * SCALE * PQPINV

               BATCH (M+1) = PQX * F1
               BATCH (M+2) = PQY * F1
               BATCH (M+3) = PQZ * F1
               M = M + 3
 2220       CONTINUE
 2200    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 3-center spss and
C                psss (AB|CC) type integrals.
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
            XP = PXVAL - PXSUB
            YP = PYVAL - PYSUB
            ZP = PZVAL - PZSUB
            PSCALE = SCALEP (IJ)
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

                   F0 = SCALE * F0
                   F1 = SCALE * F1
               ELSE
                   TINV = ONE / T
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = HALF * TINV * F0
               END IF

               F1 = - F1 * QVAL * PQPINV

               BATCH (M+1) = XP * F0 + PQX * F1
               BATCH (M+2) = YP * F0 + PQY * F1
               BATCH (M+3) = ZP * F0 + PQZ * F1
               M = M + 3
 2330       CONTINUE
 2300    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 3-center spss and
C                psss (AA|CD) type integrals.
C
C
   24    M = 0
         DO 2400 IJ = 1,MIJ
            PVAL = P (IJ)
            PSCALE = SCALEP (IJ)
            DO 2440 KL = 1,MKL
               QVAL = Q (KL)
               PQMULT = PVAL * QVAL
               PQPLUS = PVAL + QVAL
               PQPINV = ONE / PQPLUS
               PQX = X1 - QX (KL)
               PQY = Y1 - QY (KL)
               PQZ = Z1 - QZ (KL)
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
               ELSE
                   TINV = ONE / T
                   F1 = TINV * FOURTH * DSQRT (PI*TINV)
               END IF

               F1 = - F1 * QVAL * SCALE * PQPINV

               BATCH (M+1) = PQX * F1
               BATCH (M+2) = PQY * F1
               BATCH (M+3) = PQZ * F1
               M = M + 3
 2440       CONTINUE
 2400    CONTINUE

         RETURN
C
C
C             ...the 'K4' loop over all ij- and kl-exponent pairs
C                in present ij and kl block for 4-center spss and
C                psss (AB|CD) type integrals.
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

         M = 0
         DO 2500 IJ = 1,MIJ
            PVAL = P (IJ)
            PXVAL = PX (IJ)
            PYVAL = PY (IJ)
            PZVAL = PZ (IJ)
            XP = PXVAL - PXSUB
            YP = PYVAL - PYSUB
            ZP = PZVAL - PZSUB
            PSCALE = SCALEP (IJ)
            DO 2550 KL = 1,MKL
               QVAL = Q (KL)
               PQMULT = PVAL * QVAL
               PQPLUS = PVAL + QVAL
               PQPINV = ONE / PQPLUS
               PQX = PXVAL - QX (KL)
               PQY = PYVAL - QY (KL)
               PQZ = PZVAL - QZ (KL)
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

                   F0 = SCALE * F0
                   F1 = SCALE * F1
               ELSE
                   TINV = ONE / T
                   F0 = SCALE * HALF * DSQRT (PI*TINV)
                   F1 = HALF * TINV * F0
               END IF

               F1 = - F1 * QVAL * PQPINV

               BATCH (M+1) = XP * F0 + PQX * F1
               BATCH (M+2) = YP * F0 + PQY * F1
               BATCH (M+3) = ZP * F0 + PQZ * F1
               M = M + 3
 2550       CONTINUE
 2500    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
