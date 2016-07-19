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
         SUBROUTINE  ERD__SSSS_PCGTO_BLOCK
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
C  OPERATION   : ERD__SSSS_PCGTO_BLOCK
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation is designed to provide ultrafast block
C                evaluation of a batch of normalized electron repulsion
C                integrals between s-shell primitive gaussian type
C                orbitals.
C
C                A batch is defined here as containing all possible
C                integrals, that is its dimension is determined by
C                the total number of primitive functions (here = 1)
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
C                                     ssss
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
C                                    ssss integral batch
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
C                                    cartesian ssss integrals
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

         DOUBLE PRECISION  DELTA
         DOUBLE PRECISION  EXP1,EXP2,EXP3,EXP4
         DOUBLE PRECISION  F0
         DOUBLE PRECISION  PVAL,QVAL
         DOUBLE PRECISION  PQPLUS,PQMULT
         DOUBLE PRECISION  PQX,PQY,PQZ
         DOUBLE PRECISION  PSCALE
         DOUBLE PRECISION  PXVAL,PYVAL,PZVAL
         DOUBLE PRECISION  RNPQSQ
         DOUBLE PRECISION  SCALE
         DOUBLE PRECISION  T,TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4
         DOUBLE PRECISION  X12,Y12,Z12,X34,Y34,Z34
         DOUBLE PRECISION  SIXTH,FIFTH,FOURTH,THIRD,HALF,ONE,PI

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

         PARAMETER  (SIXTH   =  0.166666666666667D0)
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
               M = M + 1
               BATCH (M) = PSCALE * SCALEQ (KL)
     +                            / (PVAL * QVAL * DSQRT (PVAL+QVAL))
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
         RNPQSQ = PQX * PQX + PQY * PQY + PQZ * PQZ

         M = 0
         DO 2000 IJ = 1,MIJ
            PVAL = P (IJ)
            PSCALE = SCALEP (IJ)
            DO 2200 KL = 1,MKL
               QVAL = Q (KL)
               PQMULT = PVAL * QVAL
               PQPLUS = PVAL + QVAL
               T = RNPQSQ * PQMULT / PQPLUS
               SCALE = PSCALE * SCALEQ (KL) / (PQMULT * DSQRT (PQPLUS))

               IF (T.LE.TMAX) THEN
                   TGRID = INT (T * TVSTEP + HALF)
                   DELTA = TGRID * TSTEP - T

                   F0 = ((((( FTABLE (6,TGRID)  * DELTA * SIXTH
     +                      + FTABLE (5,TGRID)) * DELTA * FIFTH
     +                      + FTABLE (4,TGRID)) * DELTA * FOURTH
     +                      + FTABLE (3,TGRID)) * DELTA * THIRD
     +                      + FTABLE (2,TGRID)) * DELTA * HALF
     +                      + FTABLE (1,TGRID)) * DELTA
     +                      + FTABLE (0,TGRID)
               ELSE
                   F0 = HALF * DSQRT (PI/T)
               END IF

               M = M + 1
               BATCH (M) = SCALE * F0
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
            PQX = PX (IJ) - X3
            PQY = PY (IJ) - Y3
            PQZ = PZ (IJ) - Z3
            RNPQSQ = PQX * PQX + PQY * PQY + PQZ * PQZ
            PSCALE = SCALEP (IJ)
            DO 3300 KL = 1,MKL
               QVAL = Q (KL)
               PQMULT = PVAL * QVAL
               PQPLUS = PVAL + QVAL
               T = RNPQSQ * PQMULT / PQPLUS
               SCALE = PSCALE * SCALEQ (KL) / (PQMULT * DSQRT (PQPLUS))

               IF (T.LE.TMAX) THEN
                   TGRID = INT (T * TVSTEP + HALF)
                   DELTA = TGRID * TSTEP - T

                   F0 = ((((( FTABLE (6,TGRID)  * DELTA * SIXTH
     +                      + FTABLE (5,TGRID)) * DELTA * FIFTH
     +                      + FTABLE (4,TGRID)) * DELTA * FOURTH
     +                      + FTABLE (3,TGRID)) * DELTA * THIRD
     +                      + FTABLE (2,TGRID)) * DELTA * HALF
     +                      + FTABLE (1,TGRID)) * DELTA
     +                      + FTABLE (0,TGRID)
               ELSE
                   F0 = HALF * DSQRT (PI/T)
               END IF

               M = M + 1
               BATCH (M) = SCALE * F0
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
            DO 4400 KL = 1,MKL
               QVAL = Q (KL)
               PQMULT = PVAL * QVAL
               PQPLUS = PVAL + QVAL
               PQX = X1 - QX (KL)
               PQY = Y1 - QY (KL)
               PQZ = Z1 - QZ (KL)
               T = (PQX*PQX + PQY*PQY + PQZ*PQZ) * PQMULT / PQPLUS
               SCALE = PSCALE * SCALEQ (KL) / (PQMULT * DSQRT (PQPLUS))

               IF (T.LE.TMAX) THEN
                   TGRID = INT (T * TVSTEP + HALF)
                   DELTA = TGRID * TSTEP - T

                   F0 = ((((( FTABLE (6,TGRID)  * DELTA * SIXTH
     +                      + FTABLE (5,TGRID)) * DELTA * FIFTH
     +                      + FTABLE (4,TGRID)) * DELTA * FOURTH
     +                      + FTABLE (3,TGRID)) * DELTA * THIRD
     +                      + FTABLE (2,TGRID)) * DELTA * HALF
     +                      + FTABLE (1,TGRID)) * DELTA
     +                      + FTABLE (0,TGRID)
               ELSE
                   F0 = HALF * DSQRT (PI/T)
               END IF

               M = M + 1
               BATCH (M) = SCALE * F0
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
            DO 5500 KL = 1,MKL
               QVAL = Q (KL)
               PQMULT = PVAL * QVAL
               PQPLUS = PVAL + QVAL
               PQX = PXVAL - QX (KL)
               PQY = PYVAL - QY (KL)
               PQZ = PZVAL - QZ (KL)
               T = (PQX*PQX + PQY*PQY + PQZ*PQZ) * PQMULT / PQPLUS
               SCALE = PSCALE * SCALEQ (KL) / (PQMULT * DSQRT (PQPLUS))

               IF (T.LE.TMAX) THEN
                   TGRID = INT (T * TVSTEP + HALF)
                   DELTA = TGRID * TSTEP - T

                   F0 = ((((( FTABLE (6,TGRID)  * DELTA * SIXTH
     +                      + FTABLE (5,TGRID)) * DELTA * FIFTH
     +                      + FTABLE (4,TGRID)) * DELTA * FOURTH
     +                      + FTABLE (3,TGRID)) * DELTA * THIRD
     +                      + FTABLE (2,TGRID)) * DELTA * HALF
     +                      + FTABLE (1,TGRID)) * DELTA
     +                      + FTABLE (0,TGRID)
               ELSE
                   F0 = HALF * DSQRT (PI/T)
               END IF

               M = M + 1
               BATCH (M) = SCALE * F0
 5500       CONTINUE
 5000    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
