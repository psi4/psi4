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
         SUBROUTINE  ERD__1111_CSGTO
     +
     +                    ( IMAX,ZMAX,
     +                      NALPHA,NCOEFF,NCSUM,
     +                      NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                      NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                      SHELL1,SHELL2,SHELL3,SHELL4,
     +                      X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                      ALPHA,CC,CCBEG,CCEND,
     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                      L1CACHE,TILE,NCTROW,
     +                      SCREEN,
     +                      ICORE,
     +
     +                                NBATCH,
     +                                NFIRST,
     +                                ZCORE )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__1111_CSGTO
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : ERD__SET_IJ_KL_PAIRS
C                ERD__1111_DEF_BLOCKS
C                ERD__PREPARE_CTR
C                ERD__SSSS_PCGTO_BLOCK
C                ERD__SSSP_PCGTO_BLOCK
C                ERD__SSPP_PCGTO_BLOCK
C                ERD__SPPP_PCGTO_BLOCK
C                ERD__PPPP_PCGTO_BLOCK
C                ERD__CTR_4INDEX_BLOCK
C                ERD__CTR_RS_EXPAND
C                ERD__CTR_TU_EXPAND
C                ERD__CTR_4INDEX_REORDER
C                ERD__MAP_IJKL_TO_IKJL
C  DESCRIPTION : This operation calculates a batch of contracted
C                electron repulsion integrals on up to four different
C                centers between spherical gaussian type shells.
C
C                Special fast routine for integrals involving s- and
C                p-type shells only!
C
C
C                  Input (x = 1,2,3 and 4):
C
C                    IMAX,ZMAX    =  maximum integer,flp memory
C                    NALPHA       =  total # of exponents
C                    NCOEFF       =  total # of contraction coeffs
C                    NCSUM        =  total # of contractions
C                    NCGTOx       =  # of contractions for csh x
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x
C                    SHELLx       =  the shell type for csh x
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = 1,2,3 and 4
C                    ALPHA        =  primitive exponents for csh
C                                    1,2,3,4 in that order
C                    CC           =  full set (including zeros) of
C                                    contraction coefficients for csh
C                                    1,2,3,4 in that order, for each
C                                    csh individually such that an
C                                    (I,J) element corresponds to the
C                                    I-th primitive and J-th contraction.
C                    CC(BEG)END   =  (lowest)highest nonzero primitive
C                                    index for contractions for csh
C                                    1,2,3,4 in that order. They are
C                                    different from (1)NPGTOx only for
C                                    segmented contractions
C                    FTABLE       =  Fm (T) table for interpolation
C                                    in low T region
C                    MGRID        =  maximum m in Fm (T) table
C                    NGRID        =  # of T's for which Fm (T) table
C                                    was set up
C                    TMAX         =  maximum T in Fm (T) table
C                    TSTEP        =  difference between two consecutive
C                                    T's in Fm (T) table
C                    TVSTEP       =  Inverse of TSTEP
C                    L1CACHE      =  Size of level 1 cache in units of
C                                    8 Byte
C                    TILE         =  Number of rows and columns in
C                                    units of 8 Byte of level 1 cache
C                                    square tile array used for
C                                    performing optimum matrix
C                                    transpositions
C                    NCTROW       =  minimum # of rows that are
C                                    accepted for blocked contractions
C                    SCREEN       =  is true, if screening will be
C                                    done at primitive integral level
C                    ICORE        =  integer scratch space
C                    ZCORE (part) =  flp scratch space
C
C                  Output:
C
C                    NBATCH       =  # of integrals in batch
C                    NFIRST       =  first address location inside the
C                                    ZCORE array containing the first
C                                    integral
C                    ZCORE        =  full batch of contracted (12|34)
C                                    integrals over spherical gaussians
C                                    starting at ZCORE (NFIRST)
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

         LOGICAL     ATOMIC,ATOM12,ATOM34,ATOM23
         LOGICAL     BLOCKED
         LOGICAL     EMPTY
         LOGICAL     EQUAL12,EQUAL34
         LOGICAL     MEMORY
         LOGICAL     REORDER
         LOGICAL     SCREEN
         LOGICAL     SWAPRS,SWAPTU

         INTEGER     CTMOVE,XTMOVE
         INTEGER     I,J,K,L
         INTEGER     IMAX,ZMAX
         INTEGER     IN,OUT
         INTEGER     INDEXR,INDEXS,INDEXT,INDEXU
         INTEGER     IPRIM1,IPRIM2,IPRIM3,IPRIM4
         INTEGER     IPUSED,IPSAVE,IPPAIR
         INTEGER     L1CACHE,TILE,NCTROW
         INTEGER     LCC1,LCC2,LCC3,LCC4
         INTEGER     LCCSEG1,LCCSEG2,LCCSEG3,LCCSEG4
         INTEGER     LEXP1,LEXP2,LEXP3,LEXP4
         INTEGER     MIJ,MKL,MIJKL
         INTEGER     MGRID,NGRID
         INTEGER     MXPRIM,MNPRIM
         INTEGER     NALPHA,NCOEFF,NCSUM
         INTEGER     NBATCH,NFIRST
         INTEGER     NCGTO1,NCGTO2,NCGTO3,NCGTO4,NCGTO12,NCGTO34
         INTEGER     NCGTOR,NCGTOS,NCGTOT,NCGTOU
         INTEGER     NCTR
         INTEGER     NIJ,NKL
         INTEGER     NIJBLK,NKLBLK,NIJBEG,NKLBEG,NIJEND,NKLEND
         INTEGER     NPGTO1,NPGTO2,NPGTO3,NPGTO4,NPGTO12,NPGTO34
         INTEGER     NPSIZE,NCSIZE,NWSIZE
         INTEGER     NXYZ1,NXYZ2,NXYZ3,NXYZ4,NXYZT
         INTEGER     SHELL1,SHELL2,SHELL3,SHELL4,SHELLP,SHELLT
         INTEGER     ZCBATCH,ZPBATCH,ZWORK,
     +               ZNORM1,ZNORM2,ZNORM3,ZNORM4,
     +               ZRHO12,ZRHO34,
     +               ZP,ZPX,ZPY,ZPZ,ZSCPK2,
     +               ZQ,ZQX,ZQY,ZQZ,ZSCQK2

         INTEGER     CCBEG (1:NCSUM)
         INTEGER     CCEND (1:NCSUM)
         INTEGER     ICORE (1:IMAX)
         INTEGER     IXOFF (1:4)

         DOUBLE PRECISION    PREFACT
         DOUBLE PRECISION    RN12SQ,RN34SQ
         DOUBLE PRECISION    SPNORM
         DOUBLE PRECISION    TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION    X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4
         DOUBLE PRECISION    X12,Y12,Z12,X34,Y34,Z34
         DOUBLE PRECISION    ZERO,ONE

         DOUBLE PRECISION    ALPHA (1:NALPHA)
         DOUBLE PRECISION    CC    (1:NCOEFF)
         DOUBLE PRECISION    ZCORE (1:ZMAX)

         DOUBLE PRECISION    FTABLE (0:MGRID,0:NGRID)

         PARAMETER  (ZERO    = 0.D0)
         PARAMETER  (ONE     = 1.D0)
         PARAMETER  (PREFACT = 9.027033336764101D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...decide as early as possible if a zero batch of
C                integrals is expected.
C
C
         SHELLP = SHELL1 + SHELL2
         SHELLT = SHELLP + SHELL3 + SHELL4

         ATOM12 = (X1.EQ.X2) .AND. (Y1.EQ.Y2) .AND. (Z1.EQ.Z2)
         ATOM23 = (X2.EQ.X3) .AND. (Y2.EQ.Y3) .AND. (Z2.EQ.Z3)
         ATOM34 = (X3.EQ.X4) .AND. (Y3.EQ.Y4) .AND. (Z3.EQ.Z4)

         ATOMIC = (ATOM12 .AND. ATOM34 .AND. ATOM23)

         IF (ATOMIC .AND. (MOD(SHELLT,2).EQ.1)) THEN
             NBATCH = 0
             RETURN
         END IF
C
C
C             ...set the pointers to the alpha exponents, contraction
C                coefficients and segmented contraction boundaries.
C
C
         LEXP1 = 1
         LEXP2 = LEXP1 + NPGTO1
         LEXP3 = LEXP2 + NPGTO2
         LEXP4 = LEXP3 + NPGTO3

         LCC1 = 1
         LCC2 = LCC1 + NPGTO1 * NCGTO1
         LCC3 = LCC2 + NPGTO2 * NCGTO2
         LCC4 = LCC3 + NPGTO3 * NCGTO3

         LCCSEG1 = 1
         LCCSEG2 = LCCSEG1 + NCGTO1
         LCCSEG3 = LCCSEG2 + NCGTO2
         LCCSEG4 = LCCSEG3 + NCGTO3
C
C
C             ...determine csh equality between center pairs 1,2
C                and 3,4 in increasing order of complexity:
C
C                 centers -> shells -> exponents -> ctr coefficients
C
C
         EQUAL12 = ATOM12

         IF (EQUAL12) THEN
             EQUAL12 =     (SHELL1 .EQ. SHELL2)
     +               .AND. (NPGTO1 .EQ. NPGTO2)
     +               .AND. (NCGTO1 .EQ. NCGTO2)
             IF (EQUAL12) THEN
               K = LEXP1 - 1
               L = LEXP2 - 1
               DO 100 I = 1,NPGTO1
                  EQUAL12 = EQUAL12 .AND. (ALPHA (K+I).EQ.ALPHA (L+I))
  100          CONTINUE
               IF (EQUAL12) THEN
                 K = LCC1 - 1
                 L = LCC2 - 1
                 DO 110 J = 1,NCGTO1
                    IF (EQUAL12) THEN
                      DO 120 I = 1,NPGTO1
                         EQUAL12 = EQUAL12 .AND. (CC (K+I).EQ.CC (L+I))
  120                 CONTINUE
                      K = K + NPGTO1
                      L = L + NPGTO1
                    END IF
  110            CONTINUE
               END IF
             END IF
         END IF

         EQUAL34 = ATOM34

         IF (EQUAL34) THEN
             EQUAL34 =     (SHELL3 .EQ. SHELL4)
     +               .AND. (NPGTO3 .EQ. NPGTO4)
     +               .AND. (NCGTO3 .EQ. NCGTO4)
             IF (EQUAL34) THEN
               K = LEXP3 - 1
               L = LEXP4 - 1
               DO 130 I = 1,NPGTO3
                  EQUAL34 = EQUAL34 .AND. (ALPHA (K+I).EQ.ALPHA (L+I))
  130          CONTINUE
               IF (EQUAL34) THEN
                 K = LCC3 - 1
                 L = LCC4 - 1
                 DO 140 J = 1,NCGTO3
                    IF (EQUAL34) THEN
                      DO 150 I = 1,NPGTO3
                         EQUAL34 = EQUAL34 .AND. (CC (K+I).EQ.CC (L+I))
  150                 CONTINUE
                      K = K + NPGTO3
                      L = L + NPGTO3
                    END IF
  140            CONTINUE
               END IF
             END IF
         END IF
C
C
C             ...calculate relevant data for the [12|34] batch of
C                integrals, such as dimensions, total # of integrals
C                to be expected, relevant ij and kl primitive exponent
C                pairs, etc... The integral prefactor PREFACT has been
C                set as a parameter, its value being = 16 / sqrt(pi).
C                Calculate here also the overall norm factor SPNORM due
C                to presence of s- or p-type shells. The contribution
C                to SPNORM is very simple: each s-type shell -> * 1.0,
C                each p-type shell -> * 2.0.
C
C
         NXYZ1  = SHELL1 + SHELL1 + 1
         NXYZ2  = SHELL2 + SHELL2 + 1
         NXYZ3  = SHELL3 + SHELL3 + 1
         NXYZ4  = SHELL4 + SHELL4 + 1

         NXYZT  = NXYZ1 * NXYZ2 * NXYZ3 * NXYZ4

         IF (.NOT.ATOM12) THEN
             X12 = X1 - X2
             Y12 = Y1 - Y2
             Z12 = Z1 - Z2
             RN12SQ = X12 * X12 + Y12 * Y12 + Z12 * Z12
         ELSE
             X12 = ZERO
             Y12 = ZERO
             Z12 = ZERO
             RN12SQ = ZERO
         END IF

         IF (.NOT.ATOM34) THEN
             X34 = X3 - X4
             Y34 = Y3 - Y4
             Z34 = Z3 - Z4
             RN34SQ = X34 * X34 + Y34 * Y34 + Z34 * Z34
         ELSE
             X34 = ZERO
             Y34 = ZERO
             Z34 = ZERO
             RN34SQ = ZERO
         END IF

         IF (EQUAL12) THEN
             NPGTO12 = (NPGTO1*(NPGTO1+1))/2
             NCGTO12 = (NCGTO1*(NCGTO1+1))/2
         ELSE
             NPGTO12 = NPGTO1 * NPGTO2
             NCGTO12 = NCGTO1 * NCGTO2
         END IF

         IF (EQUAL34) THEN
             NPGTO34 = (NPGTO3*(NPGTO3+1))/2
             NCGTO34 = (NCGTO3*(NCGTO3+1))/2
         ELSE
             NPGTO34 = NPGTO3 * NPGTO4
             NCGTO34 = NCGTO3 * NCGTO4
         END IF

         NCTR = NCGTO12 * NCGTO34

         IPRIM1 = 1
         IPRIM2 = IPRIM1 + NPGTO12
         IPRIM3 = IPRIM2 + NPGTO12
         IPRIM4 = IPRIM3 + NPGTO34

         SWAPRS  = NPGTO1 .GT. NPGTO2
         SWAPTU  = NPGTO3 .GT. NPGTO4

         SPNORM = ONE
         IF (SHELL1.EQ.1) THEN
             SPNORM = SPNORM + SPNORM
         END IF
         IF (SHELL2.EQ.1) THEN
             SPNORM = SPNORM + SPNORM
         END IF
         IF (SHELL3.EQ.1) THEN
             SPNORM = SPNORM + SPNORM
         END IF
         IF (SHELL4.EQ.1) THEN
             SPNORM = SPNORM + SPNORM
         END IF

         CALL  ERD__SET_IJ_KL_PAIRS
     +
     +              ( NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                NPGTO12,NPGTO34,
     +                ATOM12,ATOM34,
     +                EQUAL12,EQUAL34,
     +                SWAPRS,SWAPTU,
     +                X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                RN12SQ,RN34SQ,
     +                PREFACT,
     +                ALPHA (LEXP1),ALPHA (LEXP2),
     +                ALPHA (LEXP3),ALPHA (LEXP4),
     +                FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                SCREEN,
     +
     +                         EMPTY,
     +                         NIJ,NKL,
     +                         ICORE (IPRIM1),ICORE (IPRIM2),
     +                         ICORE (IPRIM3),ICORE (IPRIM4),
     +                         ZCORE (1) )
     +
     +
         IF (EMPTY) THEN
             NBATCH = 0
             RETURN
         END IF
C
C
C             ...decide on the primitive [12|34] block size and
C                return array sizes and pointers for the primitive
C                [12|34] generation. Perform also some preparation
C                steps for contraction.
C
C
         MEMORY = .FALSE.

         CALL  ERD__1111_DEF_BLOCKS
     +
     +              ( ZMAX,
     +                NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                NIJ,NKL,
     +                NCGTO12,NCGTO34,NCTR,
     +                NXYZT,
     +                L1CACHE,NCTROW,
     +                MEMORY,
     +
     +                          NIJBLK,NKLBLK,
     +                          NPSIZE,NCSIZE,NWSIZE,
     +                          MXPRIM,MNPRIM,
     +                          ZCBATCH,ZPBATCH,ZWORK,
     +                          ZNORM1,ZNORM2,ZNORM3,ZNORM4,
     +                          ZRHO12,ZRHO34,
     +                          ZP,ZPX,ZPY,ZPZ,ZSCPK2,
     +                          ZQ,ZQX,ZQY,ZQZ,ZSCQK2 )
     +
     +
         BLOCKED = (NIJBLK.LT.NIJ) .OR. (NKLBLK.LT.NKL)

         CALL  ERD__PREPARE_CTR
     +
     +              ( NCSIZE,
     +                NIJ,NKL,
     +                NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                SHELL1,SHELL2,SHELL3,SHELL4,
     +                ALPHA (LEXP1),ALPHA (LEXP2),
     +                ALPHA (LEXP3),ALPHA (LEXP4),
     +                PREFACT,SPNORM,
     +                EQUAL12,EQUAL34,
     +                BLOCKED,
     +                ZCORE (1),
     +
     +                          ZCORE (ZNORM1),ZCORE (ZNORM2),
     +                          ZCORE (ZNORM3),ZCORE (ZNORM4),
     +                          ZCORE (ZRHO12),ZCORE (ZRHO34),
     +                          ZCORE (ZCBATCH) )
     +
     +
         IPUSED = IPRIM4 + NPGTO34
         IPSAVE = IPUSED + MNPRIM
         IPPAIR = IPSAVE + MXPRIM
C
C
C             ...evaluate [12|34] in blocks over ij and kl pairs
C                and add to final contracted (12|34) with full
C                contracted index ranges r,s,t,u. The keyword REORDER
C                indicates, if the primitive [12|34] blocks needs to
C                be transposed before being contracted.
C
C
         REORDER = .FALSE.

         IF (SHELLT.EQ.0) THEN

             DO 1000 NIJBEG = 1,NIJ,NIJBLK
                NIJEND = MIN0 (NIJBEG+NIJBLK-1,NIJ)
                MIJ = NIJEND - NIJBEG + 1
                DO 1100 NKLBEG = 1,NKL,NKLBLK
                   NKLEND = MIN0 (NKLBEG+NKLBLK-1,NKL)
                   MKL = NKLEND - NKLBEG + 1
                   MIJKL = MIJ * MKL
                   NPSIZE = NXYZT * MIJKL

                   CALL  ERD__SSSS_PCGTO_BLOCK
     +
     +                        ( NPSIZE,
     +                          ATOMIC,ATOM12,ATOM34,
     +                          MIJ,MKL,
     +                          NIJ,NIJBEG,NIJEND,NKL,NKLBEG,NKLEND,
     +                          NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                          X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                          X12,Y12,Z12,X34,Y34,Z34,
     +                          ALPHA (LEXP1),ALPHA (LEXP2),
     +                          ALPHA (LEXP3),ALPHA (LEXP4),
     +                          FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                          ICORE (IPRIM1+NIJBEG-1),
     +                          ICORE (IPRIM2+NIJBEG-1),
     +                          ICORE (IPRIM3+NKLBEG-1),
     +                          ICORE (IPRIM4+NKLBEG-1),
     +                          ZCORE (ZNORM1),ZCORE (ZNORM2),
     +                          ZCORE (ZNORM3),ZCORE (ZNORM4),
     +                          ZCORE (ZRHO12),ZCORE (ZRHO34),
     +                          ZCORE (ZP),
     +                          ZCORE (ZPX),ZCORE (ZPY),ZCORE (ZPZ),
     +                          ZCORE (ZSCPK2),
     +                          ZCORE (ZQ),
     +                          ZCORE (ZQX),ZCORE (ZQY),ZCORE (ZQZ),
     +                          ZCORE (ZSCQK2),
     +
     +                               ZCORE (ZPBATCH) )
     +
     +
                   CALL  ERD__CTR_4INDEX_BLOCK
     +
     +                        ( NPSIZE,NCSIZE,NWSIZE,
     +                          NXYZT,MIJKL,
     +                          MIJ,MKL,NCGTO12,NCGTO34,
     +                          NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                          NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                          MXPRIM,MNPRIM,
     +                          CC (LCC1),CC (LCC2),CC (LCC3),CC (LCC4),
     +                          CCBEG (LCCSEG1),CCBEG (LCCSEG2),
     +                          CCBEG (LCCSEG3),CCBEG (LCCSEG4),
     +                          CCEND (LCCSEG1),CCEND (LCCSEG2),
     +                          CCEND (LCCSEG3),CCEND (LCCSEG4),
     +                          ICORE (IPRIM1+NIJBEG-1),
     +                          ICORE (IPRIM2+NIJBEG-1),
     +                          ICORE (IPRIM3+NKLBEG-1),
     +                          ICORE (IPRIM4+NKLBEG-1),
     +                          L1CACHE,TILE,NCTROW,
     +                          EQUAL12,EQUAL34,
     +                          SWAPRS,SWAPTU,
     +                          REORDER,
     +                          BLOCKED,
     +                          ICORE (IPUSED),
     +                          ICORE (IPSAVE),
     +                          ICORE (IPPAIR),
     +                          ZCORE (ZPBATCH),
     +                          ZCORE (ZWORK),
     +
     +                                    ZCORE (ZCBATCH) )
     +
     +
 1100           CONTINUE
 1000        CONTINUE

         ELSE IF (SHELLT.EQ.1) THEN

             DO 2000 NIJBEG = 1,NIJ,NIJBLK
                NIJEND = MIN0 (NIJBEG+NIJBLK-1,NIJ)
                MIJ = NIJEND - NIJBEG + 1
                DO 2200 NKLBEG = 1,NKL,NKLBLK
                   NKLEND = MIN0 (NKLBEG+NKLBLK-1,NKL)
                   MKL = NKLEND - NKLBEG + 1
                   MIJKL = MIJ * MKL
                   NPSIZE = NXYZT * MIJKL

                   CALL  ERD__SSSP_PCGTO_BLOCK
     +
     +                        ( NPSIZE,
     +                          ATOMIC,ATOM12,ATOM34,
     +                          MIJ,MKL,
     +                          NIJ,NIJBEG,NIJEND,NKL,NKLBEG,NKLEND,
     +                          NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                          SHELL1,SHELL3,SHELLP,
     +                          X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                          X12,Y12,Z12,X34,Y34,Z34,
     +                          ALPHA (LEXP1),ALPHA (LEXP2),
     +                          ALPHA (LEXP3),ALPHA (LEXP4),
     +                          FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                          ICORE (IPRIM1+NIJBEG-1),
     +                          ICORE (IPRIM2+NIJBEG-1),
     +                          ICORE (IPRIM3+NKLBEG-1),
     +                          ICORE (IPRIM4+NKLBEG-1),
     +                          ZCORE (ZNORM1),ZCORE (ZNORM2),
     +                          ZCORE (ZNORM3),ZCORE (ZNORM4),
     +                          ZCORE (ZRHO12),ZCORE (ZRHO34),
     +                          ZCORE (ZP),
     +                          ZCORE (ZPX),ZCORE (ZPY),ZCORE (ZPZ),
     +                          ZCORE (ZSCPK2),
     +                          ZCORE (ZQ),
     +                          ZCORE (ZQX),ZCORE (ZQY),ZCORE (ZQZ),
     +                          ZCORE (ZSCQK2),
     +
     +                               ZCORE (ZPBATCH) )
     +
     +
                   CALL  ERD__CTR_4INDEX_BLOCK
     +
     +                        ( NPSIZE,NCSIZE,NWSIZE,
     +                          NXYZT,MIJKL,
     +                          MIJ,MKL,NCGTO12,NCGTO34,
     +                          NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                          NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                          MXPRIM,MNPRIM,
     +                          CC (LCC1),CC (LCC2),CC (LCC3),CC (LCC4),
     +                          CCBEG (LCCSEG1),CCBEG (LCCSEG2),
     +                          CCBEG (LCCSEG3),CCBEG (LCCSEG4),
     +                          CCEND (LCCSEG1),CCEND (LCCSEG2),
     +                          CCEND (LCCSEG3),CCEND (LCCSEG4),
     +                          ICORE (IPRIM1+NIJBEG-1),
     +                          ICORE (IPRIM2+NIJBEG-1),
     +                          ICORE (IPRIM3+NKLBEG-1),
     +                          ICORE (IPRIM4+NKLBEG-1),
     +                          L1CACHE,TILE,NCTROW,
     +                          EQUAL12,EQUAL34,
     +                          SWAPRS,SWAPTU,
     +                          REORDER,
     +                          BLOCKED,
     +                          ICORE (IPUSED),
     +                          ICORE (IPSAVE),
     +                          ICORE (IPPAIR),
     +                          ZCORE (ZPBATCH),
     +                          ZCORE (ZWORK),
     +
     +                                    ZCORE (ZCBATCH) )
     +
     +
 2200           CONTINUE
 2000        CONTINUE

         ELSE IF (SHELLT.EQ.2) THEN

             DO 3000 NIJBEG = 1,NIJ,NIJBLK
                NIJEND = MIN0 (NIJBEG+NIJBLK-1,NIJ)
                MIJ = NIJEND - NIJBEG + 1
                DO 3300 NKLBEG = 1,NKL,NKLBLK
                   NKLEND = MIN0 (NKLBEG+NKLBLK-1,NKL)
                   MKL = NKLEND - NKLBEG + 1
                   MIJKL = MIJ * MKL
                   NPSIZE = NXYZT * MIJKL

                   CALL  ERD__SSPP_PCGTO_BLOCK
     +
     +                        ( NPSIZE,
     +                          ATOMIC,ATOM12,ATOM34,
     +                          MIJ,MKL,
     +                          NIJ,NIJBEG,NIJEND,NKL,NKLBEG,NKLEND,
     +                          NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                          SHELL1,SHELL3,SHELLP,
     +                          X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                          X12,Y12,Z12,X34,Y34,Z34,
     +                          ALPHA (LEXP1),ALPHA (LEXP2),
     +                          ALPHA (LEXP3),ALPHA (LEXP4),
     +                          FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                          ICORE (IPRIM1+NIJBEG-1),
     +                          ICORE (IPRIM2+NIJBEG-1),
     +                          ICORE (IPRIM3+NKLBEG-1),
     +                          ICORE (IPRIM4+NKLBEG-1),
     +                          ZCORE (ZNORM1),ZCORE (ZNORM2),
     +                          ZCORE (ZNORM3),ZCORE (ZNORM4),
     +                          ZCORE (ZRHO12),ZCORE (ZRHO34),
     +                          ZCORE (ZP),
     +                          ZCORE (ZPX),ZCORE (ZPY),ZCORE (ZPZ),
     +                          ZCORE (ZSCPK2),
     +                          ZCORE (ZQ),
     +                          ZCORE (ZQX),ZCORE (ZQY),ZCORE (ZQZ),
     +                          ZCORE (ZSCQK2),
     +
     +                               ZCORE (ZPBATCH) )
     +
     +
                   CALL  ERD__CTR_4INDEX_BLOCK
     +
     +                        ( NPSIZE,NCSIZE,NWSIZE,
     +                          NXYZT,MIJKL,
     +                          MIJ,MKL,NCGTO12,NCGTO34,
     +                          NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                          NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                          MXPRIM,MNPRIM,
     +                          CC (LCC1),CC (LCC2),CC (LCC3),CC (LCC4),
     +                          CCBEG (LCCSEG1),CCBEG (LCCSEG2),
     +                          CCBEG (LCCSEG3),CCBEG (LCCSEG4),
     +                          CCEND (LCCSEG1),CCEND (LCCSEG2),
     +                          CCEND (LCCSEG3),CCEND (LCCSEG4),
     +                          ICORE (IPRIM1+NIJBEG-1),
     +                          ICORE (IPRIM2+NIJBEG-1),
     +                          ICORE (IPRIM3+NKLBEG-1),
     +                          ICORE (IPRIM4+NKLBEG-1),
     +                          L1CACHE,TILE,NCTROW,
     +                          EQUAL12,EQUAL34,
     +                          SWAPRS,SWAPTU,
     +                          REORDER,
     +                          BLOCKED,
     +                          ICORE (IPUSED),
     +                          ICORE (IPSAVE),
     +                          ICORE (IPPAIR),
     +                          ZCORE (ZPBATCH),
     +                          ZCORE (ZWORK),
     +
     +                                    ZCORE (ZCBATCH) )
     +
     +
 3300           CONTINUE
 3000        CONTINUE

         ELSE IF (SHELLT.EQ.3) THEN

             DO 4000 NIJBEG = 1,NIJ,NIJBLK
                NIJEND = MIN0 (NIJBEG+NIJBLK-1,NIJ)
                MIJ = NIJEND - NIJBEG + 1
                DO 4400 NKLBEG = 1,NKL,NKLBLK
                   NKLEND = MIN0 (NKLBEG+NKLBLK-1,NKL)
                   MKL = NKLEND - NKLBEG + 1
                   MIJKL = MIJ * MKL
                   NPSIZE = NXYZT * MIJKL

                   CALL  ERD__SPPP_PCGTO_BLOCK
     +
     +                        ( NPSIZE,
     +                          ATOMIC,ATOM12,ATOM34,
     +                          MIJ,MKL,
     +                          NIJ,NIJBEG,NIJEND,NKL,NKLBEG,NKLEND,
     +                          NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                          SHELL1,SHELL3,SHELLP,
     +                          X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                          X12,Y12,Z12,X34,Y34,Z34,
     +                          ALPHA (LEXP1),ALPHA (LEXP2),
     +                          ALPHA (LEXP3),ALPHA (LEXP4),
     +                          FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                          ICORE (IPRIM1+NIJBEG-1),
     +                          ICORE (IPRIM2+NIJBEG-1),
     +                          ICORE (IPRIM3+NKLBEG-1),
     +                          ICORE (IPRIM4+NKLBEG-1),
     +                          ZCORE (ZNORM1),ZCORE (ZNORM2),
     +                          ZCORE (ZNORM3),ZCORE (ZNORM4),
     +                          ZCORE (ZRHO12),ZCORE (ZRHO34),
     +                          ZCORE (ZP),
     +                          ZCORE (ZPX),ZCORE (ZPY),ZCORE (ZPZ),
     +                          ZCORE (ZSCPK2),
     +                          ZCORE (ZQ),
     +                          ZCORE (ZQX),ZCORE (ZQY),ZCORE (ZQZ),
     +                          ZCORE (ZSCQK2),
     +
     +                               ZCORE (ZPBATCH) )
     +
     +
                   CALL  ERD__CTR_4INDEX_BLOCK
     +
     +                        ( NPSIZE,NCSIZE,NWSIZE,
     +                          NXYZT,MIJKL,
     +                          MIJ,MKL,NCGTO12,NCGTO34,
     +                          NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                          NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                          MXPRIM,MNPRIM,
     +                          CC (LCC1),CC (LCC2),CC (LCC3),CC (LCC4),
     +                          CCBEG (LCCSEG1),CCBEG (LCCSEG2),
     +                          CCBEG (LCCSEG3),CCBEG (LCCSEG4),
     +                          CCEND (LCCSEG1),CCEND (LCCSEG2),
     +                          CCEND (LCCSEG3),CCEND (LCCSEG4),
     +                          ICORE (IPRIM1+NIJBEG-1),
     +                          ICORE (IPRIM2+NIJBEG-1),
     +                          ICORE (IPRIM3+NKLBEG-1),
     +                          ICORE (IPRIM4+NKLBEG-1),
     +                          L1CACHE,TILE,NCTROW,
     +                          EQUAL12,EQUAL34,
     +                          SWAPRS,SWAPTU,
     +                          REORDER,
     +                          BLOCKED,
     +                          ICORE (IPUSED),
     +                          ICORE (IPSAVE),
     +                          ICORE (IPPAIR),
     +                          ZCORE (ZPBATCH),
     +                          ZCORE (ZWORK),
     +
     +                                    ZCORE (ZCBATCH) )
     +
     +
 4400           CONTINUE
 4000        CONTINUE

         ELSE

             DO 5000 NIJBEG = 1,NIJ,NIJBLK
                NIJEND = MIN0 (NIJBEG+NIJBLK-1,NIJ)
                MIJ = NIJEND - NIJBEG + 1
                DO 5500 NKLBEG = 1,NKL,NKLBLK
                   NKLEND = MIN0 (NKLBEG+NKLBLK-1,NKL)
                   MKL = NKLEND - NKLBEG + 1
                   MIJKL = MIJ * MKL
                   NPSIZE = NXYZT * MIJKL

                   CALL  ERD__PPPP_PCGTO_BLOCK
     +
     +                        ( NPSIZE,
     +                          ATOMIC,ATOM12,ATOM34,
     +                          MIJ,MKL,
     +                          NIJ,NIJBEG,NIJEND,NKL,NKLBEG,NKLEND,
     +                          NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                          X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                          X12,Y12,Z12,X34,Y34,Z34,
     +                          ALPHA (LEXP1),ALPHA (LEXP2),
     +                          ALPHA (LEXP3),ALPHA (LEXP4),
     +                          FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                          ICORE (IPRIM1+NIJBEG-1),
     +                          ICORE (IPRIM2+NIJBEG-1),
     +                          ICORE (IPRIM3+NKLBEG-1),
     +                          ICORE (IPRIM4+NKLBEG-1),
     +                          ZCORE (ZNORM1),ZCORE (ZNORM2),
     +                          ZCORE (ZNORM3),ZCORE (ZNORM4),
     +                          ZCORE (ZRHO12),ZCORE (ZRHO34),
     +                          ZCORE (ZP),
     +                          ZCORE (ZPX),ZCORE (ZPY),ZCORE (ZPZ),
     +                          ZCORE (ZSCPK2),
     +                          ZCORE (ZQ),
     +                          ZCORE (ZQX),ZCORE (ZQY),ZCORE (ZQZ),
     +                          ZCORE (ZSCQK2),
     +
     +                               ZCORE (ZPBATCH) )
     +
     +
                   CALL  ERD__CTR_4INDEX_BLOCK
     +
     +                        ( NPSIZE,NCSIZE,NWSIZE,
     +                          NXYZT,MIJKL,
     +                          MIJ,MKL,NCGTO12,NCGTO34,
     +                          NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                          NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                          MXPRIM,MNPRIM,
     +                          CC (LCC1),CC (LCC2),CC (LCC3),CC (LCC4),
     +                          CCBEG (LCCSEG1),CCBEG (LCCSEG2),
     +                          CCBEG (LCCSEG3),CCBEG (LCCSEG4),
     +                          CCEND (LCCSEG1),CCEND (LCCSEG2),
     +                          CCEND (LCCSEG3),CCEND (LCCSEG4),
     +                          ICORE (IPRIM1+NIJBEG-1),
     +                          ICORE (IPRIM2+NIJBEG-1),
     +                          ICORE (IPRIM3+NKLBEG-1),
     +                          ICORE (IPRIM4+NKLBEG-1),
     +                          L1CACHE,TILE,NCTROW,
     +                          EQUAL12,EQUAL34,
     +                          SWAPRS,SWAPTU,
     +                          REORDER,
     +                          BLOCKED,
     +                          ICORE (IPUSED),
     +                          ICORE (IPSAVE),
     +                          ICORE (IPPAIR),
     +                          ZCORE (ZPBATCH),
     +                          ZCORE (ZWORK),
     +
     +                                    ZCORE (ZCBATCH) )
     +
     +
 5500           CONTINUE
 5000        CONTINUE

         END IF
C
C
C             ...expand the contraction indices (if necessary):
C
C                   batch (nxyzt,r>=s,t>=u) --> batch (nxyzt,r,s,t,u)
C
C                and reorder the contraction indices (if necessary):
C
C                   batch (nxyzt,r,s,t,u) --> batch (nxyzt,i,j,k,l)
C
C                such that they are in final correspondence:
C
C                                    i --> 1
C                                    j --> 2
C                                    k --> 3
C                                    l --> 4
C
         NCTR = NCGTO1 * NCGTO2 * NCGTO3 * NCGTO4
         NBATCH = NXYZT * NCTR

         IN = ZCBATCH
         OUT = ZCBATCH + NBATCH

         IF (EQUAL12 .AND. NCGTO12.GT.1) THEN
             CALL  ERD__CTR_RS_EXPAND
     +
     +                  ( NXYZT,NCGTO12,NCGTO34,
     +                    NCGTO1,NCGTO2,
     +                    ZCORE (IN),
     +
     +                            ZCORE (OUT) )
     +
     +
C             WRITE (*,*) ' Finished ctr rs expansion '

             I = IN
             IN = OUT
             OUT = I
         END IF

         NCGTO12 = NCGTO1 * NCGTO2

         IF (EQUAL34 .AND. NCGTO34.GT.1) THEN
             CALL  ERD__CTR_TU_EXPAND
     +
     +                  ( NXYZT*NCGTO12,NCGTO34,
     +                    NCGTO3,NCGTO4,
     +                    ZCORE (IN),
     +
     +                            ZCORE (OUT) )
     +
     +
C             WRITE (*,*) ' Finished ctr tu expansion '

             I = IN
             IN = OUT
             OUT = I
         END IF

         NCGTO34 = NCGTO3 * NCGTO4

         IF ((SWAPRS.OR.SWAPTU) .AND. NCTR.GT.1) THEN

              IXOFF (1) = 1
              IXOFF (2) = NCGTO1
              IXOFF (3) = NCGTO1 * NCGTO2
              IXOFF (4) = NCGTO1 * NCGTO2 * NCGTO3

              IF (SWAPRS) THEN
                  INDEXR = 2
                  INDEXS = 1
                  NCGTOR = NCGTO2
                  NCGTOS = NCGTO1
              ELSE
                  INDEXR = 1
                  INDEXS = 2
                  NCGTOR = NCGTO1
                  NCGTOS = NCGTO2
              END IF

              IF (SWAPTU) THEN
                  INDEXT = 4
                  INDEXU = 3
                  NCGTOT = NCGTO4
                  NCGTOU = NCGTO3
              ELSE
                  INDEXT = 3
                  INDEXU = 4
                  NCGTOT = NCGTO3
                  NCGTOU = NCGTO4
              END IF

              CALL  ERD__CTR_4INDEX_REORDER
     +
     +                   ( NXYZT,NCTR,
     +                     NCGTOR,NCGTOS,NCGTOT,NCGTOU,
     +                     IXOFF (INDEXR),IXOFF (INDEXS),
     +                     IXOFF (INDEXT),IXOFF (INDEXU),
     +                     ZCORE (IN),
     +
     +                             ZCORE (OUT) )
     +
     +
C              WRITE (*,*) ' Finished ctr reorder '

              I = IN
              IN = OUT
              OUT = I
         END IF
C
C
C             ...reorder contracted (12|34) batch:
C
C                      batch (nxyz1,nxyz2,nxyz3,nxyz4,rstu) -->
C                               batch (nxyz1,r,nxyz2,s,nxyz3,t,nxyz4,u)
C
C                Do this in three steps (if necessary):
C
C                   i) batch (nxyz1,nxyz2,nxyz3,nxyz4,rstu) -->
C                               batch (nxyz1,nxyz2,nxyz3,rst,nxyz4,u)
C
C                  ii) batch (nxyz1,nxyz2,nxyz3,rst,nxyz4,u) -->
C                               batch (nxyz1,nxyz2,rs,nxyz3,t,nxyz4,u)
C
C                 iii) batch (nxyz1,nxyz2,rs,nxyz3,t,nxyz4,u) -->
C                               batch (nxyz1,r,nxyz2,s,nxyz3,t,nxyz4,u)
C
C
         XTMOVE = NXYZ2 * NXYZ3 * NXYZ4
         CTMOVE = NCGTO12 * NCGTO3

         IF (XTMOVE.GT.1 .OR. CTMOVE.GT.1) THEN

             IF (NXYZ4.GT.1) THEN
                 CALL  ERD__MAP_IJKL_TO_IKJL
     +
     +                      ( NXYZ1*NXYZ2*NXYZ3,
     +                        NXYZ4,
     +                        CTMOVE,
     +                        NCGTO4,
     +                        TILE,
     +                        ZCORE (IN),
     +
     +                                 ZCORE (OUT) )
     +
     +
                 I = IN
                 IN = OUT
                 OUT = I
             END IF

             IF (NXYZ3.GT.1 .AND. NCGTO12.GT.1) THEN
                 CALL  ERD__MAP_IJKL_TO_IKJL
     +
     +                      ( NXYZ1*NXYZ2,
     +                        NXYZ3,
     +                        NCGTO12,
     +                        NXYZ4*NCGTO34,
     +                        TILE,
     +                        ZCORE (IN),
     +
     +                                 ZCORE (OUT) )
     +
     +
                 I = IN
                 IN = OUT
                 OUT = I
             END IF

             IF (NXYZ2.GT.1 .AND. NCGTO1.GT.1) THEN
                 CALL  ERD__MAP_IJKL_TO_IKJL
     +
     +                      ( NXYZ1,
     +                        NXYZ2,
     +                        NCGTO1,
     +                        NXYZ3*NXYZ4*NCGTO2*NCGTO34,
     +                        TILE,
     +                        ZCORE (IN),
     +
     +                                 ZCORE (OUT) )
     +
     +
                 I = IN
                 IN = OUT
                 OUT = I
             END IF

         END IF
C
C
C             ...set final pointer to integrals in ZCORE array.
C
C
         NFIRST = IN
C
C
C             ...ready!
C
C
         RETURN
         END
