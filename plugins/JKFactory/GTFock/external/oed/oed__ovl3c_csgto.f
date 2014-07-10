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
         SUBROUTINE  OED__OVL3C_CSGTO
     +
     +                    ( IMAX,ZMAX,
     +                      NALPHA,NCOEFF,NCSUM,
     +                      NCGTO1,NCGTO2,NCGTO3,
     +                      NPGTO1,NPGTO2,NPGTO3,
     +                      SHELL1,SHELL2,SHELL3,
     +                      X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,
     +                      ALPHA,CC,CCBEG,CCEND,
     +                      L1CACHE,TILE,NCTROW,
     +                      SPHERIC,SCREEN,
     +                      ICORE,
     +
     +                                NBATCH,
     +                                NFIRST,
     +                                ZCORE )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__OVL3C_CSGTO
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__OVL3C_SET_ABC
C                OED__OVL3C_SET_IJK_EXPS
C                OED__OVL3C_F00_DEF_BLOCKS
C                OED__OVL3C_PREPARE_CTR
C                OED__OVL3C_F00_PCGTO_BLOCK
C                OED__CTR_3INDEX_BLOCK
C                OED__CTR_RS_EXPAND
C                OED__CTR_3INDEX_REORDER
C                OED__TRANSPOSE_BATCH
C                OED__XYZ_TO_RY_ABC
C                OED__CARTESIAN_NORMS
C                OED__HRR_MATRIX
C                OED__HRR_TRANSFORM
C                OED__SPHERICAL_TRANSFORM
C                OED__NORMALIZE_CARTESIAN
C                OED__MOVE_RY
C  DESCRIPTION : This operation calculates a batch of contracted
C                3-center overlap integrals on up to three different
C                centers between spherical or cartesian gaussian type
C                shells.
C
C
C                  Input (x = 1,2 and 3):
C
C                    IMAX,ZMAX    =  maximum integer + flp memory
C                    NALPHA       =  total # of exponents
C                    NCOEFF       =  total # of contraction coeffs
C                    NCSUM        =  total # of contractions
C                    NCGTOx       =  # of contractions for csh x
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x
C                    SHELLx       =  the shell type for csh x
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = 1,2 and 3
C                    ALPHA        =  primitive exponents for csh
C                                    1,2 and 3 in that order
C                    CC           =  full set (including zeros) of
C                                    contraction coefficients for csh
C                                    1,2 and 3 in that order, for each
C                                    csh individually such that an
C                                    (I,J) element corresponds to the
C                                    I-th primitive and J-th contraction.
C                    CC(BEG)END   =  (lowest)highest nonzero primitive
C                                    index for contractions for csh
C                                    1,2 and 3 in that order. They are
C                                    different from (1)NPGTOx only for
C                                    segmented contractions
C                    L1CACHE      =  Size of level 1 cache in units of
C                                    8 Byte
C                    TILE         =  Number of rows and columns in
C                                    units of 8 Byte of level 1 cache
C                                    square tile array used for
C                                    performing optimum matrix
C                                    transpositions
C                    NCTROW       =  minimum # of rows that are
C                                    accepted for blocked contractions
C                    SPHERIC      =  is true, if spherical integrals
C                                    are wanted, false if cartesian
C                                    ones are wanted
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
C                    ZCORE        =  full batch of contracted (123)
C                                    3-center overlap integrals over
C                                    cartesian or spherical gaussians
C                                    starting at ZCORE (NFIRST)
C
C
C
C      --- NOTES ABOUT THE OVERALL 3-CENTER OVLP INTEGRAL PREFACTOR ---
C
C                The overal 3-center overlap integral prefactor is
C                defined here as follows. Consider the normalization
C                factors for a primitive cartesian GTO and for a
C                spherical GTO belonging to the angular momentum
C                L = l+m+n:
C
C
C                    lmn                        l m n           2
C                 GTO  (x,y,z) = N (l,m,n,a) * x y z  * exp (-ar )
C                    a
C
C
C                    LM                     L    LM                2
C                 GTO  (r,t,p) = N (L,a) * r  * Y  (t,p) * exp (-ar )
C                    a
C
C
C                where a = alpha exponent, t = theta and p = phi and
C                N (l,m,n,a) and N (L,a) denote the respective
C                cartesian and spherical normalization factors such
C                that:
C
C
C                              lmn            lmn
C                 integral {GTO  (x,y,z) * GTO   (x,y,z) dx dy dz} = 1
C                              a              a
C
C
C                              LM             LM
C                 integral {GTO  (r,t,p) * GTO  (r,t,p) dr dt dp} = 1
C                              a              a
C
C
C                The normalization constants have then the following
C                values, assuming the spherical harmonics are
C                normalized:
C
C                              _____________________________________
C                             /      2^(2L+1+1/2) * a^((2L+3)/2)
C            N (l,m,n,a) =   / ----------------------------------------
C                          \/ (2l-1)!!(2m-1)!!(2n-1)!! * pi * sqrt (pi)
C
C
C                                   ____________________________
C                                  / 2^(2L+3+1/2) * a^((2L+3)/2)
C                     N (L,a) =   / -----------------------------
C                               \/     (2L+1)!! * sqrt (pi)
C
C
C                Note, that the extra pi under the square root in
C                N (l,m,n,a) belongs to the normalization of the
C                spherical harmonic functions and therefore does not
C                appear in the expression for N (L,a). The common
C                L-,l-,m-,n- and a-independent part of the cartesian
C                norm is a scalar quantity needed for all integrals
C                no matter what L-,l-,m-,n- and a-values they have:
C
C                                       _____________
C                                      /  2^(1+1/2)
C                     N (0,0,0,0) =   / --------------
C                                   \/  pi * sqrt (pi)
C
C
C                Also every 3-center overlap integral has a factor of
C                pi**(3/2) associated with it, hence the overall common
C                factor for all 3-center overlap integrals will be
C                N(0,0,0,0)**3 times pi**(3/2), which is equal to:
C
C                                       _____________
C                                      /   sqrt (2)
C                    PREFACT  =  4 *  / --------------
C                                   \/  pi * sqrt (pi)
C
C
C                and is set as a parameter inside the present routine.
C                The alpha exponent dependent part of the norms:
C
C                                   ______________
C                                 \/ a^((2L+3)/2)
C
C                will be calculated separately (see below) and their
C                inclusion in evaluating the primitive cartesian
C                3-center overlap [F00] integrals will be essential
C                for numerical stability during contraction.
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     ATOMIC
         LOGICAL     BLOCKED
         LOGICAL     EMPTY
         LOGICAL     EQUALAB,EQUALAC,EQUALBC
         LOGICAL     EQUALIJ,EQUALIK,EQUALJK
         LOGICAL     MEMORY
         LOGICAL     REORDER
         LOGICAL     SCREEN
         LOGICAL     SPHERIC
         LOGICAL     SWAP12,SWAP13,SWAP23
         LOGICAL     SWAPAC,SWAPBC
         LOGICAL     SWAPRS

         INTEGER     IHNROW,IHROW,IHSCR
         INTEGER     IMAX,ZMAX
         INTEGER     IN,OUT
         INTEGER     INDEXA,INDEXB,INDEXC
         INTEGER     INDEXI,INDEXJ,INDEXK
         INTEGER     INDEXR,INDEXS,INDEXT
         INTEGER     IPRIMI,IPRIMJ,IPRIMK
         INTEGER     IPUSED,IPSAVE,IPPAIR
         INTEGER     ISNROWA,ISNROWB,ISNROWC
         INTEGER     ISROWA,ISROWB,ISROWC
         INTEGER     IUSED,ZUSED
         INTEGER     L1CACHE,TILE,NCTROW
         INTEGER     LCC1,LCC2,LCC3
         INTEGER     LCCA,LCCB,LCCC
         INTEGER     LCCI,LCCJ,LCCK
         INTEGER     LCCSEGA,LCCSEGB,LCCSEGC
         INTEGER     LCCSEGI,LCCSEGJ,LCCSEGK
         INTEGER     LEXP1,LEXP2,LEXP3
         INTEGER     LEXPA,LEXPB,LEXPC
         INTEGER     LEXPI,LEXPJ,LEXPK
         INTEGER     MIJ,MIJK
         INTEGER     MOVE,NOTMOVE
         INTEGER     MXPRIM,MNPRIM
         INTEGER     MXSHELL,MXSIZE
         INTEGER     NABCOOR,NACCOOR
         INTEGER     NALPHA,NCOEFF,NCSUM
         INTEGER     NBATCH,NFIRST
         INTEGER     NCGTO1,NCGTO2,NCGTO3
         INTEGER     NCGTOA,NCGTOB,NCGTOC
         INTEGER     NCGTOI,NCGTOJ,NCGTOK,NCGTOIJ
         INTEGER     NCGTOR,NCGTOS,NCGTOT
         INTEGER     NCOLHRR,NROWHRR,NROTHRR,NXYZHRR
         INTEGER     NCTR
         INTEGER     NIJ,NIJBLK,NIJBEG,NIJEND,NK,NIJK
         INTEGER     NINT1D
         INTEGER     NPGTO1,NPGTO2,NPGTO3
         INTEGER     NPGTOA,NPGTOB,NPGTOC
         INTEGER     NPGTOI,NPGTOJ,NPGTOK,NPGTOIJ
         INTEGER     NPSIZE,NCSIZE,NWSIZE
         INTEGER     NROTA,NROTB,NROTC
         INTEGER     NROWA,NROWB,NROWC
         INTEGER     NRYA,NRYB,NRYC
         INTEGER     NXYZA,NXYZB,NXYZC
         INTEGER     NXYZET,NXYZFT,NXYZP,NXYZQ
         INTEGER     POS1,POS2
         INTEGER     SHELL1,SHELL2,SHELL3
         INTEGER     SHELLA,SHELLB,SHELLC
         INTEGER     SHELLI,SHELLJ,SHELLK
         INTEGER     SHELLP,SHELLQ
         INTEGER     TEMP

         INTEGER     ZCBATCH,ZPBATCH,ZWORK,
     +               ZNORMI,ZNORMJ,ZNORMK,
     +               ZBASE,ZCNORM,
     +               ZRHOIJK,
     +               ZQAX,ZQAY,ZQAZ,ZQINVHF,ZSCALE,
     +               ZINT1DX,ZINT1DY,ZINT1DZ,
     +               ZSROTA,ZSROTB,ZSROTC,
     +               ZHROT

         INTEGER     CCBEG (1:NCSUM)
         INTEGER     CCEND (1:NCSUM)
         INTEGER     ICORE (1:IMAX)
         INTEGER     IXOFF (1:3)

         DOUBLE PRECISION  ABX,ABY,ABZ,ACX,ACY,ACZ
         DOUBLE PRECISION  PREFACT
         DOUBLE PRECISION  RNABSQ,RNACSQ,RNBCSQ
         DOUBLE PRECISION  RNIJSQ,RNIKSQ,RNJKSQ
         DOUBLE PRECISION  SPNORM
         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB,XC,YC,ZC
         DOUBLE PRECISION  XI,YI,ZI,XJ,YJ,ZJ,XK,YK,ZK

         DOUBLE PRECISION  ALPHA (1:NALPHA)
         DOUBLE PRECISION  CC    (1:NCOEFF)
         DOUBLE PRECISION  ZCORE (1:ZMAX)

         PARAMETER  (PREFACT = 2.015835484307046D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...fix the A,B,C labels from the 1,2,3 ones. Calculate
C                the relevant data for the A,B,C batch of 3-center
C                overlap integrals.
C
C
         LEXP1 = 1
         LEXP2 = LEXP1 + NPGTO1
         LEXP3 = LEXP2 + NPGTO2

         LCC1  = 1
         LCC2  = LCC1 + NPGTO1 * NCGTO1
         LCC3  = LCC2 + NPGTO2 * NCGTO2

         CALL  OED__OVL3C_SET_ABC
     +
     +              ( NCGTO1,NCGTO2,NCGTO3,
     +                NPGTO1,NPGTO2,NPGTO3,
     +                SHELL1,SHELL2,SHELL3,
     +                X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,
     +                ALPHA (LEXP1),ALPHA (LEXP2),ALPHA (LEXP3),
     +                CC (LCC1),CC (LCC2),CC (LCC3),
     +                SPHERIC,
     +
     +                            NCGTOA,NCGTOB,NCGTOC,
     +                            NPGTOA,NPGTOB,NPGTOC,
     +                            SHELLA,SHELLB,SHELLC,SHELLP,SHELLQ,
     +                            MXSHELL,
     +                            XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,
     +                            ATOMIC,EQUALAB,EQUALBC,EQUALAC,
     +                            ABX,ABY,ABZ,ACX,ACY,ACZ,
     +                            NABCOOR,NACCOOR,
     +                            RNABSQ,RNACSQ,RNBCSQ,
     +                            SPNORM,
     +                            NXYZA,NXYZB,NXYZC,
     +                            NXYZET,NXYZFT,NXYZP,NXYZQ,
     +                            NRYA,NRYB,NRYC,
     +                            INDEXA,INDEXB,INDEXC,
     +                            SWAP12,SWAP23,SWAP13,
     +                            LEXPA,LEXPB,LEXPC,
     +                            LCCA,LCCB,LCCC,
     +                            LCCSEGA,LCCSEGB,LCCSEGC,
     +                            NXYZHRR,NCOLHRR,NROTHRR,
     +                            EMPTY )
     +
     +
C         WRITE (*,*) ' Finished set abc '
C         WRITE (*,*) ' Index A,B,C = ',INDEXA,INDEXB,INDEXC

         IF (EMPTY) THEN
             NBATCH = 0
             RETURN
         END IF
C
C
C             ...enter the cartesian contracted (f00) 3-center
C                overlap batch generation. Set the ij and k
C                primitive exponent sets and the corresponding
C                exponential prefactors.
C
C
         CALL  OED__OVL3C_SET_IJK
     +
     +              ( NCGTOA,NCGTOB,NCGTOC,
     +                NPGTOA,NPGTOB,NPGTOC,
     +                SHELLA,SHELLB,SHELLC,
     +                EQUALAB,EQUALAC,EQUALBC,
     +                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,
     +                RNABSQ,RNACSQ,RNBCSQ,
     +                INDEXA,INDEXB,INDEXC,
     +                LEXPA,LEXPB,LEXPC,
     +                LCCA,LCCB,LCCC,
     +                LCCSEGA,LCCSEGB,LCCSEGC,
     +
     +                         NCGTOI,NCGTOJ,NCGTOK,NCGTOIJ,
     +                         NPGTOI,NPGTOJ,NPGTOK,NPGTOIJ,
     +                         SHELLI,SHELLJ,SHELLK,
     +                         EQUALIJ,EQUALIK,EQUALJK,
     +                         XI,YI,ZI,XJ,YJ,ZJ,XK,YK,ZK,
     +                         RNIJSQ,RNIKSQ,RNJKSQ,
     +                         SWAPAC,SWAPBC,SWAPRS,
     +                         INDEXI,INDEXJ,INDEXK,
     +                         LEXPI,LEXPJ,LEXPK,
     +                         LCCI,LCCJ,LCCK,
     +                         LCCSEGI,LCCSEGJ,LCCSEGK )
     +
     +
         IPRIMI = 1
         IPRIMJ = IPRIMI + NPGTOIJ
         IPRIMK = IPRIMJ + NPGTOIJ

         CALL  OED__OVL3C_SET_IJK_TRIPLES
     +
     +              ( NPGTOI,NPGTOJ,NPGTOK,NPGTOIJ,
     +                ATOMIC,EQUALIJ,
     +                SWAPRS,
     +                RNIJSQ,RNIKSQ,RNJKSQ,
     +                PREFACT,
     +                ALPHA (LEXPI),ALPHA (LEXPJ),ALPHA (LEXPK),
     +                SCREEN,
     +
     +                         EMPTY,
     +                         NIJ,NK,NIJK,
     +                         ICORE (IPRIMI),
     +                         ICORE (IPRIMJ),
     +                         ICORE (IPRIMK),
     +                         ZCORE (1) )
     +
     +
         IF (EMPTY) THEN
             NBATCH = 0
             RETURN
         END IF
C
C
C             ...decide on the primitive [f00] block size and
C                return array sizes and pointers for the primitive
C                [f00] generation. Perform also some preparation
C                steps for contraction.
C
C
         NCTR = NCGTOIJ * NCGTOK

         MEMORY = .FALSE.

         CALL  OED__OVL3C_F00_DEF_BLOCKS
     +
     +              ( ZMAX,
     +                NPGTOI,NPGTOJ,NPGTOK,
     +                SHELLQ,
     +                NIJ,NK,NIJK,
     +                NCGTOIJ,NCGTOK,NCTR,
     +                NXYZFT,
     +                L1CACHE,NCTROW,
     +                MEMORY,
     +
     +                        NIJBLK,
     +                        NPSIZE,NCSIZE,NWSIZE,
     +                        NINT1D,
     +                        MXPRIM,MNPRIM,
     +                        ZCBATCH,ZPBATCH,ZWORK,
     +                        ZNORMI,ZNORMJ,ZNORMK,
     +                        ZRHOIJK,
     +                        ZQAX,ZQAY,ZQAZ,ZQINVHF,ZSCALE,
     +                        ZINT1DX,ZINT1DY,ZINT1DZ )
     +
     +
         BLOCKED = NIJBLK .LT. NIJ

         CALL  OED__OVL3C_PREPARE_CTR
     +
     +              ( NCSIZE,
     +                NIJK,
     +                NPGTOI,NPGTOJ,NPGTOK,
     +                SHELLI,SHELLJ,SHELLK,
     +                ALPHA (LEXPI),ALPHA (LEXPJ),ALPHA (LEXPK),
     +                PREFACT,SPNORM,
     +                EQUALIJ,EQUALIK,EQUALJK,
     +                BLOCKED,
     +                ZCORE (1),
     +
     +                          ZCORE (ZNORMI),
     +                          ZCORE (ZNORMJ),
     +                          ZCORE (ZNORMK),
     +                          ZCORE (ZRHOIJK),
     +                          ZCORE (ZCBATCH) )
     +
     +
         IPUSED = IPRIMK + NPGTOK
         IPSAVE = IPUSED + MNPRIM
         IPPAIR = IPSAVE + MXPRIM
C
C
C             ...evaluate unnormalized rescaled [f00] 3-center
C                overlap integrals in blocks over ij pairs and
C                add to final contracted (f00) 3-center overlap
C                integrals. The keyword REORDER indicates, if the
C                primitive [f00] 3-center overlap integrals need
C                to be transposed before being contracted.
C
C
         REORDER = .TRUE.

         DO 1000 NIJBEG = 1,NIJ,NIJBLK
            NIJEND = MIN0 (NIJBEG+NIJBLK-1,NIJ)
            MIJ = NIJEND - NIJBEG + 1
            MIJK = MIJ * NK

            CALL  OED__OVL3C_F00_PCGTO_BLOCK
     +
     +                 ( NPSIZE,NINT1D,
     +                   ATOMIC,
     +                   MIJ,NK,MIJK,
     +                   NIJ,NIJBEG,NIJEND,NK,1,NK,
     +                   NPGTOI,NPGTOJ,NPGTOK,
     +                   NXYZFT,NXYZQ,
     +                   SHELLA,SHELLQ,
     +                   XA,YA,ZA,XI,YI,ZI,XJ,YJ,ZJ,XK,YK,ZK,
     +                   ALPHA (LEXPI),ALPHA (LEXPJ),ALPHA (LEXPK),
     +                   ICORE (IPRIMI),ICORE (IPRIMJ),ICORE (IPRIMK),
     +                   ZCORE (ZNORMI),ZCORE (ZNORMJ),ZCORE (ZNORMK),
     +                   ZCORE (ZRHOIJK),
     +                   ZCORE (ZQAX),ZCORE (ZQAY),ZCORE (ZQAZ),
     +                   ZCORE (ZQINVHF),ZCORE (ZSCALE),
     +                   ZCORE (ZINT1DX),
     +                   ZCORE (ZINT1DY),
     +                   ZCORE (ZINT1DZ),
     +
     +                             ZCORE (ZPBATCH) )
     +
     +
C            WRITE (*,*) ' Finished f00 3-center overlap pcgto block '

            CALL  OED__CTR_3INDEX_BLOCK
     +
     +                 ( NPSIZE,NCSIZE,NWSIZE,
     +                   NXYZFT,
     +                   MIJ,NK,MIJK,NCGTOIJ,
     +                   NPGTOI,NPGTOJ,NPGTOK,
     +                   NCGTOI,NCGTOJ,NCGTOK,
     +                   MXPRIM,MNPRIM,
     +                   CC (LCCI),CC (LCCJ),CC (LCCK),
     +                   CCBEG (LCCSEGI),CCBEG (LCCSEGJ),
     +                   CCBEG (LCCSEGK),
     +                   CCEND (LCCSEGI),CCEND (LCCSEGJ),
     +                   CCEND (LCCSEGK),
     +                   ICORE (IPRIMI+NIJBEG-1),
     +                   ICORE (IPRIMJ+NIJBEG-1),
     +                   ICORE (IPRIMK),
     +                   L1CACHE,TILE,NCTROW,
     +                   EQUALIJ,
     +                   SWAPRS,
     +                   REORDER,
     +                   BLOCKED,
     +                   ICORE (IPUSED),
     +                   ICORE (IPSAVE),
     +                   ICORE (IPPAIR),
     +                   ZCORE (ZPBATCH),
     +                   ZCORE (ZWORK),
     +
     +                             ZCORE (ZCBATCH) )
     +
     +
C            WRITE (*,*) ' Finished 3 index ctr block '

 1000    CONTINUE

C         CALL  OED__PRINT_BATCH
C     +
C     +              ( NXYZFT,1,1,1,UNITID,  ZCORE (ZCBATCH) )
C     +
C         STOP
C
C
C             ...the unnormalized cartesian (f00) contracted 3-center
C                overlap batch is ready. Expand the contraction indices
C                (if necessary):
C
C                   batch (nxyzft,k,i>=j) --> batch (nxyzft,k,i,j)
C
C                and reorder the contraction index part (if necessary):
C
C                   batch (nxyzft,k,i,j) --> batch (nxyzet,1,2,3)
C
C                The array IXOFF (x) indicates the total # of indices
C                to the left of x without including the nxyzft-part.
C                For the left most IXOFF value it is convenient to
C                set it equal to 1 instead of 0. Note, that the IXOFF
C                array indicates the true # of indices to the left
C                after! the batch has been transposed (see below) and
C                can be used as initial values when moving the
C                ry-components later on during the HRR and cartesian ->
C                spherical transformation procedure.
C                
C                For efficient application of the HRR contraction
C                scheme we need the batch elements ordered as:
C
C                            batch (1,2,3,nxyzft)
C
C                hence we transpose the batch after the reordering.
C
C                The space partitioning of the flp array for all of
C                these steps will be as follows:
C
C
C                          |  Zone 1  |  Zone 2  |
C
C                in which Zone 1 and 2 are 2 batches of HRR maximum
C                size. This can be done because we have always the
C                following dimension inequality:
C
C                              NXYZFT =< NXYZHRR
C
C
         IXOFF (1) = 1
         IXOFF (2) = NCGTO1
         IXOFF (3) = NCGTO1 * NCGTO2

         NCTR = IXOFF (3) * NCGTO3
         MXSIZE = NCTR * NXYZHRR

         IN = ZCBATCH
         OUT = IN + MXSIZE

         IF (EQUALIJ .AND. NCGTOIJ.GT.1) THEN
             CALL  OED__CTR_RS_EXPAND
     +
     +                  ( NXYZFT*NCGTOK,NCGTOIJ,
     +                    NCGTOI,NCGTOJ,
     +                    ZCORE (IN),
     +
     +                            ZCORE (OUT) )
     +
     +
C             WRITE (*,*) ' Finished ctr rs expansion '

             TEMP = IN
             IN = OUT
             OUT = TEMP
         END IF

         IF (SWAPRS) THEN
             INDEXR = INDEXK
             INDEXS = INDEXJ
             INDEXT = INDEXI
             NCGTOR = NCGTOK
             NCGTOS = NCGTOJ
             NCGTOT = NCGTOI
         ELSE
             INDEXR = INDEXK
             INDEXS = INDEXI
             INDEXT = INDEXJ
             NCGTOR = NCGTOK
             NCGTOS = NCGTOI
             NCGTOT = NCGTOJ
         END IF

         REORDER = .NOT. ((INDEXR.EQ.1).AND.(INDEXS.EQ.2))

         IF (REORDER .AND. NCTR.GT.1) THEN

             CALL  OED__CTR_3INDEX_REORDER
     +
     +                  ( NXYZFT,NCTR,
     +                    NCGTOR,NCGTOS,NCGTOT,
     +                    IXOFF (INDEXR),
     +                    IXOFF (INDEXS),
     +                    IXOFF (INDEXT),
     +                    ZCORE (IN),
     +
     +                            ZCORE (OUT) )
     +
     +
C             WRITE (*,*) ' Finished 3 index ctr reorder '

             TEMP = IN
             IN = OUT
             OUT = TEMP
         END IF

         IF (NXYZFT.GT.1 .AND. NCTR.GT.1) THEN

             CALL  OED__TRANSPOSE_BATCH
     +
     +                  ( NXYZFT,NCTR,
     +                    TILE,
     +                    ZCORE (IN),
     +
     +                            ZCORE (OUT) )
     +
     +
C             WRITE (*,*) ' Finished transpose batch '

             TEMP = IN
             IN = OUT
             OUT = TEMP
         END IF
C
C
C             ...enter the HRR contraction and cartesian -> spherical
C                transformation / cartesian normalization section.
C                The sequence of events is to apply the HRR followed
C                by cartesian -> spherical transformations or cartesian
C                normalizations and immediate final positioning of
C                the finished parts to correspond with the contraction
C                indices 1,2,3. The sequence of events is (where ' means
C                spherical or cartesian normalization and [] means the
C                indices are in correspondence):
C
C                       batch (123,f00) --> batch (123,e0,c)
C                       batch (123,e0,c) --> batch (123,e0,c')
C                       batch (123,e0,c') --> batch (123[c'],e0)
C                       batch (123[c'],e0) --> batch (123[c'],a,b)
C                       batch (123[c'],a,b) --> batch (123[c'],a,b')
C                       batch (123[c'],a,b') --> batch (123[b'c'],a)
C                       batch (123[b'c'],a) --> batch (123[b'c'],a')
C                       batch (123[b'c'],a') --> batch (123[a'b'c'])
C
C
C                The space partitioning of the flp array will be
C                as follows:
C
C
C                 |  Zone 1  |  Zone 2  |  Zone 3  |  Zone 4  |
C
C
C                 Zone 1 and 2:  2 batches of MXSIZE maximum size
C                                (set previously)
C
C                       Zone 3:  cart -> spher transformation data
C                                             or
C                                cartesian normalization factors
C
C                       Zone 4:  HRR contraction data
C
C
C                Determine memory allocation offsets for the entire HRR
C                procedure + cartesian -> spherical transformations or
C                cartesian normalizations and generate the transformation
C                matrices + associated data for those shells > p-shell.
C                The offsets are as follows (x=A,B,C):
C
C                    IN = offset for input HRR batch
C                   OUT = offset for output HRR batch
C
C                ZSROTx = offset for x-part transformation matrix
C               ISNROWx = offset for # of non-zero XYZ contribution row
C                         labels for x-part transformation matrix
C                ISROWx = offset for non-zero XYZ contribution row
C                         labels for x-part transformation matrix
C
C                 ZHROT = offset for HRR transformation matrix
C                IHNROW = offset for # of nonzero row labels for
C                         each HRR matrix column
C                 IHROW = offset for nonzero row labels for the HRR
C                         matrix
C                 IHSCR = integer scratch space for HRR matrix
C                         assembly
C
C
C                In case of s- or p-shells no transformation matrix is
C                generated, hence if we have s- and/or p-shells, then
C                no call to the cartesian -> spherical transformation
C                or cartesian normalization routines needs to be done.
C                Instead all integrals have to be multiplied by a factor
C                SPNORM, which has the following value for each s- and
C                p-shell:
C
C                       For s-type shell  =  1
C
C                       For p-type shell  =  2 * norm for s-type
C
C                This factor was introduced together with the overall
C                prefactor during evaluation of the primitive integrals
C                in order to save multiplications.
C
C
         ZBASE = MAX (IN,OUT) + MXSIZE

         IF (SPHERIC) THEN
             IF (MXSHELL.GT.1) THEN
                 CALL  OED__XYZ_TO_RY_ABC
     +
     +                      ( NXYZA,NXYZB,NXYZC,
     +                        NRYA,NRYB,NRYC,
     +                        SHELLA,SHELLB,SHELLC,
     +                        1,ZBASE,
     +
     +                                 NROWA,NROWB,NROWC,
     +                                 NROTA,NROTB,NROTC,
     +                                 ZSROTA,ZSROTB,ZSROTC,
     +                                 ISNROWA,ISNROWB,ISNROWC,
     +                                 ISROWA,ISROWB,ISROWC,
     +                                 IUSED,ZUSED,
     +                                 ICORE,ZCORE )
     +
     +
C                 WRITE (*,*) ' Finished xyz to ry abc '
             ELSE
                 IUSED = 0
                 ZUSED = 0
             END IF
         ELSE
             IF (MXSHELL.GT.1) THEN
                 ZCNORM = ZBASE
                 CALL  OED__CARTESIAN_NORMS
     +
     +                      ( MXSHELL,
     +
     +                                 ZCORE (ZCNORM))
     +
     +
C                 WRITE (*,*) ' Finished cartesian partial norms '
                 IUSED = 0
                 ZUSED = MXSHELL + 1
             ELSE
                 IUSED = 0
                 ZUSED = 0
             END IF
         END IF

         IHNROW = 1 + IUSED
         IHROW  = IHNROW + NCOLHRR + NCOLHRR
         IHSCR  = IHROW + NROTHRR + NROTHRR
         ZHROT  = ZBASE + ZUSED
C
C
C             ...do first HRR step (if c > 0):
C
C                    (123,f00) --> (123,e=a,a+1,...,a+b,0,c)
C
C
         IF (SHELLC.NE.0) THEN

             CALL  OED__HRR_MATRIX
     +
     +                  ( NROTHRR,NCOLHRR,
     +                    NXYZFT,NXYZA,NXYZQ,
     +                    SHELLA,SHELLC,SHELLQ,
     +                    NACCOOR,
     +                    ACX,ACY,ACZ,
     +                    ICORE (IHSCR),
     +
     +                             POS1,POS2,
     +                             NROWHRR,
     +                             ICORE (IHNROW),
     +                             ICORE (IHROW),
     +                             ZCORE (ZHROT) )
     +
     +
C             WRITE (*,*) ' Finished HRR f00 matrix '

             CALL  OED__HRR_TRANSFORM
     +
     +                  ( NCTR,
     +                    NROWHRR,NXYZFT,NXYZET*NXYZC,
     +                    NXYZET,NXYZC,
     +                    ICORE (IHNROW+POS1-1),
     +                    ICORE (IHROW+POS2-1),
     +                    ZCORE (ZHROT+POS2-1),
     +                    ZCORE (IN),
     +
     +                             ZCORE (OUT) )
     +
     +
C             WRITE (*,*) ' Finished HRR f00 -> e0c '
             TEMP = IN
             IN = OUT
             OUT = TEMP
C
C
C             ...do cart -> spher transformation / cart normalization
C                (if c > 1):
C
C                        (123,e0,c) --> (123,e0,c')
C
C
             IF (SHELLC.GT.1) THEN
                 IF (SPHERIC) THEN
                     CALL  OED__SPHERICAL_TRANSFORM
     +
     +                          ( NCTR*NXYZET,
     +                            NROWC,NXYZC,NRYC,
     +                            ICORE (ISNROWC),
     +                            ICORE (ISROWC),
     +                            ZCORE (ZSROTC),
     +                            ZCORE (IN),
     +
     +                                    ZCORE (OUT) )
     +
     +
C                     WRITE (*,*) ' Finished sph quart c '
                     TEMP = IN
                     IN = OUT
                     OUT = TEMP
                 ELSE
                     CALL  OED__NORMALIZE_CARTESIAN
     +
     +                          ( NCTR*NXYZET,
     +                            NXYZC,
     +                            SHELLC,
     +                            ZCORE (ZCNORM),
     +
     +                                    ZCORE (IN) )
     +
     +
C                     WRITE (*,*) ' Finished normalized cart c '
                 END IF
             END IF
         END IF
C
C
C             ...move transformed c shell (if size > 1):
C
C                        (123,e0,c') --> (123[c'],e0)
C
C
         IF (NRYC.GT.1) THEN
             NBATCH = NCTR * NXYZET * NRYC
             NOTMOVE = IXOFF (INDEXC)
             MOVE = NBATCH / (NOTMOVE * NRYC)
             IF (MOVE.GT.1) THEN

                 CALL  OED__MOVE_RY
     +
     +                      ( NBATCH,3,
     +                        NOTMOVE,MOVE,NRYC,
     +                        INDEXC,
     +                        TILE,
     +                        ZCORE (IN),
     +
     +                                IXOFF,
     +                                ZCORE (OUT) )
     +
     +
C                 WRITE (*,*) ' Finished move ry c '
                 TEMP = IN
                 IN = OUT
                 OUT = TEMP
             END IF
         END IF
C
C
C             ...do second HRR step (if b > 0):
C
C                    (123[c'],e0) --> (123[c'],a,b)
C
C
         IF (SHELLB.NE.0) THEN

             CALL  OED__HRR_MATRIX
     +
     +                  ( NROTHRR,NCOLHRR,
     +                    NXYZET,NXYZA,NXYZP,
     +                    SHELLA,SHELLB,SHELLP,
     +                    NABCOOR,
     +                    ABX,ABY,ABZ,
     +                    ICORE (IHSCR),
     +
     +                             POS1,POS2,
     +                             NROWHRR,
     +                             ICORE (IHNROW),
     +                             ICORE (IHROW),
     +                             ZCORE (ZHROT) )
     +
     +
C             WRITE (*,*) ' Finished HRR e0 matrix '

             CALL  OED__HRR_TRANSFORM
     +
     +                  ( NCTR*NRYC,
     +                    NROWHRR,NXYZET,NXYZA*NXYZB,
     +                    NXYZA,NXYZB,
     +                    ICORE (IHNROW+POS1-1),
     +                    ICORE (IHROW+POS2-1),
     +                    ZCORE (ZHROT+POS2-1),
     +                    ZCORE (IN),
     +
     +                             ZCORE (OUT) )
     +
     +
C             WRITE (*,*) ' Finished HRR e0 -> ab '
             TEMP = IN
             IN = OUT
             OUT = TEMP
C
C
C             ...do cart -> spher transformation / cart normalization
C                (if b > 1):
C
C                        (123[c'],a,b) --> (123[c'],a,b')
C
C
             IF (SHELLB.GT.1) THEN
                 IF (SPHERIC) THEN
                     CALL  OED__SPHERICAL_TRANSFORM
     +
     +                          ( NCTR*NRYC*NXYZA,
     +                            NROWB,NXYZB,NRYB,
     +                            ICORE (ISNROWB),
     +                            ICORE (ISROWB),
     +                            ZCORE (ZSROTB),
     +                            ZCORE (IN),
     +
     +                                    ZCORE (OUT) )
     +
     +
C                     WRITE (*,*) ' Finished sph quart b '
                     TEMP = IN
                     IN = OUT
                     OUT = TEMP
                 ELSE
                     CALL  OED__NORMALIZE_CARTESIAN
     +
     +                          ( NCTR*NRYC*NXYZA,
     +                            NXYZB,
     +                            SHELLB,
     +                            ZCORE (ZCNORM),
     +
     +                                    ZCORE (IN) )
     +
     +
C                     WRITE (*,*) ' Finished normalized cart b '
                 END IF
             END IF
         END IF
C
C
C             ...move transformed b shell (if size > 1):
C
C                        (123[c'],a,b') --> (123[b'c'],a)
C
C
         IF (NRYB.GT.1) THEN
             NBATCH = NCTR * NRYC * NXYZA * NRYB
             NOTMOVE = IXOFF (INDEXB)
             MOVE = NBATCH / (NOTMOVE * NRYB)
             IF (MOVE.GT.1) THEN

                 CALL  OED__MOVE_RY
     +
     +                      ( NBATCH,3,
     +                        NOTMOVE,MOVE,NRYB,
     +                        INDEXB,
     +                        TILE,
     +                        ZCORE (IN),
     +
     +                                IXOFF,
     +                                ZCORE (OUT) )
     +
     +
C                 WRITE (*,*) ' Finished move ry b '
                 TEMP = IN
                 IN = OUT
                 OUT = TEMP
             END IF
         END IF
C
C
C             ...do cart -> spher transformation / cart normalization
C                (if a > 1):
C
C                        (123[b'c'],a) --> (123[b'c'],a')
C
C
         IF (SHELLA.GT.1) THEN
             IF (SPHERIC) THEN
                 CALL  OED__SPHERICAL_TRANSFORM
     +
     +                      ( NCTR*NRYC*NRYB,
     +                        NROWA,NXYZA,NRYA,
     +                        ICORE (ISNROWA),
     +                        ICORE (ISROWA),
     +                        ZCORE (ZSROTA),
     +                        ZCORE (IN),
     +
     +                                ZCORE (OUT) )
     +
     +
C                 WRITE (*,*) ' Finished sph quart a '
                 TEMP = IN
                 IN = OUT
                 OUT = TEMP
             ELSE
                 CALL  OED__NORMALIZE_CARTESIAN
     +
     +                      ( NCTR*NRYC*NRYB,
     +                        NXYZA,
     +                        SHELLA,
     +                        ZCORE (ZCNORM),
     +
     +                                ZCORE (IN) )
     +
     +
C                 WRITE (*,*) ' Finished normalized cart a '
             END IF
         END IF
C
C
C             ...move transformed a shell (if size > 1):
C
C                        (123[b'c'],a') --> (123[a'b'c'])
C
C
         NBATCH = NCTR * NRYC * NRYB * NRYA

         IF (NRYA.GT.1) THEN
             NOTMOVE = IXOFF (INDEXA)
             MOVE = NBATCH / (NOTMOVE * NRYA)
             IF (MOVE.GT.1) THEN

                 CALL  OED__MOVE_RY
     +
     +                      ( NBATCH,3,
     +                        NOTMOVE,MOVE,NRYA,
     +                        INDEXA,
     +                        TILE,
     +                        ZCORE (IN),
     +
     +                                IXOFF,
     +                                ZCORE (OUT) )
     +
     +
C                 WRITE (*,*) ' Finished move ry a '
                 TEMP = IN
                 IN = OUT
                 OUT = TEMP
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
