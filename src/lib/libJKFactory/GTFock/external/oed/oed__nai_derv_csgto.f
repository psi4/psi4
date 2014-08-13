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
         SUBROUTINE  OED__NAI_DERV_CSGTO
     +
     +                    ( IMAX,ZMAX,
     +                      NALPHA,NCOEFF,NCSUM,
     +                      NCGTO1,NCGTO2,
     +                      NPGTO1,NPGTO2,
     +                      SHELL1,SHELL2,
     +                      X1,Y1,Z1,X2,Y2,Z2,
     +                      NUCLEI,
     +                      XN,YN,ZN,NCHARGE,
     +                      IXDERC,
     +                      DER1X,DER1Y,DER1Z,
     +                      DER2X,DER2Y,DER2Z,
     +                      DERCX,DERCY,DERCZ,
     +                      ALPHA,CC,CCBEG,CCEND,
     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                      L1CACHE,TILE,NCTROW,
     +                      SPHERIC,SCREEN,
     +                      ICORE,
     +
     +                                NBATCH,
     +                                NFIRST,
     +                                ZCORE )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_DERV_CSGTO
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__NAI_SET_DERV_AB
C                OED__NAI_SET_DERV_IJC_TRIPLES
C                OED__NAI_SET_DERV_SEQUENCE
C                OED__NAI_DERV_DEF_BLOCKS
C                OED__NAI_PREPARE_CTR
C                OED__NAI_DERV_PCGTO_BLOCK
C                OED__CTR_2INDEX_BLOCK
C                OED__CTR_RS_EXPAND
C                OED__CTR_2INDEX_REORDER
C                OED__TRANSPOSE_BATCH
C                OED__XYZ_TO_RY_AB
C                OED__CARTESIAN_NORMS
C                OED__SPHERICAL_TRANSFORM
C                OED__NORMALIZE_CARTESIAN
C                OED__MOVE_RY
C  DESCRIPTION : This operation calculates a batch of differentiated
C                contracted nuclear attraction integrals for a set of
C                nuclear attraction centers on up to two different
C                electronic centers between spherical or cartesian
C                gaussian type shells.
C
C
C                  Input (x = 1 and 2):
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
C                                    y = 1 and 2
C                    NUCLEI       =  # of nuclear attraction centers
C                    XN,YN,ZN     =  the x,y,z-coordinates for all
C                                    nuclear attraction centers
C                    NCHARGE      =  the nuclear charges for all
C                                    nuclear attraction centers
C                    IXDERC       =  the index of which of the nuclear
C                                    attraction centers is to be
C                                    differentiated. If that index
C                                    corresponds to one of the centers
C                                    1 and/or 2, it will already be
C                                    differentiated along with these
C                                    centers, hence values transmitted
C                                    for DERCX,DERCY,DERCZ are
C                                    irrelevant in that case. If no
C                                    nuclear attraction center is to
C                                    be differentiated besides those
C                                    which are possibly equal to
C                                    centers 1 and/or 2, set this
C                                    index value =< 0
C                    DERyp        =  the order of differentiation on
C                                    centers y = 1 and 2 with respect
C                                    to the p = x,y,z coordinates
C                    DERCp        =  the order of differentiation for
C                                    the IXDERC-th nuclear attraction
C                                    center with respect to the
C                                    p = x,y,z coordinates
C                    ALPHA        =  primitive exponents for csh
C                                    1 and 2 in that order
C                    CC           =  full set (including zeros) of
C                                    contraction coefficients for csh
C                                    1 and 2 in that order, for each
C                                    csh individually such that an
C                                    (I,J) element corresponds to the
C                                    I-th primitive and J-th contraction.
C                    CC(BEG)END   =  (lowest)highest nonzero primitive
C                                    index for contractions for csh
C                                    1 and 2 in that order. They are
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
C                    NBATCH       =  # of derivative integrals in batch
C                    NFIRST       =  first address location inside the
C                                    ZCORE array containing the first
C                                    derivative integral
C                    ZCORE        =  full batch of contracted (1|2)
C                                    derivative nuclear attraction
C                                    integrals over cartesian or
C                                    spherical gaussians starting at
C                                    ZCORE (NFIRST)
C
C
C
C              --- NOTES ABOUT THE OVERALL NAI INTEGRAL PREFACTOR ---
C
C                The overal nuclear attraction integral prefactor is
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
C                Also every nuclear attraction integral has a factor
C                of 2*pi associated with it, hence the overall common
C                factor for all nuclear attraction integrals will be
C                N(0,0,0,0)**2 times 2*pi. In order to avoid sign
C                changes of the nuclear charges, the negative sign of
C                these is incorporated also into the prefactor, which
C                is thus equal to:
C
C                                            _________
C                            PREFACT  =  - \/ 32 / pi 
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
C                derivative nuclear attraction [A|B] integrals will be
C                essential for numerical stability during contraction.
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
         LOGICAL     CASEI,CASEII,CASEIII
         LOGICAL     DIFFA,DIFFB,DIFFC
         LOGICAL     DIFFX,DIFFY,DIFFZ
         LOGICAL     EMPTY
         LOGICAL     EQUALAB
         LOGICAL     MEMORY
         LOGICAL     ONECASE
         LOGICAL     REORDER
         LOGICAL     SCREEN
         LOGICAL     SPHERIC
         LOGICAL     SWAP12,SWAPRS

         INTEGER     DER1X,DER1Y,DER1Z
         INTEGER     DER2X,DER2Y,DER2Z
         INTEGER     DERAX,DERAY,DERAZ
         INTEGER     DERBX,DERBY,DERBZ
         INTEGER     DERCX,DERCY,DERCZ
         INTEGER     IMAX,ZMAX
         INTEGER     IN,OUT
         INTEGER     INDEXA,INDEXB
         INTEGER     INDEXR,INDEXS
         INTEGER     INUCCEN
         INTEGER     IPRIMA,IPRIMB
         INTEGER     IPUSED,IPSAVE,IPPAIR
         INTEGER     ISNROWA,ISNROWB
         INTEGER     ISROWA,ISROWB
         INTEGER     IUSED,ZUSED
         INTEGER     IXAEQB,IXCEQA,IXCEQB
         INTEGER     IXDERC
         INTEGER     L1CACHE,TILE,NCTROW
         INTEGER     LCC1,LCC2
         INTEGER     LCCA,LCCB
         INTEGER     LCCSEGA,LCCSEGB
         INTEGER     LEXP1,LEXP2
         INTEGER     LEXPA,LEXPB
         INTEGER     MGRID,NGRID
         INTEGER     MIJ,MIJCEN,MGQPIJ,MGIJCEN
         INTEGER     MOVE,NOTMOVE
         INTEGER     MXPRIM,MNPRIM
         INTEGER     MXSHELL
         INTEGER     NALPHA,NCOEFF,NCSUM
         INTEGER     NBATCH,NFIRST
         INTEGER     NCENA,NCENB,NCENC
         INTEGER     NCGTO1,NCGTO2
         INTEGER     NCGTOA,NCGTOB,NCGTOAB
         INTEGER     NCGTOR,NCGTOS
         INTEGER     NCTR
         INTEGER     NDERX,NDERY,NDERZ
         INTEGER     NGQP,NMOM,NGQSCR
         INTEGER     NIJ,NIJBLK,NIJBEG,NIJEND
         INTEGER     NINT1DX,NINT1DY,NINT1DZ
         INTEGER     NPGTO1,NPGTO2
         INTEGER     NPGTOA,NPGTOB,NPGTOAB
         INTEGER     NPSIZE,NCSIZE,NWSIZE
         INTEGER     NROTA,NROTB
         INTEGER     NROWA,NROWB
         INTEGER     NRYA,NRYB
         INTEGER     NUCLEI
         INTEGER     NXYZA,NXYZB,NXYZT
         INTEGER     SHELL1,SHELL2
         INTEGER     SHELLA,SHELLB,SHELLP
         INTEGER     TEMP

         INTEGER     ZCBATCH,ZPBATCH,ZWORK,
     +               ZNORMA,ZNORMB,
     +               ZBASE,ZCNORM,
     +               ZRHOAB,
     +               ZPX,ZPY,ZPZ,ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCALE,
     +               ZRTS,ZWTS,ZGQSCR,ZTVAL,
     +               ZR1X,ZR1Y,ZR1Z,ZR2,
     +               ZEXP2A,ZEXP2B,ZEXP2AB,
     +               ZINT1DX,ZINT1DY,ZINT1DZ,
     +               ZSROTA,ZSROTB

         INTEGER     CCBEG (1:NCSUM)
         INTEGER     CCEND (1:NCSUM)
         INTEGER     ICORE (1:IMAX)
         INTEGER     IXOFF (1:2)

         DOUBLE PRECISION  ABX,ABY,ABZ
         DOUBLE PRECISION  PREFACT
         DOUBLE PRECISION  RNABSQ
         DOUBLE PRECISION  SPNORM
         DOUBLE PRECISION  TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB

         DOUBLE PRECISION  XN      (1:NUCLEI)
         DOUBLE PRECISION  YN      (1:NUCLEI)
         DOUBLE PRECISION  ZN      (1:NUCLEI)
         DOUBLE PRECISION  NCHARGE (1:NUCLEI)

         DOUBLE PRECISION  ALPHA (1:NALPHA)
         DOUBLE PRECISION  CC    (1:NCOEFF)
         DOUBLE PRECISION  ZCORE (1:ZMAX)

         DOUBLE PRECISION  FTABLE (0:MGRID,0:NGRID)

         PARAMETER  (PREFACT = -3.191538243211461D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...fix the A,B labels from the 1,2 ones. Calculate
C                the relevant data for the A,B batch of derivative
C                nuclear attraction integrals.
C
C
         LEXP1 = 1
         LEXP2 = LEXP1 + NPGTO1
         LCC1  = 1
         LCC2  = LCC1 + NPGTO1 * NCGTO1

         CALL  OED__NAI_SET_DERV_AB
     +
     +              ( NCGTO1,NCGTO2,
     +                NPGTO1,NPGTO2,
     +                SHELL1,SHELL2,
     +                X1,Y1,Z1,X2,Y2,Z2,
     +                NUCLEI,
     +                XN,YN,ZN,
     +                IXDERC,
     +                DER1X,DER1Y,DER1Z,
     +                DER2X,DER2Y,DER2Z,
     +                DERCX,DERCY,DERCZ,
     +                ALPHA (LEXP1),ALPHA (LEXP2),
     +                CC (LCC1),CC (LCC2),
     +                SPHERIC,
     +
     +                            NCGTOA,NCGTOB,
     +                            NPGTOA,NPGTOB,
     +                            SHELLA,SHELLB,SHELLP,
     +                            MXSHELL,
     +                            XA,YA,ZA,XB,YB,ZB,
     +                            NDERX,NDERY,NDERZ,
     +                            DERAX,DERAY,DERAZ,
     +                            DERBX,DERBY,DERBZ,
     +                            DIFFA,DIFFB,DIFFC,
     +                            DIFFX,DIFFY,DIFFZ,
     +                            IXAEQB,IXCEQA,IXCEQB,
     +                            ATOMIC,EQUALAB,
     +                            ABX,ABY,ABZ,RNABSQ,
     +                            SPNORM,
     +                            NXYZA,NXYZB,NXYZT,
     +                            NRYA,NRYB,
     +                            INDEXA,INDEXB,
     +                            SWAP12,SWAPRS,
     +                            LEXPA,LEXPB,
     +                            LCCA,LCCB,
     +                            LCCSEGA,LCCSEGB,
     +                            EMPTY )
     +
     +
C         WRITE (*,*) ' Finished set derivative ab '
C         WRITE (*,*) ' Index A,B = ',INDEXA,INDEXB

         IF (EMPTY) THEN
             NBATCH = 0
             RETURN
         END IF
C
C
C             ...set the i and j primitive exponent sets, the relevant
C                nuclear attraction centers and the corresponding
C                exponential prefactors.
C
C
         IF (EQUALAB) THEN
             NPGTOAB = (NPGTOA*(NPGTOA+1))/2
             NCGTOAB = (NCGTOA*(NCGTOA+1))/2
         ELSE
             NPGTOAB = NPGTOA * NPGTOB
             NCGTOAB = NCGTOA * NCGTOB
         END IF

         IPRIMA = 1
         IPRIMB = IPRIMA + NPGTOAB
         INUCCEN = IPRIMB + NPGTOAB

         CALL  OED__NAI_SET_DERV_IJC_TRIPLES
     +
     +              ( NUCLEI,
     +                XN,YN,ZN,NCHARGE,
     +                NPGTOA,NPGTOB,NPGTOAB,
     +                DIFFC,IXDERC,
     +                ATOMIC,EQUALAB,
     +                SWAPRS,
     +                XA,YA,ZA,XB,YB,ZB,
     +                RNABSQ,
     +                PREFACT,
     +                ALPHA (LEXPA),ALPHA (LEXPB),
     +                FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                SCREEN,
     +
     +                        EMPTY,
     +                        NIJ,NCENA,NCENB,NCENC,
     +                        ICORE (IPRIMA),ICORE (IPRIMB),
     +                        IXCEQA,IXCEQB,
     +                        ICORE (INUCCEN),
     +                        ZCORE (1) )
     +
     +
C         WRITE (*,*) ' Finished set derivative ijc triples '

         IF (EMPTY) THEN
             NBATCH = 0
             RETURN
         END IF
C
C
C             ...we are now ready to calculate the cartesian contracted
C                (a|b) derivative nuclear attraction batch, which is
C                obtained as a sum over all attraction centers:
C
C                           dxdydz (a|b)  =  sum dxdydz (a|C|b)
C                                             C
C
C                Three cases can arise for each basic (a|C|b) derivative
C                integral, each of which has to be dealt separately and
C                summed into (a|b):
C
C                    Case I    =  2-center integrals of type (a|A|b)
C                    Case II   =  2-center integrals of type (a|B|b)
C                    Case III  =  3-center integrals of type (a|C|b)
C
C                Note, that cases I and II cannot arise once dxdydz
C                operates on a center C different from A and B. In that
C                case we only have case III integrals.
C
C                Cases I,II and III differ mainly in the necessary
C                initial shell dimensions for the 1D integrals. We
C                have:
C
C                    Case I   : dxdydz operator will be applied on
C                               the b shell only => initial shell
C                               sizes (x-component as example):
C
C                                     a shell = SHELLA
C                                     b shell = SHELLB + DERAX + DERBX
C
C                    Case II  : dxdydz operator will be applied on
C                               the a shell only => initial shell
C                               sizes (x-component as example):
C
C                                     a shell = SHELLA + DERAX + DERBX
C                                     b shell = SHELLB
C
C                    Case III : dxdydz operator will be applied on
C                               both a and b shells => initial shell
C                               sizes (x-component as example):
C
C                                     a shell = SHELLA + DERAX
C                                     b shell = SHELLB + DERBX
C
C
C                Determine next which of the cases are present and
C                decide on the primitive [a|b] block size and return
C                array sizes and pointers for the primitive [a|b]
C                generation. Perform also some preparation steps
C                for contraction.
C
C
         CASEI   = NCENA .GT. 0
         CASEII  = NCENB .GT. 0
         CASEIII = NCENC .GT. 0

C         WRITE (*,*) ' CASEI   = ',CASEI
C         WRITE (*,*) ' CASEII  = ',CASEII
C         WRITE (*,*) ' CASEIII = ',CASEIII

         ONECASE =      (CASEI  .AND.(.NOT.CASEII).AND.(.NOT.CASEIII))
     +             .OR. (CASEII .AND.(.NOT.CASEI) .AND.(.NOT.CASEIII))
     +             .OR. (CASEIII.AND.(.NOT.CASEI) .AND.(.NOT.CASEII) )

         NGQP = 1 + (SHELLP + NDERX + NDERY + NDERZ) / 2
         NMOM = 2 * NGQP - 1
         NGQSCR = 5 * NMOM + 2 * NGQP - 2

         MEMORY = .FALSE.

         CALL  OED__NAI_DERV_DEF_BLOCKS
     +
     +              ( ZMAX,
     +                NPGTOA,NPGTOB,
     +                SHELLA,SHELLB,SHELLP,
     +                NIJ,NCGTOAB,
     +                MAX (NCENA,NCENB,NCENC),
     +                NGQP,NGQSCR,
     +                NXYZT,
     +                DERAX,DERAY,DERAZ,
     +                DERBX,DERBY,DERBZ,
     +                CASEI,CASEII,CASEIII,
     +                L1CACHE,NCTROW,
     +                MEMORY,
     +
     +                       NIJBLK,
     +                       NBATCH,
     +                       NPSIZE,NCSIZE,NWSIZE,
     +                       NINT1DX,NINT1DY,NINT1DZ,
     +                       MXPRIM,MNPRIM,
     +                       ZCBATCH,ZPBATCH,ZWORK,
     +                       ZNORMA,ZNORMB,
     +                       ZRHOAB,
     +                       ZPX,ZPY,ZPZ,ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCALE,
     +                       ZRTS,ZWTS,ZGQSCR,ZTVAL,
     +                       ZR1X,ZR1Y,ZR1Z,ZR2,
     +                       ZEXP2A,ZEXP2B,ZEXP2AB,
     +                       ZINT1DX,ZINT1DY,ZINT1DZ )
     +
     +
         BLOCKED = (.NOT.ONECASE) .OR. (NIJBLK .LT. NIJ)

C         WRITE (*,*) ' BLOCKED = ',BLOCKED

         CALL  OED__NAI_PREPARE_CTR
     +
     +              ( NCSIZE,
     +                NIJ,
     +                NPGTOA,NPGTOB,
     +                SHELLA,SHELLB,
     +                ALPHA (LEXPA),ALPHA (LEXPB),
     +                PREFACT,SPNORM,
     +                EQUALAB,
     +                BLOCKED,
     +                ZCORE (1),
     +
     +                        ZCORE (ZNORMA),ZCORE (ZNORMB),
     +                        ZCORE (ZRHOAB),
     +                        ZCORE (ZCBATCH) )
     +
     +
         IPUSED = INUCCEN + NCENC
         IPSAVE = IPUSED + MNPRIM
         IPPAIR = IPSAVE + MXPRIM
C
C
C             ...evaluate unnormalized rescaled [a|b] derivative
C                nuclear attraction integrals in blocks over ij pairs
C                and add to final contracted (a|b) derivative nuclear
C                attraction integrals. The keyword REORDER indicates,
C                if the primitive [a|b] derivative nuclear attraction
C                integrals need to be transposed before being
C                contracted.
C
C
         REORDER = .TRUE.
C
C
C             ...case I integrals first (if any). Note, that case I
C                integrals can only arise if both A and B centers
C                differ!
C
C
         IF (CASEI) THEN

             DO 1000 NIJBEG = 1,NIJ,NIJBLK
                NIJEND = MIN0 (NIJBEG+NIJBLK-1,NIJ)
                MIJ = NIJEND - NIJBEG + 1
                MGQPIJ = NGQP * MIJ

                CALL  OED__NAI_DERV_2CEN_PCGTO_BLOCK
     +
     +                     ( NBATCH,
     +                       NINT1DX,NINT1DY,NINT1DZ,
     +                       MIJ,
     +                       NIJ,NIJBEG,NIJEND,
     +                       NGQP,NMOM,NGQSCR,MGQPIJ,
     +                       NPGTOA,NPGTOB,
     +                       NXYZA,NXYZB,
     +                       SHELLA,SHELLB,SHELLP,
     +                       XA,YA,ZA,XB,YB,ZB,
     +                       ABX,ABY,ABZ,
     +                       'A',NCHARGE (IXCEQA),
     +                       NDERX,NDERY,NDERZ,
     +                       DERAX,DERAY,DERAZ,
     +                       DERBX,DERBY,DERBZ,
     +                       DIFFX,DIFFY,DIFFZ,
     +                       ALPHA (LEXPA),ALPHA (LEXPB),
     +                       FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                       ICORE (IPRIMA+NIJBEG-1),
     +                       ICORE (IPRIMB+NIJBEG-1),
     +                       ZCORE (ZNORMA),ZCORE (ZNORMB),
     +                       ZCORE (ZRHOAB),
     +                       ZCORE (ZPX),ZCORE (ZPY),ZCORE (ZPZ),
     +                       ZCORE (ZPAX),ZCORE (ZPAY),ZCORE (ZPAZ),
     +                       ZCORE (ZPINVHF),ZCORE (ZSCALE),
     +                       ZCORE (ZRTS),ZCORE (ZWTS),ZCORE (ZGQSCR),
     +                       ZCORE (ZTVAL),
     +                       ZCORE (ZR1X),ZCORE (ZR1Y),ZCORE (ZR1Z),
     +                       ZCORE (ZR2),
     +                       ZCORE (ZEXP2A),ZCORE (ZEXP2B),
     +                       ZCORE (ZINT1DX),
     +                       ZCORE (ZINT1DY),
     +                       ZCORE (ZINT1DZ),
     +
     +                                 ZCORE (ZPBATCH) )
     +
     +
C            WRITE (*,*) ' Finished case I ab nai pcgto derv block '

                CALL  OED__CTR_2INDEX_BLOCK
     +
     +                     ( NPSIZE,NCSIZE,NWSIZE,
     +                       NXYZT,
     +                       MIJ,NCGTOAB,
     +                       NPGTOA,NPGTOB,
     +                       NCGTOA,NCGTOB,
     +                       MXPRIM,MNPRIM,
     +                       CC (LCCA),CC (LCCB),
     +                       CCBEG (LCCSEGA),CCBEG (LCCSEGB),
     +                       CCEND (LCCSEGA),CCEND (LCCSEGB),
     +                       ICORE (IPRIMA+NIJBEG-1),
     +                       ICORE (IPRIMB+NIJBEG-1),
     +                       L1CACHE,TILE,NCTROW,
     +                       EQUALAB,
     +                       SWAPRS,
     +                       REORDER,
     +                       BLOCKED,
     +                       ICORE (IPUSED),
     +                       ICORE (IPSAVE),
     +                       ICORE (IPPAIR),
     +                       ZCORE (ZPBATCH),
     +                       ZCORE (ZWORK),
     +
     +                                 ZCORE (ZCBATCH) )
     +
     +
C                WRITE (*,*) ' Finished case I 2 index ctr block '

 1000        CONTINUE

         END IF
C
C
C             ...case II integrals next (if any). Note, that case II
C                integrals can only arise if both A and B centers
C                differ!
C
C
         IF (CASEII) THEN

             DO 2000 NIJBEG = 1,NIJ,NIJBLK
                NIJEND = MIN0 (NIJBEG+NIJBLK-1,NIJ)
                MIJ = NIJEND - NIJBEG + 1
                MGQPIJ = NGQP * MIJ

                CALL  OED__NAI_DERV_2CEN_PCGTO_BLOCK
     +
     +                     ( NBATCH,
     +                       NINT1DX,NINT1DY,NINT1DZ,
     +                       MIJ,
     +                       NIJ,NIJBEG,NIJEND,
     +                       NGQP,NMOM,NGQSCR,MGQPIJ,
     +                       NPGTOA,NPGTOB,
     +                       NXYZA,NXYZB,
     +                       SHELLA,SHELLB,SHELLP,
     +                       XA,YA,ZA,XB,YB,ZB,
     +                       ABX,ABY,ABZ,
     +                       'B',NCHARGE (IXCEQB),
     +                       NDERX,NDERY,NDERZ,
     +                       DERAX,DERAY,DERAZ,
     +                       DERBX,DERBY,DERBZ,
     +                       DIFFX,DIFFY,DIFFZ,
     +                       ALPHA (LEXPA),ALPHA (LEXPB),
     +                       FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                       ICORE (IPRIMA+NIJBEG-1),
     +                       ICORE (IPRIMB+NIJBEG-1),
     +                       ZCORE (ZNORMA),ZCORE (ZNORMB),
     +                       ZCORE (ZRHOAB),
     +                       ZCORE (ZPX),ZCORE (ZPY),ZCORE (ZPZ),
     +                       ZCORE (ZPAX),ZCORE (ZPAY),ZCORE (ZPAZ),
     +                       ZCORE (ZPINVHF),ZCORE (ZSCALE),
     +                       ZCORE (ZRTS),ZCORE (ZWTS),ZCORE (ZGQSCR),
     +                       ZCORE (ZTVAL),
     +                       ZCORE (ZR1X),ZCORE (ZR1Y),ZCORE (ZR1Z),
     +                       ZCORE (ZR2),
     +                       ZCORE (ZEXP2A),ZCORE (ZEXP2B),
     +                       ZCORE (ZINT1DX),
     +                       ZCORE (ZINT1DY),
     +                       ZCORE (ZINT1DZ),
     +
     +                                 ZCORE (ZPBATCH) )
     +
     +
C            WRITE (*,*) ' Finished case II ab nai pcgto derv block '

                CALL  OED__CTR_2INDEX_BLOCK
     +
     +                     ( NPSIZE,NCSIZE,NWSIZE,
     +                       NXYZT,
     +                       MIJ,NCGTOAB,
     +                       NPGTOA,NPGTOB,
     +                       NCGTOA,NCGTOB,
     +                       MXPRIM,MNPRIM,
     +                       CC (LCCA),CC (LCCB),
     +                       CCBEG (LCCSEGA),CCBEG (LCCSEGB),
     +                       CCEND (LCCSEGA),CCEND (LCCSEGB),
     +                       ICORE (IPRIMA+NIJBEG-1),
     +                       ICORE (IPRIMB+NIJBEG-1),
     +                       L1CACHE,TILE,NCTROW,
     +                       EQUALAB,
     +                       SWAPRS,
     +                       REORDER,
     +                       BLOCKED,
     +                       ICORE (IPUSED),
     +                       ICORE (IPSAVE),
     +                       ICORE (IPPAIR),
     +                       ZCORE (ZPBATCH),
     +                       ZCORE (ZWORK),
     +
     +                                 ZCORE (ZCBATCH) )
     +
     +
C                WRITE (*,*) ' Finished case II 2 index ctr block '

 2000        CONTINUE

         END IF
C
C
C             ...case III integrals last (if any).
C
C
         IF (CASEIII) THEN

             DO 3000 NIJBEG = 1,NIJ,NIJBLK
                NIJEND = MIN0 (NIJBEG+NIJBLK-1,NIJ)
                MIJ = NIJEND - NIJBEG + 1
                MIJCEN  = MIJ * NCENC
                MGIJCEN = NGQP * MIJCEN

                CALL  OED__NAI_DERV_3CEN_PCGTO_BLOCK
     +
     +                     ( NBATCH,
     +                       NINT1DX,NINT1DY,NINT1DZ,
     +                       ATOMIC,
     +                       MIJ,NCENC,MIJCEN,
     +                       NIJ,NIJBEG,NIJEND,
     +                       NGQP,NMOM,NGQSCR,MGIJCEN,
     +                       NPGTOA,NPGTOB,
     +                       NXYZA,NXYZB,
     +                       SHELLA,SHELLB,SHELLP,
     +                       XA,YA,ZA,XB,YB,ZB,
     +                       ABX,ABY,ABZ,
     +                       NUCLEI,
     +                       XN,YN,ZN,NCHARGE,
     +                       ICORE (INUCCEN),
     +                       DERAX,DERAY,DERAZ,
     +                       DERBX,DERBY,DERBZ,
     +                       DERCX,DERCY,DERCZ,
     +                       DIFFX,DIFFY,DIFFZ,
     +                       DIFFA,DIFFB,DIFFC,
     +                       IXAEQB,
     +                       ALPHA (LEXPA),ALPHA (LEXPB),
     +                       FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                       ICORE (IPRIMA+NIJBEG-1),
     +                       ICORE (IPRIMB+NIJBEG-1),
     +                       ZCORE (ZNORMA),ZCORE (ZNORMB),
     +                       ZCORE (ZRHOAB),
     +                       ZCORE (ZPX),ZCORE (ZPY),ZCORE (ZPZ),
     +                       ZCORE (ZPAX),ZCORE (ZPAY),ZCORE (ZPAZ),
     +                       ZCORE (ZPINVHF),ZCORE (ZSCALE),
     +                       ZCORE (ZRTS),ZCORE (ZWTS),ZCORE (ZGQSCR),
     +                       ZCORE (ZTVAL),
     +                       ZCORE (ZR1X),ZCORE (ZR1Y),ZCORE (ZR1Z),
     +                       ZCORE (ZR2),
     +                       ZCORE (ZEXP2A),ZCORE (ZEXP2B),
     +                       ZCORE (ZEXP2AB),
     +                       ZCORE (ZINT1DX),
     +                       ZCORE (ZINT1DY),
     +                       ZCORE (ZINT1DZ),
     +
     +                                 ZCORE (ZPBATCH) )
     +
     +
C            WRITE (*,*) ' Finished case III ab nai pcgto derv block '

                CALL  OED__CTR_2INDEX_BLOCK
     +
     +                     ( NPSIZE,NCSIZE,NWSIZE,
     +                       NXYZT,
     +                       MIJ,NCGTOAB,
     +                       NPGTOA,NPGTOB,
     +                       NCGTOA,NCGTOB,
     +                       MXPRIM,MNPRIM,
     +                       CC (LCCA),CC (LCCB),
     +                       CCBEG (LCCSEGA),CCBEG (LCCSEGB),
     +                       CCEND (LCCSEGA),CCEND (LCCSEGB),
     +                       ICORE (IPRIMA+NIJBEG-1),
     +                       ICORE (IPRIMB+NIJBEG-1),
     +                       L1CACHE,TILE,NCTROW,
     +                       EQUALAB,
     +                       SWAPRS,
     +                       REORDER,
     +                       BLOCKED,
     +                       ICORE (IPUSED),
     +                       ICORE (IPSAVE),
     +                       ICORE (IPPAIR),
     +                       ZCORE (ZPBATCH),
     +                       ZCORE (ZWORK),
     +
     +                                 ZCORE (ZCBATCH) )
     +
     +
C                WRITE (*,*) ' Finished case III 2 index ctr block '

 3000        CONTINUE

         END IF
C
C
C             ...the unnormalized cartesian (a|b) contracted nuclear
C                attraction derivative batch is ready. Expand the
C                contraction indices (if necessary):
C
C                   batch (nxyzt,r>=s) --> batch (nxyzt,r,s)
C
C                and reorder the contraction index part (if necessary):
C
C                   batch (nxyzt,r,s) --> batch (nxyzt,1,2)
C
C                The array IXOFF (x) indicates the total # of indices
C                to the left of x without including the nxyzt-part.
C                For the left most IXOFF value it is convenient to
C                set it equal to 1 instead of 0. Note, that the IXOFF
C                array indicates the true # of indices to the left
C                after! the batch has been transposed (see below) and
C                can be used as initial values when moving the
C                ry-components later during the cartesian -> spherical
C                transformation procedure.
C                
C                Efficient application of the cartesian -> spherical
C                transformation requires the batch elements to be
C                ordered as:
C
C                            batch (1,2,nxyzt)
C
C                hence we transpose the batch after the reordering.
C
C                The space partitioning of the flp array for all of
C                these steps will be as follows:
C
C
C                          |  Zone 1  |  Zone 2  |
C
C                in which Zone 1 and 2 are 2 batches of cartesian
C                (a|b) size.
C
C
         IXOFF (1) = 1
         IXOFF (2) = NCGTO1

         NCTR = NCGTO1 * NCGTO2
         NBATCH = NCTR * NXYZT

         IN = ZCBATCH
         OUT = IN + NBATCH

         IF (EQUALAB .AND. NCGTOAB.GT.1) THEN
             CALL  OED__CTR_RS_EXPAND
     +
     +                  ( NXYZT,NCGTOAB,
     +                    NCGTOA,NCGTOB,
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

         REORDER = SWAP12 .NEQV. SWAPRS

         IF (REORDER .AND. NCTR.GT.1) THEN

             IF (SWAPRS) THEN
                 INDEXR = INDEXB
                 INDEXS = INDEXA
                 NCGTOR = NCGTOB
                 NCGTOS = NCGTOA
             ELSE
                 INDEXR = INDEXA
                 INDEXS = INDEXB
                 NCGTOR = NCGTOA
                 NCGTOS = NCGTOB
             END IF

             CALL  OED__CTR_2INDEX_REORDER
     +
     +                  ( NXYZT,NCTR,
     +                    NCGTOR,NCGTOS,
     +                    IXOFF (INDEXR),IXOFF (INDEXS),
     +                    ZCORE (IN),
     +
     +                            ZCORE (OUT) )
     +
     +
C             WRITE (*,*) ' Finished ctr reorder '

             TEMP = IN
             IN = OUT
             OUT = TEMP
         END IF

         IF (NXYZT.GT.1 .AND. NCTR.GT.1) THEN

             CALL  OED__TRANSPOSE_BATCH
     +
     +                  ( NXYZT,NCTR,
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
C             ...enter the cartesian -> spherical transformation or
C                cartesian normalization section. After one shell
C                is transformed it is moved to its final position.
C                This leads to the sequence of events (where ' means
C                spherical or cartesian normalization and [] means the
C                indices are in positional correspondence):
C
C                       batch (ij,a,b) --> batch (ij,a,b')
C                       batch (ij,a,b') --> batch (ij[b'],a)
C                       batch (ij[b'],a) --> batch (ij[b'],a')
C                       batch (ij[b'],a') --> batch (ij[a'b'])
C
C
C                The space partitioning of the flp array will be
C                as follows:
C
C
C                     |  Zone 1  |  Zone 2  |  Zone 3  |
C
C
C                 Zone 1 and 2:  2 batches of cartesian (a|b) size
C                                (set previously)
C
C                       Zone 3:  cart -> spher transformation data
C                                             or
C                                cartesian normalization factors
C
C
C                Determine the memory allocation offsets for the
C                cartesian -> spherical transformations or cartesian
C                normalizations and generate the transformation
C                matrices + associated data for those shells > p-shell.
C                The offsets are as follows (x=A,B):
C
C                    IN = offset for input cartesian (a|b) batch
C                   OUT = offset for output cartesian (a|b) batch
C
C                ZSROTx = offset for x-part transformation matrix
C               ISNROWx = offset for # of non-zero XYZ contribution row
C                         labels for x-part transformation matrix
C                ISROWx = offset for non-zero XYZ contribution row
C                         labels for x-part transformation matrix
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
         ZBASE = MAX (IN,OUT) + NBATCH

         IF (SPHERIC) THEN
             IF (MXSHELL.GT.1) THEN
                 CALL  OED__XYZ_TO_RY_AB
     +
     +                      ( NXYZA,NXYZB,
     +                        NRYA,NRYB,
     +                        SHELLA,SHELLB,
     +                        1,ZBASE,
     +
     +                                 NROWA,NROWB,
     +                                 NROTA,NROTB,
     +                                 ZSROTA,ZSROTB,
     +                                 ISNROWA,ISNROWB,
     +                                 ISROWA,ISROWB,
     +                                 IUSED,ZUSED,
     +                                 ICORE,ZCORE )
     +
     +
C                 WRITE (*,*) ' Finished xyz to ry ab '
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
C
C
C             ...do cart -> spher transformation / cart normalization
C                (if b > 1):
C
C                        (ij,a,b) --> (ij,a,b')
C
C
         IF (SHELLB.GT.1) THEN
             IF (SPHERIC) THEN
                 CALL  OED__SPHERICAL_TRANSFORM
     +
     +                      ( NCTR*NXYZA,
     +                        NROWB,NXYZB,NRYB,
     +                        ICORE (ISNROWB),
     +                        ICORE (ISROWB),
     +                        ZCORE (ZSROTB),
     +                        ZCORE (IN),
     +
     +                                ZCORE (OUT) )
     +
     +
C                 WRITE (*,*) ' Finished sph quart b '
                 TEMP = IN
                 IN = OUT
                 OUT = TEMP
             ELSE
                 CALL  OED__NORMALIZE_CARTESIAN
     +
     +                      ( NCTR*NXYZA,
     +                        NXYZB,
     +                        SHELLB,
     +                        ZCORE (ZCNORM),
     +
     +                                ZCORE (IN) )
     +
     +
C                 WRITE (*,*) ' Finished normalized cart b '
             END IF
         END IF
C
C
C             ...move transformed b shell (if size > 1):
C
C                        (ij,a,b') --> (ij[b'],a)
C
C
         IF (NRYB.GT.1) THEN
             NBATCH = NCTR * NXYZA * NRYB
             NOTMOVE = IXOFF (INDEXB)
             MOVE = NBATCH / (NOTMOVE * NRYB)
             IF (MOVE.GT.1) THEN

                 CALL  OED__MOVE_RY
     +
     +                      ( NBATCH,2,
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
C             ...do cart -> spher / cart normalization (if a > 1):
C
C                        (ij[b'],a) --> (ij[b'],a')
C
C
         IF (SHELLA.GT.1) THEN
             IF (SPHERIC) THEN
                 CALL  OED__SPHERICAL_TRANSFORM
     +
     +                      ( NCTR*NRYB,
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
     +                      ( NCTR*NRYB,
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
C                        (ij[b'],a') --> (ij[a'b'])
C
C
         NBATCH = NCTR * NRYB * NRYA

         IF (NRYA.GT.1) THEN
             NOTMOVE = IXOFF (INDEXA)
             MOVE = NBATCH / (NOTMOVE * NRYA)
             IF (MOVE.GT.1) THEN

                 CALL  OED__MOVE_RY
     +
     +                      ( NBATCH,2,
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
