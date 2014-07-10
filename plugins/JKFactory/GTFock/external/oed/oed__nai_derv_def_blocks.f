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
         SUBROUTINE  OED__NAI_DERV_DEF_BLOCKS
     +
     +               ( ZMAX,
     +                 NPGTOA,NPGTOB,
     +                 SHELLA,SHELLB,SHELLP,
     +                 NIJ,NRS,NCEN,
     +                 NGQP,NGQSCR,
     +                 NXYZT,
     +                 DERAX,DERAY,DERAZ,
     +                 DERBX,DERBY,DERBZ,
     +                 CASEI,CASEII,CASEIII,
     +                 L1CACHE,NCTROW,
     +                 MEMORY,
     +
     +                     NIJBLK,
     +                     NBATCH,
     +                     NPSIZE,NCSIZE,NWSIZE,
     +                     NINT1DX,NINT1DY,NINT1DZ,
     +                     MXPRIM,MNPRIM,
     +                     ZCBATCH,ZPBATCH,ZWORK,
     +                     ZNORMA,ZNORMB,
     +                     ZRHOAB,
     +                     ZPX,ZPY,ZPZ,ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCALE,
     +                     ZRTS,ZWTS,ZGQSCR,ZTVAL,
     +                     ZR1X,ZR1Y,ZR1Z,ZR2,
     +                     ZEXP2A,ZEXP2B,ZEXP2AB,
     +                     ZINT1DX,ZINT1DY,ZINT1DZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_DERV_DEF_BLOCKS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine handles the memory partitions for the
C                generation of the entire contracted derivative (a|b)
C                nuclear attraction batch. It determines the block size
C                partitions (if any) for the ij exponent pairs of the
C                whole primitive derivative [a|b] batch and returns
C                pointers to the flp data sections needed by the (a|b)
C                generation.
C
C                The routine is also equipped with the possibility
C                of returning just the minimum / optimum flp memory
C                needed, without evaluating the (a|b) generation flp
C                pointers. This is useful for establishing just the
C                overall memory size needed for the (a|b) generation.
C                The keyword MEMORY has to be set true for this case
C                and obviously the value of ZMAX passed is irrelevant.
C
C
C                  Input:
C
C                    ZMAX         =  maximum flp memory
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x=A,B
C                    SHELLx       =  the shell type for contraction
C                                    shells x = A,B and P=A+B
C                    NIJ          =  total # of ij primitive index
C                                    pairs corresponding to the csh
C                                    pairs A,B
C                    NRS          =  total # of rs contraction index
C                                    pairs corresponding to the csh
C                                    pairs A,B
C                    NCEN         =  # of nuclear attraction centers
C                                    to be considered here
C                    NGQP         =  # of gaussian quadrature points
C                                    (roots)
C                    NGQSCR       =  size of gaussian quadrature
C                                    scratch space needed to calculate
C                                    all the quadrature roots
C                    NXYZT        =  total # of cartesian monomial
C                                    pairs
C                    DERxp        =  the order of differentiation on
C                                    centers x = A,B with respect to
C                                    the p = x,y,z coordinates
C                    CASEx        =  Type of nuclear attraction
C                                    integrals present: x = I are
C                                    (a|A|b) integrals, x = II are
C                                    (a|B|b) integrals and x = III are
C                                    (a|C|b) integrals with C = not A,B
C                    L1CACHE      =  Size of level 1 cache in units of
C                                    8 Byte
C                    NCTROW       =  minimum # of rows that are
C                                    accepted for blocked contractions
C                    MEMORY       =  if this keyword is true, the
C                                    routine will only determine the
C                                    minimum / optimum flp memory and
C                                    store these values into NPSIZE and
C                                    NCSIZE, respectively (see below)
C
C                  Output:
C
C                    NIJBLK       =  contains the block size for the
C                                    ij primitive indices in order to
C                                    perform efficient contractions.
C                    NBATCH       =  size of the integral batch array
C                                    during primitive integral block
C                                    generation. This array is also
C                                    needed for 1D integral derivative
C                                    evaluation, hence its size might
C                                    be larger than the final primitive
C                                    integral block
C                    NxSIZE       =  size of the primitive integral
C                                    block (x=P), contracted integrals
C                                    (x=C) and working (x=W) arrays.
C                                    If MEMORY is true, NPSIZE and
C                                    NCSIZE will contain respectively
C                                    the minimum and optimum flp memory
C                    NINT1Dx      =  space needed for each of the 1D
C                                    x = X,Y,Z integral arrays (they
C                                    might be different due to different
C                                    orders of differentiation for
C                                    each cartesian component)
C                    MXPRIM       =  the maximum # of primitives
C                                    between all i and j primitives,
C                                    i.e. = max (i,j)
C                    MNPRIM       =  the minimum # of primitives
C                                    between all i and j primitives,
C                                    i.e. = min (i,j)
C                    Z.....       =  pointers for space partition of
C                                    the flp array (see below)
C
C
C                The space partitioning of the flp array will be
C                as follows:
C
C
C                      |  Zone 1  |  Zone 2  |  Zone 3  |
C
C
C                   Zone 1 : final (a|b) derivative batch
C                   Zone 2 : complete set of A,B norms and complete
C                            set (after screening!) of exponential
C                            prefactors and nuclear attraction center
C                            information (coordinates + nuclear charges)
C                   Zone 3 : will be used at two different stages:
C                             i) during derivated primitive [a|b]
C                                block generation
C                            ii) during the contraction stage
C
C
C                While the first two Zones are straightforward in
C                their definition, Zone 3 requires some more detailed
C                explanation:
C
C                 Stage i) : During this stage the primitive blocks
C                            have to be assembled from the derivated
C                            1D integrals, which in turn are generated
C                            in single derivative steps. The derivation
C                            of the 1D integrals needs an extra work
C                            space for accumulation, which will be
C                            provided at the same memory location where
C                            the final primitive block will be. This
C                            can be done, because the primitive block
C                            can only be assembled after! all the
C                            derivated 1D integrals are ready.
C
C                Stage ii) : This stage needs the primitive blocks
C                            plus some workspace for performing the
C                            contraction step.
C
C
C                Memory allocation offsets for the primitive [a|b]
C                derivative batches generation + contraction to (a|b):
C
C
C                  --- Zone 1 ---
C
C                   ZCBATCH = offset for contracted (a|b) batch
C
C                  --- Zone 2 ---
C
C                   ZNORMA = offset for all A norms
C                   ZNORMB = offset for all B norms
C                   ZRHOAB = offset for all AB exponential prefactors
C
C                  --- Zone 3: for Stage i) only ---
C
C                   ZPBATCH = offset for derivated 1D integrals
C                             accumulation work space and for
C                             (blocked) derivated [a|b] batch
C
C                   ZPX = offset for (blocked) PX values
C                   ZPY = offset for (blocked) PY values
C                   ZPZ = offset for (blocked) PZ values
C                   ZPAX = offset for (blocked) PAX values
C                   ZPAY = offset for (blocked) PAY values
C                   ZPAZ = offset for (blocked) PAZ values
C                   ZPINVHF = offset for (blocked) 1/2P values
C                   ZSCALE = offset for (blocked) scale values
C
C                   ZRTS = offset for (blocked) quadrature roots 
C                   ZWTS = offset for (blocked) quadrature weights
C                   ZGQSCR = offset for (blocked) quadrature scratch
C                   ZTVAL = offset for (blocked) T exponents
C
C                   ZR1X = offset for (blocked) VRR coeff R1X
C                   ZR1Y = offset for (blocked) VRR coeff R1Y
C                   ZR1Z = offset for (blocked) VRR coeff R1Z
C                   ZR2 = offset for (blocked) VRR coeff R2
C
C                   ZEXP2A = offset for 2 x A-exponent values
C                   ZEXP2B = offset for 2 x B-exponent values
C                   ZEXP2AB = offset for 2 x (A+B)-exponent values
C
C                   ZINT1DX = offset for (blocked) derivative 1DX
C                             integrals accumulation
C                   ZINT1DY = offset for (blocked) derivative 1DY
C                             integrals accumulation
C                   ZINT1DZ = offset for (blocked) derivative 1DZ
C                             integrals accumulation
C
C                  --- Zone 3: for Stage ii) only ---
C
C                   ZPBATCH = offset for (blocked) derivated [a|b]
C                             batch
C
C                   ZWORK = offset for contraction working array
C
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     CASEI,CASEII,CASEIII
         LOGICAL     MEMORY

         INTEGER     DERAX,DERAY,DERAZ
         INTEGER     DERBX,DERBY,DERBZ
         INTEGER     IJDIV
         INTEGER     ISHELLX,ISHELLY,ISHELLZ
         INTEGER     L1CACHE,NCTROW
         INTEGER     MIJ,MIJBIG,NCEN,MIJCEN,MGIJCEN
         INTEGER     MXPRIM,MNPRIM
         INTEGER     MXSTEP
         INTEGER     N
         INTEGER     NBATCH
         INTEGER     NGQP,NGQSCR
         INTEGER     NIJ,NIJBLK
         INTEGER     NINT1DX,NINT1DY,NINT1DZ
         INTEGER     NPGTOA,NPGTOB
         INTEGER     NPSIZE,NCSIZE,NWSIZE
         INTEGER     NRS
         INTEGER     NXYZT
         INTEGER     PFACT
         INTEGER     POW2IJ
         INTEGER     PSHELL
         INTEGER     SHELLA,SHELLB,SHELLP
         INTEGER     SHELLX,SHELLY,SHELLZ
         INTEGER     SIZEAX,SIZEAY,SIZEAZ
         INTEGER     SIZEBX,SIZEBY,SIZEBZ
         INTEGER     SIZEPX,SIZEPY,SIZEPZ
         INTEGER     WCOL
         INTEGER     ZMAX
         INTEGER     ZONE1,ZONE2,ZONE3,ZONE3B,ZONE3C
         INTEGER     ZONEMIN,ZONEOPT
         INTEGER     ZCBATCH,ZPBATCH,ZWORK,
     +               ZNORMA,ZNORMB,ZRHOAB,
     +               ZPX,ZPY,ZPZ,ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCALE,
     +               ZRTS,ZWTS,ZGQSCR,ZTVAL,
     +               ZR1X,ZR1Y,ZR1Z,ZR2,
     +               ZEXP2A,ZEXP2B,ZEXP2AB,
     +               ZINT1DX,ZINT1DY,ZINT1DZ

         DOUBLE PRECISION   DLOG2

         PARAMETER  (DLOG2  =  0.6931471805599D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...determine optimum size of [a|b] blocks:
C
C                The next few lines of code are among the most
C                important for performance of the entire overlap
C                integral package. Any badly chosen sizes for the
C                primitive blocks has a severe effect on the integral
C                evaluation timing. Two opposite effects have
C                to be considered:
C
C                1) Reuse level 1 cache data as much as possible,
C                   i.e. choose the primitive blocks in such a
C                   way that most of the heavy computational work
C                   is being done reusing data in each cache line.
C
C                2) Choose the size of the primitive blocks as
C                   large as possible observing point 1), to reduce
C                   subroutine call overheads as much as possible.
C
C                After some experimenting, the following procedure
C                was found to work best:
C
C                  i) The optimum [a|b] block size is assumed to
C                     be directly proportional to the level 1 cache
C                     size.
C
C                 ii) The proportionality factor PFACT has to be
C                     determined by 'experiment', running a real
C                     molecular case with fully contracted basis sets.
C
C
         PFACT = 4

         MIJ = PFACT * L1CACHE / NXYZT
         MIJ = MAX (1,MIJ)
         MIJ = MIN (MIJ,NIJ)
C
C
C             ...the optimum block size MIJ has been determined.
C                We have to see now how it fits with the rest of
C                the needed data into the maximum memory given.
C                If necessary subdivide the optimum block size
C                into convenient subblocks that fit into memory.
C
C                The subdivision of the MIJ block is done in such
C                a way that successive powers of 2 are checked as
C                divisors of MIJ. The maximum # of division steps
C                MXSTEP can thus be calculated beforehand by knowing
C                the fact that the smallest subdivided MIJ block
C                is of size 1. The succesive divisors are checked
C                in the following order:
C
C
C                            step  |  div of MIJ
C                           ----------------------
C                              1   |       1
C                              2   |       2
C                              3   |       4
C                              4   |       8
C                              5   |      16
C                             ...         ...
C
C                The routine stops, if even a MIJ block of minimum
C                size 1 cannot be accomodated into the memory.
C
C
         SIZEPX  = (SHELLP + DERAX + DERBX + 1)
         SIZEPY  = (SHELLP + DERAY + DERBY + 1)
         SIZEPZ  = (SHELLP + DERAZ + DERBZ + 1)
         SHELLX  = 0
         SHELLY  = 0
         SHELLZ  = 0
         ISHELLX = 0
         ISHELLY = 0
         ISHELLZ = 0

         IF (CASEI) THEN
             SIZEAX  = (SHELLA + 1)
             SIZEAY  = (SHELLA + 1)
             SIZEAZ  = (SHELLA + 1)
             SIZEBX  = (SHELLB + DERAX + DERBX + 1)
             SIZEBY  = (SHELLB + DERAY + DERBY + 1)
             SIZEBZ  = (SHELLB + DERAZ + DERBZ + 1)
             SHELLX  = MAX (SHELLX  , SIZEPX * MIN (SIZEAX,SIZEBX))
             SHELLY  = MAX (SHELLY  , SIZEPY * MIN (SIZEAY,SIZEBY))
             SHELLZ  = MAX (SHELLZ  , SIZEPZ * MIN (SIZEAZ,SIZEBZ))
             ISHELLX = MAX (ISHELLX , SIZEAX * SIZEBX)
             ISHELLY = MAX (ISHELLY , SIZEAY * SIZEBY)
             ISHELLZ = MAX (ISHELLZ , SIZEAZ * SIZEBZ)
         END IF

         IF (CASEII) THEN
             SIZEAX  = (SHELLA + DERAX + DERBX + 1)
             SIZEAY  = (SHELLA + DERAY + DERBY + 1)
             SIZEAZ  = (SHELLA + DERAZ + DERBZ + 1)
             SIZEBX  = (SHELLB + 1)
             SIZEBY  = (SHELLB + 1)
             SIZEBZ  = (SHELLB + 1)
             SHELLX  = MAX (SHELLX  , SIZEPX * MIN (SIZEAX,SIZEBX))
             SHELLY  = MAX (SHELLY  , SIZEPY * MIN (SIZEAY,SIZEBY))
             SHELLZ  = MAX (SHELLZ  , SIZEPZ * MIN (SIZEAZ,SIZEBZ))
             ISHELLX = MAX (ISHELLX , SIZEAX * SIZEBX)
             ISHELLY = MAX (ISHELLY , SIZEAY * SIZEBY)
             ISHELLZ = MAX (ISHELLZ , SIZEAZ * SIZEBZ)
         END IF

         IF (CASEIII) THEN
             SIZEAX  = (SHELLA + DERAX + 1)
             SIZEAY  = (SHELLA + DERAY + 1)
             SIZEAZ  = (SHELLA + DERAZ + 1)
             SIZEBX  = (SHELLB + DERBX + 1)
             SIZEBY  = (SHELLB + DERBY + 1)
             SIZEBZ  = (SHELLB + DERBZ + 1)
             SHELLX  = MAX (SHELLX  , SIZEPX * MIN (SIZEAX,SIZEBX))
             SHELLY  = MAX (SHELLY  , SIZEPY * MIN (SIZEAY,SIZEBY))
             SHELLZ  = MAX (SHELLZ  , SIZEPZ * MIN (SIZEAZ,SIZEBZ))
             ISHELLX = MAX (ISHELLX , SIZEAX * SIZEBX)
             ISHELLY = MAX (ISHELLY , SIZEAY * SIZEBY)
             ISHELLZ = MAX (ISHELLZ , SIZEAZ * SIZEBZ)
         END IF

         PSHELL = MAX (SHELLX,SHELLY,SHELLZ)

         MXPRIM = MAX (NPGTOA,NPGTOB)
         MNPRIM = MIN (NPGTOA,NPGTOB)

         NCSIZE = NXYZT * NRS

         ZONE2 = NPGTOA + NPGTOB + NIJ
C
C
C             ...if the MEMORY keyword is activated, the routine
C                will only determine the optimum and minimum flp
C                memory (in that order), place them respectively
C                into the NCSIZE and NPSIZE variables and exit.
C
C
         IF (MEMORY) THEN

             MIJCEN  = MIJ * NCEN
             MGIJCEN = NGQP * MIJCEN
             NINT1DX = MGIJCEN * ISHELLX
             NINT1DY = MGIJCEN * ISHELLY
             NINT1DZ = MGIJCEN * ISHELLZ
             NBATCH  = MAX (MGIJCEN*PSHELL,NXYZT*MIJ)
             NPSIZE  = NXYZT * MIJ
             WCOL    = MNPRIM
             NWSIZE  = NCTROW * WCOL
             NPSIZE  = MAX (NPSIZE,NWSIZE)
             NWSIZE  = NPSIZE
             ZONE1   = NCSIZE
             ZONE3B  = NBATCH + NINT1DX + NINT1DY + NINT1DZ
     +                 + NGQSCR + MIJCEN + 8*MIJ + 9*MGIJCEN
             ZONE3C  = NPSIZE + NWSIZE
             ZONE3   = MAX (ZONE3B,ZONE3C)

             ZONEOPT = ZONE1 + ZONE2 + ZONE3

             MIJ     = 1
             MIJCEN  = MIJ * NCEN
             MGIJCEN = NGQP * MIJCEN
             NINT1DX = MGIJCEN * ISHELLX
             NINT1DY = MGIJCEN * ISHELLY
             NINT1DZ = MGIJCEN * ISHELLZ
             NBATCH  = MAX (MGIJCEN*PSHELL,NXYZT*MIJ)
             NPSIZE  = NXYZT * MIJ
             WCOL    = MNPRIM
             NWSIZE  = NCTROW * WCOL
             NPSIZE  = MAX (NPSIZE,NWSIZE)
             NWSIZE  = NPSIZE
             ZONE1   = NCSIZE
             ZONE3B  = NBATCH + NINT1DX + NINT1DY + NINT1DZ
     +                 + NGQSCR + MIJCEN + 8*MIJ + 9*MGIJCEN
             ZONE3C  = NPSIZE + NWSIZE
             ZONE3   = MAX (ZONE3B,ZONE3C)

             ZONEMIN = ZONE1 + ZONE2 + ZONE3

             NPSIZE  = ZONEMIN
             NCSIZE  = ZONEOPT

             RETURN

         ELSE
C
C
C             ...the actual fitting into the maximum memory given.
C
C
             POW2IJ = DLOG ( DFLOAT (MIJ) ) / DLOG2
             MXSTEP = POW2IJ + 1

             IJDIV = 1
             MIJBIG = MIJ

             DO 100 N = 1,MXSTEP

                MIJ = MAX (1, MIJBIG / IJDIV )

                MIJCEN  = MIJ * NCEN
                MGIJCEN = NGQP * MIJCEN
                NINT1DX = MGIJCEN * ISHELLX
                NINT1DY = MGIJCEN * ISHELLY
                NINT1DZ = MGIJCEN * ISHELLZ
                NBATCH  = MAX (MGIJCEN*PSHELL,NXYZT*MIJ)
                NPSIZE  = NXYZT * MIJ
                WCOL    = MNPRIM
                NWSIZE  = NCTROW * WCOL
                NPSIZE  = MAX (NPSIZE,NWSIZE)
                NWSIZE  = NPSIZE
                ZONE1   = NCSIZE
                ZONE3B  = NBATCH + NINT1DX + NINT1DY + NINT1DZ
     +                    + NGQSCR + MIJCEN + 8*MIJ + 9*MGIJCEN
                ZONE3C  = NPSIZE + NWSIZE
                ZONE3   = MAX (ZONE3B,ZONE3C)

                IF (ZONE1+ZONE2+ZONE3.LE.ZMAX) THEN

                    NIJBLK = MIJ
                    NWSIZE = ZMAX - ZONE1 - ZONE2 - NPSIZE
C
C
C             ...generate the memory allocation pointers.
C
C
                    ZCBATCH = 1

                    ZNORMA = ZCBATCH + NCSIZE
                    ZNORMB = ZNORMA + NPGTOA
                    ZRHOAB = ZNORMB + NPGTOB

                    ZPBATCH = ZRHOAB + NIJ

                    ZPX = ZPBATCH + NBATCH
                    ZPY = ZPX + MIJ
                    ZPZ = ZPY + MIJ
                    ZPAX = ZPZ + MIJ
                    ZPAY = ZPAX + MIJ
                    ZPAZ = ZPAY + MIJ
                    ZPINVHF = ZPAZ + MIJ
                    ZSCALE = ZPINVHF + MIJ

                    ZRTS = ZSCALE + MGIJCEN
                    ZWTS = ZRTS + MGIJCEN
                    ZGQSCR = ZWTS + MGIJCEN
                    ZTVAL = ZGQSCR + NGQSCR

                    ZR1X = ZTVAL + MIJCEN
                    ZR1Y = ZR1X + MGIJCEN
                    ZR1Z = ZR1Y + MGIJCEN
                    ZR2 = ZR1Z + MGIJCEN

                    ZEXP2A = ZR2 + MGIJCEN
                    ZEXP2B = ZEXP2A + MGIJCEN
                    ZEXP2AB = ZEXP2B + MGIJCEN

                    ZINT1DX = ZEXP2AB + MIJ
                    ZINT1DY = ZINT1DX + NINT1DX
                    ZINT1DZ = ZINT1DY + NINT1DY

                    ZWORK = ZPBATCH + NPSIZE

                    RETURN
                END IF

                IJDIV = IJDIV + IJDIV

  100        CONTINUE

             WRITE (*,*) ' Memory allocation failed for derv (a|b) ! '
             WRITE (*,*) ' NIJ,MIJ = ',NIJ,MIJ
             WRITE (*,*) ' Increase flp memory! '
             WRITE (*,*) ' (oed__nai_derv_def_blocks) '

             STOP

         END IF
C
C
C             ...ready!
C
C
         END
