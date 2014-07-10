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
         SUBROUTINE  OED__XYZ_E0_DEF_BLOCKS
     +
     +                    ( ZMAX,
     +                      MOMENTX,MOMENTY,MOMENTZ,
     +                      NPGTOA,NPGTOB,
     +                      SHELLP,
     +                      NIJ,NRS,
     +                      NXYZT,
     +                      L1CACHE,NCTROW,
     +                      MEMORY,
     +
     +                              NIJBLK,
     +                              NPSIZE,NCSIZE,NWSIZE,
     +                              NINT1D,
     +                              MXPRIM,MNPRIM,
     +                              ZCBATCH,ZPBATCH,ZWORK,
     +                              ZNORMA,ZNORMB,
     +                              ZRHOAB,ZXP,ZYP,ZZP,
     +                              ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCALE,
     +                              ZINT1DX,ZINT1DY,ZINT1DZ,
     +                              ZSCRMX,ZSCRMY,ZSCRMZ,
     +                              ZXMINUS,ZYMINUS,ZZMINUS )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__XYZ_E0_DEF_BLOCKS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine handles the memory partitions for the
C                generation of the entire contracted (e|0) moment
C                batch. It determines the block size partitions (if any)
C                for the ij exponent paris of the whole primitive
C                [e|0] batch and returns pointers to the flp data
C                sections needed by the (e|0) generation.
C
C                The routine is also equipped with the possibility
C                of returning just the minimum / optimum flp memory
C                needed, without evaluating the (e|0) generation flp
C                pointers. This is useful for establishing just the
C                overall memory size needed for the (e|0) generation.
C                The keyword MEMORY has to be set true for this case
C                and obviously the value of ZMAX passed is irrelevant.
C
C
C                  Input:
C
C                    ZMAX         =  maximum flp memory
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x=A,B
C                    SHELLP       =  the shell sum A+B
C                    NIJ          =  total # of ij primitive index
C                                    pairs corresponding to the csh
C                                    pairs A,B
C                    NRS          =  total # of rs contraction index
C                                    pairs corresponding to the csh
C                                    pairs A,B
C                    NXYZT        =  total # of cartesian monomial
C                                    pairs
C                    L1CACHE      =  Size of level 1 cache in units of
C                                    8 Byte
C                    NCTROW       =  minimum # of rows that are
C                                    accepted for blocked contractions
C                    MOMENTx      =  if MOMENTx (x = X,Y,Z) is greater
C                                    than zero, then set aside memory
C                                    in SCRMTx
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
C                    NxSIZE       =  size of the primitive integral
C                                    block (x=P), contracted integrals
C                                    (x=C) and working (x=W) arrays.
C                                    If MEMORY is true, NPSIZE and
C                                    NCSIZE will contain respectively
C                                    the minimum and optimum flp memory
C                    NINT1D       =  space needed for the 1D X,Y,Z
C                                    integral arrays
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
C                 |  Zone 1  |  Zone 2  |  Zone 3  |  Zone 4  |
C
C
C                   Zone 1: final (e|0) batch
C                   Zone 2: block [e|0]
C                   Zone 3 : complete set of A,B norms and complete
C                            set of AB exponential prefactors
C                   Zone 4 : i) scratch for block [e|0] generation
C                           ii) scratch for (e|0) generation
C
C
C                Memory allocation offsets for the primitive [e|0]
C                batches generation + contraction to (e|0):
C
C
C                  --- Zone 1 ---
C
C                   ZCBATCH = offset for contracted (e|0) batch
C
C                  --- Zone 2 ---
C
C                   ZPBATCH = offset for (blocked) [e|0] batch
C
C                  --- Zone 3 ---
C
C                   ZNORMA = offset for all A norms
C                   ZNORMB = offset for all B norms
C                   ZRHOAB = offset for all AB exponential prefactors
C
C                  --- Zone 4: for block [e|0] generation only ---
C
C                   ZPAX = offset for (blocked) PAX values
C                   ZPAY = offset for (blocked) PAY values
C                   ZPAZ = offset for (blocked) PAZ values
C                   ZPINVHF = offset for (blocked) 1/2P values
C                   ZSCALE = offset for (blocked) scale values
C
C                   ZINT1DX = offset for (blocked) 1DX integrals
C                   ZINT1DY = offset for (blocked) 1DY integrals
C                   ZINT1DZ = offset for (blocked) 1DZ integrals
C
C                   ZXP = offset for (blocked) XP values
C                   ZYP = offset for (blocked) YP values
C                   ZZP = offset for (blocked) ZP values
C
C                   ZXMINUS = offset for (blocked) previous integrals
C                   ZYMINUS = offset for (blocked) previous integrals
C                   ZZMINUS = offset for (blocked) previous integrals
C
C                   ZSCRMX = offset for (blocked) 1DX scratch integrals
C                   ZSCRMY = offset for (blocked) 1DY scratch integrals
C                   ZSCRMZ = offset for (blocked) 1DZ scratch integrals
C
C                  --- Zone 4: for contraction only ---
C
C                   ZWORK = offset for contraction working array
C
C
C
C  AUTHOR      : Norbert Flocke
C                  - Wrote original OED package
C
C  MODIFIED    : Thomas Watson Jr.                   p  q  r
C                  - Modified OED package to handle X, Y, Z integrals
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     MEMORY
         LOGICAL     CASEX,CASEY,CASEZ
         LOGICAL     CASE1,CASE2,CASE3

         INTEGER     IJDIV
         INTEGER     L1CACHE,NCTROW
         INTEGER     MOMENTX,MOMENTY,MOMENTZ
         INTEGER     N
         INTEGER     MIJ,MIJBIG
         INTEGER     MXPRIM,MNPRIM
         INTEGER     MXSTEP
         INTEGER     NIJ,NIJBLK
         INTEGER     NINT1D
         INTEGER     NPGTOA,NPGTOB
         INTEGER     NPSIZE,NCSIZE,NWSIZE
         INTEGER     NRS
         INTEGER     NXYZT
         INTEGER     PFACT
         INTEGER     POW2IJ
         INTEGER     SHELLP
         INTEGER     WCOL
         INTEGER     ZMAX
         INTEGER     ZONE12,ZONE3,ZONE4,ZONE4B,ZONE4C
         INTEGER     ZONEMIN,ZONEOPT
         INTEGER     ZCBATCH,ZPBATCH,ZWORK,
     +               ZNORMA,ZNORMB,ZRHOAB,
     +               ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCALE,
     +               ZINT1DX,ZINT1DY,ZINT1DZ,
     +               ZXP,ZYP,ZZP,ZSCRMX,ZSCRMY,ZSCRMZ,
     +               ZXMINUS,ZYMINUS,ZZMINUS,ZEND

         DOUBLE PRECISION   DLOG2

         PARAMETER  (DLOG2  =  0.6931471805599D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...determine optimum size of [e|0] blocks:
C
C                The next few lines of code are among the most
C                important for performance of the entire moment
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
C                  i) The optimum [e|0] block size is assumed to
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
C                If necessary subdivide the optimum MIJ block size
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
C             ...Watson
C                 Determine the cases to see whether memory is set
C                 aside for scratch matrices for X, Y, Z integrals.
C
C
         CASEX = (MOMENTX .GT. 0)
         CASEY = (MOMENTY .GT. 0)
         CASEZ = (MOMENTZ .GT. 0)

         CASE1 = (CASEX .NEQV. CASEY .NEQV. CASEZ)
         CASE2 = (CASEX .AND. CASEY) .NEQV. (CASEX .AND.CASEZ) .NEQV.
     +           (CASEY .AND. CASEZ)
         CASE3 = (CASEX .AND. CASEY .AND. CASEZ)


         MXPRIM = MAX (NPGTOA,NPGTOB)
         MNPRIM = MIN (NPGTOA,NPGTOB)

         NCSIZE  = NXYZT * NRS

         ZONE3 = NPGTOA + NPGTOB + NIJ
C
C
C             ...if the MEMORY keyword is activated, the routine
C                will only determine the optimum and minimum flp
C                memory (in that order), place them respectively
C                into the NCSIZE and NPSIZE variables and exit.
C
C
         IF (MEMORY) THEN

             NINT1D  = MIJ * (SHELLP+1)
             NPSIZE  = NXYZT * MIJ
             WCOL    = MNPRIM
             NWSIZE  = NCTROW * WCOL
             NPSIZE  = MAX (NPSIZE,NWSIZE)
             NWSIZE  = NPSIZE
             ZONE12  = NPSIZE + NCSIZE

             IF (CASE1)  ZONE4B = 5 * MIJ + 4 * NINT1D + 4 * NIJ
             IF (CASE2)  ZONE4B = 5 * MIJ + 5 * NINT1D + 5 * NIJ
             IF (CASE3)  ZONE4B = 5 * MIJ + 6 * NINT1D + 6 * NIJ

             ZONE4C  = NWSIZE
             ZONE4   = MAX (ZONE4B,ZONE4C)

             ZONEOPT = ZONE12 + ZONE3 + ZONE4

             MIJ     = 1
             NINT1D  = MIJ * (SHELLP+1)
             NPSIZE  = NXYZT * MIJ
             WCOL    = MNPRIM
             NWSIZE  = NCTROW * WCOL
             NPSIZE  = MAX (NPSIZE,NWSIZE)
             NWSIZE  = NPSIZE
             ZONE12  = NPSIZE + NCSIZE

             IF (CASE1)  ZONE4B = 5 * MIJ + 4 * NINT1D + 4 * NIJ
             IF (CASE2)  ZONE4B = 5 * MIJ + 5 * NINT1D + 5 * NIJ
             IF (CASE3)  ZONE4B = 5 * MIJ + 6 * NINT1D + 6 * NIJ

             ZONE4C  = NWSIZE
             ZONE4   = MAX (ZONE4B,ZONE4C)

             ZONEMIN = ZONE12 + ZONE3 + ZONE4

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

                NINT1D  = MIJ * (SHELLP+1)
                NPSIZE  = NXYZT * MIJ
                WCOL    = MNPRIM
                NWSIZE  = NCTROW * WCOL
                NPSIZE  = MAX (NPSIZE,NWSIZE)
                NWSIZE  = NPSIZE
                ZONE12  = NPSIZE + NCSIZE
C
C
C             ...Watson
C                 ZONE4B changes depending on X,Y,or Z integrals
C
C
                IF (CASE1)  ZONE4B = 5 * MIJ + 4 * NINT1D + 4 * NIJ
                IF (CASE2)  ZONE4B = 5 * MIJ + 5 * NINT1D + 5 * NIJ
                IF (CASE3)  ZONE4B = 5 * MIJ + 6 * NINT1D + 6 * NIJ

                ZONE4C  = NWSIZE
                ZONE4   = MAX (ZONE4B,ZONE4C)

                IF (ZONE12+ZONE3+ZONE4.LE.ZMAX) THEN

                    NIJBLK = MIJ
                    NWSIZE = ZMAX - ZONE12 - ZONE3
C
C
C             ...generate the memory allocation pointers.
C
C
                    ZCBATCH = 1
                    ZPBATCH = NCSIZE + 1

                    ZNORMA = ZPBATCH + NPSIZE
                    ZNORMB = ZNORMA + NPGTOA
                    ZRHOAB = ZNORMB + NPGTOB

                    ZPAX = ZRHOAB + NIJ
                    ZPAY = ZPAX + MIJ
                    ZPAZ = ZPAY + MIJ
                    ZPINVHF = ZPAZ + MIJ
                    ZSCALE = ZPINVHF + MIJ

                    ZINT1DX = ZSCALE + MIJ
                    ZINT1DY = ZINT1DX + NINT1D
                    ZINT1DZ = ZINT1DY + NINT1D

                    ZXP = ZINT1DZ + NINT1D
                    ZYP = ZXP + NIJ
                    ZZP = ZYP + NIJ
C
C
C             ...Watson
C                Begin allocating memory for possible X,Y,Z
C                polynomial integrals.
C
C
                    ZSCRMX = ZZP + NIJ
                    IF (MOMENTX .GT. 0) THEN
                        ZXMINUS = ZSCRMX  + NINT1D
                        ZSCRMY  = ZXMINUS + NIJ
                    ELSE
                        ZXMINUS = ZSCRMX
                        ZSCRMY  = ZXMINUS
                    END IF

                    IF (MOMENTY .GT. 0) THEN
                        ZYMINUS = ZSCRMY  + NINT1D
                        ZSCRMZ  = ZYMINUS + NIJ
                    ELSE
                        ZYMINUS = ZSCRMY
                        ZSCRMZ  = ZYMINUS
                    END IF

                    IF (MOMENTZ .GT. 0) THEN
                        ZZMINUS = ZSCRMZ +  NINT1D
                        ZEND    = ZZMINUS + NIJ
                    ELSE
                        ZZMINUS = ZSCRMZ
                        ZEND    = ZZMINUS
                    END IF

                    ZWORK = ZPAX

                    RETURN
                END IF

                IJDIV = IJDIV + IJDIV

  100        CONTINUE

             WRITE (*,*) ' Memory allocation failed for (e|0) ! '
             WRITE (*,*) ' NIJ,MIJ = ',NIJ,MIJ
             WRITE (*,*) ' Increase flp memory! '
             WRITE (*,*) ' (oed__ovl_e0_def_blocks) '

             STOP

         END IF
C
C
C             ...ready!
C
C
         END
