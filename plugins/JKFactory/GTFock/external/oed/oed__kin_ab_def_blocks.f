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
         SUBROUTINE  OED__KIN_AB_DEF_BLOCKS
     +
     +                    ( ZMAX,
     +                      NPGTOA,NPGTOB,
     +                      SHELLA,SHELLB,SHELLP,
     +                      NIJ,NRS,
     +                      NXYZT,
     +                      L1CACHE,NCTROW,
     +                      MEMORY,
     +
     +                              NIJBLK,
     +                              NPSIZE,NCSIZE,NWSIZE,
     +                              NKIN1D,NOVL1D,
     +                              MXPRIM,MNPRIM,
     +                              ZCBATCH,ZPBATCH,ZWORK,
     +                              ZNORMA,ZNORMB,
     +                              ZRHOAB,
     +                              ZEA,ZEB,ZE2AB,
     +                              ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCALE,
     +                              ZKIN1DX,ZKIN1DY,ZKIN1DZ,
     +                              ZOVL1DX,ZOVL1DY,ZOVL1DZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__KIN_AB_DEF_BLOCKS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine handles the memory partitions for the
C                generation of the entire contracted (a|b) kinetic
C                batch. It determines the block size partitions (if any)
C                for the ij exponent paris of the whole primitive
C                [a|b] batch and returns pointers to the flp data
C                sections needed by the (a|b) generation.
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
C                    SHELLx       =  the shell type for csh x=A,B
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
C                    NKIN1D       =  space needed for the kinetic
C                                    1D X,Y,Z integral arrays
C                    NOVL1D       =  space needed for the overlap
C                                    1D X,Y,Z integral arrays
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
C                   Zone 1: final (a|b) batch
C                   Zone 2: block [a|b]
C                   Zone 3 : complete set of A,B norms and complete
C                            set of AB exponential prefactors 
C                   Zone 4 : i) scratch for block [a|b] generation
C                           ii) scratch for (a|b) generation
C
C
C                Memory allocation offsets for the primitive [a|b]
C                batches generation + contraction to (a|b):
C
C
C                  --- Zone 1 ---
C
C                   ZCBATCH = offset for contracted (a|b) batch
C
C                  --- Zone 2 ---
C
C                   ZPBATCH = offset for (blocked) [a|b] batch
C
C                  --- Zone 3 ---
C
C                   ZNORMA = offset for all A norms
C                   ZNORMB = offset for all B norms
C                   ZRHOAB = offset for all AB exponential prefactors
C
C                  --- Zone 4: for block [a|b] generation only ---
C
C                   ZEA = offset for (blocked) A exponent values
C                   ZEB = offset for (blocked) B exponent values
C                   ZE2AB = offset for (blocked) 2*(A*B exponent) values
C                   ZPAX = offset for (blocked) PAX values
C                   ZPAY = offset for (blocked) PAY values
C                   ZPAZ = offset for (blocked) PAZ values
C                   ZPINVHF = offset for (blocked) 1/2P values
C                   ZSCALE = offset for (blocked) scale values
C
C                   ZKIN1DX = offset for (blocked) kinetic 1DX integrals
C                   ZKIN1DY = offset for (blocked) kinetic 1DY integrals
C                   ZKIN1DZ = offset for (blocked) kinetic 1DZ integrals
C                   ZOVL1DX = offset for (blocked) overlap 1DX integrals
C                   ZOVL1DY = offset for (blocked) overlap 1DY integrals
C                   ZOVL1DZ = offset for (blocked) overlap 1DZ integrals
C
C                  --- Zone 4: for contraction only ---
C
C                   ZWORK = offset for contraction working array
C
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

         LOGICAL     MEMORY

         INTEGER     IJDIV
         INTEGER     L1CACHE,NCTROW
         INTEGER     N
         INTEGER     MIJ,MIJBIG
         INTEGER     MXPRIM,MNPRIM
         INTEGER     MXSTEP
         INTEGER     NIJ,NIJBLK
         INTEGER     NKIN1D,NOVL1D
         INTEGER     NPGTOA,NPGTOB
         INTEGER     NPSIZE,NCSIZE,NWSIZE
         INTEGER     NRS
         INTEGER     NXYZT
         INTEGER     PFACT
         INTEGER     POW2IJ
         INTEGER     SHELLA,SHELLB,SHELLP
         INTEGER     SKIN1D,SOVL1D
         INTEGER     WCOL
         INTEGER     ZMAX
         INTEGER     ZONE12,ZONE3,ZONE4,ZONE4B,ZONE4C
         INTEGER     ZONEMIN,ZONEOPT
         INTEGER     ZCBATCH,ZPBATCH,ZWORK,
     +               ZNORMA,ZNORMB,ZRHOAB,
     +               ZEA,ZEB,ZE2AB,ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCALE,
     +               ZKIN1DX,ZKIN1DY,ZKIN1DZ,ZOVL1DX,ZOVL1DY,ZOVL1DZ

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
C                important for performance of the entire kinetic
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
         MXPRIM = MAX (NPGTOA,NPGTOB)
         MNPRIM = MIN (NPGTOA,NPGTOB)

         NCSIZE  = NXYZT * NRS

         ZONE3 = NPGTOA + NPGTOB + NIJ

         SKIN1D = (SHELLA + 1) * (SHELLB + 1)
         SOVL1D = (SHELLP + 3) * (SHELLB + 2)
C
C
C             ...if the MEMORY keyword is activated, the routine
C                will only determine the optimum and minimum flp
C                memory (in that order), place them respectively
C                into the NCSIZE and NPSIZE variables and exit.
C
C
         IF (MEMORY) THEN

             NKIN1D  = MIJ * SKIN1D
             NOVL1D  = MIJ * SOVL1D
             NPSIZE  = NXYZT * MIJ
             WCOL    = MNPRIM
             NWSIZE  = NCTROW * WCOL
             NPSIZE  = MAX (NPSIZE,NWSIZE)
             NWSIZE  = NPSIZE
             ZONE12  = NPSIZE + NCSIZE
             ZONE4B  = 8 * MIJ + 3 * (NKIN1D + NOVL1D)
             ZONE4C  = NWSIZE
             ZONE4   = MAX (ZONE4B,ZONE4C)

             ZONEOPT = ZONE12 + ZONE3 + ZONE4

             MIJ     = 1
             NKIN1D  = MIJ * SKIN1D
             NOVL1D  = MIJ * SOVL1D
             NPSIZE  = NXYZT * MIJ
             WCOL    = MNPRIM
             NWSIZE  = NCTROW * WCOL
             NPSIZE  = MAX (NPSIZE,NWSIZE)
             NWSIZE  = NPSIZE
             ZONE12  = NPSIZE + NCSIZE
             ZONE4B  = 8 * MIJ + 3 * (NKIN1D + NOVL1D)
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

                MIJ = MAX0 (1, MIJBIG / IJDIV )

                NKIN1D  = MIJ * SKIN1D
                NOVL1D  = MIJ * SOVL1D
                NPSIZE  = NXYZT * MIJ
                WCOL    = MNPRIM
                NWSIZE  = NCTROW * WCOL
                NPSIZE  = MAX (NPSIZE,NWSIZE)
                NWSIZE  = NPSIZE
                ZONE12  = NPSIZE + NCSIZE
                ZONE4B  = 8 * MIJ + 3 * (NKIN1D + NOVL1D)
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

                    ZEA = ZRHOAB + NIJ
                    ZEB = ZEA + MIJ
                    ZE2AB = ZEB + MIJ
                    ZPAX = ZE2AB + MIJ
                    ZPAY = ZPAX + MIJ
                    ZPAZ = ZPAY + MIJ
                    ZPINVHF = ZPAZ + MIJ
                    ZSCALE = ZPINVHF + MIJ

                    ZKIN1DX = ZSCALE + MIJ
                    ZKIN1DY = ZKIN1DX + NKIN1D
                    ZKIN1DZ = ZKIN1DY + NKIN1D
                    ZOVL1DX = ZKIN1DZ + NKIN1D
                    ZOVL1DY = ZOVL1DX + NOVL1D
                    ZOVL1DZ = ZOVL1DY + NOVL1D

                    ZWORK = ZEA

                    RETURN
                END IF

                IJDIV = IJDIV + IJDIV

  100        CONTINUE

             WRITE (*,*) ' Memory allocation failed for (a|b) ! '
             WRITE (*,*) ' NIJ,MIJ = ',NIJ,MIJ
             WRITE (*,*) ' Increase flp memory! '
             WRITE (*,*) ' (oed__kin_ab_def_blocks) '

             STOP

         END IF
C
C
C             ...ready!
C
C
         END
