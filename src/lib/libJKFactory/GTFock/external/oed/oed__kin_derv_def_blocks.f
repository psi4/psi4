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
         SUBROUTINE  OED__KIN_DERV_DEF_BLOCKS
     +
     +                    ( ZMAX,
     +                      NPGTOA,NPGTOB,
     +                      SHELLA,SHELLB,SHELLP,
     +                      NIJ,NRS,
     +                      NXYZT,
     +                      DERAX,DERAY,DERAZ,
     +                      DERBX,DERBY,DERBZ,
     +                      L1CACHE,NCTROW,
     +                      MEMORY,
     +
     +                              NIJBLK,
     +                              NBATCH,
     +                              NPSIZE,NCSIZE,NWSIZE,
     +                              NKIN1D,
     +                              NOVL1DX,NOVL1DY,NOVL1DZ,
     +                              MXPRIM,MNPRIM,
     +                              ZCBATCH,ZPBATCH,ZWORK,
     +                              ZNORMA,ZNORMB,
     +                              ZRHOAB,
     +                              ZEA,ZEB,ZE2A,ZE2B,ZE2AB,
     +                              ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCALE,
     +                              ZKIN1DX,ZKIN1DY,ZKIN1DZ,
     +                              ZOVL1DX,ZOVL1DY,ZOVL1DZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__KIN_DERV_DEF_BLOCKS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine handles the memory partitions for the
C                generation of the entire contracted derivative (a|b)
C                kinetic batch. It determines the block size partitions
C                (if any) for the ij exponent paris of the whole
C                primitive [a|b] batch and returns pointers to the flp
C                data sections needed by the (a|b) generation.
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
C                    NXYZT        =  total # of cartesian monomial
C                                    pairs
C                    DERxp        =  the order of differentiation on
C                                    centers x = A,B with respect to
C                                    the p = x,y,z coordinates
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
C                                    needed for 1D overlap integral
C                                    derivative evaluation, hence its
C                                    size might be larger than the
C                                    final primitive integral block
C                    NxSIZE       =  size of the primitive integral
C                                    block (x=P), contracted integrals
C                                    (x=C) and working (x=W) arrays.
C                                    If MEMORY is true, NPSIZE and
C                                    NCSIZE will contain respectively
C                                    the minimum and optimum flp memory
C                    NKIN1D       =  space needed for each of the
C                                    kinetic 1D X,Y,Z integral arrays
C                    NOVL1Dx      =  space needed for each of the 1D
C                                    overlap x = X,Y,Z integral arrays
C                                    (they might be different due to
C                                    different orders of differentiation
C                                    for each cartesian component)
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
C                            prefactors
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
C                   ZEA = offset for A-exponent values
C                   ZEB = offset for B-exponent values
C                   ZE2A = offset for 2 x A-exponent values
C                   ZE2B = offset for 2 x B-exponent values
C                   ZE2AB = offset for 2 x (A*B) exponent product
C                           values
C
C                   ZPAX = offset for (blocked) PAX values
C                   ZPAY = offset for (blocked) PAY values
C                   ZPAZ = offset for (blocked) PAZ values
C                   ZPINVHF = offset for (blocked) 1/2P values
C                   ZSCALE = offset for (blocked) scale values
C
C                   ZKIN1DX = offset for (blocked) derivative 1DX
C                             kinetic integrals accumulation
C                   ZKIN1DY = offset for (blocked) derivative 1DY
C                             kinetic integrals accumulation
C                   ZKIN1DZ = offset for (blocked) derivative 1DZ
C                             kinetic integrals accumulation
C                   ZOVL1DX = offset for (blocked) derivative 1DX
C                             overlap integrals accumulation
C                   ZOVL1DY = offset for (blocked) derivative 1DY
C                             overlap accumulation
C                   ZOVL1DZ = offset for (blocked) derivative 1DZ
C                             overlap accumulation
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
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     MEMORY

         INTEGER     DERAX,DERAY,DERAZ
         INTEGER     DERBX,DERBY,DERBZ
         INTEGER     IJDIV
         INTEGER     L1CACHE,NCTROW
         INTEGER     MIJ,MIJBIG
         INTEGER     MXPRIM,MNPRIM
         INTEGER     MXSTEP
         INTEGER     N
         INTEGER     NBATCH
         INTEGER     NIJ,NIJBLK
         INTEGER     NKIN1D
         INTEGER     NOVL1DX,NOVL1DY,NOVL1DZ
         INTEGER     NPGTOA,NPGTOB
         INTEGER     NPSIZE,NCSIZE,NWSIZE
         INTEGER     NRS
         INTEGER     NXYZT
         INTEGER     PFACT
         INTEGER     POW2IJ
         INTEGER     PSHELL
         INTEGER     SHELLA,SHELLB,SHELLP
         INTEGER     SHELLX,SHELLY,SHELLZ
         INTEGER     SKIN1D
         INTEGER     SOVL1DX,SOVL1DY,SOVL1DZ
         INTEGER     WCOL
         INTEGER     ZMAX
         INTEGER     ZONE1,ZONE2,ZONE3,ZONE3B,ZONE3C
         INTEGER     ZONEMIN,ZONEOPT
         INTEGER     ZCBATCH,ZPBATCH,ZWORK,
     +               ZNORMA,ZNORMB,ZRHOAB,
     +               ZEA,ZEB,ZE2A,ZE2B,ZE2AB,
     +               ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCALE,
     +               ZKIN1DX,ZKIN1DY,ZKIN1DZ,
     +               ZOVL1DX,ZOVL1DY,ZOVL1DZ

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
         SHELLX =   (SHELLP + DERAX + DERBX + 3)
     +            * MIN ((SHELLA + DERAX + 2),(SHELLB + DERBX + 2))
         SHELLY =   (SHELLP + DERAY + DERBY + 3)
     +            * MIN ((SHELLA + DERAY + 2),(SHELLB + DERBY + 2))
         SHELLZ =   (SHELLP + DERAZ + DERBZ + 3)
     +            * MIN ((SHELLA + DERAZ + 2),(SHELLB + DERBZ + 2))

         PSHELL = MAX (SHELLX,SHELLY,SHELLZ)

         SKIN1D  = (SHELLA         + 1) * (SHELLB         + 1)
         SOVL1DX = (SHELLA + DERAX + 2) * (SHELLB + DERBX + 2)
         SOVL1DY = (SHELLA + DERAY + 2) * (SHELLB + DERBY + 2)
         SOVL1DZ = (SHELLA + DERAZ + 2) * (SHELLB + DERBZ + 2)

         MXPRIM = MAX (NPGTOA,NPGTOB)
         MNPRIM = MIN (NPGTOA,NPGTOB)

         NCSIZE  = NXYZT * NRS

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

             NKIN1D  = MIJ * SKIN1D
             NOVL1DX = MIJ * SOVL1DX
             NOVL1DY = MIJ * SOVL1DY
             NOVL1DZ = MIJ * SOVL1DZ
             NBATCH  = MIJ * MAX (PSHELL,NXYZT)
             NPSIZE  = NXYZT * MIJ
             WCOL    = MNPRIM
             NWSIZE  = NCTROW * WCOL
             NPSIZE  = MAX (NPSIZE,NWSIZE)
             NWSIZE  = NPSIZE
             ZONE1   = NCSIZE
             ZONE3B  = NBATCH + 3*NKIN1D + 10*MIJ
     +                 + NOVL1DX + NOVL1DY + NOVL1DZ
             ZONE3C  = NPSIZE + NWSIZE
             ZONE3   = MAX (ZONE3B,ZONE3C)

             ZONEOPT = ZONE1 + ZONE2 + ZONE3

             MIJ     = 1
             NKIN1D  = MIJ * SKIN1D
             NOVL1DX = MIJ * SOVL1DX
             NOVL1DY = MIJ * SOVL1DY
             NOVL1DZ = MIJ * SOVL1DZ
             NBATCH  = MIJ * MAX (PSHELL,NXYZT)
             NPSIZE  = NXYZT * MIJ
             WCOL    = MNPRIM
             NWSIZE  = NCTROW * WCOL
             NPSIZE  = MAX (NPSIZE,NWSIZE)
             NWSIZE  = NPSIZE
             ZONE1   = NCSIZE
             ZONE3B  = NBATCH + 3*NKIN1D + 10*MIJ
     +                 + NOVL1DX + NOVL1DY + NOVL1DZ
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

                NKIN1D  = MIJ * SKIN1D
                NOVL1DX = MIJ * SOVL1DX
                NOVL1DY = MIJ * SOVL1DY
                NOVL1DZ = MIJ * SOVL1DZ
                NBATCH  = MIJ * MAX (PSHELL,NXYZT)
                NPSIZE  = NXYZT * MIJ
                WCOL    = MNPRIM
                NWSIZE  = NCTROW * WCOL
                NPSIZE  = MAX (NPSIZE,NWSIZE)
                NWSIZE  = NPSIZE
                ZONE1   = NCSIZE
                ZONE3B  = NBATCH + 3*NKIN1D + 10*MIJ
     +                    + NOVL1DX + NOVL1DY + NOVL1DZ
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

                    ZEA = ZPBATCH + NBATCH
                    ZEB = ZEA + MIJ
                    ZE2A = ZEB + MIJ
                    ZE2B = ZE2A + MIJ
                    ZE2AB = ZE2B + MIJ

                    ZPAX = ZE2AB + MIJ
                    ZPAY = ZPAX + MIJ
                    ZPAZ = ZPAY + MIJ
                    ZPINVHF = ZPAZ + MIJ
                    ZSCALE = ZPINVHF + MIJ

                    ZKIN1DX = ZSCALE + MIJ
                    ZKIN1DY = ZKIN1DX + NKIN1D
                    ZKIN1DZ = ZKIN1DY + NKIN1D
                    ZOVL1DX = ZKIN1DZ + NKIN1D
                    ZOVL1DY = ZOVL1DX + NOVL1DX
                    ZOVL1DZ = ZOVL1DY + NOVL1DY

                    ZWORK = ZPBATCH + NPSIZE

                    RETURN
                END IF

                IJDIV = IJDIV + IJDIV

  100        CONTINUE

             WRITE (*,*) ' Memory allocation failed for derv (a|b) ! '
             WRITE (*,*) ' NIJ,MIJ = ',NIJ,MIJ
             WRITE (*,*) ' Increase flp memory! '
             WRITE (*,*) ' (oed__kin_derv_def_blocks) '

             STOP

         END IF
C
C
C             ...ready!
C
C
         END
