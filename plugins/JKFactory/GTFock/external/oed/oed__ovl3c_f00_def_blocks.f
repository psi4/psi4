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
         SUBROUTINE  OED__OVL3C_F00_DEF_BLOCKS
     +
     +                    ( ZMAX,
     +                      NPGTOI,NPGTOJ,NPGTOK,
     +                      SHELLQ,
     +                      NIJ,NK,NIJK,
     +                      NRS,NT,NRST,
     +                      NXYZT,
     +                      L1CACHE,NCTROW,
     +                      MEMORY,
     +
     +                              NIJBLK,
     +                              NPSIZE,NCSIZE,NWSIZE,
     +                              NINT1D,
     +                              MXPRIM,MNPRIM,
     +                              ZCBATCH,ZPBATCH,ZWORK,
     +                              ZNORMI,ZNORMJ,ZNORMK,
     +                              ZRHOIJK,
     +                              ZQAX,ZQAY,ZQAZ,ZQINVHF,ZSCALE,
     +                              ZINT1DX,ZINT1DY,ZINT1DZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__OVL3C_F00_DEF_BLOCKS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine handles the memory partitions for the
C                generation of the entire contracted (f00) 3-center
C                overlap batch. It determines the block size partitions
C                (if any) for the ij exponent pairs of the whole
C                primitive [f00] batch and returns pointers to the
C                flp data sections needed by the (f00) generation.
C                The single k exponent set will be treated as one
C                entity and will not be blocked.
C
C                The routine is also equipped with the possibility
C                of returning just the minimum / optimum flp memory
C                needed, without evaluating the (f00) generation flp
C                pointers. This is useful for establishing just the
C                overall memory size needed for the (f00) generation.
C                The keyword MEMORY has to be set true for this case
C                and obviously the value of ZMAX passed is irrelevant.
C
C
C                  Input:
C
C                    ZMAX         =  maximum flp memory
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x=I,J,K
C                    SHELLQ       =  the shell sum I+J+K
C                    NIJ          =  total # of ij primitive index
C                                    pairs corresponding to the csh
C                                    pairs I,J
C                    NK           =  total # of k exponent indices
C                                    corresponding to the csh K
C                    NIJK         =  # of ij primitive indices times
C                                    # of k exponent indices
C                    NRS          =  total # of rs contraction index
C                                    pairs corresponding to the csh
C                                    pairs I,J
C                    NT           =  total # of t contraction indices
C                                    corresponding to the csh K
C                    NRST         =  # of rs contraction indices times
C                                    # of t contraction indices
C                    NXYZT        =  total # of cartesian monomial
C                                    triples
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
C                    NINT1D       =  space needed for the 1D X,Y,Z
C                                    integral arrays
C                    MXPRIM       =  the maximum # of primitives
C                                    between all i,j and k primitives,
C                                    i.e. = max (i,j,k)
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
C                   Zone 1: final (f00) batch
C                   Zone 2: block [f00]
C                   Zone 3 : complete set of I,J,K norms and complete
C                            set of IJK exponential prefactors 
C                   Zone 4 : i) scratch for block [f00] generation
C                           ii) scratch for (f00) generation
C
C
C                Memory allocation offsets for the primitive [f00]
C                batches generation + contraction to (f00):
C
C
C                  --- Zone 1 ---
C
C                   ZCBATCH = offset for contracted (f00) batch
C
C                  --- Zone 2 ---
C
C                   ZPBATCH = offset for (blocked) [f00] batch
C
C                  --- Zone 3 ---
C
C                   ZNORMI = offset for all I norms
C                   ZNORMJ = offset for all J norms
C                   ZNORMK = offset for all K norms
C                   ZRHOIJK = offset for all IJK exponential prefactors
C
C                  --- Zone 4: for block [f00] generation only ---
C
C                   ZQAX = offset for (blocked) QAX values
C                   ZQAY = offset for (blocked) QAY values
C                   ZQAZ = offset for (blocked) QAZ values
C                   ZQINVHF = offset for (blocked) 1/2Q values
C                   ZSCALE = offset for (blocked) scale values
C
C                   ZINT1DX = offset for (blocked) 1DX integrals
C                   ZINT1DY = offset for (blocked) 1DY integrals
C                   ZINT1DZ = offset for (blocked) 1DZ integrals
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
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     MEMORY

         INTEGER     IJDIV
         INTEGER     L1CACHE,NCTROW
         INTEGER     N
         INTEGER     MIJ,MIJK,MIJT,MIJBIG
         INTEGER     MXPRIM,MNPRIM
         INTEGER     MXSTEP
         INTEGER     NIJ,NIJBLK,NK,NIJK
         INTEGER     NINT1D
         INTEGER     NPGTOI,NPGTOJ,NPGTOK
         INTEGER     NPSIZE,NCSIZE,NWSIZE
         INTEGER     NRS,NT,NRST
         INTEGER     NXYZT
         INTEGER     PFACT
         INTEGER     POW2IJ
         INTEGER     SHELLQ
         INTEGER     WCOL
         INTEGER     ZMAX
         INTEGER     ZONE12,ZONE3,ZONE4,ZONE4B,ZONE4C
         INTEGER     ZONEMIN,ZONEOPT
         INTEGER     ZCBATCH,ZPBATCH,ZWORK,
     +               ZNORMI,ZNORMJ,ZNORMK,ZRHOIJK,
     +               ZQAX,ZQAY,ZQAZ,ZQINVHF,ZSCALE,
     +               ZINT1DX,ZINT1DY,ZINT1DZ

         DOUBLE PRECISION   DLOG2

         PARAMETER  (DLOG2  =  0.6931471805599D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...determine optimum size of [f00] blocks:
C
C                The next few lines of code are among the most
C                important for performance of the entire 3-center
C                overlap integral package. Any badly chosen sizes
C                for the primitive blocks has a severe effect on the
C                integral evaluation timing. Two opposite effects have
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
C                  i) The optimum [f00] block size is assumed to
C                     be directly proportional to the level 1 cache
C                     size.
C
C                 ii) The proportionality factor PFACT has to be
C                     determined by 'experiment', running a real
C                     molecular case with fully contracted basis sets.
C
C
         PFACT = 4

         MIJ = (PFACT * L1CACHE) / (NXYZT * NK)
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
         MXPRIM = MAX (NPGTOI,NPGTOJ,NPGTOK)
         MNPRIM = MIN (NPGTOI,NPGTOJ)

         NCSIZE = NXYZT * NRST

         ZONE3 = NPGTOI + NPGTOJ + NPGTOK + NIJK
C
C
C             ...if the MEMORY keyword is activated, the routine
C                will only determine the optimum and minimum flp
C                memory (in that order), place them respectively
C                into the NCSIZE and NPSIZE variables and exit.
C
C
         IF (MEMORY) THEN

             MIJK    = MIJ * NK
             MIJT    = MIJ * NT
             NINT1D  = MIJK * (SHELLQ+1)
             NPSIZE  = NXYZT * MAX (MIJK,MIJT)
             WCOL    = MNPRIM
             NWSIZE  = NCTROW * WCOL
             NPSIZE  = MAX (NPSIZE,NWSIZE)
             NWSIZE  = NPSIZE
             ZONE12  = NPSIZE + NCSIZE
             ZONE4B  = 5 * MIJK + 3 * NINT1D
             ZONE4C  = NWSIZE
             ZONE4   = MAX (ZONE4B,ZONE4C)

             ZONEOPT = ZONE12 + ZONE3 + ZONE4

             MIJ     = 1
             MIJK    = MIJ * NK
             MIJT    = MIJ * NT
             NINT1D  = MIJK * (SHELLQ+1)
             NPSIZE  = NXYZT * MAX (MIJK,MIJT)
             WCOL    = MNPRIM
             NWSIZE  = NCTROW * WCOL
             NPSIZE  = MAX (NPSIZE,NWSIZE)
             NWSIZE  = NPSIZE
             ZONE12  = NPSIZE + NCSIZE
             ZONE4B  = 5 * MIJK + 3 * NINT1D
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
                MIJK = MIJ * NK
                MIJT = MIJ * NT

                NINT1D  = MIJK * (SHELLQ+1)
                NPSIZE  = NXYZT * MAX (MIJK,MIJT)
                WCOL    = MNPRIM
                NWSIZE  = NCTROW * WCOL
                NPSIZE  = MAX (NPSIZE,NWSIZE)
                NWSIZE  = NPSIZE
                ZONE12  = NPSIZE + NCSIZE
                ZONE4B  = 5 * MIJK + 3 * NINT1D
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

                    ZNORMI = ZPBATCH + NPSIZE
                    ZNORMJ = ZNORMI + NPGTOI
                    ZNORMK = ZNORMJ + NPGTOJ
                    ZRHOIJK = ZNORMK + NPGTOK

                    ZQAX = ZRHOIJK + NIJK
                    ZQAY = ZQAX + MIJK
                    ZQAZ = ZQAY + MIJK
                    ZQINVHF = ZQAZ + MIJK
                    ZSCALE = ZQINVHF + MIJK

                    ZINT1DX = ZSCALE + MIJK
                    ZINT1DY = ZINT1DX + NINT1D
                    ZINT1DZ = ZINT1DY + NINT1D

                    ZWORK = ZQAX

                    RETURN
                END IF

                IJDIV = IJDIV + IJDIV

  100        CONTINUE

             WRITE (*,*) ' Memory allocation failed for (f00) ! '
             WRITE (*,*) ' NIJ,MIJ = ',NIJ,MIJ
             WRITE (*,*) ' Increase flp memory! '
             WRITE (*,*) ' (oed__ovl3c_f00_def_blocks) '

             STOP

         END IF
C
C
C             ...ready!
C
C
         END
