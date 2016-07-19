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
        SUBROUTINE  ERD__1111_DEF_BLOCKS
     +
     +                    ( ZMAX,
     +                      NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                      NIJ,NKL,
     +                      NRS,NTU,NRSTU,
     +                      NXYZT,
     +                      L1CACHE,NCTROW,
     +                      MEMORY,
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
C------------------------------------------------------------------------
C  OPERATION   : ERD__1111_DEF_BLOCKS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine handles the memory partitions for the
C                generation of the entire contracted (12|34) batches
C                for integrals over s- and/or p-functions only.
C                It determines the block size partitions (if any) for
C                the ij and kl exponent pairs of the whole primitive
C                [12|34] batch and returns pointers to the flp data
C                sections needed by the (12|34) generation.
C
C                The routine is also equipped with the possibility
C                of returning just the minimum / optimum flp memory
C                needed, without evaluating the (12|34) generation flp
C                pointers. This is useful for establishing just the
C                overall memory size needed for the (12|34) generation.
C                The keyword MEMORY has to be set true for this case
C                and obviously the value of ZMAX passed is irrelevant.
C
C
C                  Input:
C
C                    ZMAX         =  maximum flp memory
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x=1,2,3,4
C                    NIJ(KL)      =  total # of ij(kl) primitive
C                                    index pairs corresponding to
C                                    the csh pairs 1,2(3,4)
C                    NRS(TU)      =  total # of rs(tu) contraction
C                                    index pairs corresponding to
C                                    the csh pairs 1,2(3,4)
C                    NRSTU        =  total # of rstu contraction
C                                    index quadruplets (= NRS*NTU)
C                    NXYZT        =  total # of cartesian monomial
C                                    quadruplets
C                    L1CACHE      =  Size of level 1 cache in units of
C                                    8 Byte
C                    NCTROW       =  minimum # of rows that are
C                                    accepted for blocked contractions
C                    MEMORY       =  if this keyword is true, the
C                                    routine will only determine the
C                                    minimum / optimum flp memory and
C                                    store these values into NIJBLK and
C                                    NKLBLK, respectively (see below)
C
C                  Output:
C
C                    NIJ(KL)BLK   =  if MEMORY is false, they contain
C                                    the block sizes for the ij(kl)
C                                    primitive indices in order to
C                                    perform efficient contractions.
C                                    if MEMORY is true, they will have
C                                    the values for the minimum /
C                                    optimum flp memory
C                    NxSIZE       =  size of the primitive integral
C                                    block (x=P), contracted integrals
C                                    (x=C) and working (x=W) arrays
C                    MXPRIM       =  the maximum # of primitives
C                                    between all i,j,k,l primitives,
C                                    i.e. = max (i,j,k,l)
C                    MNPRIM       =  the minimum # of primitives
C                                    between i and j primitives and
C                                    k and l primitives and form the
C                                    maximum between these two values,
C                                    i.e. = max (min(i,j),min(k,l))
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
C                   Zone 1: final (12|34) batch
C                   Zone 2: block [12|34] and half transformed (12|34)
C                   Zone 3 : complete set of 1,2,3,4 norms and
C                            complete set (after screening!) of
C                            exponential prefactors 
C                   Zone 4 :  i) scratch for block [12|34] generation
C                            ii) scratch for partial (12|34) generation
C
C
C                Memory allocation offsets for the primitive [12|34]
C                batches generation + contraction to (12|34):
C
C
C                  --- Zone 1 ---
C
C                   ZCBATCH = offset for contracted (12|34) batch
C
C                  --- Zone 2 ---
C
C                   ZPBATCH = offset for (blocked) [12|34] batch (K4)
C
C                  --- Zone 3 ---
C
C                   ZNORM1 = offset for all 1 norms
C                   ZNORM2 = offset for all 2 norms
C                   ZNORM3 = offset for all 3 norms
C                   ZNORM4 = offset for all 4 norms
C
C                   ZRHO12 = offset for all 12 exponential prefactors
C                   ZRHO34 = offset for all 34 exponential prefactors
C
C                  --- Zone 4: for block [12|34] only ---
C
C                   ZP = offset for (blocked) P values (K2)
C                   ZPX = offset for (blocked) PX values (K2)
C                   ZPY = offset for (blocked) PY values (K2)
C                   ZPZ = offset for (blocked) PZ values (K2)
C                   ZSCPK2 = offset for (blocked) P scale values (K2)
C
C                   ZQ = offset for (blocked) Q values (K2)
C                   ZQX = offset for (blocked) QX values (K2)
C                   ZQY = offset for (blocked) QY values (K2)
C                   ZQZ = offset for (blocked) QZ values (K2)
C                   ZSCQK2 = offset for (blocked) Q scale values (K2)
C
C                  --- Zone 4: for contraction only ---
C
C                   ZWORK = offset for contraction working array
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

         INTEGER     IJDIV,KLDIV
         INTEGER     L1CACHE,NCTROW
         INTEGER     M,N
         INTEGER     MIJ,MKL,MIJKL,MRSKL
         INTEGER     MIJBIG,MKLBIG
         INTEGER     MXPRIM,MNPRIM
         INTEGER     MXSTEP
         INTEGER     NIJ,NKL
         INTEGER     NIJBLK,NKLBLK
         INTEGER     NPGTO1,NPGTO2,NPGTO3,NPGTO4
         INTEGER     NPMINRS,NPMINTU
         INTEGER     NPSIZE,NCSIZE,NWSIZE
         INTEGER     NRS,NTU,NRSTU
         INTEGER     NXYZT
         INTEGER     PFACT
         INTEGER     POW2IJ,POW2KL
         INTEGER     WCOL
         INTEGER     ZONE12,ZONE3,ZONE4,ZONE4B,ZONE4C
         INTEGER     ZMAX
         INTEGER     ZCBATCH,ZPBATCH,ZWORK,
     +               ZNORM1,ZNORM2,ZNORM3,ZNORM4,
     +               ZRHO12,ZRHO34,
     +               ZP,ZPX,ZPY,ZPZ,ZSCPK2,
     +               ZQ,ZQX,ZQY,ZQZ,ZSCQK2

         DOUBLE PRECISION   DLOG2

         PARAMETER  (DLOG2  =  0.6931471805599D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...determine optimum size of [12|34] blocks:
C
C                The next few lines of code are among the most
C                important for performance of the entire integral
C                package. Any badly chosen sizes for the primitive
C                blocks has a severe effect on the integral
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
C                  i) The optimum [12|34] block size is assumed to
C                     be directly proportional to the level 1 cache
C                     size.
C
C                 ii) The proportionality factor PFACT has to be
C                     determined by 'experiment', running a real
C                     molecular case with fully contracted basis sets.
C
C                iii) After having found the optimum # of exponent
C                     quadruplets per [12|34] block, it is best
C                     to choose the IJ-part as large as possible,
C                     keeping the KL-part size to a minimum. The
C                     reason behind this is that a KL-part size
C                     of 1 does not trigger an intermediate IJ <-> KL
C                     transposition during the contraction procedure.
C
C
         PFACT = 4

         MIJKL = PFACT * L1CACHE / NXYZT
         MIJKL = MAX (1,MIJKL)

         MIJ = MIN (MIJKL,    NIJ)
         MKL = MIN (MIJKL/MIJ,NKL)
C
C
C             ...the optimum block sizes MIJ x MKL have been found.
C                We have to see now how it fits with the rest of
C                the needed data into the maximum memory given.
C                If necessary subdivide the optimum block size
C                into convenient subblocks that fit into memory.
C                Note, that the primitive array must also hold the
C                intermediate halftransformed contracted integrals.
C
C                The subdivision of the MIJ x MKL block is done in
C                such a way that successive powers of 2 are checked
C                as divisors of MIJ and MKL. Only one of the divisors
C                is incremented in each trial step. The maximum #
C                of steps MXSTEP can thus be calculated beforehand
C                by knowing the fact that the smallest subdivided
C                MIJ x MKL block is of size 1 x 1. The succesive
C                divisors are checked in the following order:
C
C
C                        step  |  div of MIJ  |  div of MKL
C                       ------------------------------------
C                          1   |       1      |      1
C                          2   |       1      |      2
C                          3   |       2      |      2
C                          4   |       2      |      4
C                          5   |       4      |      4
C                          6   |       4      |      8
C                          7   |       8      |      8
C                         ...         ...           ...
C
C                The routine stops with a diagnostic, if even a
C                MIJ x MKL block of minimum size 1 x 1 cannot be
C                accomodated into the memory.
C
C
         NPMINRS = MIN (NPGTO1,NPGTO2)
         NPMINTU = MIN (NPGTO3,NPGTO4)

         NCSIZE  = NXYZT * NRSTU

         MXPRIM = MAX (NPGTO1,NPGTO2,NPGTO3,NPGTO4)
         MNPRIM = MAX (NPMINRS,NPMINTU)

         ZONE3 = NPGTO1 + NPGTO2 + NPGTO3 + NPGTO4 + NIJ + NKL
C
C
C             ...if the MEMORY keyword is activated, the routine
C                will only determine the optimum and minimum flp
C                memory (in that order), place them respectively
C                into the NKLBLK and NIJBLK variables and exit.
C
C
         IF (MEMORY) THEN

             MIJKL   = MIJ * MKL
             MRSKL   = NRS * MKL
             NPSIZE  = NXYZT * MAX (MIJKL,MRSKL)
             NWSIZE  = NPSIZE
             ZONE12  = NPSIZE + NCSIZE
             ZONE4B  = 5*(MIJ+MKL)
             WCOL    = MNPRIM
             ZONE4C  = NWSIZE + NCTROW * WCOL
             ZONE4   = MAX (ZONE4B,ZONE4C)

             NKLBLK  = ZONE12 + ZONE3 + ZONE4

             MIJ     = 1
             MKL     = 1
             MIJKL   = MIJ * MKL
             MRSKL   = NRS * MKL
             NPSIZE  = NXYZT * MAX (MIJKL,MRSKL)
             NWSIZE  = NPSIZE
             ZONE12  = NPSIZE + NCSIZE
             ZONE4B  = 5*(MIJ+MKL)
             WCOL    = MNPRIM
             ZONE4C  = NWSIZE + NCTROW * WCOL
             ZONE4   = MAX (ZONE4B,ZONE4C)

             NIJBLK  = ZONE12 + ZONE3 + ZONE4

             RETURN

         ELSE
C
C
C             ...the actual fitting into the maximum memory given.
C
C
             POW2IJ = DLOG ( DFLOAT (MIJ) ) / DLOG2
             POW2KL = DLOG ( DFLOAT (MKL) ) / DLOG2

             MXSTEP = 2*MIN (POW2IJ,POW2KL) + ABS (POW2IJ - POW2KL) + 1

             MIJBIG = MIJ
             MKLBIG = MKL

             IJDIV = 1
             KLDIV = 1

             DO 20 N = 1,MXSTEP

                MIJ = MAX (1, MIJBIG / IJDIV )
                MKL = MAX (1, MKLBIG / KLDIV )

                MIJKL   = MIJ * MKL
                MRSKL   = NRS * MKL
                NPSIZE  = NXYZT * MAX (MIJKL,MRSKL)
                NWSIZE  = NPSIZE
                ZONE12  = NPSIZE + NCSIZE
                ZONE4B  = 5*(MIJ+MKL)
                WCOL    = MNPRIM
                ZONE4C  = NWSIZE + NCTROW * WCOL
                ZONE4   = MAX (ZONE4B,ZONE4C)

                IF (ZONE12+ZONE3+ZONE4.LE.ZMAX) THEN
                    NIJBLK = MIJ
                    NKLBLK = MKL
                    NWSIZE = ZMAX - ZONE12 - ZONE3
C
C
C             ...generate the memory allocation pointers.
C
C
                    ZCBATCH = 1
                    ZPBATCH = NCSIZE + 1

                    ZNORM1 = ZPBATCH + NPSIZE
                    ZNORM2 = ZNORM1 + NPGTO1
                    ZNORM3 = ZNORM2 + NPGTO2
                    ZNORM4 = ZNORM3 + NPGTO3

                    ZRHO12 = ZNORM4 + NPGTO4
                    ZRHO34 = ZRHO12 + NIJ

                    ZP = ZRHO34 + NKL
                    ZPX = ZP + MIJ
                    ZPY = ZPX + MIJ
                    ZPZ = ZPY + MIJ
                    ZSCPK2 = ZPZ + MIJ

                    ZQ = ZSCPK2 + MIJ
                    ZQX = ZQ + MKL
                    ZQY = ZQX + MKL
                    ZQZ = ZQY + MKL
                    ZSCQK2 = ZQZ + MKL

                    ZWORK = ZP

                    RETURN
                END IF

                IF (MIJ.EQ.1) THEN
                    KLDIV = KLDIV * 2
                ELSE IF (MKL.EQ.1) THEN
                    IJDIV = IJDIV * 2
                ELSE
                    M = MOD (N+1,2)
                    IJDIV = IJDIV * (M+1)
                    KLDIV = KLDIV * (2-M)
                END IF

   20        CONTINUE

             WRITE (*,*) ' Memory allocation failed for (12|34) ! '
             WRITE (*,*) ' NIJ,NKL,MIJ,MKL = ',NIJ,NKL,MIJ,MKL
             WRITE (*,*) ' Increase flp memory! '
             WRITE (*,*) ' (erd__1111_def_blocks) '

             STOP

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
