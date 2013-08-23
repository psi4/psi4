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
         SUBROUTINE  ERD__E0F0_DEF_BLOCKS
     +
     +               ( ZMAX,
     +                 NPGTOA,NPGTOB,NPGTOC,NPGTOD,
     +                 SHELLP,SHELLQ,
     +                 NIJ,NKL,
     +                 NRS,NTU,NRSTU,
     +                 NGQP,NGQSCR,
     +                 NXYZT,
     +                 L1CACHE,NCTROW,
     +                 MEMORY,
     +
     +                     NIJBLK,NKLBLK,
     +                     NPSIZE,NCSIZE,NWSIZE,
     +                     NINT2D,
     +                     MXPRIM,MNPRIM,
     +                     ZCBATCH,ZPBATCH,ZWORK,
     +                     ZNORMA,ZNORMB,ZNORMC,ZNORMD,
     +                     ZRHOAB,ZRHOCD,
     +                     ZP,ZPX,ZPY,ZPZ,ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCPK2,
     +                     ZQ,ZQX,ZQY,ZQZ,ZQCX,ZQCY,ZQCZ,ZQINVHF,ZSCQK2,
     +                     ZRTS,ZWTS,ZGQSCR,ZTVAL,ZPQPINV,ZSCPQK4,
     +                     ZB00,ZB01,ZB10,
     +                     ZC00X,ZC00Y,ZC00Z,ZD00X,ZD00Y,ZD00Z,
     +                     ZINT2DX,ZINT2DY,ZINT2DZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__E0F0_DEF_BLOCKS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine handles the memory partitions for the
C                generation of the entire contracted (e0|f0) batch.
C                It determines the block size partitions (if any) for
C                the ij and kl exponent paris of the whole primitive
C                [e0|f0] batch and returns pointers to the flp data
C                sections needed by the (e0|f0) generation.
C
C                The routine is also equipped with the possibility
C                of returning just the minimum / optimum flp memory
C                needed, without evaluating the (e0|f0) generation flp
C                pointers. This is useful for establishing just the
C                overall memory size needed for the (e0|f0) generation.
C                The keyword MEMORY has to be set true for this case
C                and obviously the value of ZMAX passed is irrelevant.
C
C
C                  Input:
C
C                    ZMAX         =  maximum flp memory
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x=A,B,C,D
C                    SHELLP(Q)    =  the shell sums for csh: A+B (C+D)
C                    NIJ(KL)      =  total # of ij(kl) primitive
C                                    index pairs corresponding to
C                                    the csh pairs A,B(C,D)
C                    NRS(TU)      =  total # of rs(tu) contraction
C                                    index pairs corresponding to
C                                    the csh pairs A,B(C,D)
C                    NRSTU        =  total # of rstu contraction
C                                    index quadruplets (= NRS*NTU)
C                    NGQP         =  # of gaussian quadrature points
C                                    (roots)
C                    NGQSCR       =  size of gaussian quadrature
C                                    scratch space needed to calculate
C                                    all the quadrature roots
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
C                    NINT2D       =  space needed for the 2D X,Y,Z
C                                    integral arrays
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
C                   Zone 1 : final (e0|f0) batch
C                   Zone 2 : block [e0|f0] and half transformed (e0|f0)
C                   Zone 3 : complete set of A,B,C,D norms and
C                            complete set of AB and CD exponential
C                            prefactors 
C                   Zone 4 : i) scratch for block [e0|f0] generation
C                           ii) scratch for partial (e0|f0) generation
C
C
C                Memory allocation offsets for the primitive [e0|f0]
C                batches generation + contraction to (e0|f0):
C
C
C                  --- Zone 1 ---
C
C                   ZCBATCH = offset for contracted (e0|f0) batch
C
C                  --- Zone 2 ---
C
C                   ZPBATCH = offset for (blocked) [e0|f0] batch (K4)
C
C                  --- Zone 3 ---
C
C                   ZNORMA = offset for all A norms
C                   ZNORMB = offset for all B norms
C                   ZNORMC = offset for all C norms
C                   ZNORMD = offset for all D norms
C
C                   ZRHOAB = offset for all AB exponential prefactors
C                   ZRHOCD = offset for all CD exponential prefactors
C
C                  --- Zone 4: for block [e0|f0] generation only ---
C
C                   ZP = offset for (blocked) P values (K2)
C                   ZPX = offset for (blocked) PX values (K2)
C                   ZPY = offset for (blocked) PY values (K2)
C                   ZPZ = offset for (blocked) PZ values (K2)
C                   ZPAX = offset for (blocked) PAX values (K2)
C                   ZPAY = offset for (blocked) PAY values (K2
C                   ZPAZ = offset for (blocked) PAZ values (K2)
C                   ZPINVHF = offset for (blocked) 1/2P values (K2)
C                   ZSCPK2 = offset for (blocked) P scale values (K2)
C
C                   ZQ = offset for (blocked) Q values (K2)
C                   ZQX = offset for (blocked) QX values (K2)
C                   ZQY = offset for (blocked) QY values (K2)
C                   ZQZ = offset for (blocked) QZ values (K2)
C                   ZQCX = offset for (blocked) QCX values (K2)
C                   ZQCY = offset for (blocked) QCY values (K2
C                   ZQCZ = offset for (blocked) QCZ values (K2)
C                   ZQINVHF = offset for (blocked) 1/2Q values (K2)
C                   ZSCQK2 = offset for (blocked) Q scale values (K2)
C
C                   ZRTS = offset for (blocked) quad roots (K4)
C                   ZWTS = offset for (blocked) quad weights (K4)
C                   ZGQSCR = offset for (blocked) quad scratch (K4)
C                   ZTVAL = offset for (blocked) T exponents (K4)
C                   ZPQPINV = offset for (blocked) 1/(P+Q) values (K4)
C                   ZSCPQK4 = offset for (blocked) PQ scale values (K4)
C
C                   ZB00 = offset for (blocked) VRR coeff B00 (K4)
C                   ZB01 = offset for (blocked) VRR coeff B01 (K4)
C                   ZB10 = offset for (blocked) VRR coeff B10 (K4)
C                   ZC00X = offset for (blocked) VRR coeff C00X (K4)
C                   ZC00Y = offset for (blocked) VRR coeff C00Y (K4)
C                   ZC00Z = offset for (blocked) VRR coeff C00Z (K4)
C                   ZD00X = offset for (blocked) VRR coeff D00X (K4)
C                   ZD00Y = offset for (blocked) VRR coeff D00Y (K4)
C                   ZD00Z = offset for (blocked) VRR coeff D00Z (K4)
C
C                   ZINT2DX = offset for (blocked) 2DX integrals (K4)
C                   ZINT2DY = offset for (blocked) 2DY integrals (K4)
C                   ZINT2DZ = offset for (blocked) 2DZ integrals (K4)
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

         INTEGER     IJDIV,KLDIV
         INTEGER     L1CACHE,NCTROW
         INTEGER     M,N
         INTEGER     MIJ,MKL,MIJKL,MRSKL
         INTEGER     MIJBIG,MKLBIG
         INTEGER     MXPRIM,MNPRIM
         INTEGER     MXSTEP
         INTEGER     NIJ,NKL
         INTEGER     NIJBLK,NKLBLK
         INTEGER     NINT2D
         INTEGER     NGQP,NGQSCR,MGQIJKL
         INTEGER     NPGTOA,NPGTOB,NPGTOC,NPGTOD
         INTEGER     NPMINRS,NPMINTU
         INTEGER     NPSIZE,NCSIZE,NWSIZE
         INTEGER     NRS,NTU,NRSTU
         INTEGER     NXYZT
         INTEGER     PFACT
         INTEGER     POW2IJ,POW2KL
         INTEGER     SHELLP,SHELLQ
         INTEGER     WCOL
         INTEGER     ZONE12,ZONE3,ZONE4,ZONE4B,ZONE4C
         INTEGER     ZMAX
         INTEGER     ZCBATCH,ZPBATCH,ZWORK,
     +               ZNORMA,ZNORMB,ZNORMC,ZNORMD,
     +               ZRHOAB,ZRHOCD,
     +               ZP,ZPX,ZPY,ZPZ,ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCPK2,
     +               ZQ,ZQX,ZQY,ZQZ,ZQCX,ZQCY,ZQCZ,ZQINVHF,ZSCQK2,
     +               ZRTS,ZWTS,ZGQSCR,ZTVAL,ZPQPINV,ZSCPQK4,
     +               ZB00,ZB01,ZB10,
     +               ZC00X,ZC00Y,ZC00Z,ZD00X,ZD00Y,ZD00Z,
     +               ZINT2DX,ZINT2DY,ZINT2DZ

         DOUBLE PRECISION   DLOG2

         PARAMETER  (DLOG2  =  0.6931471805599D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...determine optimum size of [e0|f0] blocks:
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
C                  i) The optimum [e0|f0] block size is assumed to
C                     be directly proportional to the level 1 cache
C                     size.
C
C                 ii) The proportionality factor PFACT has to be
C                     determined by 'experiment', running a real
C                     molecular case with fully contracted basis sets.
C
C                iii) After having found the optimum # of exponent
C                     quadruplets per [e0|f0] block, it is best
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
C
C             ...the optimum block sizes MIJ x MKL have been found.
C                We have to see now how it fits with the rest of
C                the needed data into the maximum memory given.
C                If necessary subdivide the optimum MIJ x MKL block
C                size into convenient subblocks that fit into memory.
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
C
         NPMINRS = MIN (NPGTOA,NPGTOB)
         NPMINTU = MIN (NPGTOC,NPGTOD)

         NCSIZE = NXYZT * NRSTU

         MXPRIM = MAX (NPGTOA,NPGTOB,NPGTOC,NPGTOD)
         MNPRIM = MAX (NPMINRS,NPMINTU)

         ZONE3 = NPGTOA + NPGTOB + NPGTOC + NPGTOD + NIJ + NKL
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
             MGQIJKL = NGQP * MIJKL
             NINT2D  = MGQIJKL * (SHELLP+1) * (SHELLQ+1)
             MRSKL   = NRS * MKL
             NPSIZE  = NXYZT * MAX (MIJKL,MRSKL)
             NWSIZE  = NPSIZE
             ZONE12  = NPSIZE + NCSIZE
             ZONE4B  = NGQSCR + 2*MIJKL+9*(MIJ+MKL)+12*MGQIJKL+3*NINT2D
             WCOL    = MNPRIM
             ZONE4C  = NWSIZE + NCTROW * WCOL
             ZONE4   = MAX (ZONE4B,ZONE4C)

             NKLBLK  = ZONE12 + ZONE3 + ZONE4

             MIJ     = 1
             MKL     = 1
             MIJKL   = MIJ * MKL
             MGQIJKL = NGQP * MIJKL
             NINT2D  = MGQIJKL * (SHELLP+1) * (SHELLQ+1)
             MRSKL   = NRS * MKL
             NPSIZE  = NXYZT * MAX (MIJKL,MRSKL)
             NWSIZE  = NPSIZE
             ZONE12  = NPSIZE + NCSIZE
             ZONE4B  = NGQSCR + 2*MIJKL+9*(MIJ+MKL)+12*MGQIJKL+3*NINT2D
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
                MGQIJKL = NGQP * MIJKL
                NINT2D  = MGQIJKL * (SHELLP+1) * (SHELLQ+1)
                MRSKL   = NRS * MKL
                NPSIZE  = NXYZT * MAX (MIJKL,MRSKL)
                NWSIZE  = NPSIZE
                ZONE12  = NPSIZE + NCSIZE
                ZONE4B  = NGQSCR + 2*MIJKL + 9*(MIJ+MKL) + 12*MGQIJKL
     +                           + 3*NINT2D
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

                    ZNORMA = ZPBATCH + NPSIZE
                    ZNORMB = ZNORMA + NPGTOA
                    ZNORMC = ZNORMB + NPGTOB
                    ZNORMD = ZNORMC + NPGTOC

                    ZRHOAB = ZNORMD + NPGTOD
                    ZRHOCD = ZRHOAB + NIJ

                    ZP = ZRHOCD + NKL
                    ZPX = ZP + MIJ
                    ZPY = ZPX + MIJ
                    ZPZ = ZPY + MIJ
                    ZPAX = ZPZ + MIJ
                    ZPAY = ZPAX + MIJ
                    ZPAZ = ZPAY + MIJ
                    ZPINVHF = ZPAZ + MIJ
                    ZSCPK2 = ZPINVHF + MIJ

                    ZQ = ZSCPK2 + MIJ
                    ZQX = ZQ + MKL
                    ZQY = ZQX + MKL
                    ZQZ = ZQY + MKL
                    ZQCX = ZQZ + MKL
                    ZQCY = ZQCX + MKL
                    ZQCZ = ZQCY + MKL
                    ZQINVHF = ZQCZ + MKL
                    ZSCQK2 = ZQINVHF + MKL

                    ZRTS = ZSCQK2 + MKL
                    ZWTS = ZRTS + MGQIJKL
                    ZGQSCR = ZWTS + MGQIJKL
                    ZTVAL = ZGQSCR + NGQSCR
                    ZPQPINV = ZTVAL + MIJKL
                    ZSCPQK4 = ZPQPINV + MIJKL

                    ZB00 = ZSCPQK4 + MGQIJKL
                    ZB01 = ZB00 + MGQIJKL
                    ZB10 = ZB01 + MGQIJKL
                    ZC00X = ZB10 + MGQIJKL
                    ZC00Y = ZC00X + MGQIJKL
                    ZC00Z = ZC00Y + MGQIJKL
                    ZD00X = ZC00Z + MGQIJKL
                    ZD00Y = ZD00X + MGQIJKL
                    ZD00Z = ZD00Y + MGQIJKL

                    ZINT2DX = ZD00Z + MGQIJKL
                    ZINT2DY = ZINT2DX + NINT2D
                    ZINT2DZ = ZINT2DY + NINT2D

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

             WRITE (*,*) ' Memory allocation failed for (e0|f0) ! '
             WRITE (*,*) ' NIJ,NKL,MIJ,MKL = ',NIJ,NKL,MIJ,MKL
             WRITE (*,*) ' Increase flp memory! '
             WRITE (*,*) ' (erd__e0f0_def_blocks) '

             STOP

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
