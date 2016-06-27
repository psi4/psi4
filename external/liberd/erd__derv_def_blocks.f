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
         SUBROUTINE  ERD__DERV_DEF_BLOCKS
     +
     +               ( ZMAX,
     +                 NPGTOA,NPGTOB,NPGTOC,NPGTOD,
     +                 SHELLA,SHELLB,SHELLC,SHELLD,
     +                 SHELLP,SHELLQ,
     +                 NIJ,NKL,
     +                 NRS,NTU,NRSTU,
     +                 NGQP,NGQSCR,
     +                 NXYZT,
     +                 DERAX,DERAY,DERAZ,
     +                 DERBX,DERBY,DERBZ,
     +                 DERCX,DERCY,DERCZ,
     +                 DERDX,DERDY,DERDZ,
     +                 DERPX,DERPY,DERPZ,
     +                 DERQX,DERQY,DERQZ,
     +                 PRIMTYP,
     +                 L1CACHE,NCTROW,
     +                 MEMORY,
     +
     +                     NIJBLK,NKLBLK,
     +                     NBATCH,
     +                     NPSIZE,NCSIZE,NWSIZE,
     +                     NINT2DX,NINT2DY,NINT2DZ,
     +                     MXPRIM,MNPRIM,
     +                     ZCBATCH,ZPBATCH,ZWORK,
     +                     ZNORMA,ZNORMB,ZNORMC,ZNORMD,
     +                     ZRHOAB,ZRHOCD,
     +                     ZP,ZPX,ZPY,ZPZ,ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCPK2,
     +                     ZQ,ZQX,ZQY,ZQZ,ZQCX,ZQCY,ZQCZ,ZQINVHF,ZSCQK2,
     +                     ZRTS,ZWTS,ZGQSCR,ZTVAL,ZPQPINV,ZSCPQK4,
     +                     ZB00,ZB01,ZB10,
     +                     ZC00X,ZC00Y,ZC00Z,ZD00X,ZD00Y,ZD00Z,
     +                     ZEXP2A,ZEXP2B,ZEXP2C,ZEXP2D,
     +                     ZINT2DX,ZINT2DY,ZINT2DZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__DERV_DEF_BLOCKS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine handles the memory partitions for the
C                generation of the entire contracted derivated (e0|cd)
C                or (ab|cd) batch. It determines the block size
C                partitions (if any) for the ij and kl exponent paris
C                of the whole primitive derivated [e0|cd] or [ab|cd]
C                batch and returns pointers to the flp data sections
C                needed by the (e0|cd) or (ab|cd) generation.
C
C                The routine is also equipped with the possibility
C                of returning just the minimum / optimum flp memory
C                needed, without evaluating the (e0|cd) or (ab|cd)
C                generation flp pointers. This is useful for
C                establishing just the overall memory size needed for
C                the (e0|cd) or (ab|cd) generation. The keyword MEMORY
C                has to be set true for this case and obviously the
C                value of ZMAX passed is irrelevant.
C
C
C                  Input:
C
C                    ZMAX         =  maximum flp memory
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x=A,B,C,D
C                    SHELLx       =  the shell type for csh x=A,B,C,D
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
C                    DERyp        =  the order of differentiation on
C                                    centers y = A,B,C,D with respect
C                                    to the p = x,y,z coordinates
C                                    and their sums y=P,Q=A+B,C+D
C                    PRIMTYP      =  character variable, indicating
C                                    which type of primitives will
C                                    be generated and contracted. Can
C                                    be only 'E0CD' or 'ABCD'
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
C                    NIJ(KL)BLK   =  block sizes for the ij(kl)
C                                    primitive indices in order to
C                                    perform efficient contractions.
C                    NBATCH       =  size of the integral batch array
C                                    during primitive integral block
C                                    generation. This array is also
C                                    needed for 2D integral derivative
C                                    evaluation, hence its size might
C                                    be larger than the final primitive
C                                    integral block
C                    NxSIZE       =  size of the primitive integral
C                                    block (x=P), contracted integrals
C                                    (x=C) and working (x=W) arrays
C                    NINT2Dx      =  space needed for the 2D x=X,Y,Z
C                                    integral arrays. They are usually
C                                    not equal due to derivations
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
C                    |  Zone 1  |  Zone 2  |  Zone 3  |
C
C
C                   Zone 1: final contracted (e0|cd) or (ab|cd) batch
C                   Zone 2: complete set of A,B,C,D norms and
C                           complete set (after screening!) of
C                           exponential prefactors 
C                   Zone 3: will be used at two different stages:
C                            i) during derivated primitive [e0|cd]
C                               or [ab|cd] block generation
C                           ii) during the contraction stage
C
C
C                While the first two Zones are straightforward in
C                their definition, Zone 3 requires some more detailed
C                explanation:
C
C                 Stage i) : During this stage the primitive blocks
C                            have to be assembled from the derivated
C                            2D integrals, which in turn are generated
C                            in single derivative steps. The derivation
C                            of the 2D integrals needs an extra work
C                            space for accumulation, which will be
C                            provided at the same memory location where
C                            the final primitive block will be. This
C                            can be done, because the primitive block
C                            can only be assembled after! all the
C                            derivated 2D integrals are ready.
C
C                Stage ii) : This stage needs the primitive blocks
C                            plus some workspace for performing the
C                            two partial contraction steps.
C
C
C                Memory allocation offsets for the primitive [e0|cd]/
C                [ab|cd] batches generation + contraction to (e0|cd)/
C                (ab|cd):
C
C
C                  --- Zone 1 ---
C
C                   ZCBATCH = offset for contracted (e0|cd) or (ab|cd)
C                             batch
C
C                  --- Zone 2 ---
C
C                   ZNORMA = offset for all A norms
C                   ZNORMB = offset for all B norms
C                   ZNORMC = offset for all C norms
C                   ZNORMD = offset for all D norms
C
C                   ZRHOAB = offset for all AB exponential prefactors
C                   ZRHOCD = offset for all CD exponential prefactors
C
C                  --- Zone 3: for Stage i) only ---
C
C                   ZPBATCH = offset for derivated 2D integrals
C                             accumulation work space and for
C                             (blocked) derivated [e0|cd] or [ab|cd]
C                             batch (K4)
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
C                   ZEXP2A = offset for 2 x A-exponent values (K4)
C                   ZEXP2B = offset for 2 x B-exponent values (K4)
C                   ZEXP2C = offset for 2 x C-exponent values (K4)
C                   ZEXP2D = offset for 2 x D-exponent values (K4)
C
C                   ZINT2DX = offset for (blocked) derivative 2DX
C                             integrals accumulation (K4)
C                   ZINT2DY = offset for (blocked) derivative 2DY
C                             integrals accumulation (K4)
C                   ZINT2DZ = offset for (blocked) derivative 2DZ
C                             integrals accumulation (K4)
C
C                  --- Zone 3: for Stage ii) only ---
C
C                   ZPBATCH = offset for (blocked) derivated [e0|cd]
C                             or [ab|cd] batch (K4)
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

         CHARACTER*4   PRIMTYP

         INTEGER     DERAX,DERBX,DERCX,DERDX
         INTEGER     DERAY,DERBY,DERCY,DERDY
         INTEGER     DERAZ,DERBZ,DERCZ,DERDZ
         INTEGER     DERPX,DERPY,DERPZ
         INTEGER     DERQX,DERQY,DERQZ
         INTEGER     IJDIV,KLDIV
         INTEGER     ISHELLX,ISHELLY,ISHELLZ
         INTEGER     L1CACHE,NCTROW
         INTEGER     M,N
         INTEGER     MIJ,MKL,MIJKL,MRSKL
         INTEGER     MIJBIG,MKLBIG
         INTEGER     MXPRIM,MNPRIM
         INTEGER     MXSTEP
         INTEGER     NBATCH
         INTEGER     NIJ,NKL
         INTEGER     NIJBLK,NKLBLK
         INTEGER     NINT2DX,NINT2DY,NINT2DZ
         INTEGER     NGQP,NGQSCR,MGQIJKL
         INTEGER     NPGTOA,NPGTOB,NPGTOC,NPGTOD
         INTEGER     NPMINRS,NPMINTU
         INTEGER     NPSIZE,NCSIZE,NWSIZE
         INTEGER     NRS,NTU,NRSTU
         INTEGER     NXYZT
         INTEGER     PFACT
         INTEGER     POW2IJ,POW2KL
         INTEGER     PSHELL
         INTEGER     SHELLA,SHELLB,SHELLC,SHELLD
         INTEGER     SHELLP,SHELLQ
         INTEGER     SHELLX,SHELLY,SHELLZ
         INTEGER     WCOL
         INTEGER     ZONE1,ZONE2,ZONE3,ZONE3B,ZONE3C
         INTEGER     ZMAX
         INTEGER     ZCBATCH,ZPBATCH,ZWORK,
     +               ZNORMA,ZNORMB,ZNORMC,ZNORMD,
     +               ZRHOAB,ZRHOCD,
     +               ZP,ZPX,ZPY,ZPZ,ZPAX,ZPAY,ZPAZ,ZPINVHF,ZSCPK2,
     +               ZQ,ZQX,ZQY,ZQZ,ZQCX,ZQCY,ZQCZ,ZQINVHF,ZSCQK2,
     +               ZRTS,ZWTS,ZGQSCR,ZTVAL,ZPQPINV,ZSCPQK4,
     +               ZB00,ZB01,ZB10,ZC00X,ZC00Y,ZC00Z,ZD00X,ZD00Y,ZD00Z,
     +               ZEXP2A,ZEXP2B,ZEXP2C,ZEXP2D,
     +               ZINT2DX,ZINT2DY,ZINT2DZ

         DOUBLE PRECISION   DLOG2

         PARAMETER  (DLOG2  =  0.6931471805599D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...determine optimum size of [e0|cd] or [ab|cd] blocks:
C
C                The next few lines of code are among the most
C                important for performance of the entire integral
C                derivative package. Any badly chosen sizes for the
C                primitive blocks has a severe effect on the integral
C                derivative evaluation timing. Two opposite effects
C                have to be considered:
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
C                  i) The optimum [e0|cd]/[ab|cd] block size is
C                     assumed to be directly proportional to the
C                     level 1 cache size.
C
C                 ii) The proportionality factor PFACT has to be
C                     determined by 'experiment', running a real
C                     molecular case with fully contracted basis sets.
C
C                iii) After having found the optimum # of exponent
C                     quadruplets per [e0|cd]/[ab|cd] block, it is best
C                     to choose the IJ-part as large as possible,
C                     keeping the KL-part size to a minimum. The
C                     reason behind this is that a KL-part size
C                     of 1 does not trigger an intermediate IJ <-> KL
C                     transposition during the contraction procedure.
C
C                Note, that the difference between the [e0|cd] or
C                [ab|cd] block case is in the size of the # of
C                monomials quadruplets, which is given by the variable
C                NXYZT.
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
C                The routine stops, if even a MIJ x MKL block of
C                minimum size 1 x 1 cannot be accomodated into the
C                memory.
C
C
C
         IF (PRIMTYP.EQ.'E0CD') THEN

             SHELLX =   (SHELLP + 1)
     +                * (SHELLQ + DERQX + 1)
     +                * MIN ((SHELLC + DERCX + 1),(SHELLD + DERDX + 1))
             SHELLY =   (SHELLP + 1)
     +                * (SHELLQ + DERQY + 1)
     +                * MIN ((SHELLC + DERCY + 1),(SHELLD + DERDY + 1))
             SHELLZ =   (SHELLP + 1)
     +                * (SHELLQ + DERQZ + 1)
     +                * MIN ((SHELLC + DERCZ + 1),(SHELLD + DERDZ + 1))

             PSHELL = MAX (SHELLX,SHELLY,SHELLZ)

             ISHELLX =   (SHELLP + 1)
     +                 * (SHELLC + DERCX + 1)
     +                 * (SHELLD + DERDX + 1)
             ISHELLY =   (SHELLP + 1)
     +                 * (SHELLC + DERCY + 1)
     +                 * (SHELLD + DERDY + 1)
             ISHELLZ =   (SHELLP + 1)
     +                 * (SHELLC + DERCZ + 1)
     +                 * (SHELLD + DERDZ + 1)
         ELSE
             SHELLX =   (SHELLP + DERPX + 1)
     +                * (SHELLQ + DERQX + 1)
     +                * MIN ((SHELLA + DERAX + 1),(SHELLB + DERBX + 1))
     +                * MIN ((SHELLC + DERCX + 1),(SHELLD + DERDX + 1))
             SHELLY =   (SHELLP + DERPY + 1)
     +                * (SHELLQ + DERQY + 1)
     +                * MIN ((SHELLA + DERAY + 1),(SHELLB + DERBY + 1))
     +                * MIN ((SHELLC + DERCY + 1),(SHELLD + DERDY + 1))
             SHELLZ =   (SHELLP + DERPZ + 1)
     +                * (SHELLQ + DERQZ + 1)
     +                * MIN ((SHELLA + DERAZ + 1),(SHELLB + DERBZ + 1))
     +                * MIN ((SHELLC + DERCZ + 1),(SHELLD + DERDZ + 1))

             PSHELL = MAX (SHELLX,SHELLY,SHELLZ)

             ISHELLX =   (SHELLA + DERAX + 1)
     +                 * (SHELLB + DERBX + 1)
     +                 * (SHELLC + DERCX + 1)
     +                 * (SHELLD + DERDX + 1)
             ISHELLY =   (SHELLA + DERAY + 1)
     +                 * (SHELLB + DERBY + 1)
     +                 * (SHELLC + DERCY + 1)
     +                 * (SHELLD + DERDY + 1)
             ISHELLZ =   (SHELLA + DERAZ + 1)
     +                 * (SHELLB + DERBZ + 1)
     +                 * (SHELLC + DERCZ + 1)
     +                 * (SHELLD + DERDZ + 1)
         END IF

         NPMINRS = MIN (NPGTOA,NPGTOB)
         NPMINTU = MIN (NPGTOC,NPGTOD)

         NCSIZE  = NXYZT * NRSTU

         MXPRIM = MAX (NPGTOA,NPGTOB,NPGTOC,NPGTOD)
         MNPRIM = MAX (NPMINRS,NPMINTU)

         ZONE2 = NPGTOA + NPGTOB + NPGTOC + NPGTOD + NIJ + NKL
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
             NINT2DX = MGQIJKL * ISHELLX
             NINT2DY = MGQIJKL * ISHELLY
             NINT2DZ = MGQIJKL * ISHELLZ
             MRSKL   = NRS * MKL
             NBATCH  = MAX (MGQIJKL*PSHELL,NXYZT*MIJKL)
             NPSIZE  = NXYZT * MAX (MIJKL,MRSKL)
             WCOL    = MNPRIM
             NWSIZE  = NPSIZE + NCTROW * WCOL
             ZONE1   = NCSIZE
             ZONE3B  = NBATCH + NINT2DX + NINT2DY + NINT2DZ
     +                 + NGQSCR + 2*MIJKL + 9*(MIJ+MKL) + 16*MGQIJKL
             ZONE3C  = NPSIZE + NWSIZE
             ZONE3   = MAX (ZONE3B,ZONE3C)

             NKLBLK  = ZONE1 + ZONE2 + ZONE3

             MIJ     = 1
             MKL     = 1
             MIJKL   = MIJ * MKL
             MGQIJKL = NGQP * MIJKL
             NINT2DX = MGQIJKL * ISHELLX
             NINT2DY = MGQIJKL * ISHELLY
             NINT2DZ = MGQIJKL * ISHELLZ
             MRSKL   = NRS * MKL
             NBATCH  = MAX (MGQIJKL*PSHELL,NXYZT*MIJKL)
             NPSIZE  = NXYZT * MAX (MIJKL,MRSKL)
             WCOL    = MNPRIM
             NWSIZE  = NPSIZE + NCTROW * WCOL
             ZONE1   = NCSIZE
             ZONE3B  = NBATCH + NINT2DX + NINT2DY + NINT2DZ
     +                 + NGQSCR + 2*MIJKL + 9*(MIJ+MKL) + 16*MGQIJKL
             ZONE3C  = NPSIZE + NWSIZE
             ZONE3   = MAX (ZONE3B,ZONE3C)

             NIJBLK  = ZONE1 + ZONE2 + ZONE3

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

             DO 100 N = 1,MXSTEP

                MIJ = MAX (1, MIJBIG / IJDIV )
                MKL = MAX (1, MKLBIG / KLDIV )

                MIJKL   = MIJ * MKL
                MGQIJKL = NGQP * MIJKL
                NINT2DX = MGQIJKL * ISHELLX
                NINT2DY = MGQIJKL * ISHELLY
                NINT2DZ = MGQIJKL * ISHELLZ
                MRSKL   = NRS * MKL
                NBATCH  = MAX (MGQIJKL*PSHELL,NXYZT*MIJKL)
                NPSIZE  = NXYZT * MAX (MIJKL,MRSKL)
                WCOL    = MNPRIM
                NWSIZE  = NPSIZE + NCTROW * WCOL
                ZONE1   = NCSIZE
                ZONE3B  = NBATCH + NINT2DX + NINT2DY + NINT2DZ
     +                    + NGQSCR + 2*MIJKL + 9*(MIJ+MKL) + 16*MGQIJKL
                ZONE3C  = NPSIZE + NWSIZE
                ZONE3   = MAX (ZONE3B,ZONE3C)

                IF (ZONE1+ZONE2+ZONE3.LE.ZMAX) THEN

                    NIJBLK = MIJ
                    NKLBLK = MKL
                    NWSIZE = ZMAX - ZONE1 - ZONE2 - NPSIZE
C
C
C             ...generate the memory allocation pointers.
C
C
                    ZCBATCH = 1

                    ZNORMA = ZCBATCH + NCSIZE
                    ZNORMB = ZNORMA + NPGTOA
                    ZNORMC = ZNORMB + NPGTOB
                    ZNORMD = ZNORMC + NPGTOC

                    ZRHOAB = ZNORMD + NPGTOD
                    ZRHOCD = ZRHOAB + NIJ

                    ZPBATCH = ZRHOCD + NKL

                    ZP = ZPBATCH + NBATCH
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

                    ZEXP2A = ZD00Z + MGQIJKL
                    ZEXP2B = ZEXP2A + MGQIJKL
                    ZEXP2C = ZEXP2B + MGQIJKL
                    ZEXP2D = ZEXP2C + MGQIJKL

                    ZINT2DX = ZEXP2D + MGQIJKL
                    ZINT2DY = ZINT2DX + NINT2DX
                    ZINT2DZ = ZINT2DY + NINT2DY

                    ZWORK = ZPBATCH + NPSIZE

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

  100        CONTINUE

             WRITE (*,*) ' Memory allocation failed for derv blocks ! '
             WRITE (*,*) ' NIJ,NKL,MIJ,MKL = ',NIJ,NKL,MIJ,MKL
             WRITE (*,*) ' Increase flp memory! '
             WRITE (*,*) ' (erd__derv_def_blocks) '

             STOP

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
