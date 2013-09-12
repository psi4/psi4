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
         SUBROUTINE  ERD__CTR_4INDEX_BLOCK
     +
     +                    ( NPSIZE,NCSIZE,NWSIZE,
     +                      NXYZT,MIJKL,
     +                      MIJ,MKL,NRS,NTU,
     +                      NPR,NPS,NPT,NPU,
     +                      NCR,NCS,NCT,NCU,
     +                      MXPRIM,MNPRIM,
     +                      CCR,CCS,CCT,CCU,
     +                      CCBEGR,CCBEGS,CCBEGT,CCBEGU,
     +                      CCENDR,CCENDS,CCENDT,CCENDU,
     +                      PRIMR,PRIMS,PRIMT,PRIMU,
     +                      L1CACHE,TILE,NCTROW,
     +                      EQUALRS,EQUALTU,
     +                      SWAPRS,SWAPTU,
     +                      PTRANS,
     +                      BLOCKED,
     +                      PUSED,PSAVE,PPAIR,
     +                      PBATCH,
     +                      WORK,
     +
     +                                CBATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__CTR_4INDEX_BLOCK
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : ERD__TRANSPOSE_BATCH
C                ERD__CTR_1ST_HALF
C                ERD__MAP_IJKL_TO_IKJL
C                ERD__CTR_2ND_HALF_NEW
C                ERD__CTR_2ND_HALF_UPDATE
C  DESCRIPTION : This operation performs a four-indexed contraction
C                on a block of primitive integrals and updates all
C                contracted integrals:
C
C
C                  sum  [ab|cd]      -->  (ab|cd)
C                 ijkl         ijkl              rstu
C
C                                             r = 1,NCR
C                                             s = 1,NCS
C                                             t = 1,NCT
C                                             u = 1,NCU
C
C                For optimum performance of the contraction procedure,
C                the contraction routine is written to operate on the
C                block of primitive integrals ordered in the following
C                form (i.e. primitive indices to the far right):
C
C                              pbatch (nxyzt,kl,ij)
C
C                where ij and kl are the partial primitive index pairs
C                being treated here. If on entry to this routine the
C                primitive integrals are ordered in the form:
C
C                              pbatch (kl,ij,nxyzt)
C
C                a transposition has to preceed the actual contraction
C                procedure. This is triggered by the keyword PTRANS,
C                which has to be set true, if such a transposition is
C                needed. The complete set of cartesian monomial
C                quadruplets is kept together for efficiency through
C                the entire contraction sequence.
C
C                STRATEGY:
C
C                Perform overall contraction in 2 half contraction
C                steps, using an intermediate reordering if necessary:
C
C
C                1) Contract over ij primitives:
C
C                  (nxyzt,kl,rs) = ccr (r,i) * ccs (s,j) * (nxyzt,kl,ij)
C
C                2) Reorder kl <-> rs (if needed):
C
C                  (nxyzt,kl,rs) --> (nxyzt,rs,kl)
C
C                3) Contract over kl primitives and add result to the
C                   final contraction vector:
C
C                  (nxyzt,rs,tu) = (nxyzt,rs,tu)
C                                + cct (t,k) * ccu (u,l) * (nxyzt,rs,kl)
C
C
C
C                  Input:
C
C                    NxSIZE       =  size of the primitive integral
C                                    block (x=P), contracted integrals
C                                    (x=C) and working (x=W) arrays
C                    NXYZT        =  total # of cartesian monomial
C                                    quadruplets
C                    MIJKL        =  total # of partial primitive
C                                    index quadruplets
C                    MIJ(KL)      =  # of partial ij(kl) primitive
C                                    pairs to be transformed
C                    NRS(TU)      =  # of rs(tu) contraction pairs
C                    NPx          =  # of respective i,j,k,l primitives
C                                    for contractions x=R,S,T,U
C                    NCx          =  # of contractions for x=R,S,T,U
C                    MXPRIM       =  the maximum # of primitives
C                                    between all i,j,k,l primitives,
C                                    i.e. = max (i,j,k,l)
C                    MNPRIM       =  the minimum # of primitives
C                                    between i and j primitives and
C                                    k and l primitives and form the
C                                    maximum between these two values,
C                                    i.e. = max (min(i,j),min(k,l))
C                    CCx          =  full set (including zeros) of
C                                    contraction coefficients for
C                                    x=R,S,T,U contractions
C                    CCBEGx       =  lowest nonzero primitive i,j,k,l
C                                    index for x=R,S,T,U contractions
C                    CCENDx       =  highest nonzero primitive i,j,k,l
C                                    index for x=R,S,T,U contractions
C                    PRIMx        =  primitive i,j,k,l indices for
C                                    the x=R,S,T,U contractions
C                    L1CACHE      =  Size of level 1 cache in units of
C                                    8 Byte
C                    TILE         =  Number of rows and columns in
C                                    units of 8 Byte of level 1 cache
C                                    square tile array used for
C                                    performing optimum matrix
C                                    transpositions
C                    NCTROW       =  minimum # of rows that are
C                                    accepted for blocked contractions
C                    EQUALRS(TU)  =  is true, if only the lower
C                                    triangle of ij(kl) primitive
C                                    indices is present and consequently
C                                    only the lower triangle of rs(tu)
C                                    contractions needs to be evaluated
C                    SWAPRS(TU)   =  if this is true, the 1st quarter
C                                    transformation will be over R(T)
C                                    followed by the 2nd over S(U).
C                                    If false, the order will be
C                                    reversed: 1st over S(U) then 2nd
C                                    over R(T)
C                    PTRANS       =  if true, a necessary primitive
C                                    integral transposition needs to
C                                    be done in order to bring the
C                                    ijkl primitive indices to the
C                                    far right position
C                    BLOCKED      =  if false, only one call will be
C                                    made to the present contraction
C                                    routine. The contraction batch
C                                    has not! been initialized to zero
C                                    and there is no need to perform
C                                    an update of the contracted batch.
C                    Pxxxx        =  intermediate storage arrays for
C                                    primitive labels to bundle
C                                    contraction steps in do loops
C                                    (xxxx = USED,SAVE,PAIR)
C                    PBATCH       =  the batch of primitive integrals
C                                    to be contracted
C                    WORK         =  the working array for intermediate
C                                    storage
C                    CBATCH       =  the batch of contracted integrals
C                                    before contraction update
C
C                  Output:
C
C                    CBATCH       =  the update batch of contracted
C                                    integrals after contraction
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

         LOGICAL     BLOCKED
         LOGICAL     EQUALRS,EQUALTU
         LOGICAL     INWORK
         LOGICAL     PTRANS
         LOGICAL     SWAPRS,SWAPTU

         INTEGER     L1CACHE,L1FREE,L1USED
         INTEGER     MIJ,MKL,NRS,NTU
         INTEGER     MIJKL
         INTEGER     MXPRIM,MNPRIM
         INTEGER     N
         INTEGER     NBLOCK
         INTEGER     NCOL,NROW
         INTEGER     NCR,NCS,NCT,NCU
         INTEGER     NCTROW
         INTEGER     NPMAX,NPMIN
         INTEGER     NPR,NPS,NPT,NPU
         INTEGER     NPSIZE,NCSIZE,NWSIZE
         INTEGER     NXYZT
         INTEGER     TILE
         INTEGER     WROW,WUSED

         INTEGER     CCBEGR (1:NCR)
         INTEGER     CCBEGS (1:NCS)
         INTEGER     CCBEGT (1:NCT)
         INTEGER     CCBEGU (1:NCU)

         INTEGER     CCENDR (1:NCR)
         INTEGER     CCENDS (1:NCS)
         INTEGER     CCENDT (1:NCT)
         INTEGER     CCENDU (1:NCU)

         INTEGER     PRIMR  (1:MIJ)
         INTEGER     PRIMS  (1:MIJ)
         INTEGER     PRIMT  (1:MKL)
         INTEGER     PRIMU  (1:MKL)

         INTEGER     PPAIR  (1:MXPRIM)
         INTEGER     PSAVE  (1:MXPRIM)
         INTEGER     PUSED  (1:MNPRIM)

         DOUBLE PRECISION    CCR (1:NPR,1:NCR)
         DOUBLE PRECISION    CCS (1:NPS,1:NCS)
         DOUBLE PRECISION    CCT (1:NPT,1:NCT)
         DOUBLE PRECISION    CCU (1:NPU,1:NCU)

         DOUBLE PRECISION    CBATCH (1:NCSIZE)
         DOUBLE PRECISION    PBATCH (1:NPSIZE)
         DOUBLE PRECISION    WORK   (1:NWSIZE)
C
C
C------------------------------------------------------------------------
C
C
C             ...check, if a primitive batch transposition is needed:
C
C                        (kl,ij,nxyzt) --> (nxyzt,kl,ij)
C
C
         IF (PTRANS .AND. (MIJKL.GT.1) .AND. (NXYZT.GT.1)) THEN

             CALL  ERD__TRANSPOSE_BATCH
     +
     +                  ( MIJKL,NXYZT,
     +                    TILE,
     +                    PBATCH,
     +
     +                            WORK )
     +
     +
             INWORK = .TRUE.
         ELSE
             INWORK = .FALSE.
         END IF
C
C
C             ...prepare for contraction over ij. Logical variable
C                INWORK controls where the actual significant data
C                (i.e. the integrals to be contracted) is residing.
C                If INWORK is true, then array WORK contains the
C                integrals, if false, they are in array PBATCH.
C                One of the key variables to be determined here
C                is the blocking size of the invariant indices N
C                such that the cache is used efficiently. The
C                blocking size has to be adapted also to the size
C                of the working space available. L1USED contains
C                the amount of data that will occupy the cache
C                besides the three main big arrays containing the
C                initial, quarter transformed and halftransformed
C                integrals. This extra data is the contraction
C                coefficients, their segmentation limits and the
C                primitive indices. L1FREE indicates the size of
C                the level 1 cache that is finally available for
C                the three big integral arrays.
C
C
         NPMAX = MAX (NPR,NPS)
         NPMIN = MIN (NPR,NPS)

         L1USED = MIJ + NCR + NCS + NPR + NPS + NPR*NCR + NPS*NCS
         L1FREE = (3*L1CACHE/4) - L1USED

         NCOL = NPMIN + NRS + MIJ
         NROW = L1FREE / NCOL
C
C
C             ...do contraction over ij.
C
C
         N = NXYZT * MKL

         IF (INWORK) THEN

             WUSED = N * MIJ
             WROW = (NWSIZE - WUSED) / NPMIN
             NBLOCK = MIN (N,NROW,WROW)
             NBLOCK = MAX (NBLOCK,NCTROW)

             CALL  ERD__CTR_1ST_HALF
     +
     +                  ( N,
     +                    NPMAX,NPMIN,
     +                    MIJ,NRS,
     +                    NBLOCK,
     +                    NCR,NCS,
     +                    NPR,NPS,
     +                    CCR,CCS,
     +                    CCBEGR,CCBEGS,
     +                    CCENDR,CCENDS,
     +                    PRIMR,PRIMS,
     +                    EQUALRS,
     +                    SWAPRS,
     +                    PUSED,PSAVE,PPAIR,
     +                    WORK (1),WORK (WUSED+1),
     +
     +                            PBATCH )
     +
     +
             INWORK = .FALSE.

         ELSE

             WUSED = N * NRS
             WROW = (NWSIZE - WUSED) / NPMIN
             NBLOCK = MIN (N,NROW,WROW)
             NBLOCK = MAX (NBLOCK,NCTROW)

             CALL  ERD__CTR_1ST_HALF
     +
     +                  ( N,
     +                    NPMAX,NPMIN,
     +                    MIJ,NRS,
     +                    NBLOCK,
     +                    NCR,NCS,
     +                    NPR,NPS,
     +                    CCR,CCS,
     +                    CCBEGR,CCBEGS,
     +                    CCENDR,CCENDS,
     +                    PRIMR,PRIMS,
     +                    EQUALRS,
     +                    SWAPRS,
     +                    PUSED,PSAVE,PPAIR,
     +                    PBATCH,WORK (WUSED+1),
     +
     +                            WORK (1) )
     +
     +
             INWORK = .TRUE.

         END IF
C
C
C             ...reorder rs <-> kl , if needed.
C
C
         IF (MKL.GT.1 .AND. NRS.GT.1) THEN

             IF (INWORK) THEN

                 CALL  ERD__MAP_IJKL_TO_IKJL
     +
     +                  ( NXYZT,MKL,NRS,1,
     +                    TILE,
     +                    WORK,
     +
     +                            PBATCH )
     +
     +
                 INWORK = .FALSE.

             ELSE

                 CALL  ERD__MAP_IJKL_TO_IKJL
     +
     +                  ( NXYZT,MKL,NRS,1,
     +                    TILE,
     +                    PBATCH,
     +
     +                            WORK )
     +
     +
                 INWORK = .TRUE.

             END IF
         END IF
C
C
C             ...prepare for contraction over kl. Same procedure
C                as for ij (see comments above). The contracted
C                result will be added (updated) directly to the
C                final contracted CBATCH array. Do not perform a
C                contracted batch update, if no contraction blocking
C                is necessary.
C
C
         NPMAX = MAX (NPT,NPU)
         NPMIN = MIN (NPT,NPU)

         L1USED = MKL + NCT + NCU + NPT + NPU + NPT*NCT + NPU*NCU
         L1FREE = (3*L1CACHE/4) - L1USED

         NCOL = NPMIN + NTU + MKL
         NROW = L1FREE / NCOL
C
C
C             ...do contraction over kl.
C
C
         N = NXYZT * NRS

         IF (INWORK) THEN

             WUSED = N * MKL
             WROW = (NWSIZE - WUSED) / NPMIN
             NBLOCK = MIN (N,NROW,WROW)
             NBLOCK = MAX (NBLOCK,NCTROW)

             IF (BLOCKED) THEN
                 CALL  ERD__CTR_2ND_HALF_UPDATE
     +
     +                      ( N,
     +                        NPMAX,NPMIN,
     +                        MKL,NTU,
     +                        NBLOCK,
     +                        NCT,NCU,
     +                        NPT,NPU,
     +                        CCT,CCU,
     +                        CCBEGT,CCBEGU,
     +                        CCENDT,CCENDU,
     +                        PRIMT,PRIMU,
     +                        EQUALTU,
     +                        SWAPTU,
     +                        PUSED,PSAVE,PPAIR,
     +                        WORK (1),WORK (WUSED+1),
     +
     +                                 CBATCH )
     +
     +
             ELSE
                 CALL  ERD__CTR_2ND_HALF_NEW
     +
     +                      ( N,
     +                        NPMAX,NPMIN,
     +                        MKL,NTU,
     +                        NBLOCK,
     +                        NCT,NCU,
     +                        NPT,NPU,
     +                        CCT,CCU,
     +                        CCBEGT,CCBEGU,
     +                        CCENDT,CCENDU,
     +                        PRIMT,PRIMU,
     +                        EQUALTU,
     +                        SWAPTU,
     +                        PUSED,PSAVE,PPAIR,
     +                        WORK (1),WORK (WUSED+1),
     +
     +                                 CBATCH )
     +
     +
             END IF

         ELSE

             WUSED = 0
             WROW = (NWSIZE - WUSED) / NPMIN
             NBLOCK = MIN (N,NROW,WROW)
             NBLOCK = MAX (NBLOCK,NCTROW)

             IF (BLOCKED) THEN
                 CALL  ERD__CTR_2ND_HALF_UPDATE
     +
     +                      ( N,
     +                        NPMAX,NPMIN,
     +                        MKL,NTU,
     +                        NBLOCK,
     +                        NCT,NCU,
     +                        NPT,NPU,
     +                        CCT,CCU,
     +                        CCBEGT,CCBEGU,
     +                        CCENDT,CCENDU,
     +                        PRIMT,PRIMU,
     +                        EQUALTU,
     +                        SWAPTU,
     +                        PUSED,PSAVE,PPAIR,
     +                        PBATCH,WORK (1),
     +
     +                               CBATCH )
     +
     +
             ELSE
                 CALL  ERD__CTR_2ND_HALF_NEW
     +
     +                      ( N,
     +                        NPMAX,NPMIN,
     +                        MKL,NTU,
     +                        NBLOCK,
     +                        NCT,NCU,
     +                        NPT,NPU,
     +                        CCT,CCU,
     +                        CCBEGT,CCBEGU,
     +                        CCENDT,CCENDU,
     +                        PRIMT,PRIMU,
     +                        EQUALTU,
     +                        SWAPTU,
     +                        PUSED,PSAVE,PPAIR,
     +                        PBATCH,WORK (1),
     +
     +                               CBATCH )
     +
     +
             END IF

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
