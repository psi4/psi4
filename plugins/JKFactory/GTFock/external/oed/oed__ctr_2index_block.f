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
         SUBROUTINE  OED__CTR_2INDEX_BLOCK
     +
     +                    ( NPSIZE,NCSIZE,NWSIZE,
     +                      NXYZT,
     +                      MIJ,NRS,
     +                      NPR,NPS,
     +                      NCR,NCS,
     +                      MXPRIM,MNPRIM,
     +                      CCR,CCS,
     +                      CCBEGR,CCBEGS,
     +                      CCENDR,CCENDS,
     +                      PRIMR,PRIMS,
     +                      L1CACHE,TILE,NCTROW,
     +                      EQUALRS,
     +                      SWAPRS,
     +                      PTRANS,
     +                      BLOCKED,
     +                      PUSED,PSAVE,PPAIR,
     +                      PBATCH,
     +                      WORK,
     +
     +                                CBATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__CTR_2INDEX_BLOCK
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__TRANSPOSE_BATCH
C                OED__CTR_PAIR_NEW
C                OED__CTR_PAIR_UPDATE
C  DESCRIPTION : This operation performs a two-indexed contraction from
C                a block of one electron primitives to contracted
C                integrals:
C
C
C                      sum  [a|b]    -->  (a|b)
C                       ij       ij            rs
C
C                                                r = 1,NCR
C                                                s = 1,NCS
C
C                For optimum performance of the contraction procedure,
C                the contraction routine is written to operate on the
C                batch block of primitive integrals ordered in the
C                following form:
C
C                               pbatch (nxyzt,ij)
C
C                where ij are the partial primitive index pairs being
C                treated here. Their total numbers are MIJ and the
C                individual i,j indices are sitting in PRIMR,PRIMS.
C                If on entry of this routine the primitive integrals
C                are ordered in the form (ij,nxyzt), a transposition
C                has to preceed the actual contraction procedure.
C                This is triggered by the keyword PTRANS, which is
C                set true, if such a transposition is needed. The
C                complete set of cartesian monomial quadruplets is
C                kept together for efficiency through the entire
C                contraction sequence.
C
C                STRATEGY:
C
C                Perform overall contraction in one step:
C
C                  (nxyzt,rs) = (nxyzt,rs)
C                                + ccr (r,i) * ccs (s,j) * (nxyzt,ij)
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
C                    MIJ          =  # of ij primitive pairs to be
C                                    transformed
C                    NRS          =  # of rs contraction pairs
C                    NPx          =  # of respective i,j primitives
C                                    for contractions x=R,S
C                    NCx          =  # of contractions for x=R,S
C                    MXPRIM       =  the maximum # of primitives
C                                    between all i and j primitives,
C                                    i.e. = max (i,j)
C                    MNPRIM       =  the minimum # of primitives
C                                    between all i and j primitives,
C                                    i.e. = min (i,j)
C                    CCx          =  full set (including zeros) of
C                                    contraction coefficients for
C                                    x=R,S contractions
C                    CCBEGx       =  lowest nonzero primitive i,j
C                                    index for x=R,S contractions
C                    CCENDx       =  highest nonzero primitive i,j
C                                    index for x=R,S contractions
C                    PRIMx        =  primitive i,j indices for the
C                                    x=R,S contractions
C                    L1CACHE      =  Size of level 1 cache in units of
C                                    8 Byte
C                    TILE         =  Number of rows and columns in
C                                    units of 8 Byte of level 1 cache
C                                    square tile array used for
C                                    performing optimum matrix
C                                    transpositions
C                    NCTROW       =  minimum # of rows that are
C                                    accepted for blocked contractions
C                    EQUALRS      =  is true, if only the lower
C                                    triangle of ij primitive indices
C                                    is present and consequently only
C                                    the lower triangle of the rs
C                                    contractions needs to be evaluated
C                    SWAPRS       =  if this is true, the 1st half
C                                    transformation will be over R
C                                    followed by the 2nd over S.
C                                    If false, the order will be
C                                    reversed: 1st over S then 2nd
C                                    over R
C                    PTRANS       =  if true, a necessary primitive
C                                    integral transposition needs to
C                                    be done in order to bring the ij
C                                    primitive indices to the far right
C                                    position
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
         LOGICAL     EQUALRS
         LOGICAL     INWORK
         LOGICAL     PTRANS
         LOGICAL     SWAPRS

         INTEGER     L1CACHE,L1FREE,L1USED
         INTEGER     MIJ,NRS
         INTEGER     MXPRIM,MNPRIM
         INTEGER     NBLOCK
         INTEGER     NCOL,NROW
         INTEGER     NCR,NCS
         INTEGER     NCTROW
         INTEGER     NPR,NPS
         INTEGER     NPSIZE,NCSIZE,NWSIZE
         INTEGER     NXYZT
         INTEGER     TILE
         INTEGER     WROW

         INTEGER     CCBEGR (1:NCR)
         INTEGER     CCBEGS (1:NCS)

         INTEGER     CCENDR (1:NCR)
         INTEGER     CCENDS (1:NCS)

         INTEGER     PRIMR  (1:MIJ)
         INTEGER     PRIMS  (1:MIJ)

         INTEGER     PPAIR  (1:MXPRIM)
         INTEGER     PSAVE  (1:MXPRIM)
         INTEGER     PUSED  (1:MNPRIM)

         DOUBLE PRECISION    CCR (1:NPR,1:NCR)
         DOUBLE PRECISION    CCS (1:NPS,1:NCS)

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
C                        (ij,nxyzt) --> (nxyzt,ij)
C
C
         IF (PTRANS .AND. (MIJ.GT.1) .AND. (NXYZT.GT.1)) THEN

             CALL  OED__TRANSPOSE_BATCH
     +
     +                  ( MIJ,NXYZT,
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
C                is the blocking size of the invariant indices
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
         L1USED = MIJ + NCR + NCS + NPR + NPS + NPR*NCR + NPS*NCS
         L1FREE = (3*L1CACHE/4) - L1USED

         NCOL = MNPRIM + NRS + MIJ
         NROW = L1FREE / NCOL
C
C
C             ...do contraction over ij. If primitive integrals are
C                in working array, use original PBATCH array as
C                working array! Do not perform a contracted batch
C                update, if no contraction blocking is necessary.
C
C
         IF (INWORK) THEN

             WROW = NPSIZE / MNPRIM
             NBLOCK = MIN (NXYZT,NROW,WROW)
             NBLOCK = MAX (NBLOCK,NCTROW)

             IF (BLOCKED) THEN
                 CALL  OED__CTR_PAIR_UPDATE
     +
     +                      ( NXYZT,
     +                        MXPRIM,MNPRIM,
     +                        MIJ,NRS,
     +                        NBLOCK,
     +                        NCR,NCS,
     +                        NPR,NPS,
     +                        CCR,CCS,
     +                        CCBEGR,CCBEGS,
     +                        CCENDR,CCENDS,
     +                        PRIMR,PRIMS,
     +                        EQUALRS,
     +                        SWAPRS,
     +                        PUSED,PSAVE,PPAIR,
     +                        WORK,PBATCH,
     +
     +                                CBATCH )
     +
     +
             ELSE
                 CALL  OED__CTR_PAIR_NEW
     +
     +                      ( NXYZT,
     +                        MXPRIM,MNPRIM,
     +                        MIJ,NRS,
     +                        NBLOCK,
     +                        NCR,NCS,
     +                        NPR,NPS,
     +                        CCR,CCS,
     +                        CCBEGR,CCBEGS,
     +                        CCENDR,CCENDS,
     +                        PRIMR,PRIMS,
     +                        EQUALRS,
     +                        SWAPRS,
     +                        PUSED,PSAVE,PPAIR,
     +                        WORK,PBATCH,
     +
     +                                CBATCH )
     +
     +
             END IF

         ELSE

             WROW = NWSIZE / MNPRIM
             NBLOCK = MIN (NXYZT,NROW,WROW)
             NBLOCK = MAX (NBLOCK,NCTROW)

             IF (BLOCKED) THEN
                 CALL  OED__CTR_PAIR_UPDATE
     +
     +                      ( NXYZT,
     +                        MXPRIM,MNPRIM,
     +                        MIJ,NRS,
     +                        NBLOCK,
     +                        NCR,NCS,
     +                        NPR,NPS,
     +                        CCR,CCS,
     +                        CCBEGR,CCBEGS,
     +                        CCENDR,CCENDS,
     +                        PRIMR,PRIMS,
     +                        EQUALRS,
     +                        SWAPRS,
     +                        PUSED,PSAVE,PPAIR,
     +                        PBATCH,WORK,
     +
     +                                CBATCH )
     +
     +
             ELSE
                 CALL  OED__CTR_PAIR_NEW
     +
     +                      ( NXYZT,
     +                        MXPRIM,MNPRIM,
     +                        MIJ,NRS,
     +                        NBLOCK,
     +                        NCR,NCS,
     +                        NPR,NPS,
     +                        CCR,CCS,
     +                        CCBEGR,CCBEGS,
     +                        CCENDR,CCENDS,
     +                        PRIMR,PRIMS,
     +                        EQUALRS,
     +                        SWAPRS,
     +                        PUSED,PSAVE,PPAIR,
     +                        PBATCH,WORK,
     +
     +                                CBATCH )
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
