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
         SUBROUTINE  OED__MEMORY_OVL3C_CSGTO
     +
     +                    ( NALPHA,NCOEFF,
     +                      NCGTO1,NCGTO2,NCGTO3,
     +                      NPGTO1,NPGTO2,NPGTO3,
     +                      SHELL1,SHELL2,SHELL3,
     +                      X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,
     +                      ALPHA,CC,
     +                      L1CACHE,NCTROW,
     +                      SPHERIC,
     +
     +                                IMIN,IOPT,
     +                                ZMIN,ZOPT )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__MEMORY_OVL3C_CSGTO
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__OVL3C_SET_ABC
C                OED__OVL3C_F00_DEF_BLOCKS
C  DESCRIPTION : This operation calculates the minimum and optimum
C                integer/flp memory needed for evaluating a batch
C                of contracted 3-center overlap integrals on up to
C                three different centers between cartesian or spherical
C                gaussian type shells.
C
C
C                  Input (x = 1,2 and 3):
C
C                    NALPHA       =  total # of exponents
C                    NCOEFF       =  total # of contraction coeffs
C                    NCGTOx       =  # of contractions for csh x
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x
C                    SHELLx       =  the shell type for csh x
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = 1,2 and 3
C                    ALPHA        =  primitive exponents for csh
C                                    1,2 and 3 in that order
C                    CC           =  contraction coefficient for csh
C                                    1,2 and 3 in that order, for each
C                                    csh individually such that an
C                                    (I,J) element corresponds to the
C                                    I-th primitive and J-th contraction.
C                    L1CACHE      =  Size of level 1 cache in units of
C                                    8 Byte
C                    NCTROW       =  minimum # of rows that are
C                                    accepted for blocked contractions
C                    SPHERIC      =  is true, if spherical integrals
C                                    are wanted, false if cartesian
C                                    ones are wanted
C
C                  Output:
C
C                    IMIN,IOPT    =  minimum/optimum integer memory
C                    ZMIN,ZOPT    =  minimum/optimum flp memory
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

         LOGICAL     EMPTY
         LOGICAL     EQUALAB,EQUALAC,EQUALBC
         LOGICAL     MEMORY
         LOGICAL     SPHERIC

         LOGICAL     DUMMYL (1:6)

         INTEGER     IMIN,INEED,IOPT
         INTEGER     L1CACHE,NCTROW
         INTEGER     LCC1,LCC2,LCC3
         INTEGER     LEXP1,LEXP2,LEXP3
         INTEGER     MXPRIM,MNPRIM
         INTEGER     MXSHELL,MXSIZE
         INTEGER     NALPHA,NCOEFF
         INTEGER     NCGTO1,NCGTO2,NCGTO3
         INTEGER     NCGTOA,NCGTOB,NCGTOC
         INTEGER     NCGTOK,NCGTOIJ
         INTEGER     NCOLHRR,NROTHRR,NXYZHRR
         INTEGER     NCTR
         INTEGER     NIJ,NK,NIJK
         INTEGER     NPGTO1,NPGTO2,NPGTO3
         INTEGER     NPGTOA,NPGTOB,NPGTOC
         INTEGER     NPGTOI,NPGTOJ,NPGTOK,NPGTOIJ
         INTEGER     NROTA,NROTB,NROTC
         INTEGER     NROWA,NROWB,NROWC
         INTEGER     NRYA,NRYB,NRYC
         INTEGER     NXYZA,NXYZB,NXYZC,NXYZFT
         INTEGER     SHELL1,SHELL2,SHELL3
         INTEGER     SHELLA,SHELLB,SHELLC
         INTEGER     SHELLQ
         INTEGER     ZMIN,ZNEED,ZOPT,ZOUT1,ZOUT2

         INTEGER     DUMMYI (1:29)

         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3

         DOUBLE PRECISION  ALPHA (1:NALPHA)
         DOUBLE PRECISION  CC    (1:NCOEFF)

         DOUBLE PRECISION  DUMMYR (1:24)
C
C
C------------------------------------------------------------------------
C
C
C             ...set initial memory values.
C
C
         IMIN = 0
         IOPT = 0
         ZMIN = 0
         ZOPT = 0
C
C
C             ...fix the A,B,C labels from the 1,2,3 ones. Calculate
C                the relevant data for the A,B,C batch of 3-center
C                overlap integrals. Most of this data is not needed to
C                evaluate the memory requirements and is dumped into
C                the DUMMYx arrays with x=I,R,L standing for integer,
C                real and logical, respectively.
C
C
         LEXP1 = 1
         LEXP2 = LEXP1 + NPGTO1
         LEXP3 = LEXP2 + NPGTO2

         LCC1  = 1
         LCC2  = LCC1 + NPGTO1 * NCGTO1
         LCC3  = LCC2 + NPGTO2 * NCGTO2

         CALL  OED__OVL3C_SET_ABC
     +
     +              ( NCGTO1,NCGTO2,NCGTO3,
     +                NPGTO1,NPGTO2,NPGTO3,
     +                SHELL1,SHELL2,SHELL3,
     +                X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,
     +                ALPHA (LEXP1),ALPHA (LEXP2),ALPHA (LEXP3),
     +                CC (LCC1),CC (LCC2),CC (LCC3),
     +                SPHERIC,
     +
     +                          NCGTOA,NCGTOB,NCGTOC,
     +                          NPGTOA,NPGTOB,NPGTOC,
     +                          SHELLA,SHELLB,SHELLC,
     +                          DUMMYI (1),
     +                          SHELLQ,
     +                          MXSHELL,
     +                          DUMMYR (1), DUMMYR (2), DUMMYR (3),
     +                          DUMMYR (4), DUMMYR (5), DUMMYR (6),
     +                          DUMMYR (7), DUMMYR (8), DUMMYR (9),
     +                          DUMMYL (1),
     +                          EQUALAB,EQUALBC,EQUALAC,
     +                          DUMMYR (10),DUMMYR (11),DUMMYR (12),
     +                          DUMMYR (13),DUMMYR (14),DUMMYR (15),
     +                          DUMMYI (2), DUMMYI (3),
     +                          DUMMYR (16),DUMMYR (17),DUMMYR (18),
     +                          DUMMYR (19),
     +                          NXYZA,NXYZB,NXYZC,
     +                          DUMMYI (4), NXYZFT,
     +                          DUMMYI (5), DUMMYI (6),
     +                          NRYA,NRYB,NRYC,
     +                          DUMMYI (7), DUMMYI (8), DUMMYI (9),
     +                          DUMMYL (2), DUMMYL (3), DUMMYL (4),
     +                          DUMMYI (10),DUMMYI (11),DUMMYI (12),
     +                          DUMMYI (13),DUMMYI (14),DUMMYI (15),
     +                          DUMMYI (16),DUMMYI (17),DUMMYI (18),
     +                          NXYZHRR,NCOLHRR,NROTHRR,
     +                          EMPTY )
     +
     +
         IF (EMPTY) THEN
             RETURN
         END IF
C
C
C             ...simulate the cartesian contracted (f00) batch
C                generation. Set the ij and k primitive data. Again,
C                most of this data is not needed to evaluate the
C                memory requirements, so add it to the DUMMYx arrays.
C
C
         CALL  OED__OVL3C_SET_IJK
     +
     +              ( NCGTOA,NCGTOB,NCGTOC,
     +                NPGTOA,NPGTOB,NPGTOC,
     +                SHELLA,SHELLB,SHELLC,
     +                EQUALAB,EQUALAC,EQUALBC,
     +                DUMMYR (1), DUMMYR (2), DUMMYR (3),
     +                DUMMYR (4), DUMMYR (5), DUMMYR (6),
     +                DUMMYR (7), DUMMYR (8), DUMMYR (9),
     +                DUMMYR (10),DUMMYR (11),DUMMYR (12),
     +                DUMMYI (1), DUMMYI (2), DUMMYI (3),
     +                DUMMYI (4), DUMMYI (5), DUMMYI (6),
     +                DUMMYI (7), DUMMYI (8), DUMMYI (9),
     +                DUMMYI (10),DUMMYI (11),DUMMYI (12),
     +
     +                         DUMMYI (13),DUMMYI (14),
     +                         NCGTOK,NCGTOIJ,
     +                         NPGTOI,NPGTOJ,NPGTOK,NPGTOIJ,
     +                         DUMMYI (15),DUMMYI (16),DUMMYI (17),
     +                         DUMMYL (1), DUMMYL (2), DUMMYL (3),
     +                         DUMMYR (13),DUMMYR (14),DUMMYR (15),
     +                         DUMMYR (16),DUMMYR (17),DUMMYR (18),
     +                         DUMMYR (19),DUMMYR (20),DUMMYR (21),
     +                         DUMMYR (22),DUMMYR (23),DUMMYR (24),
     +                         DUMMYL (4), DUMMYL (5), DUMMYL (6),
     +                         DUMMYI (18),DUMMYI (19),DUMMYI (20),
     +                         DUMMYI (21),DUMMYI (22),DUMMYI (23),
     +                         DUMMYI (24),DUMMYI (25),DUMMYI (26),
     +                         DUMMYI (27),DUMMYI (28),DUMMYI (29) )
     +
     +
C
C
C             ...at this point we would determine the IJ and K
C                exponent sets necessay to evaluate the cartesian
C                contracted (f00) batch after a possible screening
C                of the primitives. Since at the moment we do not
C                apply screening for memory determination, we use the
C                complete set of IJ and K exponents. If this will be
C                changed in the future, this is where a memory routine
C                handling the screening of the IJ and K sets should
C                be placed.
C
C
         NK = NPGTOK
         NIJ = NPGTOIJ
         NIJK = NIJ * NK

         ZNEED = NIJK
         INEED = NPGTOIJ + NPGTOIJ + NPGTOK

         IMIN = MAX0 (IMIN,INEED)
         IOPT = MAX0 (IOPT,INEED)
         ZMIN = MAX0 (ZMIN,ZNEED)
         ZOPT = MAX0 (ZOPT,ZNEED)
C
C
C             ...determine minimum and optimum flp needs for the
C                unnormalized cartesian (f00) contracted 3-center
C                overlap batch generation.
C
C
         NCTR = NCGTOIJ * NCGTOK

         MEMORY = .TRUE.

         CALL  OED__OVL3C_F00_DEF_BLOCKS
     +
     +              ( 0,
     +                NPGTOI,NPGTOJ,NPGTOK,
     +                SHELLQ,
     +                NIJ,NK,NIJK,
     +                NCGTOIJ,NCGTOK,NCTR,
     +                NXYZFT,
     +                L1CACHE,NCTROW,
     +                MEMORY,
     +
     +                        DUMMYI (1),
     +                        ZOUT1,ZOUT2,
     +                        DUMMYI (2),DUMMYI (3),
     +                        MXPRIM,MNPRIM,
     +                        DUMMYI (4) ,DUMMYI (5) ,DUMMYI (6),
     +                        DUMMYI (7) ,DUMMYI (8) ,DUMMYI (9),
     +                        DUMMYI (10),DUMMYI (11),DUMMYI (12),
     +                        DUMMYI (13),DUMMYI (14),DUMMYI (15),
     +                        DUMMYI (16),DUMMYI (17),DUMMYI (18) )
     +
     +
         INEED = INEED + MXPRIM + MXPRIM + MNPRIM

         IMIN = MAX (IMIN,INEED)
         IOPT = MAX (IOPT,INEED)
         ZMIN = MAX (ZMIN,ZOUT1)
         ZOPT = MAX (ZOPT,ZOUT2)
C
C
C             ...determine the integer/flp memory needs for the next
C                steps:
C
C
C                1) expanding the contraction indices (if any)
C                2) reordering the contraction indices (if any)
C                3) transposing the contraction indices (if any)
C                4) HRR contraction
C                5) cartesian -> spherical transformation or
C                   cartesian normalization
C
C
C                The space partitioning of the flp array will be
C                as follows:
C
C
C                 |  Zone 1  |  Zone 2  |  Zone 3  |  Zone 4  |
C
C
C                 Zone 1 and 2:  2 batches of HRR maximum size
C
C                       Zone 3:  cart -> spher transformation data
C                                             or
C                                cartesian normalization factors
C
C                       Zone 4:  HRR contraction data
C
C
         NCTR = NCGTO1 * NCGTO2 * NCGTO3
         MXSIZE = NCTR * NXYZHRR
C
C
C             ...memory for Zone 1 and 2.
C
C
         INEED = 0
         ZNEED = MXSIZE + MXSIZE
C
C
C             ...memory for Zone 3.
C
C
         IF (SPHERIC) THEN
             IF (MXSHELL.GT.1) THEN

                 IF (SHELLC.GT.1) THEN
                     NROWC = (SHELLC/2+1)*(SHELLC/2+2)/2
                     NROTC = NROWC * NRYC
                     INEED = INEED + NRYC + NROTC
                     ZNEED = ZNEED + NROTC + NXYZC
                 END IF

                 IF ((SHELLB.GT.1) .AND. (SHELLB.NE.SHELLC)) THEN
                     NROWB = (SHELLB/2+1)*(SHELLB/2+2)/2
                     NROTB = NROWB * NRYB
                     INEED = INEED + NRYB + NROTB
                     ZNEED = ZNEED + NROTB + NXYZB
                 END IF

                 IF ((SHELLA.GT.1) .AND. (SHELLA.NE.SHELLB)
     +                             .AND. (SHELLA.NE.SHELLC)) THEN
                     NROWA = (SHELLA/2+1)*(SHELLA/2+2)/2
                     NROTA = NROWA * NRYA
                     INEED = INEED + NRYA + NROTA
                     ZNEED = ZNEED + NROTA + NXYZA
                 END IF

             END IF
         ELSE
             IF (MXSHELL.GT.1) THEN
                 ZNEED = ZNEED + MXSHELL + 1
             END IF
         END IF
C
C
C             ...memory for Zone 4.
C
C
         INEED = INEED + 4 * NCOLHRR + 2 * NROTHRR
         ZNEED = ZNEED + 2 * NROTHRR

         IMIN = MAX0 (IMIN,INEED)
         IOPT = MAX0 (IOPT,INEED)
         ZMIN = MAX0 (ZMIN,ZNEED)
         ZOPT = MAX0 (ZOPT,ZNEED)
C
C
C             ...ready!
C
C
         RETURN
         END
