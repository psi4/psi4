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
         SUBROUTINE  OED__MEMORY_NAI_CSGTO
     +
     +                    ( NALPHA,NCOEFF,
     +                      NCGTO1,NCGTO2,
     +                      NPGTO1,NPGTO2,
     +                      SHELL1,SHELL2,
     +                      X1,Y1,Z1,X2,Y2,Z2,
     +                      NUCLEI,
     +                      ALPHA,CC,
     +                      L1CACHE,NCTROW,
     +                      SPHERIC,
     +
     +                                IMIN,IOPT,
     +                                ZMIN,ZOPT )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__MEMORY_NAI_CSGTO
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__NAI_SET_AB
C                OED__NAI_E0_DEF_BLOCKS
C  DESCRIPTION : This operation calculates the minimum and optimum
C                integer/flp memory needed for evaluating a batch
C                of contracted electron nuclear attraction integrals
C                between cartesian or spherical gaussian type shells.
C
C
C                  Input (x = 1 and 2):
C
C                    NALPHA       =  total # of exponents
C                    NCOEFF       =  total # of contraction coeffs
C                    NCGTOx       =  # of contractions for csh x
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x
C                    SHELLx       =  the shell type for csh x
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = 1 and 2
C                    NUCLEI       =  # of nuclear attraction centers
C                    ALPHA        =  primitive exponents for csh
C                                    1 and 2 in that order
C                    CC           =  contraction coefficient for csh
C                                    1 and 2 in that order, for each
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
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     EMPTY
         LOGICAL     EQUALAB
         LOGICAL     MEMORY
         LOGICAL     SPHERIC

         LOGICAL     DUMMYL (1:3)

         INTEGER     IMIN,INEED,IOPT
         INTEGER     L1CACHE,NCTROW
         INTEGER     LCC1,LCC2
         INTEGER     LEXP1,LEXP2
         INTEGER     MXPRIM,MNPRIM
         INTEGER     MXSHELL,MXSIZE
         INTEGER     NALPHA,NCOEFF
         INTEGER     NCEN
         INTEGER     NCGTO1,NCGTO2
         INTEGER     NCGTOA,NCGTOB,NCGTOAB
         INTEGER     NCOLHRR,NROTHRR,NXYZHRR
         INTEGER     NCTR
         INTEGER     NGQP,NMOM,NGQSCR
         INTEGER     NIJ
         INTEGER     NPGTO1,NPGTO2
         INTEGER     NPGTOA,NPGTOB,NPGTOAB
         INTEGER     NROTA,NROTB
         INTEGER     NROWA,NROWB
         INTEGER     NRYA,NRYB
         INTEGER     NUCLEI
         INTEGER     NXYZA,NXYZB,NXYZET
         INTEGER     SHELL1,SHELL2
         INTEGER     SHELLA,SHELLB
         INTEGER     SHELLP
         INTEGER     ZMIN,ZNEED,ZOPT,ZOUT1,ZOUT2

         INTEGER     DUMMYI (1:32)

         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2

         DOUBLE PRECISION  ALPHA (1:NALPHA)
         DOUBLE PRECISION  CC    (1:NCOEFF)

         DOUBLE PRECISION  DUMMYR (1:11)
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
C             ...fix the A,B labels from the 1,2 ones. Calculate the
C                relevant data for the A,B batch of nuclear attraction
C                integrals. Most of this data is not needed to evaluate
C                the memory requirements and is dumped into the DUMMYx
C                arrays with x=I,R,L standing for integer, real and
C                logical, respectively.
C
C
         LEXP1 = 1
         LEXP2 = LEXP1 + NPGTO1
         LCC1 = 1
         LCC2 = LCC1 + NPGTO1 * NCGTO1

         CALL  OED__NAI_SET_AB
     +
     +              ( NCGTO1,NCGTO2,
     +                NPGTO1,NPGTO2,
     +                SHELL1,SHELL2,
     +                X1,Y1,Z1,X2,Y2,Z2,
     +                ALPHA (LEXP1),ALPHA (LEXP2),
     +                CC (LCC1),CC (LCC2),
     +                SPHERIC,
     +
     +                            NCGTOA,NCGTOB,
     +                            NPGTOA,NPGTOB,
     +                            SHELLA,SHELLB,SHELLP,
     +                            MXSHELL,
     +                            DUMMYR (1), DUMMYR (2), DUMMYR (3),
     +                            DUMMYR (4), DUMMYR (5), DUMMYR (6),
     +                            DUMMYL (1),
     +                            EQUALAB,
     +                            DUMMYR (7), DUMMYR (8), DUMMYR (9),
     +                            DUMMYI (1),
     +                            DUMMYR (10),DUMMYR (11),
     +                            NXYZA,NXYZB,NXYZET,
     +                            DUMMYI (4),
     +                            NRYA,NRYB,
     +                            DUMMYI (2), DUMMYI (3),
     +                            DUMMYL (2), DUMMYL (3),
     +                            DUMMYI (4), DUMMYI (5),
     +                            DUMMYI (6), DUMMYI (7),
     +                            DUMMYI (8), DUMMYI (9),
     +                            NXYZHRR,NCOLHRR,NROTHRR,
     +                            EMPTY )
     +
     +
         IF (EMPTY) THEN
             RETURN
         END IF
C
C
C             ...simulate the cartesian contracted (e|0) nuclear
C                attraction batch generation.
C
C
         IF (EQUALAB) THEN
             NPGTOAB = (NPGTOA*(NPGTOA+1))/2
             NCGTOAB = (NCGTOA*(NCGTOA+1))/2
         ELSE
             NPGTOAB = NPGTOA * NPGTOB
             NCGTOAB = NCGTOA * NCGTOB
         END IF
C
C
C             ...at this point we would determine the IJ exponent
C                pairs and the # of relevant nuclear attraction
C                centers necessay to evaluate the cartesian contracted
C                (e|0) nuclear attraction batch after a possible
C                screening of the primitives. Since at the moment we
C                do not apply screening for memory determination,
C                we use the complete set of IJ pairs and the complete
C                set of nuclear attraction centers. In any event
C                this will be changed in the future, this is where a
C                memory routine handling the IJ pair determination
C                should be placed.
C
C
         NIJ = NPGTOAB
         NCEN = NUCLEI

         ZNEED = NIJ + 4 * NCEN
         INEED = NPGTOAB + NPGTOAB + NUCLEI

         IMIN = MAX (IMIN,INEED)
         IOPT = MAX (IOPT,INEED)
         ZMIN = MAX (ZMIN,ZNEED)
         ZOPT = MAX (ZOPT,ZNEED)
C
C
C             ...determine minimum and optimum flp needs for the
C                unnormalized cartesian (e|0) contracted nuclear
C                attraction batch generation.
C
C
         NGQP = 1 + SHELLP / 2
         NMOM = 2 * NGQP - 1
         NGQSCR = 5 * NMOM + 2 * NGQP - 2

         MEMORY = .TRUE.

         CALL  OED__NAI_E0_DEF_BLOCKS
     +
     +              ( 0,
     +                NPGTOA,NPGTOB,
     +                SHELLP,
     +                NIJ,NCGTOAB,NCEN,
     +                NGQP,NGQSCR,
     +                NXYZET,
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
     +                        DUMMYI (16),DUMMYI (17),DUMMYI (18),
     +                        DUMMYI (19),DUMMYI (20),DUMMYI (21),
     +                        DUMMYI (22),DUMMYI (23),DUMMYI (24),
     +                        DUMMYI (25),DUMMYI (26),DUMMYI (27),
     +                        DUMMYI (28),DUMMYI (29),DUMMYI (30),
     +                        DUMMYI (31),DUMMYI (32) )
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
         NCTR = NCGTO1 * NCGTO2
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

                 IF (SHELLB.GT.1) THEN
                     NROWB = (SHELLB/2+1)*(SHELLB/2+2)/2
                     NROTB = NROWB * NRYB
                     INEED = INEED + NRYB + NROTB
                     ZNEED = ZNEED + NROTB + NXYZB
                 END IF

                 IF ((SHELLA.GT.1) .AND. (SHELLA.NE.SHELLB)) THEN
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
