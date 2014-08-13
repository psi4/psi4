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
         SUBROUTINE  ERD__MEMORY_CSGTO
     +
     +                    ( NALPHA,NCOEFF,
     +                      NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                      NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                      SHELL1,SHELL2,SHELL3,SHELL4,
     +                      X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                      ALPHA,CC,
     +                      L1CACHE,NCTROW,
     +                      SPHERIC,
     +
     +                                IMIN,IOPT,
     +                                ZMIN,ZOPT )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__MEMORY_CSGTO
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : ERD__SET_ABCD
C                ERD__E0F0_DEF_BLOCKS
C  DESCRIPTION : This operation calculates the minimum and optimum
C                integer/flp memory needed for evaluating a batch
C                of contracted electron repulsion integrals on up to
C                four different centers between cartesian or spherical
C                gaussian type shells.
C
C
C                  Input (x = 1,2,3 and 4):
C
C                    NALPHA       =  total # of exponents
C                    NCOEFF       =  total # of contraction coeffs
C                    NCGTOx       =  # of contractions for csh x
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x
C                    SHELLx       =  the shell type for csh x
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = 1,2,3 and 4
C                    ALPHA        =  primitive exponents for csh
C                                    1,2,3,4 in that order
C                    CC           =  contraction coefficient for csh
C                                    1,2,3,4 in that order, for each
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
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     EMPTY
         LOGICAL     EQUALAB,EQUALCD
         LOGICAL     MEMORY
         LOGICAL     SPHERIC

         LOGICAL     DUMMYL (1:8)

         INTEGER     IMIN,INEED,IOPT
         INTEGER     L1CACHE,NCTROW
         INTEGER     LCC1,LCC2,LCC3,LCC4
         INTEGER     LEXP1,LEXP2,LEXP3,LEXP4
         INTEGER     MXPRIM,MNPRIM
         INTEGER     MXSHELL,MXSIZE
         INTEGER     NALPHA,NCOEFF
         INTEGER     NCGTO1,NCGTO2,NCGTO3,NCGTO4
         INTEGER     NCGTOA,NCGTOB,NCGTOC,NCGTOD,NCGTOAB,NCGTOCD
         INTEGER     NCOLHRR,NROTHRR,NXYZHRR
         INTEGER     NCTR
         INTEGER     NGQP,NMOM,NGQSCR
         INTEGER     NIJ,NKL
         INTEGER     NPGTO1,NPGTO2,NPGTO3,NPGTO4
         INTEGER     NPGTOA,NPGTOB,NPGTOC,NPGTOD,NPGTOAB,NPGTOCD
         INTEGER     NROTA,NROTB,NROTC,NROTD
         INTEGER     NROWA,NROWB,NROWC,NROWD
         INTEGER     NRYA,NRYB,NRYC,NRYD
         INTEGER     NXYZA,NXYZB,NXYZC,NXYZD
         INTEGER     NXYZT,NXYZET,NXYZFT
         INTEGER     SHELL1,SHELL2,SHELL3,SHELL4
         INTEGER     SHELLA,SHELLB,SHELLC,SHELLD
         INTEGER     SHELLP,SHELLQ,SHELLT
         INTEGER     ZMIN,ZNEED,ZOPT,ZOUT1,ZOUT2

         INTEGER     DUMMYI (1:49)

         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4

         DOUBLE PRECISION  ALPHA (1:NALPHA)
         DOUBLE PRECISION  CC    (1:NCOEFF)

         DOUBLE PRECISION  DUMMYR (1:21)
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
C             ...fix the A,B,C,D labels from the 1,2,3,4 ones.
C                Calculate the relevant data for the A,B,C,D batch of
C                integrals. Most of this data is not needed to
C                evaluate the memory requirements and is dumped into
C                the DUMMYx arrays with x=I,R,L standing for integer,
C                real and logical, respectively.
C
C
         LEXP1 = 1
         LEXP2 = LEXP1 + NPGTO1
         LEXP3 = LEXP2 + NPGTO2
         LEXP4 = LEXP3 + NPGTO3

         LCC1 = 1
         LCC2 = LCC1 + NPGTO1 * NCGTO1
         LCC3 = LCC2 + NPGTO2 * NCGTO2
         LCC4 = LCC3 + NPGTO3 * NCGTO3

         CALL  ERD__SET_ABCD
     +
     +              ( NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                SHELL1,SHELL2,SHELL3,SHELL4,
     +                X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                ALPHA (LEXP1),ALPHA (LEXP2),
     +                ALPHA (LEXP3),ALPHA (LEXP4),
     +                CC (LCC1),CC (LCC2),CC (LCC3),CC (LCC4),
     +                SPHERIC,
     +
     +                            NCGTOA,NCGTOB,NCGTOC,NCGTOD,
     +                            NPGTOA,NPGTOB,NPGTOC,NPGTOD,
     +                            SHELLA,SHELLB,SHELLC,SHELLD,
     +                            SHELLP,SHELLQ,SHELLT,
     +                            MXSHELL,
     +                            DUMMYR (1), DUMMYR (2), DUMMYR (3),
     +                            DUMMYR (4), DUMMYR (5), DUMMYR (6),
     +                            DUMMYR (7), DUMMYR (8), DUMMYR (9),
     +                            DUMMYR (10),DUMMYR (11),DUMMYR (12),
     +                            DUMMYL (1), DUMMYL (2), DUMMYL (3),
     +                            EQUALAB,EQUALCD,
     +                            DUMMYR (13),DUMMYR (14),DUMMYR (15),
     +                            DUMMYR (16),DUMMYR (17),DUMMYR (18),
     +                            DUMMYI (1), DUMMYI (2),
     +                            DUMMYR (19),DUMMYR (20),
     +                            DUMMYR (21),
     +                            NXYZA,NXYZB,NXYZC,NXYZD,
     +                            NXYZET,NXYZFT,
     +                            DUMMYI (3), DUMMYI (4),
     +                            NRYA,NRYB,NRYC,NRYD,
     +                            DUMMYI (5), DUMMYI (6),
     +                            DUMMYI (7), DUMMYI (8),
     +                            DUMMYL (4), DUMMYL (5), DUMMYL (6),
     +                            DUMMYL (7), DUMMYL (8),
     +                            DUMMYI (9), DUMMYI (10),
     +                            DUMMYI (11),DUMMYI (12),
     +                            DUMMYI (13),DUMMYI (14),
     +                            DUMMYI (15),DUMMYI (16),
     +                            DUMMYI (17),DUMMYI (18),
     +                            DUMMYI (19),DUMMYI (20),
     +                            NXYZHRR,NCOLHRR,NROTHRR,
     +                            EMPTY )
     +
     +
         IF (EMPTY) THEN
             RETURN
         END IF
C
C
C             ...simulate the cartesian contracted (e0|f0) batch
C                generation.
C
C
         IF (EQUALAB) THEN
             NPGTOAB = (NPGTOA*(NPGTOA+1))/2
             NCGTOAB = (NCGTOA*(NCGTOA+1))/2
         ELSE
             NPGTOAB = NPGTOA * NPGTOB
             NCGTOAB = NCGTOA * NCGTOB
         END IF

         IF (EQUALCD) THEN
             NPGTOCD = (NPGTOC*(NPGTOC+1))/2
             NCGTOCD = (NCGTOC*(NCGTOC+1))/2
         ELSE
             NPGTOCD = NPGTOC * NPGTOD
             NCGTOCD = NCGTOC * NCGTOD
         END IF

         NCTR = NCGTOAB * NCGTOCD
         NXYZT = NXYZET * NXYZFT
C
C
C             ...at this point we would determine the IJ and KL
C                exponent pairs necessay to evaluate the cartesian
C                contracted (e0|f0) batch after a possible screening
C                of the primitives. Since at the moment we do not
C                apply screening for memory determination, we use the
C                complete set of IJ and KL pairs. In any event this
C                will be changed in the future, this is where a memory
C                routine handling the IJ and KL pair determination
C                should be placed.
C
C
         NIJ = NPGTOAB
         NKL = NPGTOCD

         ZNEED = NPGTOAB + NPGTOCD
         INEED = 2 * ZNEED

         IMIN = MAX (IMIN,INEED)
         IOPT = MAX (IOPT,INEED)
         ZMIN = MAX (ZMIN,ZNEED)
         ZOPT = MAX (ZOPT,ZNEED)
C
C
C             ...determine minimum and optimum flp needs for the
C                unnormalized cartesian (e0|f0) contracted batch
C                generation.
C
C
         NGQP = 1 + SHELLT / 2
         NMOM = 2 * NGQP - 1
         NGQSCR = 5 * NMOM + 2 * NGQP - 2

         MEMORY = .TRUE.

         CALL  ERD__E0F0_DEF_BLOCKS
     +
     +              ( 0,
     +                NPGTOA,NPGTOB,NPGTOC,NPGTOD,
     +                SHELLP,SHELLQ,
     +                NIJ,NKL,
     +                NCGTOAB,NCGTOCD,NCTR,
     +                NGQP,NGQSCR,
     +                NXYZT,
     +                L1CACHE,NCTROW,
     +                MEMORY,
     +
     +                     ZOUT1,ZOUT2,
     +                     DUMMYI (1),DUMMYI (2),DUMMYI (3),
     +                     DUMMYI (4),
     +                     MXPRIM,MNPRIM,
     +                     DUMMYI (5) ,DUMMYI (6) ,DUMMYI (7),
     +                     DUMMYI (8) ,DUMMYI (9) ,DUMMYI (10),
     +                     DUMMYI (11),DUMMYI (12),DUMMYI (13),
     +                     DUMMYI (14),DUMMYI (15),DUMMYI (16),
     +                     DUMMYI (17),DUMMYI (18),DUMMYI (19),
     +                     DUMMYI (20),DUMMYI (21),DUMMYI (22),
     +                     DUMMYI (23),DUMMYI (24),DUMMYI (25),
     +                     DUMMYI (26),DUMMYI (27),DUMMYI (28),
     +                     DUMMYI (29),DUMMYI (30),DUMMYI (31),
     +                     DUMMYI (32),DUMMYI (33),DUMMYI (34),
     +                     DUMMYI (35),DUMMYI (36),DUMMYI (37),
     +                     DUMMYI (38),DUMMYI (39),DUMMYI (40),
     +                     DUMMYI (41),DUMMYI (42),DUMMYI (43),
     +                     DUMMYI (44),DUMMYI (45),DUMMYI (46),
     +                     DUMMYI (47),DUMMYI (48),DUMMYI (49) )
     +
     +
         INEED = INEED + 2 * MXPRIM + MNPRIM

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
         NCTR = NCGTO1 * NCGTO2 * NCGTO3 * NCGTO4
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

                 IF (SHELLD.GT.1) THEN
                     NROWD = (SHELLD/2+1)*(SHELLD/2+2)/2
                     NROTD = NROWD * NRYD
                     INEED = INEED + NRYD + NROTD
                     ZNEED = ZNEED + NROTD + NXYZD
                 END IF

                 IF ((SHELLC.GT.1) .AND. (SHELLC.NE.SHELLD)) THEN
                     NROWC = (SHELLC/2+1)*(SHELLC/2+2)/2
                     NROTC = NROWC * NRYC
                     INEED = INEED + NRYC + NROTC
                     ZNEED = ZNEED + NROTC + NXYZC
                 END IF

                 IF ((SHELLB.GT.1) .AND. (SHELLB.NE.SHELLC)
     +                             .AND. (SHELLB.NE.SHELLD)) THEN
                     NROWB = (SHELLB/2+1)*(SHELLB/2+2)/2
                     NROTB = NROWB * NRYB
                     INEED = INEED + NRYB + NROTB
                     ZNEED = ZNEED + NROTB + NXYZB
                 END IF

                 IF ((SHELLA.GT.1) .AND. (SHELLA.NE.SHELLB)
     +                             .AND. (SHELLA.NE.SHELLC)
     +                             .AND. (SHELLA.NE.SHELLD)) THEN
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
