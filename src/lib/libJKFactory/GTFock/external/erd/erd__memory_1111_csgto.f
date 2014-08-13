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
         SUBROUTINE  ERD__MEMORY_1111_CSGTO
     +
     +                    ( NALPHA,NCOEFF,
     +                      NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                      NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                      SHELL1,SHELL2,SHELL3,SHELL4,
     +                      X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                      ALPHA,CC,
     +                      L1CACHE,NCTROW,
     +
     +                                IMIN,IOPT,
     +                                ZMIN,ZOPT )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__MEMORY_1111_CSGTO
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : ERD__MEMORY_1111_BLOCKS
C  DESCRIPTION : This operation calculates the minimum and optimum
C                integer/flp memory needed for evaluating a batch
C                of contracted electron repulsion integrals on up to
C                four different centers involving s- and p-type shells
C                only!
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

         LOGICAL     ATOMIC,ATOM12,ATOM34,ATOM23
         LOGICAL     EQUAL12,EQUAL34
         LOGICAL     MEMORY

         INTEGER     I,J,K,L
         INTEGER     IMIN,INEED,IOPT
         INTEGER     L1CACHE,NCTROW
         INTEGER     LCC1,LCC2,LCC3,LCC4
         INTEGER     LEXP1,LEXP2,LEXP3,LEXP4
         INTEGER     MXPRIM,MNPRIM
         INTEGER     NALPHA,NCOEFF
         INTEGER     NCGTO1,NCGTO2,NCGTO3,NCGTO4,NCGTO12,NCGTO34
         INTEGER     NCTR
         INTEGER     NIJ,NKL
         INTEGER     NPGTO1,NPGTO2,NPGTO3,NPGTO4,NPGTO12,NPGTO34
         INTEGER     NXYZ1,NXYZ2,NXYZ3,NXYZ4,NXYZT
         INTEGER     SHELL1,SHELL2,SHELL3,SHELL4,SHELLP,SHELLT
         INTEGER     ZMIN,ZNEED,ZOPT,ZOUT1,ZOUT2

         INTEGER     DUMMYI (1:22)

         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4

         DOUBLE PRECISION  ALPHA (1:NALPHA)
         DOUBLE PRECISION  CC    (1:NCOEFF)
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
C             ...simulate the cartesian contracted (12|34) batch
C                generation.
C
C
         SHELLP = SHELL1 + SHELL2
         SHELLT = SHELLP + SHELL3 + SHELL4

         ATOM12 = (X1.EQ.X2) .AND. (Y1.EQ.Y2) .AND. (Z1.EQ.Z2)
         ATOM23 = (X2.EQ.X3) .AND. (Y2.EQ.Y3) .AND. (Z2.EQ.Z3)
         ATOM34 = (X3.EQ.X4) .AND. (Y3.EQ.Y4) .AND. (Z3.EQ.Z4)

         ATOMIC = (ATOM12 .AND. ATOM34 .AND. ATOM23)

         IF (ATOMIC .AND. (MOD(SHELLT,2).EQ.1)) THEN
             RETURN
         END IF

         LEXP1 = 1
         LEXP2 = LEXP1 + NPGTO1
         LEXP3 = LEXP2 + NPGTO2
         LEXP4 = LEXP3 + NPGTO3

         LCC1 = 1
         LCC2 = LCC1 + NPGTO1 * NCGTO1
         LCC3 = LCC2 + NPGTO2 * NCGTO2
         LCC4 = LCC3 + NPGTO3 * NCGTO3
C
C
C             ...determine csh equality between center pairs 1,2
C                and 3,4 in increasing order of complexity:
C
C                 centers -> shells -> exponents -> ctr coefficients
C
C
         EQUAL12 = ATOM12

         IF (EQUAL12) THEN
             EQUAL12 =     (SHELL1 .EQ. SHELL2)
     +               .AND. (NPGTO1 .EQ. NPGTO2)
     +               .AND. (NCGTO1 .EQ. NCGTO2)
             IF (EQUAL12) THEN
               K = LEXP1 - 1
               L = LEXP2 - 1
               DO 100 I = 1,NPGTO1
                  EQUAL12 = EQUAL12 .AND. (ALPHA (K+I).EQ.ALPHA (L+I))
  100          CONTINUE
               IF (EQUAL12) THEN
                 K = LCC1 - 1
                 L = LCC2 - 1
                 DO 110 J = 1,NCGTO1
                    IF (EQUAL12) THEN
                      DO 120 I = 1,NPGTO1
                         EQUAL12 = EQUAL12 .AND. (CC (K+I).EQ.CC (L+I))
  120                 CONTINUE
                      K = K + NPGTO1
                      L = L + NPGTO1
                    END IF
  110            CONTINUE
               END IF
             END IF
         END IF

         EQUAL34 = ATOM34

         IF (EQUAL34) THEN
             EQUAL34 =     (SHELL3 .EQ. SHELL4)
     +               .AND. (NPGTO3 .EQ. NPGTO4)
     +               .AND. (NCGTO3 .EQ. NCGTO4)
             IF (EQUAL34) THEN
               K = LEXP3 - 1
               L = LEXP4 - 1
               DO 130 I = 1,NPGTO3
                  EQUAL34 = EQUAL34 .AND. (ALPHA (K+I).EQ.ALPHA (L+I))
  130          CONTINUE
               IF (EQUAL34) THEN
                 K = LCC3 - 1
                 L = LCC4 - 1
                 DO 140 J = 1,NCGTO3
                    IF (EQUAL34) THEN
                      DO 150 I = 1,NPGTO3
                         EQUAL34 = EQUAL34 .AND. (CC (K+I).EQ.CC (L+I))
  150                 CONTINUE
                      K = K + NPGTO3
                      L = L + NPGTO3
                    END IF
  140            CONTINUE
               END IF
             END IF
         END IF
C
C
C             ...calculate relevant data for the [12|34] batch of
C                integrals, such as dimensions, total # of integrals
C                to be expected, etc...
C
C
         NXYZ1  = SHELL1 + SHELL1 + 1
         NXYZ2  = SHELL2 + SHELL2 + 1
         NXYZ3  = SHELL3 + SHELL3 + 1
         NXYZ4  = SHELL4 + SHELL4 + 1

         NXYZT  = NXYZ1 * NXYZ2 * NXYZ3 * NXYZ4

         IF (EQUAL12) THEN
             NPGTO12 = (NPGTO1*(NPGTO1+1))/2
             NCGTO12 = (NCGTO1*(NCGTO1+1))/2
         ELSE
             NPGTO12 = NPGTO1 * NPGTO2
             NCGTO12 = NCGTO1 * NCGTO2
         END IF

         IF (EQUAL34) THEN
             NPGTO34 = (NPGTO3*(NPGTO3+1))/2
             NCGTO34 = (NCGTO3*(NCGTO3+1))/2
         ELSE
             NPGTO34 = NPGTO3 * NPGTO4
             NCGTO34 = NCGTO3 * NCGTO4
         END IF

         NCTR = NCGTO12 * NCGTO34
C
C
C             ...at this point we would determine the IJ and KL
C                exponent pairs necessay to evaluate the cartesian
C                contracted (12|34) batch after a possible screening
C                of the primitives. Since at the moment we do not
C                apply screening for memory determination, we use the
C                complete set of IJ and KL pairs. In any event this
C                will be changed in the future, this is where a memory
C                routine handling the IJ and KL pair determination
C                should be placed.
C
C
         NIJ = NPGTO12
         NKL = NPGTO34

         ZNEED = NPGTO12 + NPGTO34
         INEED = 2 * ZNEED

         IMIN = MAX0 (IMIN,INEED)
         IOPT = MAX0 (IOPT,INEED)
         ZMIN = MAX0 (ZMIN,ZNEED)
         ZOPT = MAX0 (ZOPT,ZNEED)
C
C
C             ...determine minimum and optimum flp needs for the
C                unnormalized cartesian (12|34) contracted batch
C                generation.
C
C
         MEMORY = .TRUE.

         CALL  ERD__1111_DEF_BLOCKS
     +
     +              ( 0,
     +                NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                NIJ,NKL,
     +                NCGTO12,NCGTO34,NCTR,
     +                NXYZT,
     +                L1CACHE,NCTROW,
     +                MEMORY,
     +
     +                          ZOUT1,ZOUT2,
     +                          DUMMYI (1),DUMMYI (2),DUMMYI (3),
     +                          MXPRIM,MNPRIM,
     +                          DUMMYI (4) ,DUMMYI (5) ,DUMMYI (6) ,
     +                          DUMMYI (7) ,DUMMYI (8) ,DUMMYI (9) ,
     +                          DUMMYI (10),DUMMYI (11),DUMMYI (12),
     +                          DUMMYI (13),DUMMYI (14),DUMMYI (15),
     +                          DUMMYI (16),DUMMYI (17),DUMMYI (18),
     +                          DUMMYI (19),DUMMYI (20),DUMMYI (21),
     +                          DUMMYI (22) )
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
C                1) expanding the contraction indices (if any)
C                2) reordering the contraction indices (if any)
C
C                The space partitioning of the flp array will be
C                as follows:
C
C
C                         |  Zone 1  |  Zone 2  |
C
C
C                 Zone 1 and 2:  2 batches of final (12|34) size
C
C
         NCTR = NCGTO1 * NCGTO2 * NCGTO3 * NCGTO4

         ZNEED = 2 * NCTR * NXYZT

         ZMIN = MAX0 (ZMIN,ZNEED)
         ZOPT = MAX0 (ZOPT,ZNEED)
C
C
C             ...ready!
C
C
         RETURN
         END
