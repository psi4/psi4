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
         SUBROUTINE  OED__XYZ_TO_RY_ABC
     +
     +                    ( NXYZA,NXYZB,NXYZC,
     +                      NRYA,NRYB,NRYC,
     +                      SHELLA,SHELLB,SHELLC,
     +                      ISTART,ZSTART,
     +
     +                            NROWA,NROWB,NROWC,
     +                            NROTA,NROTB,NROTC,
     +                            Z00A,Z00B,Z00C,
     +                            I0A1,I0B1,I0C1,
     +                            I0A2,I0B2,I0C2,
     +                            IUSED,ZUSED,
     +                            ICORE,ZCORE )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__XYZ_TO_RY_ABC
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__XYZ_TO_RY_MATRIX
C  DESCRIPTION : This operation generates all the info needed for
C                cartesian -> spherical transformation of a batch
C                of cartesian integrals (ABC). No duplicate info
C                is generated, i.e. the above transmitted pointers
C                of the type Z00x, I0x1 and I0x2 with x=A,B,C may
C                coincide if the shell types among the A,B,C are
C                equal (see below). The transformation data is placed
C                in two arrays (integer and flp) at the appropriate
C                locations.
C
C                Input:
C
C                    NXYZx = dimension of xyz-monomial basis
C                            corresponding to the shell SHELLx
C                            for x=A,B,C.
C
C                     NRYx = dimension of ry-spherical basis
C                            corresponding to the shell SHELLx
C                            for x=A,B,C.
C
C                   SHELLx = the shell types for A,B,C.
C
C                I(Z)START = Starting location for the integer (flp)
C                            data.
C
C
C                Output:
C
C                    NROWx = maximum # of xyz-monomials contributing
C                            to the ry-components for x=A,B,C.
C
C                    NROTx = maximum # of elements in transformation
C                            matrix for x=A,B,C. This is equal to
C                            NROWx times NRYx.
C
C                     Z00x = pointer for the transformation matrix
C                            elements for x=A,B,C.
C
C                     I0x1 = pointer for the # of row indices leading
C                            to non-zero elements in the transformation
C                            matrix for x=A,B,C.
C
C                     I0x2 = pointer for the row index labels of the
C                            non-zero elements in the transformation
C                            matrix for x=A,B,C.
C
C                 I(Z)USED = # of integer (flp) words used.
C
C                 I(Z)CORE = The integer (flp) arrays holding the
C                            transformation data.
C
C
C                If any of the A,B,C shell labels are equal, their
C                offsets will be set equal in the following sequence,
C                governed by their order of usage:
C
C                    1) If shell B = C, then:
C
C                              Z00B = Z00C
C                              I0B1 = I0C1
C                              I0B2 = I0C2
C
C                    2) If shell A = B or C, then:
C
C                              Z00A = Z00B or Z00C
C                              I0A1 = I0B1 or I0C1
C                              I0A2 = I0B2 or I0C2
C
C                Only mutually different transformation matrices +
C                associated data are generated.
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         LOGICAL    BDATA

         INTEGER    I0A1,I0B1,I0C1,I0A2,I0B2,I0C2
         INTEGER    ISTART,ZSTART
         INTEGER    IUSED,ZUSED
         INTEGER    NROTA,NROTB,NROTC
         INTEGER    NROWA,NROWB,NROWC
         INTEGER    NRYA,NRYB,NRYC
         INTEGER    NXYZA,NXYZB,NXYZC
         INTEGER    SHELLA,SHELLB,SHELLC
         INTEGER    Z00A,Z00B,Z00C,Z0AT,Z0BT,Z0CT

         INTEGER    ICORE (*)

         DOUBLE PRECISION   ZCORE (*)
C
C
C------------------------------------------------------------------------
C
C
C             ...shell C data.
C
C
         IUSED = 0
         ZUSED = 0

         IF (SHELLC.GT.1) THEN
             NROWC = (SHELLC/2+1)*(SHELLC/2+2)/2
             NROTC = NROWC * NRYC
             Z00C = ZSTART
             Z0CT = Z00C + NROTC
             I0C1 = ISTART
             I0C2 = I0C1 + NRYC

             CALL  OED__XYZ_TO_RY_MATRIX
     +
     +                  ( NXYZC,NRYC,NROWC,
     +                    SHELLC,
     +                    ZCORE (Z0CT),
     +
     +                            ICORE (I0C1),
     +                            ICORE (I0C2),
     +                            ZCORE (Z00C) )
     +
     +
             IUSED = NRYC + NROTC
             ZUSED = NROTC
         ELSE
             NROWC = 0
             NROTC = 0
             Z00C = ZSTART
             I0C2 = ISTART
         END IF
C
C
C             ...shell B data.
C
C
         BDATA = .FALSE.

         IF (SHELLB.GT.1) THEN
             IF (SHELLB.EQ.SHELLC) THEN
                 Z00B = Z00C
                 I0B1 = I0C1
                 I0B2 = I0C2
                 NROWB = NROWC
                 NROTB = NROTC
             ELSE
                 NROWB = (SHELLB/2+1)*(SHELLB/2+2)/2
                 NROTB = NROWB * NRYB
                 Z00B = Z00C + NROTC
                 Z0BT = Z00B + NROTB
                 I0B1 = I0C2 + NROTC
                 I0B2 = I0B1 + NRYB

                 CALL  OED__XYZ_TO_RY_MATRIX
     +
     +                      ( NXYZB,NRYB,NROWB,
     +                        SHELLB,
     +                        ZCORE (Z0BT),
     +
     +                                ICORE (I0B1),
     +                                ICORE (I0B2),
     +                                ZCORE (Z00B) )
     +
     +
                 IUSED = IUSED + NRYB + NROTB
                 ZUSED = ZUSED + NROTB
                 BDATA = .TRUE.
             END IF
         ELSE
             NROWB = 0
             NROTB = 0
             Z00B = Z00C
             I0B2 = I0C2
         END IF
C
C
C             ...shell A data (being careful, using SHELLC data if
C                SHELLB data is not present!).
C
C
         IF (SHELLA.GT.1) THEN
             IF (SHELLA.EQ.SHELLB) THEN
                 Z00A = Z00B
                 I0A1 = I0B1
                 I0A2 = I0B2
                 NROWA = NROWB
                 NROTA = NROTB
             ELSE IF (SHELLA.EQ.SHELLC) THEN
                 Z00A = Z00C
                 I0A1 = I0C1
                 I0A2 = I0C2
                 NROWA = NROWC
                 NROTA = NROTC
             ELSE
                 NROWA = (SHELLA/2+1)*(SHELLA/2+2)/2
                 NROTA = NROWA * NRYA
                 IF (BDATA) THEN
                     Z00A = Z00B + NROTB
                     Z0AT = Z00A + NROTA
                     I0A1 = I0B2 + NROTB
                     I0A2 = I0A1 + NRYA
                 ELSE
                     Z00A = Z00C + NROTC
                     Z0AT = Z00A + NROTA
                     I0A1 = I0C2 + NROTC
                     I0A2 = I0A1 + NRYA
                 END IF

                 CALL  OED__XYZ_TO_RY_MATRIX
     +
     +                      ( NXYZA,NRYA,NROWA,
     +                        SHELLA,
     +                        ZCORE (Z0AT),
     +
     +                                ICORE (I0A1),
     +                                ICORE (I0A2),
     +                                ZCORE (Z00A) )
     +
     +
                 IUSED = IUSED + NRYA + NROTA
                 ZUSED = ZUSED + NROTA
             END IF
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
