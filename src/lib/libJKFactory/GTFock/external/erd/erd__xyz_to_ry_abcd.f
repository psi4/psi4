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
         SUBROUTINE  ERD__XYZ_TO_RY_ABCD
     +
     +                    ( NXYZA,NXYZB,NXYZC,NXYZD,
     +                      NRYA,NRYB,NRYC,NRYD,
     +                      SHELLA,SHELLB,SHELLC,SHELLD,
     +                      ISTART,ZSTART,
     +
     +                            NROWA,NROWB,NROWC,NROWD,
     +                            NROTA,NROTB,NROTC,NROTD,
     +                            Z00A,Z00B,Z00C,Z00D,
     +                            I0A1,I0B1,I0C1,I0D1,
     +                            I0A2,I0B2,I0C2,I0D2,
     +                            IUSED,ZUSED,
     +                            ICORE,ZCORE )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__XYZ_TO_RY_ABCD
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : ERD__XYZ_TO_RY_MATRIX
C  DESCRIPTION : This operation generates all the info needed for
C                cartesian -> spherical transformation of a batch
C                of cartesian integrals (AB|CD). No duplicate info
C                is generated, i.e. the above transmitted pointers
C                of the type Z00x, I0x1 and I0x2 with x=A,B,C,D may
C                coincide if the shell types among the A,B,C,D are
C                equal (see below). The transformation data is placed
C                in two arrays (integer and flp) at the appropriate
C                locations.
C
C                Input:
C
C                    NXYZx = dimension of xyz-monomial basis
C                            corresponding to the shell SHELLx
C                            for x=A,B,C,D.
C
C                     NRYx = dimension of ry-spherical basis
C                            corresponding to the shell SHELLx
C                            for x=A,B,C,D.
C
C                   SHELLx = the shell types for A,B,C,D.
C
C                I(Z)START = Starting location for the integer (flp)
C                            data.
C
C
C                Output:
C
C                    NROWx = maximum # of xyz-monomials contributing
C                            to the ry-components for x=A,B,C,D.
C
C                    NROTx = maximum # of elements in transformation
C                            matrix for x=A,B,C,D. This is equal to
C                            NROWx times NRYx.
C
C                     Z00x = pointer for the transformation matrix
C                            elements for x=A,B,C,D.
C
C                     I0x1 = pointer for the # of row indices leading
C                            to non-zero elements in the transformation
C                            matrix for x=A,B,C,D.
C
C                     I0x2 = pointer for the row index labels of the
C                            non-zero elements in the transformation
C                            matrix for x=A,B,C,D.
C
C                 I(Z)USED = # of integer (flp) words used.
C
C                 I(Z)CORE = The integer (flp) arrays holding the
C                            transformation data.
C
C
C                If any of the A,B,C,D shell labels are equal, their
C                offsets will be set equal in the following sequence,
C                governed by their order of usage:
C
C                    1) If shell C = D, then:
C
C                              Z00C = Z00D
C                              I0C1 = I0D1
C                              I0C2 = I0D2
C
C                    2) If shell B = C or D, then:
C
C                              Z00B = Z00C or Z00D
C                              I0B1 = I0C1 or I0D1
C                              I0B2 = I0C2 or I0D2
C
C                    3) If shell A = B or C or D, then:
C
C                              Z00A = Z00B or Z00C or Z00D
C                              I0A1 = I0B1 or I0C1 or I0D1
C                              I0A2 = I0B2 or I0C2 or I0D2
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

         LOGICAL    BDATA,CDATA

         INTEGER    I0A1,I0B1,I0C1,I0D1,I0A2,I0B2,I0C2,I0D2
         INTEGER    ISTART,ZSTART
         INTEGER    IUSED,ZUSED
         INTEGER    NROTA,NROTB,NROTC,NROTD
         INTEGER    NROWA,NROWB,NROWC,NROWD
         INTEGER    NRYA,NRYB,NRYC,NRYD
         INTEGER    NXYZA,NXYZB,NXYZC,NXYZD
         INTEGER    SHELLA,SHELLB,SHELLC,SHELLD
         INTEGER    Z00A,Z00B,Z00C,Z00D,Z0AT,Z0BT,Z0CT,Z0DT

         INTEGER    ICORE (*)

         DOUBLE PRECISION   ZCORE (*)
C
C
C------------------------------------------------------------------------
C
C
C             ...shell D data.
C
C
         IUSED = 0
         ZUSED = 0

         IF (SHELLD.GT.1) THEN
             NROWD = (SHELLD/2+1)*(SHELLD/2+2)/2
             NROTD = NROWD * NRYD
             Z00D = ZSTART
             Z0DT = Z00D + NROTD
             I0D1 = ISTART
             I0D2 = I0D1 + NRYD

             CALL  ERD__XYZ_TO_RY_MATRIX
     +
     +                  ( NXYZD,NRYD,NROWD,
     +                    SHELLD,
     +                    ZCORE (Z0DT),
     +
     +                            ICORE (I0D1),
     +                            ICORE (I0D2),
     +                            ZCORE (Z00D) )
     +
     +
             IUSED = NRYD + NROTD
             ZUSED = NROTD
         ELSE
             NROWD = 0
             NROTD = 0
             Z00D = ZSTART
             I0D2 = ISTART
         END IF
C
C
C             ...shell C data.
C
C
         CDATA = .FALSE.

         IF (SHELLC.GT.1) THEN
             IF (SHELLC.EQ.SHELLD) THEN
                 Z00C = Z00D
                 I0C1 = I0D1
                 I0C2 = I0D2
                 NROWC = NROWD
                 NROTC = NROTD
             ELSE
                 NROWC = (SHELLC/2+1)*(SHELLC/2+2)/2
                 NROTC = NROWC * NRYC
                 Z00C = Z00D + NROTD
                 Z0CT = Z00C + NROTC
                 I0C1 = I0D2 + NROTD
                 I0C2 = I0C1 + NRYC

                 CALL  ERD__XYZ_TO_RY_MATRIX
     +
     +                      ( NXYZC,NRYC,NROWC,
     +                        SHELLC,
     +                        ZCORE (Z0CT),
     +
     +                                ICORE (I0C1),
     +                                ICORE (I0C2),
     +                                ZCORE (Z00C) )
     +
     +
                 IUSED = IUSED + NRYC + NROTC
                 ZUSED = ZUSED + NROTC
                 CDATA = .TRUE.
             END IF
         ELSE
             NROWC = 0
             NROTC = 0
             Z00C = Z00D
             I0C2 = I0D2
         END IF
C
C
C             ...shell B data (being careful, using SHELLD data if
C                SHELLC data is not present!).
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
             ELSE IF (SHELLB.EQ.SHELLD) THEN
                 Z00B = Z00D
                 I0B1 = I0D1
                 I0B2 = I0D2
                 NROWB = NROWD
                 NROTB = NROTD
             ELSE
                 NROWB = (SHELLB/2+1)*(SHELLB/2+2)/2
                 NROTB = NROWB * NRYB
                 IF (CDATA) THEN
                     Z00B = Z00C + NROTC
                     Z0BT = Z00B + NROTB
                     I0B1 = I0C2 + NROTC
                     I0B2 = I0B1 + NRYB
                 ELSE
                     Z00B = Z00D + NROTD
                     Z0BT = Z00B + NROTB
                     I0B1 = I0D2 + NROTD
                     I0B2 = I0B1 + NRYB
                 END IF

                 CALL  ERD__XYZ_TO_RY_MATRIX
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
C                SHELLB data is not present or using SHELLD data if
C                also SHELLC data is not present!).
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
             ELSE IF (SHELLA.EQ.SHELLD) THEN
                 Z00A = Z00D
                 I0A1 = I0D1
                 I0A2 = I0D2
                 NROWA = NROWD
                 NROTA = NROTD
             ELSE
                 NROWA = (SHELLA/2+1)*(SHELLA/2+2)/2
                 NROTA = NROWA * NRYA
                 IF (BDATA) THEN
                     Z00A = Z00B + NROTB
                     Z0AT = Z00A + NROTA
                     I0A1 = I0B2 + NROTB
                     I0A2 = I0A1 + NRYA
                 ELSE IF (CDATA) THEN
                     Z00A = Z00C + NROTC
                     Z0AT = Z00A + NROTA
                     I0A1 = I0C2 + NROTC
                     I0A2 = I0A1 + NRYA
                 ELSE
                     Z00A = Z00D + NROTD
                     Z0AT = Z00A + NROTA
                     I0A1 = I0D2 + NROTD
                     I0A2 = I0A1 + NRYA
                 END IF

                 CALL  ERD__XYZ_TO_RY_MATRIX
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
