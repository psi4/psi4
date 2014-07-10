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
         SUBROUTINE  OED__OVL3C_SET_ABC
     +
     +                    ( NCGTO1,NCGTO2,NCGTO3,
     +                      NPGTO1,NPGTO2,NPGTO3,
     +                      SHELL1,SHELL2,SHELL3,
     +                      X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,
     +                      EXP1,EXP2,EXP3,
     +                      CC1,CC2,CC3,
     +                      SPHERIC,
     +
     +                            NCGTOA,NCGTOB,NCGTOC,
     +                            NPGTOA,NPGTOB,NPGTOC,
     +                            SHELLA,SHELLB,SHELLC,SHELLP,SHELLQ,
     +                            MXSHELL,
     +                            XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,
     +                            ATOMIC,EQUALAB,EQUALBC,EQUALAC,
     +                            ABX,ABY,ABZ,ACX,ACY,ACZ,
     +                            NABCOOR,NACCOOR,
     +                            RNABSQ,RNACSQ,RNBCSQ,
     +                            SPNORM,
     +                            NXYZA,NXYZB,NXYZC,
     +                            NXYZET,NXYZFT,NXYZP,NXYZQ,
     +                            NRYA,NRYB,NRYC,
     +                            INDEXA,INDEXB,INDEXC,
     +                            SWAP12,SWAP23,SWAP13,
     +                            LEXPA,LEXPB,LEXPC,
     +                            LCCA,LCCB,LCCC,
     +                            LCCSEGA,LCCSEGB,LCCSEGC,
     +                            NXYZHRR,NCOLHRR,NROTHRR,
     +                            EMPTY )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__OVL3C_SET_ABC
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine handles the logistics on how to evaluate
C                the 3-center (123) overlap integral batch in the most
C                efficient way. It performs the label map:
C
C                                 (123) --> (ABC)
C
C                according to certain criteria that have to be met for
C                efficiency. The freedom we have in making the internal
C                association 1,2,3 -> A,B,C follows from the 6-fold
C                permutational symmetry of the 3-center overlap
C                integrals:
C
C                  (123) = (132) = (213) = (231) = (312) = (321)
C
C                Which kind of these permutations have to be applied
C                is governed by the least amount of work (# of steps)
C                for the combined two HRR procedures:
C
C                     (A+B+C,0,0) -> (A+B,C,0) -> (A,B,C)
C
C                The final A and B shell labels have thus to obey the
C                relation A>=B>=C.
C
C
C                  Input (x = 1,2 and 3):
C
C                   NCGTOx        =  # of contractions for csh x
C                   NPGTOx        =  # of primitives per contraction
C                                    for csh x
C                   SHELLx        =  the shell type for csh x
C                   Xy,Yy,Zy      =  the x,y,z-coordinates for centers
C                                    y = 1,2 and 3
C                   EXPx          =  primitive exponents for csh x
C                   CCx           =  contraction coeffs for csh x
C                   SPHERIC       =  is true, if spherical integrals
C                                    are wanted, false if cartesian
C                                    ones are wanted
C
C                  Output (x = A,B and C):
C
C                   NCGTOx        =  # of contractions for csh x
C                   NPGTOx        =  # of primitives per contraction
C                                    for csh x
C                   SHELLx        =  the shell type for csh x
C                   SHELLy        =  the shell sums y = P = A+B and
C                                    y = Q = A+B+C
C                   MXSHELL       =  the largest (maximum) shell type
C                   Xy,Yy,Zy      =  the x,y,z-coordinates for centers
C                                    y = A,B and C
C                   ATOMIC        =  indicates, if purely atomic
C                                    integrals will be evaluated
C                   EQUALxy       =  indicates, if csh x and csh y are
C                                    considered to be equal for the
C                                    pairs AB,BC and AC
C                   xyX,xyY,xyZ   =  the x,y,z-coordinate differences
C                                    between centers xy = AB and xy = AC
C                   NxyCOOR       =  # of non-zero x,y,z-coordinate
C                                    differences between centers xy = AB
C                                    and xy = AC
C                   RNxySQ        =  square of the magnitude of the
C                                    distance between centers xy = AB,AC
C                                    and BC
C                   SPNORM        =  normalization factor due to
C                                    presence of s- and p-type shells.
C                                    For each s-type shell there is a
C                                    factor of 1 and for each p-type
C                                    shell a factor of 2
C                   NXYZx         =  # of cartesian monomials for csh x
C                   NXYZET        =  sum of # of cartesian monomials
C                                    for all shells in the range
C                                    E = A,...,A+B
C                   NXYZFT        =  sum of # of cartesian monomials
C                                    for all shells in the range
C                                    F = A,...,A+B+C
C                   NXYZy         =  # of cartesian monomials for
C                                    the y = P = A+B and y = Q = A+B+C
C                                    shells
C                   NRYx          =  # of spherical functions for csh x
C                   INDEXx        =  index A,B,C -> 1,2,3 map
C                   SWAPxy        =  is .true., if a swap x <-> y has
C                                    been performed for xy = 12,23,13
C                   LEXPx         =  pointers to locate appropriate
C                                    section of the exponent array
C                                    corresponding to csh x
C                   LCCx          =  pointers to locate appropriate
C                                    section of the contraction coeff
C                                    array corresponding to csh x
C                   LCCSEGx       =  pointers to locate appropriate
C                                    section of the lowest and highest
C                                    primitive index array defining
C                                    segmented contraction boundaries
C                                    for csh x
C                   NXYZHRR       =  maximum dimension of cartesian and
C                                    spherical components during the
C                                    entire HRR contraction procedure
C                   NCOLHRR       =  maximum # of HRR rotation matrix
C                                    columns needed to generate the
C                                    final HRR rotation matrix
C                   NROTHRR       =  maximum # of HRR rotation matrix
C                                    elements needed to generate the
C                                    final HRR rotation matrix
C                   EMPTY         =  logical flag, indicating if an
C                                    empty batch of integrals is
C                                    expected.
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     ATOMIC
         LOGICAL     ATOM12,ATOM13,ATOM23
         LOGICAL     CASE1,CASE2
         LOGICAL     EMPTY
         LOGICAL     EQUAL12,EQUAL13,EQUAL23
         LOGICAL     EQUALAB,EQUALAC,EQUALBC
         LOGICAL     SPHERIC
         LOGICAL     SWAP12,SWAP23,SWAP13

         INTEGER     I,J,M
         INTEGER     INDEXA,INDEXB,INDEXC
         INTEGER     LCCSEGA,LCCSEGB,LCCSEGC
         INTEGER     LCCA,LCCB,LCCC
         INTEGER     LEXPA,LEXPB,LEXPC
         INTEGER     MXSHELL
         INTEGER     NABCOOR,NACCOOR
         INTEGER     NCGTO1,NCGTO2,NCGTO3
         INTEGER     NCGTOA,NCGTOB,NCGTOC
         INTEGER     NGH,NGHO
         INTEGER     NPGTO1,NPGTO2,NPGTO3
         INTEGER     NPGTOA,NPGTOB,NPGTOC
         INTEGER     NROW,NCOL,NROT
         INTEGER     NRYA,NRYB,NRYC
         INTEGER     NXYZA,NXYZB,NXYZC
         INTEGER     NXYZET,NXYZFT,NXYZP,NXYZQ
         INTEGER     NXYZG,NXYZH,NXYZI,NXYZGO,NXYZHO
         INTEGER     NXYZHRR,NCOLHRR,NROTHRR
         INTEGER     SHELL1,SHELL2,SHELL3
         INTEGER     SHELLA,SHELLB,SHELLC
         INTEGER     SHELLG,SHELLH
         INTEGER     SHELLP,SHELLQ

         INTEGER     ADD (0:2)

         DOUBLE PRECISION  ABX,ABY,ABZ,ACX,ACY,ACZ,BCX,BCY,BCZ
         DOUBLE PRECISION  RNABSQ,RNACSQ,RNBCSQ
         DOUBLE PRECISION  SPNORM
         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB,XC,YC,ZC
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  EXP1 (1:NPGTO1)
         DOUBLE PRECISION  EXP2 (1:NPGTO2)
         DOUBLE PRECISION  EXP3 (1:NPGTO3)

         DOUBLE PRECISION  CC1 (1:NPGTO1,1:NCGTO1)
         DOUBLE PRECISION  CC2 (1:NPGTO2,1:NCGTO2)
         DOUBLE PRECISION  CC3 (1:NPGTO3,1:NCGTO3)

         DATA  ADD    /0,0,1/

         PARAMETER  (ZERO    =  0.D0)
         PARAMETER  (ONE     =  1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...generate all 1,2,3 data. Decide as early as possible,
C                if a zero batch of overlap integrals is expected.
C
C
         EMPTY  = .FALSE.

         ATOM12 = (X1.EQ.X2) .AND. (Y1.EQ.Y2) .AND. (Z1.EQ.Z2)
         ATOM13 = (X1.EQ.X3) .AND. (Y1.EQ.Y3) .AND. (Z1.EQ.Z3)
         ATOM23 = (X2.EQ.X3) .AND. (Y2.EQ.Y3) .AND. (Z2.EQ.Z3)
         ATOMIC = ATOM12 .AND. ATOM13
         SHELLQ = SHELL1 + SHELL2 + SHELL3

         MXSHELL = MAX0 (SHELL1,SHELL2,SHELL3)

         CASE1 = MOD (SHELLQ,2) .EQ. 1
         CASE2 = SPHERIC. AND. (MXSHELL.GT.(SHELLQ-MXSHELL))

         IF (ATOMIC .AND. (CASE1.OR.CASE2)) THEN
             EMPTY = .TRUE.
             RETURN
         END IF
C
C
C             ...determine csh equality between center pairs 1,2 ,
C                1,3 and 2,3 in increasing order of complexity:
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
               DO I = 1,NPGTO1
                  EQUAL12 = EQUAL12 .AND. (EXP1(I).EQ.EXP2(I))
               END DO
               IF (EQUAL12) THEN
                 DO J = 1,NCGTO1
                    IF (EQUAL12) THEN
                      DO I = 1,NPGTO1
                         EQUAL12 = EQUAL12 .AND. (CC1(I,J).EQ.CC2(I,J))
                      END DO
                    END IF
                 END DO
               END IF
             END IF
         END IF

         EQUAL13 = ATOM13

         IF (EQUAL13) THEN
             EQUAL13 =     (SHELL1 .EQ. SHELL3)
     +               .AND. (NPGTO1 .EQ. NPGTO3)
     +               .AND. (NCGTO1 .EQ. NCGTO3)
             IF (EQUAL13) THEN
               DO I = 1,NPGTO1
                  EQUAL13 = EQUAL13 .AND. (EXP1(I).EQ.EXP3(I))
               END DO
               IF (EQUAL13) THEN
                 DO J = 1,NCGTO1
                    IF (EQUAL13) THEN
                      DO I = 1,NPGTO1
                         EQUAL13 = EQUAL13 .AND. (CC1(I,J).EQ.CC3(I,J))
                      END DO
                    END IF
                 END DO
               END IF
             END IF
         END IF

         EQUAL23 = EQUAL12 .AND. EQUAL13

         IF (.NOT.EQUAL12 .AND. .NOT.EQUAL13) THEN
           EQUAL23 = ATOM23
           IF (EQUAL23) THEN
             EQUAL23 =     (SHELL2 .EQ. SHELL3)
     +               .AND. (NPGTO2 .EQ. NPGTO3)
     +               .AND. (NCGTO2 .EQ. NCGTO3)
             IF (EQUAL23) THEN
               DO I = 1,NPGTO2
                  EQUAL23 = EQUAL23 .AND. (EXP2(I).EQ.EXP3(I))
               END DO
               IF (EQUAL23) THEN
                 DO J = 1,NCGTO2
                    IF (EQUAL23) THEN
                      DO I = 1,NPGTO2
                         EQUAL23 = EQUAL23 .AND. (CC2(I,J).EQ.CC3(I,J))
                      END DO
                    END IF
                 END DO
               END IF
             END IF
           END IF
         END IF
C
C
C             ...decide on the 1 <-> 2 , 1 <-> 3 and 2 <-> 3 swapping.
C
C
         SWAP12 = SHELL1 .LT. SHELL2
         SWAP13 = SHELL1 .LT. SHELL3
         SWAP23 = SHELL2 .LT. SHELL3
C
C
C             ...according to the previously gathered info, set the
C                new A,B,C shells, # of primitives + contraction coeffs
C                as well as pointers to the alpha exponents and
C                contraction coefficients. Also set the info for
C                evaluation of the [f00] overlap batches and for the
C                HRR steps later on.
C
C
         IF (SWAP12) THEN
             IF (SWAP23) THEN
                 XA = X3
                 YA = Y3
                 ZA = Z3
                 XB = X2
                 YB = Y2
                 ZB = Z2
                 XC = X1
                 YC = Y1
                 ZC = Z1
                 EQUALAB = EQUAL23
                 EQUALBC = EQUAL12
                 EQUALAC = EQUAL13
                 SHELLA = SHELL3
                 SHELLB = SHELL2
                 SHELLC = SHELL1
                 NPGTOA = NPGTO3
                 NPGTOB = NPGTO2
                 NPGTOC = NPGTO1
                 NCGTOA = NCGTO3
                 NCGTOB = NCGTO2
                 NCGTOC = NCGTO1
                 INDEXA = 3
                 INDEXB = 2
                 INDEXC = 1
                 LEXPC = 1
                 LEXPB = LEXPC + NPGTO1
                 LEXPA = LEXPB + NPGTO2
                 LCCC = 1
                 LCCB = LCCC + NPGTO1 * NCGTO1
                 LCCA = LCCB + NPGTO2 * NCGTO2
                 LCCSEGC = 1
                 LCCSEGB = LCCSEGC + NCGTO1
                 LCCSEGA = LCCSEGB + NCGTO2
             ELSE
                 IF (SWAP13) THEN
                     XA = X2
                     YA = Y2
                     ZA = Z2
                     XB = X3
                     YB = Y3
                     ZB = Z3
                     XC = X1
                     YC = Y1
                     ZC = Z1
                     EQUALAB = EQUAL23
                     EQUALBC = EQUAL13
                     EQUALAC = EQUAL12
                     SHELLA = SHELL2
                     SHELLB = SHELL3
                     SHELLC = SHELL1
                     NPGTOA = NPGTO2
                     NPGTOB = NPGTO3
                     NPGTOC = NPGTO1
                     NCGTOA = NCGTO2
                     NCGTOB = NCGTO3
                     NCGTOC = NCGTO1
                     INDEXA = 2
                     INDEXB = 3
                     INDEXC = 1
                     LEXPC = 1
                     LEXPA = LEXPC + NPGTO1
                     LEXPB = LEXPA + NPGTO2
                     LCCC = 1
                     LCCA = LCCC + NPGTO1 * NCGTO1
                     LCCB = LCCA + NPGTO2 * NCGTO2
                     LCCSEGC = 1
                     LCCSEGA = LCCSEGC + NCGTO1
                     LCCSEGB = LCCSEGA + NCGTO2
                 ELSE
                     XA = X2
                     YA = Y2
                     ZA = Z2
                     XB = X1
                     YB = Y1
                     ZB = Z1
                     XC = X3
                     YC = Y3
                     ZC = Z3
                     EQUALAB = EQUAL12
                     EQUALBC = EQUAL13
                     EQUALAC = EQUAL23
                     SHELLA = SHELL2
                     SHELLB = SHELL1
                     SHELLC = SHELL3
                     NPGTOA = NPGTO2
                     NPGTOB = NPGTO1
                     NPGTOC = NPGTO3
                     NCGTOA = NCGTO2
                     NCGTOB = NCGTO1
                     NCGTOC = NCGTO3
                     INDEXA = 2
                     INDEXB = 1
                     INDEXC = 3
                     LEXPB = 1
                     LEXPA = LEXPB + NPGTO1
                     LEXPC = LEXPA + NPGTO2
                     LCCB = 1
                     LCCA = LCCB + NPGTO1 * NCGTO1
                     LCCC = LCCA + NPGTO2 * NCGTO2
                     LCCSEGB = 1
                     LCCSEGA = LCCSEGB + NCGTO1
                     LCCSEGC = LCCSEGA + NCGTO2
                 END IF
             END IF
         ELSE
             IF (SWAP23) THEN
                 IF (SWAP13) THEN
                     XA = X3
                     YA = Y3
                     ZA = Z3
                     XB = X1
                     YB = Y1
                     ZB = Z1
                     XC = X2
                     YC = Y2
                     ZC = Z2
                     EQUALAB = EQUAL13
                     EQUALBC = EQUAL12
                     EQUALAC = EQUAL23
                     SHELLA = SHELL3
                     SHELLB = SHELL1
                     SHELLC = SHELL2
                     NPGTOA = NPGTO3
                     NPGTOB = NPGTO1
                     NPGTOC = NPGTO2
                     NCGTOA = NCGTO3
                     NCGTOB = NCGTO1
                     NCGTOC = NCGTO2
                     INDEXA = 3
                     INDEXB = 1
                     INDEXC = 2
                     LEXPB = 1
                     LEXPC = LEXPB + NPGTO1
                     LEXPA = LEXPC + NPGTO2
                     LCCB = 1
                     LCCC = LCCB + NPGTO1 * NCGTO1
                     LCCA = LCCC + NPGTO2 * NCGTO2
                     LCCSEGB = 1
                     LCCSEGC = LCCSEGB + NCGTO1
                     LCCSEGA = LCCSEGC + NCGTO2
                 ELSE
                     XA = X1
                     YA = Y1
                     ZA = Z1
                     XB = X3
                     YB = Y3
                     ZB = Z3
                     XC = X2
                     YC = Y2
                     ZC = Z2
                     EQUALAB = EQUAL13
                     EQUALBC = EQUAL23
                     EQUALAC = EQUAL12
                     SHELLA = SHELL1
                     SHELLB = SHELL3
                     SHELLC = SHELL2
                     NPGTOA = NPGTO1
                     NPGTOB = NPGTO3
                     NPGTOC = NPGTO2
                     NCGTOA = NCGTO1
                     NCGTOB = NCGTO3
                     NCGTOC = NCGTO2
                     INDEXA = 1
                     INDEXB = 3
                     INDEXC = 2
                     LEXPA = 1
                     LEXPC = LEXPA + NPGTO1
                     LEXPB = LEXPC + NPGTO2
                     LCCA = 1
                     LCCC = LCCA + NPGTO1 * NCGTO1
                     LCCB = LCCC + NPGTO2 * NCGTO2
                     LCCSEGA = 1
                     LCCSEGC = LCCSEGA + NCGTO1
                     LCCSEGB = LCCSEGC + NCGTO2
                 END IF
             ELSE
                 XA = X1
                 YA = Y1
                 ZA = Z1
                 XB = X2
                 YB = Y2
                 ZB = Z2
                 XC = X3
                 YC = Y3
                 ZC = Z3
                 EQUALAB = EQUAL12
                 EQUALBC = EQUAL23
                 EQUALAC = EQUAL13
                 SHELLA = SHELL1
                 SHELLB = SHELL2
                 SHELLC = SHELL3
                 NPGTOA = NPGTO1
                 NPGTOB = NPGTO2
                 NPGTOC = NPGTO3
                 NCGTOA = NCGTO1
                 NCGTOB = NCGTO2
                 NCGTOC = NCGTO3
                 INDEXA = 1
                 INDEXB = 2
                 INDEXC = 3
                 LEXPA = 1
                 LEXPB = LEXPA + NPGTO1
                 LEXPC = LEXPB + NPGTO2
                 LCCA = 1
                 LCCB = LCCA + NPGTO1 * NCGTO1
                 LCCC = LCCB + NPGTO2 * NCGTO2
                 LCCSEGA = 1
                 LCCSEGB = LCCSEGA + NCGTO1
                 LCCSEGC = LCCSEGB + NCGTO2
             END IF
         END IF
C
C
C             ...the new A,B and C shells are set. Calculate the
C                following info:
C
C                1) total shell sum A + B
C                2) cartesian monomial dimensions
C                3) cartesian monomial dimensions for shell sums
C                4) spherical dimensions (= cartesian, if no spherical)
C                5) the overall norm factor due to s- or p-type shells
C                6) the maximum size of the cartesian HRR part
C                   There will be two HRR steps involved:
C
C                           [f00] -> [e0c]      HRR I
C                           [e0c] -> [abc]      HRR II
C
C
         SHELLP = SHELLA + SHELLB

         NXYZA  = (SHELLA+1)*(SHELLA+2)/2
         NXYZB  = (SHELLB+1)*(SHELLB+2)/2
         NXYZC  = (SHELLC+1)*(SHELLC+2)/2
         NXYZP  = (SHELLP+1)*(SHELLP+2)/2
         NXYZQ  = (SHELLQ+1)*(SHELLQ+2)/2
         NXYZET =   ((SHELLP+1)*(SHELLP+2)*(SHELLP+3)/6)
     +            - ((SHELLA  )*(SHELLA+1)*(SHELLA+2)/6)
         NXYZFT =   ((SHELLQ+1)*(SHELLQ+2)*(SHELLQ+3)/6)
     +            - ((SHELLA  )*(SHELLA+1)*(SHELLA+2)/6)

         NRYA = SHELLA + SHELLA + 1
         NRYB = SHELLB + SHELLB + 1
         NRYC = SHELLC + SHELLC + 1

         IF (.NOT.SPHERIC) THEN
             NRYA = NXYZA
             NRYB = NXYZB
             NRYC = NXYZC
         END IF

         SPNORM = ONE
         IF (SHELLA.EQ.1) THEN
             SPNORM = SPNORM + SPNORM
         END IF
         IF (SHELLB.EQ.1) THEN
             SPNORM = SPNORM + SPNORM
         END IF
         IF (SHELLC.EQ.1) THEN
             SPNORM = SPNORM + SPNORM
         END IF

         NXYZHRR = MAX (NXYZET * NXYZC , NRYC * NXYZA * NXYZB)
C
C
C             ...calculate the coordinate differences between centers
C                A and B and calculate the square of the magnitude
C                of the distance between A and B. Also determine the
C                number of non-zero coordinate differences for the
C                pair A,B. Same for centers A and C. For centers
C                B and C we only need the square of the magnitude
C                of the distance between B and C.
C
C
         IF (.NOT.EQUALAB) THEN
             ABX = XA - XB
             ABY = YA - YB
             ABZ = ZA - ZB
             RNABSQ = ABX * ABX + ABY * ABY + ABZ * ABZ
             NABCOOR = 3
             IF (DABS(ABX).EQ.ZERO) NABCOOR = NABCOOR - 1
             IF (DABS(ABY).EQ.ZERO) NABCOOR = NABCOOR - 1
             IF (DABS(ABZ).EQ.ZERO) NABCOOR = NABCOOR - 1
         ELSE
             ABX = ZERO
             ABY = ZERO
             ABZ = ZERO
             RNABSQ = ZERO
             NABCOOR = 0
         END IF

         IF (.NOT.EQUALAC) THEN
             ACX = XA - XC
             ACY = YA - YC
             ACZ = ZA - ZC
             RNACSQ = ACX * ACX + ACY * ACY + ACZ * ACZ
             NACCOOR = 3
             IF (DABS(ACX).EQ.ZERO) NACCOOR = NACCOOR - 1
             IF (DABS(ACY).EQ.ZERO) NACCOOR = NACCOOR - 1
             IF (DABS(ACZ).EQ.ZERO) NACCOOR = NACCOOR - 1
         ELSE
             ACX = ZERO
             ACY = ZERO
             ACZ = ZERO
             RNACSQ = ZERO
             NACCOOR = 0
         END IF

         IF (.NOT.EQUALBC) THEN
             BCX = XB - XC
             BCY = YB - YC
             BCZ = ZB - ZC
             RNBCSQ = BCX * BCX + BCY * BCY + BCZ * BCZ
         ELSE
             RNBCSQ = ZERO
         END IF
C
C
C             ...if HRR I and/or II contractions are to be performed,
C                calculate NCOLHRR (maximum # of HRR rotation matrix
C                columns needed to generate the final HRR rotation
C                matrices) and NROTHRR (maximum # of HRR rotation
C                matrix elements).
C
C                First find maximum values for the HRR I part.
C
C
         NCOLHRR = 0
         NROTHRR = 0

         IF (SHELLC.NE.0) THEN

             NGH = NXYZFT
             NXYZG = NXYZFT
             NXYZH = 1
             NXYZGO = NXYZFT
             NXYZHO = 1
             NXYZI = NXYZQ
             SHELLG = SHELLQ
             NROW = 1
             NCOL = NGH
             NROT = NGH

             DO 100 SHELLH = 1,SHELLC
                NXYZGO = NXYZGO - NXYZI
                NXYZHO = NXYZHO + SHELLH + 1
                NGHO = NXYZGO * NXYZHO

                IF (NACCOOR.EQ.3) THEN
                    M = 1 + SHELLH/3
                    NROW = NROW +  M*(M + ADD (MOD(SHELLH,3)))
                ELSE IF (NACCOOR.EQ.2) THEN
                    NROW = NROW + SHELLH/2 + 1
                ELSE IF (NACCOOR.EQ.1) THEN
                    NROW = NROW + 1
                END IF

                NCOL = MAX0 (NGHO,NCOL)
                NROT = MAX0 (NROW*NGHO,NROT)

                NGH = NGHO
                NXYZG = NXYZGO
                NXYZH = NXYZHO
                NXYZI = NXYZI - SHELLG - 1
                SHELLG = SHELLG - 1
  100        CONTINUE

             NCOLHRR = NCOL
             NROTHRR = NROT

         END IF
C
C
C             ...next find maximum values for the HRR II part
C                and set overall maximum values.
C
C
         IF (SHELLB.NE.0) THEN

             NGH = NXYZET
             NXYZG = NXYZET
             NXYZH = 1
             NXYZGO = NXYZET
             NXYZHO = 1
             NXYZI = NXYZP
             SHELLG = SHELLP
             NROW = 1
             NCOL = NGH
             NROT = NGH

             DO 200 SHELLH = 1,SHELLB
                NXYZGO = NXYZGO - NXYZI
                NXYZHO = NXYZHO + SHELLH + 1
                NGHO = NXYZGO * NXYZHO

                IF (NABCOOR.EQ.3) THEN
                    M = 1 + SHELLH/3
                    NROW = NROW +  M*(M + ADD (MOD(SHELLH,3)))
                ELSE IF (NABCOOR.EQ.2) THEN
                    NROW = NROW + SHELLH/2 + 1
                ELSE IF (NABCOOR.EQ.1) THEN
                    NROW = NROW + 1
                END IF

                NCOL = MAX0 (NGHO,NCOL)
                NROT = MAX0 (NROW*NGHO,NROT)

                NGH = NGHO
                NXYZG = NXYZGO
                NXYZH = NXYZHO
                NXYZI = NXYZI - SHELLG - 1
                SHELLG = SHELLG - 1
  200        CONTINUE

             NCOLHRR = MAX (NCOL,NCOLHRR)
             NROTHRR = MAX (NROT,NROTHRR)

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
