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
         SUBROUTINE  ERD__SET_ABCD
     +
     +                    ( NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                      NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                      SHELL1,SHELL2,SHELL3,SHELL4,
     +                      X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                      EXP1,EXP2,EXP3,EXP4,
     +                      CC1,CC2,CC3,CC4,
     +                      SPHERIC,
     +
     +                            NCGTOA,NCGTOB,NCGTOC,NCGTOD,
     +                            NPGTOA,NPGTOB,NPGTOC,NPGTOD,
     +                            SHELLA,SHELLB,SHELLC,SHELLD,
     +                            SHELLP,SHELLQ,SHELLT,
     +                            MXSHELL,
     +                            XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     +                            ATOMIC,ATOMAB,ATOMCD,
     +                            EQUALAB,EQUALCD,
     +                            ABX,ABY,ABZ,CDX,CDY,CDZ,
     +                            NABCOOR,NCDCOOR,
     +                            RNABSQ,RNCDSQ,
     +                            SPNORM,
     +                            NXYZA,NXYZB,NXYZC,NXYZD,
     +                            NXYZET,NXYZFT,NXYZP,NXYZQ,
     +                            NRYA,NRYB,NRYC,NRYD,
     +                            INDEXA,INDEXB,INDEXC,INDEXD,
     +                            SWAP12,SWAP34,SWAPRS,SWAPTU,TR1234,
     +                            LEXPA,LEXPB,LEXPC,LEXPD,
     +                            LCCA,LCCB,LCCC,LCCD,
     +                            LCCSEGA,LCCSEGB,LCCSEGC,LCCSEGD,
     +                            NXYZHRR,NCOLHRR,NROTHRR,
     +                            EMPTY )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__SET_ABCD
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine handles the logistics on how to evaluate
C                the (12|34) integral batch in the most efficient way.
C                It performs the label map:
C
C                             (12|34) --> (AB|CD)
C
C                according to certain criteria that have to be met for
C                efficiency. The freedom we have in making the internal
C                association 1,2,3,4 -> A,B,C,D follows from the 8-fold
C                permutational symmetry of the integrals in (12|34):
C
C                     (12|34) = (21|34) = (12|43) = (21|43) =
C                     (34|12) = (43|12) = (34|21) = (43|21)
C
C                where the first line has only 12 -> 21 and 34 -> 43
C                switches and the second line is obtained from the
C                first by bra <-> ket transpositions.
C
C                The type of switch to be applied is simply governed
C                by the demand that the final A,B,C,D shell labels
C                obey the relations A>=B and C>=D, since this means
C                the least amount of work (# of steps) for the HRR
C                procedure. The decision to perform a bra <-> ket
C                transposition comes from handling memory allocations
C                during both HRR on the bra and ket sides.
C
C
C                  Input (x = 1,2,3 and 4):
C
C                   NCGTOx        =  # of contractions for csh x
C                   NPGTOx        =  # of primitives per contraction
C                                    for csh x
C                   SHELLx        =  the shell type for csh x
C                   Xy,Yy,Zy      =  the x,y,z-coordinates for centers
C                                    y = 1,2,3 and 4
C                   EXPx          =  primitive exponents for csh x
C                   CCx           =  contraction coeffs for csh x
C                   SPHERIC       =  is true, if spherical integrals
C                                    are wanted, false if cartesian
C                                    ones are wanted
C
C                  Output (x = A,B,C and D):
C
C                   NCGTOx        =  # of contractions for csh x
C                   NPGTOx        =  # of primitives per contraction
C                                    for csh x
C                   SHELLx        =  the shell type for csh x
C                   SHELLy        =  the shell sums: y = P,Q,T =
C                                    A+B,C+D,P+Q
C                   MXSHELL       =  the largest (maximum) shell type
C                   Xy,Yy,Zy      =  the x,y,z-coordinates for centers
C                                    y = A,B,C and D
C                   ATOMIC        =  indicates, if purely atomic
C                                    integrals will be evaluated
C                   ATOMxy        =  indicates, if atoms x and y are
C                                    considered to be equal for the
C                                    pairs xy = AB and CD
C                   EQUALxy       =  indicates, if csh x and csh y are
C                                    considered to be equal for the
C                                    pairs xy = AB and CD
C                   xxX,xxY,xxZ   =  the x,y,z-coordinate differences
C                                    for centers xx = AB and CD
C                   NxxCOOR       =  # of non-zero x,y,z-coordinate
C                                    differences for centers xx = AB
C                                    and CD
C                   RNxxSQ        =  square of the magnitude of the
C                                    distance between centers xx = AB
C                                    and CD
C                   SPNORM        =  normalization factor due to
C                                    presence of s- and p-type shells.
C                                    For each s-type shell there is a
C                                    factor of 1 and for each p-type
C                                    shell a factor of 2
C                   NXYZx         =  # of cartesian monomials for csh x
C                   NXYZE(F)T     =  sum of # of cartesian monomials
C                                    for all shells in the range
C                                    E = A,...,A+B and in the range
C                                    F = C,...,C+D
C                   NXYZy         =  # of cartesian monomials for
C                                    y = P,Q shells
C                   NRYx          =  # of spherical functions for csh x
C                   INDEXx        =  index A,B,C,D -> 1,2,3,4 map
C                   SWAPxy        =  is .true. for xy = 12 and 34, if
C                                    a swap 1 <-> 2 and 3 <-> 4 has
C                                    been performed
C                   SWAPRS(TU)    =  is set .true. if the contraction
C                                    order of the primitives pair AB(CD)
C                                    will be performed in reverse order
C                                    BA(DC) for efficiency reasons
C                   TR1234        =  is .true., if a bra <-> ket
C                                    transposition has been applied
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
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         LOGICAL     ATOMIC,ATOM12,ATOM34,ATOM23,ATOMAB,ATOMCD
         LOGICAL     CASE1,CASE2
         LOGICAL     EMPTY
         LOGICAL     EQUAL12,EQUAL34,EQUALAB,EQUALCD
         LOGICAL     SPHERIC
         LOGICAL     SWAP12,SWAP34,SWAPRS,SWAPTU,TR1234

         INTEGER     I,J,M
         INTEGER     INDEXA,INDEXB,INDEXC,INDEXD
         INTEGER     LCCSEGA,LCCSEGB,LCCSEGC,LCCSEGD
         INTEGER     LCCA,LCCB,LCCC,LCCD
         INTEGER     LEXPA,LEXPB,LEXPC,LEXPD
         INTEGER     MXSHELL
         INTEGER     NABCOOR,NCDCOOR
         INTEGER     NCC1,NCC2,NCC3
         INTEGER     NCGTO1,NCGTO2,NCGTO3,NCGTO4
         INTEGER     NCGTOA,NCGTOB,NCGTOC,NCGTOD
         INTEGER     NGH,NGHO
         INTEGER     NHRR1ST,NHRR2ND
         INTEGER     NPGTO1,NPGTO2,NPGTO3,NPGTO4
         INTEGER     NPGTOA,NPGTOB,NPGTOC,NPGTOD
         INTEGER     NROW,NCOL,NROT
         INTEGER     NRY1,NRY2,NRY3,NRY4
         INTEGER     NRYA,NRYB,NRYC,NRYD
         INTEGER     NXYZ1,NXYZ2,NXYZ3,NXYZ4
         INTEGER     NXYZA,NXYZB,NXYZC,NXYZD
         INTEGER     NXYZE,NXYZF,NXYZET,NXYZFT
         INTEGER     NXYZG,NXYZH,NXYZI,NXYZGO,NXYZHO
         INTEGER     NXYZP,NXYZQ
         INTEGER     NXYZHRR,NCOLHRR,NROTHRR
         INTEGER     SHELL1,SHELL2,SHELL3,SHELL4
         INTEGER     SHELLA,SHELLB,SHELLC,SHELLD
         INTEGER     SHELLG,SHELLH
         INTEGER     SHELLP,SHELLQ,SHELLT

         INTEGER     ADD (0:2)

         DOUBLE PRECISION  ABX,ABY,ABZ,CDX,CDY,CDZ
         DOUBLE PRECISION  RNABSQ,RNCDSQ
         DOUBLE PRECISION  SPNORM
         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  EXP1 (1:NPGTO1)
         DOUBLE PRECISION  EXP2 (1:NPGTO2)
         DOUBLE PRECISION  EXP3 (1:NPGTO3)
         DOUBLE PRECISION  EXP4 (1:NPGTO4)

         DOUBLE PRECISION  CC1 (1:NPGTO1,1:NCGTO1)
         DOUBLE PRECISION  CC2 (1:NPGTO2,1:NCGTO2)
         DOUBLE PRECISION  CC3 (1:NPGTO3,1:NCGTO3)
         DOUBLE PRECISION  CC4 (1:NPGTO4,1:NCGTO4)

         DATA  ADD    /0,0,1/

         PARAMETER  (ZERO    =  0.D0)
         PARAMETER  (ONE     =  1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...generate all 1,2,3,4 data. Decide as early as
C                possible, if a zero batch of integrals is expected.
C
C
         EMPTY = .FALSE.

         ATOM12 = (X1.EQ.X2) .AND. (Y1.EQ.Y2) .AND. (Z1.EQ.Z2)
         ATOM23 = (X2.EQ.X3) .AND. (Y2.EQ.Y3) .AND. (Z2.EQ.Z3)
         ATOM34 = (X3.EQ.X4) .AND. (Y3.EQ.Y4) .AND. (Z3.EQ.Z4)

         ATOMIC = (ATOM12 .AND. ATOM34 .AND. ATOM23)

         SHELLP = SHELL1 + SHELL2
         SHELLQ = SHELL3 + SHELL4
         SHELLT = SHELLP + SHELLQ

         MXSHELL = MAX (SHELL1,SHELL2,SHELL3,SHELL4)

         CASE1 = MOD (SHELLT,2) .EQ. 1
         CASE2 = SPHERIC. AND. ((MXSHELL+MXSHELL).GT.SHELLT)

         IF (ATOMIC .AND. (CASE1.OR.CASE2)) THEN
             EMPTY = .TRUE.
             RETURN
         END IF
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

         EQUAL34 = ATOM34

         IF (EQUAL34) THEN
             EQUAL34 =     (SHELL3 .EQ. SHELL4)
     +               .AND. (NPGTO3 .EQ. NPGTO4)
     +               .AND. (NCGTO3 .EQ. NCGTO4)
             IF (EQUAL34) THEN
               DO I = 1,NPGTO3
                  EQUAL34 = EQUAL34 .AND. (EXP3(I).EQ.EXP4(I))
               END DO
               IF (EQUAL34) THEN
                 DO J = 1,NCGTO3
                    IF (EQUAL34) THEN
                      DO I = 1,NPGTO3
                         EQUAL34 = EQUAL34 .AND. (CC3(I,J).EQ.CC4(I,J))
                      END DO
                    END IF
                 END DO
               END IF
             END IF
         END IF
C
C
C             ...set the cartesian and spherical dimensions. In case
C                no spherical transformations are wanted, set the
C                corresponding dimensions equal to the cartesian ones.
C
C
         NXYZ1  = (SHELL1+1)*(SHELL1+2)/2
         NXYZ2  = (SHELL2+1)*(SHELL2+2)/2
         NXYZ3  = (SHELL3+1)*(SHELL3+2)/2
         NXYZ4  = (SHELL4+1)*(SHELL4+2)/2

         NRY1 = SHELL1 + SHELL1 + 1
         NRY2 = SHELL2 + SHELL2 + 1
         NRY3 = SHELL3 + SHELL3 + 1
         NRY4 = SHELL4 + SHELL4 + 1

         IF (.NOT.SPHERIC) THEN
             NRY1 = NXYZ1
             NRY2 = NXYZ2
             NRY3 = NXYZ3
             NRY4 = NXYZ4
         END IF
C
C
C             ...decide on the 1 <-> 2 and/or 3 <-> 4 swapping.
C
C
         SWAP12 = SHELL1.LT.SHELL2
         SWAP34 = SHELL3.LT.SHELL4
C
C
C             ...calculate NXYZHRR for the two possible HRR
C                and (if any) cartesian -> spherical transformation
C                sequences:
C
C                    i) initial dimension:  NXYZE * NXYZF
C                   ii) perform HRR on 34:  NXYZE * NXYZ3 * NXYZ4
C                  iii) cart -> sph on 34:  NXYZE * NRY3 * NRY4
C                   iv) perform HRR on 12:  NXYZ1 * NXYZ2 * NRY3 * NRY4
C                    v) cart -> sph on 12:  NRY1 * NRY2 * NRY3 * NRY4
C
C
C                    i) initial dimension:  NXYZE * NXYZF
C                   ii) perform HRR on 12:  NXYZ1 * NXYZ2 * NXYZF
C                  iii) cart -> sph on 12:  NRY1 * NRY2 * NXYZF
C                   iv) perform HRR on 34:  NRY1 * NRY2 * NXYZ3 * NXYZ4
C                    v) cart -> sph on 34:  NRY1 * NRY2 * NRY3 * NRY4
C
C
C                The only dimension increasing steps are ii) and iv)
C                in both cases. Hence we first find the maximum
C                between ii) and iv) for both sequences and then we
C                take the overall minimum of these two maxima.
C                Since the order of sequence of the HRR on the A,B,C,D
C                labels is CD followed by AB, the overall minimum
C                will decide if to perform a bra <-> ket transposition
C                on the 12/34 labels.
C
C
         SHELLA = MAX0 (SHELL1,SHELL2)
         SHELLB = MIN0 (SHELL1,SHELL2)
         SHELLC = MAX0 (SHELL3,SHELL4)
         SHELLD = MIN0 (SHELL3,SHELL4)

         NXYZE  =   ((SHELLP+1)*(SHELLP+2)*(SHELLP+3)/6)
     +            - ((SHELLA  )*(SHELLA+1)*(SHELLA+2)/6)
         NXYZF  =   ((SHELLQ+1)*(SHELLQ+2)*(SHELLQ+3)/6)
     +            - ((SHELLC  )*(SHELLC+1)*(SHELLC+2)/6)

         IF (SHELLB.EQ.0 .AND. SHELLD.EQ.0) THEN
             NXYZHRR = NXYZE * NXYZF
             TR1234 = .FALSE.
         ELSE
             NHRR1ST = MAX0 (NXYZE*NXYZ3*NXYZ4,NXYZ1*NXYZ2*NRY3*NRY4)
             NHRR2ND = MAX0 (NXYZF*NXYZ1*NXYZ2,NXYZ3*NXYZ4*NRY1*NRY2)
             NXYZHRR = MIN0 (NHRR1ST,NHRR2ND)
             TR1234  = NHRR1ST .GT. NHRR2ND
         END IF
C
C
C             ...according to the previously gathered info, set the
C                new A,B,C,D shells, # of primitives + contraction
C                coeffs as well as pointers to the alpha exponents
C                and contraction coefficients. Also set the info for
C                evaluation of the [e0|f0] batches and for the HRR
C                steps later on.
C
C
         NCC1 = NPGTO1 * NCGTO1
         NCC2 = NPGTO2 * NCGTO2
         NCC3 = NPGTO3 * NCGTO3

         IF (.NOT.TR1234) THEN

             NXYZET = NXYZE
             NXYZFT = NXYZF
             ATOMAB = ATOM12
             ATOMCD = ATOM34
             EQUALAB = EQUAL12
             EQUALCD = EQUAL34

             IF (.NOT.SWAP12) THEN
                 XA = X1
                 YA = Y1
                 ZA = Z1
                 XB = X2
                 YB = Y2
                 ZB = Z2
                 SHELLA = SHELL1
                 SHELLB = SHELL2
                 NPGTOA = NPGTO1
                 NPGTOB = NPGTO2
                 NCGTOA = NCGTO1
                 NCGTOB = NCGTO2
                 NXYZA = NXYZ1
                 NXYZB = NXYZ2
                 NRYA = NRY1
                 NRYB = NRY2
                 INDEXA = 1
                 INDEXB = 2
                 LEXPA = 1
                 LEXPB = LEXPA + NPGTO1
                 LCCA = 1
                 LCCB = LCCA + NCC1
                 LCCSEGA = 1
                 LCCSEGB = LCCSEGA + NCGTO1
             ELSE
                 XA = X2
                 YA = Y2
                 ZA = Z2
                 XB = X1
                 YB = Y1
                 ZB = Z1
                 SHELLA = SHELL2
                 SHELLB = SHELL1
                 NPGTOA = NPGTO2
                 NPGTOB = NPGTO1
                 NCGTOA = NCGTO2
                 NCGTOB = NCGTO1
                 NXYZA = NXYZ2
                 NXYZB = NXYZ1
                 NRYA = NRY2
                 NRYB = NRY1
                 INDEXA = 2
                 INDEXB = 1
                 LEXPB = 1
                 LEXPA = LEXPB + NPGTO1
                 LCCB = 1
                 LCCA = LCCB + NCC1
                 LCCSEGB = 1
                 LCCSEGA = LCCSEGB + NCGTO1
             END IF

             IF (.NOT.SWAP34) THEN
                 XC = X3
                 YC = Y3
                 ZC = Z3
                 XD = X4
                 YD = Y4
                 ZD = Z4
                 SHELLC = SHELL3
                 SHELLD = SHELL4
                 NPGTOC = NPGTO3
                 NPGTOD = NPGTO4
                 NCGTOC = NCGTO3
                 NCGTOD = NCGTO4
                 NXYZC = NXYZ3
                 NXYZD = NXYZ4
                 NRYC = NRY3
                 NRYD = NRY4
                 INDEXC = 3
                 INDEXD = 4
                 LEXPC = 1 + NPGTO1 + NPGTO2
                 LEXPD = LEXPC + NPGTO3
                 LCCC = 1 + NCC1 + NCC2
                 LCCD = LCCC + NCC3
                 LCCSEGC = 1 + NCGTO1 + NCGTO2
                 LCCSEGD = LCCSEGC + NCGTO3
             ELSE
                 XC = X4
                 YC = Y4
                 ZC = Z4
                 XD = X3
                 YD = Y3
                 ZD = Z3
                 SHELLC = SHELL4
                 SHELLD = SHELL3
                 NPGTOC = NPGTO4
                 NPGTOD = NPGTO3
                 NCGTOC = NCGTO4
                 NCGTOD = NCGTO3
                 NXYZC = NXYZ4
                 NXYZD = NXYZ3
                 NRYC = NRY4
                 NRYD = NRY3
                 INDEXC = 4
                 INDEXD = 3
                 LEXPD = 1 + NPGTO1 + NPGTO2
                 LEXPC = LEXPD + NPGTO3
                 LCCD = 1 + NCC1 + NCC2
                 LCCC = LCCD + NCC3
                 LCCSEGD = 1 + NCGTO1 + NCGTO2
                 LCCSEGC = LCCSEGD + NCGTO3
             END IF

         ELSE

             NXYZET = NXYZF
             NXYZFT = NXYZE
             ATOMAB = ATOM34
             ATOMCD = ATOM12
             EQUALAB = EQUAL34
             EQUALCD = EQUAL12

             IF (.NOT.SWAP12) THEN
                 XC = X1
                 YC = Y1
                 ZC = Z1
                 XD = X2
                 YD = Y2
                 ZD = Z2
                 SHELLC = SHELL1
                 SHELLD = SHELL2
                 NPGTOC = NPGTO1
                 NPGTOD = NPGTO2
                 NCGTOC = NCGTO1
                 NCGTOD = NCGTO2
                 NXYZC = NXYZ1
                 NXYZD = NXYZ2
                 NRYC = NRY1
                 NRYD = NRY2
                 INDEXC = 1
                 INDEXD = 2
                 LEXPC = 1
                 LEXPD = LEXPC + NPGTO1
                 LCCC = 1
                 LCCD = LCCC + NCC1
                 LCCSEGC = 1
                 LCCSEGD = LCCSEGC + NCGTO1
             ELSE
                 XC = X2
                 YC = Y2
                 ZC = Z2
                 XD = X1
                 YD = Y1
                 ZD = Z1
                 SHELLC = SHELL2
                 SHELLD = SHELL1
                 NPGTOC = NPGTO2
                 NPGTOD = NPGTO1
                 NCGTOC = NCGTO2
                 NCGTOD = NCGTO1
                 NXYZC = NXYZ2
                 NXYZD = NXYZ1
                 NRYC = NRY2
                 NRYD = NRY1
                 INDEXC = 2
                 INDEXD = 1
                 LEXPD = 1
                 LEXPC = LEXPD + NPGTO1
                 LCCD = 1
                 LCCC = LCCD + NCC1
                 LCCSEGD = 1
                 LCCSEGC = LCCSEGD + NCGTO1
             END IF

             IF (.NOT.SWAP34) THEN
                 XA = X3
                 YA = Y3
                 ZA = Z3
                 XB = X4
                 YB = Y4
                 ZB = Z4
                 SHELLA = SHELL3
                 SHELLB = SHELL4
                 NPGTOA = NPGTO3
                 NPGTOB = NPGTO4
                 NCGTOA = NCGTO3
                 NCGTOB = NCGTO4
                 NXYZA = NXYZ3
                 NXYZB = NXYZ4
                 NRYA = NRY3
                 NRYB = NRY4
                 INDEXA = 3
                 INDEXB = 4
                 LEXPA = 1 + NPGTO1 + NPGTO2
                 LEXPB = LEXPA + NPGTO3
                 LCCA = 1 + NCC1 + NCC2
                 LCCB = LCCA + NCC3
                 LCCSEGA = 1 + NCGTO1 + NCGTO2
                 LCCSEGB = LCCSEGA + NCGTO3
             ELSE
                 XA = X4
                 YA = Y4
                 ZA = Z4
                 XB = X3
                 YB = Y3
                 ZB = Z3
                 SHELLA = SHELL4
                 SHELLB = SHELL3
                 NPGTOA = NPGTO4
                 NPGTOB = NPGTO3
                 NCGTOA = NCGTO4
                 NCGTOB = NCGTO3
                 NXYZA = NXYZ4
                 NXYZB = NXYZ3
                 NRYA = NRY4
                 NRYB = NRY3
                 INDEXA = 4
                 INDEXB = 3
                 LEXPB = 1 + NPGTO1 + NPGTO2
                 LEXPA = LEXPB + NPGTO3
                 LCCB = 1 + NCC1 + NCC2
                 LCCA = LCCB + NCC3
                 LCCSEGB = 1 + NCGTO1 + NCGTO2
                 LCCSEGA = LCCSEGB + NCGTO3
             END IF
         END IF
C
C
C             ...the new A,B,C,D shells are set. Calculate the
C                following info: 1) control variables to be used
C                during contraction, 2) total shell values for
C                electrons 1 and 2 in [AB|CD], 3) their corresponding
C                cartesian monomial sizes and 4) the overall norm
C                factor SPNORM due to presence of s- or p-type shells.
C                The latter is necessary, because for such shells
C                there will be no calls to the cartesian normalization
C                or spherical transformation routines. The contribution
C                to SPNORM is very simple: each s-type shell -> * 1.0,
C                each p-type shell -> * 2.0.
C
C
         SWAPRS = NPGTOA .GT. NPGTOB
         SWAPTU = NPGTOC .GT. NPGTOD

         SHELLP = SHELLA + SHELLB
         SHELLQ = SHELLC + SHELLD

         NXYZP  = (SHELLP+1)*(SHELLP+2)/2
         NXYZQ  = (SHELLQ+1)*(SHELLQ+2)/2

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
         IF (SHELLD.EQ.1) THEN
             SPNORM = SPNORM + SPNORM
         END IF
C
C
C             ...calculate the coordinate differences between centers
C                A and B and between centers C and D and calculate the
C                square of the magnitude of the distances. Also
C                determine the number of non-zero coordinate differences
C                for each pair A,B and C,D.
C
C
         IF (.NOT.ATOMAB) THEN
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

         IF (.NOT.ATOMCD) THEN
             CDX = XC - XD
             CDY = YC - YD
             CDZ = ZC - ZD
             RNCDSQ = CDX * CDX + CDY * CDY + CDZ * CDZ
             NCDCOOR = 3
             IF (DABS(CDX).EQ.ZERO) NCDCOOR = NCDCOOR - 1
             IF (DABS(CDY).EQ.ZERO) NCDCOOR = NCDCOOR - 1
             IF (DABS(CDZ).EQ.ZERO) NCDCOOR = NCDCOOR - 1
         ELSE
             CDX = ZERO
             CDY = ZERO
             CDZ = ZERO
             RNCDSQ = ZERO
             NCDCOOR = 0
         END IF
C
C
C             ...if HRR contractions are to be performed, calculate
C                NCOLHRR (maximum # of HRR rotation matrix columns
C                needed to generate the final HRR rotation matrices) and
C                NROTHRR (maximum # of HRR rotation matrix elements).
C
C                First find maximum values for the HRR on the AB-part.
C
C
         NCOLHRR = 0
         NROTHRR = 0

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

             DO 100 SHELLH = 1,SHELLB
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
  100        CONTINUE

             NCOLHRR = NCOL
             NROTHRR = NROT

         END IF
C
C
C             ...next find maximum values for the HRR on the CD-part
C                and set overall maximum values.
C
C
         IF (SHELLD.NE.0) THEN

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

             DO 200 SHELLH = 1,SHELLD
                NXYZGO = NXYZGO - NXYZI
                NXYZHO = NXYZHO + SHELLH + 1
                NGHO = NXYZGO * NXYZHO

                IF (NCDCOOR.EQ.3) THEN
                    M = 1 + SHELLH/3
                    NROW = NROW +  M*(M + ADD (MOD(SHELLH,3)))
                ELSE IF (NCDCOOR.EQ.2) THEN
                    NROW = NROW + SHELLH/2 + 1
                ELSE IF (NCDCOOR.EQ.1) THEN
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

             NCOLHRR = MAX0 (NCOL,NCOLHRR)
             NROTHRR = MAX0 (NROT,NROTHRR)

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
