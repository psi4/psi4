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
         SUBROUTINE  ERD__SET_DERV_ABCD
     +
     +                    ( NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                      NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                      SHELL1,SHELL2,SHELL3,SHELL4,
     +                      X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                      EXP1,EXP2,EXP3,EXP4,
     +                      CC1,CC2,CC3,CC4,
     +                      DER1X,DER1Y,DER1Z,
     +                      DER2X,DER2Y,DER2Z,
     +                      DER3X,DER3Y,DER3Z,
     +                      DER4X,DER4Y,DER4Z,
     +                      SPHERIC,
     +
     +                            NCGTOA,NCGTOB,NCGTOC,NCGTOD,
     +                            NPGTOA,NPGTOB,NPGTOC,NPGTOD,
     +                            SHELLA,SHELLB,SHELLC,SHELLD,
     +                            SHELLP,SHELLQ,SHELLT,
     +                            MXSHELL,
     +                            XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     +                            NDERX,NDERY,NDERZ,
     +                            DERAX,DERAY,DERAZ,
     +                            DERBX,DERBY,DERBZ,
     +                            DERCX,DERCY,DERCZ,
     +                            DERDX,DERDY,DERDZ,
     +                            DERPX,DERPY,DERPZ,
     +                            DERQX,DERQY,DERQZ,
     +                            DIFFA,DIFFB,DIFFC,DIFFD,
     +                            DIFFX,DIFFY,DIFFZ,
     +                            CENEQS,
     +                            NZSHELL,NZNXYZ,
     +                            PRIMTYP,ANGMTYP,
     +                            ATOMAB,ATOMCD,
     +                            EQUALAB,EQUALCD,
     +                            ABX,ABY,ABZ,CDX,CDY,CDZ,
     +                            NABCOOR,
     +                            RNABSQ,RNCDSQ,
     +                            SPNORM,
     +                            NXYZET,NXYZP,NXYZBRA,NXYZT,
     +                            NXYZA,NXYZB,NXYZC,NXYZD,
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
C  OPERATION   : ERD__SET_DERV_ABCD
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine handles the logistics on how to evaluate
C                the (12|34) derivative integral batch in the most
C                efficient way. It performs the label map:
C
C                             (12|34) --> (AB|CD)
C
C                and returns some additional data of crucial importance
C                for evaluation of the (AB|CD) derivative integrals.
C
C                The freedom in making the internal association 1,2,3,4
C                -> A,B,C,D follows from the 8-fold permutational
C                symmetry of the derivative integrals in (12|34):
C
C                     (12|34) = (21|34) = (12|43) = (21|43) =
C                     (34|12) = (43|12) = (34|21) = (43|21)
C
C                where the first line has only 12 -> 21 and 34 -> 43
C                switches and the second line is obtained from the
C                first by bra <-> ket transpositions. Note, that this
C                8-fold permutational symmetry only holds as long
C                as we permute the differential operators with it.
C
C                The type of switch to be applied is governed by
C                the following rules:
C
C                  1) The final A,B,C,D shell labels should obey
C                     the relations A>=B and C>=D, since this means
C                     the least amount of work (# of steps) for the
C                     HRR procedure, either at primitive or at
C                     contracted level. Notice here, that if the
C                     HRR is to be applied at the primitive level,
C                     then it will be done before! any derivations
C                     are performed on either A,B,C or D shells.
C                     Hence in that case the size of these shells
C                     must include the derivative operator sizes
C                     too when figuring out the relations A>=B and
C                     C>=D.
C
C                  2) Whenever possible, try to apply the HRR at
C                     contracted level on the bra side. The HRR
C                     at contracted level can only be applied, if
C                     no differential operators are associated with
C                     the corresponding centers. Hence, if the
C                     differential operators are associated with
C                     only the bra in (12|34), then we interchange
C                     bra <-> ket, so that we can apply the HRR at
C                     contracted level to the undifferentiated side
C                     AB in (AB|CD).
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
C                   DERyp         =  the order of differentiation on
C                                    centers y = 1,2,3,4 with respect
C                                    to the p = x,y,z coordinates
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
C                   Xx,Yx,Zx      =  the x,y,z-coordinates for centers
C                                    x = A,B,C and D
C                   NDERp         =  the total order of differentiation
C                                    with respect to the p = x,y,z
C                                    coordinates
C                   DERyp         =  the order of differentiation on
C                                    centers y = A,B,C,D,P=A+B,Q=C+D
C                                    with respect to the p = x,y,z
C                                    coordinates
C                   DIFFy         =  is true, if differentiation will be
C                                    performed on centers y = A,B,C,D
C                                    involving the y = x,y,z coordinates
C                   CENEQS (I,J)  =  center equality indicator of size
C                                    4 x 4. The matrix is defined as
C                                    follows: if the centers indexed by
C                                    I and J are equal => value = 1,
C                                    if not => value = 0. The indices
C                                    correspond to the A,B,C,D ordering,
C                                    i.e. 1st index -> A, 2nd -> B, etc
C                   NZSHELL (I)   =  the I-th nonzero shell type within
C                                    the shell type center sequence
C                                    A,B,C,D
C                   NZNXYZ (I)    =  the I-th # of cartesian monomials
C                                    corresponding to the I-th nonzero
C                                    shell in NZSHELL (I)
C                   PRIMTYP       =  character variable, indicating
C                                    which type of primitives will
C                                    be generated and contracted. Can
C                                    be only 'E0CD' or 'ABCD'
C                   ANGMTYP       =  character variable, indicating
C                                    the overall angular momentum type
C                                    combination without observing
C                                    the order. Can be only 'SSSS',
C                                    'SSSX','SSXX','SXXX' or 'XXXX',
C                                    where the S-symbol indicates the
C                                    presence of an s-shell and the
C                                    X-symbol the presence of a shell
C                                    >= p-shell
C                   ATOMxy        =  indicates, if atoms x and y are
C                                    considered to be equal for the
C                                    pairs xy = AB and CD
C                   EQUALxy       =  indicates, if csh x and csh y are
C                                    considered to be equal for the
C                                    pairs xy = AB and CD
C                   xxX,xxY,xxZ   =  the x,y,z-coordinate differences
C                                    for centers xx = AB and CD
C                   NABCOOR       =  # of non-zero x,y,z-coordinate
C                                    differences for centers A and B
C                   RNxxSQ        =  square of the magnitude of the
C                                    distance between centers xx = AB
C                                    and CD
C                   SPNORM        =  normalization factor due to
C                                    presence of s- and p-type shells.
C                                    For each s-type shell there is a
C                                    factor of 1 and for each p-type
C                                    shell a factor of 2
C                   NXYZET        =  sum of # of cartesian monomials
C                                    for all shells in the range
C                                    E = A,...,A+B
C                   NXYZP         =  # of cartesian monomials for the
C                                    P = A + B shell
C                   NXYZBRA       =  total # of cartesian monomials
C                                    for the bra side of the primitive
C                                    derivative integral batch 
C                   NXYZT         =  total # of cartesian monomials
C                                    for the primitive derivative
C                                    integral batch 
C                   NXYZx         =  # of cartesian monomials for csh x
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
C                   NXYZHRR       =  maximum dimension of monomial part
C                                    after contraction of primitives.
C                                    Takes into consideration eventual
C                                    HRR at contracted level on E -> AB
C                                    part and cartesian -> spherical
C                                    transformations
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
C
C             ...include files and declare variables.
C
C
         IMPLICIT  NONE

         LOGICAL   ALERT
         LOGICAL   ATOM12,ATOM13,ATOM14,ATOM23,ATOM24,ATOM34
         LOGICAL   ATOMIC,ATOMAB,ATOMCD
         LOGICAL   DIFF1,DIFF2,DIFF3,DIFF4
         LOGICAL   DIFFA,DIFFB,DIFFC,DIFFD
         LOGICAL   DIFFBRA,DIFFKET
         LOGICAL   DIFFX,DIFFY,DIFFZ
         LOGICAL   EMPTY
         LOGICAL   EQUAL12,EQUAL34
         LOGICAL   EQUALAB,EQUALCD
         LOGICAL   SAMEAB,SAMEAC,SAMEAD,SAMEBC,SAMEBD,SAMECD
         LOGICAL   SPHERIC
         LOGICAL   SWAP12,SWAP34,SWAPRS,SWAPTU,TR1234
         LOGICAL   SWAPX,SWAPY,SWAPZ

         CHARACTER*4   ANGMTYP
         CHARACTER*4   PRIMTYP

         INTEGER     ATOM1,ATOM2,ATOM3,ATOM4
         INTEGER     ATOMA,ATOMB,ATOMC,ATOMD
         INTEGER     DER1X,DER2X,DER3X,DER4X
         INTEGER     DER1Y,DER2Y,DER3Y,DER4Y
         INTEGER     DER1Z,DER2Z,DER3Z,DER4Z
         INTEGER     DERAX,DERBX,DERCX,DERDX
         INTEGER     DERAY,DERBY,DERCY,DERDY
         INTEGER     DERAZ,DERBZ,DERCZ,DERDZ
         INTEGER     DERPX,DERPY,DERPZ
         INTEGER     DERQX,DERQY,DERQZ
         INTEGER     I,J,M
         INTEGER     INDEXA,INDEXB,INDEXC,INDEXD
         INTEGER     LCCSEGA,LCCSEGB,LCCSEGC,LCCSEGD
         INTEGER     LCCA,LCCB,LCCC,LCCD
         INTEGER     LEXPA,LEXPB,LEXPC,LEXPD
         INTEGER     MXSHELL
         INTEGER     NABCOOR
         INTEGER     NCC1,NCC2,NCC3
         INTEGER     NCGTO1,NCGTO2,NCGTO3,NCGTO4
         INTEGER     NCGTOA,NCGTOB,NCGTOC,NCGTOD
         INTEGER     NDERX,NDERY,NDERZ
         INTEGER     NGH,NGHO
         INTEGER     NPGTO1,NPGTO2,NPGTO3,NPGTO4
         INTEGER     NPGTOA,NPGTOB,NPGTOC,NPGTOD
         INTEGER     NROW,NCOL,NROT
         INTEGER     NRY1,NRY2,NRY3,NRY4
         INTEGER     NRYA,NRYB,NRYC,NRYD
         INTEGER     NXYZ1,NXYZ2,NXYZ3,NXYZ4
         INTEGER     NXYZA,NXYZB,NXYZC,NXYZD
         INTEGER     NXYZET,NXYZP,NXYZBRA,NXYZT
         INTEGER     NXYZG,NXYZH,NXYZI,NXYZGO,NXYZHO
         INTEGER     NXYZHRR,NCOLHRR,NROTHRR
         INTEGER     SHELL1,SHELL2,SHELL3,SHELL4
         INTEGER     SHELLA,SHELLB,SHELLC,SHELLD
         INTEGER     SHELLG,SHELLH
         INTEGER     SHELLP,SHELLQ,SHELLT

         INTEGER     ADD  (0:2)

         INTEGER     NZNXYZ  (1:4)
         INTEGER     NZSHELL (1:4)

         INTEGER     CENEQS (1:4,1:4)

         DOUBLE PRECISION  ABX,ABY,ABZ,CDX,CDY,CDZ
         DOUBLE PRECISION  IVSUM1,IVSUM2,IVSUM3,IVSUM4
         DOUBLE PRECISION  RNABSQ,RNCDSQ
         DOUBLE PRECISION  SPNORM
         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
         DOUBLE PRECISION  ZERO,HALF,ONE

         DOUBLE PRECISION  EXP1 (1:NPGTO1)
         DOUBLE PRECISION  EXP2 (1:NPGTO2)
         DOUBLE PRECISION  EXP3 (1:NPGTO3)
         DOUBLE PRECISION  EXP4 (1:NPGTO4)

         DOUBLE PRECISION  CC1 (1:NPGTO1,1:NCGTO1)
         DOUBLE PRECISION  CC2 (1:NPGTO2,1:NCGTO2)
         DOUBLE PRECISION  CC3 (1:NPGTO3,1:NCGTO3)
         DOUBLE PRECISION  CC4 (1:NPGTO4,1:NCGTO4)

         DATA  ADD    /0,0,1/

         PARAMETER  (ZERO  =  0.D0)
         PARAMETER  (HALF  =  0.5D0)
         PARAMETER  (ONE   =  1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...generate all 1,2,3,4 data. Decide as early as
C                possible, if a zero batch of derivative integrals
C                is expected.
C
C
         EMPTY = .FALSE.

         ATOM12 = (X1.EQ.X2) .AND. (Y1.EQ.Y2) .AND. (Z1.EQ.Z2)
         ATOM23 = (X2.EQ.X3) .AND. (Y2.EQ.Y3) .AND. (Z2.EQ.Z3)
         ATOM34 = (X3.EQ.X4) .AND. (Y3.EQ.Y4) .AND. (Z3.EQ.Z4)

         ATOMIC = (ATOM12 .AND. ATOM34 .AND. ATOM23)

         IF (ATOMIC) THEN
             EMPTY = .TRUE.
             RETURN
         END IF

         MXSHELL = MAX0 (SHELL1,SHELL2,SHELL3,SHELL4)

         ATOM13 = (X1.EQ.X3) .AND. (Y1.EQ.Y3) .AND. (Z1.EQ.Z3)
         ATOM14 = (X1.EQ.X4) .AND. (Y1.EQ.Y4) .AND. (Z1.EQ.Z4)
         ATOM24 = (X2.EQ.X4) .AND. (Y2.EQ.Y4) .AND. (Z2.EQ.Z4)

         ATOM1 = 1
         ATOM2 = 2
         ATOM3 = 3
         ATOM4 = 4

         IF (ATOM12) ATOM1 = ATOM2
         IF (ATOM13) ATOM1 = ATOM3
         IF (ATOM14) ATOM1 = ATOM4
         IF (ATOM23) ATOM2 = ATOM3
         IF (ATOM24) ATOM2 = ATOM4
         IF (ATOM34) ATOM3 = ATOM4
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
               DO 10 I = 1,NPGTO1
                  EQUAL12 = EQUAL12 .AND. (EXP1(I).EQ.EXP2(I))
   10          CONTINUE
               IF (EQUAL12) THEN
                 DO 12 J = 1,NCGTO1
                    IF (EQUAL12) THEN
                      DO 14 I = 1,NPGTO1
                         EQUAL12 = EQUAL12 .AND. (CC1(I,J).EQ.CC2(I,J))
   14                 CONTINUE
                    END IF
   12            CONTINUE
               END IF
             END IF
         END IF

         EQUAL34 = ATOM34

         IF (EQUAL34) THEN
             EQUAL34 =     (SHELL3 .EQ. SHELL4)
     +               .AND. (NPGTO3 .EQ. NPGTO4)
     +               .AND. (NCGTO3 .EQ. NCGTO4)
             IF (EQUAL34) THEN
               DO 30 I = 1,NPGTO3
                  EQUAL34 = EQUAL34 .AND. (EXP3(I).EQ.EXP4(I))
   30          CONTINUE
               IF (EQUAL34) THEN
                 DO 32 J = 1,NCGTO3
                    IF (EQUAL34) THEN
                      DO 34 I = 1,NPGTO3
                         EQUAL34 = EQUAL34 .AND. (CC3(I,J).EQ.CC4(I,J))
   34                 CONTINUE
                    END IF
   32            CONTINUE
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
C                Be careful to include the derivative operators!
C                This is a bit trickey, since we have three possible
C                situations due to the x,y,z components. The rule
C                right now is to make the swapping, if at least
C                two of the components are in favor to do so.
C
C
         SWAPX = (SHELL1 + DER1X) .LT. (SHELL2 + DER2X)
         SWAPY = (SHELL1 + DER1Y) .LT. (SHELL2 + DER2Y)
         SWAPZ = (SHELL1 + DER1Z) .LT. (SHELL2 + DER2Z)

         SWAP12 =      (SWAPX .AND. SWAPY)
     +            .OR. (SWAPX .AND. SWAPZ)
     +            .OR. (SWAPY .AND. SWAPZ)

         SWAPX = (SHELL3 + DER3X) .LT. (SHELL4 + DER4X)
         SWAPY = (SHELL3 + DER3Y) .LT. (SHELL4 + DER4Y)
         SWAPZ = (SHELL3 + DER3Z) .LT. (SHELL4 + DER4Z)

         SWAP34 =      (SWAPX .AND. SWAPY)
     +            .OR. (SWAPX .AND. SWAPZ)
     +            .OR. (SWAPY .AND. SWAPZ)
C
C
C             ...analyze the situation (center location + coordinate
C                type) of the derivative operators.
C
C
         DIFFX = (DER1X + DER2X + DER3X  + DER4X) .NE. 0
         DIFFY = (DER1Y + DER2Y + DER3Y  + DER4Y) .NE. 0
         DIFFZ = (DER1Z + DER2Z + DER3Z  + DER4Z) .NE. 0

         DIFF1 = (DER1X + DER1Y + DER1Z) .NE. 0
         DIFF2 = (DER2X + DER2Y + DER2Z) .NE. 0
         DIFF3 = (DER3X + DER3Y + DER3Z) .NE. 0
         DIFF4 = (DER4X + DER4Y + DER4Z) .NE. 0

         DIFFBRA = DIFF1 .OR. DIFF2
         DIFFKET = DIFF3 .OR. DIFF4

         IF (DIFFBRA .AND. DIFFKET) THEN
             PRIMTYP = 'ABCD'
             TR1234 = .FALSE.
         ELSE IF (DIFFBRA .AND. .NOT.DIFFKET) THEN
             PRIMTYP = 'E0CD'
             TR1234 = .TRUE.
         ELSE IF (DIFFKET .AND. .NOT.DIFFBRA) THEN
             PRIMTYP = 'E0CD'
             TR1234 = .FALSE.
         ELSE
             WRITE (*,*) ' No differential operator! '
             WRITE (*,*) ' DIFFBRA,DIFFKET = ',DIFFBRA,DIFFKET
             WRITE (*,*) ' erd__set_derv_abcd '
             STOP
         END IF
C
C
C             ...according to the previously gathered info, set the
C                new A,B,C,D shells, # of primitives + contraction
C                coeffs as well as pointers to the alpha exponents
C                and contraction coefficients.
C
C
         NCC1 = NPGTO1 * NCGTO1
         NCC2 = NPGTO2 * NCGTO2
         NCC3 = NPGTO3 * NCGTO3

         IF (.NOT.TR1234) THEN

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
                 ATOMA = ATOM1
                 ATOMB = ATOM2
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
                 DIFFA = DIFF1
                 DIFFB = DIFF2
                 DERAX = DER1X
                 DERAY = DER1Y
                 DERAZ = DER1Z
                 DERBX = DER2X
                 DERBY = DER2Y
                 DERBZ = DER2Z
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
                 ATOMA = ATOM2
                 ATOMB = ATOM1
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
                 DIFFA = DIFF2
                 DIFFB = DIFF1
                 DERAX = DER2X
                 DERAY = DER2Y
                 DERAZ = DER2Z
                 DERBX = DER1X
                 DERBY = DER1Y
                 DERBZ = DER1Z
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
                 ATOMC = ATOM3
                 ATOMD = ATOM4
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
                 DIFFC = DIFF3
                 DIFFD = DIFF4
                 DERCX = DER3X
                 DERCY = DER3Y
                 DERCZ = DER3Z
                 DERDX = DER4X
                 DERDY = DER4Y
                 DERDZ = DER4Z
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
                 ATOMC = ATOM4
                 ATOMD = ATOM3
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
                 DIFFC = DIFF4
                 DIFFD = DIFF3
                 DERCX = DER4X
                 DERCY = DER4Y
                 DERCZ = DER4Z
                 DERDX = DER3X
                 DERDY = DER3Y
                 DERDZ = DER3Z
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
                 ATOMC = ATOM1
                 ATOMD = ATOM2
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
                 DIFFC = DIFF1
                 DIFFD = DIFF2
                 DERCX = DER1X
                 DERCY = DER1Y
                 DERCZ = DER1Z
                 DERDX = DER2X
                 DERDY = DER2Y
                 DERDZ = DER2Z
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
                 ATOMC = ATOM2
                 ATOMD = ATOM1
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
                 DIFFC = DIFF2
                 DIFFD = DIFF1
                 DERCX = DER2X
                 DERCY = DER2Y
                 DERCZ = DER2Z
                 DERDX = DER1X
                 DERDY = DER1Y
                 DERDZ = DER1Z
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
                 ATOMA = ATOM3
                 ATOMB = ATOM4
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
                 DIFFA = DIFF3
                 DIFFB = DIFF4
                 DERAX = DER3X
                 DERAY = DER3Y
                 DERAZ = DER3Z
                 DERBX = DER4X
                 DERBY = DER4Y
                 DERBZ = DER4Z
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
                 ATOMA = ATOM4
                 ATOMB = ATOM3
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
                 DIFFA = DIFF4
                 DIFFB = DIFF3
                 DERAX = DER4X
                 DERAY = DER4Y
                 DERAZ = DER4Z
                 DERBX = DER3X
                 DERBY = DER3Y
                 DERBZ = DER3Z
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
C                P=A+B,Q=C+D and T=A+B+C+D, 3) the total derivative
C                orders per cartesian component for the shell sums
C                P=A+B and Q=C+D and 4) the overall norm factor SPNORM
C                due to presence of s- or p-type shells. The latter
C                is necessary, because for such shells there will be
C                no calls to the cartesian normalization or spherical
C                transformation routines. The contribution to SPNORM
C                is very simple: each s-type shell -> * 1.0, each
C                p-type shell -> * 2.0.
C
C
         SWAPRS = NPGTOA .GT. NPGTOB
         SWAPTU = NPGTOC .GT. NPGTOD

         SHELLP = SHELLA + SHELLB
         SHELLQ = SHELLC + SHELLD
         SHELLT = SHELLP + SHELLQ

         DERPX = DERAX + DERBX
         DERPY = DERAY + DERBY
         DERPZ = DERAZ + DERBZ
         DERQX = DERCX + DERDX
         DERQY = DERCY + DERDY
         DERQZ = DERCZ + DERDZ

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
C                square of the magnitude of the distances.
C
C
         IF (.NOT.ATOMAB) THEN
             ABX = XA - XB
             ABY = YA - YB
             ABZ = ZA - ZB
             RNABSQ = ABX * ABX + ABY * ABY + ABZ * ABZ
         ELSE
             ABX = ZERO
             ABY = ZERO
             ABZ = ZERO
             RNABSQ = ZERO
         END IF

         IF (.NOT.ATOMCD) THEN
             CDX = XC - XD
             CDY = YC - YD
             CDZ = ZC - ZD
             RNCDSQ = CDX * CDX + CDY * CDY + CDZ * CDZ
         ELSE
             CDX = ZERO
             CDY = ZERO
             CDZ = ZERO
             RNCDSQ = ZERO
         END IF
C
C
C             ...set the A,B,C,D center equality map. This map
C                is a 4 x 4 dimensional array, which contains
C                a 1 in (i,j) position if centers indexed by
C                i and j are equal and a 0 otherwise. The indices
C                of the map are such that they have the following
C                correspondence:
C
C                             Index i value: 1 2 3 4
C                     corresponds to center: A B C D
C
C
         CENEQS (1,1) = 1
         CENEQS (2,1) = 1
         CENEQS (3,1) = 1
         CENEQS (4,1) = 1
         CENEQS (2,2) = 1
         CENEQS (3,2) = 1
         CENEQS (4,2) = 1
         CENEQS (3,3) = 1
         CENEQS (4,3) = 1
         CENEQS (4,4) = 1

         IF (ATOMA .NE. ATOMB) CENEQS (2,1) = 0
         IF (ATOMA .NE. ATOMC) CENEQS (3,1) = 0
         IF (ATOMA .NE. ATOMD) CENEQS (4,1) = 0
         IF (ATOMB .NE. ATOMC) CENEQS (3,2) = 0
         IF (ATOMB .NE. ATOMD) CENEQS (4,2) = 0
         IF (ATOMC .NE. ATOMD) CENEQS (4,3) = 0

         CENEQS (1,2) = CENEQS (2,1)
         CENEQS (1,3) = CENEQS (3,1)
         CENEQS (1,4) = CENEQS (4,1)
         CENEQS (2,3) = CENEQS (3,2)
         CENEQS (2,4) = CENEQS (4,2)
         CENEQS (3,4) = CENEQS (4,3)
C
C
C             ...check, if the derivative orders for each center
C                and coordinate match the center equality pattern.
C                If not, stop the derivative integral calculation.
C                If yes, determine the order of the derivation for
C                each coordinate.
C
C
         SAMEAB = CENEQS (2,1) .EQ. 1
         SAMEAC = CENEQS (3,1) .EQ. 1
         SAMEAD = CENEQS (4,1) .EQ. 1
         SAMEBC = CENEQS (3,2) .EQ. 1
         SAMEBD = CENEQS (4,2) .EQ. 1
         SAMECD = CENEQS (4,3) .EQ. 1

         IVSUM1 = ONE / DFLOAT (  CENEQS(1,1)
     +                          + CENEQS(2,1)
     +                          + CENEQS(3,1)
     +                          + CENEQS(4,1) )
         IVSUM2 = ONE / DFLOAT (  CENEQS(1,2)
     +                          + CENEQS(2,2)
     +                          + CENEQS(3,2)
     +                          + CENEQS(4,2) )
         IVSUM3 = ONE / DFLOAT (  CENEQS(1,3)
     +                          + CENEQS(2,3)
     +                          + CENEQS(3,3)
     +                          + CENEQS(4,3) )
         IVSUM4 = ONE / DFLOAT (  CENEQS(1,4)
     +                          + CENEQS(2,4)
     +                          + CENEQS(3,4)
     +                          + CENEQS(4,4) )
C
C
C             ...x coordinate derivatives (if any).
C
C
         IF (DIFFX) THEN

             ALERT =      (SAMEAB .AND. DERAX.NE.DERBX)
     +               .OR. (SAMEAC .AND. DERAX.NE.DERCX)
     +               .OR. (SAMEAD .AND. DERAX.NE.DERDX)
     +               .OR. (SAMEBC .AND. DERBX.NE.DERCX)
     +               .OR. (SAMEBD .AND. DERBX.NE.DERDX)
     +               .OR. (SAMECD .AND. DERCX.NE.DERDX)

             IF (ALERT) THEN
                 WRITE (*,*) ' Center equality / x derv mismatch! '
                 WRITE (*,*) ' DERAX,DERBX,DERCX,DERDX = ',
     +                         DERAX,DERBX,DERCX,DERDX
                 WRITE (*,*) ' erd__set_derv_abcd '
                 STOP
             END IF

             NDERX = INT (  DFLOAT (DERAX) * IVSUM1
     +                    + DFLOAT (DERBX) * IVSUM2
     +                    + DFLOAT (DERCX) * IVSUM3
     +                    + DFLOAT (DERDX) * IVSUM4
     +                    + HALF )
         ELSE
             NDERX = 0
         END IF
C
C
C             ...y coordinate derivatives (if any).
C
C
         IF (DIFFY) THEN

             ALERT =      (SAMEAB .AND. DERAY.NE.DERBY)
     +               .OR. (SAMEAC .AND. DERAY.NE.DERCY)
     +               .OR. (SAMEAD .AND. DERAY.NE.DERDY)
     +               .OR. (SAMEBC .AND. DERBY.NE.DERCY)
     +               .OR. (SAMEBD .AND. DERBY.NE.DERDY)
     +               .OR. (SAMECD .AND. DERCY.NE.DERDY)

             IF (ALERT) THEN
                 WRITE (*,*) ' Center equality / y derv mismatch! '
                 WRITE (*,*) ' DERAY,DERBY,DERCY,DERDY = ',
     +                         DERAY,DERBY,DERCY,DERDY
                 WRITE (*,*) ' erd__set_derv_abcd '
                 STOP
             END IF

             NDERY = INT (  DFLOAT (DERAY) * IVSUM1
     +                    + DFLOAT (DERBY) * IVSUM2
     +                    + DFLOAT (DERCY) * IVSUM3
     +                    + DFLOAT (DERDY) * IVSUM4
     +                    + HALF )
         ELSE
             NDERY = 0
         END IF
C
C
C             ...z coordinate derivatives (if any).
C
C
         IF (DIFFZ) THEN

             ALERT =      (SAMEAB .AND. DERAZ.NE.DERBZ)
     +               .OR. (SAMEAC .AND. DERAZ.NE.DERCZ)
     +               .OR. (SAMEAD .AND. DERAZ.NE.DERDZ)
     +               .OR. (SAMEBC .AND. DERBZ.NE.DERCZ)
     +               .OR. (SAMEBD .AND. DERBZ.NE.DERDZ)
     +               .OR. (SAMECD .AND. DERCZ.NE.DERDZ)

             IF (ALERT) THEN
                 WRITE (*,*) ' Center equality / z derv mismatch! '
                 WRITE (*,*) ' DERAZ,DERBZ,DERCZ,DERDZ = ',
     +                         DERAZ,DERBZ,DERCZ,DERDZ
                 WRITE (*,*) ' erd__set_derv_abcd '
                 STOP
             END IF

             NDERZ = INT (  DFLOAT (DERAZ) * IVSUM1
     +                    + DFLOAT (DERBZ) * IVSUM2
     +                    + DFLOAT (DERCZ) * IVSUM3
     +                    + DFLOAT (DERDZ) * IVSUM4
     +                    + HALF )
         ELSE
             NDERZ = 0
         END IF
C
C
C             ...calculate data relevant for the two possible cases:
C
C                 1) A HRR at contracted level will be applied
C                    to the E -> AB part.
C
C                 2) No HRR at contracted level.
C
C                Explanation of variable NXYZHRR:
C
C                  The value of NXYZHRR, which will contain the maximum
C                  dimension that will be encountered for the monomial
C                  transformation part after contraction. The monomial
C                  transformations consist of a possible HRR at
C                  contracted level and (if any) cartesian -> spherical
C                  transformation sequences. For the two cases above
C                  we have:
C
C                 Case 1):
C
C                  i) initial dimension:  NXYZET * NXYZC * NXYZD
C                 ii) cart -> sph on CD:  NXYZET * NRYC * NRYD
C                iii) perform HRR on AB:  NXYZA * NXYZB * NRYC * NRYD
C                 iv) cart -> sph on AB:  NRYA * NRYB * NRYC * NRYD
C
C                 Case 2):
C
C                  i) initial dimension:  NXYZA * NXYZB * NXYZC * NXYZD
C                 ii) cart -> sph on CD:  NXYZA * NXYZB * NRYC * NRYD
C                iii) cart -> sph on AB:  NRYA * NRYB * NRYC * NRYD
C
C
C                The only dimension peaks are steps i) and iii) in
C                Case 1) and step i) in Case 2).
C
C                Other info that needs to be calculated for Case 1):
C
C                  a) total monomial dimension for E = A,...,A+B
C                  b) # of nonzero coordinate differences between
C                     centers A and B
C                  c) values NCOLHRR (maximum # of HRR rotation
C                     matrix columns needed to generate the final
C                     HRR rotation matrices) and NROTHRR (maximum
C                     # of HRR rotation matrix elements)
C
C                We also need to evaluate the total # of monomial
C                quadruplets expected for the primitive derivative
C                integral batch for both cases.
C
C
         IF (PRIMTYP.EQ.'E0CD') THEN

             NXYZP   =    ((SHELLP+1)*(SHELLP+2))/2
             NXYZET  =   (((SHELLP+1)*(SHELLP+2)*(SHELLP+3))/6)
     +                 - (((SHELLA  )*(SHELLA+1)*(SHELLA+2))/6)

             NXYZBRA = NXYZET
             NXYZT   = NXYZBRA * NXYZC * NXYZD
             NXYZHRR = MAX0 (NXYZT,NXYZA*NXYZB*NRYC*NRYD)

             NABCOOR = 3
             IF (DABS(ABX).EQ.ZERO) NABCOOR = NABCOOR - 1
             IF (DABS(ABY).EQ.ZERO) NABCOOR = NABCOOR - 1
             IF (DABS(ABZ).EQ.ZERO) NABCOOR = NABCOOR - 1

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
  100            CONTINUE

                 NCOLHRR = NCOL
                 NROTHRR = NROT

             END IF

             M = 0

             IF (SHELLC.GT.0) THEN
                 M = M + 1
                 NZNXYZ (M) = NXYZC
                 NZSHELL (M) = SHELLC
             END IF

             IF (SHELLD.GT.0) THEN
                 M = M + 1
                 NZNXYZ (M) = NXYZD
                 NZSHELL (M) = SHELLD
             END IF

             IF (M.EQ.0) THEN
                 ANGMTYP = 'SSXX'
             ELSE IF (M.EQ.1) THEN
                 ANGMTYP = 'SXXX'
             ELSE
                 ANGMTYP = 'XXXX'
             END IF

         ELSE

             NXYZBRA = NXYZA * NXYZB
             NXYZT   = NXYZBRA * NXYZC * NXYZD
             NXYZHRR = NXYZT

             M = 0

             IF (SHELLA.GT.0) THEN
                 M = M + 1
                 NZNXYZ (M) = NXYZA
                 NZSHELL (M) = SHELLA
             END IF

             IF (SHELLB.GT.0) THEN
                 M = M + 1
                 NZNXYZ (M) = NXYZB
                 NZSHELL (M) = SHELLB
             END IF

             IF (SHELLC.GT.0) THEN
                 M = M + 1
                 NZNXYZ (M) = NXYZC
                 NZSHELL (M) = SHELLC
             END IF

             IF (SHELLD.GT.0) THEN
                 M = M + 1
                 NZNXYZ (M) = NXYZD
                 NZSHELL (M) = SHELLD
             END IF

             IF (M.EQ.0) THEN
                 ANGMTYP = 'SSSS'
             ELSE IF (M.EQ.1) THEN
                 ANGMTYP = 'SSSX'
             ELSE IF (M.EQ.2) THEN
                 ANGMTYP = 'SSXX'
             ELSE IF (M.EQ.3) THEN
                 ANGMTYP = 'SXXX'
             ELSE
                 ANGMTYP = 'XXXX'
             END IF

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
