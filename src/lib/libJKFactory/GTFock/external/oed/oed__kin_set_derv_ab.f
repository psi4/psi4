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
         SUBROUTINE  OED__KIN_SET_DERV_AB
     +
     +                    ( NCGTO1,NCGTO2,
     +                      NPGTO1,NPGTO2,
     +                      SHELL1,SHELL2,
     +                      X1,Y1,Z1,X2,Y2,Z2,
     +                      EXP1,EXP2,
     +                      CC1,CC2,
     +                      DER1X,DER1Y,DER1Z,
     +                      DER2X,DER2Y,DER2Z,
     +                      SPHERIC,
     +
     +                                 NCGTOA,NCGTOB,
     +                                 NPGTOA,NPGTOB,
     +                                 SHELLA,SHELLB,SHELLP,
     +                                 MXSHELL,
     +                                 XA,YA,ZA,XB,YB,ZB,
     +                                 NDERX,NDERY,NDERZ,
     +                                 DERAX,DERAY,DERAZ,
     +                                 DERBX,DERBY,DERBZ,
     +                                 DIFFA,DIFFB,
     +                                 DIFFX,DIFFY,DIFFZ,
     +                                 ABX,ABY,ABZ,RNABSQ,
     +                                 SPNORM,
     +                                 NXYZA,NXYZB,NXYZT,
     +                                 NRYA,NRYB,
     +                                 INDEXA,INDEXB,
     +                                 SWAP12,SWAPRS,
     +                                 LEXPA,LEXPB,
     +                                 LCCA,LCCB,
     +                                 LCCSEGA,LCCSEGB,
     +                                 EMPTY )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__KIN_SET_DERV_AB
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine handles the logistics on how to evaluate
C                the (1|2) derivative kinetic integral batch in the
C                most efficient way. It performs the label map:
C
C                               (1|2) --> (A|B)
C
C                and returns some additional data of crucial importance
C                for evaluation of the (A|B) derivative kinetic
C                integrals.
C
C                The freedom we have in making the internal association
C                1,2 -> A,B follows from the 2-fold permutational
C                symmetry of the kinetic integrals in (1|2):
C
C                                 (1|2) = (2|1)
C
C                Note, that this permutational symmetry only holds as
C                long as we permute the differential operators with it.
C
C                Application of the 1 <-> 2 switch is pretty irrelevant
C                here, due to the fact the routine which generates the
C                initial 1D integrals has the possibility to either
C                apply the HRR steps on the A or B shell. Still we
C                choose here to have A>=B for convenience and possibly
C                slightly better performance.
C
C
C                  Input (x = 1 and 2):
C
C                    NCGTOx       =  # of contractions for csh x
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x
C                    SHELLx       =  the shell type for csh x
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = 1 and 2
C                    EXPx         =  primitive exponents for csh x
C                    CCx          =  contraction coeffs for csh x
C                    DERyp        =  the order of differentiation on
C                                    centers y = 1 and 2 with respect
C                                    to the p = x,y,z coordinates
C                    SPHERIC      =  is true, if spherical integrals
C                                    are wanted, false if cartesian
C                                    ones are wanted
C
C                  Output (x = A and B):
C
C                    NCGTOx       =  # of contractions for csh x
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x
C                    SHELLx       =  the shell type for csh x
C                    SHELLP       =  the shell sum A+B
C                    MXSHELL      =  the largest (maximum) shell type
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = A and B
C                    NDERp        =  the total order of differentiation
C                                    with respect to the p = x,y,z
C                                    coordinates
C                    DERyp        =  the order of differentiation on
C                                    centers y = A and B with respect
C                                    to the p = x,y,z coordinates
C                    DIFFy        =  is true, if differentiation will
C                                    be performed on centers y = A and B
C                                    involving the y = x,y,z coordinates
C                    ABX,ABY,ABZ  =  the x,y,z-coordinate differences
C                                    between centers A and B
C                    RNABSQ       =  square of the magnitude of the
C                                    distance between centers A and B
C                    SPNORM       =  normalization factor due to
C                                    presence of s- and p-type shells.
C                                    For each s-type shell there is a
C                                    factor of 1 and for each p-type
C                                    shell a factor of 2
C                    NXYZx        =  # of cartesian monomials for csh x
C                    NXYZT        =  total # cartesian monomial
C                                    combinations, that is NXYZA * NXYZB
C                    NRYx         =  # of spherical functions for csh x
C                    INDEXx       =  index A,B -> 1,2 map
C                    SWAP12       =  is .true., if a swap 1 <-> 2 has
C                                    been performed
C                    SWAPRS       =  is set .true. if the contraction
C                                    order of the primitives pair AB
C                                    will be performed in reverse order
C                                    BA for efficiency reasons
C                    LEXPx        =  pointers to locate appropriate
C                                    section of the exponent array
C                                    corresponding to csh x
C                    LCCx         =  pointers to locate appropriate
C                                    section of the contraction coeff
C                                    array corresponding to csh x
C                    LCCSEGx      =  pointers to locate appropriate
C                                    section of the lowest and highest
C                                    primitive index array defining
C                                    segmented contraction boundaries
C                                    for csh x
C                    EMPTY        =  logical flag, indicating if an
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

         LOGICAL     ATOMIC
         LOGICAL     DIFFA,DIFFB
         LOGICAL     DIFFX,DIFFY,DIFFZ
         LOGICAL     EMPTY
         LOGICAL     SPHERIC
         LOGICAL     SWAP12,SWAPRS

         INTEGER     DER1X,DER1Y,DER1Z
         INTEGER     DER2X,DER2Y,DER2Z
         INTEGER     DERAX,DERAY,DERAZ
         INTEGER     DERBX,DERBY,DERBZ
         INTEGER     INDEXA,INDEXB
         INTEGER     LCCSEGA,LCCSEGB
         INTEGER     LCCA,LCCB
         INTEGER     LEXPA,LEXPB
         INTEGER     MXSHELL
         INTEGER     NCGTO1,NCGTO2
         INTEGER     NCGTOA,NCGTOB
         INTEGER     NDERX,NDERY,NDERZ
         INTEGER     NPGTO1,NPGTO2
         INTEGER     NPGTOA,NPGTOB
         INTEGER     NRYA,NRYB
         INTEGER     NXYZA,NXYZB,NXYZT
         INTEGER     SHELL1,SHELL2
         INTEGER     SHELLA,SHELLB
         INTEGER     SHELLP

         DOUBLE PRECISION  ABX,ABY,ABZ
         DOUBLE PRECISION  RNABSQ
         DOUBLE PRECISION  SPNORM
         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  EXP1 (1:NPGTO1)
         DOUBLE PRECISION  EXP2 (1:NPGTO2)

         DOUBLE PRECISION  CC1 (1:NPGTO1,1:NCGTO1)
         DOUBLE PRECISION  CC2 (1:NPGTO2,1:NCGTO2)

         PARAMETER  (ZERO    =  0.D0)
         PARAMETER  (ONE     =  1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...generate all 1,2 data. Decide as early as possible,
C                if a zero batch of kinetic derivative integrals is
C                expected.
C
C
         ATOMIC = (X1.EQ.X2) .AND. (Y1.EQ.Y2) .AND. (Z1.EQ.Z2)

         IF (ATOMIC) THEN
             EMPTY = .TRUE.
             RETURN
         END IF

         EMPTY  = .FALSE.
         SHELLP = SHELL1 + SHELL2
         MXSHELL = MAX0 (SHELL1,SHELL2)
C
C
C             ...decide on the 1 <-> 2 swapping.
C
C
         SWAP12 = SHELL1 .LT. SHELL2
C
C
C             ...according to the previously gathered info, set the
C                new A,B shells, # of primitives + contraction coeffs
C                as well as pointers to the alpha exponents and
C                contraction coefficients.
C
C
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
             LCCB = LCCA + NPGTO1 * NCGTO1
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
             LCCA = LCCB + NPGTO1 * NCGTO1
             LCCSEGB = 1
             LCCSEGA = LCCSEGB + NCGTO1
         END IF
C
C
C             ...the new A and B shells are set. Calculate the
C                following info:
C
C                1) control variable to determine contraction order
C                2) cartesian monomial dimensions
C                3) spherical dimensions (= cartesian, if no spherical)
C                4) the overall norm factor due to s- or p-type shells
C
C
         SWAPRS = NPGTOA .GT. NPGTOB

         NXYZA  = (SHELLA+1)*(SHELLA+2)/2
         NXYZB  = (SHELLB+1)*(SHELLB+2)/2
         NXYZT  = NXYZA * NXYZB

         NRYA = SHELLA + SHELLA + 1
         NRYB = SHELLB + SHELLB + 1

         IF (.NOT.SPHERIC) THEN
             NRYA = NXYZA
             NRYB = NXYZB
         END IF

         SPNORM = ONE
         IF (SHELLA.EQ.1) THEN
             SPNORM = SPNORM + SPNORM
         END IF
         IF (SHELLB.EQ.1) THEN
             SPNORM = SPNORM + SPNORM
         END IF
C
C
C             ...calculate the coordinate differences between centers
C                A and B and calculate the square of the magnitude
C                of the distance between A and B. Set also the A,B
C                equality index and determine if differentiation is
C                is to be performed on centers A and/or B and x,y,z-
C                component. Total differentiation orders for each
C                component are also determined.
C
C
         ABX = XA - XB
         ABY = YA - YB
         ABZ = ZA - ZB
         RNABSQ = ABX * ABX + ABY * ABY + ABZ * ABZ

         NDERX = DERAX + DERBX
         NDERY = DERAY + DERBY
         NDERZ = DERAZ + DERBZ

         DIFFA = (DERAX + DERAY + DERAZ) .NE. 0
         DIFFB = (DERBX + DERBY + DERBZ) .NE. 0

         DIFFX = NDERX .NE. 0
         DIFFY = NDERY .NE. 0
         DIFFZ = NDERZ .NE. 0
C
C
C             ...ready!
C
C
         RETURN
         END
