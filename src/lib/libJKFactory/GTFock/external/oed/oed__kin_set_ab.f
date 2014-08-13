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
         SUBROUTINE  OED__KIN_SET_AB
     +
     +                    ( NCGTO1,NCGTO2,
     +                      NPGTO1,NPGTO2,
     +                      SHELL1,SHELL2,
     +                      X1,Y1,Z1,X2,Y2,Z2,
     +                      EXP1,EXP2,
     +                      CC1,CC2,
     +                      SPHERIC,
     +
     +                            NCGTOA,NCGTOB,
     +                            NPGTOA,NPGTOB,
     +                            SHELLA,SHELLB,SHELLP,
     +                            MXSHELL,
     +                            XA,YA,ZA,XB,YB,ZB,
     +                            ATOMIC,EQUALAB,
     +                            ABX,ABY,ABZ,RNABSQ,
     +                            SPNORM,
     +                            NXYZA,NXYZB,NXYZT,
     +                            NRYA,NRYB,
     +                            INDEXA,INDEXB,
     +                            SWAP12,SWAPRS,
     +                            LEXPA,LEXPB,
     +                            LCCA,LCCB,
     +                            LCCSEGA,LCCSEGB,
     +                            EMPTY )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__KIN_SET_AB
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine handles the logistics on how to evaluate
C                the (1|2) kinetic integral batch in the most efficient
C                way. It performs the label map:
C
C                               (1|2) --> (A|B)
C
C                according to certain criteria that have to be met for
C                efficiency. The freedom we have in making the internal
C                association 1,2 -> A,B follows from the 2-fold
C                permutational symmetry of the kinetic integrals in
C                (1|2):
C
C                                 (1|2) = (2|1)
C
C                If the 1 <-> 2 switch has to be applied is simply
C                governed by the demand that the final A and B shell
C                labels have to obey the relation A>=B, since this means
C                the least amount of space and work needed to generate
C                the necessary 1D overlap and kinetic integrals.
C
C
C                  Input (x = 1 and 2):
C
C                   NCGTOx        =  # of contractions for csh x
C                   NPGTOx        =  # of primitives per contraction
C                                    for csh x
C                   SHELLx        =  the shell type for csh x
C                   Xy,Yy,Zy      =  the x,y,z-coordinates for centers
C                                    y = 1 and 2
C                   EXPx          =  primitive exponents for csh x
C                   CCx           =  contraction coeffs for csh x
C                   SPHERIC       =  is true, if spherical integrals
C                                    are wanted, false if cartesian
C                                    ones are wanted
C
C                  Output (x = A and B):
C
C                   NCGTOx        =  # of contractions for csh x
C                   NPGTOx        =  # of primitives per contraction
C                                    for csh x
C                   SHELLx        =  the shell type for csh x
C                   SHELLP        =  the shell sum A+B
C                   MXSHELL       =  the largest (maximum) shell type
C                   Xy,Yy,Zy      =  the x,y,z-coordinates for centers
C                                    y = A and B
C                   ATOMIC        =  indicates, if purely atomic
C                                    integrals will be evaluated
C                   EQUALAB       =  indicates, if csh x and csh y are
C                                    considered to be equal for the
C                                    pair AB
C                   ABX,ABY,ABZ   =  the x,y,z-coordinate differences
C                                    between centers A and B
C                   RNABSQ        =  square of the magnitude of the
C                                    distance between centers A and B
C                   SPNORM        =  normalization factor due to
C                                    presence of s- and p-type shells.
C                                    For each s-type shell there is a
C                                    factor of 1 and for each p-type
C                                    shell a factor of 2
C                   NXYZx         =  # of cartesian monomials for csh x
C                   NXYZT         =  total # cartesian monomial
C                                    combinations, that is NXYZA * NXYZB
C                   NRYx          =  # of spherical functions for csh x
C                   INDEXx        =  index A,B -> 1,2 map
C                   SWAP12        =  is .true., if a swap 1 <-> 2 has
C                                    been performed
C                   SWAPRS        =  is set .true. if the contraction
C                                    order of the primitives pair AB
C                                    will be performed in reverse order
C                                    BA for efficiency reasons
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

         LOGICAL     ATOMIC
         LOGICAL     CASE1,CASE2
         LOGICAL     EMPTY
         LOGICAL     EQUALAB
         LOGICAL     SPHERIC
         LOGICAL     SWAP12,SWAPRS

         INTEGER     I,J
         INTEGER     INDEXA,INDEXB
         INTEGER     LCCSEGA,LCCSEGB
         INTEGER     LCCA,LCCB
         INTEGER     LEXPA,LEXPB
         INTEGER     MXSHELL
         INTEGER     NCGTO1,NCGTO2
         INTEGER     NCGTOA,NCGTOB
         INTEGER     NPGTO1,NPGTO2
         INTEGER     NPGTOA,NPGTOB
         INTEGER     NRYA,NRYB
         INTEGER     NXYZA,NXYZB,NXYZT
         INTEGER     SHELL1,SHELL2
         INTEGER     SHELLA,SHELLB,SHELLP

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
C                if a zero batch of kinetic integrals is expected.
C
C
         EMPTY  = .FALSE.
         ATOMIC = (X1.EQ.X2) .AND. (Y1.EQ.Y2) .AND. (Z1.EQ.Z2)
         SHELLP = SHELL1 + SHELL2

         MXSHELL = MAX0 (SHELL1,SHELL2)

         CASE1 = MOD (SHELLP,2) .EQ. 1
         CASE2 = SPHERIC. AND. (SHELL1.NE.SHELL2)

         IF (ATOMIC .AND. (CASE1.OR.CASE2)) THEN
             EMPTY = .TRUE.
             RETURN
         END IF
C
C
C             ...determine csh equality between centers 1 and 2
C                in increasing order of complexity:
C
C                 centers -> shells -> exponents -> ctr coefficients
C
C
         EQUALAB = ATOMIC

         IF (EQUALAB) THEN
             EQUALAB =     (SHELL1 .EQ. SHELL2)
     +               .AND. (NPGTO1 .EQ. NPGTO2)
     +               .AND. (NCGTO1 .EQ. NCGTO2)
             IF (EQUALAB) THEN
               DO I = 1,NPGTO1
                  EQUALAB = EQUALAB .AND. (EXP1(I).EQ.EXP2(I))
               END DO
               IF (EQUALAB) THEN
                 DO J = 1,NCGTO1
                    IF (EQUALAB) THEN
                      DO I = 1,NPGTO1
                         EQUALAB = EQUALAB .AND. (CC1(I,J).EQ.CC2(I,J))
                      END DO
                    END IF
                 END DO
               END IF
             END IF
         END IF
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
C                contraction coefficients. Also set the info for
C                evaluation of the [a|b] kinetic batches.
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

         NXYZA = (SHELLA+1)*(SHELLA+2)/2
         NXYZB = (SHELLB+1)*(SHELLB+2)/2
         NXYZT = NXYZA * NXYZB

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
C                of the distance between A and B.
C
C
         IF (.NOT.ATOMIC) THEN
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
C
C
C             ...ready!
C
C
         RETURN
         END
