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
         SUBROUTINE  OED__XYZ_TO_RY_MATRIX
     +
     +                    ( NXYZ,NRY,NROWMX,
     +                      L,
     +                      TEMP,
     +
     +                            NROW,
     +                            ROW,
     +                            TMAT )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__XYZ_TO_RY_MATRIX
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation generates the transformation matrix for
C                performing a primitive gaussian integral batch
C                transformation from cartesian to spherical gaussian
C                type orbitals.
C
C                Input:
C
C                     NXYZ = dimension of xyz-monomial basis
C                            corresponding to the specified shell
C                            quantum number L.
C                            Must at least be equal to (L+1)*(L+2)/2.
C
C                      NRY = dimension of ry-spherical basis
C                            corresponding to the specified shell
C                            quantum number L.
C                            Must at least be equal to 2L+1.
C
C                   NROWMX = maximum number of xyz-monomials in one
C                            ry-solid harmonic function corresponding
C                            to the specified shell quantum number L.
C                            Must at least be equal to (L/2+1)*(L/2+2)/2.
C
C                        L = Shell quantum number.
C
C                     TEMP = temporary scratch array.
C
C
C                Output:
C
C                     NROW (I) = Number of xyz-monomials contributing
C                                to the I-th ry-component.
C
C                    ROW (K,I) = K-th xyz-monomial position number
C                                containing nonzero contribution to
C                                the I-th ry-component.
C
C                   TMAT (K,I) = K-th nonzero xyz-monomial coefficient
C                                to the I-th ry-component.
C
C
C                Each ry-component for a specific shell type corresponds
C                to a sherical harmonic function Y (L,M) with M ranging
C                from -L to +L. The TMAT matrix ry-component indices are
C                ordered such that +L corresponds to the 1st index and
C                -L to the last one.
C
C                For the ordering of the xyz-monomials used in the
C                expansion of the ry-components we use the usual
C                preference ordering p>q>r in the x-,y- and z-exponents
C                of x^p y^q z^r. For L=3 we would thus have the 10
C                monomials arranged as follows:
C
C                                #  |  x  y  z
C                              -----------------
C                                1  |  3  0  0
C                                2  |  2  1  0
C                                3  |  2  0  1
C                                4  |  1  2  0
C                                5  |  1  1  1
C                                6  |  1  0  2
C                                7  |  0  3  0
C                                8  |  0  2  1
C                                9  |  0  1  2
C                               10  |  0  0  3
C
C                This ordering is also assumed when evaluating batches
C                of cartesian gaussian integrals, hence the present
C                routine and the cartesian gaussian integral generation
C                routines go hand in hand and must be viewed together
C                when performing xyz-monomial basis ordering changes!
C
C                The ultimate goal of using the present routine is
C                to achieve achieve a transformation matrix T relating
C                cartesian gaussian functions to spherical ones.
C                At his stage it is convenient to introduce all
C                factors independent! of the gaussian exponents that
C                enter the final normalization constant for each
C                spherical gaussian function. Such a function is given
C                by the following expression:
C
C                    LM                     L    LM                2
C                 GTO  (r,t,p) = N (L,a) * r  * Y  (t,p) * exp (-ar )
C                    a
C
C                where t = theta and p = phi and N (L,a) denotes a
C                normalization factor such that
C
C                              LM       *    LM
C                 integral {GTO  (r,t,p)  GTO  (r,t,p) dr dt dp} = 1
C                              a             a
C
C                Since the spherical harmonics are normalized, we can
C                integrate out over the angles t and p and arrive at
C                the following defining expression for N (L,a):
C
C                        2  r=oo  2L+2          2
C                 N (L,a)   int  r     exp (-2ar ) dr  =  1
C                           r=0
C
C                which leads to:
C
C                               ____________________________
C                              / 2^(2L+3+1/2) * a^((2L+3)/2)
C                 N (L,a) =   / -----------------------------
C                           \/     (2L+1)!! * sqrt (pi)
C
C
C                The following part of N (L,a) has already been
C                taken care of during the cartesian integrals
C                generation:
C
C                               ____________________________
C                              / 2^(1+1/2) * a^((2L+3)/2)
C                 N (L,a) =   / -----------------------------
C                           \/          sqrt (pi)
C
C
C                hence we are left only with the following part
C                of N (L,a) to deal with here:
C
C                                   __________
C                                  / 2^(2L+2)
C                                 / -----------
C                               \/   (2L+1)!!
C
C                which is designated RNORM (to reflect its origin)
C                and is being calculated in logarithmic form to avoid
C                integer overflow for the double factorial. The
C                normalization constant for the spherical harmonic
C                part, depending on both L and M quantum numbers,
C                is designated YNORM and is evaluated stepwise
C                from its base value at M = -L or L downwards to
C                -1 or 0. Note, that YNORM should also include a factor
C                in terms of pi, namely 1/sqrt(pi), but this has also
C                already been incorporated into the cartesian integrals
C                evaluation routine.
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

         INTEGER    I,K,L,M,P,T,V
         INTEGER    KEND,PEND
         INTEGER    PP,TT,VV
         INTEGER    LMM,LPM
         INTEGER    NRY,NXYZ,NROWMX
         INTEGER    RY,XYZ

         INTEGER    NROW (1:NRY)

         INTEGER    ROW (1:NROWMX,1:NRY)

         DOUBLE PRECISION   A,B,C,D,X
         DOUBLE PRECISION   RNORM,YNORM,RYNORM
         DOUBLE PRECISION   XCOL
         DOUBLE PRECISION   XM,XK,XL,XT,XV
         DOUBLE PRECISION   XKK,XLL,XPP
         DOUBLE PRECISION   XLMM,XLPM
         DOUBLE PRECISION   ZERO,HALF,ONE,TWO,TEN

         DOUBLE PRECISION   TEMP (1:NXYZ)

         DOUBLE PRECISION   TMAT (1:NROWMX,1:NRY)

         PARAMETER  (ZERO    =  0.D0)
         PARAMETER  (HALF    =  0.5D0)
         PARAMETER  (ONE     =  1.D0)
         PARAMETER  (TWO     =  2.D0)
         PARAMETER  (TEN     =  10.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...calculate the RNORM and the initial value of
C                YNORM at M = -L. Form initial combined norm RYNORM.
C
C
         RNORM = DFLOAT (2*L+2) * DLOG10 (TWO)
         DO 10 I = 1,L
            RNORM = RNORM - DLOG10 (DFLOAT (I+I+1))
   10    CONTINUE
         RNORM = TEN ** (HALF*RNORM)

         YNORM = HALF
         DO 20 I = 1,L
            YNORM = YNORM * DFLOAT (I+I+1) / DFLOAT (I+I)
   20    CONTINUE
         YNORM = DSQRT (YNORM)

         RYNORM = RNORM * YNORM
C
C
C             ...evaluate monomial expansion for spherical harmonics
C                Y(L,M) for the cases M < 0. Note that the loop below
C                uses M, which goes from L down to 1. This M represents
C                in fact the absolute value of M, i.e. |M|.
C
C                There are two ways of coding this:
C
C                    1) symbolically, in which the nominators and
C                       denominators of all fractions are handled
C                       separately, performing common divisor cuts
C                       at every single step. This procedure is
C                       slower than the direct evaluation but more
C                       accurate. Also the symbollic procedure suffers
C                       from integer overflow at about L = 16, so
C                       there is a limit on the shell size one can
C                       safely handle.
C
C                    2) the direct mode, in which all fractions are
C                       multiplied together as flp's.
C
C                This routine uses the direct mode.
C
C
         XL = DFLOAT (L)
         XLL = XL + XL

         DO 100 M = L,1,-1

            DO 110 XYZ = 1,NXYZ
               TEMP (XYZ) = ZERO
  110       CONTINUE

            LMM = L - M
            LPM = L + M
            XM = DFLOAT (M)
            XLPM = DFLOAT (LPM)
            XLMM = DFLOAT (LMM)
            KEND = LMM / 2
            PEND = (M - 1) / 2
            
            IF (M.EQ.L) THEN
                A = XL
            ELSE
                IF (MOD(LMM,2).EQ.0) THEN
                    X = XLMM / (XLPM + ONE)
                    RYNORM = RYNORM * DSQRT (X)
                    A = A * (XM/(XM+ONE)) / X
                ELSE
                    RYNORM = RYNORM * DSQRT (XLMM * (XLPM + ONE))
                    A = A * (XM/(XM+ONE)) / XLMM
                END IF
            END IF

            B = A
            DO 120 K = 0,KEND
               IF (K.NE.0) THEN
                   XK = DFLOAT (K)
                   XKK = XK + XK
                   X = XLMM - XKK + ONE
                   B = - B * (X/(XLL-XKK+ONE)) * ((X+ONE)/TWO)
               END IF

               C = B
               DO 130 P = 0,PEND
                  PP = P + P
                  IF (P.NE.0) THEN
                      XPP = DFLOAT (PP)
                      X = XM - XPP
                      C = - C * (X/XPP) * ((X+ONE)/(XPP+ONE))
                  END IF

                  D = C
                  DO 140 T = K,0,-1
                     TT = T + T
                     XT = DFLOAT (T)
                     IF (T.NE.K) THEN
                         D = D / (XK - XT)
                     END IF

                     XCOL = D
                     DO 150 V = 0,T
                        VV = V + V
                        IF (V.EQ.0) THEN
                            DO 160 I = 2,T
                               XCOL = XCOL / DFLOAT (I)
  160                       CONTINUE
                        ELSE
                            XV = DFLOAT (V)
                            XCOL = XCOL * ((XT-XV+ONE)/XV)
                        END IF

                        XYZ = (LMM-VV+PP+2)*(LMM-VV+PP+1)/2 + LMM-TT+1

                        TEMP (XYZ) = TEMP (XYZ) + XCOL

  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE

            RY = LPM+1

            I = 0
            DO 170 XYZ = 1,NXYZ
               IF (TEMP (XYZ).NE.ZERO)  THEN
                   I = I + 1
                   ROW (I,RY) = XYZ
                   TMAT (I,RY) = RYNORM * TEMP (XYZ)
               END IF
  170       CONTINUE
            NROW (RY) = I
  100    CONTINUE
C
C
C             ...evaluate monomial expansion for spherical harmonics
C                Y(L,M) for the cases M >= 0.
C
C
         RYNORM = RNORM * YNORM

         DO 200 M = L,0,-1

            DO 210 XYZ = 1,NXYZ
               TEMP (XYZ) = ZERO
  210       CONTINUE

            LMM = L - M
            LPM = L + M
            XM = DFLOAT (M)
            XLPM = DFLOAT (LPM)
            XLMM = DFLOAT (LMM)
            KEND = LMM / 2
            PEND = M / 2
            
            IF (M.EQ.0) THEN
                RYNORM = RYNORM / DSQRT (TWO)
            END IF

            IF (M.EQ.L) THEN
                A = ONE
            ELSE
                IF (MOD(LMM,2).EQ.0) THEN
                    X = XLMM / (XLPM + ONE)
                    RYNORM = RYNORM * DSQRT (X)
                    A = A / X
                ELSE
                    RYNORM = RYNORM * DSQRT (XLMM * (XLPM + ONE))
                    A = A / XLMM
                END IF
            END IF

            B = A
            DO 220 K = 0,KEND
               IF (K.NE.0) THEN
                   XK = DFLOAT (K)
                   XKK = XK + XK
                   X = XLMM - XKK + ONE
                   B = - B * (X/(XLL-XKK+ONE)) * ((X+ONE)/TWO)
               END IF

               C = B
               DO 230 P = 0,PEND
                  PP = P + P
                  IF (P.NE.0) THEN
                      XPP = DFLOAT (PP)
                      X = XM - XPP + ONE
                      C = - C * (X/XPP) * ((X+ONE)/(XPP-ONE))
                  END IF

                  D = C
                  DO 240 T = K,0,-1
                     TT = T + T
                     XT = DFLOAT (T)
                     IF (T.NE.K) THEN
                         D = D / (XK - XT)
                     END IF

                     XCOL = D
                     DO 250 V = 0,T
                        VV = V + V
                        IF (V.EQ.0) THEN
                            DO 260 I = 2,T
                               XCOL = XCOL / DFLOAT (I)
  260                       CONTINUE
                        ELSE
                            XV = DFLOAT (V)
                            XCOL = XCOL * ((XT-XV+ONE)/XV)
                        END IF

                        XYZ = (LMM-VV+PP+1)*(LMM-VV+PP)/2 + LMM-TT+1
 
                        TEMP (XYZ) = TEMP (XYZ) + XCOL

  250                CONTINUE
  240             CONTINUE
  230          CONTINUE
  220       CONTINUE

            RY = LMM+1

            I = 0
            DO 270 XYZ = 1,NXYZ
               IF (TEMP (XYZ).NE.0)  THEN
                   I = I + 1
                   ROW (I,RY) = XYZ
                   TMAT (I,RY) = RYNORM * TEMP (XYZ)
               END IF
  270       CONTINUE
            NROW (RY) = I
  200    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
