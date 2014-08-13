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
         SUBROUTINE  OED__RYS_X_ROOTS_WEIGHTS
     +
     +                    ( NT,NTGQP,
     +                      NGQP,NMOM,
     +                      TVAL,RYSZERO,
     +                      A,B,
     +                      MOM,
     +                      DIA,OFF,
     +                      ROW1,ROW2,
     +
     +                                RTS,
     +                                WTS )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__RYS_X_ROOTS_WEIGHTS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : Master routine to evaluate roots and weights in the
C                interval [0,1] over the Rys weight function:
C
C
C                                       exp(-T*x)
C                             W   (x) = ---------
C                              Rys      2*sqrt(x)
C
C
C                using the general Gaussian Quadrature technique
C                consisting of the following basic steps:
C
C                   1) Calculate the auxillary polynomial moments
C                   2) Calculate the auxillary polynomial coefficients
C                   3) Set up the tridiagonal symmetric terminal
C                      matrix
C                   4) Solve the tridiagonal symmetric terminal
C                      matrix for roots and weights.
C
C                In this routine all T-values are treated all at once.
C                The value of the T-parameter dictates which type of
C                auxillary set of polynomials is to be used for the
C                modified Chebyshev algorithm. The auxillary polynomials
C                are often chosen (as is the case here) to be orthogonal
C                relative to some classical weight function.
C
C                The classical weights used to establish the auxillary
C                polynomials are the following (notation according to
C                M.Abramowitz and I.A.Stegun, Handbook of Mathematical
C                Functions, 1964):
C
C
C                i) Range of validity: 0 =< T =< 30
C
C                   Here we use shifted Jacobi weights and polynomials:
C
C
C                        W      (p,q,x)  =  (1-x)^(p-q) * x^(q-1)
C                         Jacobi
C
C                        i-th Jacobi polynomial  =  G (p,q,x)
C                                                    i
C
C                   with conditions: x in interval [0,1]
C                                    p-q > -1 , q > 0
C                                    i = 0,1,2,...
C
C
C                ii) Range of validity: 1 < T =< oo (infinity)
C
C                   Here we use generalized Laguerre weights and
C                   polynomials:
C
C
C                        W        (a,x)  =  exp^(-x) * x^a
C                         Laguerre
C
C                        i-th Laguerre polynomial  =  L (a,x)
C                                                      i
C
C                   with conditions: x in interval [0,inf]
C                                    a > -1
C                                    i = 0,1,2,...
C
C
C                Range of validity means that for all the specified
C                T's within the range the resulting moment integrals
C                are accurate to within 1.D-16.
C
C
C                  Input:
C
C                    NT           =  # of T-values
C                    NTGQP        =  # of roots times # of T-values
C                    NGQP         =  # of gaussian quadrature points
C                                    (roots)
C                    NMOM         =  # of necessary moment integrals
C                                    to calculate the quadrature roots
C                    TVAL         =  the T-values
C                    RYSZERO      =  the zeroth Rys moments for all
C                                    T-values
C                    A,B          =  will contain the recurrence
C                                    coefficients for the auxillary
C                                    polynomials
C                    MOM          =  will contain the normed auxillary
C                                    polynomial modified moments
C                    DIA,OFF      =  will contain the diagonal and
C                                    offdiagonal elements of the
C                                    tridiagonal symmetric terminal
C                                    matrix
C                    ROW1,ROW2    =  will be used to evaluate the
C                                    tridiagonal elements of the
C                                    symmetric terminal matrix in an
C                                    efficient way using Sack and
C                                    Donovan's method
C
C                  Output:
C
C                    RTS          =  the roots array
C                    WTS          =  the weights array
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

         INCLUDE    'oed__jacobi.inc'

         INTEGER    I,J,M,N
         INTEGER    IMAX,JMAX
         INTEGER    IP1,JP1
         INTEGER    ITER
         INTEGER    MOMMAX
         INTEGER    MXITER
         INTEGER    NGQP
         INTEGER    NMOM
         INTEGER    NRTS
         INTEGER    NT,NTGQP

         DOUBLE PRECISION   C,D,F,G,P,R,S
         DOUBLE PRECISION   BINC,SINC
         DOUBLE PRECISION   LIM1,LIM2,LIM3
         DOUBLE PRECISION   MOMI,MOMIM1,MOMIP1,MOMZERO
         DOUBLE PRECISION   R1
         DOUBLE PRECISION   ROOT
         DOUBLE PRECISION   SCALE
         DOUBLE PRECISION   SIGMA,THETA
         DOUBLE PRECISION   T,TINV,TINV2,TINVHF,TINVSQ,TEXP,TPOWER
         DOUBLE PRECISION   TEST1,TEST2
         DOUBLE PRECISION   ZMOM,ZINV
         DOUBLE PRECISION   ZERO,EXTREM,VSMALL,SMALL,TENTH,
     +                      HALF,ONE,ONEP5,TWO,THREE,TEN,TLIMIT

         DOUBLE PRECISION   A       (1:NMOM)
         DOUBLE PRECISION   B       (2:NMOM)
         DOUBLE PRECISION   MOM     (1:NMOM)
         DOUBLE PRECISION   DIA     (1:NGQP)
         DOUBLE PRECISION   OFF     (1:NGQP)
         DOUBLE PRECISION   ROW1    (1:NMOM)
         DOUBLE PRECISION   ROW2    (1:NMOM-1)

         DOUBLE PRECISION   RYSZERO (1:NT)
         DOUBLE PRECISION   TVAL    (1:NT)

         DOUBLE PRECISION   RTS     (1:NTGQP)
         DOUBLE PRECISION   WTS     (1:NTGQP)

         PARAMETER  (ZERO    =  0.D0)
         PARAMETER  (EXTREM  =  1.0D-300)
         PARAMETER  (VSMALL  =  1.0D-16)
         PARAMETER  (SMALL   =  1.0D-6)
         PARAMETER  (TENTH   =  0.1D0)
         PARAMETER  (HALF    =  0.5D0)
         PARAMETER  (ONE     =  1.D0)
         PARAMETER  (ONEP5   =  1.5D0)
         PARAMETER  (TWO     =  2.D0)
         PARAMETER  (THREE   =  3.D0)
         PARAMETER  (TEN     =  10.D0)
         PARAMETER  (TLIMIT  =  15.D0)

         PARAMETER  (MOMMAX  =  30)
         PARAMETER  (MXITER  =  30)
C
C
C------------------------------------------------------------------------
C
C
C             ...main loop over all T values. Check which T case
C                applies.
C
C
         NRTS = 0

         DO 1000 N = 1,NT

            T = TVAL (N)
            MOMZERO = RYSZERO (N)

            IF (T.LE.TLIMIT) THEN
C
C
C             ...The Jacobi section. Check first, if the number of
C                moments wanted is within limits for the Jacobi case.
C                If ok, we proceed by calculating two steps at the
C                same time:
C
C                1) normed shifted Jacobi modified moments defined as
C                   integrals on [0,1] over the Rys weight function:
C
C                                        exp(-Tx)
C                             W   (x) = ---------
C                              Rys      2*sqrt(x)
C
C                   and shifted Jacobi polynomials:
C
C                                G (0.5,0.5,x)
C                                 i
C
C                   A 3-term recursion relation can be given for these
C                   moments as follows:
C
C                       MOM (i+1) = R * MOM (i)  +  S * MOM (i-1)
C
C                   where
C
C                      R  =  (2i+1)/2T  +  (2i+1)/(4i-1)(4i+3)
C                      S  =  2i(2i+1)(2i-1)**2 / (4i-3)(4i+1)(4i-1)**2
C
C                   However, in order to evaluate the moments to the
C                   required accuracy, a downward recursion must be
C                   applied starting with some tiny seed values for
C                   MOM (i+1) and MOM (i):
C
C                       MOM (i-1) = 1/S * ( MOM (i+1) - R * MOM (i))
C
C                   The sequence of moments thus obtained is not
C                   normalized. Hence the downward recursion is
C                   performed down to the zeroth moment ZMOM and the
C                   relevant moments are normalized by subsequent
C                   division with ZMOM.
C
C                   The above downward recursion technique is extremely
C                   unstable if T is very small, in which case even the
C                   small seed values are not small enough to prevent
C                   overflow of the lowest not normalized moments.
C                   In such a case it is better to build the normed
C                   sequence starting with the first moment after the
C                   following considerations. An expansion of MOM (i)
C                   in terms of powers of T shows that the leading
C                   term is in the i-th power of T, i.e.
C
C                                    inf           k
C                       MOM (i)  =   sum  c (i,k) T
C                                   k = i
C
C                   with leading coefficient c (i,i) getting smaller
C                   with increasing i. Hence for very small T it is
C                   sufficient to set:
C
C                                             i
C                       MOM (i)  =   c (i,i) T     (for T very small)
C
C
C                   an equation which is valid in terms of computer
C                   accuray if T is less or equal to the minimum
C                   possible nonzero number within the precision of
C                   the mantissa allowed (here in the present code
C                   this is double precision, hence T less or equal
C                   1.D-16). Rather than evaluating each leading
C                   coefficient c (i,i) individualy (they are given
C                   by very! complicated expressions in terms of sums
C                   of double factorials), the procedure adopted here
C                   was to predetermine each c (i,i) numerically by
C                   running the routine with T = 1.D-16 and presetting
C                   the resulting coefficients in a data array CSMALL
C                   to be used whenever T =< 1.D-16.
C
C                   If T is in the range 1.D-16 to TMAX = 30.D0, then
C                   the routine was calibrated such that the mantissa
C                   of the resulting moments are accurate to the 16th
C                   decimal place. Calibration in this context means the
C                   predetermination of the maximum moment (controled
C                   in the code by its index IMAX) that has to be
C                   generated to perform the downward recurrence
C                   relation starting with the tiny seeds. Hence the
C                   calibration depends on four things:
C
C                       i) The total maximum number MOMMAX of normed
C                          moments wanted.
C
C                      ii) The mantissa accuray wanted for the normed
C                          moments.
C
C                     iii) The tiny values of the two seeds (the nonzero
C                          seed should correpond to the smallest
C                          possible nonzero number representable on
C                          the present computer).
C
C                      iv) T range. Each T range requires a different
C                          maximum moment (i.e. a different IMAX value)
C                          to start with.
C
C                   In order to perform a calibration for a different
C                   machine, one has to write a separate small program
C                   using this routine with different seeds, T- and IMAX
C                   values. Once the seeds have been set, this little
C                   program should recalculate the MOMMAX normed
C                   moments for increasing IMAX values starting with
C                   IMAX= MOMMAX. The IMAX value finally taken for a
C                   particular T-value should then obey the following
C                   inequality:
C
C
C                   | MOM (i,IMAX) - MOM (i,IMAX-1)|
C                   | -----------------------------| < mantissa accuracy
C                   |         MOM (i,IMAX-1)       |
C
C
C                       COMMENTS FOR FAST EVALUATION OF THE MOMENTS
C                      ---------------------------------------------
C                   As can be seen from the downward recursion formula
C                   for the moments, all we need are the values of the R
C                   and 1/S recursion parameters, which are conveniently
C                   precomputed (with R being split as R = R1/2T + R2)
C                   and supplied via an include table 'oed__jacobi.inc'.
C                   Hence 'oed__jacobi.inc' will have the values of the
C                   complicated R2 and SINV = 1/S expressions in the
C                   range from i = 1 to 100. The R1 values are simply
C                   calculated by using the appropriate decrement value
C                   of 2. The include file 'oed__jacobi.inc' will also
C                   contain the CSMALL array for the moments of very
C                   small T's.
C
C
C                2) recurrence coefficients for the shifted Jacobi
C                   polynomials G (0.5,0.5,x), denoted simply by G (x):
C
C
C                         G (x)  =  1
C                          0
C
C                         G (x)  =  (x - A )
C                          1              1
C
C                         G   (x)  =  (x - A   ) * G (x) - B   G   (x)
C                          i+1              i+1     i       i+1 i-1
C
C
C                   The result consists of the recurrence coefficients:
C
C                               A (I) , I = 1,NMOM
C                               B (I) , I = 2,NMOM
C
C                   whose values are given by the following expressions:
C
C
C                        A (i+1) = 4i*(2i+1)-1 / (4i+3)(4i-1)
C
C                        B (i+1) = (4i**2)((4i**2)-4i+1) /
C                                  (4i-3)(4i+1)((4i-1)**2)
C
C                   Since these are complicated and time consuming
C                   expressions, they are precalculated and included
C                   into the include file 'oed__jacobi.inc'.
C
C
C
                IF (NMOM.GT.MOMMAX) THEN
                    WRITE (*,*) ' # of Jacobi moments exceeded! '
                    WRITE (*,*) ' NMOM = ',NMOM
                    WRITE (*,*) ' Maximum value = ',MOMMAX
                    WRITE (*,*) ' oed__rys_x_roots_weights '
                    STOP
                END IF
C
C
C             ...the very small T case.
C
C
                IF (T .LE. VSMALL) THEN

                    IMAX = MIN0 (NMOM,16)
                    A (1) = AJAC (1)
                    MOM (1) = CSMALL (1) * T

                    TPOWER = T
                    DO I = 2,IMAX
                       TPOWER = TPOWER * T
                       A (I) = AJAC (I)
                       B (I) = BJAC (I)
                       MOM (I) = CSMALL (I) * TPOWER
                    END DO

                    DO I = IMAX+1,NMOM
                       A (I) = AJAC (I)
                       B (I) = BJAC (I)
                       MOM (I) = ZERO
                    END DO

                ELSE
C
C
C             ...the general Jacobi case. Set maximum number of moments
C                necessary to get required moments to an accuracy of
C                at least 1.D-16. See above in the routine description
C                for details and calibration of the setting.
C
C
                    IF (NMOM.LE.5) THEN
                        IF (T.LT.SMALL) THEN
                            IMAX = NMOM + 1
                        ELSE IF (T.LT.TENTH) THEN
                            IMAX = NMOM + 3
                        ELSE IF (T.LT.TWO) THEN
                            IMAX = NMOM + 7
                        ELSE IF (T.LT.TEN) THEN
                            IMAX = NMOM + 13
                        ELSE
                            IMAX = NMOM + 22
                        END IF
                    ELSE
                        IF (T.LT.SMALL) THEN
                            IMAX = NMOM
                        ELSE IF (T.LT.TENTH) THEN
                            IMAX = NMOM + 2
                        ELSE IF (T.LT.TWO) THEN
                            IMAX = NMOM + 4
                        ELSE IF (T.LT.TEN) THEN
                            IMAX = NMOM + 8
                        ELSE
                            IMAX = NMOM + 16
                        END IF
                    END IF
C
C
C             ...proceed by setting seed values for downward recursion
C                to obtain minimal solution and start recursive
C                evaluation for all Jacobi moments down to first moment
C                (two loops are necessary here due to calculation of
C                higher moments than actually needed later on).
C
C
                    MOMI = EXTREM
                    MOMIP1 = ZERO

                    TINVHF = HALF / T
                    R1 = DFLOAT (2*IMAX+5)

                    DO I = IMAX+1,NMOM+2,-1
                       R1 = R1 - TWO
                       R = R1 * TINVHF  +  R2 (I)
                       MOMIM1  =  SINV (I) * (MOMIP1 - R * MOMI)
                       MOMIP1 = MOMI
                       MOMI = MOMIM1
                    END DO

                    DO I = NMOM+1,2,-1
                       R1 = R1 - TWO
                       R = R1 * TINVHF  +  R2 (I)
                       MOMIM1  =  SINV (I) * (MOMIP1 - R * MOMI)
                       MOM (I-1) = MOMIM1
                       MOMIP1 = MOMI
                       MOMI = MOMIM1
                    END DO
C
C
C             ...evaluate zeroth moment and normalize sequence.
C                If the absolute zeroth moment is less than the
C                approximate absolute nonzero minimum (set here
C                equal to 1.D-300), the normalization looses its
C                meaning and the calculation must be stopped.
C                Set also here the recurrence relation coefficients
C                A and B for the shifted Jacobi polynomials.
C
C
                    R = THREE * TINVHF  +  R2 (1)
                    ZMOM  =  SINV (1) * (MOMIP1 - R * MOMI)

                    IF (DABS (ZMOM) .LT. EXTREM) THEN
                        WRITE (*,*) ' 0th Jacobi moment found to be 0! '
                        WRITE (*,*) ' ZMOM = ',ZMOM
                        WRITE (*,*) ' Recalibrate maximum settings! '
                        WRITE (*,*) ' oed__rys_x_roots_weights '
                        STOP
                    END IF

                    A (1) = AJAC (1)
                    ZINV = ONE / ZMOM
                    MOM (1) = MOM (1) * ZINV

                    DO I = 2,NMOM
                       A (I) = AJAC (I)
                       B (I) = BJAC (I)
                       MOM (I) = MOM (I) * ZINV
                    END DO

                END IF

            ELSE
C
C
C             ...The Laguerre section. As for the Jacobi case,
C                we calculate two steps at the same time:
C
C                1) normed Laguerre modified moments defined as
C                   integrals on [0,1] over the Rys weight function:
C
C                                        exp(-Tx)
C                             W   (x) = ---------
C                              Rys      2*sqrt(x)
C
C                   and monic generalized 'scaled' Laguerre polynomials:
C
C                                L (-0.5,Tx)
C                                 i
C
C                   A closed formula can be given for these moments in
C                   terms of monic generalized 'scaled' Laguerre
C                   polynomials with generalization +0.5:
C
C                       MOM (i) = SCALE * L   (+0.5,T)
C                                          i-1
C                   where
C
C                               - exp(-T)
C                       SCALE = ---------   ;  F0(T) = Rys zero moment
C                                2T*F0(T)
C
C                   The recursion relation for the +0.5 polynomials is:
C
C
C                     L   (+0.5,T) = R * L   (+0.5,T) - S * L   (+0.5,T)
C                      i-1                i-2                i-3
C
C                   where
C
C                          R  =  (T - 2i + 5/2) / T
C                          S  =  (i - 2)(i - 3/2) / T*T
C
C
C                   All moments MOM (i);i=1,2,...NMOM using the above
C                   outlined algorithm are evaluated.
C
C
C                2) recurrence coefficients for T-scaled generalized
C                   monic Laguerre polynomials L (-0.5,Tx), denoted
C                   simply by L (Tx):
C
C
C                       L (Tx)    =   1
C                        0
C
C                       L (Tx)    =   (x - A )
C                        1                  1
C
C                       L   (Tx)  =   (x - A   ) * L (Tx) - B   L   (Tx)
C                        i+1                i+1     i        i+1 i-1
C
C
C                   The result consists of the recurrence coefficients:
C
C                             A (I) , I = 1,NMOM
C                             B (I) , I = 2,NMOM
C
C                   whose values are given by the following expressions:
C
C                           A (1) = 1 / 2T
C
C                         A (i+1) = (2i+1/2) / T
C
C                         B (i+1) = i*(i-1/2) / (T*T)
C
C
C
                TEXP = DEXP (-T)
                TINV = ONE / T
                TINV2 = TWO * TINV
                TINVHF = HALF * TINV
                TINVSQ = TINV * TINV

                SCALE  =  - TINVHF * TEXP / MOMZERO

                IF (NMOM.EQ.1) THEN
                    A (1) = TINVHF
                    MOM (1) = SCALE
                ELSE
                    A (1) = TINVHF
                    A (2) = TINVHF + TINV2
                    B (2) = HALF * TINVSQ
                    MOM (1) = SCALE
                    R = ONE - (ONEP5 * TINV)
                    MOM (2) = SCALE * R

                    S = ZERO
                    BINC = HALF
                    SINC = - HALF
                    LIM2 = R
                    LIM3 = ONE

                    DO I = 3,NMOM
                       BINC = BINC + TWO
                       A (I) = A (I-1) + TINV2
                       B (I) = B (I-1) + BINC * TINVSQ
                       SINC = SINC + TWO
                       R = R - TINV2
                       S = S + SINC * TINVSQ
                       LIM1 = R * LIM2 - S * LIM3
                       MOM (I) = SCALE * LIM1
                       LIM3 = LIM2
                       LIM2 = LIM1
                    END DO

                END IF

            END IF
C
C
C             ...This section calculates a symmetric terminal matrix
C                using the normed modified moments MOM and related
C                recurrence coefficients A and B of monic polynomials
C                established in the previous section. The algorithm is
C                a transcription of the LQMD algorithm described by Sack
C                and Donovan in Numer. Mathematik 18, 465-478 (1972).
C
C                The needed data is as follows (NMOM = 2*NGQP-1):
C
C                   1) Normed modified moments: MOM (I),I=1,NMOM
C
C                   2) Recurrence coefficients of the corresponding
C                      monic polynomials: A (I),I=1,NMOM and B (I),
C                      I=2,NMOM. The recurrence relation for the monic
C                      polynomials is defined as:
C
C                         P (x)  =  1
C                          0
C
C                         P (x)  =  (x - A )
C                          1              1
C
C                         P   (x)  =  (x - A   ) * P (x) - B   P   (x)
C                          i+1              i+1     i       i+1 i-1
C
C                The result will consist in the diagonal and offdiagonal
C                parts of the symmetric terminal matrix:
C
C                            DIA (I) , I = 1,NGQP
C                            OFF (I) , I = 1,NGQP-1
C
C                Proceed now with Sack and Donovan's algorithm.
C                Handle low number (1 or 2) of quadrature points first.
C
C
C
            IF (NGQP.EQ.1) THEN
                DIA (1) = MOM (1) + A (1)
            ELSE IF (NGQP.EQ.2) THEN
                SIGMA = MOM (1) + A (1)
                DIA (1) = SIGMA
                THETA = (A (2) - SIGMA) * MOM (1) + MOM (2) + B (2)
                OFF (1) = DSQRT (THETA)
                DIA (2) = (((A(3)-SIGMA)*MOM(2)+MOM(3)+B(3)*MOM(1))
     +                    /THETA) - MOM(1)+A(2)
            ELSE
C
C
C             ...Handle case for number of quadrature points > 2.
C                Set maximum values for I and J and evaluate first
C                diagonal element.
C
C
                IMAX = NGQP - 1
                JMAX = NGQP + IMAX
                DO J = 1,JMAX
                   ROW1 (J) = MOM (J)
                END DO
                SIGMA = ROW1 (1) + A (1)
                DIA (1) = SIGMA
C
C
C             ...evaluate 2nd row of terminal matrix.
C
C
                ROW2 (1) = (A (2) - SIGMA) * ROW1 (1) + ROW1 (2) + B (2)
                THETA = ROW2 (1)
                OFF (1) = DSQRT (THETA)
                JMAX = JMAX - 1
                DO J = 2,JMAX
                   JP1 = J + 1
                   ROW2 (J) = (A (JP1) - SIGMA) * ROW1 (J) + ROW1 (JP1)
     +                        + B (JP1) * ROW1 (J-1)
                END DO
                SIGMA = (ROW2 (2) / THETA) - ROW1 (1) + A (2)
                DIA (2) = SIGMA
C
C
C             ...proceed with higher rows.
C
C
                DO 200 I = 2,IMAX
                   IP1 = I + 1
                   JMAX = JMAX - 1

                   IF ( MOD (I,2) .EQ. 0 ) THEN
                        DO J = I,JMAX
                           JP1 = J + 1
                           ROW1 (J) = (A (JP1) - SIGMA) * ROW2 (J)
     +                              + ROW2 (JP1) + B (JP1) * ROW2 (J-1)
     +                              - THETA * ROW1 (J)
                        END DO
                        SIGMA = A (IP1) - (ROW2 (I) / ROW2 (I-1))
     +                                  + (ROW1 (IP1) / ROW1 (I))
                        THETA = ROW1 (I) / ROW2 (I-1)
                   ELSE
                        DO J = I,JMAX
                           JP1 = J + 1
                           ROW2 (J) = (A (JP1) - SIGMA) * ROW1 (J)
     +                              + ROW1 (JP1) + B (JP1) * ROW1 (J-1)
     +                              - THETA * ROW2 (J)
                        END DO
                        SIGMA = A (IP1) - (ROW1 (I) / ROW1 (I-1))
     +                                  + (ROW2 (IP1) / ROW2 (I))
                        THETA = ROW2 (I) / ROW1 (I-1)
                   END IF

                   DIA (IP1) = SIGMA
                   OFF (I) = DSQRT (THETA)

  200           CONTINUE
            END IF
C
C
C             ...The last section computes the gaussian quadrature roots
C                and weights from a previously established symmetric
C                terminal matrix by the Golub-Welsch algorithm (see
C                G.H. Golub and J.H. Welsch, Math. of Computation 23,
C                p. 221-230 and A1-A10, 1969), which is based on
C                a result of Wilf (see H.S. Wilf, Mathematics for the
C                Physical Sciences, New York: Wiley, Problem 9, p. 80).
C
C                Wilf has shown that if Z (K,I) is the K-th element of
C                the I-th normalized eigenvector of the terminal matrix
C                corresponding to the I-th eigenvalue D (I), then the
C                roots RTS (i.e. the zeros of the N-th orthogonal
C                monic polynomial) and weights WTS to be used for the
C                gaussian quadrature are given by:
C
C                           RTS (I) = D (I)
C                           WTS (I) = MOMZERO * (Z (1,I)**2)
C
C                where MOMZERO is the value of the definite integral
C                over the weight function W (x) alone. In our case
C                it is equal to the value of the zeroth Rys moment.
C
C                The present section performs hence a diagonalization
C                of the tridiagonal symmetric terminal matrix keeping
C                only the first component of the eigenvectors and sets
C                the roots and weights equal to the above relations.
C                The diagonalization code was derived from the routine
C                IMTQL2 in the EISPACK collection and uses the implicit
C                QL method.
C
C                Note, that the original diagonals DIA and offdiagonals
C                OFF of the terminal matrix are destroyed during the
C                diagonalization process !!!
C
C                Handle special case, if order of terminal matrix is 1.
C
C
C
            IF (NGQP.EQ.1) THEN
                NRTS = NRTS + 1
                RTS (NRTS) = DIA (1)
                WTS (NRTS) = MOMZERO
            ELSE
C
C
C             ...initialize vector for collecting first component
C                of eigenvectors. To save space, array A is used,
C                which can be done safely, since its dimension NMOM
C                (# of moments) is always >= # of roots NGQP.
C
C
                A (1) = ONE
                DO J = 2,NGQP
                   A (J) = ZERO
                END DO
C
C
C             ...QL iterations.
C
C
                OFF (NGQP) = ZERO

                DO 300 J = 1,NGQP

                   ITER = 0

 3000              DO M = J,NGQP
                      IF (M.EQ.NGQP) GOTO 3300
                      TEST1 = DABS (DIA (M))  +  DABS (DIA (M+1))
                      TEST2 = TEST1 + DABS (OFF (M))
                      IF (TEST2 .EQ. TEST1) GOTO 3300
                   END DO

 3300              P = DIA (J)
                   IF (M.EQ.J) GOTO 300
                   IF (ITER.EQ.MXITER) THEN
                       WRITE (*,*) ' Root number ',J,' not converged! '
                       WRITE (*,*) ' Program bug? More iterations? '
                       WRITE (*,*) ' oed__rys_x_roots_weights '
                       STOP
                   END IF
                   ITER = ITER + 1

                   G = (DIA (J+1) - P) / (TWO * OFF (J))
                   R = DSQRT (G*G + ONE)
                   G = DIA (M) - P + OFF (J) / (G + DSIGN (R,G))
                   S = ONE
                   C = ONE
                   P = ZERO

                   DO I = M-1,J,-1
                      F = S * OFF (I)
                      D = C * OFF (I)
                      R = DSQRT (F*F + G*G)
                      OFF (I+1) = R
                      IF (R .EQ. ZERO) THEN
                          DIA (I+1) = DIA (I+1) - P
                          OFF (M) = ZERO
                          GOTO 3000
                      END IF
                      S = F / R
                      C = G / R
                      G = DIA (I+1) - P
                      R = (DIA (I) - G) * S + TWO * C * D
                      P = S * R
                      DIA (I+1) = G + P
                      G = C * R - D
                      F = A (I+1)
                      A (I+1) = S * A (I) + C * F
                      A (I) = C * A (I) - S * F
                   END DO

                   DIA (J) = DIA (J) - P
                   OFF (J) = G
                   OFF (M) = ZERO
                   GOTO 3000

  300           CONTINUE
C
C
C             ...calculate roots and weights. Since it is known that
C                the roots must lay between 0 and 1, a check is
C                made on them to see if they are actually within this
C                range.
C
C
                DO I = 1,NGQP
                   ROOT = DIA (I)
                   IF (ROOT.LT.ZERO .OR. ROOT.GT.ONE) THEN
                       WRITE (*,*) ' Quadrature root not in range 0-1! '
                       WRITE (*,*) ' I,I-th Root = ',I,ROOT
                       WRITE (*,*) ' oed__rys_x_roots_weights '
                   STOP
                   END IF
                   RTS (NRTS+I) = ROOT
                   WTS (NRTS+I) = MOMZERO * (A (I)**2)
                END DO

                NRTS = NRTS + NGQP

            END IF
C
C
C             ...next T-value
C
C
 1000    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
