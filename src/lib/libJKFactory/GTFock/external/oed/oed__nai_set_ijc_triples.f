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
         SUBROUTINE  OED__NAI_SET_IJC_TRIPLES
     +
     +                    ( NUCLEI,
     +                      XN,YN,ZN,NCHARGE,
     +                      NPGTOA,NPGTOB,NPGTOAB,
     +                      ATOMIC,EQUALAB,
     +                      SWAPRS,
     +                      XA,YA,ZA,XB,YB,ZB,
     +                      RNABSQ,
     +                      PREFACT,
     +                      ALPHAA,ALPHAB,
     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                      SCREEN,
     +
     +                               EMPTY,
     +                               NIJ,NCEN,
     +                               PRIMA,PRIMB,
     +                               NUCCEN,
     +                               ZCORE )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_SET_IJC_TRIPLES
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation determines the ij primitive exponent
C                pairs which will actually be considered for
C                contraction and the nuclear attraction centers that
C                will be processed. It is here where the internal
C                prescreening (if any) of nuclear attraction primitives
C                must be applied.
C
C                There are two ways the ij pairs can be ordered:
C                1) i runs fastest or 2) j runs fastest. The ordering
C                to be chosen depends on the order the contraction will
C                be performed. We have the following rule:
C
C                     If # primitives on A > # primitives on B
C                     (i.e. the SWAPRS = .true. case), the contraction
C                     will be performed first on the A primitives
C                     followed by the B primitives in order to save
C                     intermediate storage space for the 1st half
C                     transformation. Hence each B primitive will have
C                     associated a set of A primitives for the 1st
C                     half contraction and thus the order for the
C                     primitive pairs ij must be such that i runs
C                     fastest. If # primitives on A =< # primitives
C                     on B we have to let j run fastest.
C
C                The prescreening will be based on the fact that
C                for normalized cartesian integrals the primitive
C                (s|s) nuclear attraction integrals will be largest
C                in magnitude for all exponent combinations and
C                nuclear attraction centers. The normalized (s|s)
C                nuclear attraction integral value expression is:
C
C                 Sqrt (32/pi) * (-Z) * (ab)**(3/4) * 1/(a+b)
C                              * rho (a,b) * F0 ((a+b)*R(PC)**2)
C
C                where a and b are the individual gaussian exponents,
C                R(PC) is the distance between points P, the location of
C                the gaussian product function and C, the nuclear
C                attraction center, Z is the nuclear charge at center
C                C, and rho (a,b) is the exponential prefactor given
C                by the expression:
C
C                                      - (ab/(a+b)) * R(AB)**2
C                       rho (a,b) = exp
C
C                where R(AB)**2 is the square of the distance between
C                the atomic centers A and B. The exponential prefactor
C                measures to which extent the two gaussians overlap
C                in space and thus is a direct measure of the magnitude
C                of their product. The Boys function F0 argument is
C                a measure of the coulombic interaction between the
C                electronic cloud and the nuclear attraction center
C                and is a measure on how far apart the electronic cloud
C                and the nuclear attraction center are in space. 
C
C                There are two kinds of prescreening applied here:
C
C                    1) Prescreening of ij exponents
C                    2) Prescreening of nuclear attraction centers
C
C                These prescreenings will be applied according to the
C                following scheme:
C
C                    Find exponents a(min) and b(min), giving the
C                    largest electronic cloud distribution in space.
C
C                    Loop over all nuclear attraction centers C
C
C                         Calculate minimum distance d(C,min) between
C                         line segment AB and C
C
C                         Calculate (s|s) NAI integral using a(min),
C                         b(min) and R(PC) = d(C,min)
C
C                         Skip center C, if |(s|s)| < TOL
C
C                    continue
C
C                    From the reduced set of centers C, choose the
C                    one that has given the largest (s|s) integral
C                    previously. From that chosen center use its
C                    R(PC) value for the Boys function argument
C                    calculations and the corresponding nuclear charge
C                    Z.
C
C                    Loop over all ij exponent pairs
C
C                         Calculate (s|s) NAI integral using i-th and
C                         j-th exponent and use R(PC) = d(C(min),min)
C
C                         Skip ij pair, if |(s|s)| < TOL
C
C                    continue
C
C                This prescreening scheme will require as many
C                evaluations of the Boys function as there are nuclear
C                attraction centers plus ij exponent pairs.
C
C                Since evaluation of the exponential prefactors is
C                rather costly, the ones surviving the prescreening
C                are transmitted back to the calling routine for
C                further use.
C
C
C                  Input:
C
C                    NUCLEI       =  the original # of nuclear
C                                    attraction centers
C                    XN,YN,ZN     =  the original x,y,z-coordinates
C                                    of the nuclear attraction centers
C                    NCHARGE      =  the original charges of the
C                                    nuclear attraction centers
C                    NPGTOx       =  # of primitives per contraction
C                                    for contraction shells x = A,B
C                    NPGTOAB      =  # of primitive pairs for
C                                    contraction shell pair AB
C                    ATOMIC       =  indicates, if atomic integrals
C                                    are evaluated
C                    EQUALAB      =  indicates, if csh A and csh B are
C                                    considered to be equal
C                    SWAPRS       =  is .true. if the contraction order
C                                    of the primitive pairs AB will
C                                    be performed in reverse order
C                                    BA for efficiency reasons
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = A and B
C                    RNABSQ       =  the square distance R(AB)**2
C                                    between centers A and B
C                    PREFACT      =  the value of sqrt (8)
C                    ALPHAx       =  primitive exponents for centers
C                                    x=A,B
C                    FTABLE       =  Fm (T) table for interpolation
C                                    in low T region
C                    MGRID        =  maximum m in Fm (T) table
C                    NGRID        =  # of T's for which Fm (T) table
C                                    was set up
C                    TMAX         =  maximum T in Fm (T) table
C                    TSTEP        =  difference between two consecutive
C                                    T's in Fm (T) table
C                    TVSTEP       =  Inverse of TSTEP
C                    SCREEN       =  is true, if screening will be
C                                    performed
C
C                  Output:
C
C                    EMPTY        =  is true, if no ij pairs were
C                                    established
C                    NIJ          =  # of ij pairs after screening
C                    NCEN         =  # of ij nuclear attraction
C                                    centers after screening
C                    PRIMx        =  i,j labels of primitives for x=A,B
C                    NUCCEN       =  will hold the index labels of
C                                    the nuclear attraction centers
C                                    surviving the screening process
C                    ZCORE        =  flp array holding NIJ exponential
C                                    prefactors rho (a,b)
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT  NONE

         LOGICAL   ATOMIC
         LOGICAL   EMPTY
         LOGICAL   EQUALAB
         LOGICAL   SCREEN
         LOGICAL   SWAPRS

         INTEGER   I,J,N
         INTEGER   MGRID,NGRID
         INTEGER   NIJ,NCEN
         INTEGER   NPGTOA,NPGTOB,NPGTOAB
         INTEGER   NUCLEI
         INTEGER   TGRID

         INTEGER   NUCCEN (1:NUCLEI)
         INTEGER   PRIMA  (1:NPGTOAB)
         INTEGER   PRIMB  (1:NPGTOAB)

         DOUBLE PRECISION  A,B
         DOUBLE PRECISION  AB,ABMIN
         DOUBLE PRECISION  DELTA
         DOUBLE PRECISION  F0
         DOUBLE PRECISION  OED__DSQMIN_LINE_SEGMENTS
         DOUBLE PRECISION  P,PINV,PMIN
         DOUBLE PRECISION  PREFACT
         DOUBLE PRECISION  RHOAB
         DOUBLE PRECISION  RSQ,RNABSQ,RMINSQ
         DOUBLE PRECISION  SMAXAB
         DOUBLE PRECISION  SSINT,SSMAX
         DOUBLE PRECISION  T,TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB,XC,YC,ZC
         DOUBLE PRECISION  ZATOM,ZFORIJ
         DOUBLE PRECISION  ZERO,SIXTH,FIFTH,FOURTH,THIRD,HALF,
     +                     ZP75,ONE,PI,TOL

         DOUBLE PRECISION  ALPHAA  (1:NPGTOA)
         DOUBLE PRECISION  ALPHAB  (1:NPGTOB)

         DOUBLE PRECISION  XN      (1:NUCLEI)
         DOUBLE PRECISION  YN      (1:NUCLEI)
         DOUBLE PRECISION  ZN      (1:NUCLEI)
         DOUBLE PRECISION  NCHARGE (1:NUCLEI)

         DOUBLE PRECISION  ZCORE (*)

         DOUBLE PRECISION  FTABLE (0:MGRID,0:NGRID)

         PARAMETER  (ZERO    =  0.D0)
         PARAMETER  (SIXTH   =  0.166666666666667D0)
         PARAMETER  (FIFTH   =  0.2D0)
         PARAMETER  (FOURTH  =  0.25D0)
         PARAMETER  (THIRD   =  0.333333333333333D0)
         PARAMETER  (HALF    =  0.5D0)
         PARAMETER  (ZP75    =  0.75D0)
         PARAMETER  (ONE     =  1.D0)
         PARAMETER  (PI      =  3.141592653589793D0)
         PARAMETER  (TOL     =  1.D-14)
C
C
C------------------------------------------------------------------------
C
C
C             ...the no screening section.
C
C
         EMPTY = .FALSE.

         IF (.NOT.SCREEN) THEN

             NIJ = 0
             IF (EQUALAB) THEN
                 DO I = 1,NPGTOA
                 DO J = 1,I
                    NIJ = NIJ + 1
                    ZCORE (NIJ) = ONE
                    PRIMA (NIJ) = I
                    PRIMB (NIJ) = J
                 END DO
                 END DO
             ELSE
                 IF (SWAPRS) THEN
                     IF (ATOMIC) THEN
                         DO J = 1,NPGTOB
                         DO I = 1,NPGTOA
                            NIJ = NIJ + 1
                            ZCORE (NIJ) = ONE
                            PRIMA (NIJ) = I
                            PRIMB (NIJ) = J
                         END DO
                         END DO
                     ELSE
                         DO J = 1,NPGTOB
                            B = ALPHAB (J)
                            DO I = 1,NPGTOA
                               A = ALPHAA (I)
                               NIJ = NIJ + 1
                               ZCORE (NIJ) = DEXP (- A*B*RNABSQ / (A+B))
                               PRIMA (NIJ) = I
                               PRIMB (NIJ) = J
                            END DO
                         END DO
                     END IF
                 ELSE
                     IF (ATOMIC) THEN
                         DO I = 1,NPGTOA
                         DO J = 1,NPGTOB
                            NIJ = NIJ + 1
                            ZCORE (NIJ) = ONE
                            PRIMA (NIJ) = I
                            PRIMB (NIJ) = J
                         END DO
                         END DO
                     ELSE
                         DO I = 1,NPGTOA
                            A = ALPHAA (I)
                            DO J = 1,NPGTOB
                               B = ALPHAB (J)
                               NIJ = NIJ + 1
                               ZCORE (NIJ) = DEXP (- A*B*RNABSQ / (A+B))
                               PRIMA (NIJ) = I
                               PRIMB (NIJ) = J
                            END DO
                         END DO
                     END IF
                 END IF
             END IF

             NCEN = NUCLEI

             DO N = 1,NUCLEI
                NUCCEN (N) = N
             END DO

             RETURN

         END IF
C
C
C             ...enter screening. Determine overall maxima and minima.
C
C
         A = ALPHAA (1)
         DO I = 2,NPGTOA
            A = DMIN1 (A,ALPHAA (I))
         END DO

         B = ALPHAB (1)
         DO I = 2,NPGTOB
            B = DMIN1 (B,ALPHAB (I))
         END DO

         PMIN = A + B
         ABMIN = A * B
         PINV = ONE / PMIN

         SMAXAB = PREFACT * (ABMIN ** ZP75) * PINV
     +                    * DEXP (-ABMIN * RNABSQ * PINV)

         SSMAX = ZERO
C
C
C             ...perform screening on nuclear attraction centers.
C
C
         NCEN = 0

         DO N = 1,NUCLEI

            XC = XN (N)
            YC = YN (N)
            ZC = ZN (N)
            ZATOM = NCHARGE (N)

            RSQ = OED__DSQMIN_LINE_SEGMENTS
     +
     +                 ( XA,YA,ZA,
     +                   XB,YB,ZB,
     +                   XC,YC,ZC,
     +                   XC,YC,ZC )
     +
     +
            T = PMIN * RSQ

            IF (T.EQ.ZERO) THEN
                F0 = ONE
            ELSE IF (T.LE.TMAX) THEN
                TGRID = INT (T * TVSTEP + HALF)
                DELTA = TGRID * TSTEP - T

                F0 = ((((( FTABLE (6,TGRID)  * DELTA * SIXTH
     +                   + FTABLE (5,TGRID)) * DELTA * FIFTH
     +                   + FTABLE (4,TGRID)) * DELTA * FOURTH
     +                   + FTABLE (3,TGRID)) * DELTA * THIRD
     +                   + FTABLE (2,TGRID)) * DELTA * HALF
     +                   + FTABLE (1,TGRID)) * DELTA
     +                   + FTABLE (0,TGRID)
            ELSE
                F0 = HALF * DSQRT (PI/T)
            END IF

            SSINT = SMAXAB * F0 * ZATOM

            IF (DABS (SSINT).GE.TOL) THEN
                NCEN = NCEN + 1
                NUCCEN (NCEN) = N
                IF (DABS (SSINT).GT.SSMAX) THEN
                    RMINSQ = RSQ
                    ZFORIJ = ZATOM
                    SSMAX = SSINT
                END IF
            ELSE
C                WRITE (*,*) ' Skip nai center C = ',N
            END IF

         END DO

         IF (NCEN.EQ.0) THEN
             EMPTY = .TRUE.
             RETURN
         END IF
C
C
C             ...perform screening on ij exponent pairs.
C
C
         NIJ = 0

         IF (EQUALAB) THEN

             DO I = 1,NPGTOA
                A = ALPHAA (I)
                DO J = 1,I
                   B = ALPHAB (J)
                   P = A + B
                   AB = A * B
                   PINV = ONE / P
                   T = P * RMINSQ

                   IF (T.EQ.ZERO) THEN
                       F0 = ONE
                   ELSE IF (T.LE.TMAX) THEN
                       TGRID = INT (T * TVSTEP + HALF)
                       DELTA = TGRID * TSTEP - T

                       F0 = ((((( FTABLE (6,TGRID)  * DELTA * SIXTH
     +                          + FTABLE (5,TGRID)) * DELTA * FIFTH
     +                          + FTABLE (4,TGRID)) * DELTA * FOURTH
     +                          + FTABLE (3,TGRID)) * DELTA * THIRD
     +                          + FTABLE (2,TGRID)) * DELTA * HALF
     +                          + FTABLE (1,TGRID)) * DELTA
     +                          + FTABLE (0,TGRID)
                   ELSE
                       F0 = HALF * DSQRT (PI/T)
                   END IF

                   SSINT = PREFACT * (AB**ZP75) * PINV * F0 * ZFORIJ

                   IF (DABS (SSINT).GE.TOL) THEN
                       NIJ = NIJ + 1
                       ZCORE (NIJ) = ONE
                       PRIMA (NIJ) = I
                       PRIMB (NIJ) = J
                   ELSE
C                       WRITE (*,*) ' Skip atom nai I J = ',I,J
                   END IF

                END DO
             END DO

         ELSE

             IF (SWAPRS) THEN
                 IF (ATOMIC) THEN
                     DO J = 1,NPGTOB
                        B = ALPHAB (J)
                        DO I = 1,NPGTOA
                           A = ALPHAA (I)
                           P = A + B
                           AB = A * B
                           PINV = ONE / P
                           T = P * RMINSQ

                           IF (T.EQ.ZERO) THEN
                               F0 = ONE
                           ELSE IF (T.LE.TMAX) THEN
                               TGRID = INT (T * TVSTEP + HALF)
                               DELTA = TGRID * TSTEP - T

                               F0 = (((((FTABLE (6,TGRID) *DELTA*SIXTH
     +                                 + FTABLE (5,TGRID))*DELTA*FIFTH
     +                                 + FTABLE (4,TGRID))*DELTA*FOURTH
     +                                 + FTABLE (3,TGRID))*DELTA*THIRD
     +                                 + FTABLE (2,TGRID))*DELTA*HALF
     +                                 + FTABLE (1,TGRID))*DELTA
     +                                 + FTABLE (0,TGRID)
                           ELSE
                               F0 = HALF * DSQRT (PI/T)
                           END IF

                           SSINT = PREFACT * (AB**ZP75) * PINV
     +                                                  * F0 * ZFORIJ
                           IF (DABS (SSINT).GE.TOL) THEN
                               NIJ = NIJ + 1
                               ZCORE (NIJ) = ONE
                               PRIMA (NIJ) = I
                               PRIMB (NIJ) = J
                           ELSE
C                               WRITE (*,*) ' Skip nai I J = ',I,J
                           END IF

                        END DO
                     END DO
                 ELSE
                     DO J = 1,NPGTOB
                        B = ALPHAB (J)
                        DO I = 1,NPGTOA
                           A = ALPHAA (I)
                           P = A + B
                           AB = A * B
                           PINV = ONE / P
                           RHOAB = DEXP (- AB * RNABSQ * PINV)
                           T = P * RMINSQ

                           IF (T.EQ.ZERO) THEN
                               F0 = ONE
                           ELSE IF (T.LE.TMAX) THEN
                               TGRID = INT (T * TVSTEP + HALF)
                               DELTA = TGRID * TSTEP - T

                               F0 = (((((FTABLE (6,TGRID) *DELTA*SIXTH
     +                                 + FTABLE (5,TGRID))*DELTA*FIFTH
     +                                 + FTABLE (4,TGRID))*DELTA*FOURTH
     +                                 + FTABLE (3,TGRID))*DELTA*THIRD
     +                                 + FTABLE (2,TGRID))*DELTA*HALF
     +                                 + FTABLE (1,TGRID))*DELTA
     +                                 + FTABLE (0,TGRID)
                           ELSE
                               F0 = HALF * DSQRT (PI/T)
                           END IF

                           SSINT = PREFACT * (AB**ZP75) * PINV * RHOAB
     +                                                  * F0 * ZFORIJ
                           IF (DABS (SSINT).GE.TOL) THEN
                               NIJ = NIJ + 1
                               ZCORE (NIJ) = RHOAB
                               PRIMA (NIJ) = I
                               PRIMB (NIJ) = J
                           ELSE
C                               WRITE (*,*) ' Skip nai I J = ',I,J
                           END IF

                        END DO
                     END DO
                 END IF
             ELSE
                 IF (ATOMIC) THEN
                     DO I = 1,NPGTOA
                        A = ALPHAA (I)
                        DO J = 1,NPGTOB
                           B = ALPHAB (J)
                           P = A + B
                           AB = A * B
                           PINV = ONE / P
                           T = P * RMINSQ

                           IF (T.EQ.ZERO) THEN
                               F0 = ONE
                           ELSE IF (T.LE.TMAX) THEN
                               TGRID = INT (T * TVSTEP + HALF)
                               DELTA = TGRID * TSTEP - T

                               F0 = (((((FTABLE (6,TGRID) *DELTA*SIXTH
     +                                 + FTABLE (5,TGRID))*DELTA*FIFTH
     +                                 + FTABLE (4,TGRID))*DELTA*FOURTH
     +                                 + FTABLE (3,TGRID))*DELTA*THIRD
     +                                 + FTABLE (2,TGRID))*DELTA*HALF
     +                                 + FTABLE (1,TGRID))*DELTA
     +                                 + FTABLE (0,TGRID)
                           ELSE
                               F0 = HALF * DSQRT (PI/T)
                           END IF

                           SSINT = PREFACT * (AB**ZP75) * PINV
     +                                                  * F0 * ZFORIJ
                           IF (DABS (SSINT).GE.TOL) THEN
                               NIJ = NIJ + 1
                               ZCORE (NIJ) = ONE
                               PRIMA (NIJ) = I
                               PRIMB (NIJ) = J
                           ELSE
C                               WRITE (*,*) ' Skip nai I J = ',I,J
                           END IF

                        END DO
                     END DO
                 ELSE
                     DO I = 1,NPGTOA
                        A = ALPHAA (I)
                        DO J = 1,NPGTOB
                           B = ALPHAB (J)
                           P = A + B
                           AB = A * B
                           PINV = ONE / P
                           RHOAB = DEXP (- AB * RNABSQ * PINV)
                           T = P * RMINSQ

                           IF (T.EQ.ZERO) THEN
                               F0 = ONE
                           ELSE IF (T.LE.TMAX) THEN
                               TGRID = INT (T * TVSTEP + HALF)
                               DELTA = TGRID * TSTEP - T

                               F0 = (((((FTABLE (6,TGRID) *DELTA*SIXTH
     +                                 + FTABLE (5,TGRID))*DELTA*FIFTH
     +                                 + FTABLE (4,TGRID))*DELTA*FOURTH
     +                                 + FTABLE (3,TGRID))*DELTA*THIRD
     +                                 + FTABLE (2,TGRID))*DELTA*HALF
     +                                 + FTABLE (1,TGRID))*DELTA
     +                                 + FTABLE (0,TGRID)
                           ELSE
                               F0 = HALF * DSQRT (PI/T)
                           END IF

                           SSINT = PREFACT * (AB**ZP75) * PINV * RHOAB
     +                                                  * F0 * ZFORIJ
                           IF (DABS (SSINT).GE.TOL) THEN
                               NIJ = NIJ + 1
                               ZCORE (NIJ) = RHOAB
                               PRIMA (NIJ) = I
                               PRIMB (NIJ) = J
                           ELSE
C                               WRITE (*,*) ' Skip nai I J = ',I,J
                           END IF

                        END DO
                     END DO
                 END IF
             END IF

         END IF

         IF (NIJ.EQ.0) THEN
             EMPTY = .TRUE.
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
