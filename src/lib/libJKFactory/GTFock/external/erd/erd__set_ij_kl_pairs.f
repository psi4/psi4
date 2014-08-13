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
         SUBROUTINE  ERD__SET_IJ_KL_PAIRS
     +
     +                    ( NPGTOA,NPGTOB,NPGTOC,NPGTOD,
     +                      NPGTOAB,NPGTOCD,
     +                      ATOMAB,ATOMCD,
     +                      EQUALAB,EQUALCD,
     +                      SWAPRS,SWAPTU,
     +                      XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     +                      RNABSQ,RNCDSQ,
     +                      PREFACT,
     +                      ALPHAA,ALPHAB,ALPHAC,ALPHAD,
     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                      SCREEN,
     +
     +                               EMPTY,
     +                               NIJ,NKL,
     +                               PRIMA,PRIMB,PRIMC,PRIMD,
     +                               RHO )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__SET_IJ_KL_PAIRS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation determines the ij and kl primitive
C                exponent pairs which will actually be considered
C                for contraction. It is here where the internal
C                prescreening (if any) of primitives must be applied.
C
C                There are two ways the individual ij and kl pairs
C                can be ordered: 1) i,k run fastest or 2) j,l run
C                fastest. The ordering to be chosen depends on the
C                order the contraction will be performed. We have
C                the following rule, explained here for the ij pairs,
C                with the same being applied for the kl pairs:
C
C                     If # primitives on A > # primitives on B
C                     (i.e. the SWAPRS = .true. case), the contraction
C                     will be performed first on the A primitives
C                     followed by the B primitives in order to save
C                     intermediate storage space for the 1st quarter
C                     transformation. Hence each B primitive will have
C                     associated a set of A primitives for the 1st
C                     quarter contraction and thus the order for
C                     the primitive pairs ij must be such that i
C                     runs fastest. If # primitives on A =< # primitives
C                     on B we have to let j run fastest.
C
C                The prescreening will be based on the fact that
C                for normalized cartesian integrals the primitive
C                (ss|ss) integrals will be largest in magnitude for
C                all exponent pair combinations. The normalized
C                (ss|ss) integral expression is:
C
C                16/Sqrt(Pi) * (abcd)**(3/4) * rho(a,b) * rho(c,d)
C                            * F0 ((pq/(p+q))*R(PQ)**2) / (pq*Sqrt(p+q))
C
C                where a,b,c,d are the individual gaussian exponents,
C                p = a+b and q = c+d their partial sums, R(PQ) is the
C                length between points P and Q where the gaussian
C                product function is located and rho(a,b) and rho(c,d)
C                are exponential prefactors given by the expressions:
C
C                                     - (ab/(a+b)) * R(AB)**2
C                       rho(a,b) = exp
C
C                                     - (cd/(c+d)) * R(CD)**2
C                       rho(c,d) = exp
C
C                where a and b are the i-th and j-th exponent and
C                R(AB)**2 is the square of the distance between the
C                i-th and j-th atomic centers A and B. The exponential
C                prefactor measures to which extent the two i-th and
C                j-th primitive gaussians overlap in space and thus
C                is a direct measure of the magnitude of their product.
C                The Boys function F0 argument on the other hand is
C                a measure of the coulombic interaction between the
C                two electronic clouds and is thus a measure on how
C                far apart the ij and kl gaussian pairs are in space. 
C
C                The prescreening will be done in a 'K2' fashion,
C                fixing first the maximum scaling factor for the CD
C                part and testing all AB pairs and then vice versa.
C                The maximum scaling factors are given by:
C
C
C                   For the AB part:
C
C                         SMAXAB = (ab)**(3/4) * rho(a,b) / (a+b)
C                                  using smallest a and b exponents
C
C                   For the CD part:
C
C                         SMAXCD = (cd)**(3/4) * rho(c,d) / (c+d)
C                                  using smallest c and d exponents
C
C
C                and are precalculated before entering both K2 loops
C                including a multiplication by the prefactor 16/Sqrt(Pi)
C                for convenience purposes.
C
C                As far as the F0 arguments go, we fix p and q to be
C                the respective minima and for R(PQ) we use the minimum
C                possible distance between the two line segments
C                joining atomic centers A with B and C with D.
C
C                An ij or kl pair will be rejected, if the corresponding
C                upper bound of the (ss|ss) integral is less than a
C                certain tolerance value TOL.
C
C                Since evaluation of the exponential prefactors is
C                rather costly, the ones surviving the prescreening
C                are transmitted back to the calling routine for
C                further use.
C
C
C                  Input:
C
C                    NPGTOx       =  # of primitives per contraction
C                                    for contraction shells x = A,B,C,D
C                    NPGTOxy      =  # of primitive pairs for
C                                    contraction shell pairs xy = AB
C                                    and CD
C                    ATOMxy       =  indicates, if atoms x and y are
C                                    considered to be equal for the
C                                    pairs xy = AB and CD
C                    EQUALxy      =  indicates, if csh x and csh y are
C                                    considered to be equal for the
C                                    pairs xy = AB and CD
C                    SWAPRS(TU)   =  is .true. if the contraction order
C                                    of the primitive pairs AB(CD) will
C                                    be performed in reverse order
C                                    BA(DC) for efficiency reasons
C                    RNxySQ       =  the square distances R(xy)**2
C                                    between centers x and y. We have
C                                    xy = AB and CD
C                    PREFACT      =  the value of 16/sqrt(pi)
C                    ALPHAx       =  primitive exponents for centers
C                                    x=A,B,C,D
C                    FTABLE       =  Fm (T) table for interpolation
C                                    for low T's
C                    MGRID        =  maximum m in Fm (T) table
C                    NGRID        =  # of T's for which Fm (T) table
C                                    was set up
C                    TMAX         =  maximum T in Fm (T) table
C                    TSTEP        =  difference between two consecutive
C                                    T's in Fm (T) table
C                    TVSTEP       =  inverse of TSTEP
C                    SCREEN       =  is true, if screening will be
C                                    performed
C
C                  Output:
C
C                    EMPTY        =  is true, if no ij and/or kl pairs
C                                    were established
C                    NIJ,NKL      =  # of ij and kl pairs
C                    PRIMx        =  i,j,k,l labels of primitives
C                                    for x=A,B,C,D
C                    RHO          =  NIJ exponential prefactors rho(a,b)
C                                    + NKL exponential prefactors
C                                    rho(c,d), in that order
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

         LOGICAL     ATOMAB,ATOMCD
         LOGICAL     EMPTY
         LOGICAL     EQUALAB,EQUALCD
         LOGICAL     SCREEN
         LOGICAL     SWAPRS,SWAPTU

         INTEGER     I,J,K,L
         INTEGER     MGRID,NGRID,TGRID
         INTEGER     NIJ,NKL
         INTEGER     NPGTOA,NPGTOB,NPGTOC,NPGTOD
         INTEGER     NPGTOAB,NPGTOCD

         INTEGER     PRIMA (1:NPGTOAB)
         INTEGER     PRIMB (1:NPGTOAB)
         INTEGER     PRIMC (1:NPGTOCD)
         INTEGER     PRIMD (1:NPGTOCD)

         DOUBLE PRECISION  A,B,C,D
         DOUBLE PRECISION  AB,CD
         DOUBLE PRECISION  ABMIN,CDMIN
         DOUBLE PRECISION  DELTA
         DOUBLE PRECISION  ERD__DSQMIN_LINE_SEGMENTS
         DOUBLE PRECISION  F0
         DOUBLE PRECISION  P,Q,PMIN,QMIN,PINV,QINV,PQPINV
         DOUBLE PRECISION  PREFACT
         DOUBLE PRECISION  RHOAB,RHOCD
         DOUBLE PRECISION  RNABSQ,RNCDSQ,RMINSQ
         DOUBLE PRECISION  SMAXAB,SMAXCD
         DOUBLE PRECISION  SSSSMX
         DOUBLE PRECISION  T,TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION  XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
         DOUBLE PRECISION  ZERO,SIXTH,FIFTH,FOURTH,THIRD,HALF,
     +                     ZP75,ONE,PI,TOL

         DOUBLE PRECISION  ALPHAA  (1:NPGTOA)
         DOUBLE PRECISION  ALPHAB  (1:NPGTOB)
         DOUBLE PRECISION  ALPHAC  (1:NPGTOC)
         DOUBLE PRECISION  ALPHAD  (1:NPGTOD)

         DOUBLE PRECISION  RHO (1:NPGTOAB+NPGTOCD)

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
                    RHO (NIJ) = ONE
                    PRIMA (NIJ) = I
                    PRIMB (NIJ) = J
                 END DO
                 END DO
             ELSE
                 IF (SWAPRS) THEN
                     IF (ATOMAB) THEN
                         DO J = 1,NPGTOB
                         DO I = 1,NPGTOA
                            NIJ = NIJ + 1
                            RHO (NIJ) = ONE
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
                               RHO (NIJ) = DEXP (- A*B*RNABSQ / (A+B))
                               PRIMA (NIJ) = I
                               PRIMB (NIJ) = J
                            END DO
                         END DO
                     END IF
                 ELSE
                     IF (ATOMAB) THEN
                         DO I = 1,NPGTOA
                         DO J = 1,NPGTOB
                            NIJ = NIJ + 1
                            RHO (NIJ) = ONE
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
                               RHO (NIJ) = DEXP (- A*B*RNABSQ / (A+B))
                               PRIMA (NIJ) = I
                               PRIMB (NIJ) = J
                            END DO
                         END DO
                     END IF
                 END IF
             END IF

             NKL = 0
             IF (EQUALCD) THEN
                 DO K = 1,NPGTOC
                 DO L = 1,K
                    NKL = NKL + 1
                    RHO (NIJ+NKL) = ONE
                    PRIMC (NKL) = K
                    PRIMD (NKL) = L
                 END DO
                 END DO
             ELSE
                 IF (SWAPTU) THEN
                     IF (ATOMCD) THEN
                         DO L = 1,NPGTOD
                         DO K = 1,NPGTOC
                            NKL = NKL + 1
                            RHO (NIJ+NKL) = ONE
                            PRIMC (NKL) = K
                            PRIMD (NKL) = L
                         END DO
                         END DO
                     ELSE
                         DO L = 1,NPGTOD
                            D = ALPHAD (L)
                            DO K = 1,NPGTOC
                               C = ALPHAC (K)
                               NKL = NKL + 1
                               RHO (NIJ+NKL) = DEXP (-C*D*RNCDSQ/(C+D))
                               PRIMC (NKL) = K
                               PRIMD (NKL) = L
                            END DO
                         END DO
                     END IF
                 ELSE
                     IF (ATOMCD) THEN
                         DO K = 1,NPGTOC
                         DO L = 1,NPGTOD
                            NKL = NKL + 1
                            RHO (NIJ+NKL) = ONE
                            PRIMC (NKL) = K
                            PRIMD (NKL) = L
                         END DO
                         END DO
                     ELSE
                         DO K = 1,NPGTOC
                            C = ALPHAC (K)
                            DO L = 1,NPGTOD
                               D = ALPHAD (L)
                               NKL = NKL + 1
                               RHO (NIJ+NKL) = DEXP (-C*D*RNCDSQ/(C+D))
                               PRIMC (NKL) = K
                               PRIMD (NKL) = L
                            END DO
                         END DO
                     END IF
                 END IF
             END IF

             RETURN

         END IF
C
C
C             ...the screening section: determine overall maxima
C                and minima.
C
C
         RMINSQ = ERD__DSQMIN_LINE_SEGMENTS
     +
     +                 ( XA,YA,ZA,
     +                   XB,YB,ZB,
     +                   XC,YC,ZC,
     +                   XD,YD,ZD )
     +
     +
         A = ALPHAA (1)
         DO I = 2,NPGTOA
            A = DMIN1 (A,ALPHAA (I))
         END DO

         B = ALPHAB (1)
         DO I = 2,NPGTOB
            B = DMIN1 (B,ALPHAB (I))
         END DO

         C = ALPHAC (1)
         DO I = 2,NPGTOC
            C = DMIN1 (C,ALPHAC (I))
         END DO

         D = ALPHAD (1)
         DO I = 2,NPGTOD
            D = DMIN1 (D,ALPHAD (I))
         END DO

         PMIN = A + B
         QMIN = C + D
         ABMIN = A * B
         CDMIN = C * D
         PINV = ONE / PMIN
         QINV = ONE / QMIN

         SMAXAB = PREFACT * (ABMIN ** ZP75)
     +                    * DEXP (-ABMIN * RNABSQ * PINV) * PINV
         SMAXCD = PREFACT * (CDMIN ** ZP75)
     +                    * DEXP (-CDMIN * RNCDSQ * QINV) * QINV
C
C
C             ...perform K2 primitive screening on A,B part.
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
                   PQPINV = ONE / (P + QMIN)
                   T = RMINSQ * P * QMIN * PQPINV

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

                   SSSSMX = (AB ** ZP75) * SMAXCD * F0
     +                                   * PINV * DSQRT (PQPINV)

                   IF (SSSSMX.GE.TOL) THEN
                       NIJ = NIJ + 1
                       RHO (NIJ) = ONE
                       PRIMA (NIJ) = I
                       PRIMB (NIJ) = J
                   ELSE
C                       WRITE (*,*) ' Skip atom I J = ',I,J
                   END IF

                END DO
             END DO

         ELSE

             IF (SWAPRS) THEN
                 IF (ATOMAB) THEN
                     DO J = 1,NPGTOB
                        B = ALPHAB (J)
                        DO I = 1,NPGTOA
                           A = ALPHAA (I)
                           P = A + B
                           AB = A * B
                           PINV = ONE / P
                           PQPINV = ONE / (P + QMIN)
                           T = RMINSQ * P * QMIN * PQPINV

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

                           SSSSMX = (AB ** ZP75) * SMAXCD * F0
     +                                    * PINV * DSQRT (PQPINV)

                           IF (SSSSMX.GE.TOL) THEN
                               NIJ = NIJ + 1
                               RHO (NIJ) = ONE
                               PRIMA (NIJ) = I
                               PRIMB (NIJ) = J
                           ELSE
C                               WRITE (*,*) ' Skip I J = ',I,J
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
                           PQPINV = ONE / (P + QMIN)
                           RHOAB = DEXP (- AB * RNABSQ * PINV)
                           T = RMINSQ * P * QMIN * PQPINV

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

                           SSSSMX = (AB ** ZP75) * RHOAB * SMAXCD * F0
     +                                    * PINV * DSQRT (PQPINV)

                           IF (SSSSMX.GE.TOL) THEN
                               NIJ = NIJ + 1
                               RHO (NIJ) = RHOAB
                               PRIMA (NIJ) = I
                               PRIMB (NIJ) = J
                           ELSE
C                               WRITE (*,*) ' Skip I J = ',I,J
                           END IF

                        END DO
                     END DO
                 END IF
             ELSE
                 IF (ATOMAB) THEN
                     DO I = 1,NPGTOA
                        A = ALPHAA (I)
                        DO J = 1,NPGTOB
                           B = ALPHAB (J)
                           P = A + B
                           AB = A * B
                           PINV = ONE / P
                           PQPINV = ONE / (P + QMIN)
                           T = RMINSQ * P * QMIN * PQPINV

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

                           SSSSMX = (AB ** ZP75) * SMAXCD * F0
     +                                    * PINV * DSQRT (PQPINV)

                           IF (SSSSMX.GE.TOL) THEN
                               NIJ = NIJ + 1
                               RHO (NIJ) = ONE
                               PRIMA (NIJ) = I
                               PRIMB (NIJ) = J
                           ELSE
C                               WRITE (*,*) ' Skip I J = ',I,J
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
                           PQPINV = ONE / (P + QMIN)
                           RHOAB = DEXP (- AB * RNABSQ * PINV)
                           T = RMINSQ * P * QMIN * PQPINV

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

                           SSSSMX = (AB ** ZP75) * RHOAB * SMAXCD * F0
     +                                    * PINV * DSQRT (PQPINV)

                           IF (SSSSMX.GE.TOL) THEN
                               NIJ = NIJ + 1
                               RHO (NIJ) = RHOAB
                               PRIMA (NIJ) = I
                               PRIMB (NIJ) = J
                           ELSE
C                               WRITE (*,*) ' Skip I J = ',I,J
                           END IF

                        END DO
                     END DO
                 END IF
             END IF

         END IF

         IF (NIJ.EQ.0) THEN
             EMPTY = .TRUE.
             RETURN
         END IF
C
C
C             ...perform K2 primitive screening on C,D part.
C
C
         NKL = 0

         IF (EQUALCD) THEN

             DO K = 1,NPGTOC
                C = ALPHAC (K)
                DO L = 1,K
                   D = ALPHAD (L)
                   Q = C + D
                   CD = C * D
                   QINV = ONE / Q
                   PQPINV = ONE / (PMIN + Q)
                   T = RMINSQ * PMIN * Q * PQPINV

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

                   SSSSMX = (CD ** ZP75) * SMAXAB * F0
     +                                   * QINV * DSQRT (PQPINV)

                   IF (SSSSMX.GE.TOL) THEN
                       NKL = NKL + 1
                       RHO (NIJ+NKL) = ONE
                       PRIMC (NKL) = K
                       PRIMD (NKL) = L
                   ELSE
C                       WRITE (*,*) ' Skip atom K L = ',K,L
                   END IF

                END DO
             END DO

         ELSE

             IF (SWAPTU) THEN
                 IF (ATOMCD) THEN
                     DO L = 1,NPGTOD
                        D = ALPHAD (L)
                        DO K = 1,NPGTOC
                           C = ALPHAC (K)
                           Q = C + D
                           CD = C * D
                           QINV = ONE / Q
                           PQPINV = ONE / (PMIN + Q)
                           T = RMINSQ * PMIN * Q * PQPINV

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

                           SSSSMX = (CD ** ZP75) * SMAXAB * F0
     +                                    * QINV * DSQRT (PQPINV)

                           IF (SSSSMX.GE.TOL) THEN
                               NKL = NKL + 1
                               RHO (NIJ+NKL) = ONE
                               PRIMC (NKL) = K
                               PRIMD (NKL) = L
                           ELSE
C                               WRITE (*,*) ' Skip K L = ',K,L
                           END IF

                        END DO
                     END DO
                 ELSE
                     DO L = 1,NPGTOD
                        D = ALPHAD (L)
                        DO K = 1,NPGTOC
                           C = ALPHAC (K)
                           Q = C + D
                           CD = C * D
                           QINV = ONE / Q
                           PQPINV = ONE / (PMIN + Q)
                           RHOCD = DEXP (- CD * RNCDSQ * QINV)
                           T = RMINSQ * PMIN * Q * PQPINV

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

                           SSSSMX = (CD ** ZP75) * RHOCD * SMAXAB * F0
     +                                    * QINV * DSQRT (PQPINV)

                           IF (SSSSMX.GE.TOL) THEN
                               NKL = NKL + 1
                               RHO (NIJ+NKL) = RHOCD
                               PRIMC (NKL) = K
                               PRIMD (NKL) = L
                           ELSE
C                               WRITE (*,*) ' Skip K L = ',K,L
                           END IF

                        END DO
                     END DO
                 END IF
             ELSE
                 IF (ATOMCD) THEN
                     DO K = 1,NPGTOC
                        C = ALPHAC (K)
                        DO L = 1,NPGTOD
                           D = ALPHAD (L)
                           Q = C + D
                           CD = C * D
                           QINV = ONE / Q
                           PQPINV = ONE / (PMIN + Q)
                           T = RMINSQ * PMIN * Q * PQPINV

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

                           SSSSMX = (CD ** ZP75) * SMAXAB * F0
     +                                    * QINV * DSQRT (PQPINV)

                           IF (SSSSMX.GE.TOL) THEN
                               NKL = NKL + 1
                               RHO (NIJ+NKL) = ONE
                               PRIMC (NKL) = K
                               PRIMD (NKL) = L
                           ELSE
C                               WRITE (*,*) ' Skip K L = ',K,L
                           END IF

                        END DO
                     END DO
                 ELSE
                     DO K = 1,NPGTOC
                        C = ALPHAC (K)
                        DO L = 1,NPGTOD
                           D = ALPHAD (L)
                           Q = C + D
                           CD = C * D
                           QINV = ONE / Q
                           PQPINV = ONE / (PMIN + Q)
                           RHOCD = DEXP (- CD * RNCDSQ * QINV)
                           T = RMINSQ * PMIN * Q * PQPINV

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

                           SSSSMX = (CD ** ZP75) * RHOCD * SMAXAB * F0
     +                                    * QINV * DSQRT (PQPINV)

                           IF (SSSSMX.GE.TOL) THEN
                               NKL = NKL + 1
                               RHO (NIJ+NKL) = RHOCD
                               PRIMC (NKL) = K
                               PRIMD (NKL) = L
                           ELSE
C                               WRITE (*,*) ' Skip K L = ',K,L
                           END IF

                        END DO
                     END DO
                 END IF
             END IF

         END IF

         IF (NKL.EQ.0) THEN
             EMPTY = .TRUE.
             RETURN
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
