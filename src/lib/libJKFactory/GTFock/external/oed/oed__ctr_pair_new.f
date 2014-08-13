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
         SUBROUTINE  OED__CTR_PAIR_NEW
     +
     +                    ( N,
     +                      NPMAX,NPMIN,
     +                      MIJ,NRS,
     +                      NBLOCK,
     +                      NCR,NCS,
     +                      NPR,NPS,
     +                      CCR,CCS,
     +                      CCBEGR,CCBEGS,
     +                      CCENDR,CCENDS,
     +                      PRIMR,PRIMS,
     +                      EQUALRS,
     +                      SWAPRS,
     +                      PUSED,PSAVE,PPAIR,
     +                      X,W,
     +
     +                             Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__CTR_PAIR_NEW
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation performs a pair contraction step on
C                the incomming integrals over primitives in blocked
C                form over invariant indices n and generates new
C                transformed integrals:
C
C                    y (n,rs) = sum  ccr (r,i) * ccs (s,j) * x (n,ij)
C                                ij
C
C                where ccr and ccs are the arrays containing the
C                contraction coefficients. The sum is over the i and
C                j primitives which are transmitted in arrays PRIMR
C                and PRIMS, respectively, and may constitute only
C                a subset of the full range of primitives.
C
C                The contraction is split into two quarter steps,
C                the order of which is determined by the # of i and
C                j primitives:
C
C                   a) w (n,j/i) = sum ccr/s (r/s,i/j) * x (n,ij)
C                                  i/j
C
C                   b) y (n,rs)  = sum ccs/r (s/r,j/i) * w (n,j/i)
C                                  j/i
C
C                Size of the intermediate w (n,i) or w (n,j) array
C                has to be kept to a minimum and blocking over the
C                invariant indices n has to be performed such that
C                the intermediate w array does not get kicked out
C                from the cache lines after the first quarter
C                transformation.
C
C                In case of csh equality (EQUALRS = .true), we only
C                have to consider the lower triangle of the primitive
C                integrals, which, with the exception of the diagonals,
C                have to be used twice.
C
C                        --- SEGMENTED CONTRACTIONS ---
C
C                Segmented contractions are those defined to be
C                within a certain consecutive i- and j-range of
C                primitives. The segmented limits for each contraction
C                are sitting respectively in CCBEGR and CCBEGS (lowest
C                limit) and CCENDR and CCENDS (highest limit) and they
C                determine which of the actual i's and j's from the
C                PRIMR and PRIMS lists have to be considered for each
C                contraction index.
C
C                The code also allows efficient contractions in case
C                there is only one contraction coefficient present
C                in a certain contraction and its value is equal to 1.
C                In such cases we can save lots of multiplications by
C                1 and since these cases are quite common for certain
C                types of basis functions it is worth including some
C                IF's inside the contraction loops to gain speed.
C
C
C                  Input:
C
C                    N            =  # of invariant indices
C                    NPMAX(MIN)   =  the maximum (minimum) # of
C                                    primitives between both primitive
C                                    sets i,j
C                    MIJ          =  # of ij primitive products to
C                                    be transformed
C                    NRS          =  # of rs contractions to be done
C                    NBLOCK       =  blocking size for invariant
C                                    indices n
C                    NCR(S)       =  # of contractions for the i -> R
C                                    (j -> S) primitives
C                    NPR(S)       =  # of i(j) primitives
C                    CCR(S)       =  full set (including zeros) of
C                                    contraction coefficients for
C                                    R(S) contractions
C                    CCBEGR(S)    =  lowest nonzero primitive i(j)
C                                    index for R(S) contractions
C                    CCENDR(S)    =  highest nonzero primitive i(j)
C                                    index for R(S) contractions
C                    PRIMR(S)     =  primitive i(j) indices
C                    EQUALRS      =  is true, if only the lower
C                                    triangle of ij primitive indices
C                                    is present and consequently
C                                    only the lower triangle of rs
C                                    contractions needs to be evaluated
C                    SWAPRS       =  if this is true, the 1st quarter
C                                    transformation is over R followed
C                                    by the 2nd over S. If false, the
C                                    order is reversed: 1st over S then
C                                    2nd over R
C                    Pxxxx        =  intermediate storage arrays for
C                                    primitive labels to bundle
C                                    contraction steps in do loops
C                                    (xxxx = USED,SAVE,PAIR)
C                    X            =  input array containing the
C                                    primitive integrals
C                    W            =  intermediate storage array
C                                    containing 1st half transformed
C                                    primitive integrals
C
C                  Output:
C
C                    Y            =  contains the new transformed
C                                    integrals
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
         IMPLICIT    NONE

         LOGICAL     CNTRCT
         LOGICAL     EQUALRS
         LOGICAL     IRANGE,JRANGE
         LOGICAL     SWAPRS

         INTEGER     I,J,K,L,M,N,R,S
         INTEGER     I1,I2,I3,I4,I5,I6,I7,I8
         INTEGER     IJ,RS
         INTEGER     IJ1,IJ2,IJ3,IJ4,IJ5,IJ6,IJ7,IJ8
         INTEGER     INEXT,JNEXT
         INTEGER     J1,J2,J3,J4,J5,J6,J7,J8
         INTEGER     MIJ,NRS
         INTEGER     NBLOCK,NSUB
         INTEGER     NCR,NCS
         INTEGER     NI,NJ
         INTEGER     NIBASE,NJBASE
         INTEGER     NILEFT,NJLEFT
         INTEGER     NIREST,NJREST
         INTEGER     NISTEP,NJSTEP
         INTEGER     NPMAX,NPMIN
         INTEGER     NPR,NPS
         INTEGER     PMIN,PMAX

         INTEGER     CCBEGR (1:NCR)
         INTEGER     CCBEGS (1:NCS)
         INTEGER     CCENDR (1:NCR)
         INTEGER     CCENDS (1:NCS)

         INTEGER     PRIMR  (1:MIJ)
         INTEGER     PRIMS  (1:MIJ)

         INTEGER     PPAIR  (1:NPMAX)
         INTEGER     PSAVE  (1:NPMAX)
         INTEGER     PUSED  (1:NPMIN)

         DOUBLE PRECISION    C1,C2,C3,C4,C5,C6,C7,C8
         DOUBLE PRECISION    ZERO,ONE

         DOUBLE PRECISION    CCR  (1:NPR,1:NCR)
         DOUBLE PRECISION    CCS  (1:NPS,1:NCS)

         DOUBLE PRECISION    W (1:NBLOCK,1:NPMIN)

         DOUBLE PRECISION    Y (1:N,1:NRS)
         DOUBLE PRECISION    X (1:N,1:MIJ)

         PARAMETER  (ZERO  =  0.D0)
         PARAMETER  (ONE   =  1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...check case R >= S or full RS range.
C
C
         IF (EQUALRS) THEN
C
C
C             ...the R >= S case. Here we always have equal # of
C                primitives I and J for both R and S. The primitives
C                I >= J are ordered such that J varies fastest.
C                Outer contraction is over S, inner over R.
C
C
           DO 1000 K = 0,N-1,NBLOCK
              NSUB = MIN0 (NBLOCK,N-K)

              RS = 0
              DO 100 S = 1,NCS
                 PMIN = CCBEGS (S)
                 PMAX = CCENDS (S)

                 DO 110 I = 1,NPR
                    PUSED (I) = 0
  110            CONTINUE

                 IF (PMIN.EQ.PMAX .AND. CCS (PMIN,S).EQ.ONE) THEN
C
C
C             ...we reach this point, if only one S contraction
C                coefficient equal to 1.0 is present.
C
C
                     DO 120 IJ = 1,MIJ
                        J = PRIMS (IJ)
                        IF (J.EQ.PMIN) THEN
                            I = PRIMR (IJ)
                            DO L = 1,NSUB
                               W (L,I) = X (K+L,IJ)
                            END DO
                            PUSED (I) = 1
                        END IF
  120                CONTINUE
                 ELSE
C
C
C             ...the general S contraction case with many coefficients.
C                For each I determine all allowed J's and perform the
C                present S contraction for each I over all J's
C                simultaneously.
C
C
                     NJ = 0
                     DO 130 IJ = 1,MIJ
                        I = PRIMR (IJ)
                        J = PRIMS (IJ)
                        JRANGE = J.GE.PMIN .AND. J.LE.PMAX

                        IF (JRANGE) THEN
                            NJ = NJ + 1
                            PSAVE (NJ) = J
                            PPAIR (NJ) = IJ
                        END IF

                        IF (IJ.EQ.MIJ) THEN
                            INEXT = 0
                        ELSE
                            INEXT = PRIMR (IJ+1)
                        END IF

                        CNTRCT = (INEXT.NE.I) .AND. (NJ.GT.0)

                        IF (CNTRCT) THEN
                            IF (NJ.EQ.1) THEN
                                IJ1 = PPAIR (1)
                                C1 = CCS (PSAVE (1),S)
                                DO L = 1,NSUB
                                   W (L,I) = C1 * X (K+L,IJ1)
                                END DO
                            ELSE IF (NJ.EQ.2) THEN
                                IJ1 = PPAIR (1)
                                IJ2 = PPAIR (2)
                                C1 = CCS (PSAVE (1),S)
                                C2 = CCS (PSAVE (2),S)
                                DO L = 1,NSUB
                                   W (L,I) =   C1 * X (K+L,IJ1)
     +                                       + C2 * X (K+L,IJ2)
                                END DO
                            ELSE IF (NJ.EQ.3) THEN
                                IJ1 = PPAIR (1)
                                IJ2 = PPAIR (2)
                                IJ3 = PPAIR (3)
                                C1 = CCS (PSAVE (1),S)
                                C2 = CCS (PSAVE (2),S)
                                C3 = CCS (PSAVE (3),S)
                                DO L = 1,NSUB
                                   W (L,I) =   C1 * X (K+L,IJ1)
     +                                       + C2 * X (K+L,IJ2)
     +                                       + C3 * X (K+L,IJ3)
                                END DO
                            ELSE IF (NJ.EQ.4) THEN
                                IJ1 = PPAIR (1)
                                IJ2 = PPAIR (2)
                                IJ3 = PPAIR (3)
                                IJ4 = PPAIR (4)
                                C1 = CCS (PSAVE (1),S)
                                C2 = CCS (PSAVE (2),S)
                                C3 = CCS (PSAVE (3),S)
                                C4 = CCS (PSAVE (4),S)
                                DO L = 1,NSUB
                                   W (L,I) =   C1 * X (K+L,IJ1)
     +                                       + C2 * X (K+L,IJ2)
     +                                       + C3 * X (K+L,IJ3)
     +                                       + C4 * X (K+L,IJ4)
                                END DO
                            ELSE IF (NJ.EQ.5) THEN
                                IJ1 = PPAIR (1)
                                IJ2 = PPAIR (2)
                                IJ3 = PPAIR (3)
                                IJ4 = PPAIR (4)
                                IJ5 = PPAIR (5)
                                C1 = CCS (PSAVE (1),S)
                                C2 = CCS (PSAVE (2),S)
                                C3 = CCS (PSAVE (3),S)
                                C4 = CCS (PSAVE (4),S)
                                C5 = CCS (PSAVE (5),S)
                                DO L = 1,NSUB
                                   W (L,I) =   C1 * X (K+L,IJ1)
     +                                       + C2 * X (K+L,IJ2)
     +                                       + C3 * X (K+L,IJ3)
     +                                       + C4 * X (K+L,IJ4)
     +                                       + C5 * X (K+L,IJ5)
                                END DO
                            ELSE IF (NJ.EQ.6) THEN
                                IJ1 = PPAIR (1)
                                IJ2 = PPAIR (2)
                                IJ3 = PPAIR (3)
                                IJ4 = PPAIR (4)
                                IJ5 = PPAIR (5)
                                IJ6 = PPAIR (6)
                                C1 = CCS (PSAVE (1),S)
                                C2 = CCS (PSAVE (2),S)
                                C3 = CCS (PSAVE (3),S)
                                C4 = CCS (PSAVE (4),S)
                                C5 = CCS (PSAVE (5),S)
                                C6 = CCS (PSAVE (6),S)
                                DO L = 1,NSUB
                                   W (L,I) =   C1 * X (K+L,IJ1)
     +                                       + C2 * X (K+L,IJ2)
     +                                       + C3 * X (K+L,IJ3)
     +                                       + C4 * X (K+L,IJ4)
     +                                       + C5 * X (K+L,IJ5)
     +                                       + C6 * X (K+L,IJ6)
                                END DO
                            ELSE IF (NJ.EQ.7) THEN
                                IJ1 = PPAIR (1)
                                IJ2 = PPAIR (2)
                                IJ3 = PPAIR (3)
                                IJ4 = PPAIR (4)
                                IJ5 = PPAIR (5)
                                IJ6 = PPAIR (6)
                                IJ7 = PPAIR (7)
                                C1 = CCS (PSAVE (1),S)
                                C2 = CCS (PSAVE (2),S)
                                C3 = CCS (PSAVE (3),S)
                                C4 = CCS (PSAVE (4),S)
                                C5 = CCS (PSAVE (5),S)
                                C6 = CCS (PSAVE (6),S)
                                C7 = CCS (PSAVE (7),S)
                                DO L = 1,NSUB
                                   W (L,I) =   C1 * X (K+L,IJ1)
     +                                       + C2 * X (K+L,IJ2)
     +                                       + C3 * X (K+L,IJ3)
     +                                       + C4 * X (K+L,IJ4)
     +                                       + C5 * X (K+L,IJ5)
     +                                       + C6 * X (K+L,IJ6)
     +                                       + C7 * X (K+L,IJ7)
                                END DO
                            ELSE IF (NJ.EQ.8) THEN
                                IJ1 = PPAIR (1)
                                IJ2 = PPAIR (2)
                                IJ3 = PPAIR (3)
                                IJ4 = PPAIR (4)
                                IJ5 = PPAIR (5)
                                IJ6 = PPAIR (6)
                                IJ7 = PPAIR (7)
                                IJ8 = PPAIR (8)
                                C1 = CCS (PSAVE (1),S)
                                C2 = CCS (PSAVE (2),S)
                                C3 = CCS (PSAVE (3),S)
                                C4 = CCS (PSAVE (4),S)
                                C5 = CCS (PSAVE (5),S)
                                C6 = CCS (PSAVE (6),S)
                                C7 = CCS (PSAVE (7),S)
                                C8 = CCS (PSAVE (8),S)
                                DO L = 1,NSUB
                                   W (L,I) =   C1 * X (K+L,IJ1)
     +                                       + C2 * X (K+L,IJ2)
     +                                       + C3 * X (K+L,IJ3)
     +                                       + C4 * X (K+L,IJ4)
     +                                       + C5 * X (K+L,IJ5)
     +                                       + C6 * X (K+L,IJ6)
     +                                       + C7 * X (K+L,IJ7)
     +                                       + C8 * X (K+L,IJ8)
                                END DO
                            ELSE
                                IJ1 = PPAIR (1)
                                IJ2 = PPAIR (2)
                                IJ3 = PPAIR (3)
                                IJ4 = PPAIR (4)
                                IJ5 = PPAIR (5)
                                IJ6 = PPAIR (6)
                                IJ7 = PPAIR (7)
                                IJ8 = PPAIR (8)
                                C1 = CCS (PSAVE (1),S)
                                C2 = CCS (PSAVE (2),S)
                                C3 = CCS (PSAVE (3),S)
                                C4 = CCS (PSAVE (4),S)
                                C5 = CCS (PSAVE (5),S)
                                C6 = CCS (PSAVE (6),S)
                                C7 = CCS (PSAVE (7),S)
                                C8 = CCS (PSAVE (8),S)
                                DO L = 1,NSUB
                                   W (L,I) =   C1 * X (K+L,IJ1)
     +                                       + C2 * X (K+L,IJ2)
     +                                       + C3 * X (K+L,IJ3)
     +                                       + C4 * X (K+L,IJ4)
     +                                       + C5 * X (K+L,IJ5)
     +                                       + C6 * X (K+L,IJ6)
     +                                       + C7 * X (K+L,IJ7)
     +                                       + C8 * X (K+L,IJ8)
                                END DO

                                NJBASE = 8
                                NJLEFT = NJ - 8
                                NJSTEP = NJLEFT / 8
                                NJREST = MOD (NJLEFT,8)

                                DO M = 1,NJSTEP
                                   IJ1 = PPAIR (NJBASE+1)
                                   IJ2 = PPAIR (NJBASE+2)
                                   IJ3 = PPAIR (NJBASE+3)
                                   IJ4 = PPAIR (NJBASE+4)
                                   IJ5 = PPAIR (NJBASE+5)
                                   IJ6 = PPAIR (NJBASE+6)
                                   IJ7 = PPAIR (NJBASE+7)
                                   IJ8 = PPAIR (NJBASE+8)
                                   C1 = CCS (PSAVE (NJBASE+1),S)
                                   C2 = CCS (PSAVE (NJBASE+2),S)
                                   C3 = CCS (PSAVE (NJBASE+3),S)
                                   C4 = CCS (PSAVE (NJBASE+4),S)
                                   C5 = CCS (PSAVE (NJBASE+5),S)
                                   C6 = CCS (PSAVE (NJBASE+6),S)
                                   C7 = CCS (PSAVE (NJBASE+7),S)
                                   C8 = CCS (PSAVE (NJBASE+8),S)
                                   DO L = 1,NSUB
                                      W (L,I) = W (L,I)
     +                                             + C1 * X (K+L,IJ1)
     +                                             + C2 * X (K+L,IJ2)
     +                                             + C3 * X (K+L,IJ3)
     +                                             + C4 * X (K+L,IJ4)
     +                                             + C5 * X (K+L,IJ5)
     +                                             + C6 * X (K+L,IJ6)
     +                                             + C7 * X (K+L,IJ7)
     +                                             + C8 * X (K+L,IJ8)
                                   END DO
                                   NJBASE = NJBASE + 8
                                END DO

                                DO M = 1,NJREST
                                   IJ1 = PPAIR (NJBASE+M)
                                   C1 = CCS (PSAVE (NJBASE+M),S)
                                   DO L = 1,NSUB
                                      W (L,I) = W (L,I)
     +                                             + C1 * X (K+L,IJ1)
                                   END DO
                                END DO

                            END IF
                            PUSED (I) = 1
C
C
C             ...don't forget to add into W the contributions due
C                to the J indices due to triangularization. Observe
C                that there is no possibility of grouping the J's
C                together, unless we resort the I and J such that
C                I varies fastest. Also observe that only J's
C                will be addressed that are =< I, hence we do
C                not run into the danger of addressing W columns
C                prematurely which have J > I, which saves us
C                from checking addressing of W columns during the
C                above simultaneous J contractions for each I. 
C
C
                            IF (I.GE.PMIN .AND. I.LE.PMAX) THEN
                                C1 = CCS (I,S)
                                DO 135 M = 1,NJ
                                   J = PSAVE (M)
                                   IF (J.NE.I) THEN
                                       IJ1 = PPAIR (M)
                                       IF (PUSED (J).EQ.1) THEN
                                           DO L = 1,NSUB
                                              W (L,J) = W (L,J)
     +                                                + C1 * X (K+L,IJ1)
                                           END DO
                                       ELSE
                                           DO L = 1,NSUB
                                              W (L,J) = C1 * X (K+L,IJ1)
                                           END DO
                                           PUSED (J) = 1
                                       END IF
                                   END IF
  135                           CONTINUE
                            END IF

                            NJ = 0
                        END IF
  130                CONTINUE
                 END IF
C
C
C             ...inner contraction over all R >= S.
C
C
                 DO 140 R = S,NCR
                    RS = RS + 1
                    PMIN = CCBEGR (R)
                    PMAX = CCENDR (R)

                    IF (PMIN.EQ.PMAX .AND. CCR (PMIN,R).EQ.ONE
     +                               .AND. PUSED (PMIN).EQ.1  ) THEN
C
C
C             ...we reach this point, if only one R contraction
C                coefficient equal to 1.0 is present.
C
C
                        DO L = 1,NSUB
                           Y (K+L,RS) = W (L,PMIN)
                        END DO
                    ELSE
C
C
C             ...the general R contraction case with many coefficients.
C                Group all relevant I's together and perform the
C                present R contraction over all I's simultaneously.
C
C
                        NI = 0
                        DO I = PMIN,PMAX
                           IF (PUSED (I).EQ.1) THEN
                               NI = NI + 1
                               PSAVE (NI) = I
                           END IF
                        END DO

                        IF (NI.EQ.0) THEN
                            DO L = 1,NSUB
                               Y (K+L,RS) = ZERO
                            END DO
                        ELSE IF (NI.EQ.1) THEN
                            I1 = PSAVE (1)
                            C1 = CCR (I1,R)
                            DO L = 1,NSUB
                               Y (K+L,RS) = C1 * W (L,I1)
                            END DO
                        ELSE IF (NI.EQ.2) THEN
                            I1 = PSAVE (1)
                            I2 = PSAVE (2)
                            C1 = CCR (I1,R)
                            C2 = CCR (I2,R)
                            DO L = 1,NSUB
                               Y (K+L,RS) =   C1 * W (L,I1)
     +                                      + C2 * W (L,I2)
                            END DO
                        ELSE IF (NI.EQ.3) THEN
                            I1 = PSAVE (1)
                            I2 = PSAVE (2)
                            I3 = PSAVE (3)
                            C1 = CCR (I1,R)
                            C2 = CCR (I2,R)
                            C3 = CCR (I3,R)
                            DO L = 1,NSUB
                               Y (K+L,RS) =   C1 * W (L,I1)
     +                                      + C2 * W (L,I2)
     +                                      + C3 * W (L,I3)
                            END DO
                        ELSE IF (NI.EQ.4) THEN
                            I1 = PSAVE (1)
                            I2 = PSAVE (2)
                            I3 = PSAVE (3)
                            I4 = PSAVE (4)
                            C1 = CCR (I1,R)
                            C2 = CCR (I2,R)
                            C3 = CCR (I3,R)
                            C4 = CCR (I4,R)
                            DO L = 1,NSUB
                               Y (K+L,RS) =   C1 * W (L,I1)
     +                                      + C2 * W (L,I2)
     +                                      + C3 * W (L,I3)
     +                                      + C4 * W (L,I4)
                            END DO
                        ELSE IF (NI.EQ.5) THEN
                            I1 = PSAVE (1)
                            I2 = PSAVE (2)
                            I3 = PSAVE (3)
                            I4 = PSAVE (4)
                            I5 = PSAVE (5)
                            C1 = CCR (I1,R)
                            C2 = CCR (I2,R)
                            C3 = CCR (I3,R)
                            C4 = CCR (I4,R)
                            C5 = CCR (I5,R)
                            DO L = 1,NSUB
                               Y (K+L,RS) =   C1 * W (L,I1)
     +                                      + C2 * W (L,I2)
     +                                      + C3 * W (L,I3)
     +                                      + C4 * W (L,I4)
     +                                      + C5 * W (L,I5)
                            END DO
                        ELSE IF (NI.EQ.6) THEN
                            I1 = PSAVE (1)
                            I2 = PSAVE (2)
                            I3 = PSAVE (3)
                            I4 = PSAVE (4)
                            I5 = PSAVE (5)
                            I6 = PSAVE (6)
                            C1 = CCR (I1,R)
                            C2 = CCR (I2,R)
                            C3 = CCR (I3,R)
                            C4 = CCR (I4,R)
                            C5 = CCR (I5,R)
                            C6 = CCR (I6,R)
                            DO L = 1,NSUB
                               Y (K+L,RS) =   C1 * W (L,I1)
     +                                      + C2 * W (L,I2)
     +                                      + C3 * W (L,I3)
     +                                      + C4 * W (L,I4)
     +                                      + C5 * W (L,I5)
     +                                      + C6 * W (L,I6)
                            END DO
                        ELSE IF (NI.EQ.7) THEN
                            I1 = PSAVE (1)
                            I2 = PSAVE (2)
                            I3 = PSAVE (3)
                            I4 = PSAVE (4)
                            I5 = PSAVE (5)
                            I6 = PSAVE (6)
                            I7 = PSAVE (7)
                            C1 = CCR (I1,R)
                            C2 = CCR (I2,R)
                            C3 = CCR (I3,R)
                            C4 = CCR (I4,R)
                            C5 = CCR (I5,R)
                            C6 = CCR (I6,R)
                            C7 = CCR (I7,R)
                            DO L = 1,NSUB
                               Y (K+L,RS) =   C1 * W (L,I1)
     +                                      + C2 * W (L,I2)
     +                                      + C3 * W (L,I3)
     +                                      + C4 * W (L,I4)
     +                                      + C5 * W (L,I5)
     +                                      + C6 * W (L,I6)
     +                                      + C7 * W (L,I7)
                            END DO
                        ELSE IF (NI.EQ.8) THEN
                            I1 = PSAVE (1)
                            I2 = PSAVE (2)
                            I3 = PSAVE (3)
                            I4 = PSAVE (4)
                            I5 = PSAVE (5)
                            I6 = PSAVE (6)
                            I7 = PSAVE (7)
                            I8 = PSAVE (8)
                            C1 = CCR (I1,R)
                            C2 = CCR (I2,R)
                            C3 = CCR (I3,R)
                            C4 = CCR (I4,R)
                            C5 = CCR (I5,R)
                            C6 = CCR (I6,R)
                            C7 = CCR (I7,R)
                            C8 = CCR (I8,R)
                            DO L = 1,NSUB
                               Y (K+L,RS) =   C1 * W (L,I1)
     +                                      + C2 * W (L,I2)
     +                                      + C3 * W (L,I3)
     +                                      + C4 * W (L,I4)
     +                                      + C5 * W (L,I5)
     +                                      + C6 * W (L,I6)
     +                                      + C7 * W (L,I7)
     +                                      + C8 * W (L,I8)
                            END DO
                        ELSE
                            I1 = PSAVE (1)
                            I2 = PSAVE (2)
                            I3 = PSAVE (3)
                            I4 = PSAVE (4)
                            I5 = PSAVE (5)
                            I6 = PSAVE (6)
                            I7 = PSAVE (7)
                            I8 = PSAVE (8)
                            C1 = CCR (I1,R)
                            C2 = CCR (I2,R)
                            C3 = CCR (I3,R)
                            C4 = CCR (I4,R)
                            C5 = CCR (I5,R)
                            C6 = CCR (I6,R)
                            C7 = CCR (I7,R)
                            C8 = CCR (I8,R)
                            DO L = 1,NSUB
                               Y (K+L,RS) =   C1 * W (L,I1)
     +                                      + C2 * W (L,I2)
     +                                      + C3 * W (L,I3)
     +                                      + C4 * W (L,I4)
     +                                      + C5 * W (L,I5)
     +                                      + C6 * W (L,I6)
     +                                      + C7 * W (L,I7)
     +                                      + C8 * W (L,I8)
                            END DO

                            NIBASE = 8
                            NILEFT = NI - 8
                            NISTEP = NILEFT / 8
                            NIREST = MOD (NILEFT,8)

                            DO M = 1,NISTEP
                               I1 = PSAVE (NIBASE+1)
                               I2 = PSAVE (NIBASE+2)
                               I3 = PSAVE (NIBASE+3)
                               I4 = PSAVE (NIBASE+4)
                               I5 = PSAVE (NIBASE+5)
                               I6 = PSAVE (NIBASE+6)
                               I7 = PSAVE (NIBASE+7)
                               I8 = PSAVE (NIBASE+8)
                               C1 = CCR (I1,R)
                               C2 = CCR (I2,R)
                               C3 = CCR (I3,R)
                               C4 = CCR (I4,R)
                               C5 = CCR (I5,R)
                               C6 = CCR (I6,R)
                               C7 = CCR (I7,R)
                               C8 = CCR (I8,R)
                               DO L = 1,NSUB
                                  Y (K+L,RS) = Y (K+L,RS)
     +                                         + C1 * W (L,I1)
     +                                         + C2 * W (L,I2)
     +                                         + C3 * W (L,I3)
     +                                         + C4 * W (L,I4)
     +                                         + C5 * W (L,I5)
     +                                         + C6 * W (L,I6)
     +                                         + C7 * W (L,I7)
     +                                         + C8 * W (L,I8)
                               END DO
                               NIBASE = NIBASE + 8
                            END DO

                            DO M = 1,NIREST
                               I1 = PSAVE (NIBASE+M)
                               C1 = CCR (I1,R)
                               DO L = 1,NSUB
                                  Y (K+L,RS) = Y (K+L,RS)
     +                                         + C1 * W (L,I1)
                               END DO
                            END DO

                        END IF
                    END IF
  140            CONTINUE
C
C
C             ...next outer S contraction for present block.
C
C
  100         CONTINUE
C
C
C             ...next block (if any).
C
C
 1000      CONTINUE

         ELSE
C
C
C             ...the full RS case. Check the order of the two quarter
C                transformations and proceed accordingly.
C
C                The case: # of I primitives R > # of J primitives S
C                The primitives I and J are ordered such that I varies
C                fastest. Outer contraction is over R, inner over S.
C
C
           IF (SWAPRS) THEN

            DO 2000 K = 0,N-1,NBLOCK
               NSUB = MIN0 (NBLOCK,N-K)

               RS = 0
               DO 200 R = 1,NCR
                  PMIN = CCBEGR (R)
                  PMAX = CCENDR (R)

                  DO 210 J = 1,NPS
                     PUSED (J) = 0
  210             CONTINUE

                  IF (PMIN.EQ.PMAX .AND. CCR (PMIN,R).EQ.ONE) THEN
C
C
C             ...we reach this point, if only one R contraction
C                coefficient equal to 1.0 is present.
C
C
                      DO 220 IJ = 1,MIJ
                         I = PRIMR (IJ)
                         IF (I.EQ.PMIN) THEN
                             J = PRIMS (IJ)
                             DO L = 1,NSUB
                                W (L,J) = X (K+L,IJ)
                             END DO
                             PUSED (J) = 1
                         END IF
  220                 CONTINUE
                  ELSE
C
C
C             ...the general R contraction case with many coefficients.
C                For each J determine all allowed I's and perform the
C                present R contraction for each J over all I's
C                simultaneously.
C
C
                      NI = 0
                      DO 230 IJ = 1,MIJ
                         I = PRIMR (IJ)
                         J = PRIMS (IJ)
                         IRANGE = I.GE.PMIN .AND. I.LE.PMAX

                         IF (IRANGE) THEN
                             NI = NI + 1
                             PSAVE (NI) = I
                             PPAIR (NI) = IJ
                         END IF

                         IF (IJ.EQ.MIJ) THEN
                             JNEXT = 0
                         ELSE
                             JNEXT = PRIMS (IJ+1)
                         END IF

                         CNTRCT = (JNEXT.NE.J) .AND. (NI.GT.0)

                         IF (CNTRCT) THEN
                             IF (NI.EQ.1) THEN
                                 IJ1 = PPAIR (1)
                                 C1 = CCR (PSAVE (1),R)
                                 DO L = 1,NSUB
                                    W (L,J) = C1 * X (K+L,IJ1)
                                 END DO
                             ELSE IF (NI.EQ.2) THEN
                                 IJ1 = PPAIR (1)
                                 IJ2 = PPAIR (2)
                                 C1 = CCR (PSAVE (1),R)
                                 C2 = CCR (PSAVE (2),R)
                                 DO L = 1,NSUB
                                    W (L,J) =   C1 * X (K+L,IJ1)
     +                                        + C2 * X (K+L,IJ2)
                                 END DO
                             ELSE IF (NI.EQ.3) THEN
                                 IJ1 = PPAIR (1)
                                 IJ2 = PPAIR (2)
                                 IJ3 = PPAIR (3)
                                 C1 = CCR (PSAVE (1),R)
                                 C2 = CCR (PSAVE (2),R)
                                 C3 = CCR (PSAVE (3),R)
                                 DO L = 1,NSUB
                                    W (L,J) =   C1 * X (K+L,IJ1)
     +                                        + C2 * X (K+L,IJ2)
     +                                        + C3 * X (K+L,IJ3)
                                 END DO
                             ELSE IF (NI.EQ.4) THEN
                                 IJ1 = PPAIR (1)
                                 IJ2 = PPAIR (2)
                                 IJ3 = PPAIR (3)
                                 IJ4 = PPAIR (4)
                                 C1 = CCR (PSAVE (1),R)
                                 C2 = CCR (PSAVE (2),R)
                                 C3 = CCR (PSAVE (3),R)
                                 C4 = CCR (PSAVE (4),R)
                                 DO L = 1,NSUB
                                    W (L,J) =   C1 * X (K+L,IJ1)
     +                                        + C2 * X (K+L,IJ2)
     +                                        + C3 * X (K+L,IJ3)
     +                                        + C4 * X (K+L,IJ4)
                                 END DO
                             ELSE IF (NI.EQ.5) THEN
                                 IJ1 = PPAIR (1)
                                 IJ2 = PPAIR (2)
                                 IJ3 = PPAIR (3)
                                 IJ4 = PPAIR (4)
                                 IJ5 = PPAIR (5)
                                 C1 = CCR (PSAVE (1),R)
                                 C2 = CCR (PSAVE (2),R)
                                 C3 = CCR (PSAVE (3),R)
                                 C4 = CCR (PSAVE (4),R)
                                 C5 = CCR (PSAVE (5),R)
                                 DO L = 1,NSUB
                                    W (L,J) =   C1 * X (K+L,IJ1)
     +                                        + C2 * X (K+L,IJ2)
     +                                        + C3 * X (K+L,IJ3)
     +                                        + C4 * X (K+L,IJ4)
     +                                        + C5 * X (K+L,IJ5)
                                 END DO
                             ELSE IF (NI.EQ.6) THEN
                                 IJ1 = PPAIR (1)
                                 IJ2 = PPAIR (2)
                                 IJ3 = PPAIR (3)
                                 IJ4 = PPAIR (4)
                                 IJ5 = PPAIR (5)
                                 IJ6 = PPAIR (6)
                                 C1 = CCR (PSAVE (1),R)
                                 C2 = CCR (PSAVE (2),R)
                                 C3 = CCR (PSAVE (3),R)
                                 C4 = CCR (PSAVE (4),R)
                                 C5 = CCR (PSAVE (5),R)
                                 C6 = CCR (PSAVE (6),R)
                                 DO L = 1,NSUB
                                    W (L,J) =   C1 * X (K+L,IJ1)
     +                                        + C2 * X (K+L,IJ2)
     +                                        + C3 * X (K+L,IJ3)
     +                                        + C4 * X (K+L,IJ4)
     +                                        + C5 * X (K+L,IJ5)
     +                                        + C6 * X (K+L,IJ6)
                                 END DO
                             ELSE IF (NI.EQ.7) THEN
                                 IJ1 = PPAIR (1)
                                 IJ2 = PPAIR (2)
                                 IJ3 = PPAIR (3)
                                 IJ4 = PPAIR (4)
                                 IJ5 = PPAIR (5)
                                 IJ6 = PPAIR (6)
                                 IJ7 = PPAIR (7)
                                 C1 = CCR (PSAVE (1),R)
                                 C2 = CCR (PSAVE (2),R)
                                 C3 = CCR (PSAVE (3),R)
                                 C4 = CCR (PSAVE (4),R)
                                 C5 = CCR (PSAVE (5),R)
                                 C6 = CCR (PSAVE (6),R)
                                 C7 = CCR (PSAVE (7),R)
                                 DO L = 1,NSUB
                                    W (L,J) =   C1 * X (K+L,IJ1)
     +                                        + C2 * X (K+L,IJ2)
     +                                        + C3 * X (K+L,IJ3)
     +                                        + C4 * X (K+L,IJ4)
     +                                        + C5 * X (K+L,IJ5)
     +                                        + C6 * X (K+L,IJ6)
     +                                        + C7 * X (K+L,IJ7)
                                 END DO
                             ELSE IF (NI.EQ.8) THEN
                                 IJ1 = PPAIR (1)
                                 IJ2 = PPAIR (2)
                                 IJ3 = PPAIR (3)
                                 IJ4 = PPAIR (4)
                                 IJ5 = PPAIR (5)
                                 IJ6 = PPAIR (6)
                                 IJ7 = PPAIR (7)
                                 IJ8 = PPAIR (8)
                                 C1 = CCR (PSAVE (1),R)
                                 C2 = CCR (PSAVE (2),R)
                                 C3 = CCR (PSAVE (3),R)
                                 C4 = CCR (PSAVE (4),R)
                                 C5 = CCR (PSAVE (5),R)
                                 C6 = CCR (PSAVE (6),R)
                                 C7 = CCR (PSAVE (7),R)
                                 C8 = CCR (PSAVE (8),R)
                                 DO L = 1,NSUB
                                    W (L,J) =   C1 * X (K+L,IJ1)
     +                                        + C2 * X (K+L,IJ2)
     +                                        + C3 * X (K+L,IJ3)
     +                                        + C4 * X (K+L,IJ4)
     +                                        + C5 * X (K+L,IJ5)
     +                                        + C6 * X (K+L,IJ6)
     +                                        + C7 * X (K+L,IJ7)
     +                                        + C8 * X (K+L,IJ8)
                                 END DO
                             ELSE
                                 IJ1 = PPAIR (1)
                                 IJ2 = PPAIR (2)
                                 IJ3 = PPAIR (3)
                                 IJ4 = PPAIR (4)
                                 IJ5 = PPAIR (5)
                                 IJ6 = PPAIR (6)
                                 IJ7 = PPAIR (7)
                                 IJ8 = PPAIR (8)
                                 C1 = CCR (PSAVE (1),R)
                                 C2 = CCR (PSAVE (2),R)
                                 C3 = CCR (PSAVE (3),R)
                                 C4 = CCR (PSAVE (4),R)
                                 C5 = CCR (PSAVE (5),R)
                                 C6 = CCR (PSAVE (6),R)
                                 C7 = CCR (PSAVE (7),R)
                                 C8 = CCR (PSAVE (8),R)
                                 DO L = 1,NSUB
                                    W (L,J) =   C1 * X (K+L,IJ1)
     +                                        + C2 * X (K+L,IJ2)
     +                                        + C3 * X (K+L,IJ3)
     +                                        + C4 * X (K+L,IJ4)
     +                                        + C5 * X (K+L,IJ5)
     +                                        + C6 * X (K+L,IJ6)
     +                                        + C7 * X (K+L,IJ7)
     +                                        + C8 * X (K+L,IJ8)
                                 END DO

                                 NIBASE = 8
                                 NILEFT = NI - 8
                                 NISTEP = NILEFT / 8
                                 NIREST = MOD (NILEFT,8)

                                 DO M = 1,NISTEP
                                    IJ1 = PPAIR (NIBASE+1)
                                    IJ2 = PPAIR (NIBASE+2)
                                    IJ3 = PPAIR (NIBASE+3)
                                    IJ4 = PPAIR (NIBASE+4)
                                    IJ5 = PPAIR (NIBASE+5)
                                    IJ6 = PPAIR (NIBASE+6)
                                    IJ7 = PPAIR (NIBASE+7)
                                    IJ8 = PPAIR (NIBASE+8)
                                    C1 = CCR (PSAVE (NIBASE+1),R)
                                    C2 = CCR (PSAVE (NIBASE+2),R)
                                    C3 = CCR (PSAVE (NIBASE+3),R)
                                    C4 = CCR (PSAVE (NIBASE+4),R)
                                    C5 = CCR (PSAVE (NIBASE+5),R)
                                    C6 = CCR (PSAVE (NIBASE+6),R)
                                    C7 = CCR (PSAVE (NIBASE+7),R)
                                    C8 = CCR (PSAVE (NIBASE+8),R)
                                    DO L = 1,NSUB
                                       W (L,J) = W (L,J)
     +                                              + C1 * X (K+L,IJ1)
     +                                              + C2 * X (K+L,IJ2)
     +                                              + C3 * X (K+L,IJ3)
     +                                              + C4 * X (K+L,IJ4)
     +                                              + C5 * X (K+L,IJ5)
     +                                              + C6 * X (K+L,IJ6)
     +                                              + C7 * X (K+L,IJ7)
     +                                              + C8 * X (K+L,IJ8)
                                    END DO
                                    NIBASE = NIBASE + 8
                                 END DO

                                 DO M = 1,NIREST
                                    IJ1 = PPAIR (NIBASE+M)
                                    C1 = CCR (PSAVE (NIBASE+M),R)
                                    DO L = 1,NSUB
                                       W (L,J) = W (L,J)
     +                                              + C1 * X (K+L,IJ1)
                                    END DO
                                 END DO

                             END IF
                             PUSED (J) = 1
                             NI = 0
                         END IF

  230                 CONTINUE
                  END IF
C
C
C             ...inner contraction over all S.
C
C
                  DO 240 S = 1,NCS
                     RS = RS + 1
                     PMIN = CCBEGS (S)
                     PMAX = CCENDS (S)

                     IF (PMIN.EQ.PMAX .AND. CCS (PMIN,S).EQ.ONE
     +                                .AND. PUSED (PMIN).EQ.1  ) THEN
C
C
C             ...we reach this point, if only one S contraction
C                coefficient equal to 1.0 is present.
C
C
                         DO L = 1,NSUB
                            Y (K+L,RS) = W (L,PMIN)
                         END DO
                     ELSE
C
C
C             ...the general S contraction case with many coefficients.
C                Group all relevant J's together and perform the
C                present S contraction over all J's simultaneously.
C
C
                         NJ = 0
                         DO J = PMIN,PMAX
                            IF (PUSED (J).EQ.1) THEN
                                NJ = NJ + 1
                                PSAVE (NJ) = J
                            END IF
                         END DO

                         IF (NJ.EQ.0) THEN
                             DO L = 1,NSUB
                                Y (K+L,RS) = ZERO
                             END DO
                         ELSE IF (NJ.EQ.1) THEN
                             J1 = PSAVE (1)
                             C1 = CCS (J1,S)
                             DO L = 1,NSUB
                                Y (K+L,RS) = C1 * W (L,J1)
                             END DO
                         ELSE IF (NJ.EQ.2) THEN
                             J1 = PSAVE (1)
                             J2 = PSAVE (2)
                             C1 = CCS (J1,S)
                             C2 = CCS (J2,S)
                             DO L = 1,NSUB
                                Y (K+L,RS) =   C1 * W (L,J1)
     +                                       + C2 * W (L,J2)
                             END DO
                         ELSE IF (NJ.EQ.3) THEN
                             J1 = PSAVE (1)
                             J2 = PSAVE (2)
                             J3 = PSAVE (3)
                             C1 = CCS (J1,S)
                             C2 = CCS (J2,S)
                             C3 = CCS (J3,S)
                             DO L = 1,NSUB
                                Y (K+L,RS) =   C1 * W (L,J1)
     +                                       + C2 * W (L,J2)
     +                                       + C3 * W (L,J3)
                             END DO
                         ELSE IF (NJ.EQ.4) THEN
                             J1 = PSAVE (1)
                             J2 = PSAVE (2)
                             J3 = PSAVE (3)
                             J4 = PSAVE (4)
                             C1 = CCS (J1,S)
                             C2 = CCS (J2,S)
                             C3 = CCS (J3,S)
                             C4 = CCS (J4,S)
                             DO L = 1,NSUB
                                Y (K+L,RS) =   C1 * W (L,J1)
     +                                       + C2 * W (L,J2)
     +                                       + C3 * W (L,J3)
     +                                       + C4 * W (L,J4)
                             END DO
                         ELSE IF (NJ.EQ.5) THEN
                             J1 = PSAVE (1)
                             J2 = PSAVE (2)
                             J3 = PSAVE (3)
                             J4 = PSAVE (4)
                             J5 = PSAVE (5)
                             C1 = CCS (J1,S)
                             C2 = CCS (J2,S)
                             C3 = CCS (J3,S)
                             C4 = CCS (J4,S)
                             C5 = CCS (J5,S)
                             DO L = 1,NSUB
                                Y (K+L,RS) =   C1 * W (L,J1)
     +                                       + C2 * W (L,J2)
     +                                       + C3 * W (L,J3)
     +                                       + C4 * W (L,J4)
     +                                       + C5 * W (L,J5)
                             END DO
                         ELSE IF (NJ.EQ.6) THEN
                             J1 = PSAVE (1)
                             J2 = PSAVE (2)
                             J3 = PSAVE (3)
                             J4 = PSAVE (4)
                             J5 = PSAVE (5)
                             J6 = PSAVE (6)
                             C1 = CCS (J1,S)
                             C2 = CCS (J2,S)
                             C3 = CCS (J3,S)
                             C4 = CCS (J4,S)
                             C5 = CCS (J5,S)
                             C6 = CCS (J6,S)
                             DO L = 1,NSUB
                                Y (K+L,RS) =   C1 * W (L,J1)
     +                                       + C2 * W (L,J2)
     +                                       + C3 * W (L,J3)
     +                                       + C4 * W (L,J4)
     +                                       + C5 * W (L,J5)
     +                                       + C6 * W (L,J6)
                             END DO
                         ELSE IF (NJ.EQ.7) THEN
                             J1 = PSAVE (1)
                             J2 = PSAVE (2)
                             J3 = PSAVE (3)
                             J4 = PSAVE (4)
                             J5 = PSAVE (5)
                             J6 = PSAVE (6)
                             J7 = PSAVE (7)
                             C1 = CCS (J1,S)
                             C2 = CCS (J2,S)
                             C3 = CCS (J3,S)
                             C4 = CCS (J4,S)
                             C5 = CCS (J5,S)
                             C6 = CCS (J6,S)
                             C7 = CCS (J7,S)
                             DO L = 1,NSUB
                                Y (K+L,RS) =   C1 * W (L,J1)
     +                                       + C2 * W (L,J2)
     +                                       + C3 * W (L,J3)
     +                                       + C4 * W (L,J4)
     +                                       + C5 * W (L,J5)
     +                                       + C6 * W (L,J6)
     +                                       + C7 * W (L,J7)
                             END DO
                         ELSE IF (NJ.EQ.8) THEN
                             J1 = PSAVE (1)
                             J2 = PSAVE (2)
                             J3 = PSAVE (3)
                             J4 = PSAVE (4)
                             J5 = PSAVE (5)
                             J6 = PSAVE (6)
                             J7 = PSAVE (7)
                             J8 = PSAVE (8)
                             C1 = CCS (J1,S)
                             C2 = CCS (J2,S)
                             C3 = CCS (J3,S)
                             C4 = CCS (J4,S)
                             C5 = CCS (J5,S)
                             C6 = CCS (J6,S)
                             C7 = CCS (J7,S)
                             C8 = CCS (J8,S)
                             DO L = 1,NSUB
                                Y (K+L,RS) =   C1 * W (L,J1)
     +                                       + C2 * W (L,J2)
     +                                       + C3 * W (L,J3)
     +                                       + C4 * W (L,J4)
     +                                       + C5 * W (L,J5)
     +                                       + C6 * W (L,J6)
     +                                       + C7 * W (L,J7)
     +                                       + C8 * W (L,J8)
                             END DO
                         ELSE
                             J1 = PSAVE (1)
                             J2 = PSAVE (2)
                             J3 = PSAVE (3)
                             J4 = PSAVE (4)
                             J5 = PSAVE (5)
                             J6 = PSAVE (6)
                             J7 = PSAVE (7)
                             J8 = PSAVE (8)
                             C1 = CCS (J1,S)
                             C2 = CCS (J2,S)
                             C3 = CCS (J3,S)
                             C4 = CCS (J4,S)
                             C5 = CCS (J5,S)
                             C6 = CCS (J6,S)
                             C7 = CCS (J7,S)
                             C8 = CCS (J8,S)
                             DO L = 1,NSUB
                                Y (K+L,RS) =   C1 * W (L,J1)
     +                                       + C2 * W (L,J2)
     +                                       + C3 * W (L,J3)
     +                                       + C4 * W (L,J4)
     +                                       + C5 * W (L,J5)
     +                                       + C6 * W (L,J6)
     +                                       + C7 * W (L,J7)
     +                                       + C8 * W (L,J8)
                             END DO

                             NJBASE = 8
                             NJLEFT = NJ - 8
                             NJSTEP = NJLEFT / 8
                             NJREST = MOD (NJLEFT,8)

                             DO M = 1,NJSTEP
                                J1 = PSAVE (NJBASE+1)
                                J2 = PSAVE (NJBASE+2)
                                J3 = PSAVE (NJBASE+3)
                                J4 = PSAVE (NJBASE+4)
                                J5 = PSAVE (NJBASE+5)
                                J6 = PSAVE (NJBASE+6)
                                J7 = PSAVE (NJBASE+7)
                                J8 = PSAVE (NJBASE+8)
                                C1 = CCS (J1,S)
                                C2 = CCS (J2,S)
                                C3 = CCS (J3,S)
                                C4 = CCS (J4,S)
                                C5 = CCS (J5,S)
                                C6 = CCS (J6,S)
                                C7 = CCS (J7,S)
                                C8 = CCS (J8,S)
                                DO L = 1,NSUB
                                   Y (K+L,RS) = Y (K+L,RS)
     +                                          + C1 * W (L,J1)
     +                                          + C2 * W (L,J2)
     +                                          + C3 * W (L,J3)
     +                                          + C4 * W (L,J4)
     +                                          + C5 * W (L,J5)
     +                                          + C6 * W (L,J6)
     +                                          + C7 * W (L,J7)
     +                                          + C8 * W (L,J8)
                                END DO
                                NJBASE = NJBASE + 8
                             END DO

                             DO M = 1,NJREST
                                J1 = PSAVE (NJBASE+M)
                                C1 = CCS (J1,S)
                                DO L = 1,NSUB
                                   Y (K+L,RS) = Y (K+L,RS)
     +                                          + C1 * W (L,J1)
                                END DO
                             END DO

                         END IF
                     END IF
  240             CONTINUE
C
C
C             ...next outer R contraction for present block.
C
C
  200          CONTINUE
C
C
C             ...next block (if any).
C
C
 2000       CONTINUE

           ELSE
C
C
C             ...the case: # of I primitives R =< # of J primitives S.
C                The primitives I and J are ordered such that J varies
C                fastest. Outer contraction is over S, inner over R.
C
C
            DO 3000 K = 0,N-1,NBLOCK
               NSUB = MIN0 (NBLOCK,N-K)

               RS = 0
               DO 300 S = 1,NCS
                  PMIN = CCBEGS (S)
                  PMAX = CCENDS (S)

                  DO 310 I = 1,NPR
                     PUSED (I) = 0
  310             CONTINUE

                  IF (PMIN.EQ.PMAX .AND. CCS (PMIN,S).EQ.ONE) THEN
C
C
C             ...we reach this point, if only one S contraction
C                coefficient equal to 1.0 is present.
C
C
                      DO 320 IJ = 1,MIJ
                         J = PRIMS (IJ)
                         IF (J.EQ.PMIN) THEN
                             I = PRIMR (IJ)
                             DO L = 1,NSUB
                                W (L,I) = X (K+L,IJ)
                             END DO
                             PUSED (I) = 1
                         END IF
  320                 CONTINUE
                  ELSE
C
C
C             ...the general S contraction case with many coefficients.
C                For each I determine all allowed J's and perform the
C                present S contraction for each I over all J's
C                simultaneously.
C
C
                      NJ = 0
                      DO 330 IJ = 1,MIJ
                         I = PRIMR (IJ)
                         J = PRIMS (IJ)
                         JRANGE = J.GE.PMIN .AND. J.LE.PMAX

                         IF (JRANGE) THEN
                             NJ = NJ + 1
                             PSAVE (NJ) = J
                             PPAIR (NJ) = IJ
                         END IF

                         IF (IJ.EQ.MIJ) THEN
                             INEXT = 0
                         ELSE
                             INEXT = PRIMR (IJ+1)
                         END IF

                         CNTRCT = (INEXT.NE.I) .AND. (NJ.GT.0)

                         IF (CNTRCT) THEN
                             IF (NJ.EQ.1) THEN
                                 IJ1 = PPAIR (1)
                                 C1 = CCS (PSAVE (1),S)
                                 DO L = 1,NSUB
                                    W (L,I) = C1 * X (K+L,IJ1)
                                 END DO
                             ELSE IF (NJ.EQ.2) THEN
                                 IJ1 = PPAIR (1)
                                 IJ2 = PPAIR (2)
                                 C1 = CCS (PSAVE (1),S)
                                 C2 = CCS (PSAVE (2),S)
                                 DO L = 1,NSUB
                                    W (L,I) =   C1 * X (K+L,IJ1)
     +                                        + C2 * X (K+L,IJ2)
                                 END DO
                             ELSE IF (NJ.EQ.3) THEN
                                 IJ1 = PPAIR (1)
                                 IJ2 = PPAIR (2)
                                 IJ3 = PPAIR (3)
                                 C1 = CCS (PSAVE (1),S)
                                 C2 = CCS (PSAVE (2),S)
                                 C3 = CCS (PSAVE (3),S)
                                 DO L = 1,NSUB
                                    W (L,I) =   C1 * X (K+L,IJ1)
     +                                        + C2 * X (K+L,IJ2)
     +                                        + C3 * X (K+L,IJ3)
                                 END DO
                             ELSE IF (NJ.EQ.4) THEN
                                 IJ1 = PPAIR (1)
                                 IJ2 = PPAIR (2)
                                 IJ3 = PPAIR (3)
                                 IJ4 = PPAIR (4)
                                 C1 = CCS (PSAVE (1),S)
                                 C2 = CCS (PSAVE (2),S)
                                 C3 = CCS (PSAVE (3),S)
                                 C4 = CCS (PSAVE (4),S)
                                 DO L = 1,NSUB
                                    W (L,I) =   C1 * X (K+L,IJ1)
     +                                        + C2 * X (K+L,IJ2)
     +                                        + C3 * X (K+L,IJ3)
     +                                        + C4 * X (K+L,IJ4)
                                 END DO
                             ELSE IF (NJ.EQ.5) THEN
                                 IJ1 = PPAIR (1)
                                 IJ2 = PPAIR (2)
                                 IJ3 = PPAIR (3)
                                 IJ4 = PPAIR (4)
                                 IJ5 = PPAIR (5)
                                 C1 = CCS (PSAVE (1),S)
                                 C2 = CCS (PSAVE (2),S)
                                 C3 = CCS (PSAVE (3),S)
                                 C4 = CCS (PSAVE (4),S)
                                 C5 = CCS (PSAVE (5),S)
                                 DO L = 1,NSUB
                                    W (L,I) =   C1 * X (K+L,IJ1)
     +                                        + C2 * X (K+L,IJ2)
     +                                        + C3 * X (K+L,IJ3)
     +                                        + C4 * X (K+L,IJ4)
     +                                        + C5 * X (K+L,IJ5)
                                 END DO
                             ELSE IF (NJ.EQ.6) THEN
                                 IJ1 = PPAIR (1)
                                 IJ2 = PPAIR (2)
                                 IJ3 = PPAIR (3)
                                 IJ4 = PPAIR (4)
                                 IJ5 = PPAIR (5)
                                 IJ6 = PPAIR (6)
                                 C1 = CCS (PSAVE (1),S)
                                 C2 = CCS (PSAVE (2),S)
                                 C3 = CCS (PSAVE (3),S)
                                 C4 = CCS (PSAVE (4),S)
                                 C5 = CCS (PSAVE (5),S)
                                 C6 = CCS (PSAVE (6),S)
                                 DO L = 1,NSUB
                                    W (L,I) =   C1 * X (K+L,IJ1)
     +                                        + C2 * X (K+L,IJ2)
     +                                        + C3 * X (K+L,IJ3)
     +                                        + C4 * X (K+L,IJ4)
     +                                        + C5 * X (K+L,IJ5)
     +                                        + C6 * X (K+L,IJ6)
                                 END DO
                             ELSE IF (NJ.EQ.7) THEN
                                 IJ1 = PPAIR (1)
                                 IJ2 = PPAIR (2)
                                 IJ3 = PPAIR (3)
                                 IJ4 = PPAIR (4)
                                 IJ5 = PPAIR (5)
                                 IJ6 = PPAIR (6)
                                 IJ7 = PPAIR (7)
                                 C1 = CCS (PSAVE (1),S)
                                 C2 = CCS (PSAVE (2),S)
                                 C3 = CCS (PSAVE (3),S)
                                 C4 = CCS (PSAVE (4),S)
                                 C5 = CCS (PSAVE (5),S)
                                 C6 = CCS (PSAVE (6),S)
                                 C7 = CCS (PSAVE (7),S)
                                 DO L = 1,NSUB
                                    W (L,I) =   C1 * X (K+L,IJ1)
     +                                        + C2 * X (K+L,IJ2)
     +                                        + C3 * X (K+L,IJ3)
     +                                        + C4 * X (K+L,IJ4)
     +                                        + C5 * X (K+L,IJ5)
     +                                        + C6 * X (K+L,IJ6)
     +                                        + C7 * X (K+L,IJ7)
                                 END DO
                             ELSE IF (NJ.EQ.8) THEN
                                 IJ1 = PPAIR (1)
                                 IJ2 = PPAIR (2)
                                 IJ3 = PPAIR (3)
                                 IJ4 = PPAIR (4)
                                 IJ5 = PPAIR (5)
                                 IJ6 = PPAIR (6)
                                 IJ7 = PPAIR (7)
                                 IJ8 = PPAIR (8)
                                 C1 = CCS (PSAVE (1),S)
                                 C2 = CCS (PSAVE (2),S)
                                 C3 = CCS (PSAVE (3),S)
                                 C4 = CCS (PSAVE (4),S)
                                 C5 = CCS (PSAVE (5),S)
                                 C6 = CCS (PSAVE (6),S)
                                 C7 = CCS (PSAVE (7),S)
                                 C8 = CCS (PSAVE (8),S)
                                 DO L = 1,NSUB
                                    W (L,I) =   C1 * X (K+L,IJ1)
     +                                        + C2 * X (K+L,IJ2)
     +                                        + C3 * X (K+L,IJ3)
     +                                        + C4 * X (K+L,IJ4)
     +                                        + C5 * X (K+L,IJ5)
     +                                        + C6 * X (K+L,IJ6)
     +                                        + C7 * X (K+L,IJ7)
     +                                        + C8 * X (K+L,IJ8)
                                 END DO
                             ELSE
                                 IJ1 = PPAIR (1)
                                 IJ2 = PPAIR (2)
                                 IJ3 = PPAIR (3)
                                 IJ4 = PPAIR (4)
                                 IJ5 = PPAIR (5)
                                 IJ6 = PPAIR (6)
                                 IJ7 = PPAIR (7)
                                 IJ8 = PPAIR (8)
                                 C1 = CCS (PSAVE (1),S)
                                 C2 = CCS (PSAVE (2),S)
                                 C3 = CCS (PSAVE (3),S)
                                 C4 = CCS (PSAVE (4),S)
                                 C5 = CCS (PSAVE (5),S)
                                 C6 = CCS (PSAVE (6),S)
                                 C7 = CCS (PSAVE (7),S)
                                 C8 = CCS (PSAVE (8),S)
                                 DO L = 1,NSUB
                                    W (L,I) =   C1 * X (K+L,IJ1)
     +                                        + C2 * X (K+L,IJ2)
     +                                        + C3 * X (K+L,IJ3)
     +                                        + C4 * X (K+L,IJ4)
     +                                        + C5 * X (K+L,IJ5)
     +                                        + C6 * X (K+L,IJ6)
     +                                        + C7 * X (K+L,IJ7)
     +                                        + C8 * X (K+L,IJ8)
                                 END DO

                                 NJBASE = 8
                                 NJLEFT = NJ - 8
                                 NJSTEP = NJLEFT / 8
                                 NJREST = MOD (NJLEFT,8)

                                 DO M = 1,NJSTEP
                                    IJ1 = PPAIR (NJBASE+1)
                                    IJ2 = PPAIR (NJBASE+2)
                                    IJ3 = PPAIR (NJBASE+3)
                                    IJ4 = PPAIR (NJBASE+4)
                                    IJ5 = PPAIR (NJBASE+5)
                                    IJ6 = PPAIR (NJBASE+6)
                                    IJ7 = PPAIR (NJBASE+7)
                                    IJ8 = PPAIR (NJBASE+8)
                                    C1 = CCS (PSAVE (NJBASE+1),S)
                                    C2 = CCS (PSAVE (NJBASE+2),S)
                                    C3 = CCS (PSAVE (NJBASE+3),S)
                                    C4 = CCS (PSAVE (NJBASE+4),S)
                                    C5 = CCS (PSAVE (NJBASE+5),S)
                                    C6 = CCS (PSAVE (NJBASE+6),S)
                                    C7 = CCS (PSAVE (NJBASE+7),S)
                                    C8 = CCS (PSAVE (NJBASE+8),S)
                                    DO L = 1,NSUB
                                       W (L,I) = W (L,I)
     +                                              + C1 * X (K+L,IJ1)
     +                                              + C2 * X (K+L,IJ2)
     +                                              + C3 * X (K+L,IJ3)
     +                                              + C4 * X (K+L,IJ4)
     +                                              + C5 * X (K+L,IJ5)
     +                                              + C6 * X (K+L,IJ6)
     +                                              + C7 * X (K+L,IJ7)
     +                                              + C8 * X (K+L,IJ8)
                                    END DO
                                    NJBASE = NJBASE + 8
                                 END DO

                                 DO M = 1,NJREST
                                    IJ1 = PPAIR (NJBASE+M)
                                    C1 = CCS (PSAVE (NJBASE+M),S)
                                    DO L = 1,NSUB
                                       W (L,I) = W (L,I)
     +                                              + C1 * X (K+L,IJ1)
                                    END DO
                                 END DO

                             END IF
                             PUSED (I) = 1
                             NJ = 0
                         END IF

  330                 CONTINUE
                  END IF
C
C
C             ...inner contraction over all R.
C
C
                  DO 340 R = 1,NCR
                     RS = RS + 1
                     PMIN = CCBEGR (R)
                     PMAX = CCENDR (R)

                     IF (PMIN.EQ.PMAX .AND. CCR (PMIN,R).EQ.ONE
     +                                .AND. PUSED (PMIN).EQ.1  ) THEN
C
C
C             ...we reach this point, if only one R contraction
C                coefficient equal to 1.0 is present.
C
C
                         DO L = 1,NSUB
                            Y (K+L,RS) = W (L,PMIN)
                         END DO
                     ELSE
C
C
C             ...the general R contraction case with many coefficients.
C                Group all relevant I's together and perform the
C                present R contraction over all I's simultaneously.
C
C
                         NI = 0
                         DO I = PMIN,PMAX
                            IF (PUSED (I).EQ.1) THEN
                                NI = NI + 1
                                PSAVE (NI) = I
                            END IF
                         END DO

                         IF (NI.EQ.0) THEN
                             DO L = 1,NSUB
                                Y (K+L,RS) = ZERO
                             END DO
                         ELSE IF (NI.EQ.1) THEN
                             I1 = PSAVE (1)
                             C1 = CCR (I1,R)
                             DO L = 1,NSUB
                                Y (K+L,RS) = C1 * W (L,I1)
                             END DO
                         ELSE IF (NI.EQ.2) THEN
                             I1 = PSAVE (1)
                             I2 = PSAVE (2)
                             C1 = CCR (I1,R)
                             C2 = CCR (I2,R)
                             DO L = 1,NSUB
                                Y (K+L,RS) =   C1 * W (L,I1)
     +                                       + C2 * W (L,I2)
                             END DO
                         ELSE IF (NI.EQ.3) THEN
                             I1 = PSAVE (1)
                             I2 = PSAVE (2)
                             I3 = PSAVE (3)
                             C1 = CCR (I1,R)
                             C2 = CCR (I2,R)
                             C3 = CCR (I3,R)
                             DO L = 1,NSUB
                                Y (K+L,RS) =   C1 * W (L,I1)
     +                                       + C2 * W (L,I2)
     +                                       + C3 * W (L,I3)
                             END DO
                         ELSE IF (NI.EQ.4) THEN
                             I1 = PSAVE (1)
                             I2 = PSAVE (2)
                             I3 = PSAVE (3)
                             I4 = PSAVE (4)
                             C1 = CCR (I1,R)
                             C2 = CCR (I2,R)
                             C3 = CCR (I3,R)
                             C4 = CCR (I4,R)
                             DO L = 1,NSUB
                                Y (K+L,RS) =   C1 * W (L,I1)
     +                                       + C2 * W (L,I2)
     +                                       + C3 * W (L,I3)
     +                                       + C4 * W (L,I4)
                             END DO
                         ELSE IF (NI.EQ.5) THEN
                             I1 = PSAVE (1)
                             I2 = PSAVE (2)
                             I3 = PSAVE (3)
                             I4 = PSAVE (4)
                             I5 = PSAVE (5)
                             C1 = CCR (I1,R)
                             C2 = CCR (I2,R)
                             C3 = CCR (I3,R)
                             C4 = CCR (I4,R)
                             C5 = CCR (I5,R)
                             DO L = 1,NSUB
                                Y (K+L,RS) =   C1 * W (L,I1)
     +                                       + C2 * W (L,I2)
     +                                       + C3 * W (L,I3)
     +                                       + C4 * W (L,I4)
     +                                       + C5 * W (L,I5)
                             END DO
                         ELSE IF (NI.EQ.6) THEN
                             I1 = PSAVE (1)
                             I2 = PSAVE (2)
                             I3 = PSAVE (3)
                             I4 = PSAVE (4)
                             I5 = PSAVE (5)
                             I6 = PSAVE (6)
                             C1 = CCR (I1,R)
                             C2 = CCR (I2,R)
                             C3 = CCR (I3,R)
                             C4 = CCR (I4,R)
                             C5 = CCR (I5,R)
                             C6 = CCR (I6,R)
                             DO L = 1,NSUB
                                Y (K+L,RS) =   C1 * W (L,I1)
     +                                       + C2 * W (L,I2)
     +                                       + C3 * W (L,I3)
     +                                       + C4 * W (L,I4)
     +                                       + C5 * W (L,I5)
     +                                       + C6 * W (L,I6)
                             END DO
                         ELSE IF (NI.EQ.7) THEN
                             I1 = PSAVE (1)
                             I2 = PSAVE (2)
                             I3 = PSAVE (3)
                             I4 = PSAVE (4)
                             I5 = PSAVE (5)
                             I6 = PSAVE (6)
                             I7 = PSAVE (7)
                             C1 = CCR (I1,R)
                             C2 = CCR (I2,R)
                             C3 = CCR (I3,R)
                             C4 = CCR (I4,R)
                             C5 = CCR (I5,R)
                             C6 = CCR (I6,R)
                             C7 = CCR (I7,R)
                             DO L = 1,NSUB
                                Y (K+L,RS) =   C1 * W (L,I1)
     +                                       + C2 * W (L,I2)
     +                                       + C3 * W (L,I3)
     +                                       + C4 * W (L,I4)
     +                                       + C5 * W (L,I5)
     +                                       + C6 * W (L,I6)
     +                                       + C7 * W (L,I7)
                             END DO
                         ELSE IF (NI.EQ.8) THEN
                             I1 = PSAVE (1)
                             I2 = PSAVE (2)
                             I3 = PSAVE (3)
                             I4 = PSAVE (4)
                             I5 = PSAVE (5)
                             I6 = PSAVE (6)
                             I7 = PSAVE (7)
                             I8 = PSAVE (8)
                             C1 = CCR (I1,R)
                             C2 = CCR (I2,R)
                             C3 = CCR (I3,R)
                             C4 = CCR (I4,R)
                             C5 = CCR (I5,R)
                             C6 = CCR (I6,R)
                             C7 = CCR (I7,R)
                             C8 = CCR (I8,R)
                             DO L = 1,NSUB
                                Y (K+L,RS) =   C1 * W (L,I1)
     +                                       + C2 * W (L,I2)
     +                                       + C3 * W (L,I3)
     +                                       + C4 * W (L,I4)
     +                                       + C5 * W (L,I5)
     +                                       + C6 * W (L,I6)
     +                                       + C7 * W (L,I7)
     +                                       + C8 * W (L,I8)
                             END DO
                         ELSE
                             I1 = PSAVE (1)
                             I2 = PSAVE (2)
                             I3 = PSAVE (3)
                             I4 = PSAVE (4)
                             I5 = PSAVE (5)
                             I6 = PSAVE (6)
                             I7 = PSAVE (7)
                             I8 = PSAVE (8)
                             C1 = CCR (I1,R)
                             C2 = CCR (I2,R)
                             C3 = CCR (I3,R)
                             C4 = CCR (I4,R)
                             C5 = CCR (I5,R)
                             C6 = CCR (I6,R)
                             C7 = CCR (I7,R)
                             C8 = CCR (I8,R)
                             DO L = 1,NSUB
                                Y (K+L,RS) =   C1 * W (L,I1)
     +                                       + C2 * W (L,I2)
     +                                       + C3 * W (L,I3)
     +                                       + C4 * W (L,I4)
     +                                       + C5 * W (L,I5)
     +                                       + C6 * W (L,I6)
     +                                       + C7 * W (L,I7)
     +                                       + C8 * W (L,I8)
                             END DO

                             NIBASE = 8
                             NILEFT = NI - 8
                             NISTEP = NILEFT / 8
                             NIREST = MOD (NILEFT,8)

                             DO M = 1,NISTEP
                                I1 = PSAVE (NIBASE+1)
                                I2 = PSAVE (NIBASE+2)
                                I3 = PSAVE (NIBASE+3)
                                I4 = PSAVE (NIBASE+4)
                                I5 = PSAVE (NIBASE+5)
                                I6 = PSAVE (NIBASE+6)
                                I7 = PSAVE (NIBASE+7)
                                I8 = PSAVE (NIBASE+8)
                                C1 = CCR (I1,R)
                                C2 = CCR (I2,R)
                                C3 = CCR (I3,R)
                                C4 = CCR (I4,R)
                                C5 = CCR (I5,R)
                                C6 = CCR (I6,R)
                                C7 = CCR (I7,R)
                                C8 = CCR (I8,R)
                                DO L = 1,NSUB
                                   Y (K+L,RS) = Y (K+L,RS)
     +                                          + C1 * W (L,I1)
     +                                          + C2 * W (L,I2)
     +                                          + C3 * W (L,I3)
     +                                          + C4 * W (L,I4)
     +                                          + C5 * W (L,I5)
     +                                          + C6 * W (L,I6)
     +                                          + C7 * W (L,I7)
     +                                          + C8 * W (L,I8)
                                END DO
                                NIBASE = NIBASE + 8
                             END DO

                             DO M = 1,NIREST
                                I1 = PSAVE (NIBASE+M)
                                C1 = CCR (I1,R)
                                DO L = 1,NSUB
                                   Y (K+L,RS) = Y (K+L,RS)
     +                                          + C1 * W (L,I1)
                                END DO
                             END DO

                         END IF
                     END IF
  340             CONTINUE
C
C
C             ...next outer S contraction for present block.
C
C
  300          CONTINUE
C
C
C             ...next block (if any).
C
C
 3000       CONTINUE

           END IF

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
