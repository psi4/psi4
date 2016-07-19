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
         SUBROUTINE  ERD__CTR_2ND_HALF_NEW
     +
     +                    ( N,
     +                      NPMAX,NPMIN,
     +                      MKL,NTU,
     +                      NBLOCK,
     +                      NCT,NCU,
     +                      NPT,NPU,
     +                      CCT,CCU,
     +                      CCBEGT,CCBEGU,
     +                      CCENDT,CCENDU,
     +                      PRIMT,PRIMU,
     +                      EQUALTU,
     +                      SWAPTU,
     +                      PUSED,PSAVE,PPAIR,
     +                      X,W,
     +
     +                             Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__CTR_2ND_HALF_NEW
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation performs the second half contraction
C                step on the incomming half transformed integrals
C                in blocked form over invariant indices n and generates
C                new fully transformed contracted batch of integrals:
C
C                   y (n,tu) = sum  cct (t,k) * ccu (u,l) * x (n,kl)
C                               kl
C
C                where cct and ccu are the arrays containing the
C                contraction coefficients. The sum is over the k and
C                l primitives which are transmitted in arrays PRIMT
C                and PRIMU, respectively, and may constitute only
C                a subset of the full range of primitives.
C
C                The contraction is split into two quarter steps,
C                the order of which is determined by the # of k and
C                l primitives:
C
C                   a) w (n,l/k) = sum cct/u (t/u,k/l) * x (n,kl)
C                                  k/l
C
C                   b) y (n,tu)  = sum ccu/t (u/t,l/k) * w (n,l/k)
C                                  l/k
C
C                Size of the intermediate w (n,k) or w (n,l) array
C                has to be kept to a minimum and blocking over the
C                invariant indices n has to be performed such that
C                the intermediate w array does not get kicked out
C                from the cache lines after the first quarter
C                transformation.
C
C                In case of csh equality (EQUALTU = .true), we only
C                have to consider the lower triangle of the primitive
C                integrals, which, with the exception of the diagonals,
C                have to be used twice.
C
C                        --- SEGMENTED CONTRACTIONS ---
C
C                Segmented contractions are those defined to be
C                within a certain consecutive k- and l-range of
C                primitives. The segmented limits for each contraction
C                are sitting respectively in CCBEGT and CCBEGU (lowest
C                limit) and CCENDT and CCENDU (highest limit) and they
C                determine which of the actual k's and l's from the
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
C                                    sets k,l
C                    MKL          =  # of kl primitive products to
C                                    be transformed
C                    NTU          =  # of tu contractions to be done
C                    NBLOCK       =  blocking size for invariant
C                                    indices n
C                    NCT(U)       =  # of contractions for the k -> T
C                                    (l -> U) primitives
C                    NPT(U)       =  # of k(l) primitives
C                    CCT(U)       =  full set (including zeros) of
C                                    contraction coefficients for
C                                    T(U) contractions
C                    CCBEGT(U)    =  lowest nonzero primitive k(l)
C                                    index for T(U) contractions
C                    CCENDT(U)    =  highest nonzero primitive k(l)
C                                    index for T(U) contractions
C                    PRIMT(U)     =  primitive k(l) indices
C                    EQUALTU      =  is true, if only the lower
C                                    triangle of kl primitive indices
C                                    is present and consequently
C                                    only the lower triangle of tu
C                                    contractions needs to be evaluated
C                    SWAPTU       =  if this is true, the 1st quarter
C                                    transformation is over T followed
C                                    by the 2nd over U. If false, the
C                                    order is reversed: 1st over U then
C                                    2nd over T
C                    Pxxxx        =  intermediate storage arrays for
C                                    primitive labels to bundle
C                                    contraction steps in do loops
C                                    (xxxx = USED,SAVE,PAIR)
C                    X            =  array containing the half
C                                    transformed integrals
C                    W            =  intermediate storage array
C                                    containing 3/4 transformed
C                                    integrals
C                    Y            =  original batch of fully contracted
C                                    integrals
C
C                  Output:
C
C                    Y            =  updated batch of fully contracted
C                                    integrals
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

         LOGICAL     CNTRCT
         LOGICAL     EQUALTU
         LOGICAL     KRANGE,LRANGE
         LOGICAL     SWAPTU

         INTEGER     I,J,K,L,M,N,T,U
         INTEGER     K1,K2,K3,K4,K5,K6,K7,K8
         INTEGER     KL,TU
         INTEGER     KL1,KL2,KL3,KL4,KL5,KL6,KL7,KL8
         INTEGER     KNEXT,LNEXT
         INTEGER     L1,L2,L3,L4,L5,L6,L7,L8
         INTEGER     MKL,NTU
         INTEGER     NBLOCK,NSUB
         INTEGER     NCT,NCU
         INTEGER     NK,NL
         INTEGER     NKBASE,NLBASE
         INTEGER     NKLEFT,NLLEFT
         INTEGER     NKREST,NLREST
         INTEGER     NKSTEP,NLSTEP
         INTEGER     NPMAX,NPMIN
         INTEGER     NPT,NPU
         INTEGER     PMIN,PMAX

         INTEGER     CCBEGT (1:NCT)
         INTEGER     CCBEGU (1:NCU)
         INTEGER     CCENDT (1:NCT)
         INTEGER     CCENDU (1:NCU)

         INTEGER     PRIMT  (1:MKL)
         INTEGER     PRIMU  (1:MKL)

         INTEGER     PPAIR  (1:NPMAX)
         INTEGER     PSAVE  (1:NPMAX)
         INTEGER     PUSED  (1:NPMIN)

         DOUBLE PRECISION    C1,C2,C3,C4,C5,C6,C7,C8
         DOUBLE PRECISION    ZERO,ONE

         DOUBLE PRECISION    CCT  (1:NPT,1:NCT)
         DOUBLE PRECISION    CCU  (1:NPU,1:NCU)

         DOUBLE PRECISION    W (1:NBLOCK,1:NPMIN)

         DOUBLE PRECISION    X (1:N,1:MKL)
         DOUBLE PRECISION    Y (1:N,1:NTU)

         PARAMETER  (ZERO  =  0.D0)
         PARAMETER  (ONE   =  1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...check case T >= U or full TU range.
C
C
         IF (EQUALTU) THEN
C
C
C             ...the T >= U case. Here we always have equal # of
C                primitives K and L for both T and U. The primitives
C                K >= L are ordered such that L varies fastest.
C                Outer contraction is over U, inner over T.
C
C
           DO 1000 I = 0,N-1,NBLOCK
              NSUB = MIN0 (NBLOCK,N-I)

              TU = 0
              DO 100 U = 1,NCU
                 PMIN = CCBEGU (U)
                 PMAX = CCENDU (U)

                 DO 110 K = 1,NPT
                    PUSED (K) = 0
  110            CONTINUE

                 IF (PMIN.EQ.PMAX .AND. CCU (PMIN,U).EQ.ONE) THEN
C
C
C             ...we reach this point, if only one U contraction
C                coefficient equal to 1.0 is present.
C
C
                     DO 120 KL = 1,MKL
                        L = PRIMU (KL)
                        IF (L.EQ.PMIN) THEN
                            K = PRIMT (KL)
                            DO J = 1,NSUB
                               W (J,K) = X (I+J,KL)
                            END DO
                            PUSED (K) = 1
                        END IF
  120                CONTINUE
                 ELSE
C
C
C             ...the general U contraction case with many coefficients.
C                For each K determine all allowed L's and perform the
C                present U contraction for each K over all L's
C                simultaneously.
C
C
                     NL = 0
                     DO 130 KL = 1,MKL
                        K = PRIMT (KL)
                        L = PRIMU (KL)
                        LRANGE = L.GE.PMIN .AND. L.LE.PMAX

                        IF (LRANGE) THEN
                            NL = NL + 1
                            PSAVE (NL) = L
                            PPAIR (NL) = KL
                        END IF

                        IF (KL.EQ.MKL) THEN
                            KNEXT = 0
                        ELSE
                            KNEXT = PRIMT (KL+1)
                        END IF

                        CNTRCT = (KNEXT.NE.K) .AND. (NL.GT.0)

                        IF (CNTRCT) THEN
                            IF (NL.EQ.1) THEN
                                KL1 = PPAIR (1)
                                C1 = CCU (PSAVE (1),U)
                                DO J = 1,NSUB
                                   W (J,K) = C1 * X (I+J,KL1)
                                END DO
                            ELSE IF (NL.EQ.2) THEN
                                KL1 = PPAIR (1)
                                KL2 = PPAIR (2)
                                C1 = CCU (PSAVE (1),U)
                                C2 = CCU (PSAVE (2),U)
                                DO J = 1,NSUB
                                   W (J,K) =   C1 * X (I+J,KL1)
     +                                       + C2 * X (I+J,KL2)
                                END DO
                            ELSE IF (NL.EQ.3) THEN
                                KL1 = PPAIR (1)
                                KL2 = PPAIR (2)
                                KL3 = PPAIR (3)
                                C1 = CCU (PSAVE (1),U)
                                C2 = CCU (PSAVE (2),U)
                                C3 = CCU (PSAVE (3),U)
                                DO J = 1,NSUB
                                   W (J,K) =   C1 * X (I+J,KL1)
     +                                       + C2 * X (I+J,KL2)
     +                                       + C3 * X (I+J,KL3)
                                END DO
                            ELSE IF (NL.EQ.4) THEN
                                KL1 = PPAIR (1)
                                KL2 = PPAIR (2)
                                KL3 = PPAIR (3)
                                KL4 = PPAIR (4)
                                C1 = CCU (PSAVE (1),U)
                                C2 = CCU (PSAVE (2),U)
                                C3 = CCU (PSAVE (3),U)
                                C4 = CCU (PSAVE (4),U)
                                DO J = 1,NSUB
                                   W (J,K) =   C1 * X (I+J,KL1)
     +                                       + C2 * X (I+J,KL2)
     +                                       + C3 * X (I+J,KL3)
     +                                       + C4 * X (I+J,KL4)
                                END DO
                            ELSE IF (NL.EQ.5) THEN
                                KL1 = PPAIR (1)
                                KL2 = PPAIR (2)
                                KL3 = PPAIR (3)
                                KL4 = PPAIR (4)
                                KL5 = PPAIR (5)
                                C1 = CCU (PSAVE (1),U)
                                C2 = CCU (PSAVE (2),U)
                                C3 = CCU (PSAVE (3),U)
                                C4 = CCU (PSAVE (4),U)
                                C5 = CCU (PSAVE (5),U)
                                DO J = 1,NSUB
                                   W (J,K) =   C1 * X (I+J,KL1)
     +                                       + C2 * X (I+J,KL2)
     +                                       + C3 * X (I+J,KL3)
     +                                       + C4 * X (I+J,KL4)
     +                                       + C5 * X (I+J,KL5)
                                END DO
                            ELSE IF (NL.EQ.6) THEN
                                KL1 = PPAIR (1)
                                KL2 = PPAIR (2)
                                KL3 = PPAIR (3)
                                KL4 = PPAIR (4)
                                KL5 = PPAIR (5)
                                KL6 = PPAIR (6)
                                C1 = CCU (PSAVE (1),U)
                                C2 = CCU (PSAVE (2),U)
                                C3 = CCU (PSAVE (3),U)
                                C4 = CCU (PSAVE (4),U)
                                C5 = CCU (PSAVE (5),U)
                                C6 = CCU (PSAVE (6),U)
                                DO J = 1,NSUB
                                   W (J,K) =   C1 * X (I+J,KL1)
     +                                       + C2 * X (I+J,KL2)
     +                                       + C3 * X (I+J,KL3)
     +                                       + C4 * X (I+J,KL4)
     +                                       + C5 * X (I+J,KL5)
     +                                       + C6 * X (I+J,KL6)
                                END DO
                            ELSE IF (NL.EQ.7) THEN
                                KL1 = PPAIR (1)
                                KL2 = PPAIR (2)
                                KL3 = PPAIR (3)
                                KL4 = PPAIR (4)
                                KL5 = PPAIR (5)
                                KL6 = PPAIR (6)
                                KL7 = PPAIR (7)
                                C1 = CCU (PSAVE (1),U)
                                C2 = CCU (PSAVE (2),U)
                                C3 = CCU (PSAVE (3),U)
                                C4 = CCU (PSAVE (4),U)
                                C5 = CCU (PSAVE (5),U)
                                C6 = CCU (PSAVE (6),U)
                                C7 = CCU (PSAVE (7),U)
                                DO J = 1,NSUB
                                   W (J,K) =   C1 * X (I+J,KL1)
     +                                       + C2 * X (I+J,KL2)
     +                                       + C3 * X (I+J,KL3)
     +                                       + C4 * X (I+J,KL4)
     +                                       + C5 * X (I+J,KL5)
     +                                       + C6 * X (I+J,KL6)
     +                                       + C7 * X (I+J,KL7)
                                END DO
                            ELSE IF (NL.EQ.8) THEN
                                KL1 = PPAIR (1)
                                KL2 = PPAIR (2)
                                KL3 = PPAIR (3)
                                KL4 = PPAIR (4)
                                KL5 = PPAIR (5)
                                KL6 = PPAIR (6)
                                KL7 = PPAIR (7)
                                KL8 = PPAIR (8)
                                C1 = CCU (PSAVE (1),U)
                                C2 = CCU (PSAVE (2),U)
                                C3 = CCU (PSAVE (3),U)
                                C4 = CCU (PSAVE (4),U)
                                C5 = CCU (PSAVE (5),U)
                                C6 = CCU (PSAVE (6),U)
                                C7 = CCU (PSAVE (7),U)
                                C8 = CCU (PSAVE (8),U)
                                DO J = 1,NSUB
                                   W (J,K) =   C1 * X (I+J,KL1)
     +                                       + C2 * X (I+J,KL2)
     +                                       + C3 * X (I+J,KL3)
     +                                       + C4 * X (I+J,KL4)
     +                                       + C5 * X (I+J,KL5)
     +                                       + C6 * X (I+J,KL6)
     +                                       + C7 * X (I+J,KL7)
     +                                       + C8 * X (I+J,KL8)
                                END DO
                            ELSE
                                KL1 = PPAIR (1)
                                KL2 = PPAIR (2)
                                KL3 = PPAIR (3)
                                KL4 = PPAIR (4)
                                KL5 = PPAIR (5)
                                KL6 = PPAIR (6)
                                KL7 = PPAIR (7)
                                KL8 = PPAIR (8)
                                C1 = CCU (PSAVE (1),U)
                                C2 = CCU (PSAVE (2),U)
                                C3 = CCU (PSAVE (3),U)
                                C4 = CCU (PSAVE (4),U)
                                C5 = CCU (PSAVE (5),U)
                                C6 = CCU (PSAVE (6),U)
                                C7 = CCU (PSAVE (7),U)
                                C8 = CCU (PSAVE (8),U)
                                DO J = 1,NSUB
                                   W (J,K) =   C1 * X (I+J,KL1)
     +                                       + C2 * X (I+J,KL2)
     +                                       + C3 * X (I+J,KL3)
     +                                       + C4 * X (I+J,KL4)
     +                                       + C5 * X (I+J,KL5)
     +                                       + C6 * X (I+J,KL6)
     +                                       + C7 * X (I+J,KL7)
     +                                       + C8 * X (I+J,KL8)
                                END DO

                                NLBASE = 8
                                NLLEFT = NL - 8
                                NLSTEP = NLLEFT / 8
                                NLREST = MOD (NLLEFT,8)

                                DO M = 1,NLSTEP
                                   KL1 = PPAIR (NLBASE+1)
                                   KL2 = PPAIR (NLBASE+2)
                                   KL3 = PPAIR (NLBASE+3)
                                   KL4 = PPAIR (NLBASE+4)
                                   KL5 = PPAIR (NLBASE+5)
                                   KL6 = PPAIR (NLBASE+6)
                                   KL7 = PPAIR (NLBASE+7)
                                   KL8 = PPAIR (NLBASE+8)
                                   C1 = CCU (PSAVE (NLBASE+1),U)
                                   C2 = CCU (PSAVE (NLBASE+2),U)
                                   C3 = CCU (PSAVE (NLBASE+3),U)
                                   C4 = CCU (PSAVE (NLBASE+4),U)
                                   C5 = CCU (PSAVE (NLBASE+5),U)
                                   C6 = CCU (PSAVE (NLBASE+6),U)
                                   C7 = CCU (PSAVE (NLBASE+7),U)
                                   C8 = CCU (PSAVE (NLBASE+8),U)
                                   DO J = 1,NSUB
                                      W (J,K) = W (J,K)
     +                                             + C1 * X (I+J,KL1)
     +                                             + C2 * X (I+J,KL2)
     +                                             + C3 * X (I+J,KL3)
     +                                             + C4 * X (I+J,KL4)
     +                                             + C5 * X (I+J,KL5)
     +                                             + C6 * X (I+J,KL6)
     +                                             + C7 * X (I+J,KL7)
     +                                             + C8 * X (I+J,KL8)
                                   END DO
                                   NLBASE = NLBASE + 8
                                END DO

                                DO M = 1,NLREST
                                   KL1 = PPAIR (NLBASE+M)
                                   C1 = CCU (PSAVE (NLBASE+M),U)
                                   DO J = 1,NSUB
                                      W (J,K) = W (J,K)
     +                                             + C1 * X (I+J,KL1)
                                   END DO
                                END DO

                            END IF
                            PUSED (K) = 1
C
C
C             ...don't forget to add into W the contributions due
C                to the L indices due to triangularization. Observe
C                that there is no possibility of grouping the L's
C                together, unless we resort the K and L such that
C                K varies fastest. Also observe that only L's
C                will be addressed that are =< K, hence we do
C                not run into the danger of addressing W columns
C                prematurely which have L > K, which saves us
C                from checking addressing of W columns during the
C                above simultaneous L contractions for each K. 
C
C
                            IF (K.GE.PMIN .AND. K.LE.PMAX) THEN
                                C1 = CCU (K,U)
                                DO 135 M = 1,NL
                                   L = PSAVE (M)
                                   IF (L.NE.K) THEN
                                       KL1 = PPAIR (M)
                                       IF (PUSED (L).EQ.1) THEN
                                           DO J = 1,NSUB
                                              W (J,L) = W (J,L)
     +                                                + C1 * X (I+J,KL1)
                                           END DO
                                       ELSE
                                           DO J = 1,NSUB
                                              W (J,L) = C1 * X (I+J,KL1)
                                           END DO
                                           PUSED (L) = 1
                                       END IF
                                   END IF
  135                           CONTINUE
                            END IF

                            NL = 0
                        END IF
  130                CONTINUE
                 END IF
C
C
C             ...inner contraction over all T >= U.
C
C
                 DO 140 T = U,NCT
                    TU = TU + 1
                    PMIN = CCBEGT (T)
                    PMAX = CCENDT (T)

                    IF (PMIN.EQ.PMAX .AND. CCT (PMIN,T).EQ.ONE
     +                               .AND. PUSED (PMIN).EQ.1  ) THEN
C
C
C             ...we reach this point, if only one T contraction
C                coefficient equal to 1.0 is present.
C
C
                        DO J = 1,NSUB
                           Y (I+J,TU) = W (J,PMIN)
                        END DO
                    ELSE
C
C
C             ...the general T contraction case with many coefficients.
C                Group all relevant K's together and perform the
C                present T contraction over all K's simultaneously.
C
C
                        NK = 0
                        DO K = PMIN,PMAX
                           IF (PUSED (K).EQ.1) THEN
                               NK = NK + 1
                               PSAVE (NK) = K
                           END IF
                        END DO

                        IF (NK.EQ.0) THEN
                            DO J = 1,NSUB
                               Y (I+J,TU) = ZERO
                            END DO
                        ELSE IF (NK.EQ.1) THEN
                            K1 = PSAVE (1)
                            C1 = CCT (K1,T)
                            DO J = 1,NSUB
                               Y (I+J,TU) = C1 * W (J,K1)
                            END DO
                        ELSE IF (NK.EQ.2) THEN
                            K1 = PSAVE (1)
                            K2 = PSAVE (2)
                            C1 = CCT (K1,T)
                            C2 = CCT (K2,T)
                            DO J = 1,NSUB
                               Y (I+J,TU) =   C1 * W (J,K1)
     +                                      + C2 * W (J,K2)
                            END DO
                        ELSE IF (NK.EQ.3) THEN
                            K1 = PSAVE (1)
                            K2 = PSAVE (2)
                            K3 = PSAVE (3)
                            C1 = CCT (K1,T)
                            C2 = CCT (K2,T)
                            C3 = CCT (K3,T)
                            DO J = 1,NSUB
                               Y (I+J,TU) =   C1 * W (J,K1)
     +                                      + C2 * W (J,K2)
     +                                      + C3 * W (J,K3)
                            END DO
                        ELSE IF (NK.EQ.4) THEN
                            K1 = PSAVE (1)
                            K2 = PSAVE (2)
                            K3 = PSAVE (3)
                            K4 = PSAVE (4)
                            C1 = CCT (K1,T)
                            C2 = CCT (K2,T)
                            C3 = CCT (K3,T)
                            C4 = CCT (K4,T)
                            DO J = 1,NSUB
                               Y (I+J,TU) =   C1 * W (J,K1)
     +                                      + C2 * W (J,K2)
     +                                      + C3 * W (J,K3)
     +                                      + C4 * W (J,K4)
                            END DO
                        ELSE IF (NK.EQ.5) THEN
                            K1 = PSAVE (1)
                            K2 = PSAVE (2)
                            K3 = PSAVE (3)
                            K4 = PSAVE (4)
                            K5 = PSAVE (5)
                            C1 = CCT (K1,T)
                            C2 = CCT (K2,T)
                            C3 = CCT (K3,T)
                            C4 = CCT (K4,T)
                            C5 = CCT (K5,T)
                            DO J = 1,NSUB
                               Y (I+J,TU) =   C1 * W (J,K1)
     +                                      + C2 * W (J,K2)
     +                                      + C3 * W (J,K3)
     +                                      + C4 * W (J,K4)
     +                                      + C5 * W (J,K5)
                            END DO
                        ELSE IF (NK.EQ.6) THEN
                            K1 = PSAVE (1)
                            K2 = PSAVE (2)
                            K3 = PSAVE (3)
                            K4 = PSAVE (4)
                            K5 = PSAVE (5)
                            K6 = PSAVE (6)
                            C1 = CCT (K1,T)
                            C2 = CCT (K2,T)
                            C3 = CCT (K3,T)
                            C4 = CCT (K4,T)
                            C5 = CCT (K5,T)
                            C6 = CCT (K6,T)
                            DO J = 1,NSUB
                               Y (I+J,TU) =   C1 * W (J,K1)
     +                                      + C2 * W (J,K2)
     +                                      + C3 * W (J,K3)
     +                                      + C4 * W (J,K4)
     +                                      + C5 * W (J,K5)
     +                                      + C6 * W (J,K6)
                            END DO
                        ELSE IF (NK.EQ.7) THEN
                            K1 = PSAVE (1)
                            K2 = PSAVE (2)
                            K3 = PSAVE (3)
                            K4 = PSAVE (4)
                            K5 = PSAVE (5)
                            K6 = PSAVE (6)
                            K7 = PSAVE (7)
                            C1 = CCT (K1,T)
                            C2 = CCT (K2,T)
                            C3 = CCT (K3,T)
                            C4 = CCT (K4,T)
                            C5 = CCT (K5,T)
                            C6 = CCT (K6,T)
                            C7 = CCT (K7,T)
                            DO J = 1,NSUB
                               Y (I+J,TU) =   C1 * W (J,K1)
     +                                      + C2 * W (J,K2)
     +                                      + C3 * W (J,K3)
     +                                      + C4 * W (J,K4)
     +                                      + C5 * W (J,K5)
     +                                      + C6 * W (J,K6)
     +                                      + C7 * W (J,K7)
                            END DO
                        ELSE IF (NK.EQ.8) THEN
                            K1 = PSAVE (1)
                            K2 = PSAVE (2)
                            K3 = PSAVE (3)
                            K4 = PSAVE (4)
                            K5 = PSAVE (5)
                            K6 = PSAVE (6)
                            K7 = PSAVE (7)
                            K8 = PSAVE (8)
                            C1 = CCT (K1,T)
                            C2 = CCT (K2,T)
                            C3 = CCT (K3,T)
                            C4 = CCT (K4,T)
                            C5 = CCT (K5,T)
                            C6 = CCT (K6,T)
                            C7 = CCT (K7,T)
                            C8 = CCT (K8,T)
                            DO J = 1,NSUB
                               Y (I+J,TU) =   C1 * W (J,K1)
     +                                      + C2 * W (J,K2)
     +                                      + C3 * W (J,K3)
     +                                      + C4 * W (J,K4)
     +                                      + C5 * W (J,K5)
     +                                      + C6 * W (J,K6)
     +                                      + C7 * W (J,K7)
     +                                      + C8 * W (J,K8)
                            END DO
                        ELSE
                            K1 = PSAVE (1)
                            K2 = PSAVE (2)
                            K3 = PSAVE (3)
                            K4 = PSAVE (4)
                            K5 = PSAVE (5)
                            K6 = PSAVE (6)
                            K7 = PSAVE (7)
                            K8 = PSAVE (8)
                            C1 = CCT (K1,T)
                            C2 = CCT (K2,T)
                            C3 = CCT (K3,T)
                            C4 = CCT (K4,T)
                            C5 = CCT (K5,T)
                            C6 = CCT (K6,T)
                            C7 = CCT (K7,T)
                            C8 = CCT (K8,T)
                            DO J = 1,NSUB
                               Y (I+J,TU) =   C1 * W (J,K1)
     +                                      + C2 * W (J,K2)
     +                                      + C3 * W (J,K3)
     +                                      + C4 * W (J,K4)
     +                                      + C5 * W (J,K5)
     +                                      + C6 * W (J,K6)
     +                                      + C7 * W (J,K7)
     +                                      + C8 * W (J,K8)
                            END DO

                            NKBASE = 8
                            NKLEFT = NK - 8
                            NKSTEP = NKLEFT / 8
                            NKREST = MOD (NKLEFT,8)

                            DO M = 1,NKSTEP
                               K1 = PSAVE (NKBASE+1)
                               K2 = PSAVE (NKBASE+2)
                               K3 = PSAVE (NKBASE+3)
                               K4 = PSAVE (NKBASE+4)
                               K5 = PSAVE (NKBASE+5)
                               K6 = PSAVE (NKBASE+6)
                               K7 = PSAVE (NKBASE+7)
                               K8 = PSAVE (NKBASE+8)
                               C1 = CCT (K1,T)
                               C2 = CCT (K2,T)
                               C3 = CCT (K3,T)
                               C4 = CCT (K4,T)
                               C5 = CCT (K5,T)
                               C6 = CCT (K6,T)
                               C7 = CCT (K7,T)
                               C8 = CCT (K8,T)
                               DO J = 1,NSUB
                                  Y (I+J,TU) = Y (I+J,TU)
     +                                         + C1 * W (J,K1)
     +                                         + C2 * W (J,K2)
     +                                         + C3 * W (J,K3)
     +                                         + C4 * W (J,K4)
     +                                         + C5 * W (J,K5)
     +                                         + C6 * W (J,K6)
     +                                         + C7 * W (J,K7)
     +                                         + C8 * W (J,K8)
                               END DO
                               NKBASE = NKBASE + 8
                            END DO

                            DO M = 1,NKREST
                               K1 = PSAVE (NKBASE+M)
                               C1 = CCT (K1,T)
                               DO J = 1,NSUB
                                  Y (I+J,TU) = Y (I+J,TU)
     +                                         + C1 * W (J,K1)
                               END DO
                            END DO

                        END IF
                    END IF
  140            CONTINUE
C
C
C             ...next outer U contraction for present block.
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
C             ...the full TU case. Check the order of the two quarter
C                transformations and proceed accordingly.
C
C                The case: # of K primitives T > # of L primitives U
C                The primitives K and L are ordered such that K varies
C                fastest. Outer contraction is over T, inner over U.
C
C
          IF (SWAPTU) THEN

            DO 2000 I = 0,N-1,NBLOCK
               NSUB = MIN0 (NBLOCK,N-I)

               TU = 0
               DO 200 T = 1,NCT
                  PMIN = CCBEGT (T)
                  PMAX = CCENDT (T)

                  DO 210 L = 1,NPU
                     PUSED (L) = 0
  210             CONTINUE

                  IF (PMIN.EQ.PMAX .AND. CCT (PMIN,T).EQ.ONE) THEN
C
C
C             ...we reach this point, if only one T contraction
C                coefficient equal to 1.0 is present.
C
C
                      DO 220 KL = 1,MKL
                         K = PRIMT (KL)
                         IF (K.EQ.PMIN) THEN
                             L = PRIMU (KL)
                             DO J = 1,NSUB
                                W (J,L) = X (I+J,KL)
                             END DO
                             PUSED (L) = 1
                         END IF
  220                 CONTINUE
                  ELSE
C
C
C             ...the general T contraction case with many coefficients.
C                For each L determine all allowed K's and perform the
C                present T contraction for each L over all K's
C                simultaneously.
C
C
                      NK = 0
                      DO 230 KL = 1,MKL
                         K = PRIMT (KL)
                         L = PRIMU (KL)
                         KRANGE = K.GE.PMIN .AND. K.LE.PMAX

                         IF (KRANGE) THEN
                             NK = NK + 1
                             PSAVE (NK) = K
                             PPAIR (NK) = KL
                         END IF

                         IF (KL.EQ.MKL) THEN
                             LNEXT = 0
                         ELSE
                             LNEXT = PRIMU (KL+1)
                         END IF

                         CNTRCT = (LNEXT.NE.L) .AND. (NK.GT.0)

                         IF (CNTRCT) THEN
                             IF (NK.EQ.1) THEN
                                 KL1 = PPAIR (1)
                                 C1 = CCT (PSAVE (1),T)
                                 DO J = 1,NSUB
                                    W (J,L) = C1 * X (I+J,KL1)
                                 END DO
                             ELSE IF (NK.EQ.2) THEN
                                 KL1 = PPAIR (1)
                                 KL2 = PPAIR (2)
                                 C1 = CCT (PSAVE (1),T)
                                 C2 = CCT (PSAVE (2),T)
                                 DO J = 1,NSUB
                                    W (J,L) =   C1 * X (I+J,KL1)
     +                                        + C2 * X (I+J,KL2)
                                 END DO
                             ELSE IF (NK.EQ.3) THEN
                                 KL1 = PPAIR (1)
                                 KL2 = PPAIR (2)
                                 KL3 = PPAIR (3)
                                 C1 = CCT (PSAVE (1),T)
                                 C2 = CCT (PSAVE (2),T)
                                 C3 = CCT (PSAVE (3),T)
                                 DO J = 1,NSUB
                                    W (J,L) =   C1 * X (I+J,KL1)
     +                                        + C2 * X (I+J,KL2)
     +                                        + C3 * X (I+J,KL3)
                                 END DO
                             ELSE IF (NK.EQ.4) THEN
                                 KL1 = PPAIR (1)
                                 KL2 = PPAIR (2)
                                 KL3 = PPAIR (3)
                                 KL4 = PPAIR (4)
                                 C1 = CCT (PSAVE (1),T)
                                 C2 = CCT (PSAVE (2),T)
                                 C3 = CCT (PSAVE (3),T)
                                 C4 = CCT (PSAVE (4),T)
                                 DO J = 1,NSUB
                                    W (J,L) =   C1 * X (I+J,KL1)
     +                                        + C2 * X (I+J,KL2)
     +                                        + C3 * X (I+J,KL3)
     +                                        + C4 * X (I+J,KL4)
                                 END DO
                             ELSE IF (NK.EQ.5) THEN
                                 KL1 = PPAIR (1)
                                 KL2 = PPAIR (2)
                                 KL3 = PPAIR (3)
                                 KL4 = PPAIR (4)
                                 KL5 = PPAIR (5)
                                 C1 = CCT (PSAVE (1),T)
                                 C2 = CCT (PSAVE (2),T)
                                 C3 = CCT (PSAVE (3),T)
                                 C4 = CCT (PSAVE (4),T)
                                 C5 = CCT (PSAVE (5),T)
                                 DO J = 1,NSUB
                                    W (J,L) =   C1 * X (I+J,KL1)
     +                                        + C2 * X (I+J,KL2)
     +                                        + C3 * X (I+J,KL3)
     +                                        + C4 * X (I+J,KL4)
     +                                        + C5 * X (I+J,KL5)
                                 END DO
                             ELSE IF (NK.EQ.6) THEN
                                 KL1 = PPAIR (1)
                                 KL2 = PPAIR (2)
                                 KL3 = PPAIR (3)
                                 KL4 = PPAIR (4)
                                 KL5 = PPAIR (5)
                                 KL6 = PPAIR (6)
                                 C1 = CCT (PSAVE (1),T)
                                 C2 = CCT (PSAVE (2),T)
                                 C3 = CCT (PSAVE (3),T)
                                 C4 = CCT (PSAVE (4),T)
                                 C5 = CCT (PSAVE (5),T)
                                 C6 = CCT (PSAVE (6),T)
                                 DO J = 1,NSUB
                                    W (J,L) =   C1 * X (I+J,KL1)
     +                                        + C2 * X (I+J,KL2)
     +                                        + C3 * X (I+J,KL3)
     +                                        + C4 * X (I+J,KL4)
     +                                        + C5 * X (I+J,KL5)
     +                                        + C6 * X (I+J,KL6)
                                 END DO
                             ELSE IF (NK.EQ.7) THEN
                                 KL1 = PPAIR (1)
                                 KL2 = PPAIR (2)
                                 KL3 = PPAIR (3)
                                 KL4 = PPAIR (4)
                                 KL5 = PPAIR (5)
                                 KL6 = PPAIR (6)
                                 KL7 = PPAIR (7)
                                 C1 = CCT (PSAVE (1),T)
                                 C2 = CCT (PSAVE (2),T)
                                 C3 = CCT (PSAVE (3),T)
                                 C4 = CCT (PSAVE (4),T)
                                 C5 = CCT (PSAVE (5),T)
                                 C6 = CCT (PSAVE (6),T)
                                 C7 = CCT (PSAVE (7),T)
                                 DO J = 1,NSUB
                                    W (J,L) =   C1 * X (I+J,KL1)
     +                                        + C2 * X (I+J,KL2)
     +                                        + C3 * X (I+J,KL3)
     +                                        + C4 * X (I+J,KL4)
     +                                        + C5 * X (I+J,KL5)
     +                                        + C6 * X (I+J,KL6)
     +                                        + C7 * X (I+J,KL7)
                                 END DO
                             ELSE IF (NK.EQ.8) THEN
                                 KL1 = PPAIR (1)
                                 KL2 = PPAIR (2)
                                 KL3 = PPAIR (3)
                                 KL4 = PPAIR (4)
                                 KL5 = PPAIR (5)
                                 KL6 = PPAIR (6)
                                 KL7 = PPAIR (7)
                                 KL8 = PPAIR (8)
                                 C1 = CCT (PSAVE (1),T)
                                 C2 = CCT (PSAVE (2),T)
                                 C3 = CCT (PSAVE (3),T)
                                 C4 = CCT (PSAVE (4),T)
                                 C5 = CCT (PSAVE (5),T)
                                 C6 = CCT (PSAVE (6),T)
                                 C7 = CCT (PSAVE (7),T)
                                 C8 = CCT (PSAVE (8),T)
                                 DO J = 1,NSUB
                                    W (J,L) =   C1 * X (I+J,KL1)
     +                                        + C2 * X (I+J,KL2)
     +                                        + C3 * X (I+J,KL3)
     +                                        + C4 * X (I+J,KL4)
     +                                        + C5 * X (I+J,KL5)
     +                                        + C6 * X (I+J,KL6)
     +                                        + C7 * X (I+J,KL7)
     +                                        + C8 * X (I+J,KL8)
                                 END DO
                             ELSE
                                 KL1 = PPAIR (1)
                                 KL2 = PPAIR (2)
                                 KL3 = PPAIR (3)
                                 KL4 = PPAIR (4)
                                 KL5 = PPAIR (5)
                                 KL6 = PPAIR (6)
                                 KL7 = PPAIR (7)
                                 KL8 = PPAIR (8)
                                 C1 = CCT (PSAVE (1),T)
                                 C2 = CCT (PSAVE (2),T)
                                 C3 = CCT (PSAVE (3),T)
                                 C4 = CCT (PSAVE (4),T)
                                 C5 = CCT (PSAVE (5),T)
                                 C6 = CCT (PSAVE (6),T)
                                 C7 = CCT (PSAVE (7),T)
                                 C8 = CCT (PSAVE (8),T)
                                 DO J = 1,NSUB
                                    W (J,L) =   C1 * X (I+J,KL1)
     +                                        + C2 * X (I+J,KL2)
     +                                        + C3 * X (I+J,KL3)
     +                                        + C4 * X (I+J,KL4)
     +                                        + C5 * X (I+J,KL5)
     +                                        + C6 * X (I+J,KL6)
     +                                        + C7 * X (I+J,KL7)
     +                                        + C8 * X (I+J,KL8)
                                 END DO

                                 NKBASE = 8
                                 NKLEFT = NK - 8
                                 NKSTEP = NKLEFT / 8
                                 NKREST = MOD (NKLEFT,8)

                                 DO M = 1,NKSTEP
                                    KL1 = PPAIR (NKBASE+1)
                                    KL2 = PPAIR (NKBASE+2)
                                    KL3 = PPAIR (NKBASE+3)
                                    KL4 = PPAIR (NKBASE+4)
                                    KL5 = PPAIR (NKBASE+5)
                                    KL6 = PPAIR (NKBASE+6)
                                    KL7 = PPAIR (NKBASE+7)
                                    KL8 = PPAIR (NKBASE+8)
                                    C1 = CCT (PSAVE (NKBASE+1),T)
                                    C2 = CCT (PSAVE (NKBASE+2),T)
                                    C3 = CCT (PSAVE (NKBASE+3),T)
                                    C4 = CCT (PSAVE (NKBASE+4),T)
                                    C5 = CCT (PSAVE (NKBASE+5),T)
                                    C6 = CCT (PSAVE (NKBASE+6),T)
                                    C7 = CCT (PSAVE (NKBASE+7),T)
                                    C8 = CCT (PSAVE (NKBASE+8),T)
                                    DO J = 1,NSUB
                                       W (J,L) = W (J,L)
     +                                              + C1 * X (I+J,KL1)
     +                                              + C2 * X (I+J,KL2)
     +                                              + C3 * X (I+J,KL3)
     +                                              + C4 * X (I+J,KL4)
     +                                              + C5 * X (I+J,KL5)
     +                                              + C6 * X (I+J,KL6)
     +                                              + C7 * X (I+J,KL7)
     +                                              + C8 * X (I+J,KL8)
                                    END DO
                                    NKBASE = NKBASE + 8
                                 END DO

                                 DO M = 1,NKREST
                                    KL1 = PPAIR (NKBASE+M)
                                    C1 = CCT (PSAVE (NKBASE+M),T)
                                    DO J = 1,NSUB
                                       W (J,L) = W (J,L)
     +                                              + C1 * X (I+J,KL1)
                                    END DO
                                 END DO

                             END IF
                             PUSED (L) = 1
                             NK = 0
                         END IF

  230                 CONTINUE
                  END IF
C
C
C             ...inner contraction over all U.
C
C
                  DO 240 U = 1,NCU
                     TU = TU + 1
                     PMIN = CCBEGU (U)
                     PMAX = CCENDU (U)

                     IF (PMIN.EQ.PMAX .AND. CCU (PMIN,U).EQ.ONE
     +                                .AND. PUSED (PMIN).EQ.1  ) THEN
C
C
C             ...we reach this point, if only one U contraction
C                coefficient equal to 1.0 is present.
C
C
                         DO J = 1,NSUB
                            Y (I+J,TU) = W (J,PMIN)
                         END DO
                     ELSE
C
C
C             ...the general U contraction case with many coefficients.
C                Group all relevant L's together and perform the
C                present U contraction over all L's simultaneously.
C
C
                         NL = 0
                         DO L = PMIN,PMAX
                            IF (PUSED (L).EQ.1) THEN
                                NL = NL + 1
                                PSAVE (NL) = L
                            END IF
                         END DO

                         IF (NL.EQ.0) THEN
                             DO J = 1,NSUB
                                Y (I+J,TU) = ZERO
                             END DO
                         ELSE IF (NL.EQ.1) THEN
                             L1 = PSAVE (1)
                             C1 = CCU (L1,U)
                             DO J = 1,NSUB
                                Y (I+J,TU) = C1 * W (J,L1)
                             END DO
                         ELSE IF (NL.EQ.2) THEN
                             L1 = PSAVE (1)
                             L2 = PSAVE (2)
                             C1 = CCU (L1,U)
                             C2 = CCU (L2,U)
                             DO J = 1,NSUB
                                Y (I+J,TU) =   C1 * W (J,L1)
     +                                       + C2 * W (J,L2)
                             END DO
                         ELSE IF (NL.EQ.3) THEN
                             L1 = PSAVE (1)
                             L2 = PSAVE (2)
                             L3 = PSAVE (3)
                             C1 = CCU (L1,U)
                             C2 = CCU (L2,U)
                             C3 = CCU (L3,U)
                             DO J = 1,NSUB
                                Y (I+J,TU) =   C1 * W (J,L1)
     +                                       + C2 * W (J,L2)
     +                                       + C3 * W (J,L3)
                             END DO
                         ELSE IF (NL.EQ.4) THEN
                             L1 = PSAVE (1)
                             L2 = PSAVE (2)
                             L3 = PSAVE (3)
                             L4 = PSAVE (4)
                             C1 = CCU (L1,U)
                             C2 = CCU (L2,U)
                             C3 = CCU (L3,U)
                             C4 = CCU (L4,U)
                             DO J = 1,NSUB
                                Y (I+J,TU) =   C1 * W (J,L1)
     +                                       + C2 * W (J,L2)
     +                                       + C3 * W (J,L3)
     +                                       + C4 * W (J,L4)
                             END DO
                         ELSE IF (NL.EQ.5) THEN
                             L1 = PSAVE (1)
                             L2 = PSAVE (2)
                             L3 = PSAVE (3)
                             L4 = PSAVE (4)
                             L5 = PSAVE (5)
                             C1 = CCU (L1,U)
                             C2 = CCU (L2,U)
                             C3 = CCU (L3,U)
                             C4 = CCU (L4,U)
                             C5 = CCU (L5,U)
                             DO J = 1,NSUB
                                Y (I+J,TU) =   C1 * W (J,L1)
     +                                       + C2 * W (J,L2)
     +                                       + C3 * W (J,L3)
     +                                       + C4 * W (J,L4)
     +                                       + C5 * W (J,L5)
                             END DO
                         ELSE IF (NL.EQ.6) THEN
                             L1 = PSAVE (1)
                             L2 = PSAVE (2)
                             L3 = PSAVE (3)
                             L4 = PSAVE (4)
                             L5 = PSAVE (5)
                             L6 = PSAVE (6)
                             C1 = CCU (L1,U)
                             C2 = CCU (L2,U)
                             C3 = CCU (L3,U)
                             C4 = CCU (L4,U)
                             C5 = CCU (L5,U)
                             C6 = CCU (L6,U)
                             DO J = 1,NSUB
                                Y (I+J,TU) =   C1 * W (J,L1)
     +                                       + C2 * W (J,L2)
     +                                       + C3 * W (J,L3)
     +                                       + C4 * W (J,L4)
     +                                       + C5 * W (J,L5)
     +                                       + C6 * W (J,L6)
                             END DO
                         ELSE IF (NL.EQ.7) THEN
                             L1 = PSAVE (1)
                             L2 = PSAVE (2)
                             L3 = PSAVE (3)
                             L4 = PSAVE (4)
                             L5 = PSAVE (5)
                             L6 = PSAVE (6)
                             L7 = PSAVE (7)
                             C1 = CCU (L1,U)
                             C2 = CCU (L2,U)
                             C3 = CCU (L3,U)
                             C4 = CCU (L4,U)
                             C5 = CCU (L5,U)
                             C6 = CCU (L6,U)
                             C7 = CCU (L7,U)
                             DO J = 1,NSUB
                                Y (I+J,TU) =   C1 * W (J,L1)
     +                                       + C2 * W (J,L2)
     +                                       + C3 * W (J,L3)
     +                                       + C4 * W (J,L4)
     +                                       + C5 * W (J,L5)
     +                                       + C6 * W (J,L6)
     +                                       + C7 * W (J,L7)
                             END DO
                         ELSE IF (NL.EQ.8) THEN
                             L1 = PSAVE (1)
                             L2 = PSAVE (2)
                             L3 = PSAVE (3)
                             L4 = PSAVE (4)
                             L5 = PSAVE (5)
                             L6 = PSAVE (6)
                             L7 = PSAVE (7)
                             L8 = PSAVE (8)
                             C1 = CCU (L1,U)
                             C2 = CCU (L2,U)
                             C3 = CCU (L3,U)
                             C4 = CCU (L4,U)
                             C5 = CCU (L5,U)
                             C6 = CCU (L6,U)
                             C7 = CCU (L7,U)
                             C8 = CCU (L8,U)
                             DO J = 1,NSUB
                                Y (I+J,TU) =   C1 * W (J,L1)
     +                                       + C2 * W (J,L2)
     +                                       + C3 * W (J,L3)
     +                                       + C4 * W (J,L4)
     +                                       + C5 * W (J,L5)
     +                                       + C6 * W (J,L6)
     +                                       + C7 * W (J,L7)
     +                                       + C8 * W (J,L8)
                             END DO
                         ELSE
                             L1 = PSAVE (1)
                             L2 = PSAVE (2)
                             L3 = PSAVE (3)
                             L4 = PSAVE (4)
                             L5 = PSAVE (5)
                             L6 = PSAVE (6)
                             L7 = PSAVE (7)
                             L8 = PSAVE (8)
                             C1 = CCU (L1,U)
                             C2 = CCU (L2,U)
                             C3 = CCU (L3,U)
                             C4 = CCU (L4,U)
                             C5 = CCU (L5,U)
                             C6 = CCU (L6,U)
                             C7 = CCU (L7,U)
                             C8 = CCU (L8,U)
                             DO J = 1,NSUB
                                Y (I+J,TU) =   C1 * W (J,L1)
     +                                       + C2 * W (J,L2)
     +                                       + C3 * W (J,L3)
     +                                       + C4 * W (J,L4)
     +                                       + C5 * W (J,L5)
     +                                       + C6 * W (J,L6)
     +                                       + C7 * W (J,L7)
     +                                       + C8 * W (J,L8)
                             END DO

                             NLBASE = 8
                             NLLEFT = NL - 8
                             NLSTEP = NLLEFT / 8
                             NLREST = MOD (NLLEFT,8)

                             DO M = 1,NLSTEP
                                L1 = PSAVE (NLBASE+1)
                                L2 = PSAVE (NLBASE+2)
                                L3 = PSAVE (NLBASE+3)
                                L4 = PSAVE (NLBASE+4)
                                L5 = PSAVE (NLBASE+5)
                                L6 = PSAVE (NLBASE+6)
                                L7 = PSAVE (NLBASE+7)
                                L8 = PSAVE (NLBASE+8)
                                C1 = CCU (L1,U)
                                C2 = CCU (L2,U)
                                C3 = CCU (L3,U)
                                C4 = CCU (L4,U)
                                C5 = CCU (L5,U)
                                C6 = CCU (L6,U)
                                C7 = CCU (L7,U)
                                C8 = CCU (L8,U)
                                DO J = 1,NSUB
                                   Y (I+J,TU) = Y (I+J,TU)
     +                                          + C1 * W (J,L1)
     +                                          + C2 * W (J,L2)
     +                                          + C3 * W (J,L3)
     +                                          + C4 * W (J,L4)
     +                                          + C5 * W (J,L5)
     +                                          + C6 * W (J,L6)
     +                                          + C7 * W (J,L7)
     +                                          + C8 * W (J,L8)
                                END DO
                                NLBASE = NLBASE + 8
                             END DO

                             DO M = 1,NLREST
                                L1 = PSAVE (NLBASE+M)
                                C1 = CCU (L1,U)
                                DO J = 1,NSUB
                                   Y (I+J,TU) = Y (I+J,TU)
     +                                          + C1 * W (J,L1)
                                END DO
                             END DO

                         END IF
                     END IF
  240             CONTINUE
C
C
C             ...next outer T contraction for present block.
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
C             ...the case: # of K primitives T =< # of L primitives U.
C                The primitives K and L are ordered such that L varies
C                fastest. Outer contraction is over U, inner over T.
C
C
            DO 3000 I = 0,N-1,NBLOCK
               NSUB = MIN0 (NBLOCK,N-I)

               TU = 0
               DO 300 U = 1,NCU
                  PMIN = CCBEGU (U)
                  PMAX = CCENDU (U)

                  DO 310 K = 1,NPT
                     PUSED (K) = 0
  310             CONTINUE

                  IF (PMIN.EQ.PMAX .AND. CCU (PMIN,U).EQ.ONE) THEN
C
C
C             ...we reach this point, if only one U contraction
C                coefficient equal to 1.0 is present.
C
C
                      DO 320 KL = 1,MKL
                         L = PRIMU (KL)
                         IF (L.EQ.PMIN) THEN
                             K = PRIMT (KL)
                             DO J = 1,NSUB
                                W (J,K) = X (I+J,KL)
                             END DO
                             PUSED (K) = 1
                         END IF
  320                 CONTINUE
                  ELSE
C
C
C             ...the general U contraction case with many coefficients.
C                For each K determine all allowed L's and perform the
C                present U contraction for each K over all L's
C                simultaneously.
C
C
                      NL = 0
                      DO 330 KL = 1,MKL
                         K = PRIMT (KL)
                         L = PRIMU (KL)
                         LRANGE = L.GE.PMIN .AND. L.LE.PMAX

                         IF (LRANGE) THEN
                             NL = NL + 1
                             PSAVE (NL) = L
                             PPAIR (NL) = KL
                         END IF

                         IF (KL.EQ.MKL) THEN
                             KNEXT = 0
                         ELSE
                             KNEXT = PRIMT (KL+1)
                         END IF

                         CNTRCT = (KNEXT.NE.K) .AND. (NL.GT.0)

                         IF (CNTRCT) THEN
                             IF (NL.EQ.1) THEN
                                 KL1 = PPAIR (1)
                                 C1 = CCU (PSAVE (1),U)
                                 DO J = 1,NSUB
                                    W (J,K) = C1 * X (I+J,KL1)
                                 END DO
                             ELSE IF (NL.EQ.2) THEN
                                 KL1 = PPAIR (1)
                                 KL2 = PPAIR (2)
                                 C1 = CCU (PSAVE (1),U)
                                 C2 = CCU (PSAVE (2),U)
                                 DO J = 1,NSUB
                                    W (J,K) =   C1 * X (I+J,KL1)
     +                                        + C2 * X (I+J,KL2)
                                 END DO
                             ELSE IF (NL.EQ.3) THEN
                                 KL1 = PPAIR (1)
                                 KL2 = PPAIR (2)
                                 KL3 = PPAIR (3)
                                 C1 = CCU (PSAVE (1),U)
                                 C2 = CCU (PSAVE (2),U)
                                 C3 = CCU (PSAVE (3),U)
                                 DO J = 1,NSUB
                                    W (J,K) =   C1 * X (I+J,KL1)
     +                                        + C2 * X (I+J,KL2)
     +                                        + C3 * X (I+J,KL3)
                                 END DO
                             ELSE IF (NL.EQ.4) THEN
                                 KL1 = PPAIR (1)
                                 KL2 = PPAIR (2)
                                 KL3 = PPAIR (3)
                                 KL4 = PPAIR (4)
                                 C1 = CCU (PSAVE (1),U)
                                 C2 = CCU (PSAVE (2),U)
                                 C3 = CCU (PSAVE (3),U)
                                 C4 = CCU (PSAVE (4),U)
                                 DO J = 1,NSUB
                                    W (J,K) =   C1 * X (I+J,KL1)
     +                                        + C2 * X (I+J,KL2)
     +                                        + C3 * X (I+J,KL3)
     +                                        + C4 * X (I+J,KL4)
                                 END DO
                             ELSE IF (NL.EQ.5) THEN
                                 KL1 = PPAIR (1)
                                 KL2 = PPAIR (2)
                                 KL3 = PPAIR (3)
                                 KL4 = PPAIR (4)
                                 KL5 = PPAIR (5)
                                 C1 = CCU (PSAVE (1),U)
                                 C2 = CCU (PSAVE (2),U)
                                 C3 = CCU (PSAVE (3),U)
                                 C4 = CCU (PSAVE (4),U)
                                 C5 = CCU (PSAVE (5),U)
                                 DO J = 1,NSUB
                                    W (J,K) =   C1 * X (I+J,KL1)
     +                                        + C2 * X (I+J,KL2)
     +                                        + C3 * X (I+J,KL3)
     +                                        + C4 * X (I+J,KL4)
     +                                        + C5 * X (I+J,KL5)
                                 END DO
                             ELSE IF (NL.EQ.6) THEN
                                 KL1 = PPAIR (1)
                                 KL2 = PPAIR (2)
                                 KL3 = PPAIR (3)
                                 KL4 = PPAIR (4)
                                 KL5 = PPAIR (5)
                                 KL6 = PPAIR (6)
                                 C1 = CCU (PSAVE (1),U)
                                 C2 = CCU (PSAVE (2),U)
                                 C3 = CCU (PSAVE (3),U)
                                 C4 = CCU (PSAVE (4),U)
                                 C5 = CCU (PSAVE (5),U)
                                 C6 = CCU (PSAVE (6),U)
                                 DO J = 1,NSUB
                                    W (J,K) =   C1 * X (I+J,KL1)
     +                                        + C2 * X (I+J,KL2)
     +                                        + C3 * X (I+J,KL3)
     +                                        + C4 * X (I+J,KL4)
     +                                        + C5 * X (I+J,KL5)
     +                                        + C6 * X (I+J,KL6)
                                 END DO
                             ELSE IF (NL.EQ.7) THEN
                                 KL1 = PPAIR (1)
                                 KL2 = PPAIR (2)
                                 KL3 = PPAIR (3)
                                 KL4 = PPAIR (4)
                                 KL5 = PPAIR (5)
                                 KL6 = PPAIR (6)
                                 KL7 = PPAIR (7)
                                 C1 = CCU (PSAVE (1),U)
                                 C2 = CCU (PSAVE (2),U)
                                 C3 = CCU (PSAVE (3),U)
                                 C4 = CCU (PSAVE (4),U)
                                 C5 = CCU (PSAVE (5),U)
                                 C6 = CCU (PSAVE (6),U)
                                 C7 = CCU (PSAVE (7),U)
                                 DO J = 1,NSUB
                                    W (J,K) =   C1 * X (I+J,KL1)
     +                                        + C2 * X (I+J,KL2)
     +                                        + C3 * X (I+J,KL3)
     +                                        + C4 * X (I+J,KL4)
     +                                        + C5 * X (I+J,KL5)
     +                                        + C6 * X (I+J,KL6)
     +                                        + C7 * X (I+J,KL7)
                                 END DO
                             ELSE IF (NL.EQ.8) THEN
                                 KL1 = PPAIR (1)
                                 KL2 = PPAIR (2)
                                 KL3 = PPAIR (3)
                                 KL4 = PPAIR (4)
                                 KL5 = PPAIR (5)
                                 KL6 = PPAIR (6)
                                 KL7 = PPAIR (7)
                                 KL8 = PPAIR (8)
                                 C1 = CCU (PSAVE (1),U)
                                 C2 = CCU (PSAVE (2),U)
                                 C3 = CCU (PSAVE (3),U)
                                 C4 = CCU (PSAVE (4),U)
                                 C5 = CCU (PSAVE (5),U)
                                 C6 = CCU (PSAVE (6),U)
                                 C7 = CCU (PSAVE (7),U)
                                 C8 = CCU (PSAVE (8),U)
                                 DO J = 1,NSUB
                                    W (J,K) =   C1 * X (I+J,KL1)
     +                                        + C2 * X (I+J,KL2)
     +                                        + C3 * X (I+J,KL3)
     +                                        + C4 * X (I+J,KL4)
     +                                        + C5 * X (I+J,KL5)
     +                                        + C6 * X (I+J,KL6)
     +                                        + C7 * X (I+J,KL7)
     +                                        + C8 * X (I+J,KL8)
                                 END DO
                             ELSE
                                 KL1 = PPAIR (1)
                                 KL2 = PPAIR (2)
                                 KL3 = PPAIR (3)
                                 KL4 = PPAIR (4)
                                 KL5 = PPAIR (5)
                                 KL6 = PPAIR (6)
                                 KL7 = PPAIR (7)
                                 KL8 = PPAIR (8)
                                 C1 = CCU (PSAVE (1),U)
                                 C2 = CCU (PSAVE (2),U)
                                 C3 = CCU (PSAVE (3),U)
                                 C4 = CCU (PSAVE (4),U)
                                 C5 = CCU (PSAVE (5),U)
                                 C6 = CCU (PSAVE (6),U)
                                 C7 = CCU (PSAVE (7),U)
                                 C8 = CCU (PSAVE (8),U)
                                 DO J = 1,NSUB
                                    W (J,K) =   C1 * X (I+J,KL1)
     +                                        + C2 * X (I+J,KL2)
     +                                        + C3 * X (I+J,KL3)
     +                                        + C4 * X (I+J,KL4)
     +                                        + C5 * X (I+J,KL5)
     +                                        + C6 * X (I+J,KL6)
     +                                        + C7 * X (I+J,KL7)
     +                                        + C8 * X (I+J,KL8)
                                 END DO

                                 NLBASE = 8
                                 NLLEFT = NL - 8
                                 NLSTEP = NLLEFT / 8
                                 NLREST = MOD (NLLEFT,8)

                                 DO M = 1,NLSTEP
                                    KL1 = PPAIR (NLBASE+1)
                                    KL2 = PPAIR (NLBASE+2)
                                    KL3 = PPAIR (NLBASE+3)
                                    KL4 = PPAIR (NLBASE+4)
                                    KL5 = PPAIR (NLBASE+5)
                                    KL6 = PPAIR (NLBASE+6)
                                    KL7 = PPAIR (NLBASE+7)
                                    KL8 = PPAIR (NLBASE+8)
                                    C1 = CCU (PSAVE (NLBASE+1),U)
                                    C2 = CCU (PSAVE (NLBASE+2),U)
                                    C3 = CCU (PSAVE (NLBASE+3),U)
                                    C4 = CCU (PSAVE (NLBASE+4),U)
                                    C5 = CCU (PSAVE (NLBASE+5),U)
                                    C6 = CCU (PSAVE (NLBASE+6),U)
                                    C7 = CCU (PSAVE (NLBASE+7),U)
                                    C8 = CCU (PSAVE (NLBASE+8),U)
                                    DO J = 1,NSUB
                                       W (J,K) = W (J,K)
     +                                              + C1 * X (I+J,KL1)
     +                                              + C2 * X (I+J,KL2)
     +                                              + C3 * X (I+J,KL3)
     +                                              + C4 * X (I+J,KL4)
     +                                              + C5 * X (I+J,KL5)
     +                                              + C6 * X (I+J,KL6)
     +                                              + C7 * X (I+J,KL7)
     +                                              + C8 * X (I+J,KL8)
                                    END DO
                                    NLBASE = NLBASE + 8
                                 END DO

                                 DO M = 1,NLREST
                                    KL1 = PPAIR (NLBASE+M)
                                    C1 = CCU (PSAVE (NLBASE+M),U)
                                    DO J = 1,NSUB
                                       W (J,K) = W (J,K)
     +                                              + C1 * X (I+J,KL1)
                                    END DO
                                 END DO

                             END IF
                             PUSED (K) = 1
                             NL = 0
                         END IF

  330                 CONTINUE
                  END IF
C
C
C             ...inner contraction over all T.
C
C
                  DO 340 T = 1,NCT
                     TU = TU + 1
                     PMIN = CCBEGT (T)
                     PMAX = CCENDT (T)

                     IF (PMIN.EQ.PMAX .AND. CCT (PMIN,T).EQ.ONE
     +                                .AND. PUSED (PMIN).EQ.1  ) THEN
C
C
C             ...we reach this point, if only one T contraction
C                coefficient equal to 1.0 is present.
C
C
                         DO J = 1,NSUB
                            Y (I+J,TU) = W (J,PMIN)
                         END DO
                     ELSE
C
C
C             ...the general T contraction case with many coefficients.
C                Group all relevant K's together and perform the
C                present T contraction over all K's simultaneously.
C
C
                         NK = 0
                         DO K = PMIN,PMAX
                            IF (PUSED (K).EQ.1) THEN
                                NK = NK + 1
                                PSAVE (NK) = K
                            END IF
                         END DO

                         IF (NK.EQ.0) THEN
                             DO J = 1,NSUB
                                Y (I+J,TU) = ZERO
                             END DO
                         ELSE IF (NK.EQ.1) THEN
                             K1 = PSAVE (1)
                             C1 = CCT (K1,T)
                             DO J = 1,NSUB
                                Y (I+J,TU) = C1 * W (J,K1)
                             END DO
                         ELSE IF (NK.EQ.2) THEN
                             K1 = PSAVE (1)
                             K2 = PSAVE (2)
                             C1 = CCT (K1,T)
                             C2 = CCT (K2,T)
                             DO J = 1,NSUB
                                Y (I+J,TU) =   C1 * W (J,K1)
     +                                       + C2 * W (J,K2)
                             END DO
                         ELSE IF (NK.EQ.3) THEN
                             K1 = PSAVE (1)
                             K2 = PSAVE (2)
                             K3 = PSAVE (3)
                             C1 = CCT (K1,T)
                             C2 = CCT (K2,T)
                             C3 = CCT (K3,T)
                             DO J = 1,NSUB
                                Y (I+J,TU) =   C1 * W (J,K1)
     +                                       + C2 * W (J,K2)
     +                                       + C3 * W (J,K3)
                             END DO
                         ELSE IF (NK.EQ.4) THEN
                             K1 = PSAVE (1)
                             K2 = PSAVE (2)
                             K3 = PSAVE (3)
                             K4 = PSAVE (4)
                             C1 = CCT (K1,T)
                             C2 = CCT (K2,T)
                             C3 = CCT (K3,T)
                             C4 = CCT (K4,T)
                             DO J = 1,NSUB
                                Y (I+J,TU) =   C1 * W (J,K1)
     +                                       + C2 * W (J,K2)
     +                                       + C3 * W (J,K3)
     +                                       + C4 * W (J,K4)
                             END DO
                         ELSE IF (NK.EQ.5) THEN
                             K1 = PSAVE (1)
                             K2 = PSAVE (2)
                             K3 = PSAVE (3)
                             K4 = PSAVE (4)
                             K5 = PSAVE (5)
                             C1 = CCT (K1,T)
                             C2 = CCT (K2,T)
                             C3 = CCT (K3,T)
                             C4 = CCT (K4,T)
                             C5 = CCT (K5,T)
                             DO J = 1,NSUB
                                Y (I+J,TU) =   C1 * W (J,K1)
     +                                       + C2 * W (J,K2)
     +                                       + C3 * W (J,K3)
     +                                       + C4 * W (J,K4)
     +                                       + C5 * W (J,K5)
                             END DO
                         ELSE IF (NK.EQ.6) THEN
                             K1 = PSAVE (1)
                             K2 = PSAVE (2)
                             K3 = PSAVE (3)
                             K4 = PSAVE (4)
                             K5 = PSAVE (5)
                             K6 = PSAVE (6)
                             C1 = CCT (K1,T)
                             C2 = CCT (K2,T)
                             C3 = CCT (K3,T)
                             C4 = CCT (K4,T)
                             C5 = CCT (K5,T)
                             C6 = CCT (K6,T)
                             DO J = 1,NSUB
                                Y (I+J,TU) =   C1 * W (J,K1)
     +                                       + C2 * W (J,K2)
     +                                       + C3 * W (J,K3)
     +                                       + C4 * W (J,K4)
     +                                       + C5 * W (J,K5)
     +                                       + C6 * W (J,K6)
                             END DO
                         ELSE IF (NK.EQ.7) THEN
                             K1 = PSAVE (1)
                             K2 = PSAVE (2)
                             K3 = PSAVE (3)
                             K4 = PSAVE (4)
                             K5 = PSAVE (5)
                             K6 = PSAVE (6)
                             K7 = PSAVE (7)
                             C1 = CCT (K1,T)
                             C2 = CCT (K2,T)
                             C3 = CCT (K3,T)
                             C4 = CCT (K4,T)
                             C5 = CCT (K5,T)
                             C6 = CCT (K6,T)
                             C7 = CCT (K7,T)
                             DO J = 1,NSUB
                                Y (I+J,TU) =   C1 * W (J,K1)
     +                                       + C2 * W (J,K2)
     +                                       + C3 * W (J,K3)
     +                                       + C4 * W (J,K4)
     +                                       + C5 * W (J,K5)
     +                                       + C6 * W (J,K6)
     +                                       + C7 * W (J,K7)
                             END DO
                         ELSE IF (NK.EQ.8) THEN
                             K1 = PSAVE (1)
                             K2 = PSAVE (2)
                             K3 = PSAVE (3)
                             K4 = PSAVE (4)
                             K5 = PSAVE (5)
                             K6 = PSAVE (6)
                             K7 = PSAVE (7)
                             K8 = PSAVE (8)
                             C1 = CCT (K1,T)
                             C2 = CCT (K2,T)
                             C3 = CCT (K3,T)
                             C4 = CCT (K4,T)
                             C5 = CCT (K5,T)
                             C6 = CCT (K6,T)
                             C7 = CCT (K7,T)
                             C8 = CCT (K8,T)
                             DO J = 1,NSUB
                                Y (I+J,TU) =   C1 * W (J,K1)
     +                                       + C2 * W (J,K2)
     +                                       + C3 * W (J,K3)
     +                                       + C4 * W (J,K4)
     +                                       + C5 * W (J,K5)
     +                                       + C6 * W (J,K6)
     +                                       + C7 * W (J,K7)
     +                                       + C8 * W (J,K8)
                             END DO
                         ELSE
                             K1 = PSAVE (1)
                             K2 = PSAVE (2)
                             K3 = PSAVE (3)
                             K4 = PSAVE (4)
                             K5 = PSAVE (5)
                             K6 = PSAVE (6)
                             K7 = PSAVE (7)
                             K8 = PSAVE (8)
                             C1 = CCT (K1,T)
                             C2 = CCT (K2,T)
                             C3 = CCT (K3,T)
                             C4 = CCT (K4,T)
                             C5 = CCT (K5,T)
                             C6 = CCT (K6,T)
                             C7 = CCT (K7,T)
                             C8 = CCT (K8,T)
                             DO J = 1,NSUB
                                Y (I+J,TU) =   C1 * W (J,K1)
     +                                       + C2 * W (J,K2)
     +                                       + C3 * W (J,K3)
     +                                       + C4 * W (J,K4)
     +                                       + C5 * W (J,K5)
     +                                       + C6 * W (J,K6)
     +                                       + C7 * W (J,K7)
     +                                       + C8 * W (J,K8)
                             END DO

                             NKBASE = 8
                             NKLEFT = NK - 8
                             NKSTEP = NKLEFT / 8
                             NKREST = MOD (NKLEFT,8)

                             DO M = 1,NKSTEP
                                K1 = PSAVE (NKBASE+1)
                                K2 = PSAVE (NKBASE+2)
                                K3 = PSAVE (NKBASE+3)
                                K4 = PSAVE (NKBASE+4)
                                K5 = PSAVE (NKBASE+5)
                                K6 = PSAVE (NKBASE+6)
                                K7 = PSAVE (NKBASE+7)
                                K8 = PSAVE (NKBASE+8)
                                C1 = CCT (K1,T)
                                C2 = CCT (K2,T)
                                C3 = CCT (K3,T)
                                C4 = CCT (K4,T)
                                C5 = CCT (K5,T)
                                C6 = CCT (K6,T)
                                C7 = CCT (K7,T)
                                C8 = CCT (K8,T)
                                DO J = 1,NSUB
                                   Y (I+J,TU) = Y (I+J,TU)
     +                                          + C1 * W (J,K1)
     +                                          + C2 * W (J,K2)
     +                                          + C3 * W (J,K3)
     +                                          + C4 * W (J,K4)
     +                                          + C5 * W (J,K5)
     +                                          + C6 * W (J,K6)
     +                                          + C7 * W (J,K7)
     +                                          + C8 * W (J,K8)
                                END DO
                                NKBASE = NKBASE + 8
                             END DO

                             DO M = 1,NKREST
                                K1 = PSAVE (NKBASE+M)
                                C1 = CCT (K1,T)
                                DO J = 1,NSUB
                                   Y (I+J,TU) = Y (I+J,TU)
     +                                          + C1 * W (J,K1)
                                END DO
                             END DO

                         END IF
                     END IF
  340             CONTINUE
C
C
C             ...next outer U contraction for present block.
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
