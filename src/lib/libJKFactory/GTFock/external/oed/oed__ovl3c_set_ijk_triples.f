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
         SUBROUTINE  OED__OVL3C_SET_IJK_TRIPLES
     +
     +                    ( NPGTOI,NPGTOJ,NPGTOK,NPGTOIJ,
     +                      ATOMIC,EQUALIJ,
     +                      SWAPRS,
     +                      RNIJSQ,RNIKSQ,RNJKSQ,
     +                      PREFACT,
     +                      ALPHAI,ALPHAJ,ALPHAK,
     +                      SCREEN,
     +
     +                               EMPTY,
     +                               NIJ,NK,NIJK,
     +                               PRIMI,PRIMJ,PRIMK,
     +                               RHO )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__OVL3C_SET_IJK_TRIPLES
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation determines the ijk primitive exponent
C                triples which will actually be considered for
C                contraction. It is here where the internal prescreening
C                (if any) of 3-center overlap primitives must be
C                applied.
C
C                There are two ways the ij pairs corresponding to the
C                contraction shells I and J can be ordered: 1) i runs
C                fastest or 2) j runs fastest. The ordering to be chosen
C                depends on the order the contraction will be performed.
C                We have the following rule:
C
C                     If # primitives on I > # primitives on J
C                     (i.e. the SWAPRS = .true. case), the contraction
C                     will be performed first on the I primitives
C                     followed by the J primitives in order to save
C                     intermediate storage space for the 1st half
C                     transformation. Hence each J primitive will have
C                     associated a set of I primitives for the 1st
C                     half contraction and thus the order for the
C                     primitive pairs ij must be such that i runs
C                     fastest. If # primitives on I =< # primitives
C                     on J we have to let j run fastest.
C
C                The prescreening will be based on the fact that
C                for normalized cartesian integrals the primitive
C                3-center (sss) overlap integrals will be largest in
C                magnitude for all exponent combinations. The normalized
C                (sss) 3-center overlap integral value expression is:
C
C                 PREFACT * (ijk)**(3/4) * (i+j+k)**(-3/2) * rho (i,j,k)
C
C                where i,j,k are the i,j,k-th gaussian exponents, and
C                and rho (i,j,k) is the exponential prefactor given by
C                the expression:
C
C                         rho (i,j,k) = IJ * IK * JK
C
C                where
C                                      - (ij/(i+j+k)) * R(IJ)**2
C                             IJ = exp
C
C                                      - (ik/(i+j+k)) * R(IK)**2
C                             IK = exp
C
C                                      - (jk/(i+j+k)) * R(JK)**2
C                             JK = exp
C
C                and R(IJ)**2 is the square of the distance between
C                the i-th and j-th atomic centers I and J. Same for
C                R(IK)**2 and R(JK)**2.
C
C                The prescreening will be done in a 'K2' fashion for
C                the ij exponent pairs, followed by a 'K1' prescreening
C                for the remaining k exponent. The following sequence
C                of events is adopted:
C
C                    Find exponents i(min), j(min) and k(min), which
C                    gives the largest contribution to the (sss)
C                    integrals.
C
C                    Prescreen all ij exponent pairs using k(min).
C
C                    Prescreen all k exponents using i(min) and j(min).
C
C
C                An ij exponent pair/ k exponent will be rejected, if
C                the corresponding (sss) integral is less than a
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
C                                    for contraction shells x = I,J,K
C                    NPGTOIJ      =  # of primitive pairs for
C                                    contraction shell pair IJ
C                    EQUALIJ      =  indicates, if csh I and csh J are
C                                    considered to be equal
C                    SWAPRS       =  is .true. if the contraction order
C                                    of the primitive pairs IJ will
C                                    be performed in reverse order
C                                    JI for efficiency reasons
C                    RNxySQ       =  the square distance R(xy)**2
C                                    between centers xy = IJ,IK,JK
C                    PREFACT      =  the value of the overall prefactor
C                    ALPHAx       =  primitive exponents for centers
C                                    x=I,J,K
C                    SCREEN       =  is true, if screening will be
C                                    performed
C
C                  Output:
C
C                    EMPTY        =  is true, if no ij pairs and/or
C                                    k exponents were established
C                    NIJ          =  # of ij pairs after screening
C                    NK           =  # of k exponents after screening
C                    NIJK         =  # of ij pairs after screening
C                                    times # of k exponents after
C                                    screening
C                    PRIMx        =  i,j,k labels of primitives for
C                                    x=I,J,K
C                    RHO          =  NIJ * NK exponential prefactors
C                                    rho (i,j,k)
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
         LOGICAL   EQUALIJ
         LOGICAL   SCREEN
         LOGICAL   SWAPRS

         INTEGER   I,J,K,M,N
         INTEGER   KBASE
         INTEGER   KSAVE
         INTEGER   NIJ,NK,NIJK
         INTEGER   NPGTOI,NPGTOJ,NPGTOK,NPGTOIJ

         INTEGER   PRIMI (1:NPGTOIJ)
         INTEGER   PRIMJ (1:NPGTOIJ)
         INTEGER   PRIMK (1:NPGTOK)

         DOUBLE PRECISION  A,B,C,D,E,P,Q
         DOUBLE PRECISION  EXPI,EXPJ,EXPK
         DOUBLE PRECISION  EXPIMIN,EXPJMIN,EXPKMIN
         DOUBLE PRECISION  PREFACT
         DOUBLE PRECISION  RHOIJK
         DOUBLE PRECISION  RNIJSQ,RNIKSQ,RNJKSQ
         DOUBLE PRECISION  SSS
         DOUBLE PRECISION  ZP75,ONE,NEG1P5,TOL

         DOUBLE PRECISION  ALPHAI (1:NPGTOI)
         DOUBLE PRECISION  ALPHAJ (1:NPGTOJ)
         DOUBLE PRECISION  ALPHAK (1:NPGTOK)

         DOUBLE PRECISION  RHO (*)

         PARAMETER  (ZP75    =  0.75D0)
         PARAMETER  (ONE     =  1.D0)
         PARAMETER  (NEG1P5  = -1.5D0)
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

             IF (EQUALIJ) THEN
                 DO I = 1,NPGTOI
                 DO J = 1,I
                    NIJ = NIJ + 1
                    PRIMI (NIJ) = I
                    PRIMJ (NIJ) = J
                 END DO
                 END DO
             ELSE
                 IF (SWAPRS) THEN
                     DO J = 1,NPGTOJ
                     DO I = 1,NPGTOI
                        NIJ = NIJ + 1
                        PRIMI (NIJ) = I
                        PRIMJ (NIJ) = J
                     END DO
                     END DO
                 ELSE
                     DO I = 1,NPGTOI
                     DO J = 1,NPGTOJ
                        NIJ = NIJ + 1
                        PRIMI (NIJ) = I
                        PRIMJ (NIJ) = J
                     END DO
                     END DO
                 END IF
             END IF

             NK = NPGTOK
             DO K = 1,NPGTOK
                PRIMK (K) = K
             END DO
C
C
C             ...calculate the rho densities for no screening section
C                and place them in linear array RHO.
C
C
             NIJK = NIJ * NK

             IF (ATOMIC) THEN
                 DO M = 1,NIJK
                    RHO (M) = ONE
                 END DO
             ELSE
                 KBASE = 0
                 DO N = 1,NK
                    K = PRIMK (N)
                    EXPK = ALPHAK (K)
                    A = EXPK * RNIKSQ
                    B = EXPK * RNJKSQ
                    DO M = 1,NIJ
                       I = PRIMI (M)
                       J = PRIMJ (M)
                       EXPI = ALPHAI (I)
                       EXPJ = ALPHAJ (J)
                       C = EXPI * (EXPJ * RNIJSQ + A) + EXPJ * B
                       Q = EXPI + EXPJ + EXPK 
                       RHO (KBASE+M) = DEXP (-C/Q)
                    END DO
                    KBASE = KBASE + NIJ
                 END DO
             END IF

             RETURN

         END IF
C
C
C             ...the screening section: determine the minimum exponents
C                for all three centers I,J,K. Save also the primitive
C                K index corresponding to the minimum K exponent. This
C                index is needed to save the exponential prefactors
C                rho (i,j,kmin) evaluated during the ij exponent pairs
C                screening.
C
C
         EXPIMIN = ALPHAI (1)
         DO I = 2,NPGTOI
            EXPIMIN = DMIN1 (EXPIMIN,ALPHAI (I))
         END DO

         EXPJMIN = ALPHAJ (1)
         DO J = 2,NPGTOJ
            EXPJMIN = DMIN1 (EXPJMIN,ALPHAJ (J))
         END DO

         KSAVE = 1
         EXPKMIN = ALPHAK (1)
         DO K = 2,NPGTOK
            EXPK = ALPHAK (K)
            IF (EXPK.LT.EXPKMIN) THEN
                KSAVE = K
                EXPKMIN = EXPK
            END IF
         END DO
C
C
C             ...perform K2 primitive screening on I,J pairs.
C                Save the exponential prefactors rho (i,j,kmin)
C                into the first NIJ positions of array RHO. No saving
C                of exponential prefactors equal to 1 (i.e. the
C                atomic case) is done, as these are set later.
C
C
         NIJ = 0

         IF (EQUALIJ) THEN
             IF (ATOMIC) THEN
                 DO I = 1,NPGTOI
                    EXPI = ALPHAI (I)
                    DO J = 1,I
                       EXPJ = ALPHAJ (J)
                       P = EXPI * EXPJ * EXPKMIN
                       Q = EXPI + EXPJ + EXPKMIN

                       SSS = PREFACT * (P**ZP75) * (Q**NEG1P5)

                       IF (SSS.GE.TOL) THEN
                           NIJ = NIJ + 1
                           PRIMI (NIJ) = I
                           PRIMJ (NIJ) = J
                       ELSE
                           WRITE (*,*) ' Skip atom ovl3c I J = ',I,J
                       END IF

                    END DO
                 END DO
             ELSE
                 A = EXPKMIN * RNIKSQ
                 B = EXPKMIN * RNJKSQ

                 DO I = 1,NPGTOI
                    EXPI = ALPHAI (I)
                    DO J = 1,I
                       EXPJ = ALPHAJ (J)
                       C = EXPI * (EXPJ * RNIJSQ + A) + EXPJ * B
                       P = EXPI * EXPJ * EXPKMIN
                       Q = EXPI + EXPJ + EXPKMIN
                       RHOIJK = DEXP (-C/Q)

                       SSS = PREFACT * (P**ZP75) * (Q**NEG1P5) * RHOIJK

                       IF (SSS.GE.TOL) THEN
                           NIJ = NIJ + 1
                           RHO (NIJ) = RHOIJK
                           PRIMI (NIJ) = I
                           PRIMJ (NIJ) = J
                       ELSE
                           WRITE (*,*) ' Skip ovl3c I J = ',I,J
                       END IF

                    END DO
                 END DO
             END IF
         ELSE
             IF (ATOMIC) THEN
                 IF (SWAPRS) THEN
                     DO J = 1,NPGTOJ
                        EXPJ = ALPHAJ (J)
                        DO I = 1,NPGTOI
                           EXPI = ALPHAI (I)
                           P = EXPI * EXPJ * EXPKMIN
                           Q = EXPI + EXPJ + EXPKMIN

                           SSS = PREFACT * (P**ZP75) * (Q**NEG1P5)

                           IF (SSS.GE.TOL) THEN
                               NIJ = NIJ + 1
                               PRIMI (NIJ) = I
                               PRIMJ (NIJ) = J
                           ELSE
                               WRITE (*,*) ' Skip atom ovl3c I J = ',I,J
                           END IF

                        END DO
                     END DO
                 ELSE
                     DO I = 1,NPGTOI
                        EXPI = ALPHAI (I)
                        DO J = 1,NPGTOJ
                           EXPJ = ALPHAJ (J)
                           P = EXPI * EXPJ * EXPKMIN
                           Q = EXPI + EXPJ + EXPKMIN

                           SSS = PREFACT*(P**ZP75)*(Q**NEG1P5)*RHOIJK

                           IF (SSS.GE.TOL) THEN
                               NIJ = NIJ + 1
                               PRIMI (NIJ) = I
                               PRIMJ (NIJ) = J
                           ELSE
                               WRITE (*,*) ' Skip atom ovl3c I J = ',I,J
                           END IF

                        END DO
                     END DO
                 END IF
             ELSE
                 A = EXPKMIN * RNIKSQ
                 B = EXPKMIN * RNJKSQ

                 IF (SWAPRS) THEN
                     DO J = 1,NPGTOJ
                        EXPJ = ALPHAJ (J)
                        DO I = 1,NPGTOI
                           EXPI = ALPHAI (I)
                           C = EXPI * (EXPJ * RNIJSQ + A) + EXPJ * B
                           P = EXPI * EXPJ * EXPKMIN
                           Q = EXPI + EXPJ + EXPKMIN
                           RHOIJK = DEXP (-C/Q)

                           SSS = PREFACT*(P**ZP75)*(Q**NEG1P5)*RHOIJK

                           IF (SSS.GE.TOL) THEN
                               NIJ = NIJ + 1
                               RHO (NIJ) = RHOIJK
                               PRIMI (NIJ) = I
                               PRIMJ (NIJ) = J
                           ELSE
                               WRITE (*,*) ' Skip ovl3c I J = ',I,J
                           END IF

                        END DO
                     END DO
                 ELSE
                     DO I = 1,NPGTOI
                        EXPI = ALPHAI (I)
                        DO J = 1,NPGTOJ
                           EXPJ = ALPHAJ (J)
                           C = EXPI * (EXPJ * RNIJSQ + A) + EXPJ * B
                           P = EXPI * EXPJ * EXPKMIN
                           Q = EXPI + EXPJ + EXPKMIN
                           RHOIJK = DEXP (-C/Q)

                           SSS = PREFACT*(P**ZP75)*(Q**NEG1P5)*RHOIJK

                           IF (SSS.GE.TOL) THEN
                               NIJ = NIJ + 1
                               RHO (NIJ) = RHOIJK
                               PRIMI (NIJ) = I
                               PRIMJ (NIJ) = J
                           ELSE
                               WRITE (*,*) ' Skip ovl3c I J = ',I,J
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
C             ...perform K1 primitive screening on K exponents.
C                No saving of calculated exponential prefactors
C                rho (imin,jmin,k).
C
C
         NK = 0

         IF (ATOMIC) THEN
             A = EXPIMIN * EXPJMIN
             B = EXPIMIN + EXPJMIN

             DO K = 1,NPGTOK
                EXPK = ALPHAK (K)
                P = A * EXPK
                Q = B + EXPK

                SSS = PREFACT * (P**ZP75) * (Q**NEG1P5)

                IF (SSS.GE.TOL) THEN
                    NK = NK + 1
                    PRIMK (NK) = K
                ELSE
                    WRITE (*,*) ' Skipping atom ovl3c K = ',K
                END IF

             END DO
         ELSE
             A = EXPIMIN * EXPJMIN
             B = EXPIMIN + EXPJMIN
             C = A * RNIJSQ
             D = EXPIMIN * RNIKSQ + EXPJMIN * RNJKSQ

             DO K = 1,NPGTOK
                EXPK = ALPHAK (K)
                E = C + EXPK * D
                P = A * EXPK
                Q = B + EXPK
                RHOIJK = DEXP (-E/Q)

                SSS = PREFACT * (P**ZP75) * (Q**NEG1P5) * RHOIJK

                IF (SSS.GE.TOL) THEN
                    NK = NK + 1
                    PRIMK (NK) = K
                ELSE
                    WRITE (*,*) ' Skipping ovl3c K = ',K
                END IF

             END DO
         END IF

         IF (NK.EQ.0) THEN
             EMPTY = .TRUE.
             RETURN
         END IF
C
C
C             ...Complete the exponential prefactor array.
C                This has to be done from the last K exponent
C                column downwards, in order to be able to place
C                those exponential prefactors sitting in first
C                column (formally):
C
C                       RHO (*,1) -------> RHO (*,KSAVE)
C
C                The final exponential prefactor array is one-
C                dimensional, thus great care has to be exercised
C                in placing its elements into the right positions.
C                For the atomic case all elements of the exponential
C                prefactor array are equal to 1.
C
C
         NIJK = NIJ * NK

         IF (ATOMIC) THEN
             DO M = 1,NIJK
                RHO (M) = ONE
             END DO
         ELSE
             KBASE = NIJK - NIJ

             DO N = NK,1,-1
                K = PRIMK (N)

                IF (K.EQ.KSAVE .AND. N.GT.1) THEN
                    DO M = 1,NIJ
                       RHO (KBASE+M) = RHO (M)
                    END DO
                ELSE
                    EXPK = ALPHAK (K)
                    A = EXPK * RNIKSQ
                    B = EXPK * RNJKSQ

                    DO M = 1,NIJ
                       I = PRIMI (M)
                       J = PRIMJ (M)
                       EXPI = ALPHAI (I)
                       EXPJ = ALPHAJ (J)
                       C = EXPI * (EXPJ * RNIJSQ + A) + EXPJ * B
                       Q = EXPI + EXPJ + EXPK 
                       RHO (KBASE+M) = DEXP (-C/Q)
                    END DO
                END IF

                KBASE = KBASE - NIJ

             END DO
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
