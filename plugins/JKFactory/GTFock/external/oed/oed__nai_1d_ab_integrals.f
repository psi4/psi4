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
         SUBROUTINE  OED__NAI_1D_AB_INTEGRALS
     +
     +                    ( SHELLP,SHELLX,
     +                      SHELLA,SHELLB,
     +                      NGEXCEN,
     +                      WTS,
     +                      R1,R2,
     +                      AB,
     +                      WTAKE,
     +                      INTSCR,
     +
     +                               INT1D )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_1D_AB_INTEGRALS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation calculates a full table of 1D AB
C                nuclear attraction integrals using the Rys vertical
C                recurrence scheme and the Rys horizontal transfer
C                scheme explained below.
C
C                Scheme:
C
C                    i) VRR => generate the (P,0) integrals
C                   ii) HRR => generate the (A,B) integrals
C
C                If WTAKE is set true, the Rys weight is being
C                multiplied to the 1D AB integrals to reduce overall
C                FLOP count. Note, that the Rys weight factor needs
C                to be introduced only two times for the starting
C                1D AB integrals for the VRR recurrence scheme,
C                namely to the (0,0) and (1,0) elements. The weight
C                factor is then automatically propagated through the
C                vertical and horizontal transfer equations (see below).
C
C
C                The recurrence schemes VRR and HRR is due to Rys,
C                Dupuis and King, J. Comp. Chem. 4, p.154-157 (1983).
C
C
C                 i) VRR Scheme:
C
C
C                   INT1D (0,0) = 1.D0    (* WEIGHT if WTAKE = .true.)
C                   INT1D (1,0) = R1      (* WEIGHT if WTAKE = .true.)
C
C                   For I = 1,...,SHELLP-1
C                       INT1D (I+1,0) = I * R2 * INT1D (I-1,0)
C                                         + R1 * INT1D (I,0)
C
C
C                ii) HRR Scheme:
C
C
C                   For J = 1,...,SHELLX  (SHELLX = Min (SHELLA,SHELLB)
C                   For I = SHELLP-J,...,0
C                       INT1D (I,J) =  +/- AB * INT1D (I,J-1)
C                                             + INT1D (I+1,J-1)
C
C
C                The 1D AB integrals are calculated for all nuclear
C                centers, all exponent pairs and all quadrature points
C                simultaneously and placed into a 3-dimensional array.
C
C                An important feature of this routine is that it can
C                be called with shell magnitudes A < B. If this case
C                happens, the VRR and HRR have to be processed
C                differently. If shell A < shell B, the HRR must be
C                defined such that P,0 -> B,A due to efficiency and
C                numerical stability reasons. This in turn needs a
C                redefinition of the R1 coefficients, which now must
C                be such that they are defined with the P shell
C                accumulated on center B instead of A. The R1
C                coefficients are (example x coordinate):
C
C                  i) P accumulated on center A:
C
C                        R1x = (Px - Ax) + term independent of Ax,Bx
C
C                 ii) P accumulated on center B:
C
C                        R1x = (Px - Bx) + term independent of Ax,Bx
C
C                When entering the present routine however, all we have
C                are the R1 values based on case i). If case ii)
C                applies, we have to form new R1 via:
C
C                          R1 -> R1 + (A - B) = R1 + AB
C
C                using the cartesian coordinate differences AB = A - B
C                between centers A and B. The constant AB is thus added
C                to all R1 values presently transmitted. The R2
C                coefficients are independent of Ax and Bx and thus
C                need not be modified.
C
C
C                  Input:
C
C                    SHELLP      =  maximum shell sum A+B
C                    SHELLX      =  minimum shell type between A and B
C                    SHELLx      =  shell type for csh x = A and B
C                    NGEXCEN     =  # of roots times # of primitive
C                                   exponent pairs times # of nuclear
C                                   attraction centers
C                    WTS         =  all quadrature Rys weights
C                    R1          =  cartesian coordinate dependent
C                                   VRR expansion coefficients based
C                                   on center A
C                    R2          =  the coordinate independent VRR
C                                   R2-coefficients
C                    AB          =  cartesian coordinate differences
C                                   A - B between sites A and B
C                    WTAKE       =  if true, the Rys weights will be
C                                   build into the 1D AB integrals;
C                                   if false, the Rys weights will not
C                                   be considered
C                    INTSCR      =  scratch array to perform the VRR
C                                   and HRR schemes
C
C
C                  Output:
C
C                    INT1D       =  all 1D AB integrals with reduced
C                                   PB -> AB dimensions
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT  NONE

         LOGICAL   HRRBA
         LOGICAL   WTAKE

         INTEGER   CASE1D
         INTEGER   I,J,N
         INTEGER   I1,I2,J1,J2,IP1,JM1
         INTEGER   NGEXCEN
         INTEGER   SHELLA,SHELLB
         INTEGER   SHELLP,SHELLX

         DOUBLE PRECISION  AB
         DOUBLE PRECISION  F
         DOUBLE PRECISION  WEIGHT
         DOUBLE PRECISION  ZERO,ONE

         DOUBLE PRECISION  R1  (1:NGEXCEN)
         DOUBLE PRECISION  R2  (1:NGEXCEN)
         DOUBLE PRECISION  WTS (1:NGEXCEN)

         DOUBLE PRECISION  INT1D  (1:NGEXCEN,0:SHELLA,0:SHELLB)
         DOUBLE PRECISION  INTSCR (1:NGEXCEN,0:SHELLP,0:SHELLX)

         PARAMETER  (ONE   =  1.D0)
         PARAMETER  (ZERO  =  0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...jump according to the 4 different cases that can arise:
C
C                  A-shell = s- or higher angular momentum
C                  B-shell = s- or higher angular momentum
C
C                each leading to specific simplifications.
C
C
         CASE1D = 2 * MIN (1,SHELLA) + MIN (1,SHELLB) + 1

         GOTO (1,2,3,4) CASE1D
C
C
C             ...the case A = s-shell and B = s-shell (no HRR!).
C                i) VRR => Evaluate: I = 0
C                                    J = 0
C
C
    1    IF (WTAKE) THEN
             DO N = 1,NGEXCEN
                INT1D (N,0,0) = WTS (N)
             END DO
         ELSE
             DO N = 1,NGEXCEN
                INT1D (N,0,0) = ONE
            END DO
         END IF

         RETURN
C
C
C             ...the cases A = s-shell and B >= p-shell (no HRR!).
C                i) VRR => Evaluate: I = 0
C                                    J = 0,SHELLB
C
C
    2    IF (AB.NE.ZERO) THEN
             DO N = 1,NGEXCEN
                R1 (N) = R1 (N) + AB
             END DO
         END IF

         IF (WTAKE) THEN
             DO N = 1,NGEXCEN
                WEIGHT = WTS (N)
                INT1D (N,0,0) = WEIGHT
                INT1D (N,0,1) = R1 (N) * WEIGHT
             END DO
         ELSE
             DO N = 1,NGEXCEN
                INT1D (N,0,0) = ONE
                INT1D (N,0,1) = R1 (N)
             END DO
         END IF

         IF (SHELLB.GT.1) THEN
             F = ONE
             DO J = 2,SHELLB
                J1 = J - 1
                J2 = J - 2
                DO N = 1,NGEXCEN
                   INT1D (N,0,J) = F * R2 (N) * INT1D (N,0,J2)
     +                               + R1 (N) * INT1D (N,0,J1)
                END DO
                F = F + ONE
             END DO
         END IF

         RETURN
C
C
C             ...the cases A >= p-shell and B = s-shell (no HRR!).
C                i) VRR => Evaluate: I = 0,SHELLA
C
C
    3    IF (WTAKE) THEN
             DO N = 1,NGEXCEN
                WEIGHT = WTS (N)
                INT1D (N,0,0) = WEIGHT
                INT1D (N,1,0) = R1 (N) * WEIGHT
             END DO
         ELSE
             DO N = 1,NGEXCEN
                INT1D (N,0,0) = ONE
                INT1D (N,1,0) = R1 (N)
             END DO
         END IF

         IF (SHELLA.GT.1) THEN
             F = ONE
             DO I = 2,SHELLA
                I1 = I - 1
                I2 = I - 2
                DO N = 1,NGEXCEN
                   INT1D (N,I,0) = F * R2 (N) * INT1D (N,I2,0)
     +                               + R1 (N) * INT1D (N,I1,0)
                END DO
                F = F + ONE
             END DO
         END IF

         RETURN
C
C
C             ...the cases A >= p-shell and B >= p-shell.
C                i) VRR => Evaluate: I = 0,SHELLP
C                                    J = 0
C
C
    4    HRRBA = SHELLA .LT. SHELLB

         IF (HRRBA .AND. AB.NE.ZERO) THEN
             DO N = 1,NGEXCEN
                R1 (N) = R1 (N) + AB
             END DO
         END IF

         IF (WTAKE) THEN
             DO N = 1,NGEXCEN
                WEIGHT = WTS (N)
                INTSCR (N,0,0) = WEIGHT
                INTSCR (N,1,0) = R1 (N) * WEIGHT
             END DO
         ELSE
             DO N = 1,NGEXCEN
                INTSCR (N,0,0) = ONE
                INTSCR (N,1,0) = R1 (N)
             END DO
         END IF

         F = ONE
         DO I = 2,SHELLP
            I1 = I - 1
            I2 = I - 2
            DO N = 1,NGEXCEN
               INTSCR (N,I,0) = F * R2 (N) * INTSCR (N,I2,0)
     +                            + R1 (N) * INTSCR (N,I1,0)
            END DO
            F = F + ONE
         END DO
C
C
C             ...ii) HRR => Evaluate: I = 0,SHELLP-SHELLA/SHELLB
C                                     J = 0,SHELLB/SHELLA
C                           and
C
C               iii) Copy INTSCR to INT1D: I = 0,SHELLA
C                                          J = 0,SHELLB
C
C
         IF (.NOT.HRRBA) THEN

             IF (AB.EQ.ZERO) THEN
                 DO J = 1,SHELLB
                    JM1 = J - 1
                    DO I = 0,SHELLP-J
                       IP1 = I + 1
                       DO N = 1,NGEXCEN
                          INTSCR (N,I,J) = INTSCR (N,IP1,JM1)
                       END DO
                    END DO
                 END DO
             ELSE
                 DO J = 1,SHELLB
                    JM1 = J - 1
                    DO I = 0,SHELLP-J
                       IP1 = I + 1
                       DO N = 1,NGEXCEN
                          INTSCR (N,I,J) = INTSCR (N,IP1,JM1)
     +                              + AB * INTSCR (N,I,JM1)
                       END DO
                    END DO
                 END DO
             END IF

             DO J = 0,SHELLB
             DO I = 0,SHELLA
             DO N = 1,NGEXCEN
                INT1D (N,I,J) = INTSCR (N,I,J)
             END DO
             END DO
             END DO

         ELSE

             IF (AB.EQ.ZERO) THEN
                 DO J = 1,SHELLA
                    JM1 = J - 1
                    DO I = 0,SHELLP-J
                       IP1 = I + 1
                       DO N = 1,NGEXCEN
                          INTSCR (N,I,J) = INTSCR (N,IP1,JM1)
                       END DO
                    END DO
                 END DO
             ELSE
                 DO J = 1,SHELLA
                    JM1 = J - 1
                    DO I = 0,SHELLP-J
                       IP1 = I + 1
                       DO N = 1,NGEXCEN
                          INTSCR (N,I,J) = INTSCR (N,IP1,JM1)
     +                              - AB * INTSCR (N,I,JM1)
                       END DO
                    END DO
                 END DO
             END IF

             DO J = 0,SHELLB
             DO I = 0,SHELLA
             DO N = 1,NGEXCEN
                INT1D (N,I,J) = INTSCR (N,J,I)
             END DO
             END DO
             END DO

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
