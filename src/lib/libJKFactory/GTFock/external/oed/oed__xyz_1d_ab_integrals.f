C  Copyright (c) 1997-1999, 2003 Massachusetts Institute of Technology
C 
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
C  USA
         SUBROUTINE  OED__XYZ_1D_AB_INTEGRALS
     +
     +                    ( SHELLP,SHELLX,
     +                      SHELLA,SHELLB,
     +                      NEXP,
     +                      PA,
     +                      PB,                            ! Watson Added
     +                      PINVHF,
     +                      AB,
     +                      INTSCR,
     +
     +                               INT1D )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__OVL_1D_AB_INTEGRALS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation calculates a full table of 1D AB
C                overlap integrals using the Rys vertical recurrence
C                scheme and the Rys horizontal transfer scheme
C                explained below.
C
C                The recurrence schemes VRR and HRR is due to Rys,
C                Dupuis and King, J. Comp. Chem. 4, p.154-157 (1983).
C
C
C                 i) VRR Scheme:
C
C
C                   INT1D (0,0) = 1.D0
C                   INT1D (1,0) = PA
C
C                   For I = 1,...,SHELLP-1
C                       INT1D (I+1,0) = I * (1/2P) * INT1D (I-1,0)
C                                             + PA * INT1D (I,0)
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
C                The 1D AB integrals are calculated for all exponent
C                pairs and placed into a 3-dimensional array.
C
C                An important feature of this routine is that it can
C                be called with shell magnitudes A < B. If this case
C                happens, the VRR and HRR have to be processed
C                differently. If shell A < shell B, the HRR must be
C                defined such that P,0 -> B,A due to efficiency and
C                numerical stability reasons. This in turn needs a
C                redefinition of the PA coefficients, which now must
C                be such that they are defined with the P shell
C                accumulated on center B instead of A. The PA
C                coefficients are (example x coordinate):
C
C                  i) P accumulated on center A:
C
C                          PA = (Px - Ax)
C
C                 ii) P accumulated on center B:
C
C                          PA = (Px - Bx)
C
C                When entering the present routine however, all we have
C                are the PA values based on case i). If case ii)
C                applies, we have to form new PA via:
C
C                          PA -> PA + (A - B) = PA + AB
C
C                using the cartesian coordinate differences AB = A - B
C                between centers A and B. The constant AB is thus added
C                to all PA values presently transmitted.
C
C
C                  Input:
C
C                    SHELLP      =  maximum shell sum A+B
C                    SHELLX      =  minimum shell type between A and B
C                    SHELLx      =  shell type for csh x = A and B
C                    NEXP        =  # of exponent pairs
C                    PA          =  current NEXP coordinate differences
C                                   P-A between centers P and A
C                    PINVHF      =  current NEXP values of 1/(2*P),
C                                   where P are the exponent sums
C                                   for contraction shells A and B
C                    AB          =  cartesian coordinate differences
C                                   A - B between sites A and B
C                    INTSCR      =  scratch array to perform the VRR
C                                   and HRR schemes
C
C
C                  Output:
C
C                    INT1D       =  all 1D AB integrals with reduced
C                                   PX -> AB dimensions
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

         INTEGER   CASE1D
         INTEGER   I,J,N
         INTEGER   I1,I2,J1,J2
         INTEGER   NEXP
         INTEGER   SHELLA,SHELLB
         INTEGER   SHELLP,SHELLX

         DOUBLE PRECISION  AB
         DOUBLE PRECISION  F,P1
         DOUBLE PRECISION  ZERO,ONE,TWO

         DOUBLE PRECISION  PA     (1:NEXP)
         DOUBLE PRECISION  PB     (1:NEXP)                 ! Watson Added
         DOUBLE PRECISION  PINVHF (1:NEXP)

         DOUBLE PRECISION  INT1D  (1:NEXP,0:SHELLA,0:SHELLB)
         DOUBLE PRECISION  INTSCR (1:NEXP,0:SHELLP,0:SHELLX)

         PARAMETER  (ZERO  =  0.D0)
         PARAMETER  (ONE   =  1.D0)
         PARAMETER  (TWO   =  2.D0)
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
    1    DO N = 1,NEXP
            INT1D (N,0,0) = ONE
         END DO

         RETURN
C
C
C             ...the cases A = s-shell and B >= p-shell (no HRR!).
C                i) VRR => Evaluate: I = 0
C                                    J = 0,SHELLB
C
C
    2    DO N = 1,NEXP
            INT1D (N,0,0) = ONE
            INT1D (N,0,1) = PB (N)
         END DO

         F = ONE
         DO J = 2,SHELLB
            J1 = J - 1
            J2 = J - 2
            DO N = 1,NEXP
               P1 = F * PINVHF (N)
               INT1D (N,0,J) =        P1 * INT1D (N,0,J2)
     +                          + PB (N) * INT1D (N,0,J1)
            END DO
            F = F + ONE
         END DO

         RETURN
C
C
C             ...the cases A >= p-shell and B = s-shell (no HRR!).
C                i) VRR => Evaluate: I = 0,SHELLA
C                                    J = 0
C
C
    3    DO N = 1,NEXP
            INT1D (N,0,0) = ONE
            INT1D (N,1,0) = PA (N)
         END DO

         F = ONE
         DO I = 2,SHELLA
            I1 = I - 1
            I2 = I - 2
            DO N = 1,NEXP
               P1 = F * PINVHF (N)
               INT1D (N,I,0) =        P1 * INT1D (N,I2,0)
     +                          + PA (N) * INT1D (N,I1,0)
            END DO
            F = F + ONE
         END DO

         RETURN
C
C
C             ...the cases A >= p-shell and B >= p-shell.
C
C                  i) VRR => Evaluate: I = 0,SHELLP
C                                      J = 0
C
C                 ii) HRR => Evaluate: I = 0,SHELLP-SHELLA/SHELLB
C                                      J = 0,SHELLB/SHELLA
C
C                iii) Copy INTSCR to INT1D: I = 0,SHELLA
C                                           J = 0,SHELLB
C
C
C
    4    HRRBA = SHELLA .LT. SHELLB

         IF (HRRBA .AND. AB.NE.ZERO) THEN
             DO N = 1,NEXP
                PA (N) = PA (N) + AB
             END DO
         END IF

         DO N = 1,NEXP
            INTSCR (N,0,0) = ONE
            INTSCR (N,1,0) = PA (N)
         END DO

         F = ONE
         DO I = 2,SHELLP
            I1 = I - 1
            I2 = I - 2
            DO N = 1,NEXP
               P1 = F * PINVHF (N)
               INTSCR (N,I,0) =        P1 * INTSCR (N,I2,0)
     +                           + PA (N) * INTSCR (N,I1,0)
            END DO
            F = F + ONE
         END DO

         IF (.NOT.HRRBA) THEN

             IF (AB.EQ.ZERO) THEN
                 DO J = 1,SHELLB
                    J1 = J - 1
                    DO I = 0,SHELLP-J
                       I1 = I + 1
                       DO N = 1,NEXP
                          INTSCR (N,I,J) = INTSCR (N,I1,J1)
                       END DO
                    END DO
                 END DO
             ELSE
                 DO J = 1,SHELLB
                    J1 = J - 1
                    DO I = 0,SHELLP-J
                       I1 = I + 1
                       DO N = 1,NEXP
                          INTSCR (N,I,J) =        INTSCR (N,I1,J1)
     +                                     + AB * INTSCR (N,I,J1)
                       END DO
                    END DO
                 END DO
             END IF

             DO J = 0,SHELLB
             DO I = 0,SHELLA
             DO N = 1,NEXP
                INT1D (N,I,J) = INTSCR (N,I,J)
             END DO
             END DO
             END DO

         ELSE

             IF (AB.EQ.ZERO) THEN
                 DO J = 1,SHELLA
                    J1 = J - 1
                    DO I = 0,SHELLP-J
                       I1 = I + 1
                       DO N = 1,NEXP
                          INTSCR (N,I,J) = INTSCR (N,I1,J1)
                       END DO
                    END DO
                 END DO
             ELSE
                 DO J = 1,SHELLA
                    J1 = J - 1
                    DO I = 0,SHELLP-J
                       I1 = I + 1
                       DO N = 1,NEXP
                          INTSCR (N,I,J) =        INTSCR (N,I1,J1)
     +                                     - AB * INTSCR (N,I,J1)
                       END DO
                    END DO
                 END DO
             END IF

             DO J = 0,SHELLB
             DO I = 0,SHELLA
             DO N = 1,NEXP
                INT1D (N,I,J) = INTSCR (N,J,I)
             END DO
             END DO
             END DO

         END IF

         IF (HRRBA .AND. AB.NE.ZERO) THEN
             DO N = 1,NEXP
                PA (N) = PA (N) - AB
             END DO
         END IF

         DO J = 0,SHELLX  
         DO I = 0,SHELLP  
         DO N = 1,NEXP 
            INTSCR (N,I,J) = ZERO
         END DO
         END DO
         END DO

C
C
C             ...ready!
C
C
         RETURN
         END
