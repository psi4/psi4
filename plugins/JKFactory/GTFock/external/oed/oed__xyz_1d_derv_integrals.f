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
         SUBROUTINE  OED__XYZ_1D_DERV_INTEGRALS
     +
     +                    ( SHELLA,SHELLB,
     +                      DERA,DERB,
     +                      NEXP,
     +                      EXP2A,EXP2B,
     +                      PINVA,PINVB,                   ! Watson Added
     +                      OVRLP,                         ! Watson Added
     +                      INT1D,
     +
     +                               OUT1D )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__OVL_1D_DERV_INTEGRALS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation performs a single derivation step on
C                the input 1D overlap integrals:
C
C            I'(n,a,b) = delta (DERA,1) *
C                              [-a*I(n,a-1,b) + 2*expa(n)*I(n,a+1,b)]
C                      + delta (DERB,1) *
C                              [-b*I(n,a,b-1) + 2*expb(n)*I(n,a,b+1)]
C
C                and returns the result in a separate array.
C
C                The derivatives of the 1D integrals are calculated for
C                all the present set of exponent pairs. The values of
C                DERA and DERB can be only 1 or 0 and only one of
C                them can be equal to 1. If both of them were equal
C                to 1 this would indicate an atomic case and the above
C                derivative expression would be exactly equal to zero.
C
C
C                  Input:
C
C                    SHELLx      =  maximum shell type for centers
C                                   x = A and B after differentiation
C                    DERx        =  indicator, if differentiation is
C                                   to be performed with respect to
C                                   centers x = A and B. Two possible
C                                   values: 0 = no differentiation,
C                                   1 = differentiate, i.e. d/dx
C                    NEXP        =  current # of exponent pairs
C                                   corresponding to the contracted
C                                   shell pairs A,B
C                    EXP2x       =  all NEXP exponents x 2 for both
C                                   centers x = A and B in the
C                                   appropriate order
C                    INT1D       =  all input 1D overlap integrals
C                                   before differentiation.
C
C
C                  Output:
C
C                    OUT1D       =  all differentiated 1D overlap
C                                   integrals
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT  NONE

         INTEGER   A,B,N
         INTEGER   AM,AP,BM,BP
         INTEGER   DERA,DERB
         INTEGER   NEXP
         INTEGER   SHELLA,SHELLB

         DOUBLE PRECISION  F
         DOUBLE PRECISION  ZERO,ONE,TWO

         DOUBLE PRECISION  EXP2A  (1:NEXP)
         DOUBLE PRECISION  EXP2B  (1:NEXP)

         DOUBLE PRECISION  PINVA  (1:NEXP)
         DOUBLE PRECISION  PINVB  (1:NEXP)

         DOUBLE PRECISION  OVRLP (1:NEXP,0:SHELLA,0:SHELLB)
         DOUBLE PRECISION  INT1D (1:NEXP,0:SHELLA,0:SHELLB)
         DOUBLE PRECISION  OUT1D (1:NEXP,0:SHELLA,0:SHELLB)

         PARAMETER  (ZERO   =  0.D0)
         PARAMETER  (ONE    =  1.D0)
         PARAMETER  (TWO    =  2.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...derivative on center A (if any).
C
C
         IF (DERA.EQ.1) THEN

             IF (SHELLA.EQ.0) THEN

                 DO B = 0,SHELLB
                 DO N = 1,NEXP
                    OUT1D (N,0,B) =   EXP2A (N) * INT1D (N,1,B)
     +                              + PINVA (N) * OVRLP (N,0,B)
                 END DO
                 END DO

             ELSE IF (SHELLA.EQ.1) THEN

                 DO B = 0,SHELLB
                 DO N = 1,NEXP
                    OUT1D (N,0,B) =   EXP2A (N) * INT1D (N,1,B)
     +                              + PINVA (N) * OVRLP (N,0,B)
                    OUT1D (N,1,B) =   EXP2A (N) * INT1D (N,2,B)
     +                                          - INT1D (N,0,B)
     +                              + PINVA (N) * OVRLP (N,1,B)
                 END DO
                 END DO

             ELSE

                 DO B = 0,SHELLB
                    DO N = 1,NEXP
                       OUT1D (N,0,B) =   EXP2A (N) * INT1D (N,1,B)
     +                                 + PINVA (N) * OVRLP (N,0,B)

                       OUT1D (N,1,B) =   EXP2A (N) * INT1D (N,2,B)
     +                                             - INT1D (N,0,B)
     +                                 + PINVA (N) * OVRLP (N,1,B)

                    END DO
                    F = TWO
                    DO A = 2,SHELLA
                       AM = A - 1
                       AP = A + 1
                       DO N = 1,NEXP
                          OUT1D (N,A,B) =   EXP2A (N) * INT1D (N,AP,B)
     +                                            - F * INT1D (N,AM,B)
     +                                    + PINVA (N) * OVRLP (N,A,B)

                       END DO
                       F = F + ONE
                    END DO
                 END DO

             END IF

         END IF
C
C
C             ...derivative on center B (if any).
C
C
         IF (DERB.EQ.1) THEN

             DO A = 0,SHELLA
             DO N = 1,NEXP
                OUT1D (N,A,0) = EXP2B (N) * INT1D (N,A,1)
     +                        + PINVB (N) * OVRLP (N,A,0)  ! Watson Added

             END DO
             END DO

             IF (SHELLB.GT.0) THEN
                 DO A = 0,SHELLA
                 DO N = 1,NEXP
                    OUT1D (N,A,1) = EXP2B (N) * INT1D (N,A,2)
     +                                        - INT1D (N,A,0)
     +                            + PINVB (N) * OVRLP (N,A,1)  ! Watson Added

                 END DO
                 END DO
             END IF

             IF (SHELLB.GT.1) THEN
                 F = TWO
                 DO B = 2,SHELLB
                    BM = B - 1
                    BP = B + 1
                    DO A = 0,SHELLA
                    DO N = 1,NEXP
                       OUT1D (N,A,B) = EXP2B (N) * INT1D (N,A,BP)
     +                                       - F * INT1D (N,A,BM)
     +                               + PINVB (N) * OVRLP (N,A,B)  ! Watson Added

                    END DO
                    END DO
                    F = F + ONE
                 END DO
             END IF

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
