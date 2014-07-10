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
         SUBROUTINE  OED__NAI_1D_SHDERV_INTEGRALS
     +
     +                    ( MGIJCEN,
     +                      SHELLA,SHELLB,
     +                      DERA,DERB,
     +                      EXP2A,EXP2B,
     +                      INT1D,
     +
     +                               OUT1D )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_1D_SHDERV_INTEGRALS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation performs a single shell derivation
C                step on the A and/or the B shell for input 1D nuclear
C                attraction integrals:
C
C
C            I'(n,a,b) = delta (DERA,1) *
C                              [-a*I(n,a-1,b) + 2*expa(n)*I(n,a+1,b)]
C                      + delta (DERB,1) *
C                              [-b*I(n,a,b-1) + 2*expb(n)*I(n,a,b+1)]
C
C
C                and returns the result in a separate array. The
C                derivatives of the 1D integrals are calculated for all
C                roots and the present set {n} of exponent pairs.
C
C
C                  Input:
C
C                    MGIJCEN     =  # of roots times # of ij primitive
C                                   index pairs times # of nuclear
C                                   attraction centers
C                    SHELLx      =  maximum shell type for centers
C                                   x = A and B after differentiation
C                    DERx        =  indicator, if differentiation is
C                                   to be performed with respect to
C                                   shells x = A and/or B. Two possible
C                                   values: 0 = no differentiation,
C                                   1 = differentiate, i.e. d/dx.
C                    EXP2x       =  the totality of all MGIJCEN
C                                   exponents x 2 for both shells
C                                   x = A and B in the appropriate
C                                   order
C                    INT1D       =  all input case I 1D nuclear
C                                   attraction integrals before
C                                   differentiation.
C
C
C                  Output:
C
C                    OUT1D       =  all differentiated case I 1D nuclear
C                                   attraction integrals
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT  NONE

         LOGICAL   FIRST

         INTEGER   A,B,N
         INTEGER   AM,AP,BM,BP
         INTEGER   DERA,DERB
         INTEGER   MGIJCEN
         INTEGER   SHELLA,SHELLB

         DOUBLE PRECISION  F
         DOUBLE PRECISION  ONE,TWO

         DOUBLE PRECISION  EXP2A  (1:MGIJCEN)
         DOUBLE PRECISION  EXP2B  (1:MGIJCEN)

         DOUBLE PRECISION  INT1D (1:MGIJCEN,0:SHELLA+DERA,
     +                                      0:SHELLB+DERB)
         DOUBLE PRECISION  OUT1D (1:MGIJCEN,0:SHELLA,
     +                                      0:SHELLB)

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
         FIRST = .TRUE.

         IF (DERA.EQ.1) THEN

             IF (SHELLA.EQ.0) THEN

                 DO B = 0,SHELLB
                 DO N = 1,MGIJCEN
                    OUT1D (N,0,B) = EXP2A (N) * INT1D (N,1,B)
                 END DO
                 END DO

             ELSE IF (SHELLA.EQ.1) THEN

                 DO B = 0,SHELLB
                 DO N = 1,MGIJCEN
                    OUT1D (N,0,B) = EXP2A (N) * INT1D (N,1,B)
                    OUT1D (N,1,B) = EXP2A (N) * INT1D (N,2,B)
     +                                        - INT1D (N,0,B)
                 END DO
                 END DO

             ELSE

                 DO B = 0,SHELLB
                    DO N = 1,MGIJCEN
                       OUT1D (N,0,B) = EXP2A (N) * INT1D (N,1,B)
                       OUT1D (N,1,B) = EXP2A (N) * INT1D (N,2,B)
     +                                           - INT1D (N,0,B)
                    END DO
                    F = TWO
                    DO A = 2,SHELLA
                       AM = A - 1
                       AP = A + 1
                       DO N = 1,MGIJCEN
                          OUT1D (N,A,B) = EXP2A (N) * INT1D (N,AP,B)
     +                                          - F * INT1D (N,AM,B)
                       END DO
                       F = F + ONE
                    END DO
                 END DO

             END IF

             FIRST = .FALSE.

         END IF
C
C
C             ...derivative on center B (if any).
C
C
         IF (DERB.EQ.1) THEN

             IF (FIRST) THEN

                 DO A = 0,SHELLA
                 DO N = 1,MGIJCEN
                    OUT1D (N,A,0) = EXP2B (N) * INT1D (N,A,1)
                 END DO
                 END DO

                 IF (SHELLB.GT.0) THEN
                     DO A = 0,SHELLA
                     DO N = 1,MGIJCEN
                        OUT1D (N,A,1) = EXP2B (N) * INT1D (N,A,2)
     +                                            - INT1D (N,A,0)
                     END DO
                     END DO
                 END IF

                 IF (SHELLB.GT.1) THEN
                     F = TWO
                     DO B = 2,SHELLB
                        BM = B - 1
                        BP = B + 1
                        DO A = 0,SHELLA
                        DO N = 1,MGIJCEN
                           OUT1D (N,A,B) = EXP2B (N) * INT1D (N,A,BP)
     +                                           - F * INT1D (N,A,BM)
                        END DO
                        END DO
                        F = F + ONE
                     END DO
                 END IF

             ELSE
 
                 DO A = 0,SHELLA
                 DO N = 1,MGIJCEN
                    OUT1D (N,A,0) = OUT1D (N,A,0) +
     +                              EXP2B (N) * INT1D (N,A,1)
                 END DO
                 END DO

                 IF (SHELLB.GT.0) THEN
                     DO A = 0,SHELLA
                     DO N = 1,MGIJCEN
                        OUT1D (N,A,1) = OUT1D (N,A,1) +
     +                                  EXP2B (N) * INT1D (N,A,2)
     +                                            - INT1D (N,A,0)
                     END DO
                     END DO
                 END IF

                 IF (SHELLB.GT.1) THEN
                     F = TWO
                     DO B = 2,SHELLB
                        BM = B - 1
                        BP = B + 1
                        DO A = 0,SHELLA
                        DO N = 1,MGIJCEN
                           OUT1D (N,A,B) = OUT1D (N,A,B) +
     +                                     EXP2B (N) * INT1D (N,A,BP)
     +                                           - F * INT1D (N,A,BM)
                        END DO
                        END DO
                        F = F + ONE
                     END DO
                 END IF

             END IF

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
