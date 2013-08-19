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
         SUBROUTINE  ERD__2D_DERV_PCD_INTEGRALS
     +
     +                    ( SHELLP,SHELLC,SHELLD,
     +                      DERC,DERD,
     +                      NGQEXQ,
     +                      EXP2C,EXP2D,
     +                      INT2D,
     +
     +                               OUT2D )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__2D_DERV_PCD_INTEGRALS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation performs a single derivation step on
C                the input 2D PCD integrals:
C
C             I'(n,p,c,d) = delta (DERC,1) *
C                           [-c*I(n,p,c-1,d) + 2*expc(n)*I(n,p,c+1,d)]
C                         + delta (DERD,1) *
C                           [-d*I(n,p,c,d-1) + 2*expd(n)*I(n,p,c,d+1)]
C
C
C                and returns the result in a separate array.
C
C                The derivatives of the 2D integrals are calculated for
C                all roots and the present set of exponent quadruplets.
C                The values of DERC and DERD can be only 1 or 0. If
C                both are 1, then centers C and D are equal.
C
C
C                  Input:
C
C                    SHELLx      =  maximum shell type for centers
C                                   x = P,C,D after differentiation
C                    DERx        =  indicator, if differentiation is
C                                   to be performed with respect to
C                                   centers x = C,D. Two possible
C                                   values: 0 = no differentiation,
C                                   1 = differentiate, i.e. d/dx
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    EXP2x       =  the NGQEXQ exponents x 2 for both
C                                   centers x = C,D in the appropriate
C                                   order.
C                    INT2D       =  all input 2D PCD integrals before
C                                   differentiation.
C
C
C                  Output:
C
C                    OUT2D       =  all differentiated 2D PCD integrals
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

         LOGICAL   FIRST

         INTEGER   C,D,N,P
         INTEGER   CM,CP,DM,DP
         INTEGER   DERC,DERD
         INTEGER   NGQEXQ
         INTEGER   SHELLP,SHELLC,SHELLD

         DOUBLE PRECISION  F
         DOUBLE PRECISION  NEGONE,NEGTWO

         DOUBLE PRECISION  EXP2C (1:NGQEXQ)
         DOUBLE PRECISION  EXP2D (1:NGQEXQ)

         DOUBLE PRECISION  INT2D (1:NGQEXQ,0:SHELLP,
     +                                     0:SHELLC+DERC,
     +                                     0:SHELLD+DERD)
         DOUBLE PRECISION  OUT2D (1:NGQEXQ,0:SHELLP,
     +                                     0:SHELLC,
     +                                     0:SHELLD)
         PARAMETER  (NEGONE = -1.D0)
         PARAMETER  (NEGTWO = -2.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...derivative on center C (if any).
C
C
         FIRST = .TRUE.

         IF (DERC.EQ.1) THEN

           IF (SHELLC.EQ.0) THEN

             DO 100 D = 0,SHELLD
             DO 100 P = 0,SHELLP
             DO 100 N = 1,NGQEXQ
                OUT2D (N,P,0,D) = EXP2C (N) * INT2D (N,P,1,D)
  100        CONTINUE

           ELSE IF (SHELLC.EQ.1) THEN

             DO 110 D = 0,SHELLD
                DO 112 P = 0,SHELLP
                DO 112 N = 1,NGQEXQ
                   OUT2D (N,P,0,D) = EXP2C (N) * INT2D (N,P,1,D)
  112           CONTINUE
                DO 114 P = 0,SHELLP
                DO 114 N = 1,NGQEXQ
                   OUT2D (N,P,1,D) = EXP2C (N) * INT2D (N,P,2,D)
     +                                         - INT2D (N,P,0,D)
  114           CONTINUE
  110        CONTINUE

           ELSE

             DO 120 D = 0,SHELLD
                DO 122 P = 0,SHELLP
                DO 122 N = 1,NGQEXQ
                   OUT2D (N,P,0,D) = EXP2C (N) * INT2D (N,P,1,D)
  122           CONTINUE
                DO 124 P = 0,SHELLP
                DO 124 N = 1,NGQEXQ
                   OUT2D (N,P,1,D) = EXP2C (N) * INT2D (N,P,2,D)
     +                                         - INT2D (N,P,0,D)
  124           CONTINUE
                F = NEGTWO
                DO 126 C = 2,SHELLC
                   CM = C - 1
                   CP = C + 1
                   DO 128 P = 0,SHELLP
                   DO 128 N = 1,NGQEXQ
                      OUT2D (N,P,C,D) = EXP2C (N) * INT2D (N,P,CP,D)
     +                                        + F * INT2D (N,P,CM,D)
  128              CONTINUE
                   F = F + NEGONE
  126           CONTINUE
  120        CONTINUE

           END IF

           FIRST = .FALSE.

         END IF
C
C
C             ...derivative on center D (if any).
C
C
         IF (DERD.EQ.1) THEN

           IF (FIRST) THEN

             IF (SHELLD.EQ.0) THEN

                 DO 200 C = 0,SHELLC
                 DO 200 P = 0,SHELLP
                 DO 200 N = 1,NGQEXQ
                    OUT2D (N,P,C,0) = EXP2D (N) * INT2D (N,P,C,1)
  200            CONTINUE

             ELSE IF (SHELLD.EQ.1) THEN

                 DO 210 C = 0,SHELLC
                 DO 210 P = 0,SHELLP
                 DO 210 N = 1,NGQEXQ
                    OUT2D (N,P,C,0) = EXP2D (N) * INT2D (N,P,C,1)
  210            CONTINUE
                 DO 212 C = 0,SHELLC
                 DO 212 P = 0,SHELLP
                 DO 212 N = 1,NGQEXQ
                    OUT2D (N,P,C,1) = EXP2D (N) * INT2D (N,P,C,2)
     +                                          - INT2D (N,P,C,0)
  212            CONTINUE

             ELSE

                 DO 220 C = 0,SHELLC
                 DO 220 P = 0,SHELLP
                 DO 220 N = 1,NGQEXQ
                    OUT2D (N,P,C,0) = EXP2D (N) * INT2D (N,P,C,1)
  220            CONTINUE
                 DO 222 C = 0,SHELLC
                 DO 222 P = 0,SHELLP
                 DO 222 N = 1,NGQEXQ
                    OUT2D (N,P,C,1) = EXP2D (N) * INT2D (N,P,C,2)
     +                                          - INT2D (N,P,C,0)
  222            CONTINUE
                 F = NEGTWO
                 DO 224 D = 2,SHELLD
                    DM = D - 1
                    DP = D + 1
                    DO 226 C = 0,SHELLC
                    DO 226 P = 0,SHELLP
                    DO 226 N = 1,NGQEXQ
                       OUT2D (N,P,C,D) = EXP2D (N) * INT2D (N,P,C,DP)
     +                                         + F * INT2D (N,P,C,DM)
  226               CONTINUE
                    F = F + NEGONE
  224            CONTINUE

             END IF

           ELSE

             IF (SHELLD.EQ.0) THEN

                 DO 230 C = 0,SHELLC
                 DO 230 P = 0,SHELLP
                 DO 230 N = 1,NGQEXQ
                    OUT2D (N,P,C,0) = OUT2D (N,P,C,0) +
     +                                EXP2D (N) * INT2D (N,P,C,1)
  230            CONTINUE

             ELSE IF (SHELLD.EQ.1) THEN

                 DO 240 C = 0,SHELLC
                 DO 240 P = 0,SHELLP
                 DO 240 N = 1,NGQEXQ
                    OUT2D (N,P,C,0) = OUT2D (N,P,C,0) +
     +                                EXP2D (N) * INT2D (N,P,C,1)
  240            CONTINUE
                 DO 242 C = 0,SHELLC
                 DO 242 P = 0,SHELLP
                 DO 242 N = 1,NGQEXQ
                       OUT2D (N,P,C,1) = OUT2D (N,P,C,1) +
     +                                   EXP2D (N) * INT2D (N,P,C,2)
     +                                             - INT2D (N,P,C,0)
  242            CONTINUE

             ELSE

                 DO 250 C = 0,SHELLC
                 DO 250 P = 0,SHELLP
                 DO 250 N = 1,NGQEXQ
                    OUT2D (N,P,C,0) = OUT2D (N,P,C,0) +
     +                                EXP2D (N) * INT2D (N,P,C,1)
  250            CONTINUE
                 DO 252 C = 0,SHELLC
                 DO 252 P = 0,SHELLP
                 DO 252 N = 1,NGQEXQ
                       OUT2D (N,P,C,1) = OUT2D (N,P,C,1) +
     +                                   EXP2D (N) * INT2D (N,P,C,2)
     +                                             - INT2D (N,P,C,0)
  252            CONTINUE
                 F = NEGTWO
                 DO 254 D = 2,SHELLD
                    DM = D - 1
                    DP = D + 1
                    DO 256 C = 0,SHELLC
                    DO 256 P = 0,SHELLP
                    DO 256 N = 1,NGQEXQ
                       OUT2D (N,P,C,D) = OUT2D (N,P,C,D) +
     +                                   EXP2D (N) * INT2D (N,P,C,DP)
     +                                         + F * INT2D (N,P,C,DM)
  256               CONTINUE
                    F = F + NEGONE
  254            CONTINUE

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
