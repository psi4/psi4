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
         SUBROUTINE  ERD__2D_DERV_ABCD_INTEGRALS
     +
     +                    ( SHELLA,SHELLB,SHELLC,SHELLD,
     +                      DERA,DERB,DERC,DERD,
     +                      NGQEXQ,
     +                      EXP2A,EXP2B,EXP2C,EXP2D,
     +                      INT2D,
     +
     +                               OUT2D )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__2D_DERV_ABCD_INTEGRALS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation performs a single derivation step on
C                the input 2D ABCD integrals:
C
C         I'(n,a,b,c,d) = delta (DERA,1) *
C                         [-a*I(n,a-1,b,c,d) + 2*expa(n)*I(n,a+1,b,c,d)]
C                       + delta (DERB,1) *
C                         [-b*I(n,a,b-1,c,d) + 2*expb(n)*I(n,a,b+1,c,d)]
C                       + delta (DERC,1) *
C                         [-c*I(n,a,b,c-1,d) + 2*expc(n)*I(n,a,b,c+1,d)]
C                       + delta (DERD,1) *
C                         [-d*I(n,a,b,c,d-1) + 2*expd(n)*I(n,a,b,c,d+1)]
C
C                and returns the result in a separate array.
C
C                The derivatives of the 2D integrals are calculated for
C                all roots and the present set of exponent quadruplets.
C                The values of DERA,DERB,DERC and DERD can be only 1
C                or 0. A set of 1's among these variables indicates
C                that these centers were identified as equal.
C
C
C                  Input:
C
C                    SHELLx      =  maximum shell type for centers
C                                   x = A,B,C,D after differentiation
C                    DERx        =  indicator, if differentiation is
C                                   to be performed with respect to
C                                   centers x = A,B,C,D. Two possible
C                                   values: 0 = no differentiation,
C                                   1 = differentiate, i.e. d/dx
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    EXP2x       =  the NGQEXQ exponents x 2 for all
C                                   centers x = A,B,C,D in the
C                                   appropriate order.
C                    INT2D       =  all input 2D ABCD integrals before
C                                   differentiation.
C
C
C                  Output:
C
C                    OUT2D       =  all differentiated 2D ABCD
C                                   integrals
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

         INTEGER   A,B,C,D,N
         INTEGER   AM,AP,BM,BP,CM,CP,DM,DP
         INTEGER   DERA,DERB,DERC,DERD
         INTEGER   NGQEXQ
         INTEGER   SHELLA,SHELLB,SHELLC,SHELLD

         DOUBLE PRECISION  F
         DOUBLE PRECISION  NEGONE,NEGTWO

         DOUBLE PRECISION  EXP2A (1:NGQEXQ)
         DOUBLE PRECISION  EXP2B (1:NGQEXQ)
         DOUBLE PRECISION  EXP2C (1:NGQEXQ)
         DOUBLE PRECISION  EXP2D (1:NGQEXQ)

         DOUBLE PRECISION  INT2D (1:NGQEXQ,0:SHELLA+DERA,
     +                                     0:SHELLB+DERB,
     +                                     0:SHELLC+DERC,
     +                                     0:SHELLD+DERD)
         DOUBLE PRECISION  OUT2D (1:NGQEXQ,0:SHELLA,
     +                                     0:SHELLB,
     +                                     0:SHELLC,
     +                                     0:SHELLD)
         PARAMETER  (NEGONE = -1.D0)
         PARAMETER  (NEGTWO = -2.D0)
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

             DO 100 D = 0,SHELLD
             DO 100 C = 0,SHELLC
             DO 100 B = 0,SHELLB
             DO 100 N = 1,NGQEXQ
                OUT2D (N,0,B,C,D) = EXP2A (N) * INT2D (N,1,B,C,D)
  100        CONTINUE

           ELSE IF (SHELLA.EQ.1) THEN

             DO 110 D = 0,SHELLD
             DO 110 C = 0,SHELLC
             DO 110 B = 0,SHELLB
             DO 110 N = 1,NGQEXQ
                OUT2D (N,0,B,C,D) = EXP2A (N) * INT2D (N,1,B,C,D)
                OUT2D (N,1,B,C,D) = EXP2A (N) * INT2D (N,2,B,C,D)
     +                                        - INT2D (N,0,B,C,D)
  110        CONTINUE

           ELSE

             DO 120 D = 0,SHELLD
             DO 120 C = 0,SHELLC
             DO 120 B = 0,SHELLB
                DO 122 N = 1,NGQEXQ
                   OUT2D (N,0,B,C,D) = EXP2A (N) * INT2D (N,1,B,C,D)
                   OUT2D (N,1,B,C,D) = EXP2A (N) * INT2D (N,2,B,C,D)
     +                                           - INT2D (N,0,B,C,D)
  122           CONTINUE
                F = NEGTWO
                DO 124 A = 2,SHELLA
                   AM = A - 1
                   AP = A + 1
                   DO 126 N = 1,NGQEXQ
                      OUT2D (N,A,B,C,D) = EXP2A (N) * INT2D (N,AP,B,C,D)
     +                                          + F * INT2D (N,AM,B,C,D)
  126              CONTINUE
                   F = F + NEGONE
  124           CONTINUE
  120        CONTINUE

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

           IF (SHELLB.EQ.0) THEN

             DO 200 D = 0,SHELLD
             DO 200 C = 0,SHELLC
             DO 200 A = 0,SHELLA
             DO 200 N = 1,NGQEXQ
                OUT2D (N,A,0,C,D) = EXP2B (N) * INT2D (N,A,1,C,D)
  200        CONTINUE

           ELSE IF (SHELLB.EQ.1) THEN

             DO 210 D = 0,SHELLD
             DO 210 C = 0,SHELLC
                DO 212 A = 0,SHELLA
                DO 212 N = 1,NGQEXQ
                   OUT2D (N,A,0,C,D) = EXP2B (N) * INT2D (N,A,1,C,D)
  212           CONTINUE
                DO 214 A = 0,SHELLA
                DO 214 N = 1,NGQEXQ
                   OUT2D (N,A,1,C,D) = EXP2B (N) * INT2D (N,A,2,C,D)
     +                                           - INT2D (N,A,0,C,D)
  214           CONTINUE
  210        CONTINUE

           ELSE

             DO 220 D = 0,SHELLD
             DO 220 C = 0,SHELLC
                DO 222 A = 0,SHELLA
                DO 222 N = 1,NGQEXQ
                   OUT2D (N,A,0,C,D) = EXP2B (N) * INT2D (N,A,1,C,D)
  222           CONTINUE
                DO 224 A = 0,SHELLA
                DO 224 N = 1,NGQEXQ
                   OUT2D (N,A,1,C,D) = EXP2B (N) * INT2D (N,A,2,C,D)
     +                                           - INT2D (N,A,0,C,D)
  224           CONTINUE
                F = NEGTWO
                DO 226 B = 2,SHELLB
                   BM = B - 1
                   BP = B + 1
                   DO 228 A = 0,SHELLA
                   DO 228 N = 1,NGQEXQ
                      OUT2D (N,A,B,C,D) = EXP2B (N) * INT2D (N,A,BP,C,D)
     +                                          + F * INT2D (N,A,BM,C,D)
  228              CONTINUE
                   F = F + NEGONE
  226           CONTINUE
  220        CONTINUE

           END IF

           FIRST = .FALSE.

          ELSE
 
           IF (SHELLB.EQ.0) THEN

             DO 230 D = 0,SHELLD
             DO 230 C = 0,SHELLC
             DO 230 A = 0,SHELLA
             DO 230 N = 1,NGQEXQ
                OUT2D (N,A,0,C,D) = OUT2D (N,A,0,C,D) +
     +                              EXP2B (N) * INT2D (N,A,1,C,D)
  230        CONTINUE

           ELSE IF (SHELLB.EQ.1) THEN

             DO 240 D = 0,SHELLD
             DO 240 C = 0,SHELLC
                DO 242 A = 0,SHELLA
                DO 242 N = 1,NGQEXQ
                   OUT2D (N,A,0,C,D) = OUT2D (N,A,0,C,D) +
     +                                 EXP2B (N) * INT2D (N,A,1,C,D)
  242           CONTINUE
                DO 244 A = 0,SHELLA
                DO 244 N = 1,NGQEXQ
                   OUT2D (N,A,1,C,D) = OUT2D (N,A,1,C,D) +
     +                                 EXP2B (N) * INT2D (N,A,2,C,D)
     +                                           - INT2D (N,A,0,C,D)
  244           CONTINUE
  240        CONTINUE

           ELSE

             DO 250 D = 0,SHELLD
             DO 250 C = 0,SHELLC
                DO 252 A = 0,SHELLA
                DO 252 N = 1,NGQEXQ
                   OUT2D (N,A,0,C,D) = OUT2D (N,A,0,C,D) +
     +                                 EXP2B (N) * INT2D (N,A,1,C,D)
  252           CONTINUE
                DO 254 A = 0,SHELLA
                DO 254 N = 1,NGQEXQ
                   OUT2D (N,A,1,C,D) = OUT2D (N,A,1,C,D) +
     +                                 EXP2B (N) * INT2D (N,A,2,C,D)
     +                                           - INT2D (N,A,0,C,D)
  254           CONTINUE
                F = NEGTWO
                DO 256 B = 2,SHELLB
                   BM = B - 1
                   BP = B + 1
                   DO 258 A = 0,SHELLA
                   DO 258 N = 1,NGQEXQ
                      OUT2D (N,A,B,C,D) = OUT2D (N,A,B,C,D) +
     +                                    EXP2B (N) * INT2D (N,A,BP,C,D)
     +                                          + F * INT2D (N,A,BM,C,D)
  258              CONTINUE
                   F = F + NEGONE
  256           CONTINUE
  250        CONTINUE

           END IF

          END IF

         END IF
C
C
C             ...derivative on center C (if any).
C
C
         IF (DERC.EQ.1) THEN

          IF (FIRST) THEN

           IF (SHELLC.EQ.0) THEN

             DO 300 D = 0,SHELLD
             DO 300 B = 0,SHELLB
             DO 300 A = 0,SHELLA
             DO 300 N = 1,NGQEXQ
                OUT2D (N,A,B,0,D) = EXP2C (N) * INT2D (N,A,B,1,D)
  300        CONTINUE

           ELSE IF (SHELLC.EQ.1) THEN

             DO 310 D = 0,SHELLD
                DO 312 B = 0,SHELLB
                DO 312 A = 0,SHELLA
                DO 312 N = 1,NGQEXQ
                   OUT2D (N,A,B,0,D) = EXP2C (N) * INT2D (N,A,B,1,D)
  312           CONTINUE
                DO 314 B = 0,SHELLB
                DO 314 A = 0,SHELLA
                DO 314 N = 1,NGQEXQ
                   OUT2D (N,A,B,1,D) = EXP2C (N) * INT2D (N,A,B,2,D)
     +                                           - INT2D (N,A,B,0,D)
  314           CONTINUE
  310        CONTINUE

           ELSE

             DO 320 D = 0,SHELLD
                DO 322 B = 0,SHELLB
                DO 322 A = 0,SHELLA
                DO 322 N = 1,NGQEXQ
                   OUT2D (N,A,B,0,D) = EXP2C (N) * INT2D (N,A,B,1,D)
  322           CONTINUE
                DO 324 B = 0,SHELLB
                DO 324 A = 0,SHELLA
                DO 324 N = 1,NGQEXQ
                   OUT2D (N,A,B,1,D) = EXP2C (N) * INT2D (N,A,B,2,D)
     +                                           - INT2D (N,A,B,0,D)
  324           CONTINUE
                F = NEGTWO
                DO 326 C = 2,SHELLC
                   CM = C - 1
                   CP = C + 1
                   DO 328 B = 0,SHELLB
                   DO 328 A = 0,SHELLA
                   DO 328 N = 1,NGQEXQ
                      OUT2D (N,A,B,C,D) = EXP2C (N) * INT2D (N,A,B,CP,D)
     +                                          + F * INT2D (N,A,B,CM,D)
  328              CONTINUE
                   F = F + NEGONE
  326           CONTINUE
  320        CONTINUE

           END IF

           FIRST = .FALSE.

          ELSE

           IF (SHELLC.EQ.0) THEN

             DO 330 D = 0,SHELLD
             DO 330 B = 0,SHELLB
             DO 330 A = 0,SHELLA
             DO 330 N = 1,NGQEXQ
                OUT2D (N,A,B,0,D) = OUT2D (N,A,B,0,D) +
     +                              EXP2C (N) * INT2D (N,A,B,1,D)
  330        CONTINUE

           ELSE IF (SHELLC.EQ.1) THEN

             DO 340 D = 0,SHELLD
                DO 342 B = 0,SHELLB
                DO 342 A = 0,SHELLA
                DO 342 N = 1,NGQEXQ
                   OUT2D (N,A,B,0,D) = OUT2D (N,A,B,0,D) +
     +                                 EXP2C (N) * INT2D (N,A,B,1,D)
  342           CONTINUE
                DO 344 B = 0,SHELLB
                DO 344 A = 0,SHELLA
                DO 344 N = 1,NGQEXQ
                   OUT2D (N,A,B,1,D) = OUT2D (N,A,B,1,D) +
     +                                 EXP2C (N) * INT2D (N,A,B,2,D)
     +                                           - INT2D (N,A,B,0,D)
  344           CONTINUE
  340        CONTINUE

           ELSE

             DO 350 D = 0,SHELLD
                DO 352 B = 0,SHELLB
                DO 352 A = 0,SHELLA
                DO 352 N = 1,NGQEXQ
                   OUT2D (N,A,B,0,D) = OUT2D (N,A,B,0,D) +
     +                                 EXP2C (N) * INT2D (N,A,B,1,D)
  352           CONTINUE
                DO 354 B = 0,SHELLB
                DO 354 A = 0,SHELLA
                DO 354 N = 1,NGQEXQ
                   OUT2D (N,A,B,1,D) = OUT2D (N,A,B,1,D) +
     +                                 EXP2C (N) * INT2D (N,A,B,2,D)
     +                                           - INT2D (N,A,B,0,D)
  354           CONTINUE
                F = NEGTWO
                DO 356 C = 2,SHELLC
                   CM = C - 1
                   CP = C + 1
                   DO 358 B = 0,SHELLB
                   DO 358 A = 0,SHELLA
                   DO 358 N = 1,NGQEXQ
                      OUT2D (N,A,B,C,D) = OUT2D (N,A,B,C,D) +
     +                                    EXP2C (N) * INT2D (N,A,B,CP,D)
     +                                          + F * INT2D (N,A,B,CM,D)
  358              CONTINUE
                   F = F + NEGONE
  356           CONTINUE
  350        CONTINUE

           END IF

          END IF

         END IF
C
C
C             ...derivative on center D (if any).
C
C
         IF (DERD.EQ.1) THEN

          IF (FIRST) THEN

           IF (SHELLD.EQ.0) THEN

             DO 400 C = 0,SHELLC
             DO 400 B = 0,SHELLB
             DO 400 A = 0,SHELLA
             DO 400 N = 1,NGQEXQ
                OUT2D (N,A,B,C,0) = EXP2D (N) * INT2D (N,A,B,C,1)
  400        CONTINUE

           ELSE IF (SHELLD.EQ.1) THEN

             DO 410 C = 0,SHELLC
             DO 410 B = 0,SHELLB
             DO 410 A = 0,SHELLA
             DO 410 N = 1,NGQEXQ
                OUT2D (N,A,B,C,0) = EXP2D (N) * INT2D (N,A,B,C,1)
  410        CONTINUE
             DO 412 C = 0,SHELLC
             DO 412 B = 0,SHELLB
             DO 412 A = 0,SHELLA
             DO 412 N = 1,NGQEXQ
                OUT2D (N,A,B,C,1) = EXP2D (N) * INT2D (N,A,B,C,2)
     +                                        - INT2D (N,A,B,C,0)
  412        CONTINUE

           ELSE

             DO 420 C = 0,SHELLC
             DO 420 B = 0,SHELLB
             DO 420 A = 0,SHELLA
             DO 420 N = 1,NGQEXQ
                OUT2D (N,A,B,C,0) = EXP2D (N) * INT2D (N,A,B,C,1)
  420        CONTINUE
             DO 422 C = 0,SHELLC
             DO 422 B = 0,SHELLB
             DO 422 A = 0,SHELLA
             DO 422 N = 1,NGQEXQ
                OUT2D (N,A,B,C,1) = EXP2D (N) * INT2D (N,A,B,C,2)
     +                                        - INT2D (N,A,B,C,0)
  422        CONTINUE
             F = NEGTWO
             DO 424 D = 2,SHELLD
                DM = D - 1
                DP = D + 1
                DO 426 C = 0,SHELLC
                DO 426 B = 0,SHELLB
                DO 426 A = 0,SHELLA
                DO 426 N = 1,NGQEXQ
                   OUT2D (N,A,B,C,D) = EXP2D (N) * INT2D (N,A,B,C,DP)
     +                                       + F * INT2D (N,A,B,C,DM)
  426           CONTINUE
                F = F + NEGONE
  424        CONTINUE

           END IF

          ELSE

           IF (SHELLD.EQ.0) THEN

             DO 430 C = 0,SHELLC
             DO 430 B = 0,SHELLB
             DO 430 A = 0,SHELLA
             DO 430 N = 1,NGQEXQ
                OUT2D (N,A,B,C,0) = OUT2D (N,A,B,C,0) +
     +                              EXP2D (N) * INT2D (N,A,B,C,1)
  430        CONTINUE

           ELSE IF (SHELLD.EQ.1) THEN

             DO 440 C = 0,SHELLC
             DO 440 B = 0,SHELLB
             DO 440 A = 0,SHELLA
             DO 440 N = 1,NGQEXQ
                OUT2D (N,A,B,C,0) = OUT2D (N,A,B,C,0) +
     +                              EXP2D (N) * INT2D (N,A,B,C,1)
  440        CONTINUE
             DO 442 C = 0,SHELLC
             DO 442 B = 0,SHELLB
             DO 442 A = 0,SHELLA
             DO 442 N = 1,NGQEXQ
                OUT2D (N,A,B,C,1) = OUT2D (N,A,B,C,1) +
     +                              EXP2D (N) * INT2D (N,A,B,C,2)
     +                                        - INT2D (N,A,B,C,0)
  442        CONTINUE

           ELSE

             DO 450 C = 0,SHELLC
             DO 450 B = 0,SHELLB
             DO 450 A = 0,SHELLA
             DO 450 N = 1,NGQEXQ
                OUT2D (N,A,B,C,0) = OUT2D (N,A,B,C,0) +
     +                              EXP2D (N) * INT2D (N,A,B,C,1)
  450        CONTINUE
             DO 452 C = 0,SHELLC
             DO 452 B = 0,SHELLB
             DO 452 A = 0,SHELLA
             DO 452 N = 1,NGQEXQ
                OUT2D (N,A,B,C,1) = OUT2D (N,A,B,C,1) +
     +                              EXP2D (N) * INT2D (N,A,B,C,2)
     +                                        - INT2D (N,A,B,C,0)
  452        CONTINUE
             F = NEGTWO
             DO 454 D = 2,SHELLD
                DM = D - 1
                DP = D + 1
                DO 456 C = 0,SHELLC
                DO 456 B = 0,SHELLB
                DO 456 A = 0,SHELLA
                DO 456 N = 1,NGQEXQ
                   OUT2D (N,A,B,C,D) = OUT2D (N,A,B,C,D) +
     +                                 EXP2D (N) * INT2D (N,A,B,C,DP)
     +                                       + F * INT2D (N,A,B,C,DM)
  456           CONTINUE
                F = F + NEGONE
  454        CONTINUE

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
