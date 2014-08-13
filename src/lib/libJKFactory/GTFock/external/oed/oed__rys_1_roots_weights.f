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
         SUBROUTINE  OED__RYS_1_ROOTS_WEIGHTS
     +
     +                    ( NT,
     +                      TVAL,
     +
     +                             RTS,
     +                             WTS )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__RYS_1_ROOTS_WEIGHTS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation returns Rys polynomial roots and weights
C                in case the number of roots and weights required
C                is = 1. All T's are treated at once so the complete
C                set of roots and weights is returned.
C
C                For the moment taken essentially unchanged from the
C                GAMESS package (routine RTS123, but removing their
C                'spaghetti' code from the 70's of unreadable
C                internested IFs and GOTOs!).
C
C                One interesting aspect of the GAMESS routines is that
C                their code returns scaled roots, i.e. their roots
C                do not ly between the range 0 and 1. To get to the
C                proper roots as needed for our package, we simply
C                set:
C
C                   root (our) = root (gamess) / (1 + root (games))
C
C
C                  Input:
C
C                    NT           =  # of T-exponents
C                    TVAL         =  the set of NT T-exponents defining
C                                    the Rys weight functions
C
C                  Output:
C
C                    RTS          =  all NT quadrature roots
C                    WTS          =  all NT quadrature weights
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         INTEGER    N
         INTEGER    NT
         INTEGER    TCASE

         INTEGER    JUMP1 (1:34)

         DOUBLE PRECISION   E,T,X
         DOUBLE PRECISION   F1
         DOUBLE PRECISION   PI4
         DOUBLE PRECISION   T1MAX
         DOUBLE PRECISION   R1,W1
         DOUBLE PRECISION   FIFTH,THIRD,HALF,ONE,TWO,FOUR

         DOUBLE PRECISION   TVAL (1:NT)
         DOUBLE PRECISION   RTS  (1:NT)
         DOUBLE PRECISION   WTS  (1:NT)

         PARAMETER  (PI4     = 7.85398163397448D-01)
         PARAMETER  (T1MAX   = 34.D0)
         PARAMETER  (FIFTH   = 0.2D0)
         PARAMETER  (THIRD   = 0.333333333333333D0)
         PARAMETER  (HALF    = 0.5D0)
         PARAMETER  (ONE     = 1.D0)
         PARAMETER  (TWO     = 2.D0)
         PARAMETER  (FOUR    = 4.D0)

         DATA  JUMP1  /1,2,2,3,3,4,4,4,4,4,
     +                 5,5,5,5,5,6,6,6,6,6,
     +                 6,6,6,6,6,6,6,6,6,6,
     +                 6,6,6,7/
C
C
C------------------------------------------------------------------------
C
C
C                 ********************************
C             ... *  # of roots and weights = 1  *
C                 ********************************
C
C
      DO 100 N = 1,NT

         T = TVAL (N)

         IF (T. LE. 3.0D-07) THEN
C
C
C             ...T-range: T essentially 0
C
C
             R1 = HALF - T * FIFTH

             WTS (N) = ONE  - T * THIRD
             RTS (N) = R1 / (ONE + R1)
             GOTO 100
         END IF

         TCASE = INT ( MIN (T+ONE,T1MAX))

         GOTO (1100,1200,1300,1400,1500,1600,1700) JUMP1 (TCASE)
C
C
C             ...T-range: 0 < T < 1
C
C
 1100    F1 = ((((((((-8.36313918003957D-08  * T
     +                +1.21222603512827D-06) * T
     +                -1.15662609053481D-05) * T
     +                +9.25197374512647D-05) * T
     +                -6.40994113129432D-04) * T
     +                +3.78787044215009D-03) * T
     +                -1.85185172458485D-02) * T
     +                +7.14285713298222D-02) * T
     +                -1.99999999997023D-01) * T
     +                +3.33333333333318D-01

         W1 = (T + T) * F1 + DEXP (-T)
         R1 = F1 / (W1 - F1)

         WTS (N) = W1
         RTS (N) = R1 / (ONE + R1)

         GOTO 100
C
C
C             ...T-range: 1 =< T < 3
C
C
 1200    X = T - TWO

         F1 = ((((((((((-1.61702782425558D-10  * X
     +                  +1.96215250865776D-09) * X
     +                  -2.14234468198419D-08) * X
     +                  +2.17216556336318D-07) * X
     +                  -1.98850171329371D-06) * X
     +                  +1.62429321438911D-05) * X
     +                  -1.16740298039895D-04) * X
     +                  +7.24888732052332D-04) * X
     +                  -3.79490003707156D-03) * X
     +                  +1.61723488664661D-02) * X
     +                  -5.29428148329736D-02) * X
     +                  +1.15702180856167D-01

         W1 = (T + T) * F1 + DEXP (-T)
         R1 = F1 / (W1 - F1)

         WTS (N) = W1
         RTS (N) = R1 / (ONE + R1)

         GOTO 100
C
C
C             ...T-range: 3 =< T < 5
C
C
 1300    X = T - FOUR

         F1 = ((((((((((-2.62453564772299D-11  * X
     +                  +3.24031041623823D-10) * X
     +                  -3.614965656163D-09)   * X
     +                  +3.760256799971D-08)   * X
     +                  -3.553558319675D-07)   * X
     +                  +3.022556449731D-06)   * X
     +                  -2.290098979647D-05)   * X
     +                  +1.526537461148D-04)   * X
     +                  -8.81947375894379D-04) * X
     +                  +4.33207949514611D-03) * X
     +                  -1.75257821619926D-02) * X
     +                  +5.28406320615584D-02

         W1 = (T + T) * F1 + DEXP (-T)
         R1 = F1 / (W1 - F1)

         WTS (N) = W1
         RTS (N) = R1 / (ONE + R1)

         GOTO 100
C
C
C             ...T-range: 5 =< T < 10
C
C
 1400    E = DEXP (-T)
         X = ONE / T

         W1 = (((((( 4.6897511375022D-01  * X
     +              -6.9955602298985D-01) * X
     +              +5.3689283271887D-01) * X
     +              -3.2883030418398D-01) * X
     +              +2.4645596956002D-01) * X
     +              -4.9984072848436D-01) * X
     +              -3.1501078774085D-06) * E + DSQRT (PI4*X)

         F1 = (W1 - E) / (T + T)
         R1 = F1 / (W1 - F1)

         WTS (N) = W1
         RTS (N) = R1 / (ONE + R1)

         GOTO 100
C
C
C             ...T-range: 10 =< T < 15
C
C
 1500    E = DEXP (-T)
         X = ONE / T

         W1 = (((-1.8784686463512D-01  * X
     +           +2.2991849164985D-01) * X
     +           -4.9893752514047D-01) * X
     +           -2.1916512131607D-05) * E + DSQRT (PI4*X)

         F1 = (W1 - E) / (T + T)
         R1 = F1 / (W1 - F1)

         WTS (N) = W1
         RTS (N) = R1 / (ONE + R1)

         GOTO 100
C
C
C             ...T-range: 15 =< T < 33
C
C
 1600    E = DEXP (-T)
         X = ONE / T

         W1 = (( 1.9623264149430D-01  * X
     +          -4.9695241464490D-01) * X
     +          -6.0156581186481D-05) * E + DSQRT (PI4*X)

         F1 = (W1 - E) / (T + T)
         R1 = F1 / (W1 - F1)

         WTS (N) = W1
         RTS (N) = R1 / (ONE + R1)

         GOTO 100
C
C
C             ...T-range: T >= 33
C
C
 1700    WTS (N) = DSQRT (PI4/T)
C         R1 = HALF / (T - HALF)
C         RTS (N) = R1 / (ONE + R1)
         RTS (N) = HALF / T

  100 CONTINUE
C
C
C             ...ready!
C
C
      RETURN
      END
