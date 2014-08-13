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
         SUBROUTINE  ERD__RYS_2_ROOTS_WEIGHTS
     +
     +                    ( NT,NTGQP,
     +                      TVAL,
     +
     +                             RTS,
     +                             WTS )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__RYS_2_ROOTS_WEIGHTS
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation returns Rys polynomial roots and weights
C                in case the number of roots and weights required
C                is = 2. All T's are treated at once so the complete
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
C                    NTGQP        =  # of roots times # of T-exponents
C                                    (= 2 * NT)
C                    TVAL         =  the set of NT T-exponents defining
C                                    the Rys weight functions
C
C                  Output:
C
C                    RTS          =  all NTGQP quadrature roots
C                    WTS          =  all NTGQP quadrature weights
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

         INTEGER    M,N
         INTEGER    NT,NTGQP
         INTEGER    TCASE

         INTEGER    JUMP2 (1:41)

         DOUBLE PRECISION   E,T,X,Y
         DOUBLE PRECISION   F1
         DOUBLE PRECISION   PI4
         DOUBLE PRECISION   R1,R2
         DOUBLE PRECISION   R12,R22
         DOUBLE PRECISION   T2MAX
         DOUBLE PRECISION   W1,W2
         DOUBLE PRECISION   W22
         DOUBLE PRECISION   HALF,ONE,TWO,FOUR,HALFP7

         DOUBLE PRECISION   TVAL (1:NT)
         DOUBLE PRECISION   RTS  (1:NTGQP)
         DOUBLE PRECISION   WTS  (1:NTGQP)

         PARAMETER  (R12 = 2.75255128608411D-01)
         PARAMETER  (R22 = 2.72474487139158D+00)
         PARAMETER  (W22 = 9.17517095361369D-02)

         PARAMETER  (PI4 = 7.85398163397448D-01)

         PARAMETER  (T2MAX = 41.D0)

         PARAMETER  (HALF    = 0.5D0)
         PARAMETER  (ONE     = 1.D0)
         PARAMETER  (TWO     = 2.D0)
         PARAMETER  (FOUR    = 4.D0)
         PARAMETER  (HALFP7  = 7.5D0)

         DATA  JUMP2  /1,2,2,3,3,4,4,4,4,4,
     +                 5,5,5,5,5,6,6,6,6,6,
     +                 6,6,6,6,6,6,6,6,6,6,
     +                 6,6,6,7,7,7,7,7,7,7,
     +                 8/
C
C
C------------------------------------------------------------------------
C
C
C                 ********************************
C             ... *  # of roots and weights = 2  *
C                 ********************************
C
C
      M = 1

      DO 200 N = 1,NT

         T = TVAL (N)

         IF (T. LE. 3.0D-07) THEN
C
C
C             ...T-range: T essentially 0
C
C
             R1 = 1.30693606237085D-01 - 2.90430236082028D-02 * T
             R2 = 2.86930639376291D+00 - 6.37623643058102D-01 * T

             WTS (M)   = 6.52145154862545D-01 - 1.22713621927067D-01 * T
             WTS (M+1) = 3.47854845137453D-01 - 2.10619711404725D-01 * T
             RTS (M)   = R1 / (ONE + R1)
             RTS (M+1) = R2 / (ONE + R2)

             M = M + 2
             GOTO 200
         END IF

         TCASE = INT ( MIN (T+ONE,T2MAX))

         GOTO (2100,2200,2300,2400,2500,2600,2700,2800) JUMP2 (TCASE)
C
C
C             ...T-range: 0 < T < 1
C
C
 2100    F1 = ((((((((-8.36313918003957D-08  * T
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

         R1 = (((((((-2.35234358048491D-09  * T
     +               +2.49173650389842D-08) * T
     +               -4.558315364581D-08)   * T
     +               -2.447252174587D-06)   * T
     +               +4.743292959463D-05)   * T
     +               -5.33184749432408D-04) * T
     +               +4.44654947116579D-03) * T
     +               -2.90430236084697D-02) * T
     +               +1.30693606237085D-01

         R2 = (((((((-2.47404902329170D-08  * T
     +               +2.36809910635906D-07) * T
     +               +1.835367736310D-06)   * T
     +               -2.066168802076D-05)   * T
     +               -1.345693393936D-04)   * T
     +               -5.88154362858038D-05) * T
     +               +5.32735082098139D-02) * T
     +               -6.37623643056745D-01) * T
     +               +2.86930639376289D+00

         W2 = ((F1 - W1) * R1 + F1) * (ONE + R2) / (R2 - R1)

         WTS (M)   = W1 - W2
         WTS (M+1) = W2
         RTS (M)   = R1 / (ONE + R1)
         RTS (M+1) = R2 / (ONE + R2)

         M = M + 2

         GOTO 200
C
C
C             ...T-range: 1 =< T < 3
C
C
 2200    X = T - TWO

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

         R1 = (((((((((-6.36859636616415D-12  * X
     +                 +8.47417064776270D-11) * X
     +                 -5.152207846962D-10)   * X
     +                 -3.846389873308D-10)   * X
     +                 +8.472253388380D-08)   * X
     +                 -1.85306035634293D-06) * X
     +                 +2.47191693238413D-05) * X
     +                 -2.49018321709815D-04) * X
     +                 +2.19173220020161D-03) * X
     +                 -1.63329339286794D-02) * X
     +                 +8.68085688285261D-02

         R2 = ((((((((( 1.45331350488343D-10  * X
     +                 +2.07111465297976D-09) * X
     +                 -1.878920917404D-08)   * X
     +                 -1.725838516261D-07)   * X
     +                 +2.247389642339D-06)   * X
     +                 +9.76783813082564D-06) * X
     +                 -1.93160765581969D-04) * X
     +                 -1.58064140671893D-03) * X
     +                 +4.85928174507904D-02) * X
     +                 -4.30761584997596D-01) * X
     +                 +1.80400974537950D+00

         W2 = ((F1 - W1) * R1 + F1) * (ONE + R2) / (R2 - R1)

         WTS (M)   = W1 - W2
         WTS (M+1) = W2
         RTS (M)   = R1 / (ONE + R1)
         RTS (M+1) = R2 / (ONE + R2)

         M = M + 2

         GOTO 200
C
C
C             ...T-range: 3 =< T < 5
C
C
 2300    X = T - FOUR

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

         R1 = ((((((((-4.11560117487296D-12  * X
     +                +7.10910223886747D-11) * X
     +                -1.73508862390291D-09) * X
     +                +5.93066856324744D-08) * X
     +                -9.76085576741771D-07) * X
     +                +1.08484384385679D-05) * X
     +                -1.12608004981982D-04) * X
     +                +1.16210907653515D-03) * X
     +                -9.89572595720351D-03) * X
     +                +6.12589701086408D-02

         R2 = (((((((((-1.80555625241001D-10  * X
     +                 +5.44072475994123D-10) * X
     +                 +1.603498045240D-08)   * X
     +                 -1.497986283037D-07)   * X
     +                 -7.017002532106D-07)   * X
     +                 +1.85882653064034D-05) * X
     +                 -2.04685420150802D-05) * X
     +                 -2.49327728643089D-03) * X
     +                 +3.56550690684281D-02) * X
     +                 -2.60417417692375D-01) * X
     +                 +1.12155283108289D+00

         W2 = ((F1 - W1) * R1 + F1) * (ONE + R2) / (R2 - R1)

         WTS (M)   = W1 - W2
         WTS (M+1) = W2
         RTS (M)   = R1 / (ONE + R1)
         RTS (M+1) = R2 / (ONE + R2)

         M = M + 2

         GOTO 200
C
C
C             ...T-range: 5 =< T < 10
C
C
 2400    E = DEXP (-T)
         X = ONE / T
         Y = T - HALFP7

         W1 = (((((( 4.6897511375022D-01  * X
     +              -6.9955602298985D-01) * X
     +              +5.3689283271887D-01) * X
     +              -3.2883030418398D-01) * X
     +              +2.4645596956002D-01) * X
     +              -4.9984072848436D-01) * X
     +              -3.1501078774085D-06) * E + DSQRT (PI4*X)

         F1 = (W1 - E) / (T + T)

         R1 = (((((((((((((-1.43632730148572D-16  * Y
     +                     +2.38198922570405D-16) * Y
     +                     +1.358319618800D-14)   * Y
     +                     -7.064522786879D-14)   * Y
     +                     -7.719300212748D-13)   * Y
     +                     +7.802544789997D-12)   * Y
     +                     +6.628721099436D-11)   * Y
     +                     -1.775564159743D-09)   * Y
     +                     +1.713828823990D-08)   * Y
     +                     -1.497500187053D-07)   * Y
     +                     +2.283485114279D-06)   * Y
     +                     -3.76953869614706D-05) * Y
     +                     +4.74791204651451D-04) * Y
     +                     -4.60448960876139D-03) * Y
     +                     +3.72458587837249D-02

         R2 = (((((((((((( 2.48791622798900D-14  * Y
     +                    -1.36113510175724D-13) * Y
     +                    -2.224334349799D-12)   * Y
     +                    +4.190559455515D-11)   * Y
     +                    -2.222722579924D-10)   * Y
     +                    -2.624183464275D-09)   * Y
     +                    +6.128153450169D-08)   * Y
     +                    -4.383376014528D-07)   * Y
     +                    -2.49952200232910D-06) * Y
     +                    +1.03236647888320D-04) * Y
     +                    -1.44614664924989D-03) * Y
     +                    +1.35094294917224D-02) * Y
     +                    -9.53478510453887D-02) * Y
     +                    +5.44765245686790D-01

         W2 = ((F1 - W1) * R1 + F1) * (ONE + R2) / (R2 - R1)

         WTS (M)   = W1 - W2
         WTS (M+1) = W2
         RTS (M)   = R1 / (ONE + R1)
         RTS (M+1) = R2 / (ONE + R2)

         M = M + 2

         GOTO 200
C
C
C             ...T-range: 10 =< T < 15
C
C
 2500    E = DEXP (-T)
         X = ONE / T

         W1 = (((-1.8784686463512D-01  * X
     +           +2.2991849164985D-01) * X
     +           -4.9893752514047D-01) * X
     +           -2.1916512131607D-05) * E + DSQRT (PI4*X)

         F1 = (W1 - E) / (T + T)

         R1 = ((((-1.01041157064226D-05  * T
     +            +1.19483054115173D-03) * T
     +            -6.73760231824074D-02) * T
     +            +1.25705571069895D+00) * T
     +            + (((-8.57609422987199D+03  * X
     +                 +5.91005939591842D+03) * X
     +                 -1.70807677109425D+03) * X
     +                 +2.64536689959503D+02) * X
     +            -2.38570496490846D+01) * E + R12 / (T - R12)

         R2 = ((( 3.39024225137123D-04  * T
     +           -9.34976436343509D-02) * T
     +           -4.22216483306320D+00) * T
     +           + (((-2.08457050986847D+03  * X
     +                -1.04999071905664D+03) * X
     +                +3.39891508992661D+02) * X
     +                -1.56184800325063D+02) * X
     +           +8.00839033297501D+00) * E + R22 / (T - R22)

         W2 = ((F1 - W1) * R1 + F1) * (ONE + R2) / (R2 - R1)

         WTS (M)   = W1 - W2
         WTS (M+1) = W2
         RTS (M)   = R1 / (ONE + R1)
         RTS (M+1) = R2 / (ONE + R2)

         M = M + 2

         GOTO 200
C
C
C             ...T-range: 15 =< T < 33
C
C
 2600    E = DEXP (-T)
         X = ONE / T

         W1 = (( 1.9623264149430D-01  * X
     +          -4.9695241464490D-01) * X
     +          -6.0156581186481D-05) * E + DSQRT (PI4*X)

         F1 = (W1 - E) / (T + T)

         R1 = ((((-1.14906395546354D-06  * T
     +            +1.76003409708332D-04) * T
     +            -1.71984023644904D-02) * T
     +            -1.37292644149838D-01) * T
     +            + (-4.75742064274859D+01  * X
     +               +9.21005186542857D+00) * X
     +            -2.31080873898939D-02) * E + R12 / (T - R12)

         R2 = ((( 3.64921633404158D-04  * T
     +           -9.71850973831558D-02) * T
     +           -4.02886174850252D+00) * T
     +           + (-1.35831002139173D+02  * X
     +              -8.66891724287962D+01) * X
     +           +2.98011277766958D+00) * E + R22 / (T - R22)

         W2 = ((F1 - W1) * R1 + F1) * (ONE + R2) / (R2 - R1)

         WTS (M)   = W1 - W2
         WTS (M+1) = W2
         RTS (M)   = R1 / (ONE + R1)
         RTS (M+1) = R2 / (ONE + R2)

         M = M + 2

         GOTO 200
C
C
C             ...T-range: 33 =< T < 40
C
C
 2700    E = DEXP (-T)

         W1 = DSQRT (PI4/T)

         W2 = ( 4.46857389308400D+00  * T
     +         -7.79250653461045D+01) * E + W22 * W1

         R1 = (-8.78947307498880D-01  * T
     +         +1.09243702330261D+01) * E + R12 / (T - R12)

         R2 = (-9.28903924275977D+00  * T
     +         +8.10642367843811D+01) * E + R22 / (T - R22)

         WTS (M)   = W1 - W2
         WTS (M+1) = W2
         RTS (M)   = R1 / (ONE + R1)
         RTS (M+1) = R2 / (ONE + R2)

         M = M + 2

         GOTO 200
C
C
C             ...T-range: T >= 40
C
C
 2800    W1 = DSQRT (PI4/T)
         W2 = W22 * W1

C         R1 = R12 / (T - R12)
C         R2 = R22 / (T - R22)
C         RTS (M)   = R1 / (ONE + R1)
C         RTS (M+1) = R2 / (ONE + R2)

         WTS (M)   = W1 - W2
         WTS (M+1) = W2
         RTS (M)   = R12 / T
         RTS (M+1) = R22 / T

         M = M + 2

  200 CONTINUE
C
C
C             ...ready!
C
C
      RETURN
      END
