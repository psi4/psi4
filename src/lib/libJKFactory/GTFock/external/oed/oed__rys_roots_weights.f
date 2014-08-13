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
         SUBROUTINE  OED__RYS_ROOTS_WEIGHTS
     +
     +                    ( NT,NTGQP,
     +                      NGQP,NMOM,
     +                      TVAL,RYSZERO,
     +                      FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                      A,B,
     +                      MOM,
     +                      DIA,OFF,
     +                      ROW1,ROW2,
     +
     +                             RTS,
     +                             WTS )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__RYS_ROOTS_WEIGHTS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : OED__RYS_1_ROOTS_WEIGHTS
C                OED__RYS_2_ROOTS_WEIGHTS
C                OED__RYS_3_ROOTS_WEIGHTS
C                OED__RYS_4_ROOTS_WEIGHTS
C                OED__RYS_5_ROOTS_WEIGHTS
C                OED__RYS_X_ROOTS_WEIGHTS
C  DESCRIPTION : This routine calculates NGQP-point Gaussian quadrature
C                rules on [0,1] over the Rys weight functions:
C
C
C                                       exp(-T*x)
C                             W   (x) = ---------
C                              Rys      2*sqrt(x)
C
C
C                for a set of NT T-exponents. Special interpolation
C                routines are provided for low number of roots and
C                weigths (NGQP < 6). On exit, NT x NGQP = NTGQP roots
C                and weights have been produced.
C
C
C                  Input:
C
C                    NT           =  # of T-exponents
C                    NTGQP        =  # of roots times # of T-exponents
C                    NGQP         =  # of gaussian quadrature points
C                                    (roots)
C                    NMOM         =  # of necessary moment integrals
C                                    to calculate the quadrature roots
C                    TVAL         =  the T-exponents
C                    RYSZERO      =  will hold the zeroth Rys moments
C                                    for all T-exponents
C                    FTABLE       =  Fm (T) table for interpolation
C                                    in low T region
C                    MGRID        =  maximum m in Fm (T) table
C                    NGRID        =  # of T's for which Fm (T) table
C                                    was set up
C                    TMAX         =  maximum T in Fm (T) table
C                    TSTEP        =  difference between two consecutive
C                                    T's in Fm (T) table
C                    TVSTEP       =  Inverse of TSTEP
C                    A,B          =  will contain the recurrence
C                                    coefficients for the auxillary
C                                    polynomials
C                    MOM          =  will contain the normed auxillary
C                                    polynomial modified moments
C                    DIA,OFF      =  will contain the diagonal and
C                                    offdiagonal elements of the
C                                    tridiagonal symmetric terminal
C                                    matrix
C                    ROW1,ROW2    =  first,second row intermediates.
C                                    Will be used to evaluate the
C                                    tridiagonal elements of the
C                                    symmetric terminal matrix in an
C                                    efficient way using Sack and
C                                    Donovan's method
C
C                  Output:
C
C                    RTS          =  the roots array
C                    WTS          =  the weights array
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

         INTEGER    MGRID
         INTEGER    N
         INTEGER    NGQP
         INTEGER    NGRID
         INTEGER    NMOM
         INTEGER    NT,NTGQP
         INTEGER    TGRID

         DOUBLE PRECISION   DELTA
         DOUBLE PRECISION   T,TMAX,TSTEP,TVSTEP
         DOUBLE PRECISION   ZERO,SIXTH,FIFTH,FOURTH,THIRD,HALF,ONE,PI

         DOUBLE PRECISION   A    (1:NMOM)
         DOUBLE PRECISION   B    (2:NMOM)
         DOUBLE PRECISION   MOM  (1:NMOM)
         DOUBLE PRECISION   DIA  (1:NGQP)
         DOUBLE PRECISION   OFF  (1:NGQP)
         DOUBLE PRECISION   ROW1 (1:NMOM)
         DOUBLE PRECISION   ROW2 (1:NMOM-1)

         DOUBLE PRECISION   RYSZERO (1:NT)
         DOUBLE PRECISION   TVAL    (1:NT)

         DOUBLE PRECISION   RTS     (1:NTGQP)
         DOUBLE PRECISION   WTS     (1:NTGQP)

         DOUBLE PRECISION   FTABLE (0:MGRID,0:NGRID)

         PARAMETER  (ZERO    =  0.D0)
         PARAMETER  (SIXTH   =  0.166666666666667D0)
         PARAMETER  (FIFTH   =  0.2D0)
         PARAMETER  (FOURTH  =  0.25D0)
         PARAMETER  (THIRD   =  0.333333333333333D0)
         PARAMETER  (HALF    =  0.5D0)
         PARAMETER  (ONE     =  1.D0)
         PARAMETER  (PI      =  3.141592653589793D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...proceed according to the number of roots and weights
C                required.
C
C
         GOTO (1000,2000,3000,4000,5000,6000)  MIN (NGQP,6)
C
C
C             ...# of roots and weights < 6
C
C
 1000    CALL  OED__RYS_1_ROOTS_WEIGHTS (NT,TVAL,          RTS,WTS)
         RETURN
 2000    CALL  OED__RYS_2_ROOTS_WEIGHTS (NT,NTGQP,TVAL,    RTS,WTS)
         RETURN
 3000    CALL  OED__RYS_3_ROOTS_WEIGHTS (NT,NTGQP,TVAL,    RTS,WTS)
         RETURN
 4000    CALL  OED__RYS_4_ROOTS_WEIGHTS (NT,NTGQP,TVAL,    RTS,WTS)
         RETURN
 5000    CALL  OED__RYS_5_ROOTS_WEIGHTS (NT,NTGQP,TVAL,    RTS,WTS)
         RETURN
C
C
C             ...# of roots and weights >= 6. Accumulate all zeroth
C                Rys moments and call the general routine.
C
C
 6000    DO N = 1,NT
            T = TVAL (N)
            IF (T.EQ.ZERO) THEN
                RYSZERO (N) = ONE
            ELSE IF (T.LE.TMAX) THEN
                TGRID = INT (T * TVSTEP + HALF)
                DELTA = TGRID * TSTEP - T

                RYSZERO (N) = ((((( FTABLE (6,TGRID)  * DELTA * SIXTH
     +                            + FTABLE (5,TGRID)) * DELTA * FIFTH
     +                            + FTABLE (4,TGRID)) * DELTA * FOURTH
     +                            + FTABLE (3,TGRID)) * DELTA * THIRD
     +                            + FTABLE (2,TGRID)) * DELTA * HALF
     +                            + FTABLE (1,TGRID)) * DELTA
     +                            + FTABLE (0,TGRID)
            ELSE
                RYSZERO (N) = HALF * DSQRT (PI/T)
            END IF
         END DO

         CALL  OED__RYS_X_ROOTS_WEIGHTS
     +
     +              ( NT,NTGQP,
     +                NGQP,NMOM,
     +                TVAL,RYSZERO,
     +                A,B,
     +                MOM,
     +                DIA,OFF,
     +                ROW1,ROW2,
     +
     +                          RTS,
     +                          WTS )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
