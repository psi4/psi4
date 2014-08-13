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
         SUBROUTINE  OED__CARTESIAN_NORMS
     +
     +                    ( L,
     +
     +                         NORM )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__CARTESIAN_NORMS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation generates all partial cartesian
C                normalization factors. The cartesian normalization
C                factors are xyz-exponent dependent and are given for
C                a specific monomial as:
C
C                                    ______________________
C                     l m n         /       2^(2L)
C                    x y z -->     / -----------------------
C                                \/ (2l-1)!!(2m-1)!!(2n-1)!!
C
C
C                where L = l+m+n. The best way to deal with these
C                factors for a complete set of monomials for fixed L
C                is to split up each factor into its l-,m- and n-
C                component:
C
C                       _______        _______        _______
C                      / 2^(2l)       / 2^(2m)       / 2^(2n)
C                     / -------  *   / -------  *   / -------
C                   \/ (2l-1)!!    \/ (2m-1)!!    \/ (2n-1)!!
C
C
C                and to precalculate each possible partial factor
C                only once for the range l = 0,1,...,L, where L denotes
C                the maximum shell quantum number that can possibly
C                arise during integral evaluation.
C
C                Note, that the first two factors for l = 0 and 1
C                are equal to 1 and 2 and in fact these are the only
C                ones that are needed for s- and p-functions.
C
C
C                  Input:
C
C                       L         =  maximum shell quantum number
C
C                  Output:
C
C                       NORM (I)  =  cartesian norms from I=0,1,...,L
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

         INTEGER    I,L

         DOUBLE PRECISION  ODD,ONE,TWO

         DOUBLE PRECISION  NORM (0:L)

         PARAMETER  (ONE = 1.D0)
         PARAMETER  (TWO = 2.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...calculate the partial norms for shell = 0 up to L.
C
C
         NORM (0) = ONE
         NORM (1) = TWO

         ODD = ONE

         DO 10 I = 2,L
            ODD = ODD + TWO
            NORM (I) = TWO * NORM (I-1) / DSQRT (ODD)
   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
