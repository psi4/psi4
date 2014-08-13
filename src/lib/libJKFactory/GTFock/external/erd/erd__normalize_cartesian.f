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
         SUBROUTINE  ERD__NORMALIZE_CARTESIAN
     +
     +                    ( M,
     +                      NXYZ,
     +                      L,
     +                      NORM,
     +
     +                             BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__NORMALIZE_CARTESIAN
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation normalizes a cartesian monomial part
C                of a batch of contracted cartesian gaussian integrals.
C                The normalization factors are xyz-exponent dependent
C                and are given for a specific monomial as:
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
C                These factors are passed in argument as NORM.
C
C
C                  Input:
C
C                    M           =  # of elements not involved in the
C                                   normalization (invariant indices)
C                    NXYZ        =  # of monomials for shell L
C                    L           =  the shell type for which the
C                                   normalization will be done
C                    NORM        =  individual normalization factors
C                                   for each monomial exponent
C                    BATCH       =  batch of unnormalized integrals
C
C                  Output:
C
C                    BATCH       =  batch of normalized integrals
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         INTEGER     I,L,M,N,X,Y,Z
         INTEGER     NXYZ
         INTEGER     YBEG

         DOUBLE PRECISION  SCALAR
         DOUBLE PRECISION  XNORM

         DOUBLE PRECISION  NORM (0:L)

         DOUBLE PRECISION  BATCH (1:M,1:NXYZ)
C
C
C------------------------------------------------------------------------
C
C
C             ...generate monomial exponents lmn in appropriate
C                sequence and normalize integrals.
C
C
         N = 0

         DO 100 X = L,0,-1
            XNORM = NORM (X)
            YBEG = L - X
            DO 110 Y = YBEG,0,-1
               Z = YBEG - Y
               SCALAR = XNORM * NORM (Y) * NORM (Z)

               N = N + 1

               DO 120 I = 1,M
                  BATCH (I,N) = SCALAR * BATCH (I,N)
  120          CONTINUE

  110       CONTINUE
  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
