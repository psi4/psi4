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
         SUBROUTINE  ERD__SPHERICAL_TRANSFORM
     +
     +                    ( M,
     +                      NROW,NXYZ,NRY,
     +                      LROW,ROW,
     +                      ROT,
     +                      X,
     +
     +                              Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__SPHERICAL_TRANSFORM
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation performs a simple cartesian to spherical
C                transformation step on a batch of contracted cartesian
C                gaussian integrals:
C
C                         y (m,r) = sum  mat (i,r) * x (m,i)
C                                    i
C
C                where i is an element of the cartesian basis, r an
C                element of the spherical basis, mat is the cartesian
C                -> spherical transformation matrix and m are the bra
C                elements not involved. Ordering of the cartesian
C                monomial basis is assumed to be:
C
C                         e f g     a b c
C                        X Y Z  >  X Y Z   , if e > a
C                                            for e=a if f > b
C                                            for e=a and f=b if g > c
C
C
C                  Input:
C
C                    M           =  # of elements not involved in the
C                                   transformation (invariant indices)
C                    NROW        =  maximum # of nonzero rows in the
C                                   transformation matrix
C                    NXYZ        =  # of cartesian monomials xyz for
C                                   the shell to be transformed
C                    NRY         =  # of spherical functions ry for
C                                   the shell to be transformed
C                    LROW (R)    =  # of xyz-monomials contributing
C                                   to the R-th ry-component
C                    ROW (I,R)   =  I-th xyz-monomial row index
C                                   containing nonzero contribution to
C                                   the R-th ry-component
C                    ROT (I,R)   =  I-th nonzero xyz-monomial to R-th
C                                   ry-component transformation matrix
C                                   element
C                    X           =  input batch of cartesian integrals
C
C                  Output:
C
C                    Y           =  output batch of spherical integrals
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

         INTEGER     I,M,N,R
         INTEGER     MROW,NROW
         INTEGER     NBASE,NLEFT,NREST,NSTEP
         INTEGER     NRY
         INTEGER     NXYZ
         INTEGER     XCOL1,XCOL2,XCOL3,XCOL4,XCOL5,XCOL6,XCOL7,XCOL8

         INTEGER     LROW (1:NRY)
         INTEGER     ROW  (1:NROW,1:NRY)

         DOUBLE PRECISION  ROT1,ROT2,ROT3,ROT4,ROT5,ROT6,ROT7,ROT8

         DOUBLE PRECISION  X (1:M,1:NXYZ)
         DOUBLE PRECISION  Y (1:M,1:NRY)

         DOUBLE PRECISION  ROT (1:NROW,1:NRY)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the cartesian -> spherical transformation.
C                Use basic row grouping of the transformation
C                to improve cache line reusing.
C
C
         DO 100 R = 1,NRY
            MROW = LROW (R)

            GOTO (1,2,3,4,5,6,7,8,9) MIN (MROW,9)
C
C
C             ...# of rows = 1
C
C
    1       XCOL1 = ROW (1,R)
            ROT1  = ROT (1,R)

            DO N = 1,M
               Y (N,R) = ROT1 * X (N,XCOL1)
            END DO

            GOTO 100
C
C
C             ...# of rows = 2
C
C
    2       XCOL1 = ROW (1,R)
            XCOL2 = ROW (2,R)
            ROT1  = ROT (1,R)
            ROT2  = ROT (2,R)

            DO N = 1,M
               Y (N,R) =   ROT1 * X (N,XCOL1)
     +                   + ROT2 * X (N,XCOL2)
            END DO

            GOTO 100
C
C
C             ...# of rows = 3
C
C
    3       XCOL1 = ROW (1,R)
            XCOL2 = ROW (2,R)
            XCOL3 = ROW (3,R)

            ROT1  = ROT (1,R)
            ROT2  = ROT (2,R)
            ROT3  = ROT (3,R)

            DO N = 1,M
               Y (N,R) =   ROT1 * X (N,XCOL1)
     +                   + ROT2 * X (N,XCOL2)
     +                   + ROT3 * X (N,XCOL3)
            END DO

            GOTO 100
C
C
C             ...# of rows = 4
C
C
    4       XCOL1 = ROW (1,R)
            XCOL2 = ROW (2,R)
            XCOL3 = ROW (3,R)
            XCOL4 = ROW (4,R)

            ROT1  = ROT (1,R)
            ROT2  = ROT (2,R)
            ROT3  = ROT (3,R)
            ROT4  = ROT (4,R)

            DO N = 1,M
               Y (N,R) =   ROT1 * X (N,XCOL1)
     +                   + ROT2 * X (N,XCOL2)
     +                   + ROT3 * X (N,XCOL3)
     +                   + ROT4 * X (N,XCOL4)
            END DO

            GOTO 100
C
C
C             ...# of rows = 5
C
C
    5       XCOL1 = ROW (1,R)
            XCOL2 = ROW (2,R)
            XCOL3 = ROW (3,R)
            XCOL4 = ROW (4,R)
            XCOL5 = ROW (5,R)

            ROT1  = ROT (1,R)
            ROT2  = ROT (2,R)
            ROT3  = ROT (3,R)
            ROT4  = ROT (4,R)
            ROT5  = ROT (5,R)

            DO N = 1,M
               Y (N,R) =   ROT1 * X (N,XCOL1)
     +                   + ROT2 * X (N,XCOL2)
     +                   + ROT3 * X (N,XCOL3)
     +                   + ROT4 * X (N,XCOL4)
     +                   + ROT5 * X (N,XCOL5)
            END DO

            GOTO 100
C
C
C             ...# of rows = 6
C
C
    6       XCOL1 = ROW (1,R)
            XCOL2 = ROW (2,R)
            XCOL3 = ROW (3,R)
            XCOL4 = ROW (4,R)
            XCOL5 = ROW (5,R)
            XCOL6 = ROW (6,R)

            ROT1  = ROT (1,R)
            ROT2  = ROT (2,R)
            ROT3  = ROT (3,R)
            ROT4  = ROT (4,R)
            ROT5  = ROT (5,R)
            ROT6  = ROT (6,R)

            DO N = 1,M
               Y (N,R) =   ROT1 * X (N,XCOL1)
     +                   + ROT2 * X (N,XCOL2)
     +                   + ROT3 * X (N,XCOL3)
     +                   + ROT4 * X (N,XCOL4)
     +                   + ROT5 * X (N,XCOL5)
     +                   + ROT6 * X (N,XCOL6)
            END DO

            GOTO 100
C
C
C             ...# of rows = 7
C
C
    7       XCOL1 = ROW (1,R)
            XCOL2 = ROW (2,R)
            XCOL3 = ROW (3,R)
            XCOL4 = ROW (4,R)
            XCOL5 = ROW (5,R)
            XCOL6 = ROW (6,R)
            XCOL7 = ROW (7,R)

            ROT1  = ROT (1,R)
            ROT2  = ROT (2,R)
            ROT3  = ROT (3,R)
            ROT4  = ROT (4,R)
            ROT5  = ROT (5,R)
            ROT6  = ROT (6,R)
            ROT7  = ROT (7,R)

            DO N = 1,M
               Y (N,R) =   ROT1 * X (N,XCOL1)
     +                   + ROT2 * X (N,XCOL2)
     +                   + ROT3 * X (N,XCOL3)
     +                   + ROT4 * X (N,XCOL4)
     +                   + ROT5 * X (N,XCOL5)
     +                   + ROT6 * X (N,XCOL6)
     +                   + ROT7 * X (N,XCOL7)
            END DO

            GOTO 100
C
C
C             ...# of rows = 8
C
C
    8       XCOL1 = ROW (1,R)
            XCOL2 = ROW (2,R)
            XCOL3 = ROW (3,R)
            XCOL4 = ROW (4,R)
            XCOL5 = ROW (5,R)
            XCOL6 = ROW (6,R)
            XCOL7 = ROW (7,R)
            XCOL8 = ROW (8,R)

            ROT1  = ROT (1,R)
            ROT2  = ROT (2,R)
            ROT3  = ROT (3,R)
            ROT4  = ROT (4,R)
            ROT5  = ROT (5,R)
            ROT6  = ROT (6,R)
            ROT7  = ROT (7,R)
            ROT8  = ROT (8,R)

            DO N = 1,M
               Y (N,R) =   ROT1 * X (N,XCOL1)
     +                   + ROT2 * X (N,XCOL2)
     +                   + ROT3 * X (N,XCOL3)
     +                   + ROT4 * X (N,XCOL4)
     +                   + ROT5 * X (N,XCOL5)
     +                   + ROT6 * X (N,XCOL6)
     +                   + ROT7 * X (N,XCOL7)
     +                   + ROT8 * X (N,XCOL8)
            END DO

            GOTO 100
C
C
C             ...# of rows > 8. Perform transformations in bundles
C                of 8 rows.
C
C
    9       XCOL1 = ROW (1,R)
            XCOL2 = ROW (2,R)
            XCOL3 = ROW (3,R)
            XCOL4 = ROW (4,R)
            XCOL5 = ROW (5,R)
            XCOL6 = ROW (6,R)
            XCOL7 = ROW (7,R)
            XCOL8 = ROW (8,R)

            ROT1  = ROT (1,R)
            ROT2  = ROT (2,R)
            ROT3  = ROT (3,R)
            ROT4  = ROT (4,R)
            ROT5  = ROT (5,R)
            ROT6  = ROT (6,R)
            ROT7  = ROT (7,R)
            ROT8  = ROT (8,R)

            DO N = 1,M
               Y (N,R) =   ROT1 * X (N,XCOL1)
     +                   + ROT2 * X (N,XCOL2)
     +                   + ROT3 * X (N,XCOL3)
     +                   + ROT4 * X (N,XCOL4)
     +                   + ROT5 * X (N,XCOL5)
     +                   + ROT6 * X (N,XCOL6)
     +                   + ROT7 * X (N,XCOL7)
     +                   + ROT8 * X (N,XCOL8)
            END DO

            NBASE = 8
            NLEFT = MROW - 8
            NSTEP = NLEFT / 8
            NREST = MOD (NLEFT,8)

            DO 90 I = 1,NSTEP

               XCOL1 = ROW (NBASE+1,R)
               XCOL2 = ROW (NBASE+2,R)
               XCOL3 = ROW (NBASE+3,R)
               XCOL4 = ROW (NBASE+4,R)
               XCOL5 = ROW (NBASE+5,R)
               XCOL6 = ROW (NBASE+6,R)
               XCOL7 = ROW (NBASE+7,R)
               XCOL8 = ROW (NBASE+8,R)

               ROT1  = ROT (NBASE+1,R)
               ROT2  = ROT (NBASE+2,R)
               ROT3  = ROT (NBASE+3,R)
               ROT4  = ROT (NBASE+4,R)
               ROT5  = ROT (NBASE+5,R)
               ROT6  = ROT (NBASE+6,R)
               ROT7  = ROT (NBASE+7,R)
               ROT8  = ROT (NBASE+8,R)

               DO N = 1,M
                  Y (N,R) = Y (N,R) + ROT1 * X (N,XCOL1)
     +                              + ROT2 * X (N,XCOL2)
     +                              + ROT3 * X (N,XCOL3)
     +                              + ROT4 * X (N,XCOL4)
     +                              + ROT5 * X (N,XCOL5)
     +                              + ROT6 * X (N,XCOL6)
     +                              + ROT7 * X (N,XCOL7)
     +                              + ROT8 * X (N,XCOL8)
               END DO

               NBASE = NBASE + 8

   90       CONTINUE

            DO 92 I = 1,NREST
               XCOL1 = ROW (NBASE+I,R)
               ROT1  = ROT (NBASE+I,R)
               DO N = 1,M
                  Y (N,R) = Y (N,R) + ROT1 * X (N,XCOL1)
               END DO
   92       CONTINUE
C
C
C             ...next spherical index.
C
C
  100    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
