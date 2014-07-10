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
         SUBROUTINE  ERD__HRR_TRANSFORM
     +
     +                    ( M,
     +                      NROW,NXYZET,NXYZAB,
     +                      NXYZA,NXYZB,
     +                      LROW,ROW,
     +                      ROT,
     +                      X,
     +
     +                              Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__HRR_TRANSFORM
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This operation performs a HRR transformation on a
C                batch of contracted cartesian gaussian integrals:
C
C                                      nxyzet
C                     y (m,nxyzab)  =   sum   x (m,i) * rot (i,nxyzab)
C                                       i=1
C
C                where rot is the HRR transformation matrix. Due to
C                the very sparse nature of this matrix, only those
C                i indices in the summation are addressed which
C                correspond to nonzero HRR transformation matrix
C                elements.
C
C
C                  Input:
C
C                     M          =  # of elements not involved in the
C                                   transformation (invariant indices)
C                     NROW       =  maximum # of nonzero row elements
C                                   per column in transformation matrix
C                     NXYZET     =  dimension of the cartesian e0-part
C                     NXYZAB     =  dimension of the resulting cartesian
C                                   ab-part
C                     NXYZA      =  dimension of the cartesian part due
C                                   to the a-shell
C                     NXYZB      =  dimension of the cartesian part due
C                                   to the b-shell
C                     LROW (N)   =  # of nonzero entries in column N of
C                                   ROT and ROW matrix. N ranges from
C                                   1 to NXYZB
C                     ROW (I,N)  =  I-th nonzero row label of column N
C                                   in ROT matrix. N ranges from 1 to
C                                   NXYZAB
C                     ROT (I,N)  =  I-th nonzero HRR transformation
C                                   matrix element of column N. N ranges
C                                   from 1 to NXYZB
C                     X          =  batch of untransformed integrals
C                                   (m,e0)
C
C                  Output:
C
C                     Y          =  batch of HRR transformed integrals
C                                   (m,ab)
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT    NONE

         INTEGER     A,B,I,J,M,N
         INTEGER     MROW,NROW
         INTEGER     NLEFT,NSTEP
         INTEGER     NXYZET,NXYZAB
         INTEGER     NXYZA,NXYZB
         INTEGER     XCOL1,XCOL2,XCOL3,XCOL4,XCOL5,XCOL6,XCOL7,XCOL8

         INTEGER     LROW (1:NXYZB)
         INTEGER     ROW  (1:NROW,1:NXYZAB)

         DOUBLE PRECISION  ROT1,ROT2,ROT3,ROT4,ROT5,ROT6,ROT7,ROT8

         DOUBLE PRECISION  X (1:M,1:NXYZET)
         DOUBLE PRECISION  Y (1:M,1:NXYZAB)

         DOUBLE PRECISION  ROT (1:NROW,1:NXYZB)
C
C
C------------------------------------------------------------------------
C
C
C             ...perform the HRR transformation. One of the main
C                properties of this transformation is that the
C                last nonzero element of the HRR transformation
C                matrix is always equal to 1. Hence we can skip
C                the multiplication with that element.
C
C                Use basic row grouping of the transformation
C                to improve cache line reusing.
C
C
         N = 1

         DO 100 B = 1,NXYZB
            MROW = LROW (B)

            GOTO (1,2,3,4,5,6,7,8,9) MIN (MROW,9)
C
C
C             ...# of rows = 1
C
C
    1       DO 10 A = 1,NXYZA
               XCOL1 = ROW (1,N)
               DO J = 1,M
                  Y (J,N) = X (J,XCOL1)
               END DO
               N = N + 1
   10       CONTINUE

            GOTO 100
C
C
C             ...# of rows = 2
C
C
    2       DO 20 A = 1,NXYZA

               XCOL1 = ROW (1,N)
               XCOL2 = ROW (2,N)
               ROT1  = ROT (1,B)

               DO J = 1,M
                  Y (J,N) = ROT1 * X (J,XCOL1)
     +                           + X (J,XCOL2)
               END DO

               N = N + 1
   20       CONTINUE

            GOTO 100
C
C
C             ...# of rows = 3
C
C
    3       DO 30 A = 1,NXYZA

               XCOL1 = ROW (1,N)
               XCOL2 = ROW (2,N)
               XCOL3 = ROW (3,N)

               ROT1  = ROT (1,B)
               ROT2  = ROT (2,B)

               DO J = 1,M
                  Y (J,N) =   ROT1 * X (J,XCOL1)
     +                      + ROT2 * X (J,XCOL2)
     +                      +        X (J,XCOL3)
               END DO

               N = N + 1
   30       CONTINUE

            GOTO 100
C
C
C             ...# of rows = 4
C
C
    4       DO 40 A = 1,NXYZA

               XCOL1 = ROW (1,N)
               XCOL2 = ROW (2,N)
               XCOL3 = ROW (3,N)
               XCOL4 = ROW (4,N)

               ROT1  = ROT (1,B)
               ROT2  = ROT (2,B)
               ROT3  = ROT (3,B)

               DO J = 1,M
                  Y (J,N) =   ROT1 * X (J,XCOL1)
     +                      + ROT2 * X (J,XCOL2)
     +                      + ROT3 * X (J,XCOL3)
     +                      +        X (J,XCOL4)
               END DO

               N = N + 1
   40       CONTINUE

            GOTO 100
C
C
C             ...# of rows = 5
C
C
    5       DO 50 A = 1,NXYZA

               XCOL1 = ROW (1,N)
               XCOL2 = ROW (2,N)
               XCOL3 = ROW (3,N)
               XCOL4 = ROW (4,N)
               XCOL5 = ROW (5,N)

               ROT1  = ROT (1,B)
               ROT2  = ROT (2,B)
               ROT3  = ROT (3,B)
               ROT4  = ROT (4,B)

               DO J = 1,M
                  Y (J,N) =   ROT1 * X (J,XCOL1)
     +                      + ROT2 * X (J,XCOL2)
     +                      + ROT3 * X (J,XCOL3)
     +                      + ROT4 * X (J,XCOL4)
     +                      +        X (J,XCOL5)
               END DO

               N = N + 1
   50       CONTINUE

            GOTO 100
C
C
C             ...# of rows = 6
C
C
    6       DO 60 A = 1,NXYZA

               XCOL1 = ROW (1,N)
               XCOL2 = ROW (2,N)
               XCOL3 = ROW (3,N)
               XCOL4 = ROW (4,N)
               XCOL5 = ROW (5,N)
               XCOL6 = ROW (6,N)

               ROT1  = ROT (1,B)
               ROT2  = ROT (2,B)
               ROT3  = ROT (3,B)
               ROT4  = ROT (4,B)
               ROT5  = ROT (5,B)

               DO J = 1,M
                  Y (J,N) =   ROT1 * X (J,XCOL1)
     +                      + ROT2 * X (J,XCOL2)
     +                      + ROT3 * X (J,XCOL3)
     +                      + ROT4 * X (J,XCOL4)
     +                      + ROT5 * X (J,XCOL5)
     +                      +        X (J,XCOL6)
               END DO

               N = N + 1
   60       CONTINUE

            GOTO 100
C
C
C             ...# of rows = 7
C
C
    7       DO 70 A = 1,NXYZA

               XCOL1 = ROW (1,N)
               XCOL2 = ROW (2,N)
               XCOL3 = ROW (3,N)
               XCOL4 = ROW (4,N)
               XCOL5 = ROW (5,N)
               XCOL6 = ROW (6,N)
               XCOL7 = ROW (7,N)

               ROT1  = ROT (1,B)
               ROT2  = ROT (2,B)
               ROT3  = ROT (3,B)
               ROT4  = ROT (4,B)
               ROT5  = ROT (5,B)
               ROT6  = ROT (6,B)

               DO J = 1,M
                  Y (J,N) =   ROT1 * X (J,XCOL1)
     +                      + ROT2 * X (J,XCOL2)
     +                      + ROT3 * X (J,XCOL3)
     +                      + ROT4 * X (J,XCOL4)
     +                      + ROT5 * X (J,XCOL5)
     +                      + ROT6 * X (J,XCOL6)
     +                      +        X (J,XCOL7)
               END DO

               N = N + 1
   70       CONTINUE

            GOTO 100
C
C
C             ...# of rows = 8
C
C
    8       DO 80 A = 1,NXYZA

               XCOL1 = ROW (1,N)
               XCOL2 = ROW (2,N)
               XCOL3 = ROW (3,N)
               XCOL4 = ROW (4,N)
               XCOL5 = ROW (5,N)
               XCOL6 = ROW (6,N)
               XCOL7 = ROW (7,N)
               XCOL8 = ROW (8,N)

               ROT1  = ROT (1,B)
               ROT2  = ROT (2,B)
               ROT3  = ROT (3,B)
               ROT4  = ROT (4,B)
               ROT5  = ROT (5,B)
               ROT6  = ROT (6,B)
               ROT7  = ROT (7,B)

               DO J = 1,M
                  Y (J,N) =   ROT1 * X (J,XCOL1)
     +                      + ROT2 * X (J,XCOL2)
     +                      + ROT3 * X (J,XCOL3)
     +                      + ROT4 * X (J,XCOL4)
     +                      + ROT5 * X (J,XCOL5)
     +                      + ROT6 * X (J,XCOL6)
     +                      + ROT7 * X (J,XCOL7)
     +                      +        X (J,XCOL8)
               END DO

               N = N + 1
   80       CONTINUE

            GOTO 100
C
C
C             ...# of rows > 8. Perform transformations in bundles
C                of 8 rows.
C
C
    9       DO 90 A = 1,NXYZA

               XCOL1 = ROW (MROW-7,N)
               XCOL2 = ROW (MROW-6,N)
               XCOL3 = ROW (MROW-5,N)
               XCOL4 = ROW (MROW-4,N)
               XCOL5 = ROW (MROW-3,N)
               XCOL6 = ROW (MROW-2,N)
               XCOL7 = ROW (MROW-1,N)
               XCOL8 = ROW (MROW  ,N)

               ROT1  = ROT (MROW-7,B)
               ROT2  = ROT (MROW-6,B)
               ROT3  = ROT (MROW-5,B)
               ROT4  = ROT (MROW-4,B)
               ROT5  = ROT (MROW-3,B)
               ROT6  = ROT (MROW-2,B)
               ROT7  = ROT (MROW-1,B)

               DO J = 1,M
                  Y (J,N) =   ROT1 * X (J,XCOL1)
     +                      + ROT2 * X (J,XCOL2)
     +                      + ROT3 * X (J,XCOL3)
     +                      + ROT4 * X (J,XCOL4)
     +                      + ROT5 * X (J,XCOL5)
     +                      + ROT6 * X (J,XCOL6)
     +                      + ROT7 * X (J,XCOL7)
     +                      +        X (J,XCOL8)
               END DO

               NLEFT = MROW - 8
               NSTEP = NLEFT / 8

               DO 92 I = 1,NSTEP

                  XCOL1 = ROW (NLEFT-7,N)
                  XCOL2 = ROW (NLEFT-6,N)
                  XCOL3 = ROW (NLEFT-5,N)
                  XCOL4 = ROW (NLEFT-4,N)
                  XCOL5 = ROW (NLEFT-3,N)
                  XCOL6 = ROW (NLEFT-2,N)
                  XCOL7 = ROW (NLEFT-1,N)
                  XCOL8 = ROW (NLEFT  ,N)

                  ROT1  = ROT (NLEFT-7,B)
                  ROT2  = ROT (NLEFT-6,B)
                  ROT3  = ROT (NLEFT-5,B)
                  ROT4  = ROT (NLEFT-4,B)
                  ROT5  = ROT (NLEFT-3,B)
                  ROT6  = ROT (NLEFT-2,B)
                  ROT7  = ROT (NLEFT-1,B)
                  ROT8  = ROT (NLEFT  ,B)

                  DO J = 1,M
                     Y (J,N) = Y (J,N) + ROT1 * X (J,XCOL1)
     +                                 + ROT2 * X (J,XCOL2)
     +                                 + ROT3 * X (J,XCOL3)
     +                                 + ROT4 * X (J,XCOL4)
     +                                 + ROT5 * X (J,XCOL5)
     +                                 + ROT6 * X (J,XCOL6)
     +                                 + ROT7 * X (J,XCOL7)
     +                                 + ROT8 * X (J,XCOL8)
                  END DO
                  NLEFT = NLEFT - 8
   92          CONTINUE

               DO 94 I = NLEFT,1,-1
                  XCOL1 = ROW (I,N)
                  ROT1  = ROT (I,B)
                  DO J = 1,M
                     Y (J,N) = Y (J,N) + ROT1 * X (J,XCOL1)
                  END DO
   94          CONTINUE

               N = N + 1
   90       CONTINUE
C
C
C             ...next b-shell monomial.
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
