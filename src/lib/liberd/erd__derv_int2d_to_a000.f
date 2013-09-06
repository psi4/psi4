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
         SUBROUTINE  ERD__DERV_INT2D_TO_A000
     +
     +                    ( SHELLA,
     +                      NGQP,NEXQ,NGQEXQ,
     +                      NXYZA,
     +                      INT2DX,INT2DY,INT2DZ,
     +                      DIFFY,DIFFZ,
     +                      TEMP1,TEMP2,
     +                      SCALE,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__DERV_INT2D_TO_A000
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of cartesian
C                derivative eris:
C
C                     [A0|00] , [0A|00] , [00|A0] or [00|0A]
C
C                adding up the contributions from all the respective
C                derivative 2D integrals:
C
C                           A000 , 0A00 , 00A0 or 000A
C
C                Simplified version of the general ABCD routine to
C                reduce loop overheads for those cases where there are
C                at least three s-shells. For comments and details see
C                the general ABCD routine.
C
C
C                  Input:
C
C                    SHELLA      =  shell type for contraction shell A
C                    NGQP        =  # of gaussian quadrature points
C                                   (roots)
C                    NEXQ        =  current # of exponent quadruplets
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    NXYZA       =  # of cartesian monomials for shell A
C                    INT2Dx      =  all current 2D A000/0A00/00A0/000A
C                                   derivative integrals for each
C                                   cartesian component (x = X,Y,Z)
C                    DIFFx       =  is true, if differentiation was
C                                   performed along the x=Y,Z direction
C                    TEMP1(2)    =  scratch arrays holding intermediate
C                                   2D A000/0A00/00A0/000A derivative
C                                   integral products
C                    SCALE       =  the NGQEXQ scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive cartesian
C                                   [A0|00]/[0A|00]/[00|A0]/[00|0A]
C                                   derivative integrals corresponding
C                                   to all current exponent quadruplets
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

         LOGICAL     DIFFY,DIFFZ

         INTEGER     I,M,N,R
         INTEGER     NGQP,NEXQ,NGQEXQ
         INTEGER     NXYZA
         INTEGER     SHELLA
         INTEGER     XA,YA,ZA
         INTEGER     YAMAX

         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  SCALE (1:NGQEXQ)
         DOUBLE PRECISION  TEMP1 (1:NGQEXQ)
         DOUBLE PRECISION  TEMP2 (1:NGQEXQ)

         DOUBLE PRECISION  BATCH  (1:NEXQ,1:NXYZA)

         DOUBLE PRECISION  INT2DX (1:NGQEXQ,0:SHELLA)
         DOUBLE PRECISION  INT2DY (1:NGQEXQ,0:SHELLA)
         DOUBLE PRECISION  INT2DZ (1:NGQEXQ,0:SHELLA)

         PARAMETER  (ZERO  =  0.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...jump according to number of roots.
C
C
         GOTO  (1,2,3,4,5,6,7,8,9,10)  MIN (NGQP,10)
C
C
C                       ********************
C                       *  # of roots = 1  *
C                       ********************
C
C
    1    I = 0
         DO 100 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 110 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO M = 1,NEXQ
                      TEMP2 (M) = TEMP1 (M)
                   END DO
               ELSE
                   DO M = 1,NEXQ
                      TEMP2 (M) = TEMP1 (M) * INT2DY (M,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   DO M = 1,NEXQ
                      BATCH (M,I) = TEMP2 (M)
                   END DO
               ELSE
                   DO M = 1,NEXQ
                      BATCH (M,I) = TEMP2 (M) * INT2DZ (M,ZA)
                   END DO
               END IF

  110       CONTINUE
  100    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 2  *
C                       ********************
C
C
    2    I = 0
         DO 200 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 210 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) =   TEMP2 (R)
     +                              + TEMP2 (R+1)
                      R = R + 2
                   END DO
               ELSE
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) =   TEMP2 (R)
     +                              * INT2DZ (R,ZA)
     +                              + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA)
                      R = R + 2
                   END DO
               END IF

  210       CONTINUE
  200    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 3  *
C                       ********************
C
C
    3    I = 0
         DO 300 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 310 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) =   TEMP2 (R)
     +                              + TEMP2 (R+1)
     +                              + TEMP2 (R+2)
                      R = R + 3
                   END DO
               ELSE
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) =   TEMP2 (R)
     +                              * INT2DZ (R,ZA)
     +                              + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA)
     +                              + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA)
                      R = R + 3
                   END DO
               END IF

  310       CONTINUE
  300    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 4  *
C                       ********************
C
C
    4    I = 0
         DO 400 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 410 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) =   TEMP2 (R)
     +                              + TEMP2 (R+1)
     +                              + TEMP2 (R+2)
     +                              + TEMP2 (R+3)
                      R = R + 4
                   END DO
               ELSE
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) =   TEMP2 (R)
     +                              * INT2DZ (R,ZA)
     +                              + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA)
     +                              + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA)
     +                              + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA)
                      R = R + 4
                   END DO
               END IF

  410       CONTINUE
  400    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 5  *
C                       ********************
C
C
    5    I = 0
         DO 500 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 510 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) =   TEMP2 (R)
     +                              + TEMP2 (R+1)
     +                              + TEMP2 (R+2)
     +                              + TEMP2 (R+3)
     +                              + TEMP2 (R+4)
                      R = R + 5
                   END DO
               ELSE
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) =   TEMP2 (R)
     +                              * INT2DZ (R,ZA)
     +                              + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA)
     +                              + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA)
     +                              + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA)
     +                              + TEMP2 (R+4)
     +                              * INT2DZ (R+4,ZA)
                      R = R + 5
                   END DO
               END IF

  510       CONTINUE
  500    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 6  *
C                       ********************
C
C
    6    I = 0
         DO 600 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 610 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) =   TEMP2 (R)
     +                              + TEMP2 (R+1)
     +                              + TEMP2 (R+2)
     +                              + TEMP2 (R+3)
     +                              + TEMP2 (R+4)
     +                              + TEMP2 (R+5)
                      R = R + 6
                   END DO
               ELSE
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) =   TEMP2 (R)
     +                              * INT2DZ (R,ZA)
     +                              + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA)
     +                              + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA)
     +                              + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA)
     +                              + TEMP2 (R+4)
     +                              * INT2DZ (R+4,ZA)
     +                              + TEMP2 (R+5)
     +                              * INT2DZ (R+5,ZA)
                      R = R + 6
                   END DO
               END IF

  610       CONTINUE
  600    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 7  *
C                       ********************
C
C
    7    I = 0
         DO 700 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 710 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) =   TEMP2 (R)
     +                              + TEMP2 (R+1)
     +                              + TEMP2 (R+2)
     +                              + TEMP2 (R+3)
     +                              + TEMP2 (R+4)
     +                              + TEMP2 (R+5)
     +                              + TEMP2 (R+6)
                      R = R + 7
                   END DO
               ELSE
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) =   TEMP2 (R)
     +                              * INT2DZ (R,ZA)
     +                              + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA)
     +                              + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA)
     +                              + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA)
     +                              + TEMP2 (R+4)
     +                              * INT2DZ (R+4,ZA)
     +                              + TEMP2 (R+5)
     +                              * INT2DZ (R+5,ZA)
     +                              + TEMP2 (R+6)
     +                              * INT2DZ (R+6,ZA)
                      R = R + 7
                   END DO
               END IF

  710       CONTINUE
  700    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 8  *
C                       ********************
C
C
    8    I = 0
         DO 800 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 810 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) =   TEMP2 (R)
     +                              + TEMP2 (R+1)
     +                              + TEMP2 (R+2)
     +                              + TEMP2 (R+3)
     +                              + TEMP2 (R+4)
     +                              + TEMP2 (R+5)
     +                              + TEMP2 (R+6)
     +                              + TEMP2 (R+7)
                      R = R + 8
                   END DO
               ELSE
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) =   TEMP2 (R)
     +                              * INT2DZ (R,ZA)
     +                              + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA)
     +                              + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA)
     +                              + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA)
     +                              + TEMP2 (R+4)
     +                              * INT2DZ (R+4,ZA)
     +                              + TEMP2 (R+5)
     +                              * INT2DZ (R+5,ZA)
     +                              + TEMP2 (R+6)
     +                              * INT2DZ (R+6,ZA)
     +                              + TEMP2 (R+7)
     +                              * INT2DZ (R+7,ZA)
                      R = R + 8
                   END DO
               END IF

  810       CONTINUE
  800    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 9  *
C                       ********************
C
C
    9    I = 0
         DO 900 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO

            DO 910 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF

               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) =   TEMP2 (R)
     +                              + TEMP2 (R+1)
     +                              + TEMP2 (R+2)
     +                              + TEMP2 (R+3)
     +                              + TEMP2 (R+4)
     +                              + TEMP2 (R+5)
     +                              + TEMP2 (R+6)
     +                              + TEMP2 (R+7)
     +                              + TEMP2 (R+8)
                      R = R + 9
                   END DO
               ELSE
                   R = 1
                   DO M = 1,NEXQ
                      BATCH (M,I) =   TEMP2 (R)
     +                              * INT2DZ (R,ZA)
     +                              + TEMP2 (R+1)
     +                              * INT2DZ (R+1,ZA)
     +                              + TEMP2 (R+2)
     +                              * INT2DZ (R+2,ZA)
     +                              + TEMP2 (R+3)
     +                              * INT2DZ (R+3,ZA)
     +                              + TEMP2 (R+4)
     +                              * INT2DZ (R+4,ZA)
     +                              + TEMP2 (R+5)
     +                              * INT2DZ (R+5,ZA)
     +                              + TEMP2 (R+6)
     +                              * INT2DZ (R+6,ZA)
     +                              + TEMP2 (R+7)
     +                              * INT2DZ (R+7,ZA)
     +                              + TEMP2 (R+8)
     +                              * INT2DZ (R+8,ZA)
                      R = R + 9
                   END DO
               END IF

  910       CONTINUE
  900    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots > 9  *
C                       ********************
C
C             ...outer loop over x-contribution. No skipping of
C                x-contribution of 0-type can be done here,
C                since the 2DX integrals carry the Rys weight!
C
C
   10    I = 0
         DO 1000 XA = SHELLA,0,-1
            YAMAX = SHELLA - XA

            DO M = 1,NGQEXQ
               TEMP1 (M) = SCALE (M) * INT2DX (M,XA)
            END DO
C
C
C             ...inner loop over y-contribution. Skip the multiplication
C                of y-contributions, if no y-coordinate derivative
C                was formed and we have a 0-contribution, as then the
C                2DY integrals are equal to 1.
C
C
            DO 1010 YA = YAMAX,0,-1
               I = I + 1
               ZA = YAMAX - YA

               IF (.NOT.DIFFY .AND. YA.EQ.0) THEN
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N)
                   END DO
               ELSE
                   DO N = 1,NGQEXQ
                      TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA)
                   END DO
               END IF
C
C
C             ...skip multiplication of z-contributions, if we
C                have a 0-type and no derivations were performed
C                on the z-coordinate, as then the 2DZ integrals
C                are equal to 1. All info concerning all three x-,
C                y- and z-contributions have been collected for all
C                exponent quadruplets at once. Sum up the 2D X,Y,Z
C                integral products to the appropriate places of the
C                batch.
C
C
               IF (.NOT.DIFFZ .AND. ZA.EQ.0) THEN
                   R = 0
                   DO M = 1,NEXQ
                      SUM = ZERO
                      DO N = 1,NGQP
                         SUM = SUM + TEMP2 (R+N)
                      END DO
                      R = R + NGQP
                      BATCH (M,I) = SUM
                   END DO
               ELSE
                   R = 0
                   DO M = 1,NEXQ
                      SUM = ZERO
                      DO N = 1,NGQP
                         SUM = SUM + TEMP2 (R+N) * INT2DZ (R+N,ZA)
                      END DO
                      R = R + NGQP
                      BATCH (M,I) = SUM
                   END DO
               END IF
C
C
C             ...next z- and y-contributions.
C
C
 1010       CONTINUE
C
C
C             ...next x-contribution.
C
C
 1000    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
