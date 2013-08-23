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
         SUBROUTINE  ERD__DERV_INT2D_TO_0000
     +
     +                    ( NGQP,NEXQ,NGQEXQ,
     +                      INT2DX,INT2DY,INT2DZ,
     +                      DIFFY,DIFFZ,
     +                      TEMP1,
     +                      SCALE,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__DERV_INT2D_TO_0000
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of
C                derivative eris [00|00], adding up the contributions
C                from all the derivative 2D 0000 integrals.
C
C                Special routine for the [00|00] eri cases.
C
C
C                  Input:
C
C                    NGQP        =  # of gaussian quadrature points
C                                   (roots)
C                    NEXQ        =  current # of exponent quadruplets
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    INT2Dx      =  all current 2D 0000 derivative
C                                   integrals for each cartesian
C                                   component (x = X,Y,Z)
C                    DIFFx       =  is true, if differentiation was
C                                   performed along the x=Y,Z direction
C                    TEMP1       =  scratch array holding intermediate
C                                   2D 0000 derivative integral products
C                    SCALE       =  the NGQEXQ scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive [00|00]
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

         INTEGER     M,N,R
         INTEGER     NGQP,NEXQ,NGQEXQ

         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  SCALE (1:NGQEXQ)
         DOUBLE PRECISION  TEMP1 (1:NGQEXQ)

         DOUBLE PRECISION  BATCH (1:NEXQ)

         DOUBLE PRECISION  INT2DX (1:NGQEXQ)
         DOUBLE PRECISION  INT2DY (1:NGQEXQ)
         DOUBLE PRECISION  INT2DZ (1:NGQEXQ)

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
    1    DO M = 1,NEXQ
            BATCH (M) = INT2DX (M) * SCALE (M)
         END DO

         IF (DIFFY) THEN
             DO M = 1,NEXQ
                BATCH (M) = BATCH (M) * INT2DY (M)
             END DO
         END IF

         IF (DIFFZ) THEN
             DO M = 1,NEXQ
                BATCH (M) = BATCH (M) * INT2DZ (M)
             END DO
         END IF

         RETURN
C
C
C                       ********************
C                       *  # of roots = 2  *
C                       ********************
C
C
    2    DO M = 1,NGQEXQ
            TEMP1 (M) = SCALE (M) * INT2DX (M)
         END DO

         IF (.NOT.DIFFY .AND. .NOT.DIFFZ) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)
     +                      + TEMP1 (R+1)
                R = R + 2
             END DO
         ELSE IF (.NOT.DIFFY) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DZ (R)
     +                      + TEMP1 (R+1) * INT2DZ (R+1)
                R = R + 2
             END DO
         ELSE IF (.NOT.DIFFZ) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DY (R)
     +                      + TEMP1 (R+1) * INT2DY (R+1)
                R = R + 2
             END DO
         ELSE
             DO N = 1,NGQEXQ
                TEMP1 (N) = TEMP1 (N) * INT2DY (N)
             END DO
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DZ (R)
     +                      + TEMP1 (R+1) * INT2DZ (R+1)
                R = R + 2
             END DO
         END IF

         RETURN
C
C
C                       ********************
C                       *  # of roots = 3  *
C                       ********************
C
C
    3    DO M = 1,NGQEXQ
            TEMP1 (M) = SCALE (M) * INT2DX (M)
         END DO

         IF (.NOT.DIFFY .AND. .NOT.DIFFZ) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)
     +                      + TEMP1 (R+1)
     +                      + TEMP1 (R+2)
                R = R + 3
             END DO
         ELSE IF (.NOT.DIFFY) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DZ (R)
     +                      + TEMP1 (R+1) * INT2DZ (R+1)
     +                      + TEMP1 (R+2) * INT2DZ (R+2)
                R = R + 3
             END DO
         ELSE IF (.NOT.DIFFZ) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DY (R)
     +                      + TEMP1 (R+1) * INT2DY (R+1)
     +                      + TEMP1 (R+2) * INT2DY (R+2)
                R = R + 3
             END DO
         ELSE
             DO N = 1,NGQEXQ
                TEMP1 (N) = TEMP1 (N) * INT2DY (N)
             END DO
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DZ (R)
     +                      + TEMP1 (R+1) * INT2DZ (R+1)
     +                      + TEMP1 (R+2) * INT2DZ (R+2)
                R = R + 3
             END DO
         END IF

         RETURN
C
C
C                       ********************
C                       *  # of roots = 4  *
C                       ********************
C
C
    4    DO M = 1,NGQEXQ
            TEMP1 (M) = SCALE (M) * INT2DX (M)
         END DO

         IF (.NOT.DIFFY .AND. .NOT.DIFFZ) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)
     +                      + TEMP1 (R+1)
     +                      + TEMP1 (R+2)
     +                      + TEMP1 (R+3)
                R = R + 4
             END DO
         ELSE IF (.NOT.DIFFY) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DZ (R)
     +                      + TEMP1 (R+1) * INT2DZ (R+1)
     +                      + TEMP1 (R+2) * INT2DZ (R+2)
     +                      + TEMP1 (R+3) * INT2DZ (R+3)
                R = R + 4
             END DO
         ELSE IF (.NOT.DIFFZ) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DY (R)
     +                      + TEMP1 (R+1) * INT2DY (R+1)
     +                      + TEMP1 (R+2) * INT2DY (R+2)
     +                      + TEMP1 (R+3) * INT2DY (R+3)
                R = R + 4
             END DO
         ELSE
             DO N = 1,NGQEXQ
                TEMP1 (N) = TEMP1 (N) * INT2DY (N)
             END DO
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DZ (R)
     +                      + TEMP1 (R+1) * INT2DZ (R+1)
     +                      + TEMP1 (R+2) * INT2DZ (R+2)
     +                      + TEMP1 (R+3) * INT2DZ (R+3)
                R = R + 4
             END DO
         END IF

         RETURN
C
C
C                       ********************
C                       *  # of roots = 5  *
C                       ********************
C
C
    5    DO M = 1,NGQEXQ
            TEMP1 (M) = SCALE (M) * INT2DX (M)
         END DO

         IF (.NOT.DIFFY .AND. .NOT.DIFFZ) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)
     +                      + TEMP1 (R+1)
     +                      + TEMP1 (R+2)
     +                      + TEMP1 (R+3)
     +                      + TEMP1 (R+4)
                R = R + 5
             END DO
         ELSE IF (.NOT.DIFFY) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DZ (R)
     +                      + TEMP1 (R+1) * INT2DZ (R+1)
     +                      + TEMP1 (R+2) * INT2DZ (R+2)
     +                      + TEMP1 (R+3) * INT2DZ (R+3)
     +                      + TEMP1 (R+4) * INT2DZ (R+4)
                R = R + 5
             END DO
         ELSE IF (.NOT.DIFFZ) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DY (R)
     +                      + TEMP1 (R+1) * INT2DY (R+1)
     +                      + TEMP1 (R+2) * INT2DY (R+2)
     +                      + TEMP1 (R+3) * INT2DY (R+3)
     +                      + TEMP1 (R+4) * INT2DY (R+4)
                R = R + 5
             END DO
         ELSE
             DO N = 1,NGQEXQ
                TEMP1 (N) = TEMP1 (N) * INT2DY (N)
             END DO
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DZ (R)
     +                      + TEMP1 (R+1) * INT2DZ (R+1)
     +                      + TEMP1 (R+2) * INT2DZ (R+2)
     +                      + TEMP1 (R+3) * INT2DZ (R+3)
     +                      + TEMP1 (R+4) * INT2DZ (R+4)
                R = R + 5
             END DO
         END IF

         RETURN
C
C
C                       ********************
C                       *  # of roots = 6  *
C                       ********************
C
C
    6    DO M = 1,NGQEXQ
            TEMP1 (M) = SCALE (M) * INT2DX (M)
         END DO

         IF (.NOT.DIFFY .AND. .NOT.DIFFZ) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)
     +                      + TEMP1 (R+1)
     +                      + TEMP1 (R+2)
     +                      + TEMP1 (R+3)
     +                      + TEMP1 (R+4)
     +                      + TEMP1 (R+5)
                R = R + 6
             END DO
         ELSE IF (.NOT.DIFFY) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DZ (R)
     +                      + TEMP1 (R+1) * INT2DZ (R+1)
     +                      + TEMP1 (R+2) * INT2DZ (R+2)
     +                      + TEMP1 (R+3) * INT2DZ (R+3)
     +                      + TEMP1 (R+4) * INT2DZ (R+4)
     +                      + TEMP1 (R+5) * INT2DZ (R+5)
                R = R + 6
             END DO
         ELSE IF (.NOT.DIFFZ) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DY (R)
     +                      + TEMP1 (R+1) * INT2DY (R+1)
     +                      + TEMP1 (R+2) * INT2DY (R+2)
     +                      + TEMP1 (R+3) * INT2DY (R+3)
     +                      + TEMP1 (R+4) * INT2DY (R+4)
     +                      + TEMP1 (R+5) * INT2DY (R+5)
                R = R + 6
             END DO
         ELSE
             DO N = 1,NGQEXQ
                TEMP1 (N) = TEMP1 (N) * INT2DY (N)
             END DO
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DZ (R)
     +                      + TEMP1 (R+1) * INT2DZ (R+1)
     +                      + TEMP1 (R+2) * INT2DZ (R+2)
     +                      + TEMP1 (R+3) * INT2DZ (R+3)
     +                      + TEMP1 (R+4) * INT2DZ (R+4)
     +                      + TEMP1 (R+5) * INT2DZ (R+5)
                R = R + 6
             END DO
         END IF

         RETURN
C
C
C                       ********************
C                       *  # of roots = 7  *
C                       ********************
C
C
    7    DO M = 1,NGQEXQ
            TEMP1 (M) = SCALE (M) * INT2DX (M)
         END DO

         IF (.NOT.DIFFY .AND. .NOT.DIFFZ) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)
     +                      + TEMP1 (R+1)
     +                      + TEMP1 (R+2)
     +                      + TEMP1 (R+3)
     +                      + TEMP1 (R+4)
     +                      + TEMP1 (R+5)
     +                      + TEMP1 (R+6)
                R = R + 7
             END DO
         ELSE IF (.NOT.DIFFY) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DZ (R)
     +                      + TEMP1 (R+1) * INT2DZ (R+1)
     +                      + TEMP1 (R+2) * INT2DZ (R+2)
     +                      + TEMP1 (R+3) * INT2DZ (R+3)
     +                      + TEMP1 (R+4) * INT2DZ (R+4)
     +                      + TEMP1 (R+5) * INT2DZ (R+5)
     +                      + TEMP1 (R+6) * INT2DZ (R+6)
                R = R + 7
             END DO
         ELSE IF (.NOT.DIFFZ) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DY (R)
     +                      + TEMP1 (R+1) * INT2DY (R+1)
     +                      + TEMP1 (R+2) * INT2DY (R+2)
     +                      + TEMP1 (R+3) * INT2DY (R+3)
     +                      + TEMP1 (R+4) * INT2DY (R+4)
     +                      + TEMP1 (R+5) * INT2DY (R+5)
     +                      + TEMP1 (R+6) * INT2DY (R+6)
                R = R + 7
             END DO
         ELSE
             DO N = 1,NGQEXQ
                TEMP1 (N) = TEMP1 (N) * INT2DY (N)
             END DO
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DZ (R)
     +                      + TEMP1 (R+1) * INT2DZ (R+1)
     +                      + TEMP1 (R+2) * INT2DZ (R+2)
     +                      + TEMP1 (R+3) * INT2DZ (R+3)
     +                      + TEMP1 (R+4) * INT2DZ (R+4)
     +                      + TEMP1 (R+5) * INT2DZ (R+5)
     +                      + TEMP1 (R+6) * INT2DZ (R+6)
                R = R + 7
             END DO
         END IF

         RETURN
C
C
C                       ********************
C                       *  # of roots = 8  *
C                       ********************
C
C
    8    DO M = 1,NGQEXQ
            TEMP1 (M) = SCALE (M) * INT2DX (M)
         END DO

         IF (.NOT.DIFFY .AND. .NOT.DIFFZ) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)
     +                      + TEMP1 (R+1)
     +                      + TEMP1 (R+2)
     +                      + TEMP1 (R+3)
     +                      + TEMP1 (R+4)
     +                      + TEMP1 (R+5)
     +                      + TEMP1 (R+6)
     +                      + TEMP1 (R+7)
                R = R + 8
             END DO
         ELSE IF (.NOT.DIFFY) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DZ (R)
     +                      + TEMP1 (R+1) * INT2DZ (R+1)
     +                      + TEMP1 (R+2) * INT2DZ (R+2)
     +                      + TEMP1 (R+3) * INT2DZ (R+3)
     +                      + TEMP1 (R+4) * INT2DZ (R+4)
     +                      + TEMP1 (R+5) * INT2DZ (R+5)
     +                      + TEMP1 (R+6) * INT2DZ (R+6)
     +                      + TEMP1 (R+7) * INT2DZ (R+7)
                R = R + 8
             END DO
         ELSE IF (.NOT.DIFFZ) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DY (R)
     +                      + TEMP1 (R+1) * INT2DY (R+1)
     +                      + TEMP1 (R+2) * INT2DY (R+2)
     +                      + TEMP1 (R+3) * INT2DY (R+3)
     +                      + TEMP1 (R+4) * INT2DY (R+4)
     +                      + TEMP1 (R+5) * INT2DY (R+5)
     +                      + TEMP1 (R+6) * INT2DY (R+6)
     +                      + TEMP1 (R+7) * INT2DY (R+7)
                R = R + 8
             END DO
         ELSE
             DO N = 1,NGQEXQ
                TEMP1 (N) = TEMP1 (N) * INT2DY (N)
             END DO
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DZ (R)
     +                      + TEMP1 (R+1) * INT2DZ (R+1)
     +                      + TEMP1 (R+2) * INT2DZ (R+2)
     +                      + TEMP1 (R+3) * INT2DZ (R+3)
     +                      + TEMP1 (R+4) * INT2DZ (R+4)
     +                      + TEMP1 (R+5) * INT2DZ (R+5)
     +                      + TEMP1 (R+6) * INT2DZ (R+6)
     +                      + TEMP1 (R+7) * INT2DZ (R+7)
                R = R + 8
             END DO
         END IF

         RETURN
C
C
C                       ********************
C                       *  # of roots = 9  *
C                       ********************
C
C
    9    DO M = 1,NGQEXQ
            TEMP1 (M) = SCALE (M) * INT2DX (M)
         END DO

         IF (.NOT.DIFFY .AND. .NOT.DIFFZ) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)
     +                      + TEMP1 (R+1)
     +                      + TEMP1 (R+2)
     +                      + TEMP1 (R+3)
     +                      + TEMP1 (R+4)
     +                      + TEMP1 (R+5)
     +                      + TEMP1 (R+6)
     +                      + TEMP1 (R+7)
     +                      + TEMP1 (R+8)
                R = R + 9
             END DO
         ELSE IF (.NOT.DIFFY) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DZ (R)
     +                      + TEMP1 (R+1) * INT2DZ (R+1)
     +                      + TEMP1 (R+2) * INT2DZ (R+2)
     +                      + TEMP1 (R+3) * INT2DZ (R+3)
     +                      + TEMP1 (R+4) * INT2DZ (R+4)
     +                      + TEMP1 (R+5) * INT2DZ (R+5)
     +                      + TEMP1 (R+6) * INT2DZ (R+6)
     +                      + TEMP1 (R+7) * INT2DZ (R+7)
     +                      + TEMP1 (R+8) * INT2DZ (R+8)
                R = R + 9
             END DO
         ELSE IF (.NOT.DIFFZ) THEN
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DY (R)
     +                      + TEMP1 (R+1) * INT2DY (R+1)
     +                      + TEMP1 (R+2) * INT2DY (R+2)
     +                      + TEMP1 (R+3) * INT2DY (R+3)
     +                      + TEMP1 (R+4) * INT2DY (R+4)
     +                      + TEMP1 (R+5) * INT2DY (R+5)
     +                      + TEMP1 (R+6) * INT2DY (R+6)
     +                      + TEMP1 (R+7) * INT2DY (R+7)
     +                      + TEMP1 (R+8) * INT2DY (R+8)
                R = R + 9
             END DO
         ELSE
             DO N = 1,NGQEXQ
                TEMP1 (N) = TEMP1 (N) * INT2DY (N)
             END DO
             R = 1
             DO M = 1,NEXQ
                BATCH (M) =   TEMP1 (R)   * INT2DZ (R)
     +                      + TEMP1 (R+1) * INT2DZ (R+1)
     +                      + TEMP1 (R+2) * INT2DZ (R+2)
     +                      + TEMP1 (R+3) * INT2DZ (R+3)
     +                      + TEMP1 (R+4) * INT2DZ (R+4)
     +                      + TEMP1 (R+5) * INT2DZ (R+5)
     +                      + TEMP1 (R+6) * INT2DZ (R+6)
     +                      + TEMP1 (R+7) * INT2DZ (R+7)
     +                      + TEMP1 (R+8) * INT2DZ (R+8)
                R = R + 9
             END DO
         END IF

         RETURN
C
C
C                       ********************
C                       *  # of roots > 9  *
C                       ********************
C
C
   10    DO M = 1,NGQEXQ
            TEMP1 (M) = SCALE (M) * INT2DX (M)
         END DO

         IF (.NOT.DIFFY .AND. .NOT.DIFFZ) THEN
             R = 0
             DO M = 1,NEXQ
                SUM = ZERO
                DO N = 1,NGQP
                   SUM = SUM + TEMP1 (R+N)
                END DO
                R = R + NGQP
                BATCH (M) = SUM
             END DO
         ELSE IF (.NOT.DIFFY) THEN
             R = 0
             DO M = 1,NEXQ
                SUM = ZERO
                DO N = 1,NGQP
                   SUM = SUM + TEMP1 (R+N) * INT2DZ (R+N)
                END DO
                R = R + NGQP
                BATCH (M) = SUM
             END DO
         ELSE IF (.NOT.DIFFZ) THEN
             R = 0
             DO M = 1,NEXQ
                SUM = ZERO
                DO N = 1,NGQP
                   SUM = SUM + TEMP1 (R+N) * INT2DY (R+N)
                END DO
                R = R + NGQP
                BATCH (M) = SUM
             END DO
         ELSE
             DO N = 1,NGQEXQ
                TEMP1 (N) = TEMP1 (N) * INT2DY (N)
             END DO
             R = 0
             DO M = 1,NEXQ
                SUM = ZERO
                DO N = 1,NGQP
                   SUM = SUM + TEMP1 (R+N) * INT2DZ (R+N)
                END DO
                R = R + NGQP
                BATCH (M) = SUM
             END DO
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
