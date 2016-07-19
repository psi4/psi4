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
         SUBROUTINE  ERD__SET_DERV_SEQUENCE
     +
     +                    ( NDERX,NDERY,NDERZ,
     +                      DERAX,DERAY,DERAZ,
     +                      DERBX,DERBY,DERBZ,
     +                      DERCX,DERCY,DERCZ,
     +                      DERDX,DERDY,DERDZ,
     +                      CENEQS,
     +
     +                               CENSQX,CENSQY,CENSQZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__SET_DERV_SEQUENCE
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine determines the sequence of centers for
C                the x,y,z coordinate derivations. The info is stored
C                in arrays CENSQp for each coordinate p=X,Y,Z. Take
C                the x coordinate for example and let CENSQX have the
C                following elements on exit from this routine:
C
C                            CENSQX (I) = 2,4,4,1
C
C                Then the sequence of single derivations is (observing
C                the index convention I=1,2,3,4 -> A,B,C,D):
C
C                     dx/dAx  dx/dDx  dx/dDx  dx/dBx  (AB|CD)
C
C                Note, that the indexing of the CENSQp arrays start
C                with 0. The reason is that the NDERp values might
C                be equal to zero.
C
C
C                  Input:
C
C                   NDERp         =  the total order of differentiation
C                                    with respect to the p = x,y,z
C                                    coordinates
C                   DERxp         =  the order of differentiation on
C                                    centers x = A,B,C,D with respect
C                                    to the p = x,y,z coordinates
C                   CENEQS (I,J)  =  center equality indicator of size
C                                    4 x 5. The first 4 columns are
C                                    defined as follows: if the
C                                    centers indexed by I and J are
C                                    equal => value = 1, if not =>
C                                    value = 0. Indexes correspond
C                                    to the A,B,C,D ordering, i.e.
C                                    1st index -> A, 2nd index -> B,
C                                    etc...
C
C                  Output:
C
C                   CENSQp        =  center sequence array for the
C                                    p = x,y,z coordinates
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

         INTEGER     DER1,DER2,DER3,DER4
         INTEGER     DERAX,DERBX,DERCX,DERDX
         INTEGER     DERAY,DERBY,DERCY,DERDY
         INTEGER     DERAZ,DERBZ,DERCZ,DERDZ
         INTEGER     I,M,N
         INTEGER     NDERX,NDERY,NDERZ

         INTEGER     CENSQX (0:NDERX)
         INTEGER     CENSQY (0:NDERY)
         INTEGER     CENSQZ (0:NDERZ)

         INTEGER     CENEQS (1:4,1:4)
C
C
C------------------------------------------------------------------------
C
C
C             ...center sequence for x coordinate.
C
C
         IF (NDERX.GT.0) THEN

             N = 0

             DER1 = DERAX
             DER2 = DERBX
             DER3 = DERCX
             DER4 = DERDX

             IF (DER1.GT.0) THEN
                 M = DER1
                 DO 110 I = 1,M
                    CENSQX (I) = 1
                    DER1 = DER1 - CENEQS (1,1)
                    DER2 = DER2 - CENEQS (2,1)
                    DER3 = DER3 - CENEQS (3,1)
                    DER4 = DER4 - CENEQS (4,1)
  110            CONTINUE
                 N = M
             END IF

             IF (DER2.GT.0) THEN
                 M = DER2
                 DO 120 I = 1,M
                    CENSQX (N+I) = 2
                    DER1 = DER1 - CENEQS (1,2)
                    DER2 = DER2 - CENEQS (2,2)
                    DER3 = DER3 - CENEQS (3,2)
                    DER4 = DER4 - CENEQS (4,2)
  120            CONTINUE
                 N = N + M
             END IF

             IF (DER3.GT.0) THEN
                 M = DER3
                 DO 130 I = 1,M
                    CENSQX (N+I) = 3
                    DER1 = DER1 - CENEQS (1,3)
                    DER2 = DER2 - CENEQS (2,3)
                    DER3 = DER3 - CENEQS (3,3)
                    DER4 = DER4 - CENEQS (4,3)
  130            CONTINUE
                 N = N + M
             END IF

             IF (DER4.GT.0) THEN
                 M = DER4
                 DO 140 I = 1,M
                    CENSQX (N+I) = 4
                    DER1 = DER1 - CENEQS (1,4)
                    DER2 = DER2 - CENEQS (2,4)
                    DER3 = DER3 - CENEQS (3,4)
                    DER4 = DER4 - CENEQS (4,4)
  140            CONTINUE
                 N = N + M
             END IF

             IF (N.NE.NDERX) THEN
                 WRITE (*,*) ' Inconsistent x derivatives! '
                 WRITE (*,*) ' N,NDERX = ',N,NDERX
                 WRITE (*,*) ' erd__set_derv_sequence '
                 STOP
             END IF

         END IF
C
C
C             ...center sequence for y coordinate.
C
C
         IF (NDERY.GT.0) THEN

             N = 0

             DER1 = DERAY
             DER2 = DERBY
             DER3 = DERCY
             DER4 = DERDY

             IF (DER1.GT.0) THEN
                 M = DER1
                 DO 210 I = 1,M
                    CENSQY (I) = 1
                    DER1 = DER1 - CENEQS (1,1)
                    DER2 = DER2 - CENEQS (2,1)
                    DER3 = DER3 - CENEQS (3,1)
                    DER4 = DER4 - CENEQS (4,1)
  210            CONTINUE
                 N = M
             END IF

             IF (DER2.GT.0) THEN
                 M = DER2
                 DO 220 I = 1,M
                    CENSQY (N+I) = 2
                    DER1 = DER1 - CENEQS (1,2)
                    DER2 = DER2 - CENEQS (2,2)
                    DER3 = DER3 - CENEQS (3,2)
                    DER4 = DER4 - CENEQS (4,2)
  220            CONTINUE
                 N = N + M
             END IF

             IF (DER3.GT.0) THEN
                 M = DER3
                 DO 230 I = 1,M
                    CENSQY (N+I) = 3
                    DER1 = DER1 - CENEQS (1,3)
                    DER2 = DER2 - CENEQS (2,3)
                    DER3 = DER3 - CENEQS (3,3)
                    DER4 = DER4 - CENEQS (4,3)
  230            CONTINUE
                 N = N + M
             END IF

             IF (DER4.GT.0) THEN
                 M = DER4
                 DO 240 I = 1,M
                    CENSQY (N+I) = 4
                    DER1 = DER1 - CENEQS (1,4)
                    DER2 = DER2 - CENEQS (2,4)
                    DER3 = DER3 - CENEQS (3,4)
                    DER4 = DER4 - CENEQS (4,4)
  240            CONTINUE
                 N = N + M
             END IF

             IF (N.NE.NDERY) THEN
                 WRITE (*,*) ' Inconsistent y derivatives! '
                 WRITE (*,*) ' N,NDERY = ',N,NDERY
                 WRITE (*,*) ' erd__set_derv_sequence '
                 STOP
             END IF

         END IF
C
C
C             ...center sequence for z coordinate.
C
C
         IF (NDERZ.GT.0) THEN

             N = 0

             DER1 = DERAZ
             DER2 = DERBZ
             DER3 = DERCZ
             DER4 = DERDZ

             IF (DER1.GT.0) THEN
                 M = DER1
                 DO 310 I = 1,M
                    CENSQZ (I) = 1
                    DER1 = DER1 - CENEQS (1,1)
                    DER2 = DER2 - CENEQS (2,1)
                    DER3 = DER3 - CENEQS (3,1)
                    DER4 = DER4 - CENEQS (4,1)
  310            CONTINUE
                 N = M
             END IF

             IF (DER2.GT.0) THEN
                 M = DER2
                 DO 320 I = 1,M
                    CENSQZ (N+I) = 2
                    DER1 = DER1 - CENEQS (1,2)
                    DER2 = DER2 - CENEQS (2,2)
                    DER3 = DER3 - CENEQS (3,2)
                    DER4 = DER4 - CENEQS (4,2)
  320            CONTINUE
                 N = N + M
             END IF

             IF (DER3.GT.0) THEN
                 M = DER3
                 DO 330 I = 1,M
                    CENSQZ (N+I) = 3
                    DER1 = DER1 - CENEQS (1,3)
                    DER2 = DER2 - CENEQS (2,3)
                    DER3 = DER3 - CENEQS (3,3)
                    DER4 = DER4 - CENEQS (4,3)
  330            CONTINUE
                 N = N + M
             END IF

             IF (DER4.GT.0) THEN
                 M = DER4
                 DO 340 I = 1,M
                    CENSQZ (N+I) = 4
                    DER1 = DER1 - CENEQS (1,4)
                    DER2 = DER2 - CENEQS (2,4)
                    DER3 = DER3 - CENEQS (3,4)
                    DER4 = DER4 - CENEQS (4,4)
  340            CONTINUE
                 N = N + M
             END IF

             IF (N.NE.NDERZ) THEN
                 WRITE (*,*) ' Inconsistent z derivatives! '
                 WRITE (*,*) ' N,NDERZ = ',N,NDERZ
                 WRITE (*,*) ' erd__set_derv_sequence '
                 STOP
             END IF

         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
