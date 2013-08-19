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
         SUBROUTINE  ERD__DERV_INT2D_TO_ABC0
     +
     +                    ( SHELLA,SHELLB,SHELLC,
     +                      NGQP,NEXQ,NGQEXQ,
     +                      NXYZA,NXYZB,NXYZC,
     +                      INT2DX,INT2DY,INT2DZ,
     +                      DIFFY,DIFFZ,
     +                      TEMP1,TEMP2,
     +                      SCALE,
     +
     +                                BATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__DERV_INT2D_TO_ABC0
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine assembles the set of batches of cartesian
C                derivative eris:
C
C                       [AB|C0] , [AB|0C] , [A0|BC] or [0A|BC]
C
C                adding up the contributions from all the respective
C                derivative 2D integrals:
C
C                          ABC0 , AB0C , A0BC or 0ABC
C
C                Simplified version of the general ABCD routine to
C                reduce loop overheads for those cases where there is
C                at least one s-shell. For comments and details see
C                the general ABCD routine.
C
C
C                  Input:
C
C                    SHELLx      =  shell types for individual csh
C                                   x = A,B,C
C                    NGQP        =  # of gaussian quadrature points
C                                   (roots)
C                    NEXQ        =  current # of exponent quadruplets
C                    NGQEXQ      =  product of # of gaussian quadrature
C                                   points times exponent quadruplets
C                    NXYZx       =  # of cartesian monomials for
C                                   x = A,B,C shells
C                    INT2Dx      =  all current 2D ABC0/AB0C/A0BC/0ABC
C                                   derivative integrals for each
C                                   cartesian component (x = X,Y,Z)
C                    DIFFx       =  is true, if differentiation was
C                                   performed along the x=Y,Z direction
C                    TEMP1(2)    =  scratch arrays holding intermediate
C                                   2D ABC0/AB0C/A0BC/0ABC derivative
C                                   integral products
C                    SCALE       =  the NGQEXQ scaling factors
C
C
C                  Output:
C
C                    BATCH       =  batch of primitive cartesian
C                                   [AB|C0]/[AB|0C]/[A0|BC]/[0A|BC]
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

         INTEGER     I,J,K,M,N,R
         INTEGER     NGQP,NEXQ,NGQEXQ
         INTEGER     NXYZA,NXYZB,NXYZC
         INTEGER     SHELLA,SHELLB,SHELLC
         INTEGER     XA,YA,ZA,XB,YB,ZB,XC,YC,ZC
         INTEGER     XAP,XBP,XCP
         INTEGER     YAMAX,YBMAX,YCMAX

         DOUBLE PRECISION  SUM
         DOUBLE PRECISION  ZERO

         DOUBLE PRECISION  SCALE (1:NGQEXQ)
         DOUBLE PRECISION  TEMP1 (1:NGQEXQ)
         DOUBLE PRECISION  TEMP2 (1:NGQEXQ)

         DOUBLE PRECISION  BATCH  (1:NEXQ,1:NXYZA,1:NXYZB,1:NXYZC)

         DOUBLE PRECISION  INT2DX (1:NGQEXQ,0:SHELLA,0:SHELLB,0:SHELLC)
         DOUBLE PRECISION  INT2DY (1:NGQEXQ,0:SHELLA,0:SHELLB,0:SHELLC)
         DOUBLE PRECISION  INT2DZ (1:NGQEXQ,0:SHELLA,0:SHELLB,0:SHELLC)

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
    1    XCP = 0
         DO 100 XC = SHELLC,0,-1
           YCMAX = SHELLC - XC
           XBP = 0
           DO 102  XB = SHELLB,0,-1
             YBMAX = SHELLB - XB
             XAP = 0
             DO 104 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC)
               END DO

               K = XCP
               DO 110 YC = YCMAX,0,-1
                 K = K + 1
                 ZC = YCMAX - YC
                 J = XBP
                 DO 112 YB = YBMAX,0,-1
                   J = J + 1
                   ZB = YBMAX - YB
                   I = XAP
                   DO 114 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB+YC.EQ.0) THEN
                         DO M = 1,NEXQ
                            TEMP2 (M) = TEMP1 (M)
                         END DO
                     ELSE
                         DO M = 1,NEXQ
                            TEMP2 (M) = TEMP1 (M) * INT2DY (M,YA,YB,YC)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB+ZC.EQ.0) THEN
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) = TEMP2 (M)
                         END DO
                     ELSE
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (M)
     +                                        * INT2DZ (M,ZA,ZB,ZC)
                         END DO
                     END IF

  114              CONTINUE
  112            CONTINUE
  110          CONTINUE
               XAP = I
  104        CONTINUE
             XBP = J
  102      CONTINUE
           XCP = K
  100    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 2  *
C                       ********************
C
C
    2    XCP = 0
         DO 200 XC = SHELLC,0,-1
           YCMAX = SHELLC - XC
           XBP = 0
           DO 202  XB = SHELLB,0,-1
             YBMAX = SHELLB - XB
             XAP = 0
             DO 204 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC)
               END DO

               K = XCP
               DO 210 YC = YCMAX,0,-1
                 K = K + 1
                 ZC = YCMAX - YC
                 J = XBP
                 DO 212 YB = YBMAX,0,-1
                   J = J + 1
                   ZB = YBMAX - YB
                   I = XAP
                   DO 214 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB,YC)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB+ZC.EQ.0) THEN
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (R)
     +                                        + TEMP2 (R+1)
                            R = R + 2
                         END DO
                     ELSE
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (R)
     +                                        * INT2DZ (R,ZA,ZB,ZC)
     +                                        + TEMP2 (R+1)
     +                                        * INT2DZ (R+1,ZA,ZB,ZC)
                            R = R + 2
                         END DO
                     END IF

  214              CONTINUE
  212            CONTINUE
  210          CONTINUE
               XAP = I
  204        CONTINUE
             XBP = J
  202      CONTINUE
           XCP = K
  200    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 3  *
C                       ********************
C
C
    3    XCP = 0
         DO 300 XC = SHELLC,0,-1
           YCMAX = SHELLC - XC
           XBP = 0
           DO 302  XB = SHELLB,0,-1
             YBMAX = SHELLB - XB
             XAP = 0
             DO 304 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC)
               END DO

               K = XCP
               DO 310 YC = YCMAX,0,-1
                 K = K + 1
                 ZC = YCMAX - YC
                 J = XBP
                 DO 312 YB = YBMAX,0,-1
                   J = J + 1
                   ZB = YBMAX - YB
                   I = XAP
                   DO 314 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB,YC)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB+ZC.EQ.0) THEN
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (R)
     +                                        + TEMP2 (R+1)
     +                                        + TEMP2 (R+2)
                            R = R + 3
                         END DO
                     ELSE
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (R)
     +                                        * INT2DZ (R,ZA,ZB,ZC)
     +                                        + TEMP2 (R+1)
     +                                        * INT2DZ (R+1,ZA,ZB,ZC)
     +                                        + TEMP2 (R+2)
     +                                        * INT2DZ (R+2,ZA,ZB,ZC)
                            R = R + 3
                         END DO
                     END IF

  314              CONTINUE
  312            CONTINUE
  310          CONTINUE
               XAP = I
  304        CONTINUE
             XBP = J
  302      CONTINUE
           XCP = K
  300    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 4  *
C                       ********************
C
C
    4    XCP = 0
         DO 400 XC = SHELLC,0,-1
           YCMAX = SHELLC - XC
           XBP = 0
           DO 402  XB = SHELLB,0,-1
             YBMAX = SHELLB - XB
             XAP = 0
             DO 404 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC)
               END DO

               K = XCP
               DO 410 YC = YCMAX,0,-1
                 K = K + 1
                 ZC = YCMAX - YC
                 J = XBP
                 DO 412 YB = YBMAX,0,-1
                   J = J + 1
                   ZB = YBMAX - YB
                   I = XAP
                   DO 414 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB,YC)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB+ZC.EQ.0) THEN
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (R)
     +                                        + TEMP2 (R+1)
     +                                        + TEMP2 (R+2)
     +                                        + TEMP2 (R+3)
                            R = R + 4
                         END DO
                     ELSE
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (R)
     +                                        * INT2DZ (R,ZA,ZB,ZC)
     +                                        + TEMP2 (R+1)
     +                                        * INT2DZ (R+1,ZA,ZB,ZC)
     +                                        + TEMP2 (R+2)
     +                                        * INT2DZ (R+2,ZA,ZB,ZC)
     +                                        + TEMP2 (R+3)
     +                                        * INT2DZ (R+3,ZA,ZB,ZC)
                            R = R + 4
                         END DO
                     END IF

  414              CONTINUE
  412            CONTINUE
  410          CONTINUE
               XAP = I
  404        CONTINUE
             XBP = J
  402      CONTINUE
           XCP = K
  400    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 5  *
C                       ********************
C
C
    5    XCP = 0
         DO 500 XC = SHELLC,0,-1
           YCMAX = SHELLC - XC
           XBP = 0
           DO 502  XB = SHELLB,0,-1
             YBMAX = SHELLB - XB
             XAP = 0
             DO 504 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC)
               END DO

               K = XCP
               DO 510 YC = YCMAX,0,-1
                 K = K + 1
                 ZC = YCMAX - YC
                 J = XBP
                 DO 512 YB = YBMAX,0,-1
                   J = J + 1
                   ZB = YBMAX - YB
                   I = XAP
                   DO 514 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB,YC)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB+ZC.EQ.0) THEN
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (R)
     +                                        + TEMP2 (R+1)
     +                                        + TEMP2 (R+2)
     +                                        + TEMP2 (R+3)
     +                                        + TEMP2 (R+4)
                            R = R + 5
                         END DO
                     ELSE
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (R)
     +                                        * INT2DZ (R,ZA,ZB,ZC)
     +                                        + TEMP2 (R+1)
     +                                        * INT2DZ (R+1,ZA,ZB,ZC)
     +                                        + TEMP2 (R+2)
     +                                        * INT2DZ (R+2,ZA,ZB,ZC)
     +                                        + TEMP2 (R+3)
     +                                        * INT2DZ (R+3,ZA,ZB,ZC)
     +                                        + TEMP2 (R+4)
     +                                        * INT2DZ (R+4,ZA,ZB,ZC)
                            R = R + 5
                         END DO
                     END IF

  514              CONTINUE
  512            CONTINUE
  510          CONTINUE
               XAP = I
  504        CONTINUE
             XBP = J
  502      CONTINUE
           XCP = K
  500    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 6  *
C                       ********************
C
C
    6    XCP = 0
         DO 600 XC = SHELLC,0,-1
           YCMAX = SHELLC - XC
           XBP = 0
           DO 602  XB = SHELLB,0,-1
             YBMAX = SHELLB - XB
             XAP = 0
             DO 604 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC)
               END DO

               K = XCP
               DO 610 YC = YCMAX,0,-1
                 K = K + 1
                 ZC = YCMAX - YC
                 J = XBP
                 DO 612 YB = YBMAX,0,-1
                   J = J + 1
                   ZB = YBMAX - YB
                   I = XAP
                   DO 614 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB,YC)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB+ZC.EQ.0) THEN
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (R)
     +                                        + TEMP2 (R+1)
     +                                        + TEMP2 (R+2)
     +                                        + TEMP2 (R+3)
     +                                        + TEMP2 (R+4)
     +                                        + TEMP2 (R+5)
                            R = R + 6
                         END DO
                     ELSE
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (R)
     +                                        * INT2DZ (R,ZA,ZB,ZC)
     +                                        + TEMP2 (R+1)
     +                                        * INT2DZ (R+1,ZA,ZB,ZC)
     +                                        + TEMP2 (R+2)
     +                                        * INT2DZ (R+2,ZA,ZB,ZC)
     +                                        + TEMP2 (R+3)
     +                                        * INT2DZ (R+3,ZA,ZB,ZC)
     +                                        + TEMP2 (R+4)
     +                                        * INT2DZ (R+4,ZA,ZB,ZC)
     +                                        + TEMP2 (R+5)
     +                                        * INT2DZ (R+5,ZA,ZB,ZC)
                            R = R + 6
                         END DO
                     END IF

  614              CONTINUE
  612            CONTINUE
  610          CONTINUE
               XAP = I
  604        CONTINUE
             XBP = J
  602      CONTINUE
           XCP = K
  600    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 7  *
C                       ********************
C
C
    7    XCP = 0
         DO 700 XC = SHELLC,0,-1
           YCMAX = SHELLC - XC
           XBP = 0
           DO 702  XB = SHELLB,0,-1
             YBMAX = SHELLB - XB
             XAP = 0
             DO 704 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC)
               END DO

               K = XCP
               DO 710 YC = YCMAX,0,-1
                 K = K + 1
                 ZC = YCMAX - YC
                 J = XBP
                 DO 712 YB = YBMAX,0,-1
                   J = J + 1
                   ZB = YBMAX - YB
                   I = XAP
                   DO 714 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB,YC)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB+ZC.EQ.0) THEN
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (R)
     +                                        + TEMP2 (R+1)
     +                                        + TEMP2 (R+2)
     +                                        + TEMP2 (R+3)
     +                                        + TEMP2 (R+4)
     +                                        + TEMP2 (R+5)
     +                                        + TEMP2 (R+6)
                            R = R + 7
                         END DO
                     ELSE
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (R)
     +                                        * INT2DZ (R,ZA,ZB,ZC)
     +                                        + TEMP2 (R+1)
     +                                        * INT2DZ (R+1,ZA,ZB,ZC)
     +                                        + TEMP2 (R+2)
     +                                        * INT2DZ (R+2,ZA,ZB,ZC)
     +                                        + TEMP2 (R+3)
     +                                        * INT2DZ (R+3,ZA,ZB,ZC)
     +                                        + TEMP2 (R+4)
     +                                        * INT2DZ (R+4,ZA,ZB,ZC)
     +                                        + TEMP2 (R+5)
     +                                        * INT2DZ (R+5,ZA,ZB,ZC)
     +                                        + TEMP2 (R+6)
     +                                        * INT2DZ (R+6,ZA,ZB,ZC)
                            R = R + 7
                         END DO
                     END IF

  714              CONTINUE
  712            CONTINUE
  710          CONTINUE
               XAP = I
  704        CONTINUE
             XBP = J
  702      CONTINUE
           XCP = K
  700    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 8  *
C                       ********************
C
C
    8    XCP = 0
         DO 800 XC = SHELLC,0,-1
           YCMAX = SHELLC - XC
           XBP = 0
           DO 802  XB = SHELLB,0,-1
             YBMAX = SHELLB - XB
             XAP = 0
             DO 804 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC)
               END DO

               K = XCP
               DO 810 YC = YCMAX,0,-1
                 K = K + 1
                 ZC = YCMAX - YC
                 J = XBP
                 DO 812 YB = YBMAX,0,-1
                   J = J + 1
                   ZB = YBMAX - YB
                   I = XAP
                   DO 814 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB,YC)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB+ZC.EQ.0) THEN
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (R)
     +                                        + TEMP2 (R+1)
     +                                        + TEMP2 (R+2)
     +                                        + TEMP2 (R+3)
     +                                        + TEMP2 (R+4)
     +                                        + TEMP2 (R+5)
     +                                        + TEMP2 (R+6)
     +                                        + TEMP2 (R+7)
                            R = R + 8
                         END DO
                     ELSE
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (R)
     +                                        * INT2DZ (R,ZA,ZB,ZC)
     +                                        + TEMP2 (R+1)
     +                                        * INT2DZ (R+1,ZA,ZB,ZC)
     +                                        + TEMP2 (R+2)
     +                                        * INT2DZ (R+2,ZA,ZB,ZC)
     +                                        + TEMP2 (R+3)
     +                                        * INT2DZ (R+3,ZA,ZB,ZC)
     +                                        + TEMP2 (R+4)
     +                                        * INT2DZ (R+4,ZA,ZB,ZC)
     +                                        + TEMP2 (R+5)
     +                                        * INT2DZ (R+5,ZA,ZB,ZC)
     +                                        + TEMP2 (R+6)
     +                                        * INT2DZ (R+6,ZA,ZB,ZC)
     +                                        + TEMP2 (R+7)
     +                                        * INT2DZ (R+7,ZA,ZB,ZC)
                            R = R + 8
                         END DO
                     END IF

  814              CONTINUE
  812            CONTINUE
  810          CONTINUE
               XAP = I
  804        CONTINUE
             XBP = J
  802      CONTINUE
           XCP = K
  800    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots = 9  *
C                       ********************
C
C
    9    XCP = 0
         DO 900 XC = SHELLC,0,-1
           YCMAX = SHELLC - XC
           XBP = 0
           DO 902  XB = SHELLB,0,-1
             YBMAX = SHELLB - XB
             XAP = 0
             DO 904 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC)
               END DO

               K = XCP
               DO 910 YC = YCMAX,0,-1
                 K = K + 1
                 ZC = YCMAX - YC
                 J = XBP
                 DO 912 YB = YBMAX,0,-1
                   J = J + 1
                   ZB = YBMAX - YB
                   I = XAP
                   DO 914 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB,YC)
                         END DO
                     END IF

                     IF (.NOT.DIFFZ .AND. ZA+ZB+ZC.EQ.0) THEN
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (R)
     +                                        + TEMP2 (R+1)
     +                                        + TEMP2 (R+2)
     +                                        + TEMP2 (R+3)
     +                                        + TEMP2 (R+4)
     +                                        + TEMP2 (R+5)
     +                                        + TEMP2 (R+6)
     +                                        + TEMP2 (R+7)
     +                                        + TEMP2 (R+8)
                            R = R + 9
                         END DO
                     ELSE
                         R = 1
                         DO M = 1,NEXQ
                            BATCH (M,I,J,K) =   TEMP2 (R)
     +                                        * INT2DZ (R,ZA,ZB,ZC)
     +                                        + TEMP2 (R+1)
     +                                        * INT2DZ (R+1,ZA,ZB,ZC)
     +                                        + TEMP2 (R+2)
     +                                        * INT2DZ (R+2,ZA,ZB,ZC)
     +                                        + TEMP2 (R+3)
     +                                        * INT2DZ (R+3,ZA,ZB,ZC)
     +                                        + TEMP2 (R+4)
     +                                        * INT2DZ (R+4,ZA,ZB,ZC)
     +                                        + TEMP2 (R+5)
     +                                        * INT2DZ (R+5,ZA,ZB,ZC)
     +                                        + TEMP2 (R+6)
     +                                        * INT2DZ (R+6,ZA,ZB,ZC)
     +                                        + TEMP2 (R+7)
     +                                        * INT2DZ (R+7,ZA,ZB,ZC)
     +                                        + TEMP2 (R+8)
     +                                        * INT2DZ (R+8,ZA,ZB,ZC)
                            R = R + 9
                         END DO
                     END IF

  914              CONTINUE
  912            CONTINUE
  910          CONTINUE
               XAP = I
  904        CONTINUE
             XBP = J
  902      CONTINUE
           XCP = K
  900    CONTINUE

         RETURN
C
C
C                       ********************
C                       *  # of roots > 9  *
C                       ********************
C
C             ...outer loops over x,x,x-triples. No skipping of
C                x,x,x-contribution of 0,0,0-type can be done here,
C                since the 2DX integrals carry the Rys weight!
C
C
   10    XCP = 0
         DO 1000 XC = SHELLC,0,-1
           YCMAX = SHELLC - XC
           XBP = 0
           DO 1002  XB = SHELLB,0,-1
             YBMAX = SHELLB - XB
             XAP = 0
             DO 1004 XA = SHELLA,0,-1
               YAMAX = SHELLA - XA

               DO M = 1,NGQEXQ
                  TEMP1 (M) = SCALE (M) * INT2DX (M,XA,XB,XC)
               END DO
C
C
C             ...inner loops over y,y,y-triples. Skip the
C                multiplication of y,y,y-contributions, if no
C                y-coordinate derivative was formed and we have a
C                0,0,0-triple, as then the 2DY integrals are
C                equal to 1.
C
C
               K = XCP
               DO 1010 YC = YCMAX,0,-1
                 K = K + 1
                 ZC = YCMAX - YC
                 J = XBP
                 DO 1012 YB = YBMAX,0,-1
                   J = J + 1
                   ZB = YBMAX - YB
                   I = XAP
                   DO 1014 YA = YAMAX,0,-1
                     I = I + 1
                     ZA = YAMAX - YA

                     IF (.NOT.DIFFY .AND. YA+YB+YC.EQ.0) THEN
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N)
                         END DO
                     ELSE
                         DO N = 1,NGQEXQ
                            TEMP2 (N) = TEMP1 (N) * INT2DY (N,YA,YB,YC)
                         END DO
                     END IF
C
C
C             ...skip multiplication of z,z,z-contributions, if we
C                have a 0,0,0-triple and no derivations were performed
C                on the z-coordinate, as then the 2DZ integrals
C                are equal to 1. All info concerning all three x,x,x-,
C                y,y,y- and z,z,z-triples have been collected for all
C                exponent quadruplets at once. Sum up the 2D X,Y,Z
C                integral products to the appropriate places of the
C                batch.
C
C
                     IF (.NOT.DIFFZ .AND. ZA+ZB+ZC.EQ.0) THEN
                         R = 0
                         DO M = 1,NEXQ
                            SUM = ZERO
                            DO N = 1,NGQP
                               SUM = SUM + TEMP2 (R+N)
                            END DO
                            R = R + NGQP
                            BATCH (M,I,J,K) = SUM
                         END DO
                     ELSE
                         R = 0
                         DO M = 1,NEXQ
                            SUM = ZERO
                            DO N = 1,NGQP
                               SUM = SUM + TEMP2 (R+N)
     +                                   * INT2DZ (R+N,ZA,ZB,ZC)
                            END DO
                            R = R + NGQP
                            BATCH (M,I,J,K) = SUM
                         END DO
                     END IF
C
C
C             ...next z,z,z-triple and y,y,y-triple.
C
C
 1014              CONTINUE
 1012            CONTINUE
 1010          CONTINUE
C
C
C             ...next x,x,x-triple.
C
C
               XAP = I
 1004        CONTINUE
             XBP = J
 1002      CONTINUE
           XCP = K
 1000    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
