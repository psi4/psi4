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
         DOUBLE PRECISION FUNCTION  OED__DSQMIN_LINE_SEGMENTS
     +
     +                     ( XP0,YP0,ZP0,
     +                       XP1,YP1,ZP1,
     +                       XQ0,YQ0,ZQ0,
     +                       XQ1,YQ1,ZQ1 )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__DSQMIN_LINE_SEGMENTS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This function returns the square of the minimum 3D
C                distance that exists between two line segments P and Q
C                in space.
C
C                The line segment P is located between points P0 and
C                P1 and likewise for line segment Q between points
C                Q0 and Q1. The routine also handles the situations
C                when the line segments P and/or Q have zero lenghts,
C                i.e. they are points in space.
C
C
C
C                        \
C                         \
C                          P0
C                           \ <-- segment P, dir vector U = P1 - P0
C                            \
C                             P1
C                              \
C                               \          /
C                                \        / <-- line L2
C                                 P(sc)  /
C                                 |\    /
C                                 | \  /
C             W = P0 - Q0         |  \/
C      C(sc,tc) = P(sc) - Q(tc)   |  /\
C                                 | /  \
C                                 |/    \
C                                 Q(tc)  \
C                                /        \
C                               /          \ <-- line L1
C                              /            \
C                             Q1             \
C                            /
C                           / <-- segment Q, dir vector V = Q1 - Q0
C                          Q0
C                         /
C
C
C
C                The two lines L1 and L2 are the extension lines on
C                which the two segments P and Q sit and are in given
C                in parametric form using the direction vectors U
C                and V as:
C
C                                 L1 = P0 + sU
C                                 L2 = Q0 + tV
C
C                with s and t parameters. C(sc,tc) is the vector
C                representing the closest approach between these two
C                lines occuring at s = sc and t = tc and is uniquely
C                perpendicular to the line direction vectors U and V:
C
C                                 U dot C(sc,tc) = 0
C                                 V dot C(sc,tc) = 0
C
C                Using C(sc,tc) = P(sc) - Q(tc) = P0 + sc*U - Q0 - tc*V
C                = W + sc*U - tc*V for these two conditions we obtain
C                two simultaneous linear equations in sc and tc:
C
C
C                     (U dot U)*sc - (U dot V)*tc = - (U dot W)
C                     (V dot U)*sc - (V dot V)*tc = - (V dot W)
C
C
C                from which we obtain:
C
C                          sc = (b*e - c*d) / (a*c - b*b)
C                          tc = (a*e - b*d) / (a*c - b*b)
C
C                where:
C
C                                 a = U dot U
C                                 b = U dot V
C                                 c = V dot V
C                                 d = U dot W
C                                 e = V dot W
C
C                If the denominator (a*c - b*b) is equal to zero, then
C                the lines L1 and L2 are parallel and the distance
C                between them is constant. We can solve for this
C                parallel distance of separation by fixing the value of
C                one parameter and using either equation to solve for
C                the other. If we select sc = 0, then we get tc = d/b
C                = e/c. The situation when one or both of the line
C                segments shrink to points is identified as a and/or c
C                being equal to zero.
C
C                The closest distance between segments may not be the
C                same as the closest distance between their extended
C                lines. The closest points on the extended infinite line
C                may be outside the range of the segments. Using the
C                parametrized line form for L1 = P0 + sU, we see that
C                the segment P is represented by all points on L1 for
C                which 0 =< s =< 1. Likewise for the segment Q with
C                0 =< t =< 1.
C
C                Minimizing the distance between the two line segments
C                P and Q means to minimize the scalar product:
C
C                  C dot C = (W0 + s*U - t*V) dot (W0 + s*U - t*V)
C
C                as a function of s and t over the restricted values
C                0 =< s =< 1 and 0 =< t =< 1 defining the segments.
C
C
C                  Input:
C
C                     XP0,YP0,ZP0 =  x,y,z coordinates of first point
C                                    defining line segment P
C                     XP1,YP1,ZP1 =  x,y,z coordinates of second point
C                                    defining line segment P
C                     XQ0,YQ0,ZQ0 =  x,y,z coordinates of first point
C                                    defining line segment Q
C                     XQ1,YQ1,ZQ1 =  x,y,z coordinates of second point
C                                    defining line segment Q
C
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         DOUBLE PRECISION    DENOM
         DOUBLE PRECISION    LINE
         DOUBLE PRECISION    PARLLEL
         DOUBLE PRECISION    SD,SN,SC,TD,TN,TC
         DOUBLE PRECISION    UDOTU,VDOTV,UDOTV,UDOTW,VDOTW
         DOUBLE PRECISION    XC,YC,ZC
         DOUBLE PRECISION    XP0,YP0,ZP0,XP1,YP1,ZP1
         DOUBLE PRECISION    XQ0,YQ0,ZQ0,XQ1,YQ1,ZQ1
         DOUBLE PRECISION    XU,YU,ZU
         DOUBLE PRECISION    XV,YV,ZV
         DOUBLE PRECISION    XW,YW,ZW
         DOUBLE PRECISION    ZERO,ONE

         PARAMETER  (PARLLEL =  1.D-12)
         PARAMETER  (LINE    =  1.D-12)
         PARAMETER  (ZERO    =  0.D0)
         PARAMETER  (ONE     =  1.D0)
C
C
C------------------------------------------------------------------------
C
C
C             ...set direction vectors U and V and reference vector W.
C
C
         XU = XP1 - XP0
         YU = YP1 - YP0
         ZU = ZP1 - ZP0
         XV = XQ1 - XQ0
         YV = YQ1 - YQ0
         ZV = ZQ1 - ZQ0
         XW = XP0 - XQ0
         YW = YP0 - YQ0
         ZW = ZP0 - ZQ0

         UDOTU = XU * XU + YU * YU + ZU * ZU
         VDOTV = XV * XV + YV * YV + ZV * ZV

         IF (UDOTU.LT.LINE .AND. VDOTV.LT.LINE) THEN
C
C
C             ...both line segments are points.
C
C
             OED__DSQMIN_LINE_SEGMENTS = XW * XW + YW * YW + ZW * ZW
             RETURN

         ELSE IF (UDOTU.LT.LINE) THEN
C
C
C             ...only P line segment is a point.
C
C
             VDOTW = XV * XW + YV * YW + ZV * ZW

             IF (VDOTW.LT.ZERO) THEN
                 OED__DSQMIN_LINE_SEGMENTS = XW * XW + YW * YW + ZW * ZW
                 RETURN
             ELSE IF (VDOTW.GT.VDOTV) THEN
                 XC = XW - XV
                 YC = YW - YV
                 ZC = ZW - ZV
             ELSE
                 TC = VDOTW / VDOTV
                 XC = XW - TC * XV
                 YC = YW - TC * YV
                 ZC = ZW - TC * ZV
             END IF

         ELSE IF (VDOTV.LT.LINE) THEN
C
C
C             ...only Q line segment is a point.
C
C
             UDOTW = - (XU * XW + YU * YW + ZU * ZW)

             IF (UDOTW.LT.ZERO) THEN
                 OED__DSQMIN_LINE_SEGMENTS = XW * XW + YW * YW + ZW * ZW
                 RETURN
             ELSE IF (UDOTW.GT.UDOTU) THEN
                 XC = XW + XU
                 YC = YW + YU
                 ZC = ZW + ZU
             ELSE
                 SC = UDOTW / UDOTU
                 XC = XW + SC * XU
                 YC = YW + SC * YU
                 ZC = ZW + SC * ZU
             END IF

         ELSE
C
C
C           ...both P and Q are line segments.
C
C
             UDOTV = XU * XV + YU * YV + ZU * ZV
             UDOTW = XU * XW + YU * YW + ZU * ZW
             VDOTW = XV * XW + YV * YW + ZV * ZW

             DENOM = UDOTU * VDOTV - UDOTV * UDOTV

             IF (DENOM.LT.PARLLEL) THEN
                 SN = ZERO
                 SD = ONE
                 TN = VDOTW
                 TD = VDOTV
             ELSE
                 SN = UDOTV * VDOTW - VDOTV * UDOTW
                 SD = DENOM
                 IF (SN.LT.ZERO) THEN
                     SN = ZERO
                     TN = VDOTW
                     TD = DENOM
                 ELSE IF (SN.GT.DENOM) THEN
                     SN = DENOM
                     SD = DENOM
                     TN = VDOTW + UDOTV
                     TD = VDOTV
                 ELSE
                     TN = UDOTU * VDOTW - UDOTV * UDOTW
                     TD = DENOM
                 END IF
             END IF

             IF (TN.LT.ZERO) THEN
                 TN = ZERO
                 SN = - UDOTW
                 IF (SN.LT.ZERO) THEN
                     SN = ZERO
                 ELSE IF (SN.GT.UDOTU) THEN
                     SN = SD
                 ELSE
                     SN = - UDOTW
                     SD = UDOTU
                 END IF
             ELSE IF (TN.GT.TD) THEN
                 TN = TD
                 SN = UDOTV - UDOTW
                 IF (SN.LT.ZERO) THEN
                     SN = ZERO
                 ELSE IF (SN.GT.UDOTU) THEN
                     SN = SD
                 ELSE
                     SD = UDOTU
                 END IF
             END IF

             SC = SN / SD
             TC = TN / TD

             XC = XW + SC * XU - TC * XV
             YC = YW + SC * YU - TC * YV
             ZC = ZW + SC * ZU - TC * ZV

         END IF

         OED__DSQMIN_LINE_SEGMENTS = XC * XC + YC * YC + ZC * ZC
C
C
C             ...ready!
C
C
         RETURN
         END
