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
         SUBROUTINE  OED__NAI_1D_COEFFICIENTS
     +
     +                    ( NGQP,NCEN,
     +                      MIJ,MGIJCEN,
     +                      ATOMIC,
     +                      NUCLEI,
     +                      XN,YN,ZN,
     +                      NUCCEN,
     +                      PX,PY,PZ,
     +                      PAX,PAY,PAZ,
     +                      PINVHF,
     +                      RTS,
     +
     +                               R1X,R1Y,R1Z,
     +                               R2 )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__NAI_1D_COEFFICIENTS
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This operation evaluates the VRR coefficients for
C                the 1D nuclear attraction integrals for the present
C                set of NGQP roots corresponding to all i,j exponent
C                pairs.
C
C
C                  Input:
C
C                    NGQP         =  # of gaussian quadrature points
C                                    (roots)
C                    NCEN         =  # of nuclear attraction centers
C                    MIJ          =  current # of ij primitive index
C                                    pairs corresponding to the
C                                    contracted shell pairs A,B
C                    MGIJCEN      =  # of roots times # of ij primitive
C                                    index pairs times # of nuclear
C                                    attraction centers
C                    ATOMIC       =  indicates, if purely atomic
C                                    integrals will be evaluated
C                    NUCLEI       =  # of nuclear attraction centers
C                    XN,YN,ZN     =  the x,y,z-coordinates for all
C                                    nuclear attraction centers
C                    NUCCEN       =  contains those index labels of
C                                    the nuclear attraction centers
C                                    that survived the screening process
C                    Px           =  current MIJ coordinates x=X,Y,Z
C                                    for the gaussian product centers
C                                    P=A+B
C                    PAx          =  current MIJ coordinate x=X,Y,Z
C                                    differences P-A between centers
C                                    P and A
C                    PINVHF       =  current MIJ values of 1/(2*P),
C                                    where P are the exponent sums for
C                                    contraction shells A and B
C                    RTS          =  all current MGIJCEN quadrature
C                                    roots
C
C                  Output:
C
C                    R1x          =  the VRR R1-coefficients (individual
C                                    cartesian components x=X,Y,Z) for
C                                    shell expansion on center P
C                    R2           =  the coordinate independent VRR
C                                    R2-coefficients
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C             ...include files and declare variables.
C
C
         IMPLICIT  NONE

         LOGICAL   ATOMIC

         INTEGER   IJ,N
         INTEGER   IXC
         INTEGER   MIJ,MGIJCEN
         INTEGER   NC,NCEN
         INTEGER   NG,NGQP
         INTEGER   NUCLEI

         INTEGER   NUCCEN (1:NCEN)

         DOUBLE PRECISION  PXIJ,PYIJ,PZIJ
         DOUBLE PRECISION  PAXIJ,PAYIJ,PAZIJ
         DOUBLE PRECISION  ROOT
         DOUBLE PRECISION  TWOP
         DOUBLE PRECISION  XC,YC,ZC

         DOUBLE PRECISION  XN     (1:NUCLEI)
         DOUBLE PRECISION  YN     (1:NUCLEI)
         DOUBLE PRECISION  ZN     (1:NUCLEI)

         DOUBLE PRECISION  PX     (1:MIJ)
         DOUBLE PRECISION  PY     (1:MIJ)
         DOUBLE PRECISION  PZ     (1:MIJ)
         DOUBLE PRECISION  PAX    (1:MIJ)
         DOUBLE PRECISION  PAY    (1:MIJ)
         DOUBLE PRECISION  PAZ    (1:MIJ)
         DOUBLE PRECISION  PINVHF (1:MIJ)

         DOUBLE PRECISION  R1X    (1:MGIJCEN)
         DOUBLE PRECISION  R1Y    (1:MGIJCEN)
         DOUBLE PRECISION  R1Z    (1:MGIJCEN)
         DOUBLE PRECISION  R2     (1:MGIJCEN)
         DOUBLE PRECISION  RTS    (1:MGIJCEN)
C
C
C------------------------------------------------------------------------
C
C
C             ...calculate the 1D coefficients.
C
C
         IF (ATOMIC) THEN
             N = 0
             DO IJ = 1,MIJ
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                TWOP = PINVHF (IJ)

                DO NC = 1,NCEN
                   IXC = NUCCEN (NC)
                   XC = XN (IXC)
                   YC = YN (IXC)
                   ZC = ZN (IXC)
                   DO NG = 1,NGQP
                      N = N + 1
                      ROOT  = RTS (N)
                      R1X (N) = ROOT * (XC - PXIJ)
                      R1Y (N) = ROOT * (YC - PYIJ)
                      R1Z (N) = ROOT * (ZC - PZIJ)
                      R2  (N) = TWOP - (TWOP * ROOT)
                   END DO
                END DO
             END DO
         ELSE
             N = 0
             DO IJ = 1,MIJ
                PXIJ = PX (IJ)
                PYIJ = PY (IJ)
                PZIJ = PZ (IJ)
                PAXIJ = PAX (IJ)
                PAYIJ = PAY (IJ)
                PAZIJ = PAZ (IJ)
                TWOP = PINVHF (IJ)

                DO NC = 1,NCEN
                   IXC = NUCCEN (NC)
                   XC = XN (IXC)
                   YC = YN (IXC)
                   ZC = ZN (IXC)
                   DO NG = 1,NGQP
                      N = N + 1
                      ROOT  = RTS (N)
                      R1X (N) = PAXIJ + ROOT * (XC - PXIJ)
                      R1Y (N) = PAYIJ + ROOT * (YC - PYIJ)
                      R1Z (N) = PAZIJ + ROOT * (ZC - PZIJ)
                      R2  (N) = TWOP - (TWOP * ROOT)
                   END DO
                END DO
             END DO
         END IF
C
C
C             ...ready!
C
C
         RETURN
         END
