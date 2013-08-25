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
         SUBROUTINE  ERD__GENER_ERI_DERV_BATCH
     +
     +                    ( IMAX,ZMAX,
     +                      NALPHA,NCOEFF,NCSUM,
     +                      NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                      NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                      SHELL1,SHELL2,SHELL3,SHELL4,
     +                      X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                      DER1X,DER1Y,DER1Z,
     +                      DER2X,DER2Y,DER2Z,
     +                      DER3X,DER3Y,DER3Z,
     +                      DER4X,DER4Y,DER4Z,
     +                      ALPHA,CC,CCBEG,CCEND,
     +                      SPHERIC,
     +                      SCREEN,
     +                      ICORE,
     +
     +                                NBATCH,
     +                                NFIRST,
     +                                ZCORE )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__GENER_ERI_DERV_BATCH
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : ERD__DERV_CSGTO
C  DESCRIPTION : Main operation that drives the calculation of a batch
C                of derivated contracted electron repulsion integrals.
C
C
C                  Input (x = 1,2,3 and 4):
C
C                    IMAX,ZMAX    =  maximum integer,flp memory
C                    NALPHA       =  total # of exponents
C                    NCOEFF       =  total # of contraction coeffs
C                    NCSUM        =  total # of contractions
C                    NCGTOx       =  # of contractions for csh x
C                    NPGTOx       =  # of primitives per contraction
C                                    for csh x
C                    SHELLx       =  the shell type for csh x
C                    Xy,Yy,Zy     =  the x,y,z-coordinates for centers
C                                    y = 1,2,3 and 4
C                    DERyp        =  the order of differentiation on
C                                    centers y = 1,2,3,4 with respect
C                                    to the p = x,y,z coordinates
C                    ALPHA        =  primitive exponents for csh
C                                    1,2,3,4 in that order
C                    CC           =  contraction coefficient for csh
C                                    1,2,3,4 in that order, for each
C                                    csh individually such that an
C                                    (I,J) element corresponds to the
C                                    I-th primitive and J-th contraction.
C                    CC(BEG)END   =  (lowest)highest nonzero primitive
C                                    index for contractions for csh
C                                    1,2,3,4 in that order. They are
C                                    different from (1)NPGTOx only for
C                                    segmented contractions
C                    SPHERIC      =  is true, if spherical integrals
C                                    are wanted, false if cartesian
C                                    ones are wanted
C                    SCREEN       =  is true, if screening will be
C                                    done at primitive integral level
C                    ICORE        =  integer scratch space
C                    ZCORE (part) =  flp scratch space
C
C
C                  Output:
C
C                    NBATCH       =  # of derivative integrals in batch
C                    NFIRST       =  first address location inside the
C                                    ZCORE array containing the first
C                                    derivative integral
C                    ZCORE        =  full batch of contracted (12|34)
C                                    derivative integrals over spherical
C                                    or cartesian gaussians starting
C                                    at ZCORE (NFIRST)
C
C
C
C                            !!! IMPORTANT !!!
C
C                For performance tuning, please see the include file
C                'erd__tuning.inc'.
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

         INCLUDE     'erd__ftable.inc'
         INCLUDE     'erd__tuning.inc'

         LOGICAL     SCREEN
         LOGICAL     SPHERIC

         INTEGER     DER1X,DER2X,DER3X,DER4X
         INTEGER     DER1Y,DER2Y,DER3Y,DER4Y
         INTEGER     DER1Z,DER2Z,DER3Z,DER4Z
         INTEGER     IMAX,ZMAX
         INTEGER     NALPHA,NCOEFF,NCSUM
         INTEGER     NBATCH,NFIRST
         INTEGER     NCGTO1,NCGTO2,NCGTO3,NCGTO4
         INTEGER     NPGTO1,NPGTO2,NPGTO3,NPGTO4
         INTEGER     SHELL1,SHELL2,SHELL3,SHELL4

         INTEGER     CCBEG (1:NCSUM)
         INTEGER     CCEND (1:NCSUM)
         INTEGER     ICORE (1:IMAX)

         DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4

         DOUBLE PRECISION  ALPHA (1:NALPHA)
         DOUBLE PRECISION  CC    (1:NCOEFF)
         DOUBLE PRECISION  ZCORE (1:ZMAX)
C
C
C------------------------------------------------------------------------
C
C
C             ...call main routine.
C
C
         CALL  ERD__DERV_CSGTO
     +
     +              ( IMAX,ZMAX,
     +                NALPHA,NCOEFF,NCSUM,
     +                NCGTO1,NCGTO2,NCGTO3,NCGTO4,
     +                NPGTO1,NPGTO2,NPGTO3,NPGTO4,
     +                SHELL1,SHELL2,SHELL3,SHELL4,
     +                X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,
     +                DER1X,DER1Y,DER1Z,
     +                DER2X,DER2Y,DER2Z,
     +                DER3X,DER3Y,DER3Z,
     +                DER4X,DER4Y,DER4Z,
     +                ALPHA,CC,CCBEG,CCEND,
     +                FTABLE,MGRID,NGRID,TMAX,TSTEP,TVSTEP,
     +                L1CACHE,TILE,NCTROW,
     +                SPHERIC,SCREEN,
     +                ICORE,
     +
     +                          NBATCH,
     +                          NFIRST,
     +                          ZCORE )
     +
     +
C
C
C             ...ready!
C
C
         RETURN
         END
