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
         SUBROUTINE  OED__MAP_IJKL_TO_IKJL
     +
     +                    ( NI,NJ,NK,NL,
     +                      TILE,
     +                      X,
     +
     +                             Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__MAP_IJKL_TO_IKJL
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine switches the medium 2 indices in the array
C                X (I,J,K,L) to produce the target array Y (I,K,J,L).
C                The routine uses tiling for optimum use of cache.
C
C
C                  Input:
C
C                    NI,NJ,NK,NL  =  size of each dimension in both
C                                    4-dimensional arrays X and Y
C                    TILE         =  tile size
C                    X            =  input 4-dimensional array
C
C                  Output:
C
C                    Y            =  output 4-dimensional array
C
C
C  AUTHOR      : Norbert Flocke
C------------------------------------------------------------------------
C
C
C             ...include files and declare variables.
C
C
         IMPLICIT   NONE

         INTEGER    I,J,K,L
         INTEGER    JJ,KK
         INTEGER    JEND,KEND
         INTEGER    NI,NJ,NK,NL
         INTEGER    TILE

         DOUBLE PRECISION   X (1:NI,1:NJ,1:NK,1:NL)
         DOUBLE PRECISION   Y (1:NI,1:NK,1:NJ,1:NL)
C
C
C------------------------------------------------------------------------
C
C
C             ...proceed.
C
C
         DO 100 L = 1,NL

            DO 110 JJ = 1,NJ,TILE
               JEND = MIN (NJ,JJ+TILE-1)
               DO 120 KK = 1,NK,TILE
                  KEND = MIN (NK,KK+TILE-1)

                  DO 130 J = JJ,JEND
                  DO 130 K = KK,KEND
                  DO 130 I = 1,NI
                     Y (I,K,J,L) = X (I,J,K,L)
  130             CONTINUE

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
