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
         SUBROUTINE  ERD__TRANSPOSE_BATCH
     +
     +                    ( NROW,NCOL,
     +                      TILE,
     +                      BATCH,
     +
     +                              OBATCH )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__TRANSPOSE_BATCH
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : none
C  DESCRIPTION : This routine transposes a batch of NROW x NCOL
C                integrals to an output batch of NCOL x NROW integrals.
C                The routine uses tiling for better use of cache.
C
C
C                  Input:
C
C                    NROW,NCOL    =  # of rows and columns in input
C                                    integral batch
C                    TILE         =  tile size
C                    BATCH        =  input integral batch
C
C                  Output:
C
C                    OBATCH       =  output integral batch
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

         INTEGER     I,J
         INTEGER     IEND,JEND
         INTEGER     II,JJ
         INTEGER     NROW,NCOL
         INTEGER     TILE

         DOUBLE PRECISION  BATCH   (1:NROW,1:NCOL)
         DOUBLE PRECISION  OBATCH  (1:NCOL,1:NROW)
C
C
C------------------------------------------------------------------------
C
C
C             ...proceed.
C
C
         DO 100 JJ = 1,NROW,TILE
            JEND = MIN (NROW,JJ+TILE-1)
            DO 110 II = 1,NCOL,TILE
               IEND = MIN (NCOL,II+TILE-1)

               DO 120 J = JJ,JEND
               DO 120 I = II,IEND
                  OBATCH (I,J) = BATCH (J,I)
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
