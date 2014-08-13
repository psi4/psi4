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
         SUBROUTINE  ERD__MOVE_RY
     +
     +                    ( NINTGRL,NINDEX,
     +                      NOTMOVE,MOVE,NRY,
     +                      INDEX,
     +                      TILE,
     +                      X,
     +
     +                             IXOFF,
     +                             Y )
     +
C------------------------------------------------------------------------
C  OPERATION   : ERD__MOVE_RY
C  MODULE      : ELECTRON REPULSION INTEGRALS DIRECT
C  MODULE-ID   : ERD
C  SUBROUTINES : ERD__TRANSPOSE_BATCH
C                ERD__MAP_IJKL_TO_IKJL
C  DESCRIPTION : This routine moves all ry-components located on the
C                far right in array X to a specific position to the left
C                in array Y:
C
C                     X (1,2,3,...,RY) ---> Y ( 1, 2, 3,...)
C                                   |          |  |  |   |
C                                    --->------^--^--^...^
C
C                The part of X that is not moved (i.e. the # of
C                invariant indices to the left in X) has been calculated
C                beforehand and is transmitted through argument.
C
C
C                  Input:
C
C                       NINTGRL    =  total # of integrals
C                       NINDEX     =  total # of possible target places
C                       NOTMOVE    =  inactive # of indices
C                       MOVE       =  # of indices that will be moved
C                       NRY        =  # of ry-components to be moved
C                       INDEX      =  place 1,2,3,... to which the
c                                     ry-components will be placed in
C                                     array Y (must be within the range
C                                     1 =< INDEX =< NINDEX)
C                       TILE       =  the level 1 cache tile
C                       X          =  initial set of integrals
C                       IXOFF (I)  =  NINDEX-fold array indicating total
C                                     # of elements preceeding place
C                                     I = 1,2,3,... before the move
C                                     (with IXOFF(1)=1)
C
C                  Output:
C
C                       IXOFF (I)  =  updated NINDEX-fold array
C                                     indicating total # of elements
C                                     preceeding place I = 1,2,3,...
C                                     after the move
C                       Y          =  final set of integrals
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

         INTEGER    I
         INTEGER    INDEX
         INTEGER    MOVE
         INTEGER    NINDEX
         INTEGER    NINTGRL
         INTEGER    NOTMOVE
         INTEGER    NRY
         INTEGER    TILE

         INTEGER    IXOFF (1:NINDEX)

         DOUBLE PRECISION   X (1:NINTGRL)
         DOUBLE PRECISION   Y (1:NINTGRL)
C
C
C------------------------------------------------------------------------
C
C
C             ...check, if the move is simply a transposition
C                or a more general ijkl -> ikjl move.
C
C
         IF (NOTMOVE.EQ.1) THEN

             CALL  ERD__TRANSPOSE_BATCH
     +
     +                  ( MOVE,NRY,
     +                    TILE,
     +                    X,
     +
     +                           Y )
     +
     +
         ELSE

             CALL  ERD__MAP_IJKL_TO_IKJL
     +
     +                  ( NOTMOVE,MOVE,NRY,1,
     +                    TILE,
     +                    X,
     +
     +                           Y )
     +
     +
         END IF
C
C
C             ...update IXOFF values.
C
C
         DO 10 I = INDEX,NINDEX
            IXOFF (I) = IXOFF (I) * NRY
   10    CONTINUE
C
C
C             ...ready!
C
C
         RETURN
         END
