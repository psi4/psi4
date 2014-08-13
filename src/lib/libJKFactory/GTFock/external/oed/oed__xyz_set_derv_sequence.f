C  Copyright (c) 1997-1999, 2003 Massachusetts Institute of Technology
C 
C  This program is free software; you can redistribute it and/or modify
C  it under the terms of the GNU General Public License as published by
C  the Free Software Foundation; either version 2 of the License, or
C  (at your option) any later version.

C  This program is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C  GNU General Public License for more details.

C  You should have received a copy of the GNU General Public License
C  along with this program; if not, write to the Free Software
C  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
C  USA
         SUBROUTINE  OED__XYZ_SET_DERV_SEQUENCE
     +
     +                    ( NDERX,NDERY,NDERZ,
     +                      DERAX,DERAY,DERAZ,
     +                      DERBX,DERBY,DERBZ,
     +                      DERMX,DERMY,DERMZ,             ! Watson Added
     +
     +                               CENSQX,CENSQY,CENSQZ )
     +
C------------------------------------------------------------------------
C  OPERATION   : OED__OVL_SET_DERV_SEQUENCE
C  MODULE      : ONE ELECTRON INTEGRALS DIRECT
C  MODULE-ID   : OED
C  SUBROUTINES : none
C  DESCRIPTION : This routine determines the sequence of centers for
C                the x,y,z-coordinate derivations. The info is stored
C                in arrays CENSQp for each coordinate p=X,Y,Z.
C
C                The indices for the centers obey the following
C                indexing scheme: index 1 -> A and index 2 -> B
C
C                Take the x-coordinate for example and let CENSQX have
C                the following elements on exit from this routine:
C
C                            CENSQX (I) = 2,1,1,2,...
C
C                Then the sequence of single derivations is:
C
C                       dx/dBx  dx/dAx  dx/dAx  dx/dBx  (A|B)
C
C                Note, that the indexing of the CENSQp arrays start
C                with 0. The reason is that the NDERp values might
C                be equal to zero.
C
C                The center sequences are determined such that those
C                on center A come first.
C
C
C                  Input:
C
C                    NDERp        =  the total order of differentiation
C                                    with respect to the p = x,y,z
C                                    coordinates
C                    DERxp        =  the order of differentiation on
C                                    centers x = A and B with respect
C                                    to the p = x,y,z coordinates
C
C                  Output:
C
C                    CENSQp       =  center sequence array for the
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

         INTEGER     DERAX,DERAY,DERAZ
         INTEGER     DERBX,DERBY,DERBZ
         INTEGER     DERMX,DERMY,DERMZ                     ! Watson Added
         INTEGER     I
         INTEGER     NDERX,NDERY,NDERZ

         INTEGER     CENSQX (0:NDERX)
         INTEGER     CENSQY (0:NDERY)
         INTEGER     CENSQZ (0:NDERZ)
C
C
C------------------------------------------------------------------------
C
C
C             ...first the center sequences due to centers A.
C
C
         DO I = 1,DERAX
            CENSQX (I) = 1
         END DO

         DO I = 1,DERAY
            CENSQY (I) = 1
         END DO

         DO I = 1,DERAZ
            CENSQZ (I) = 1
         END DO
C
C
C             ...next the center sequences due to centers B.
C
C
         DO I = 1,DERBX
            CENSQX (DERAX+I) = 2
         END DO

         DO I = 1,DERBY
            CENSQY (DERAY+I) = 2
         END DO

         DO I = 1,DERBZ
            CENSQZ (DERAZ+I) = 2
         END DO
C
C
C             ...next the center sequences due to the moments.
C
C
         DO I = 1,DERMX
            CENSQX (DERBX+I) = 3
         END DO

         DO I = 1,DERMY
            CENSQY (DERBY+I) = 3
         END DO

         DO I = 1,DERMZ
            CENSQZ (DERBZ+I) = 3
         END DO

C
C
C             ...ready!
C
C
         RETURN
         END
