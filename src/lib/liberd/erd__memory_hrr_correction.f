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
         SUBROUTINE  ERD__MEMORY_HRR_CORRECTION(
     *                     angmom, nshells, spheric,
     *                     imax_hrr, zmax_hrr)
c---------------------------------------------------------------------------
c   Uses shell angular momentum information to provide an estimate of 
c   the maximum possible memory necessary for the HRR transformation.
c   This values for integer memory and double precision memory needed 
c   are returned in the variables imax_hrr, zmax_hrr respectively.
c---------------------------------------------------------------------------
      implicit none
      integer nshells  
      integer angmom(nshells)
      integer imax_hrr, zmax_hrr
      logical spheric

      integer iunique(nshells)
      integer i, j, nunique
      integer m, n, r, s
      integer ncolhrr, nrothrr
      integer shell1, shell2, shell3, shell4
      integer shella, shellb, shellc, shelld
      integer shellp, shellq
      integer shellg, shellh
      integer ngh, ngho
      integer ncol, nrot, nrow
      integer nxyze, nxyzf
      integer nxyzhrr
      integer nxyz1, nxyz2, nxyz3, nxyz4
      integer nry1, nry2, nry3, nry4
      integer nhrr1st, nhrr2nd
      integer nxyzet, nxyzft
      integer nxyzg, nxyzh, nxyzi, nxyzp, nxyzq
      integer nxyzgo, nxyzho

      integer add(0:2)

      logical tr1234, swap12, swap34

      data add /0,0,1/

c---------------------------------------------------------------------------
c   Form an array of unique angular momenta from the shell angular momenta.
c---------------------------------------------------------------------------

      iunique(1) = angmom(1)
      nunique = 1
      do 10 i = 2, nshells
         n = nunique

c--------------------------------------------------------------------------
c   Search the current unique list to see if angmom(i) matches any.
c--------------------------------------------------------------------------

         do j = 1, n
            if (angmom(i) .eq. iunique(j)) go to 10
         enddo 

c--------------------------------------------------------------------------
c   Add this value to the unique list.
c--------------------------------------------------------------------------

         nunique = nunique + 1
         iunique(nunique) = angmom(i)
   10 continue
 
c--------------------------------------------------------------------------
c   Loop over all quadruples of unique angular momenta, calculating the
c   max. value of NCOLHRR and NROTHRR as we go.
c--------------------------------------------------------------------------

      ncolhrr = 0
      nrothrr = 0 

      do m = 1, nunique
      do n = 1, nunique
      do r = 1, nunique
      do s = 1, nunique
         shell1 = iunique(m)
         shell2 = iunique(n)
         shell3 = iunique(r)
         shell4 = iunique(s)
C
C
C             ...decide on the 1 <-> 2 and/or 3 <-> 4 swapping.
C
C
         SWAP12 = SHELL1.LT.SHELL2
         SWAP34 = SHELL3.LT.SHELL4

         SHELLP = SHELL1 + SHELL2
         SHELLQ = SHELL3 + SHELL4

         NXYZ1  = (SHELL1+1)*(SHELL1+2)/2
         NXYZ2  = (SHELL2+1)*(SHELL2+2)/2
         NXYZ3  = (SHELL3+1)*(SHELL3+2)/2
         NXYZ4  = (SHELL4+1)*(SHELL4+2)/2

         NRY1 = SHELL1 + SHELL1 + 1
         NRY2 = SHELL2 + SHELL2 + 1
         NRY3 = SHELL3 + SHELL3 + 1
         NRY4 = SHELL4 + SHELL4 + 1

         IF (.NOT.SPHERIC) THEN
             NRY1 = NXYZ1
             NRY2 = NXYZ2
             NRY3 = NXYZ3
             NRY4 = NXYZ4
         END IF

c---------------------------------------------------------------------------
c   Convert shell1, shell2, ... into shella, shellb, etc.
c---------------------------------------------------------------------------

         SHELLA = MAX0 (SHELL1,SHELL2)
         SHELLB = MIN0 (SHELL1,SHELL2)
         SHELLC = MAX0 (SHELL3,SHELL4)
         SHELLD = MIN0 (SHELL3,SHELL4)

         NXYZE  =   ((SHELLP+1)*(SHELLP+2)*(SHELLP+3)/6)
     +            - ((SHELLA  )*(SHELLA+1)*(SHELLA+2)/6)
         NXYZF  =   ((SHELLQ+1)*(SHELLQ+2)*(SHELLQ+3)/6)
     +            - ((SHELLC  )*(SHELLC+1)*(SHELLC+2)/6)

         IF (SHELLB.EQ.0 .AND. SHELLD.EQ.0) THEN
             NXYZHRR = NXYZE * NXYZF
             TR1234 = .FALSE.
         ELSE
             NHRR1ST = MAX0 (NXYZE*NXYZ3*NXYZ4,NXYZ1*NXYZ2*NRY3*NRY4)
             NHRR2ND = MAX0 (NXYZF*NXYZ1*NXYZ2,NXYZ3*NXYZ4*NRY1*NRY2)
             NXYZHRR = MIN0 (NHRR1ST,NHRR2ND)
             TR1234  = NHRR1ST .GT. NHRR2ND
         END IF

         IF (.NOT.TR1234) THEN

             NXYZET = NXYZE
             NXYZFT = NXYZF
             IF (SWAP12) THEN
                 SHELLA = SHELL2
                 SHELLB = SHELL1
             ELSE
                 SHELLA = SHELL1
                 SHELLB = SHELL2
             ENDIF 

             IF (.NOT.SWAP34) THEN
                 SHELLC = SHELL3
                 SHELLD = SHELL4
             ELSE
                 SHELLC = SHELL4
                 SHELLD = SHELL3
             ENDIF
         ELSE

             NXYZET = NXYZF
             NXYZFT = NXYZE
             IF (.NOT.SWAP12) THEN
                 SHELLC = SHELL1
                 SHELLD = SHELL2
             ELSE
                 SHELLC = SHELL2
                 SHELLD = SHELL1
             ENDIF 

             IF (.NOT.SWAP34) THEN
                 SHELLA = SHELL3
                 SHELLB = SHELL4
             ELSE
                 SHELLA = SHELL4
                 SHELLB = SHELL3
             ENDIF
         ENDIF

         SHELLP = SHELLA + SHELLB
         SHELLQ = SHELLC + SHELLD

         NXYZP  = (SHELLP+1)*(SHELLP+2)/2
         NXYZQ  = (SHELLQ+1)*(SHELLQ+2)/2

C
C             ...if HRR contractions are to be performed, calculate
C                NCOLHRR (maximum # of HRR rotation matrix columns
C                needed to generate the final HRR rotation matrices) and
C                NROTHRR (maximum # of HRR rotation matrix elements).
C
C                First find maximum values for the HRR on the AB-part.
C
C
         IF (SHELLB.NE.0) THEN

             NGH = NXYZET
             NXYZG = NXYZET
             NXYZH = 1
             NXYZGO = NXYZET
             NXYZHO = 1
             NXYZI = NXYZP
             SHELLG = SHELLP
             NROW = 1
             NCOL = NGH
             NROT = NGH

             DO 100 SHELLH = 1,SHELLB
                NXYZGO = NXYZGO - NXYZI
                NXYZHO = NXYZHO + SHELLH + 1
                NGHO = NXYZGO * NXYZHO

c---------------------------------------------------------------------------
c   Use max. possible value for NROW, since we do not know atomic 
c   coordinates of each shell.
c---------------------------------------------------------------------------

                    J = 1 + SHELLH/3
                    NROW = NROW +  J*(J + ADD (MOD(SHELLH,3)))

                NCOL = MAX0 (NGHO,NCOL)
                NROT = MAX0 (NROW*NGHO,NROT)

                NGH = NGHO
                NXYZG = NXYZGO
                NXYZH = NXYZHO
                NXYZI = NXYZI - SHELLG - 1
                SHELLG = SHELLG - 1
  100        CONTINUE

             NCOLHRR = MAX0(NCOL, NCOLHRR)
             NROTHRR = MAX0(NROT, NROTHRR)
             
         END IF
C
C
C             ...next find maximum values for the HRR on the CD-part
C                and set overall maximum values.
C
C
         IF (SHELLD.NE.0) THEN

             NGH = NXYZFT
             NXYZG = NXYZFT
             NXYZH = 1
             NXYZGO = NXYZFT
             NXYZHO = 1
             NXYZI = NXYZQ
             SHELLG = SHELLQ
             NROW = 1
             NCOL = NGH
             NROT = NGH

             DO 200 SHELLH = 1,SHELLD
                NXYZGO = NXYZGO - NXYZI
                NXYZHO = NXYZHO + SHELLH + 1
                NGHO = NXYZGO * NXYZHO

                    J = 1 + SHELLH/3
                    NROW = NROW +  J*(J + ADD (MOD(SHELLH,3)))

                NCOL = MAX0 (NGHO,NCOL)
                NROT = MAX0 (NROW*NGHO,NROT)

                NGH = NGHO
                NXYZG = NXYZGO
                NXYZH = NXYZHO
                NXYZI = NXYZI - SHELLG - 1
                SHELLG = SHELLG - 1
  200        CONTINUE

             NCOLHRR = MAX0 (NCOL,NCOLHRR)
             NROTHRR = MAX0 (NROT,NROTHRR)

         END IF

      enddo
      enddo
      enddo
      enddo

c----------------------------------------------------------------------------
c   The amount of memory for the T array in ERD__HRR_MATRIX is 2*NROTHRR.
c   The amount for the ROW array is also 2*NROTHRR.
c----------------------------------------------------------------------------

      imax_hrr = 2*nrothrr
      zmax_hrr = 2*nrothrr

      RETURN
      END
