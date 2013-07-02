*
* $Id: fci_prspnad.F,v 1.2 1999/07/28 00:23:41 d3e129 Exp $
*
      subroutine fci_prspnad(ns)
      implicit none
#include "fciinfo.fh"
      integer ns                ! [input]
c
c     Print the spin-adaption coefficents for ns open shells
c
      double precision det(ndetmx), csf(ncsfmx)
      integer icsf, idet, mdet, mcsf
c
      mdet = mnsdet(ns)
      mcsf = mnscsf(ns)
      write(6,1) ns
 1    format(/'  Spin-adaption matrix for', i3, ' open-shells '/
     $     '(determinants down columns, CSF across rows)')
      write(6,2) (icsf, icsf=1,mcsf)
 2    format(1x,1000(4x,i4:))
      do idet = 1, mdet
         call dfill(0.0d0, det, 1, mdet)
         det(idet) = 1.0d0
         call fci_spin_adapt(ns, 1, det, csf)
         write(6,3) (csf(icsf), icsf=1,mcsf)
 3       format(1x,1000f8.4)
      enddo
      write(6,*)
c
      end
