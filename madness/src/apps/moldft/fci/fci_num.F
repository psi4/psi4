*
* $Id: fci_num.F,v 1.2 1999/07/28 00:23:41 d3e129 Exp $
*
      integer function fci_num_dets(ns)
      implicit none
#include "fciinfo.fh"
c
c     return the number of determinants for ns open shell electrons
c     assuming that ns is a valid no. of open shell electrons for
c     the given FCI space and spin multiplicity
c
      integer ns                ! [input]
c
      fci_num_dets = mnsdet(ns)
c
      end
      integer function fci_num_csf(ns)
      implicit none
#include "fciinfo.fh"
c
c     return the number of configurations for ns open shell electrons
c     assuming that ns is a valid no. of open shell electrons for
c     the given FCI space and spin multiplicity
c
      integer ns                ! [input]
c
      fci_num_csf = mnscsf(ns)
c
      end
      integer function fci_num_owalks()
      implicit none
#include "fciinfo.fh"
c
c     Returns the number of orbital occupation walks in the FCI.
c     May need to be extended to accomodate symmetry.
c
      fci_num_owalks = nwalko
c
      end
      integer function fci_ndets()
      implicit none
#include "fciinfo.fh"
c
c     Returns the total number of determinants in the FCI
c     (with use of symmetry)
c
      fci_ndets = ndets
c
      end
      integer function fci_ncsf()
      implicit none
#include "fciinfo.fh"
c
c     Returns the total number of configurations in the FCI
c     (with use of symmetry)
c
      fci_ncsf = ncsf
c
      end
      integer function fci_max_num_dets()
      implicit none
#include "fciinfo.fh"
c
c     Returns the maximum number of determinants associated
c     with any orbital configuration in the full CI.
c
c     [Note that occupations K that arise as a result of
c     .factorizations such as <I|Eij|J>=<I|Eia|K><K|Eaj|J>
c     .where a is a ficitious orbital, may have more open shells.  
c     .However, these are only relevant in the internals of the 
c     .coupling coefficient generation routines]
c
      integer ns
c
      if (nelec.le.nactiv) then
         ns = nelec
      else
         ns = (2*nactiv-nelec)
      endif
c
      if (ns.le.multi .or. ns.gt.nsmax) call fci_err
     $     ('fci_max_num_dets: I am confused!',0)
c
      fci_max_num_dets = mnsdet(ns)
c
      end
