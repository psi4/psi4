*
* $Id: fci_owalk_info.F,v 1.2 1999/07/28 00:23:41 d3e129 Exp $
*
      subroutine fci_owalk_info(owalk, iocc, ns, ioff, numdet, sym)
      implicit none
#include "fciinfo.fh"
c     
c     Given the index of an orbital walk return the occupation, the 
c     offset in the CI determinant vector, the number of open shells
c     and the corresponding number of determinants, and the 
c     symmetry of the occupation.  
c     
c     The CI coefficients for this orbital occupation are addressed
c     as C(ioff+i) i=1,...,numdet.
c     
c     This routine bundles the function of several other routines
c     as a convenience.
c     
      integer owalk             ! [input]  Index of orbital occupation
      integer iocc(*)           ! [output] The orbital occupation pattern
      integer ns                ! [output] No. of open shells
      integer ioff              ! [output] Offset into CI vector
      integer numdet            ! [output] No. of determinants for this ns
      integer sym               ! [output] Symmetry 0...7
c     
      integer i
#include "bitops.fh"
c     
      if (owalk.lt.1 .or. owalk.gt.nwalko) call fci_err
     $     ('fci_owalk_info: bad owalk ', owalk)
c     
      call fci_ofind_occ(owalk, iocc)
      ns = 0
      sym = 0
      do i = 1, nactiv
         if (iocc(i) .eq. 2) then 
            ns = ns + 1
            sym = ieor(sym,orbsym(i))
         endif
      enddo
      numdet = mnsdet(ns)
      if (sym .eq. statesym) then
         ioff = indexo(owalk)
      else
         ioff = -999999999
      endif
c     
      end
