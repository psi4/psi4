*
* $Id: fci_ofndoc.F,v 1.2 1999/07/28 00:23:41 d3e129 Exp $
*
      subroutine fci_ofind_occ(iwalk,iocc)
      implicit none
#include "fciinfo.fh"
c     
c     find orbital occupation of walk with lexical index iwalk
c     
      integer iwalk             ! [input]
      integer iocc(*)           ! [output]
c
      integer iwgt, irow, iorb, k, krow, itest
c     
      iwgt = iwalk-1
      irow = nelec+1
      do 10 iorb = nactiv+1,2,-1
         do 20 k=3,1,-1
            krow = irow-(k-1)
            if(krow.gt.0) then
               if(karc(k,krow,iorb-1).ne.0) then
                  itest=iwgt-karcwt(k,krow,iorb-1)
                  if(itest.ge.0) then
                     iwgt = itest
                     irow = krow
                     iocc(iorb-1) = k
                     goto 10
                  endif
               endif
            endif
 20      continue
         call fci_err('no arc found in ofind_occ',iorb)
 10   continue
      if(iwgt.ne.0) call fci_err('invalid head in ofind_occ',iwgt)
c     
      end
