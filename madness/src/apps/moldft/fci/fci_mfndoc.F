*
* $Id: fci_mfndoc.F,v 1.2 1999/07/28 00:23:40 d3e129 Exp $
*
      subroutine fci_mfind_occ(iwalk,iocc,ns)
      implicit none
#include "fciinfo.fh"
c     
c     given the lexical index of a path on a primitive spin graph
c     for ns electrons and multi=2sz+1, return the spin occupation
c     
      integer iwalk             ! [input]
      integer iocc(*)           ! [output]
      integer ns                ! [input]
c
      integer iwgt, irow, nspt, is, k, krow, itest
c     
      iwgt = iwalk-1
      irow = multi+(ns-multi+1)/2
      nspt = ns/2 + 1
      do 10 is = ns+1,2,-1
         do 20 k=1,2
            krow = irow-(2-k)
            if(marc(k,krow,is-1,nspt).ne.0) then
               itest=iwgt-marcwt(k,krow,is-1,nspt)
               if(itest.ge.0) then
                  iwgt = itest
                  irow = krow
                  iocc(is-1) = k
                  goto 10
               endif
            endif
 20      continue
         call fci_err('no arc found in mfind_occ',is)
 10   continue
      if(iwgt.ne.0) call fci_err('invalid head in mfind_occ',iwgt)
c
      end
