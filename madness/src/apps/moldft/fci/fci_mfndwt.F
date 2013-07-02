*
* $Id: fci_mfndwt.F,v 1.2 1999/07/28 00:23:40 d3e129 Exp $
*
      subroutine fci_mfind_weight(iwalk,iocc,ns)
      implicit none
#include "fciinfo.fh"
c     
c     find the lexical index of a path on a primitive spin graph
c     for ns electrons and multi=2sz+1, given the spin occupation
c     
      integer iwalk             ! [output]
      integer iocc(*)           ! [input]
      integer ns                ! [input]
c
      integer iwgt, irow, nspt, is
c
      iwgt = 1
      irow = 1
      nspt = ns/2 + 1
      do 10 is = 1,ns
         iwgt = iwgt + marcwt(iocc(is),irow,is,nspt)
         irow = marc(iocc(is),irow,is,nspt)
         if(irow.eq.0) call fci_err('invalid arc in mfind_weight',iwgt)
 10   continue
      if(irow.ne.multi+(ns-multi+1)/2)
     &     call fci_err('invalid tail in mfind_weight',irow)
      iwalk=iwgt
c     
      end
