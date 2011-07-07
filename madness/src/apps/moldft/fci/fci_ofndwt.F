*
* $Id: fci_ofndwt.F,v 1.2 1999/07/28 00:23:41 d3e129 Exp $
*
      subroutine fci_ofind_weight(iwalk,iocc)
      implicit none
#include "fciinfo.fh"
c     
c     find the lexical index of a path on orbital graph
c     for nelec electrons and nactiv orbitals
c     
      integer iwalk             ! [output]
      integer iocc(*)           ! [input]
c
      integer iwgt, irow, iorb
c     
      iwgt = 1
      irow = 1
      do 10 iorb = 1,nactiv
         iwgt = iwgt + karcwt(iocc(iorb),irow,iorb)
         irow = karc(iocc(iorb),irow,iorb)
*        if(irow.le.0) call fci_err('invalid arc in ofind_weight',iwgt)
 10   continue
      if(irow.ne.nelec+1)
     &     call fci_err('invalid tail in ofind_weight',irow)
      iwalk=iwgt
c
      end
