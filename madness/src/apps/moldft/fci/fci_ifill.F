*
* $Id: fci_ifill.F,v 1.2 1999/07/28 00:23:40 d3e129 Exp $
*
      subroutine ifill(ival,ix,iskip,n)
      implicit none
      integer ival, ix(*), iskip, n, ipt, i
      ipt=1
      do 10 i=1,n
          ix(ipt)=ival
10        ipt=ipt+iskip
      return
      end
