*
* $Id: fci_prteop.F,v 1.2 1999/07/28 00:23:42 d3e129 Exp $
*
      subroutine fci_prteop(ionel,nonel)
      implicit integer(a-z)
      integer ionel(5,*)
c
      write(6,1)
1     format(/' one electron coupling coefficients (i j K J value)'/)
      write(6,2) ((ionel(i,j),i=1,5),j=1,nonel)
2     format(3(2i3,2i6,i3,4x))
      return
      end
