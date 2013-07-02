*
* $Id: fci_err.F,v 1.2 1999/07/28 00:23:40 d3e129 Exp $
*
      subroutine fci_err(string, code)
      implicit none
      character*(*) string
      integer code
c
      write(6,1) string, code
 1    format(/
     $     ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'/
     $     1x, a, 1x, i8/
     $     ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'/
     $     )
c
      stop 99
c
      end
