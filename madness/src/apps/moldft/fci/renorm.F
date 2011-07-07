*
* $Id: renorm.F,v 1.2 1999/07/28 00:23:43 d3e129 Exp $
*
      subroutine renorm(n,x,xnorm)
      implicit none
      double precision x(*), xnorm, s
      integer n, i
c     
      if (n.gt.0) then
         s = 0.0d0
         do 10 i = 1,n
            s = s + x(i)*x(i)
 10      continue
         if (s.ne.0.0d0) then
            xnorm = sqrt(s)
            s = 1.0d0/xnorm
            do 20 i = 1,n
               x(i) = x(i) * s
 20         continue
         else
            write(6,*) ' null vector in renorm',n
            xnorm = 0.0d0
         endif
      else
         xnorm = 0.0d0
      endif
c     
      return
      end
