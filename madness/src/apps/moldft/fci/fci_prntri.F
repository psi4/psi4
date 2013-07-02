*
* $Id: fci_prntri.F,v 1.2 1999/07/28 00:23:41 d3e129 Exp $
*
      subroutine fci_prntri(a,m)
      implicit real*8 (a-h,o-z)
c     
c     print lower triangular matrix
c     
      character*8 line
      dimension a(*)
      data line /'--------'/,ncolum/6/
    1 format(2x,1x,6(7x,i5))
    2 format(2x,4x,6(4x,a8))
    3 format(2x,i2,2x,6f12.6)
    4 format(/)
      iout=6
      ii=0
      jj=0
 200  ii=ii+1
      jj=jj+1
      kk=ncolum*jj
      nn=kk*(kk+1)/2
      mm=m
      if(m.gt.kk) mm=kk
      write(iout,1) (i,i=ii,mm)
      write(iout,2) (line,i=ii,mm)
      do 101 i=ii,m
         i1=ii+i*(i-1)/2
         i2=i*(i+1)/2
         if(i2.gt.nn) i2=i1+ncolum-1
         write(iout,3) i,(a(j),j=i1,i2)
 101  continue
      if(m.le.kk) go to 201
      write(iout,4)
      ii=kk
      go to 200
 201  return
      end
