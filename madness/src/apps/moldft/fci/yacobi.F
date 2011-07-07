*
* $Id: yacobi.F,v 1.2 1999/07/28 00:23:43 d3e129 Exp $
*
      subroutine yacobi (a,u,n,big,jb,iord)
      implicit real*8(a-h,o-z)
      dimension a(1),u(1),big(1),jb(1)
      data zero/0.0d0/
c
c     rjh ... added iord=1 (increasing evals) iord=2 (decreasing)
c
c 32 bit precision
c      data eps/1.0e-7/
c 64 bit precision
      data eps/1.0d-12/
      data one/1.0d0/
      l(i)=(i*(i-1))/2
c
      u(1)=one
      if (n.eq.1) go to 999
      j=0
      k=0
      n2=n*n
      do 10 i=1,n2
   10 u(i)=zero
      do 30 i=1,n
      j=j+i
      jj=k+i
      k=k+n
      u(jj)=one
c     section to locate biggest off-diagonal element in row i
      im1=i-1
      if (im1.eq.0) go to 30
      big(i)=zero
      do 20 jj=1,im1
      ij=l(i)+jj
      b=abs(a(ij))
      if (b.le.big(i)) go to 20
      big(i)=b
      jb(i)=jj
   20 continue
   30 continue
c     locate rotation pivot
   31 continue
      ibig=2
      do 35 i=2,n
      if (big(i).gt.big(ibig)) ibig=i
   35 continue
      biggst=big(ibig)
      jbig=jb(ibig)
      if (biggst.le.eps) go to 999
c     begin rotation
      ii=l(ibig)+ibig
      ij=l(ibig)+jbig
      jj=l(jbig)+jbig
      ab=a(ii)-a(jj)
      aa=a(ij)+a(ij)
      d=sqrt(ab*ab+aa*aa)
      if (ab.lt.zero) d=-d
      t=aa/(ab+d)
      tsq=t*t
      csq=one/(one+tsq)
      c=sqrt(csq)
      s=c*t
      ab=aa*t
      d=(a(ii)+ab+tsq*a(jj))*csq
      a(jj)=(a(jj)-ab+tsq*a(ii))*csq
      a(ii)=d
      a(ij)=zero
      iii=4
      kix=(ibig-1)*n
      kjx=(jbig-1)*n
      do 60 i=1,n
      ki=kix+i
      kj=kjx+i
      d=c*u(ki)+s*u(kj)
      u(kj)=c*u(kj)-s*u(ki)
      u(ki)=d
      if (i-jbig) 45,54,46
   45 kj=l(jbig)+i
      ki=l(ibig)+i
      go to 50
   46 kj=l(i)+jbig
      if (i-ibig)47,55,48
   47 ki=l(ibig)+i
      iii=2
      if (jb(i).eq.jbig) iii=3
      go to 50
   48 ki=l(i)+ibig
      iii=1
      if (jb(i).eq.ibig.or.jb(i).eq.jbig) iii=3
   50 d=c*a(ki)+s*a(kj)
      a(kj)=c*a(kj)-s*a(ki)
      a(ki)=d
      go to (51,52,55,60),iii
   51 b=abs(a(ki))
      if (b.le.big(i)) go to 52
      big(i)=b
      jb(i)=ibig
   52 b=abs(a(kj))
      if (b.le.big(i)) go to 60
      big(i)=b
      jb(i)=jbig
      go to 60
   54 if (i.eq.1) go to 60
   55 im1=i-1
      big(i)=zero
      do 58 j=1,im1
      ij=l(i)+j
      b=abs(a(ij))
      if (b.le.big(i)) go to 58
      big(i)=b
      jb(i)=j
   58 continue
   60 continue
c     end of rotation
      go to 31
  999 continue
c lock eigenvectors so that first coefficient is always
c positive.
      ij=1
      do 501 i=1,n
      if (u(ij).ge.zero) go to 503
      do 502 j=1,n
      u(ij+j-1)=-u(ij+j-1)
  502 continue
  503 ij=ij+n
  501 continue
c put eigen values in big and order them
      ii = 0
      do 1231 i = 1,n
          ii = ii+i
          big(i) = a(ii)
1231  continue
      call order(n,big,u,iord)
      return
      end
      subroutine order(n,e,v,iord)
      implicit real*8(a-h,o-z)
      dimension e(n),v(n,n)
      logical oswap
c
c order eigenvalues and vectors in increasing order.
c
      do 10 j=1,n
          do 20 i=j+1,n
             if (iord .eq. 1) then
                oswap = e(i).lt.e(j)
             else
                oswap = e(j).lt.e(i)
             endif
              if(oswap) then
                t=e(i)
                e(i)=e(j)
                e(j)=t
                do 30 k=1,n
                    t=v(k,i)
                    v(k,i)=v(k,j)
30                  v(k,j)=t
              endif
20        continue
10    continue
      return
      end
