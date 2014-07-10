#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      double precision function exprjh(x)
C$Id: integ.F,v 1.1 2005-03-08 23:58:03 d3g293 Exp $
      double precision x
c     
c     dumb solution to underflow problems on sun
c
      if (x.lt.-37.0d0) then
         exprjh = 0.0d0
      else
         exprjh = exp(x)
      endif
      end
      subroutine setfm
      implicit double precision (a-h,o-z)
      common/values/fm(2001,5),rdelta,delta,delo2
      dimension t(2001),et(2001)
c
c     initalize common block for computation of f0 by recursion down
c     from f200
c
      delta=28.0d0/2000.0d0
      delo2=delta*0.5d0
      rdelta=1.0d0/delta
      maxm=4
      do 10 i=1,2001
          tt=delta*dble(i-1)
          et(i)=exprjh(-tt)
          t(i)=2.0d0*tt
          fm(i,maxm+1)=0.0d0
10    continue
      do 20 i=200,maxm,-1
          rr=1.0d0/dble(2*i+1)
          do 30 ii=1,2001
              fm(ii,maxm+1)=(et(ii)+t(ii)*fm(ii,maxm+1))*rr
30        continue
20    continue
      do 40 i=maxm,1,-1
          rr=1.0d0/dble(2*i-1)
          do 50 ii=1,2001
            fm(ii,i)=(et(ii)+t(ii)*fm(ii,i+1))*rr
50        continue
40    continue
c
      end
      subroutine f0(value, t)
      implicit real*8   (a-h,o-z)
      common/values/fm(2001,5),rdelta,delta,delo2
      parameter(fac0=0.88622692545276d0,
     $          rhalf=0.5d0,rthird=0.3333333333333333d0,rquart=0.25d0)
      data t0/28.d0/
c
c     computes f0 to a relative accuracy of better than 4.e-13 for all t.
c     uses 4th order taylor expansion on grid out to t=28.0
c     asymptotic expansion accurate for t greater than 28
c
      if(t.ge.t0) then
          value = fac0 / sqrt(t)
      else
          n = idint((t+delo2)*rdelta)
          x = delta*dble(n)-t
          n = n+1
          value = fm(n,1)+x*(fm(n,2)+rhalf*x*(fm(n,3)+
     $             rthird*x*(fm(n,4)+rquart*x*fm(n,5))))
      endif
c
      end
      subroutine addin(g, i, j, k, l, fock, dens, iky)
      implicit double precision (a-h, o-z)
      dimension fock(*), dens(*), iky(*)
c
c     add (ij|kl) into the fock matrix
c
      gg = g
      g2 = gg+gg
      g4 = g2+g2
      ik = iky(i) + k
      il = iky(i) + l
      ij = iky(i) + j
      jk = iky(max(j,k)) + min(j,k)
      jl = iky(max(j,l)) + min(j,l)
      kl = iky(k) + l
      aij = g4*dens(kl)+fock(ij)
      fock(kl) = g4*dens(ij)+fock(kl)
      fock(ij) = aij
      gil=gg
      if(i.eq.k.or.j.eq.l) gg = g2
      if(j.eq.k) gil = g2
      ajk = fock(jk) - gil*dens(il)
      ail = fock(il) - gil*dens(jk)
      aik = fock(ik) - gg*dens(jl)
      fock(jl) = fock(jl) - gg*dens(ik)
      fock(jk) = ajk
      fock(il) = ail
      fock(ik) = aik
c
      end
      subroutine dfill(n,val,a,ia)
      implicit real*8 (a-h,o-z)
      dimension a(*)
c
c     initialise double precision array to scalar value
c
      if (ia.eq.1) then
         do 10 i = 1, n
            a(i) = val
 10      continue
      else
         do 20 i = 1,(n-1)*ia+1,ia
            a(i) = val
 20      continue
      endif
c
      end
