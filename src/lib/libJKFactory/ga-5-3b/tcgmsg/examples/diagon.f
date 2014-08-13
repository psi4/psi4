      subroutine orthv2(n, v, s, work1, work2)
C$Id: diagon.f,v 1.2 1995-02-02 23:24:07 d3g681 Exp $
      implicit none
      integer n
      double precision v(n,n), s(n,n), work1(n), work2(n)
c
      integer i
      double precision a, ddot
c     
c     orthonormalize vectors (v) in place over metric (s)
c     
c     v = matrix of column vectors
c     s = full square metric matrix
c
c     note ... is not suitable for high precision or if the
c              input vectors are close to linear dependence
c
c          ... assumes that the metric is symmetric
c     
c     normalize vector one
c
      call dgemv('n',n,n,1.0d0,s,n,v(1,1),1,0.0d0,work1,1)
      a = ddot(n, work1, 1, v(1,1), 1)
      call dscal(n, 1.0d0/dsqrt(a), v(1,1), 1)
c
      do 10 i = 2, n
c
c     orthog vector i to vectors j=1,...,i-1
c
         call dgemv('n',n,n,1.0d0,s,n,v(1,i),1,0.0d0,work1,1)
         call dgemv('t',n,i-1,1.0d0,v,n,work1,1,0.0d0,work2,1)
         call dgemv('n',n,i-1,-1.0d0,v,n,work2,1,1.0d0,v(1,i),1)
c
c     normalize vector i
c
         call dgemv('n',n,n,1.0d0,s,n,v(1,i),1,0.0d0,work1,1)
         a = ddot(n, work1, 1, v(1,i), 1)
         call dscal(n, 1.0d0/dsqrt(a), v(1,i), 1)
 10   continue
c
      end
      subroutine jacobi(a,iky,newbas,q,ilifq,nrow,e,iop1,
     *iop2,thresh)
      implicit real*8    (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(*),q(*),e(*),iky(*),ilifq(*)
      common/blkin/ppp(2048),omask(1024),iipt(1024),ipt(1024)
      if(iop1.eq.1) go to 88
      do 44 i=1,newbas
      ilifi=ilifq(i)
      call vclr(q(ilifi+1),1,nrow)
   44 q(ilifi+i)=1.0d0
   88 if(newbas.eq.1)goto 66
   86 te=0.0d0
      do 1 i=2,newbas
      i1=i-1
      ikyi=iky(i)
      do 1 j=1,i1
      temp= dabs(a(j+ikyi))
      if(te.lt.temp)te=temp
    1 continue
      if(te.lt.thresh)goto 99
      te=te*0.1d0
      do 22 i=2,newbas
      i1=i-1
      ip1=i+1
      ikyi=iky(i)
      itest=newbas-ip1
      ii=i+ikyi
      ilifi=ilifq(i)
      do 22 j=1,i1
      ij=j+ikyi
      vij=a(ij)
      if( dabs(vij) .lt. te) go to 22
      vii=a(ii)*0.5d0
      j1=j-1
      jp1=j+1
      ikyj=iky(j)
      jj=j+ikyj
      vjj=a(jj)*0.5d0
      temp=vii-vjj
      tem=dsqrt(temp*temp+vij*vij)
      if(temp)78,77,77
   78 tem=-tem
   77 cost=(temp+tem)/vij
      sint=dsqrt(1.0d0/(1.0d0+cost*cost))
      cost=cost*sint
      temp=vii+vjj
      a(ii)=temp+tem
      a(jj)=temp-tem
      a(ij)=0.0d0
      if(j1.le.0)go to 3
      call drot(j1,a(ikyi+1),1,a(ikyj+1),1,cost,sint)
    3 if(i1 .lt. jp1) go to 5
      do 6 k=jp1,i1
      jj=iky(k)+j
      vij=a(k+ikyi)
      a(k+ikyi)=vij*cost+a(jj)*sint
    6 a(jj)=a(jj)*cost-vij*sint
    5 if(itest)7,79,79
   79 do 8 k=ip1,newbas
      ij=iky(k)+i
      jj=j+iky(k)
      vij=a(ij)
      a(ij)=vij*cost+a(jj)*sint
    8 a(jj)=a(jj)*cost-vij*sint
    7 continue
      call drot(nrow,q(ilifi+1),1,q(ilifq(j)+1),1,cost,sint)
   22 continue
      goto 86
   99 do 11 i=1,newbas
      omask(i)=.false.
      e(i)=a(iky(i)+i)
   11 iipt(i)=i/2
      goto (67,55,55,43),iop2
c... binary sort of e.values to increasing value sequence
   55 ipt(1)=1
      do 19 j=2,newbas
      ia=1
      ib=j-1
      test=e(j)
   53 irm1=ib-ia
      if(irm1)58,50,51
   51 ibp=ia+iipt(irm1)
      if(test.lt.e(ipt(ibp)))goto 52
c...  insert into high half
      ia=ibp+1
      goto 53
c... insert into low half
   52 jj=ib
      do 54 i=ibp,ib
      ipt(jj+1)=ipt(jj)
   54 jj=jj-1
      ib=ibp-1
      goto 53
c...  end point of search
   50 jj=ipt(ia)
      if(test.ge.e(jj))goto 57
      ipt(ia+1)=jj
   58 ipt(ia)=j
      goto 19
   57  ipt(ia+1)=j
   19  continue
       goto (67,68,69,43),iop2
c...   sort by decreasing e.value(invert order)
   69 itest=newbas+1
      ip1=iipt(newbas)
      do 41    i=1,ip1
      j=itest-i
      k=ipt(i)
      ipt(i)=ipt(j)
   41 ipt(j)=k
   68 do 20 i=1,newbas
      k=ipt(i)
      iipt(k)=i
   20 ppp(i)=e(k)
   59 continue
      call dcopy(newbas,ppp(1),1,e(1),1)
c...  iipt(i)=k   means column i is to move to posn k
c...   ipt(i)=k   means column k is to move to posn i
      call sortq(q,ilifq,iipt,newbas,nrow)
      go to 67
c ... locking requested
   43 do 31 j=1,newbas
      m=ilifq(j)
      temp=0.0d0
      do 32 i=1,newbas
      vij= dabs(q(i+m))
      if(vij.lt.temp.or.omask(i))goto 32
      temp=vij
      k=i
   32 continue
      iipt(j)=k
      omask(k)=.true.
   31 ppp(k)=e(j)
      goto 59
   66 e(1)=a(1)
 67   return
      end
      subroutine sortq(q,ilifq,iipt,newbas,nrow)
      implicit real*8    (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*),ilifq(*),iipt(*)
      common/blkin/ppp(2048),omask(1024)
      juse=1
      jnext=1025
      do 10 i = 1, newbas
         omask(i) = .false.
 10   continue
      do 21 i=1,newbas
         if(omask(i))goto 21
         j=i
         call dcopy(nrow,q(ilifq(j)+1),1,ppp(juse),1)
c...  start a permutation cycle
 23      m=iipt(j)
         ilifi=ilifq(m)
         call dcopy(nrow,q(ilifi+1),1,ppp(jnext),1)
         call dcopy(nrow,ppp(juse),1,q(ilifi+1),1)
         if(m.eq.i)goto 21
         juse=jnext
         jnext=1026-jnext
         omask(m)=.true.
         j=m
         goto 23
 21   continue
      return
      end
      subroutine vclr(a,incr,n)
      implicit real*8  (a-h,o-z)
      dimension a(*)
      if(incr.eq.1) then
        do 10 loop=1,n
10      a(loop) = 0.0d0
      else
        loopi=1
        do 20 loop=1,n
           a(loopi)=0.0d0
           loopi=loopi+incr
20      continue
      endif
      return
      end
      subroutine zsqua(n, h, s)
      implicit double precision (a-h, o-z)
c
      dimension h(*), s(n,n)
c
c     put lower triangular array h into square matrix s
c     and halve the diagonal
c
      do 10 i = 1,n-1
         do 20 j = i+1,n
            s(j,i) = 0.0d0
 20      continue
 10   continue
c
      ij = 0
      do 30 i = 1,n
         do 40 j = 1,i-1
            s(j,i) = h(ij+j)
 40      continue
         ij = ij + i
         s(i,i) = h(ij) * 0.5d0
 30   continue
c
      end
      subroutine zfold(n, h, s)
      implicit double precision (a-h, o-z)
c
      dimension h(*), s(n,n)
c
c     fold square matrix s into lower triangular array h
c     ... note that the diagonal gets doubled.
c
      ij = 0
      do 10 i = 1,n
         do 20 j = 1,i
            h(ij+j) = s(i,j) + s(j,i)
 20      continue
         ij = ij + i
 10   continue
c
      end
