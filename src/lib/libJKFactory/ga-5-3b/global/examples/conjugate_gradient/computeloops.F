#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine f_matvecmul(a_mat,
     >        mult_vec,
     >        resultvec,
     >        isresultmirrored,
     >        firstrow,
     >        lastrow,
     >        myrowlist,
     >        mycollist)
      implicit none
      integer*4 isresultmirrored
      double precision mult_vec(*)
      double precision resultvec(*)
      double precision a_mat(*)
      integer*4 myrowlist(*),mycollist(*)
      integer*4 firstrow,lastrow
      integer*4 i,j,k,l
      double precision tmprowsum
      k=0
      l=1
      do i=firstrow, lastrow-1
        k=k+1
        write(6,*) k,myrowlist(k),myrowlist(k+1)
        do j=myrowlist(k), myrowlist(k+1)-1
          write(6,*) i,j,k,myrowlist(k+1),l
          tmprowsum=tmprowsum+a_mat(l)*mult_vec(mycollist(l)+1)

          l=l+1
        enddo

        resultvec(k)=tmprowsum
        tmprowsum=0

      enddo
      end


      subroutine f_computeminverse(minv,
     >         a_mat,
     >        myrowlist,
     >        mycollist,
     >        firstrow,
     >        lastrow)
      double precision a_mat(*)
      double precision minv(*)
      integer*4 myrowlist(*),mycollist(*)
      integer*4 firstrow,lastrow
      integer*4 i,j,k,l

      k=0
      l=1
      do i=firstrow, lastrow-1
        k=k+1
        do j=myrowlist(k), myrowlist(k+1)-1
          if(mycollist(l) .ge. i) then

            if(mycollist(l) .eq. i) then

              minv(k)=10.0/a_mat(l)
              if(minv(k) .le. 0) minv(k)=1.0/10.0

            endif

            if(mycollist(l)>i) minv(k)=1.0/10.0
              
            l= l+(myrowlist(k+1)-j)
            go to 100
          endif
          l=l+1
        enddo

100     continue
      enddo
      end



      subroutine f_computeminverser(minv,
     >        r,
     >        minvr,
     >        firstrow,
     >        lastrow)
      double precision minv(*)
      double precision r(*), minvr(*)
      integer*4 firstrow,lastrow
      integer*4 i,k
      k=1
      do i=firstrow, lastrow-1
        minvr(k)=minv(k)*r(k)
        k=k+1
      enddo

      end



      subroutine f_addvec(scale1,
     >         srcvec1,
     >         scale2,
     >         srcvec2,
     >         dstvec,
     >         firstrow,
     >         lastrow)
      integer*4 firstrow,lastrow
      double precision scale1,scale2
      double precision srcvec1(*),srcvec2(*)
      double precision dstvec(*)
      integer*4 i,k
      k=1
      do i=firstrow, lastrow-1
        dstvec(k)=srcvec1(k)*scale1+ srcvec2(k)*scale2
        k=k+1
      enddo

      end

      subroutine f_2addvec(scale1a,
     >         srcvec1a,
     >         scale2a,
     >         srcvec2a,
     >         dstveca,
     >         scale1b,
     >         srcvec1b,
     >         scale2b,
     >         srcvec2b,
     >         dstvecb,
     >         firstrow,
     >         lastrow)
      integer*4 firstrow,lastrow
      double precision scale1a,scale2a
      double precision srcvec1a(*),srcvec2a(*)
      double precision dstveca(*)
      double precision scale1b,scale2b
      double precision srcvec1b(*),srcvec2b(*)
      double precision dstvecb(*)
      integer*4 i,k

      k=1
      do i=firstrow, lastrow-1
        dstveca(k)=srcvec1a(k)*scale1a+ srcvec2a(k)*scale2a
        dstvecb(k)=srcvec1b(k)*scale1b+ srcvec2b(k)*scale2b
        k=k+1
      enddo

      end

