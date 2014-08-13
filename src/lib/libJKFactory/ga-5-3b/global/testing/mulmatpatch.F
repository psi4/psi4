#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c $Id: mulmatpatch.F,v 1.6 2005-11-23 10:25:18 manoj Exp $
#if (defined(CRAY) && !defined(__crayx1)) || defined(KSR)
#   define xgemm SGEMM
#   define ygemm CGEMM
#else
#   define xgemm TEST_DGEMM
#   define ygemm TEST_ZGEMM
#endif
#if defined(FUJITSU) || defined(CRAY_YMP)
# define THRESH  1.0d-10
#else
# define THRESH  1.0d-20
#endif
#define MISMATCH(x,y) abs(x-y)/max(1,abs(x)).gt.THRESH
#define NMAX 10
c
c
      program test
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      logical status
      integer bufsize, gasize
c
c***  Initialize a message passing library
c
#include "mp3.fh"
c
      call ga_initialize()
      if(ga_nodeid().eq.0)then
         print *,' GA initialized'
         call ffflush(6)
      endif
c
c     we want to force distribution of innermost loop in nga_mulmat_patch
c     by providing less buffer memory than needed
c
      if(ga_uses_ma()) then
              gasize = (NMAX**4 * 3)/ga_nnodes()
      else
              gasize =0
      endif
      bufsize = (NMAX/2 + 1)*(NMAX/3 + 1)*2 + ((NMAX/2 + 1)**2)
      bufsize = bufsize*6/7
      status = ma_init(MT_DBL, 10, gasize+bufsize+500000)
      if (.not. status)call ga_error( 'ma_init failed', -1)
      if(ga_nodeid().eq.0)then
         print *,' '
         print *,'CHECKING MATRIX MULTIPLICATION FOR PATCHES '
c        print *, 'gasize and bufsize are ', gasize, bufsize
         print *,' '
      endif
      call dpatch_test()
      call dpatch_test2()
c
      if(ga_nodeid().eq.0) print *,'All tests successful '
c      
      call ga_terminate()
c
      call MP_FINALIZE()
      end



      subroutine dpatch_test()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c
      integer n,m, dim, dim1
      parameter (n = NMAX)
      parameter (m = n*n)
      parameter (dim = 4)
      parameter (dim1 = 3)
      double precision a(n,n), b(n,n), c(n,n), buf(m), val
      double precision alpha, beta
      integer nproc, me 
      integer i, j, ailo, ajlo, bilo, bjlo, base, iran
      integer alo(dim), ahi(dim), ald(dim1), chunk(dim)
      integer blo(dim), bhi(dim)
      integer clo(dim), chi(dim)
      integer bufld(dim)
      integer g_a1, g_b1, g_c1 
      integer rows, cols, loop
      logical status, dist_same 
      iran(i) = int(drand(0)*real(i)) + 1
      dist_same = .false.
c
      me = ga_nodeid()
      nproc = ga_nnodes()
c
      do j = 1, n
         do i = 1, n
            a(i,j) = i-1 + (j-1)*n
            b(i,j) = i+j 
         enddo
      enddo
c
      ailo = iran(n/2)
      ajlo = iran(n/2)
      bilo = iran(n/2)
      bjlo = iran(n/2)
c
      do i = 1, dim
         alo(i) = 1
         ahi(i) = n
         chunk(i) = -1
         if (i.ne.dim) ald(i) = n
      enddo
c
      status = nga_create(MT_DBL, dim, ahi, 'a1', chunk, g_a1)
      if (.not. status) then
         write(6,*) ' nga_create failed'
         call ffflush(6)
         call ga_error('... exiting ',0)
      endif
c
      if(dist_same) then
         status = ga_duplicate(g_a1, g_b1, 'a1_duplicated')
         if(.not.ga_compare_distr(g_a1, g_b1))
     $           call ga_error("g_b1 distribution different",0) 
         status = ga_duplicate(g_a1, g_c1, 'a1_duplicated_again')
         if(.not.ga_compare_distr(g_a1, g_c1))
     $           call ga_error("g_c1 distribution different",0) 
      else
         chunk(dim) = n
         status = nga_create(MT_DBL, dim, ahi, 'b1', chunk, g_b1)
         if (.not. status) call ga_error('nga_create failed:b1',0) 
         chunk(dim) = 0
         chunk(dim1) = n
         status = nga_create(MT_DBL, dim, ahi, 'c1', chunk, g_c1)
         if (.not. status) call ga_error('nga_create failed:c1',0) 
      endif
c
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking nga_matmul_patch ... '
         call ffflush(6)
      endif
c
      alo(1) = 1
      ahi(1) = n
      blo(1) = 1
      bhi(1) = n
      bufld(1) = n
      do j = 2, dim-1
        alo(j) = j
        ahi(j) = j
        blo(j) = j
        bhi(j) = j
        clo(j) = j
        chi(j) = j
        bufld(j) = 1
      enddo
      do j = 1+me, n, nproc 
         alo(dim) = j
         ahi(dim) = j
         blo(dim) = j
         bhi(dim) = j
         call nga_put(g_a1,alo,ahi,a(1,j),bufld) 
         call nga_put(g_b1,blo,bhi,b(1,j),bufld)
      enddo
c
      call ga_sync()
      alo(1) = ailo
      ahi(1) = ailo+1
c     ahi(1) = ailo+n/2
      alo(dim) = ajlo
      ahi(dim) = ajlo
c     ahi(dim) = ajlo+n/3
      blo(1) = bilo
      bhi(1) = bilo
c     bhi(1) = bilo+n/3
      blo(dim) = bjlo
      bhi(dim) = bjlo
c     bhi(dim) = bjlo+n/2
      clo(1) = bilo
      chi(1) = bilo+1
c     chi(1) = bilo+n/2
      clo(dim) = bjlo
      chi(dim) = bjlo
c     chi(dim) = bjlo+n/2
c
      call nga_matmul_patch('n','n', 1d0, 0d0, 
     $                      g_a1, alo, ahi,
     $                      g_b1, blo, bhi,
     $                      g_c1, clo, chi)
      call xgemm('n','n',2,2,1,1d0,a(ailo,ajlo), n,
     $            b(bilo,bjlo),n, 0d0, c, n)
c     call xgemm('n','n',n/2+1,n/2+1,n/3+1,1d0,a(ailo,ajlo), n,
c    $            b(bilo,bjlo),n, 0d0, c, n)
c
      bufld(1) = n/2 + 1
      call nga_get(g_c1,clo,chi,buf,bufld)
      base = 0
      do j = 1, 1
         if(Mod(j,nproc).eq.me) then
            do i = 1, 2
               base = base+1
               if(ABS(buf(base)- c(i,j)).gt.1d-8) then
                  print *,me, ' error ',i,j, buf(base), c(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
c
      call ga_sync()
      if (me .eq. 0) then
         write(6,*)'  a*b: OK '
         call ffflush(6)
      endif
c
      call ga_sync()
      alo(1) = ailo
      ahi(1) = ailo+n/2
      alo(dim) = ajlo
      ahi(dim) = ajlo+n/3
      blo(1) = bilo
      bhi(1) = bilo+n/3
      blo(dim) = bjlo
      bhi(dim) = bjlo+n/2
      clo(1) = bilo
      chi(1) = bilo+n/2
      clo(dim) = bjlo
      chi(dim) = bjlo+n/2
c
      call nga_matmul_patch('n','n', 1d0, 0d0, 
     $                      g_a1, alo, ahi,
     $                      g_b1, blo, bhi,
     $                      g_c1, clo, chi)
      call xgemm('n','n',n/2+1,n/2+1,n/3+1,1d0,a(ailo,ajlo), n,
     $            b(bilo,bjlo),n, 0d0, c, n)
c
      bufld(1) = n/2 + 1
      call nga_get(g_c1,clo,chi,buf,bufld)
      base = 0
      do j = 1, 1+n/2
         if(Mod(j,nproc).eq.me) then
            do i = 1, 1+n/2
               base = base+1
               if(ABS(buf(base)- c(i,j)).gt.1d-8) then
                  print *,me, ' error ',i,j, buf(base), c(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
c
      call ga_sync()
      if (me .eq. 0) then
         write(6,*)'  a*b: OK '
         call ffflush(6)
      endif
c
      call ga_sync()
      call nga_matmul_patch('t','n', 1d0, 0d0,
     $                      g_a1, alo, ahi,
     $                      g_b1, blo, bhi,
     $                      g_c1, clo, chi)
      call xgemm('t','n',n/2+1,n/2+1,n/3+1,1d0,a(ajlo,ailo), n,
     $            b(bilo,bjlo),n, 0d0, c, n)
      call nga_get(g_c1,clo,chi,buf,bufld)
      base = 0
      do j = 1, 1+n/2
         if(Mod(j,nproc).eq.me) then
            do i = 1, 1+n/2
               base = base+1
               if(ABS(buf(base)- c(i,j)).gt.1d-8) then
                  print *,me, ' error ',i,j, buf(base), c(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
c
      call ga_sync()
      if (me .eq. 0) then
         write(6,*)'  trans(a)*b: OK '
         call ffflush(6)
      endif
c
      call ga_sync()
      call nga_matmul_patch('n','t', 1d0, 0d0,
     $                      g_a1, alo, ahi,
     $                      g_b1, blo, bhi,
     $                      g_c1, clo, chi)
      call xgemm('n','t',n/2+1,n/2+1,n/3+1,1d0,a(ailo,ajlo), n,
     $            b(bjlo,bilo),n, 0d0, c, n)
      call nga_get(g_c1,clo,chi,buf,bufld)
      base = 0
      do j = 1, 1+n/2
         if(Mod(j,nproc).eq.me) then
            do i = 1, 1+n/2
               base = base+1
               if(ABS(buf(base)- c(i,j)).gt.1d-8) then
                  print *,me, ' error ',i,j, buf(base), c(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
c
      call ga_sync()
      if (me .eq. 0) then
         write(6,*)'  a*trans(b): OK '
         call ffflush(6)
      endif
c
      call nga_matmul_patch('t','t', 1d0, 0d0,
     $                      g_a1, alo, ahi,
     $                      g_b1, blo, bhi,
     $                      g_c1, clo, chi)
      call xgemm('t','t',n/2+1,n/2+1,n/3+1,1d0,a(ajlo,ailo), n,
     $            b(bjlo,bilo),n, 0d0, c, n)
      call nga_get(g_c1,clo,chi,buf,bufld)
      base = 0
      do j = 1, 1+n/2
         if(Mod(j,nproc).eq.me) then
            do i = 1, 1+n/2
               base = base+1
               if(ABS(buf(base)- c(i,j)).gt.1d-8) then
                  print *,me, ' error ',i,j, buf(base), c(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
c
      call ga_sync()
      if (me .eq. 0) then
         write(6,*)'  trans(a)*trans(b): OK '
         call ffflush(6)
      endif
c
      status = ga_destroy(g_a1)
      status = status .and. ga_destroy(g_b1)
      status = status .and. ga_destroy(g_c1)
      if(.not. status) print *, 'ga_destroy failed'
c
      end


      subroutine dpatch_test2()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      integer n,m, dim, dim1
      parameter (n = 10)
      parameter (dim = 3)
      parameter (dim1 = 2)
      double precision alpha, beta, dot
      integer nproc, me
      integer i, j, ailo, ajlo, bilo, bjlo, base, iran
      integer alo(dim), ahi(dim), ald(dim1), chunk(dim)
      integer blo(dim), bhi(dim)
      integer clo(dim), chi(dim)
      integer g_a, g_b, g_c
      logical status
c
      me = ga_nodeid()
      nproc = ga_nnodes()
c
      do i = 1, dim
         chunk(i) = -1
         chi(i) =n
         alo(i) =1
         blo(i) =1
         clo(i) =1
         if (i.ne.dim) ald(i) = n
      enddo
c
      ahi(1)=n
      ahi(2)=1
      status = nga_create(MT_DBL, dim1, ahi, 'a', chunk, g_a)
      if(.not. status) call ga_error('create 1 failed',0)
c
      status = nga_create(MT_DBL, dim1, ahi, 'b', chunk, g_b)
      if(.not. status) call ga_error('create 2 failed',0)
c
      status = nga_create(MT_DBL, dim, chi, 'c', chunk, g_c)
      if(.not. status) call ga_error('create 2 failed',0)
c
      call ga_fill(g_a,1d0)
      call ga_fill(g_b,1d0)
      call ga_zero(g_c)

      if(me.eq.0) print *,'  [1:N,1:N,1:1]=[1,N]x[N,1]^t'
      chi(dim) = 1
      bhi(1)=1
      bhi(2)=n
      call nga_matmul_patch('n','t', 1d0, 0d0,
     $                      g_a, alo, ahi,
     $                      g_b, blo, bhi,
     $                      g_c, clo, chi)
      
c
c     call nga_print_patch(g_c,clo,chi)
      dot = nga_ddot_patch(g_c,'n', clo, chi, g_c,'n', clo, chi)
      if(abs(dot- 1d0*n*n) .gt. 0.01) then
        print *, 'error ', dot, n*n
      else
        if(me.eq.0)print *,'  OK'
      endif
      status = ga_destroy(g_a)
      status = status .and. ga_destroy(g_b)
      status = status .and. ga_destroy(g_c)
      if(.not. status) print *, 'ga_destroy failed'
c
      end


