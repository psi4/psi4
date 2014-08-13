#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c $Id: testsolve.F,v 1.11 2006/03/20 20:02:29 manoj Exp $
      program test
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer heap, stack
#define BLOCK_CYCLIC 0
c
c***  Intitialize a message passing library
c
#include "mp3.fh"
c
c     Intitialize the GA package
c
      call ga_initialize()
c     if(ga_nodeid().eq.0)print *,ga_nnodes(),' nodes'
c
c     Initialize the MA package
c
      heap = 190000
      stack= 190000
      if (.not. ma_init(MT_DBL, heap, stack))
     $    call ga_error("ma init failed",heap+stack) 
c
c
      call testit()
c
      if(ga_nodeid().eq.0) print *,'All tests successful '
c
      call ga_terminate()
c
      call MP_FINALIZE()
      end


c-----------------

      subroutine testit()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
c     
      integer n
      parameter (n = 100)
      double precision a(n,n), b(n,n), c(n,n)
      integer g_a,g_b,g_c,g_d, g_e, g_f, g_g
      integer g_aa,g_bb,g_cc,g_dd, g_ee, g_ff, g_gg
      integer  i, j
      integer ndim, dims(2), block(2), proc_grid(2), g1, g2
      integer nproc, me
      double precision dsin, sum
      logical status
c     
      nproc = ga_nnodes()
      me    = ga_nodeid()
c     
c     a() is a local copy of what the global array should start as
c
      do j = 1, n
         do i = 1, n
            a(i,j) = 1d0 * (i+j)  
            b(i,j) = DSIN(1d0* (i+j))   
        if(i.eq.j) then
               b(i,j) = 2d0 *n
               a(i,j) = i
            endif
        if(i.le.j)then
               c(i,j) = a(i,j)
            else
               c(i,j) = 0d0
            endif
         enddo
      enddo
#if 0
      if (me.eq.0) then
        open(unit=7,file='amat.dat',status='unknown')
        do  i = 1, min(20,n)
          write(7,128) (a(i,j),j=1,min(20,n))
  128     format(20f6.1)
        end do
      endif
#endif
c
c***  Create global arrays
#if BLOCK_CYCLIC
      if (me.eq.0) then
        write(6,*) '*'
        write(6,*) '* Creating Block-Cyclic Arrays'
        write(6,*) '*'
      endif
      dims(1) = n
      dims(2) = n
      block(1) = 2
      block(2) = 2
      call factor(nproc,g1,g2)
      proc_grid(1) = g1
      proc_grid(2) = g2
      g_a = ga_create_handle()
      call ga_set_data(g_a,2,dims,MT_DBL)
      call ga_set_block_cyclic_proc_grid(g_a,block,proc_grid)
      status = ga_allocate(g_a)
      g_b = ga_create_handle()
      call ga_set_data(g_b,2,dims,MT_DBL)
      call ga_set_block_cyclic_proc_grid(g_b,block,proc_grid)
      status = ga_allocate(g_b)
      g_c = ga_create_handle()
      call ga_set_data(g_c,2,dims,MT_DBL)
      call ga_set_block_cyclic_proc_grid(g_c,block,proc_grid)
      status = ga_allocate(g_c)
      g_d = ga_create_handle()
      call ga_set_data(g_d,2,dims,MT_DBL)
      call ga_set_block_cyclic_proc_grid(g_d,block,proc_grid)
      status = ga_allocate(g_d)
      dims(2) = 1
      g_e = ga_create_handle()
      call ga_set_data(g_e,2,dims,MT_DBL)
      call ga_set_block_cyclic_proc_grid(g_e,block,proc_grid)
      status = ga_allocate(g_e)
      g_f = ga_create_handle()
      call ga_set_data(g_f,2,dims,MT_DBL)
      call ga_set_block_cyclic_proc_grid(g_f,block,proc_grid)
      status = ga_allocate(g_f)
      g_g = ga_create_handle()
      call ga_set_data(g_g,2,dims,MT_DBL)
      call ga_set_block_cyclic_proc_grid(g_g,block,proc_grid)
      status = ga_allocate(g_g)
      g_aa = ga_create_handle()
      dims(2) = n
      call ga_set_data(g_aa,2,dims,MT_DBL)
      status = ga_allocate(g_aa)
      g_bb = ga_create_handle()
      call ga_set_data(g_bb,2,dims,MT_DBL)
      status = ga_allocate(g_bb)
      g_cc = ga_create_handle()
      call ga_set_data(g_cc,2,dims,MT_DBL)
      status = ga_allocate(g_cc)
      g_dd = ga_create_handle()
      call ga_set_data(g_dd,2,dims,MT_DBL)
      status = ga_allocate(g_dd)
      dims(2) = 1
      g_ee = ga_create_handle()
      call ga_set_data(g_ee,2,dims,MT_DBL)
      status = ga_allocate(g_ee)
      g_ff = ga_create_handle()
      call ga_set_data(g_ff,2,dims,MT_DBL)
      status = ga_allocate(g_ff)
      g_gg = ga_create_handle()
      call ga_set_data(g_gg,2,dims,MT_DBL)
      status = ga_allocate(g_gg)
#else
      if (.not. ga_create(MT_DBL, n, n, 'a', 1, 1, g_a))
     $     call ga_error(' ga_create failed ',2)
      if (.not. ga_create(MT_DBL, n, n, 'b', 1, 1, g_b))
     $     call ga_error(' ga_create failed ',2)
      if (.not. ga_create(MT_DBL, n, n, 'c', 1, 1, g_c))
     $     call ga_error(' ga_create failed ',2)
      if (.not. ga_create(MT_DBL, n, n, 'd', 1, 1, g_d))
     $     call ga_error(' ga_create failed ',2)
      if (.not. ga_create(MT_DBL, n, 1, 'e', 1, 1, g_e))
     $     call ga_error(' ga_create failed ',2)
      if (.not. ga_create(MT_DBL, n, 1, 'f', 1, 1, g_f))
     $     call ga_error(' ga_create failed ',2)
      if (.not. ga_create(MT_DBL, n, 1, 'g', 1, 1, g_g))
     $     call ga_error(' ga_create failed ',2)
#endif
c     
c     
c***  Fill in arrays A & B
      if (me .eq. 0) then
         print *,  ' filling in A and B  '
         call ffflush(6)
         call ga_put(g_e, 1,n, 1,1, b(1,1),n)
         call ga_put(g_f, 1,n, 1,1, b(1,1),n)
      endif
      do j = 1+me, n, nproc 
        call ga_put(g_a, 1,n, j,j, a(1,j),n)
        call ga_put(g_b, 1,n, j,j, b(1,j),n)
        call ga_put(g_c, 1,n, j,j, b(1,j),n)
      enddo
      call ga_sync()
c
c     call ga_copy(g_b,g_c)
c
      if (me .eq. 0) then
        print *,' '
        print *, '>Test of the LU-based solver with nxn rhs '
        print *,' '
        call ffflush(6)
      endif
#ifndef HAVE_SCALAPACK
      call ga_lu_solve_seq('n', g_a, g_b)
#else
#if BLOCK_CYCLIC
      call ga_copy(g_a,g_aa)
#endif
      call ga_lu_solve('n', g_a, g_b)
#if BLOCK_CYCLIC
      call ga_copy(g_aa,g_a)
#endif
#endif
#if BLOCK_CYCLIC
c      call print_block(g_b)
      call ga_copy(g_b,g_bb)
      call ga_copy(g_c,g_cc)
      call ga_copy(g_d,g_dd)
      call ga_dgemm('n','n',n,n,n, 1d0, g_aa, g_bb, 0d0, g_dd) ! d := a*b
      call ga_add(1d0, g_dd, -1d0, g_cc, g_cc) 
      sum = ga_ddot(g_cc,g_cc)
#else
c      call print_rblock(g_b)
      call ga_dgemm('n','n',n,n,n, 1d0, g_a, g_b, 0d0, g_d) ! d := a*b
      call ga_add(1d0, g_d, -1d0, g_c, g_c) 
      sum = ga_ddot(g_c,g_c)
#endif
      if (me .eq. 0) then
        print *,' '
        print *, ' dsqrt(sum) = ', dsqrt(sum)
        print *, ' n = ', n
        print *, ' norm = ', dsqrt(sum)/n
        if(dsqrt(sum)/n.lt.1d-10) then
           print *, ' test passed '
        else
           call ga_error(' test failed ',3)
        endif
        print *,' '
        call ffflush(6)
      endif
c
      if (me .eq. 0) then
        print *,' '
        print *,'>Test of the LU-based solver with a single vector rhs'
        print *,' '
        call ffflush(6)
      endif
c
#ifndef HAVE_SCALAPACK
      call ga_lu_solve_seq('n', g_a, g_e)
#else
#if BLOCK_CYCLIC
      call ga_copy(g_a,g_aa)
#endif
      call ga_lu_solve('n', g_a, g_e)
#endif
c
#if BLOCK_CYCLIC
      call ga_copy(g_e,g_ee)
      call ga_copy(g_f,g_ff)
      call ga_copy(g_g,g_gg)
      call ga_dgemm('n','n',n,1,n, 1d0, g_aa, g_ee, 0d0, g_gg) ! g := a*e
      call ga_add(1d0, g_gg, -1d0, g_ff, g_ff) 
      sum = ga_ddot(g_ff,g_ff)
#else
      call ga_dgemm('n','n',n,1,n, 1d0, g_a, g_e, 0d0, g_g) ! g := a*e
      call ga_add(1d0, g_g, -1d0, g_f, g_f) 
      sum = ga_ddot(g_f,g_f)
#endif
      if (me .eq. 0) then
        print *,' '
        print *, ' norm = ', dsqrt(sum)/n
        if(dsqrt(sum)/n.lt.1d-10) then
           print *, ' test passed '
        else
           call ga_error(' test failed ',4)
        endif           
        print *,' '
        call ffflush(6)
      endif
      end
c
      subroutine factor(p,idx,idy)
      implicit none
      integer i,j,p,idx,idy,it
      integer ip,ifac,pmax,prime(1280)
      integer fac(1280)
c
      i = 1
      ip = p
c
c    factor p completely
c    first, find all prime numbers less than or equal to p
c
      pmax = 0
      do i = 2, p
        do j = 1, pmax
          if (mod(i,prime(j)).eq.0) go to 100
        end do
        pmax = pmax + 1
        prime(pmax) = i
  100   continue
      end do
c
c    find all prime factors of p
c
      ifac = 0
      do i = 1, pmax
  200   if (mod(ip,prime(i)).eq.0) then
          ifac = ifac + 1
          fac(ifac) = prime(i)
          ip = ip/prime(i)
          go to 200
        endif
      end do
c
c    determine two factors of p of approximately the
c    same size
c
      idx = 1
      idy = 1
      do i = ifac, 1, -1
        if (idx.le.idy) then
          idx = fac(i)*idx
        else
          idy = fac(i)*idy
        endif
      end do
      return
      end
c
      subroutine print_block(g_a)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer g_a, i, j, istride, jstride
      integer idx, ld, me
      me = ga_nodeid()
      istride = 10
      if (me.eq.0) then
      jstride = 6
      else
      jstride = 4
      endif
      call nga_access_block_segment(g_a,me,idx,ld)
      do i = 1, istride
        write(6,733) me,(dbl_mb(idx+(j-1)*istride+i-1),j=1,jstride)
      end do
  733 format(i8,8f12.6)
      return
      end
c
      subroutine print_rblock(g_a)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer g_a, i, j, istride, jstride
      integer idx, ld, me, lo(2), hi(2)
      me = ga_nodeid()
      istride = 10
      jstride = 5
      call nga_distribution(g_a,me,lo,hi)
      call nga_access(g_a,lo,hi,idx,ld)
      do i = 1, istride
        write(6,733) me,(dbl_mb(idx+(j-1)*istride+i-1),j=1,jstride)
      end do
  733 format(i8,8f12.6)
      return
      end
