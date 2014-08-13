#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c $Id: testeig.F,v 1.10 2004-11-13 02:49:18 edo Exp $
      program test
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer heap, stack
c#define BLOCK_CYCLIC 0
c#define PAR_DIAG 1
c     
c
c***  Intitialize a message passing library
c
#include "mp3.fh"
c
c***  initialize PEIGS
#ifdef PAR_DIAG
#ifdef HAVE_SCALAPACK
c     something here?
#else
       call mxinit()   ! PEIGS needs mxinit
#endif
#endif
c 
c     Intitialize the GA package
c
      call ga_initialize()
c
c     Initialize the MA package
c
      heap = 190000
      stack= 190000
      if (.not. ma_init(MT_DBL, heap, stack))
     $    call ga_error("ma init failed",heap+stack)
c
      call testit()
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
      parameter (n = 10)
      double precision a(n,n), b(n,n), c(n,n), evals(n)
      double precision eva(n), evb(n)
      integer g_a,g_b,g_c,g_d
      integer  i, j, index(n),ind(n)
      integer nproc, me
      double precision dsin, sum
      logical status
#if BLOCK_CYCLIC
      integer g_aa, g_bb, g_cc, g_dd
      integer dims(2), proc_grid(2), block(2)
      integer g1, g2
#endif
c     
      nproc = ga_nnodes()
      me    = ga_nodeid()
c     
c***  a() is a local copy of what the global array should start as
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
c
c***  Create global arrays
c
#if BLOCK_CYCLIC
      if (me.eq.0) then
        write(6,*) ' '
        write(6,*) '> Creating Block-Cyclic Arrays'
        write(6,*) ' '
      endif
      dims(1) = n
      dims(2) = n
      block(1) = 16
      block(2) = 16
      call factor(nproc,g1,g2)
      proc_grid(1) = g1
      proc_grid(2) = g2
      if (me.eq.0) then
        write(6,'(a,i3,a,i3)') '  Proc grid dimensions: ',g1,' x ',g2
      endif

      g_a = ga_create_handle()
      call ga_set_data(g_a,2,dims,MT_DBL)
      call ga_set_block_cyclic_proc_grid(g_a,block,proc_grid)
      if (.not.ga_allocate(g_a))
     &     call ga_error(' ga_create a failed ', 2)

      g_b = ga_create_handle()
      call ga_set_data(g_b,2,dims,MT_DBL)
      call ga_set_block_cyclic_proc_grid(g_b,block,proc_grid)
      if (.not.ga_allocate(g_b))
     &     call ga_error(' ga_create b failed ', 2)

      g_c = ga_create_handle()
      call ga_set_data(g_c,2,dims,MT_DBL)
      call ga_set_block_cyclic_proc_grid(g_c,block,proc_grid)
      if (.not.ga_allocate(g_c))
     &     call ga_error(' ga_create c failed ', 2)

      g_d = ga_create_handle()
      call ga_set_data(g_d,2,dims,MT_DBL)
      call ga_set_block_cyclic_proc_grid(g_d,block,proc_grid)
      if (.not.ga_allocate(g_d))
     &     call ga_error(' ga_create d failed ', 2)

      if (.not. ga_create(MT_DBL, n, n, 'AA', 1, 1, g_aa))
     $     call ga_error(' ga_create aa failed ',2)
      if (.not. ga_create(MT_DBL, n, n, 'BB', 1, 1, g_bb))
     $     call ga_error(' ga_create bb failed ',2)
      if (.not. ga_create(MT_DBL, n, n, 'CC', 1, 1, g_cc))
     $     call ga_error(' ga_create cc failed ',2)
      if (.not. ga_create(MT_DBL, n, n, 'DD', 1, 1, g_dd))
     $     call ga_error(' ga_create dd failed ',2)
#else
      if (.not. ga_create(MT_DBL, n, n, 'A', 1, 1, g_a))
     $     call ga_error(' ga_create a failed ',2)
      if (.not. ga_create(MT_DBL, n, n, 'B', 1, 1, g_b))
     $     call ga_error(' ga_create b failed ',2)
      if (.not. ga_create(MT_DBL, n, n, 'C', 1, 1, g_c))
     $     call ga_error(' ga_create c failed ',2)
      if (.not. ga_create(MT_DBL, n, n, 'D', 1, 1, g_d))
     $     call ga_error(' ga_create d failed ',2)
#endif
c     
c***  Fill in arrays A & B
c
      if (me .eq. 0) then
c        write(6,21) 
 21      format(/' filling in ... ')
c        call ffflush(6)
         do j = 1, n
        call ga_put(g_a, 1,n, j,j, a(1,j),n)
        call ga_put(g_b, 1,n, j,j, b(1,j),n)
        call ga_put(g_c, 1,n, j,j, c(1,j),n)
         enddo
c    print *,'A'
        do j = 1, n
c            call GA_GET(g_a, 1,n, j,j, eva,1)
c         write(6,'(10e8.2)')(eva(i),i=1,n)
        enddo
      endif
c
c***  check symmetrization
c
      if (me .eq. 0) then
        print *,' '
        print *,'>checking ga_symmetrize '
        print *,' '
        call ffflush(6)
      endif
      call ga_symmetrize(g_c)
c
      call GA_GET(g_c,  1,n, 1,n,c,n)
      do j = ga_nodeid()+1, n, ga_nnodes()
         do i = j+1, n
            if(c(i,j).ne.c(j,i))then
                 print *, me, ' symmetrize ',i,j,c(i,j),c(j,i)
                 call ffflush(6)
                 call ga_error('exiting',-1)
            endif
         enddo
      enddo
      call ga_sync()
      if (me .eq. 0) then
        print *,' '
        print *,' ga_symmetrize is OK'
        print *,' '
        call ffflush(6)
      endif

c
c***  check symmetrization 
c
      if (me .eq. 0) then
        print *,' '
        print *,'>checking ga_transpose '
        print *,' '
        call ffflush(6)
      endif
c
      call ga_transpose(g_c,g_d)
*     call ga_print(g_c)
*     call ga_print(g_d)
      call GA_GET(g_d,  1,n, 1,n,a,n)
      do j = ga_nodeid()+1, n, ga_nnodes()
         call GA_GET(g_d,  1,n, j,j, a,n)
         do i = 1, n
            if(a(i,1).ne.c(j,i))then
                 print *, me, ' transpose ',i,j,a(i,1),c(j,i) 
                 call ffflush(6)
                 call ga_error('exiting',-1)
            endif
         enddo
      enddo
      call ga_sync()
      if (me .eq. 0) then
        print *,' '
        print *,' ga_transpose is OK'
        print *,' '
        call ffflush(6)
      endif
c
c
c***  solve the eigenproblem
      if (me .eq. 0)then
        print *,' '
        write(6,*) '>checking the generalized eigensolver ... '
        print *,' '
        call ffflush(6)
      endif
      call ga_sync()
#ifndef PAR_DIAG
      call ga_diag_seq(g_a,g_b,g_c,evals)
#else
#ifdef HAVE_SCALAPACK
#if BLOCK_CYCLIC
      call ga_copy(g_a, g_aa)
      call ga_copy(g_b, g_bb)
#endif
      call ga_pdsygv( g_a,  g_b, g_c, evals)
#if BLOCK_CYCLIC
      call ga_copy(g_aa, g_a)
      call ga_copy(g_bb, g_b)
#endif
#else
      call ga_diag(g_a,g_b,g_c,evals)
#endif
#endif
      if (me .eq. 0) then
        print *,' '
        print *,' checking multiplication'
        print *,' '
        call ffflush(6)
      endif
c
      call ga_sync()
#if BLOCK_CYCLIC
      call ga_copy(g_a, g_aa)
      call ga_copy(g_c, g_cc)
      call ga_dgemm('t','n',n,n,n, 1d0, g_cc, g_aa, 0d0, g_dd)
      call ga_dgemm('n','n',n,n,n, 1d0, g_dd, g_cc, 0d0, g_aa)
      call ga_copy(g_aa, g_a)
#else
      call ga_dgemm('t','n',n,n,n, 1d0, g_c, g_a, 0d0, g_d)
      call ga_dgemm('n','n',n,n,n, 1d0, g_d, g_c, 0d0, g_a)
#endif
c
      call ga_sync()
      if (me .eq. 0) then
         do j = 1, n
        call GA_GET(g_a, j,j, j,j, eva(j),1)
         enddo
      endif

      call ga_sync()
#if BLOCK_CYCLIC
      call ga_copy(g_b, g_bb)
      call ga_dgemm('t','n',n,n,n, 1d0, g_cc, g_bb, 0d0, g_dd)
      call ga_sync()
      call ga_dgemm('n','n',n,n,n, 1d0, g_dd, g_cc, 0d0, g_bb)
      call ga_copy(g_bb, g_b)
#else
      call ga_dgemm('t','n',n,n,n, 1d0, g_c, g_b, 0d0, g_d)
      call ga_sync()
      call ga_dgemm('n','n',n,n,n, 1d0, g_d, g_c, 0d0, g_b)
#endif
c
      call ga_sync()
      if (me .eq. 0) then
         do j = 1, n
        call GA_GET(g_b, j,j, j,j, evb(j),1)
         enddo
         write(6,*)'  j   lambda      eva      evb      eva/evb' 
         write(6,*)'  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^' 
         call ffflush(6)
         do j = 1, n
            if(ABS(evals(j) - eva(j)/evb(j)).gt.1d-5)
     &         write(6,'(i4,1h_,6(e10.3,1h,))')
     @            j,evals(j), eva(j), evb(j),eva(j)/evb(j)
         enddo
         write(6,*)'       OK         OK       OK        OK' 
         call ffflush(6)
      endif
      if (me .eq. 0) then
        print *,' '
        print *,'  eigensolver & multiplication are OK'
        print *,' '
        print *,' '
        call ffflush(6)
      endif
c
c..................................................................
c
c***  solve the std eigenproblem
      if (me .eq. 0)then
        print *,' '
        write(6,*) '>checking the standard eigensolver ... '
        print *,' '
        call ffflush(6)
      endif
      do j =1,n
         index(j)=j
         ind(j)=j
      enddo

      call ga_sync()
#ifdef PAR_DIAG
#ifdef HAVE_SCALAPACK
#if BLOCK_CYCLIC
      call ga_copy(g_a, g_aa)
#endif
      call ga_pdsyev( g_a,  g_c, evals, 16)
#if BLOCK_CYCLIC
      call ga_copy(g_aa, g_a)
#endif
#else
      call ga_diag_std(g_a,g_c,evals)
#endif
#else
      call ga_diag_std_seq(g_a,g_c,evals)
#endif
c

      call ga_zero(g_b)
      
#if BLOCK_CYCLIC
      call ga_copy(g_a, g_aa)
      call ga_copy(g_c, g_cc)
      call ga_dgemm('n','n',n,n,n, 1d0, g_aa, g_cc, 0d0, g_dd) ! d := a*c
      call ga_copy(g_dd, g_d)
#else
      call ga_dgemm('n','n',n,n,n, 1d0, g_a, g_c, 0d0, g_d) ! d := a*c
#endif
c
c  copy eigenvalues to diagonal of g_b
c
      if (me .eq. 0) call ga_scatter(g_b, evals, index, ind, n)

      call ga_sync()
#if BLOCK_CYCLIC
      call ga_copy(g_b, g_bb)
      call ga_dgemm('n','n',n,n,n, 1d0, g_cc, g_bb, 0d0, g_aa) ! a := c*b
      call ga_copy(g_aa, g_a)
#else
      call ga_dgemm('n','n',n,n,n, 1d0, g_c, g_b, 0d0, g_a) ! a := c*b
#endif
 
      call ga_sync()
      call ga_add(1d0, g_d, -1d0, g_a, g_c)
      sum = ga_ddot(g_c,g_c)
      if (me .eq. 0) then
        if(dsqrt(sum)/n.lt.1d-11)then
           print *, ' std eigensolver is OK '
        else
           print *, ' test failed: norm = ', dsqrt(sum)/n
        endif
        print *,' '
        call ffflush(6)
      endif
c     status =  MA_summarize_allocated_blocks()
      status =  ga_destroy(g_d)
      status =  ga_destroy(g_c)
      status =  ga_destroy(g_b)
      status =  ga_destroy(g_a)
      end

#if BLOCK_CYCLIC
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
#endif
