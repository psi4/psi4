#if HAVE_CONFIG_H
#   include "config.fh"
#endif
C $iD: PATch.F,v 1.21 2002/10/15 23:35:52 vinod Exp $
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
c
c#define BLOCK_CYCLIC
c#define USE_SCALAPACK_DISTR
#ifdef USE_SCALAPACK_DISTR
#define BLOCK_CYCLIC
#endif
c
      program test
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      logical status
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
      status = ma_init(MT_DBL, 500000, 900000/ga_nnodes())
      if (.not. status)call ga_error( 'ma_init failed', -1)
      if(ga_nodeid().eq.0)then
         print *,' '
         print *,'CHECKING PATCH OPERATIONS FOR DOUBLES '
         print *,' '
      endif
      call dpatch_test()
c
      if(ga_nodeid().eq.0)then
         print *,' '
         print *,'CHECKING PATCH OPERATIONS FOR DOUBLE COMPLEX'
         print *,' '
      endif
      call zpatch_test()
c 
      if(ga_nodeid().eq.0)then
         print *,' '
         print *,'CHECKING PATCH OPERATIONS FOR SINGLE PRECISION'
         print *,' '
      endif
      call spatch_test()
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
      integer n,m
      parameter (n = 128)
      parameter (m = n*n)
      double precision a(n,n), b(n,n), c(n,n), buf(m), val
      double precision alpha, beta
      integer nproc, me 
      integer i, j, ailo, ajlo, bilo, bjlo, base, iran
      integer aihi, ajhi, bihi, bjhi
      integer g_a, g_b, g_c 
      integer rows, cols, loop
      logical status, dist_same 
#ifdef BLOCK_CYCLIC
      integer ndim, dims(2)
      integer block_size(2), proc_grid(2)
#endif
      iran(i) = int(drand(1)*real(i)) + 1
      dist_same = .false.
c
      me = ga_nodeid()
      nproc = ga_nnodes()
#ifdef BLOCK_CYCLIC
      block_size(1) = 32
      block_size(2) = 32
#ifdef USE_SCALAPACK_DISTR
      if (mod(nproc,2).ne.0)
     +  call ga_error("Available procs must be divisible by 2",0)
      proc_grid(1) = 2
      proc_grid(2) = nproc/2
#endif
#endif

c
      do j = 1, n
         do i = 1, n
            a(i,j) = i-1 + (j-1)*n
            b(i,j) = i+j 
         enddo
      enddo
c
c***  Create a global array
c
#ifndef BLOCK_CYCLIC
      status = ga_create(MT_DBL, n, n, 'a', 0, 0, g_a)
#else
      g_a = ga_create_handle()
      ndim = 2
      dims(1) = n
      dims(2) = n
      call ga_set_data(g_a,ndim,dims,MT_DBL)
      call ga_set_array_name(g_a,'a')
#ifdef USE_SCALAPACK_DISTR
      call ga_set_block_cyclic_proc_grid(g_a,block_size,proc_grid)
#else
      call ga_set_block_cyclic(g_a,block_size)
#endif
      status = ga_allocate(g_a)
#endif
      if (.not. status) then
         write(6,*) ' ga_create failed'
         call ffflush(6)
         call ga_error('... exiting ',0)
      endif
c
      if(dist_same) then
         status = ga_duplicate(g_a, g_b, 'a_duplicated')
         if(.not.ga_compare_distr(g_a, g_b))
     $           call ga_error("g_b distribution different",0) 
         status = ga_duplicate(g_a, g_c, 'a_duplicated_again')
         if(.not.ga_compare_distr(g_a, g_c))
     $           call ga_error("g_c distribution different",0) 
      else
#ifndef BLOCK_CYCLIC
         status = ga_create(MT_DBL, n, n, 'b', 0, n, g_b)
#else
         g_b = ga_create_handle()
         ndim = 2
         dims(1) = n
         dims(2) = n
         call ga_set_data(g_b,ndim,dims,MT_DBL)
         call ga_set_array_name(g_b,'b')
#ifdef USE_SCALAPACK_DISTR
         call ga_set_block_cyclic_proc_grid(g_b,block_size,proc_grid)
#else
         call ga_set_block_cyclic(g_b,block_size)
#endif
         status = ga_allocate(g_b)
#endif
         if (.not. status) call ga_error('ga_create failed:b',0) 
#ifndef BLOCK_CYCLIC
         status = ga_create(MT_DBL, n, n, 'c', n, 0, g_c)
#else
         g_c = ga_create_handle()
         ndim = 2
         dims(1) = n
         dims(2) = n
         call ga_set_data(g_c,ndim,dims,MT_DBL)
         call ga_set_array_name(g_c,'c')
#ifdef USE_SCALAPACK_DISTR
         call ga_set_block_cyclic_proc_grid(g_c,block_size,proc_grid)
#else
         call ga_set_block_cyclic(g_c,block_size)
#endif
         status = ga_allocate(g_c)
#endif
         if (.not. status) call ga_error('ga_create failed:c',0) 
      endif
c
c***
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_fill_patch ... '
         call ffflush(6)
      endif
c
      val = 1d0
      call ga_fill_patch(g_a, 2,n/2, 2,n, val)
*     call ga_print(g_a, 1)
      do j = 2+me, n, nproc 
         call ga_get(g_a, 1,n/2, j,j, buf,n/2) 
         do i = 2, n/2
            if(buf(i) .ne.val ) then
               print *,me, ' error ',i,j, buf(i),val 
               call ga_error('exiting ...',0)
            endif
         enddo
      enddo
c
      call ga_sync()
      if (me .eq. 0) then
         write(6,*)'   OK '
         call ffflush(6)
      endif
c  
      do j = 1+me, n, nproc 
         call ga_put(g_a,1,n,j,j,a(1,j),n) 
      enddo
c***
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_copy_patch ... '
         call ffflush(6)
      endif
c
      do loop =1, 10
              ailo = iran(n/2)
              ajlo = iran(n/2)
              aihi = min(n, -1+ailo+n/2)
              ajhi = min(n, -1+ajlo+n/4)
              rows = aihi -ailo+1
              cols = ajhi -ajlo +1
c
              bilo = iran(n/3)
              bjlo = iran(n/3)
              bihi = bilo + rows -1
              bjhi = bjlo + cols -1
              if (me .eq. 0) then
                write(6,'(2x,1h[,4i4,1h],5h --> ,1h[,4i4,1h])')
     $                ailo,aihi,ajlo,ajhi,
     $                bilo,bihi, bjlo, bjhi
 
                call ffflush(6)
              endif

c
              call ga_copy_patch('n', g_a, ailo, aihi, ajlo, ajhi,
     &                                g_b, bilo, bihi, bjlo, bjhi)
c             call ga_print(g_a,1)
c             call ga_print(g_b,1)
              call ga_get(g_b,bilo,bihi, bjlo, bjhi, buf, rows)
              base = 0
              do j = ajlo, ajhi
                 if(Mod(j,nproc).eq.me) then
                    do i = ailo, aihi
                       base = base+1
                       if(buf(base) .ne. a(i,j)) then
                          print *,me, ' error ',i,j, buf(base), a(i,j)
                          call ga_error('exiting ...',0)
                       endif
                    enddo
                 else
                    base = base + rows
                 endif
              enddo
      enddo
c
      ailo = iran(n/2)
      ajlo = iran(n/2)
      bilo = iran(n/2)
      bjlo = iran(n/2)
c
#ifdef BLOCK_CYCLIC
      if (me .eq. 0) then
         write(6,*)'  without transpose:   OK '
         call ffflush(6)
      endif
#else
      if (me .eq. 0) then
                write(6,'(2x,1h[,4i4,1h],5h --> ,1h[,4i4,1h])')
     $                ailo,aihi,ajlo,ajhi,
     $                bilo,bihi, bjlo, bjhi
      endif

      call ga_copy_patch('n', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,1)
*     call ga_print(g_b, 1)
      call ga_get(g_b,bilo,bilo+n/3,bjlo, bjlo+n/2,buf,n/3+1) 
      base = 0
      do j = ajlo, ajlo+n/3 
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2 
               base = base+1
               if(buf(base) .ne. a(i,j)) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
      if (me .eq. 0) then
         write(6,*)'  without transpose:   OK '
         call ffflush(6)
      endif
c 
      if (me .eq. 0) then
                write(6,'(2x,1h[,4i4,2h]~,5h --> ,1h[,4i4,1h])')
     $                ailo,aihi,ajlo,ajhi,
     $                bilo,bihi, bjlo, bjhi
      endif
      call ga_copy_patch('t', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,1)
*     call ga_print(g_b, 1)
      call ga_get(g_b,bilo,bilo+n/3,bjlo, bjlo+n/2,buf,n/3+1)
      base = 0
      do i = ailo, ailo+n/2
         if(Mod(j,nproc).eq.me) then
         do j = ajlo, ajlo+n/3
               base = base+1
               if(buf(base) .ne. a(i,j)) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/3+1
         endif
      enddo
      call ga_sync()
c
      if (me .eq. 0) then
         write(6,*)'  transposed:   OK '
         call ffflush(6)
      endif
#endif
c
c***
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_scale_patch ... '
         call ffflush(6)
      endif
      call ga_copy_patch('n', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/2, bjlo, bjlo+n/3)
      val = 1.d0
      call ga_scale_patch(g_b, bilo,bilo+n/2,bjlo, bjlo+n/3,val)
      call ga_get(g_b,bilo,bilo+n/2,bjlo, bjlo+n/3,buf,n/2+1)
      base = 0
      do j = ajlo, ajlo+n/3
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2
               base = base+1
               if(buf(base) .ne. a(i,j)*val) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)*val
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
      if (me .eq. 0) then
         write(6,*)'  OK '
         call ffflush(6)
      endif
c
c***
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_add_patch ... '
         call ffflush(6)
      endif
      alpha = .1d0
      beta = .2d0
      call ga_zero(g_c)
      call ga_add_patch(alpha, g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                   beta, g_b, bilo,bilo+n/2, bjlo, bjlo+n/3,
     $                         g_c, bilo,bilo+n/2, bjlo, bjlo+n/3)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/3,buf,n/2+1)
      base = 0
      val = val*beta + alpha
      do j = ajlo, ajlo+n/3
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2
               base = base+1
               if(ABS(buf(base)- a(i,j)*val).gt.1d-5) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)*val
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
      if (me .eq. 0) then
         write(6,*)'  OK '
         call ffflush(6)
      endif
c
c***
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_ddot_patch ... '
         call ffflush(6)
      endif
      alpha= ga_ddot_patch(g_a,'n', ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                     g_c,'n', bilo,bilo+n/2, bjlo, bjlo+n/3)
      beta = 0d0
      do j = ajlo, ajlo+n/3
            do i = ailo, ailo+n/2
               beta = beta + a(i,j)**2
            enddo
      enddo

      if(ABS(beta*val- alpha).gt.1d-6*alpha) then
             print *,me, ' error ', beta*val, alpha
             call ga_error('exiting ...',0)
      endif
      call ga_sync()
c
      if (me .eq. 0) then
         write(6,*)'  OK '
         call ffflush(6)
      endif
c
c......................................................
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_matmul_patch ... '
         call ffflush(6)
      endif
      do j = 1+me, n, nproc
         call ga_put(g_b,1,n,j,j,b(1,j),n)
      enddo
      call ga_sync()
      call ga_matmul_patch('n','n', 1d0, 0d0, 
     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
      call xgemm('n','n',n/2+1,n/2+1,n/3+1,1d0,a(ailo,ajlo), n,
     $            b(bilo,bjlo),n, 0d0, c, n)
*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3)
*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print(g_c, 1)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
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
      call ga_matmul_patch('t','n', 1d0, 0d0,
     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
      call xgemm('t','n',n/2+1,n/2+1,n/3+1,1d0,a(ajlo,ailo), n,
     $            b(bilo,bjlo),n, 0d0, c, n)
*     call ga_print(g_a,1) 
*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print(g_c, 1)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
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
      call ga_matmul_patch('n','t', 1d0, 0d0,
     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
      call xgemm('n','t',n/2+1,n/2+1,n/3+1,1d0,a(ailo,ajlo), n,
     $            b(bjlo,bilo),n, 0d0, c, n)
*     call ga_print(g_a,1)
*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print(g_c, 1)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
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
      call ga_sync()
      call ga_matmul_patch('t','t', 1d0, 0d0,
     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
      call xgemm('t','t',n/2+1,n/2+1,n/3+1,1d0,a(ajlo,ailo), n,
     $            b(bjlo,bilo),n, 0d0, c, n)
*     call ga_print(g_a,1)
*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print(g_c, 1)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
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
      status = ga_destroy(g_a)
      status = status .and. ga_destroy(g_b)
      status = status .and. ga_destroy(g_c)
      if(.not. status) print *, 'ga_destroy failed'
c
      end





      subroutine zpatch_test()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c
      integer n,m
      parameter (n = 128)
      parameter (m = n*n)
      double complex a(n,n), b(n,n), c(n,n), buf(m),  val
      double complex alpha, beta
      integer nproc, me 
      integer i, j, ailo, ajlo, bilo, bjlo, base, iran
      integer aihi, ajhi, bihi, bjhi
      integer g_a, g_b, g_c 
      integer rows, cols, loop
      logical status, dist_same 
#ifdef BLOCK_CYCLIC
      integer ndim, dims(2)
      integer block_size(2), proc_grid(2)
#endif
      iran(i) = int(drand(1)*real(i)) + 1
      dist_same = .false.
c
      me = ga_nodeid()
      nproc = ga_nnodes()
#ifdef BLOCK_CYCLIC
      block_size(1) = 32
      block_size(2) = 32
#ifdef USE_SCALAPACK_DISTR
      if (mod(nproc,2).ne.0)
     +  call ga_error("Available procs must be divisible by 2",0)
      proc_grid(1) = 2
      proc_grid(2) = nproc/2
#endif
#endif
c
      do j = 1, n
         do i = 1, n
            a(i,j) = cmplx(dble(i-1), dble((j-1)*n))
            b(i,j) = cmplx(dble(i+j),1d0)
         enddo
      enddo
c
c***  Create a global array
c
#ifndef BLOCK_CYCLIC
      status = ga_create(MT_DCPL, n, n, 'a', 0, 0, g_a)
#else
      g_a = ga_create_handle()
      ndim = 2
      dims(1) = n
      dims(2) = n
      call ga_set_data(g_a,ndim,dims,MT_DCPL)
      call ga_set_array_name(g_a,'a')
#ifdef USE_SCALAPACK_DISTR
      call ga_set_block_cyclic_proc_grid(g_a,block_size,proc_grid)
#else
      call ga_set_block_cyclic(g_a,block_size)
#endif
      status = ga_allocate(g_a)
#endif
      if (.not. status) then
         write(6,*) ' ga_create failed'
         call ffflush(6)
         call ga_error('... exiting ',0)
      endif
c
      if(dist_same) then
         status = ga_duplicate(g_a, g_b, 'a_duplicated')
         if(.not.ga_compare_distr(g_a, g_b))
     $           call ga_error("g_b distribution different",0) 
         status = ga_duplicate(g_a, g_c, 'a_duplicated_again')
         if(.not.ga_compare_distr(g_a, g_c))
     $           call ga_error("g_c distribution different",0) 
      else
#ifndef BLOCK_CYCLIC
         status = ga_create(MT_DCPL, n, n, 'b', 0, n, g_b)
#else
         g_b = ga_create_handle()
         ndim = 2
         dims(1) = n
         dims(2) = n
         call ga_set_data(g_b,ndim,dims,MT_DCPL)
         call ga_set_array_name(g_b,'b')
#ifdef USE_SCALAPACK_DISTR
         call ga_set_block_cyclic_proc_grid(g_b,block_size,proc_grid)
#else
         call ga_set_block_cyclic(g_b,block_size)
#endif
         status = ga_allocate(g_b)
#endif
         if (.not. status) call ga_error('ga_create failed:b',0) 
#ifndef BLOCK_CYCLIC
         status = ga_create(MT_DCPL, n, n, 'c', n, 0, g_c)
#else
         g_c = ga_create_handle()
         ndim = 2
         dims(1) = n
         dims(2) = n
         call ga_set_data(g_c,ndim,dims,MT_DCPL)
         call ga_set_array_name(g_c,'c')
#ifdef USE_SCALAPACK_DISTR
         call ga_set_block_cyclic_proc_grid(g_c,block_size,proc_grid)
#else
         call ga_set_block_cyclic(g_c,block_size)
#endif
         status = ga_allocate(g_c)
#endif
         if (.not. status) call ga_error('ga_create failed:c',0) 
      endif
c
c***
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_fill_patch ... '
         call ffflush(6)
      endif
c
      val = (1d0,-1d0)
      call ga_fill_patch(g_a, 2,n/2, 2,n, val)
*     call ga_print(g_a, 1)
      do j = 2+me, n, nproc 
         call ga_get(g_a, 1,n/2, j,j, buf,n/2) 
         do i = 2, n/2
            if(buf(i) .ne.val ) then
               print *,me, ' error ',i,j, buf(i),val 
               call ga_error('exiting ...',0)
            endif
         enddo
      enddo
c
      call ga_sync()
      if (me .eq. 0) then
         write(6,*)'   OK '
         call ffflush(6)
      endif
c  
      do j = 1+me, n, nproc 
         call ga_put(g_a,1,n,j,j,a(1,j),n) 
      enddo
c***
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_copy_patch ... '
         call ffflush(6)
      endif
c
      do loop =1, 10
              ailo = iran(n/2)
              ajlo = iran(n/2)
              aihi = min(n, -1+ailo+n/2)
              ajhi = min(n, -1+ajlo+n/4)
              rows = aihi -ailo+1
              cols = ajhi -ajlo +1
c
              bilo = iran(n/3)
              bjlo = iran(n/3)
              bihi = bilo + rows -1
              bjhi = bjlo + cols -1
              if (me .eq. 0) then
                write(6,'(2x,1h[,4i4,1h],5h --> ,1h[,4i4,1h])')
     $                ailo,aihi,ajlo,ajhi,
     $                bilo,bihi, bjlo, bjhi
 
                call ffflush(6)
              endif

c
              call ga_copy_patch('n', g_a, ailo, aihi, ajlo, ajhi,
     &                                g_b, bilo, bihi, bjlo, bjhi)
c             call ga_print(g_a,1)
c             call ga_print(g_b,1)
              call ga_get(g_b,bilo,bihi, bjlo, bjhi, buf, rows)
              base = 0
              do j = ajlo, ajhi
                 if(Mod(j,nproc).eq.me) then
                    do i = ailo, aihi
                       base = base+1
                       if(buf(base) .ne. a(i,j)) then
                          print *,me, ' error ',i,j, buf(base), a(i,j)
                          call ga_error('exiting ...',0)
                       endif
                    enddo
                 else
                    base = base + rows
                 endif
              enddo
      enddo
c
      ailo = iran(n/2)
      ajlo = iran(n/2)
      bilo = iran(n/2)
      bjlo = iran(n/2)
c
#ifdef BLOCK_CYCLIC
      if (me .eq. 0) then
         write(6,*)'  without transpose:   OK '
         call ffflush(6)
      endif
#else
      if (me .eq. 0) then
                write(6,'(2x,1h[,4i4,1h],5h --> ,1h[,4i4,1h])')
     $                ailo,aihi,ajlo,ajhi,
     $                bilo,bihi, bjlo, bjhi
      endif

      call ga_copy_patch('n', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,1)
*     call ga_print(g_b, 1)
      call ga_get(g_b,bilo,bilo+n/3,bjlo, bjlo+n/2,buf,n/3+1) 
      base = 0
      do j = ajlo, ajlo+n/3 
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2 
               base = base+1
               if(buf(base) .ne. a(i,j)) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
      if (me .eq. 0) then
         write(6,*)'  without transpose:   OK '
         call ffflush(6)
      endif
c 
      if (me .eq. 0) then
                write(6,'(2x,1h[,4i4,2h]~,5h --> ,1h[,4i4,1h])')
     $                ailo,aihi,ajlo,ajhi,
     $                bilo,bihi, bjlo, bjhi
      endif
      call ga_copy_patch('t', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,1)
*     call ga_print(g_b, 1)
      call ga_get(g_b,bilo,bilo+n/3,bjlo, bjlo+n/2,buf,n/3+1)
      base = 0
      do i = ailo, ailo+n/2
         if(Mod(j,nproc).eq.me) then
         do j = ajlo, ajlo+n/3
               base = base+1
               if(buf(base) .ne. a(i,j)) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/3+1
         endif
      enddo
      call ga_sync()
c
      if (me .eq. 0) then
         write(6,*)'  transposed:   OK '
         call ffflush(6)
      endif
#endif
c
c***
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_scale_patch ... '
         call ffflush(6)
      endif
      call ga_copy_patch('n', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/2, bjlo, bjlo+n/3)
      val = 1.d0
      call ga_scale_patch(g_b, bilo,bilo+n/2,bjlo, bjlo+n/3,val)
      call ga_get(g_b,bilo,bilo+n/2,bjlo, bjlo+n/3,buf,n/2+1)
      base = 0
      do j = ajlo, ajlo+n/3
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2
               base = base+1
               if(buf(base) .ne. a(i,j)*val) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)*val
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
      if (me .eq. 0) then
         write(6,*)'  OK '
         call ffflush(6)
      endif
c
c***
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_add_patch ... '
         call ffflush(6)
      endif
      alpha = .1d0
      beta = .2d0
      call ga_zero(g_c)
      call ga_add_patch(alpha, g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                   beta, g_b, bilo,bilo+n/2, bjlo, bjlo+n/3,
     $                         g_c, bilo,bilo+n/2, bjlo, bjlo+n/3)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/3,buf,n/2+1)
      base = 0
      val = val*beta + alpha
      do j = ajlo, ajlo+n/3
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2
               base = base+1
               if(ABS(buf(base)- a(i,j)*val).gt.1d-5) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)*val
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
      if (me .eq. 0) then
         write(6,*)'  OK '
         call ffflush(6)
      endif
c
c***
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_zdot_patch ... '
         call ffflush(6)
      endif
      alpha= ga_zdot_patch(g_a,'n', ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                     g_c,'n', bilo,bilo+n/2, bjlo, bjlo+n/3)
      beta = (0d0,0d0)
      do j = ajlo, ajlo+n/3
            do i = ailo, ailo+n/2
               beta = beta + a(i,j)*a(i,j)
            enddo
      enddo
      if(ABS(beta*val- alpha)/abs(alpha).gt.1d-6) then
             print *,me, ' error ', beta*val, alpha
             call ga_error('exiting ...',0)
      endif
      call ga_sync()
c
      if (me .eq. 0) then
         write(6,*)'  OK '
         call ffflush(6)
      endif
c
c......................................................
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_matmul_patch ... '
         call ffflush(6)
      endif
      do j = 1+me, n, nproc
         call ga_put(g_b,1,n,j,j,b(1,j),n)
      enddo
      call ga_sync()
      call ga_matmul_patch('n','n', (1d0,0d0), (0d0,0d0),
     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
      call ygemm('n','n',n/2+1,n/2+1,n/3+1,(1d0,0d0),
     $            a(ailo,ajlo), n,
     $            b(bilo,bjlo),n, (0d0,0d0), c, n)
*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3)
*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print(g_c, 1)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
      base = 0
      do j = 1, 1+n/2
         if(Mod(j,nproc).eq.me) then
            do i = 1, 1+n/2
               base = base+1
               if(ABS(buf(base)- c(i,j))/abs(c(i,j)).gt.1d-8) then
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
      call ga_matmul_patch('t','n', (1d0,0d0), (0d0,0d0),
     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
      call ygemm('t','n',n/2+1,n/2+1,n/3+1,(1d0,0d0),
     $            a(ajlo,ailo), n,
     $            b(bilo,bjlo),n, (0d0,0d0), c, n)
*     call ga_print(g_a,1) 
*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print(g_c, 1)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
      base = 0
      do j = 1, 1+n/2
         if(Mod(j,nproc).eq.me) then
            do i = 1, 1+n/2
               base = base+1
               if(ABS(buf(base)- c(i,j))/abs(c(i,j)).gt.1d-8) then
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
      call ga_matmul_patch('n','t', (1d0,0d0), (0d0,0d0),
     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
      call ygemm('n','t',n/2+1,n/2+1,n/3+1,(1d0,0d0),
     $            a(ailo,ajlo), n,
     $            b(bjlo,bilo),n, (0d0,0d0), c, n)
*     call ga_print(g_a,1)
*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print(g_c, 1)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
      base = 0
      do j = 1, 1+n/2
         if(Mod(j,nproc).eq.me) then
            do i = 1, 1+n/2
               base = base+1
               if(ABS(buf(base)- c(i,j))/abs(c(i,j)).gt.1d-8) then
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
      call ga_sync()
      call ga_matmul_patch('t','t', (1d0,0d0), (0d0,0d0),
     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
      call ygemm('t','t',n/2+1,n/2+1,n/3+1,(1d0,0d0),
     $            a(ajlo,ailo), n,
     $            b(bilo,bjlo),n, (0d0,0d0), c, n)
*     call ga_print(g_a,1)
*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print(g_c, 1)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
      base = 0
      do j = 1, 1+n/2
         if(Mod(j,nproc).eq.me) then
            do i = 1, 1+n/2
               base = base+1
               if(ABS(buf(base)- c(i,j))/abs(c(i,j)).gt.1d-8) then
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
      status = ga_destroy(g_a)
      status = status .and. ga_destroy(g_b)
      status = status .and. ga_destroy(g_c)
      if(.not. status) print *, 'ga_destroy failed'
c
      end



      subroutine spatch_test()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c
      integer n,m
      parameter (n = 128)
      parameter (m = n*n)
      real a(n,n), b(n,n), c(n,n), buf(m), val
      real alpha, beta
      integer nproc, me
      integer i, j, ailo, ajlo, bilo, bjlo, base, iran
      integer aihi, ajhi, bihi, bjhi
      integer g_a, g_b, g_c
      integer rows, cols, loop
      logical status, dist_same
#ifdef BLOCK_CYCLIC
      integer ndim, dims(2)
      integer block_size(2), proc_grid(2)
#endif
      iran(i) = int(drand(1)*real(i)) + 1
      dist_same = .false.
c
      me = ga_nodeid()
      nproc = ga_nnodes()
#ifdef BLOCK_CYCLIC
      block_size(1) = 32
      block_size(2) = 32
#ifdef USE_SCALAPACK_DISTR
      if (mod(nproc,2).ne.0)
     +  call ga_error("Available procs must be divisible by 2",0)
      proc_grid(1) = 2
      proc_grid(2) = nproc/2
#endif
#endif
c
      do j = 1, n
         do i = 1, n
            a(i,j) = i-1 + (j-1)*n
            b(i,j) = i+j
         enddo
      enddo
c
c***  Create a global array
c
#ifndef BLOCK_CYCLIC
      status = ga_create(MT_REAL, n, n, 'a', 0, 0, g_a)
#else
      g_a = ga_create_handle()
      ndim = 2
      dims(1) = n
      dims(2) = n
      call ga_set_data(g_a,ndim,dims,MT_REAL)
      call ga_set_array_name(g_a,'a')
#ifdef USE_SCALAPACK_DISTR
      call ga_set_block_cyclic_proc_grid(g_a,block_size,proc_grid)
#else
      call ga_set_block_cyclic(g_a,block_size)
#endif
      status = ga_allocate(g_a)
#endif
      if (.not. status) then
         write(6,*) ' ga_create failed'
         call ffflush(6)
         call ga_error('... exiting ',0)
      endif
c
      if(dist_same) then
         status = ga_duplicate(g_a, g_b, 'a_duplicated')
         if(.not.ga_compare_distr(g_a, g_b))
     $           call ga_error("g_b distribution different",0)
         status = ga_duplicate(g_a, g_c, 'a_duplicated_again')
         if(.not.ga_compare_distr(g_a, g_c))
     $           call ga_error("g_c distribution different",0)
      else
#ifndef BLOCK_CYCLIC
         status = ga_create(MT_REAL, n, n, 'b', 0, n, g_b)
#else
         g_b = ga_create_handle()
         ndim = 2
         dims(1) = n
         dims(2) = n
         call ga_set_data(g_b,ndim,dims,MT_REAL)
         call ga_set_array_name(g_b,'b')
#ifdef USE_SCALAPACK_DISTR
         call ga_set_block_cyclic_proc_grid(g_b,block_size,proc_grid)
#else
         call ga_set_block_cyclic(g_b,block_size)
#endif
         status = ga_allocate(g_b)
#endif
         if (.not. status) call ga_error('ga_create failed:b',0)
#ifndef BLOCK_CYCLIC
         status = ga_create(MT_REAL, n, n, 'c', n, 0, g_c)
#else
         g_c = ga_create_handle()
         ndim = 2
         dims(1) = n
         dims(2) = n
         call ga_set_data(g_c,ndim,dims,MT_REAL)
         call ga_set_array_name(g_c,'c')
#ifdef USE_SCALAPACK_DISTR
         call ga_set_block_cyclic_proc_grid(g_c,block_size,proc_grid)
#else
         call ga_set_block_cyclic(g_c,block_size)
#endif
         status = ga_allocate(g_c)
#endif
         if (.not. status) call ga_error('ga_create failed:c',0)
      endif
c
c***
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_fill_patch ... '
         call ffflush(6)
      endif
c
      val = 1.0
      call ga_fill_patch(g_a, 2,n/2, 2,n, val)
*     call ga_print(g_a, 1)
      do j = 2+me, n, nproc
         call ga_get(g_a, 1,n/2, j,j, buf,n/2)
         do i = 2, n/2
            if(buf(i) .ne.val ) then
               print *,me, ' error ',i,j, buf(i),val
               call ga_error('exiting ...',0)
            endif
         enddo
      enddo
c
      call ga_sync()
      if (me .eq. 0) then
         write(6,*)'   OK '
         call ffflush(6)
      endif
c 
      do j = 1+me, n, nproc
         call ga_put(g_a,1,n,j,j,a(1,j),n)
      enddo
c***
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_copy_patch ... '
         call ffflush(6)
      endif
c
      do loop =1, 10
              ailo = iran(n/2)
              ajlo = iran(n/2)
              aihi = min(n, -1+ailo+n/2)
              ajhi = min(n, -1+ajlo+n/4)
              rows = aihi -ailo+1
              cols = ajhi -ajlo +1
c
              bilo = iran(n/3)
              bjlo = iran(n/3)
              bihi = bilo + rows -1
              bjhi = bjlo + cols -1
              if (me .eq. 0) then
                write(6,'(2x,1h[,4i4,1h],5h --> ,1h[,4i4,1h])')
     $                ailo,aihi,ajlo,ajhi,
     $                bilo,bihi, bjlo, bjhi

                call ffflush(6)
              endif

c
              call ga_copy_patch('n', g_a, ailo, aihi, ajlo, ajhi,
     &                                g_b, bilo, bihi, bjlo, bjhi)
c             call ga_print(g_a,1)
c             call ga_print(g_b,1)
              call ga_get(g_b,bilo,bihi, bjlo, bjhi, buf, rows)
              base = 0
              do j = ajlo, ajhi
                 if(Mod(j,nproc).eq.me) then
                    do i = ailo, aihi
                       base = base+1
                       if(buf(base) .ne. a(i,j)) then
                          print *,me, ' error ',i,j, buf(base), a(i,j)
                          call ga_error('exiting ...',0)
                       endif
                    enddo
                 else
                    base = base + rows
                 endif
              enddo
      enddo
c
      ailo = iran(n/2)
      ajlo = iran(n/2)
      bilo = iran(n/2)
      bjlo = iran(n/2)
c
#ifdef BLOCK_CYCLIC
      if (me .eq. 0) then
         write(6,*)'  without transpose:   OK '
         call ffflush(6)
      endif
#else
      if (me .eq. 0) then
                write(6,'(2x,1h[,4i4,1h],5h --> ,1h[,4i4,1h])')
     $                ailo,aihi,ajlo,ajhi,
     $                bilo,bihi, bjlo, bjhi
      endif

      call ga_copy_patch('n', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,1)
*     call ga_print(g_b, 1)
      call ga_get(g_b,bilo,bilo+n/3,bjlo, bjlo+n/2,buf,n/3+1)
      base = 0
      do j = ajlo, ajlo+n/3
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2
               base = base+1
               if(buf(base) .ne. a(i,j)) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
      if (me .eq. 0) then
         write(6,*)'  without transpose:   OK '
         call ffflush(6)
      endif
c
      if (me .eq. 0) then
                write(6,'(2x,1h[,4i4,2h]~,5h --> ,1h[,4i4,1h])')
     $                ailo,aihi,ajlo,ajhi,
     $                bilo,bihi, bjlo, bjhi
      endif
      call ga_copy_patch('t', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,1)
*     call ga_print(g_b, 1)
      call ga_get(g_b,bilo,bilo+n/3,bjlo, bjlo+n/2,buf,n/3+1)
      base = 0
      do i = ailo, ailo+n/2
         if(Mod(j,nproc).eq.me) then
         do j = ajlo, ajlo+n/3
               base = base+1
               if(buf(base) .ne. a(i,j)) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/3+1
         endif
      enddo
      call ga_sync()
c
      if (me .eq. 0) then
         write(6,*)'  transposed:   OK '
         call ffflush(6)
      endif
#endif
c
c***
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_scale_patch ... '
         call ffflush(6)
      endif

      call ga_copy_patch('n', g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                        g_b, bilo,bilo+n/2, bjlo, bjlo+n/3)
      val = 1.0
      call ga_scale_patch(g_b, bilo,bilo+n/2,bjlo, bjlo+n/3,val)
      call ga_get(g_b,bilo,bilo+n/2,bjlo, bjlo+n/3,buf,n/2+1)
      base = 0
      do j = ajlo, ajlo+n/3
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2
               base = base+1
               if(buf(base) .ne. a(i,j)*val) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)*val
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
      if (me .eq. 0) then
         write(6,*)'  OK '
         call ffflush(6)
      endif
c
c***
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_add_patch ... '
         call ffflush(6)
      endif
      alpha = 0.1
      beta = 0.2
      call ga_zero(g_c)

      call ga_add_patch(alpha, g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                   beta, g_b, bilo,bilo+n/2, bjlo, bjlo+n/3,
     $                         g_c, bilo,bilo+n/2, bjlo, bjlo+n/3)
      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/3,buf,n/2+1)
      base = 0
      val = val*beta + alpha
      do j = ajlo, ajlo+n/3
         if(Mod(j,nproc).eq.me) then
            do i = ailo, ailo+n/2
               base = base+1
               if(ABS(buf(base)- a(i,j)*val).gt.1e-3) then
                  print *,me, ' error ',i,j, buf(base), a(i,j)*val
                  call ga_error('exiting ...',0)
               endif
            enddo
         else
            base = base +n/2+1
         endif
      enddo
      call ga_sync()
c
      if (me .eq. 0) then
         write(6,*)'  OK '
         call ffflush(6)
      endif
c
c***
      if (me .eq. 0) then
         print *, ' '
         write(6,*)'> Checking ga_sdot_patch ... '
         call ffflush(6)
      endif
      alpha= ga_sdot_patch(g_a,'n', ailo,ailo+n/2, ajlo, ajlo+n/3,
     $                     g_c,'n', bilo,bilo+n/2, bjlo, bjlo+n/3)
      beta = 0.0

      do j = ajlo, ajlo+n/3
            do i = ailo, ailo+n/2
               beta = beta + a(i,j)**2
            enddo
      enddo
      if(ABS(beta*val- alpha).gt.(1.0e-2*alpha)) then
             print *,me, ' error ', beta*val, alpha
             call ga_error('exiting ...',0)
      endif
      call ga_sync()
c
      if (me .eq. 0) then
         write(6,*)'  OK '
         call ffflush(6)
      endif
c......................................................
ccc      if (me .eq. 0) then
ccc         print *, ' '
ccc         write(6,*)'> Checking ga_matmul_patch ... '
ccc         call ffflush(6)
ccc      endif
ccc      do j = 1+me, n, nproc
ccc         call ga_put(g_b,1,n,j,j,b(1,j),n)
ccc      enddo
ccc      call ga_sync()
ccc      call ga_matmul_patch('n','n', 1.0, 0.0,
ccc     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
ccc     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
ccc     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
ccc      call xgemm('n','n',n/2+1,n/2+1,n/3+1,1d0,a(ailo,ajlo), n,
ccc     $            b(bilo,bjlo),n, 0d0, c, n)
ccc*     call ga_print_patch(g_a, ailo,ailo+n/2, ajlo, ajlo+n/3)
ccc*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
ccc*     call ga_print(g_c, 1)
ccc      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
ccc      base = 0
ccc      do j = 1, 1+n/2
ccc         if(Mod(j,nproc).eq.me) then
ccc            do i = 1, 1+n/2
ccc               base = base+1
ccc               if(ABS(buf(base)- c(i,j)).gt.1e-8) then
ccc                  print *,me, ' error ',i,j, buf(base), c(i,j)
ccc                  call ga_error('exiting ...',0)
ccc               endif
ccc            enddo
ccc         else
ccc            base = base +n/2+1
ccc         endif
ccc      enddo
cccc
ccc      call ga_sync()
ccc      if (me .eq. 0) then
ccc         write(6,*)'  a*b: OK '
ccc         call ffflush(6)
ccc      endif
cccc
ccc      call ga_sync()
ccc      call ga_matmul_patch('t','n', 1d0, 0d0,
ccc     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
ccc     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
ccc     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
ccc      call xgemm('t','n',n/2+1,n/2+1,n/3+1,1.0,a(ajlo,ailo), n,
ccc     $            b(bilo,bjlo),n, 0.0, c, n)
ccc*     call ga_print(g_a,1)
ccc*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
ccc*     call ga_print(g_c, 1)
ccc      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
ccc      base = 0
ccc      do j = 1, 1+n/2
ccc         if(Mod(j,nproc).eq.me) then
ccc            do i = 1, 1+n/2
ccc               base = base+1
ccc               if(ABS(buf(base)- c(i,j)).gt.1e-8) then
ccc                  print *,me, ' error ',i,j, buf(base), c(i,j)
ccc                  call ga_error('exiting ...',0)
ccc               endif
ccc            enddo
ccc         else
ccc            base = base +n/2+1
ccc         endif
ccc      enddo
cccc
ccc      call ga_sync()
ccc      if (me .eq. 0) then
ccc         write(6,*)'  trans(a)*b: OK '
ccc         call ffflush(6)
ccc      endif
cccc
ccc      call ga_sync()
ccc      call ga_matmul_patch('n','t', 1.0, 0.0,
ccc     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
ccc     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
ccc     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
ccc      call xgemm('n','t',n/2+1,n/2+1,n/3+1,1.0,a(ailo,ajlo), n,
ccc     $            b(bjlo,bilo),n, 0.0, c, n)
ccc*     call ga_print(g_a,1)
ccc*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
ccc*     call ga_print(g_c, 1)
ccc      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
ccc      base = 0
ccc      do j = 1, 1+n/2
ccc         if(Mod(j,nproc).eq.me) then
ccc            do i = 1, 1+n/2
ccc               base = base+1
ccc               if(ABS(buf(base)- c(i,j)).gt.1e-8) then
ccc                  print *,me, ' error ',i,j, buf(base), c(i,j)
ccc                  call ga_error('exiting ...',0)
ccc               endif
ccc            enddo
ccc         else
ccc            base = base +n/2+1
ccc         endif
ccc      enddo
cccc
cccc
ccc      call ga_sync()
ccc      if (me .eq. 0) then
ccc         write(6,*)'  a*trans(b): OK '
ccc         call ffflush(6)
ccc      endif
cccc
ccc      call ga_sync()
ccc      call ga_matmul_patch('t','t', 1d0, 0d0,
ccc     $                      g_a, ailo,ailo+n/2, ajlo, ajlo+n/3,
ccc     $                      g_b, bilo,bilo+n/3, bjlo, bjlo+n/2,
ccc     $                      g_c, bilo,bilo+n/2, bjlo, bjlo+n/2)
ccc      call xgemm('t','t',n/2+1,n/2+1,n/3+1,1.0,a(ajlo,ailo), n,
ccc     $            b(bjlo,bilo),n, 0.0, c, n)
ccc*     call ga_print(g_a,1)
ccc*     call ga_print_patch(g_b, bilo,bilo+n/3, bjlo, bjlo+n/2)
ccc*     call ga_print(g_c, 1)
ccc      call ga_get(g_c,bilo,bilo+n/2,bjlo, bjlo+n/2,buf,n/2+1)
ccc      base = 0
ccc      do j = 1, 1+n/2
ccc         if(Mod(j,nproc).eq.me) then
ccc            do i = 1, 1+n/2
ccc               base = base+1
ccc               if(ABS(buf(base)- c(i,j)).gt.1e-8) then
ccc                  print *,me, ' error ',i,j, buf(base), c(i,j)
ccc                  call ga_error('exiting ...',0)
ccc               endif
ccc            enddo
ccc         else
ccc            base = base +n/2+1
ccc         endif
ccc      enddo
cccc
ccc      call ga_sync()
ccc      if (me .eq. 0) then
ccc         write(6,*)'  trans(a)*trans(b): OK '
ccc         call ffflush(6)
ccc      endif
ccc      status = ga_destroy(g_a)
ccc      status = status .and. ga_destroy(g_b)
ccc      status = status .and. ga_destroy(g_c)
ccc      if(.not. status) print *, 'ga_destroy failed'
cccc
      end


