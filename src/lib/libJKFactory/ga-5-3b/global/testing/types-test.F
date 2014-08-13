#if HAVE_CONFIG_H 
#   include "config.fh"
#endif
c $Id: test.F,v 1.64.2.11 2007-04-06 22:37:35 d3g293 Exp $
c vector boxes lack arithmetic precision 
#if defined(FUJITSU)
# define THRESH 1d-12
# define THRESHF 1e-5
#else
# define THRESH 1d-13
# define THRESHF 2e-5
#endif

#define MISMATCH(x,y) abs(x-y)/max(1d0,abs(x)).gt.THRESH
#define MISMATCHF(x,y) abs(x-y)/max(1.0,abs(x)).gt.THRESHF 

c#define NEW_API
c#define MIRROR
#define GA3
#define NGA_GATSCAT
c#define BLOCK_CYCLIC
c#define USE_SCALAPACK_DISTR
c#define USE_RESTRICTED

#ifdef USE_RESTRICTED
#  define NEW_API
#endif

#define MEM_INC 1000

#ifdef BLOCK_CYCLIC
#  define NEW_API
#  undef MIRROR
#else
#  undef USE_SCALAPAC_DISTR
#endif

      program main
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      integer heap, stack, fudge, ma_heap, me, nproc, map(4096), block
      integer g_s, ndim, dim1, i
      logical status
      parameter (heap=200*200*4, fudge=100, stack=200*200)
c     
c***  Intitialize a message passing library
c
#include "mp3.fh"
c
c***  Initialize GA
c
c     There are 2 choices: ga_initialize or ga_initialize_ltd.
c     In the first case, there is no explicit limit on memory usage.
c     In the second, user can set limit (per processor) in bytes.
c
      call ga_initialize()
      nproc = ga_nnodes()
      me = ga_nodeid()
c     we can also use GA_set_memory_limit BEFORE first ga_create call
c
      ma_heap = heap/nproc + fudge 
      ma_heap = 2*ma_heap
#ifdef USE_RESTRICTED
      ma_heap = 2*ma_heap
#endif
      call GA_set_memory_limit(util_mdtob(ma_heap))
c
      if(ga_nodeid().eq.0)then
#ifdef MIRROR
         print *,' Performing tests on Mirrored Arrays'
#endif
         print *,' GA initialized '
         call ffflush(6)
      endif
c
c***  Initialize the MA package
c     MA must be initialized before any global array is allocated
c
      status = ma_init(MT_DCPL, stack, ma_heap)
      if (.not. status) call ga_error('ma_init failed',-1) 
c
c     Uncomment the below line to register external memory allocator
c     for dynamic arrays inside GA routines.
c      call register_ext_memory()
c
      if(me.eq.(nproc-1))then
        print *, 'using ', nproc,' process(es) ', ga_cluster_nnodes(),
     $           ' cluster nodes'
        print *,'process ', me, ' is on node ',ga_cluster_nodeid(),
     $          ' with ',  ga_cluster_nprocs(-1), ' processes' 
        call ffflush(6)
      endif
c
c   create array to force staggering of memory and uneven distribution
c   of pointers
c
      dim1 = MEM_INC
      map(1) = 1
      do i = 2, nproc
        map(i) = MEM_INC*(i-1)+1
        dim1 = dim1 + MEM_INC*i
      end do
      g_s = ga_create_handle()
      ndim = 1
      call ga_set_data(g_s,ndim,dim1,MT_INT)
      call ga_set_array_name(g_s,'s')
      call ga_set_irreg_distr(g_s,map,nproc)

c
c***  Check support for single precision complex arrays
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' CHECKING SINGLE COMPLEX  '
         write(6,*)
         call ffflush(6)
      endif

      call check_complex_float()
c
c***  Check support for double precision complex arrays
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' CHECKING DOUBLE COMPLEX  '
         write(6,*)
         call ffflush(6)
      endif

      call check_complex()
      if(me.eq.0) call ga_print_stats()
      if(me.eq.0) print *,' ' 
      if(me.eq.0) print *,'All tests successful ' 
      status = ga_destroy(g_s)
c
c***  Tidy up the GA package
c
      call ga_terminate()
c
c***  Tidy up after message-passing library
c
      call MP_FINALIZE()
c
      end




      subroutine check_complex_float()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c
      integer n,m
      parameter (n = 60)
      parameter (m = 2*n)
      complex a(n,n), b(n,n), v(m),w(m)
#ifdef MIRROR
      integer ndim, dims(2), chunk(2), p_mirror
#else
#  ifdef NEW_API
      integer ndim, dims(2), chunk(2), p_mirror
#  endif
#endif
      integer iv(m), jv(m)
      logical status
      integer g_a, g_b
      integer iran, i,j, loop,nloop,maxloop, ilo, ihi, jlo, jhi, itmp
      integer nproc, me, int, ij, inc, ii, jj, nnodes
      parameter (maxloop = 100)
      integer maxproc
      parameter (maxproc = 4096)
      double precision crap, real
      double precision nwords
      complex   x, sum1, sum2, factor
      integer lprocs, inode, iproc, lproc
      integer scpl_type, istatus
#ifdef USE_RESTRICTED
      integer num_rstrctd
      integer rstrctd_list(maxproc/2)
#endif
      intrinsic int
      iran(i) = int(drand(0)*real(i)) + 1
c
      nproc = ga_nnodes()
      me    = ga_nodeid()
      inode = ga_cluster_nodeid()
      lprocs = ga_cluster_nprocs(inode)
      nnodes = ga_cluster_nnodes()
      iproc = mod(me,lprocs)
      nloop = Min(maxloop,n)
#ifdef USE_RESTRICTED
      num_rstrctd = nproc/2
      if (num_rstrctd.eq.0) num_rstrctd = 1
      do i = 1, num_rstrctd
        rstrctd_list(i) = (num_rstrctd/2) + i-1
      end do
#endif
c
c     a() is a local copy of what the global array should start as
c
      do j = 1, n
         do i = 1, n
#ifndef MIRROR
            a(i,j) = cmplx(real(i-1), real((j-1)*n))
#else
            a(i,j) = cmplx(real(inode),0.0d00)
     +             + cmplx(real(i-1), real((j-1)*n))
#endif
            b(i,j) = cmplx(-1d0,1d0)
         enddo
      enddo
c
c     Create type
c
      scpl_type = nga_register_type(8)
c
c     Create a global array
c
c     print *,ga_nodeid(), ' creating array'
      call ffflush(6)
c     call setdbg(1)
#ifdef NEW_API
      ndim = 2
      dims(1) = n
      dims(2) = n
      g_a = ga_create_handle()
      call ga_set_data(g_a,ndim,dims,scpl_type)
      call ga_set_array_name(g_a,'a')
#ifdef USE_RESTRICTED
      call ga_set_restricted(g_a, rstrctd_list, num_rstrctd)
#endif
#  ifdef MIRROR
      p_mirror = ga_pgroup_get_mirror()
      call ga_set_pgroup(g_a,p_mirror)
#  endif
      status = ga_allocate(g_a)
#else
#  ifndef MIRROR
      status = ga_create(scpl_type, n, n, 'a', 0, 0, g_a)
#  else
      ndim = 2
      dims(1) = n
      dims(2) = n
      chunk(1) = 0
      chunk(2) = 0
      p_mirror = ga_pgroup_get_mirror()
      status = nga_create_config(scpl_type, ndim, dims, 'a', chunk,
     +                           p_mirror, g_a)
#  endif
#endif
      if (.not. status) then
         write(6,*) ' ga_create failed'
         call ga_error('... exiting ',0)
      endif
#ifdef NEW_API
      g_b = ga_create_handle()
      call ga_set_data(g_b,ndim,dims,scpl_type)
      call ga_set_array_name(g_b,'b')
#  ifdef MIRROR
      call ga_set_pgroup(g_b,p_mirror)
#  endif
      if (.not.ga_allocate(g_b)) then
#else
#  ifndef MIRROR
      if (.not. ga_create(scpl_type, n, n, 'b', 0, 0, g_b)) then
#  else
      if (.not. nga_create_config(scpl_type, ndim, dims, 'b', chunk,
     _                            p_mirror, g_b)) then
#  endif
#endif
         call ga_error('ga_create failed for second array ',0)
      endif

#ifndef MIRROR
      call ga_distribution(g_a,me,ilo, ihi, jlo, jhi)
#else
      lproc = me - ga_cluster_procid(inode,0)
      call ga_distribution(g_a, lproc, ilo, ihi, jlo, jhi)
#endif
      call ga_sync()
c
c     Zero the array
c
      if (me .eq. 0) then
         write(6,21)
 21      format(/'> Checking zero ... ')
         call ffflush(6)
      endif
      call ga_zero(g_a)
c
c     Check that it is indeed zero
c
      call ga_get(g_a, 1, n, 1, n, b, n)
      call ga_sync()
      do i = 1, n
         do j = 1, n
            if(b(i,j).ne.(0d0,0d0)) then
               write(6,*) me,' zero ', i, j, b(i,j)
               call ffflush(6)
c              call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' ga_zero is OK'
         write(6,*)
      endif
      call ga_sync()
c
c     Each node fills in disjoint sections of the array
c
      if (me .eq. 0) then
         write(6,2)
 2       format(/'> Checking disjoint put ... ')
         call ffflush(6)
      endif
      call ga_sync()
c
      inc = (n-1)/4 + 1
      ij = 0
      do j = 1, n, inc
         do i = 1, n, inc
#ifndef MIRROR
            if (mod(ij,nproc) .eq. me) then
#else
            if (mod(ij,lprocs) .eq. iproc) then
#endif
               ilo = i
               ihi = min(i+inc, n)
               jlo = j
               jhi = min(j+inc, n)
               call ga_put(g_a, ilo, ihi, jlo, jhi, a(ilo, jlo), n)
            endif
            ij = ij + 1
         enddo
      enddo
      call ga_sync()
c
c     All nodes check all of a
c
      call ga_get(g_a, 1, n, 1, n, b, n)
c
      do i = 1, n
         do j = 1, n
            if (b(i,j) .ne. a(i,j)) then
               write(6,*) ' put ', me, i, j, a(i,j),b(i,j)
               call ffflush(6)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' ga_put is OK'
         write(6,*)
      endif
      call ga_sync()
c
c     Now check nloop random gets from each node
c
      if (me .eq. 0) then
         write(6,5) nloop
 5       format(/'> Checking random get (',i5,' calls)...')
         call ffflush(6)
      endif
      call ga_sync()
c
      nwords = 0
c
      crap = drand(ga_nodeid()*51 +1) !Different seed for each proc
      do loop = 1, nloop
         ilo = iran(loop)
         ihi = iran(loop)
         if (ihi.lt. ilo) then
            itmp = ihi
            ihi = ilo
            ilo = itmp
         endif
         jlo = iran(loop)
         jhi = iran(loop)
         if (jhi.lt. jlo) then
            itmp = jhi
            jhi = jlo
            jlo = itmp
         endif
c
         nwords = nwords + (ihi-ilo+1)*(jhi-jlo+1)
c
         call ga_get(g_a, ilo, ihi, jlo, jhi, b(ilo, jlo), n)
         if (me .eq. 0 .and. mod(loop-1, max(1,nloop/20)).eq.0) then
            write(6,1) loop, me, ilo, ihi, jlo, jhi, nwords
 1          format(' call ',i5, ' node ',i2,' checking get ',4i4,
     $           ' total ',d9.2)
            call ffflush(6)
         endif
         do j = jlo, jhi
            do i = ilo, ihi
               if (b(i,j) .ne. a(i,j)) then
                  write(6,*)'error:', i, j, b(i,j), a(i,j)
                  call ga_error('... exiting ',0)
               endif
            enddo
         enddo
c
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' ga_get is OK'
         write(6,*)
         call ffflush(6)
      endif
      call ga_sync()
c
c     Check the ga_copy function
c
      if (me .eq. 0) then
         write(6,*)
         write(6,*)'> Checking copy'
         write(6,*)
         call ffflush(6)
      endif
      call ga_sync()
#ifndef MIRROR
      if(me.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
#else
      if(iproc.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
#endif
      call ga_copy(g_a, g_b)
      call ga_get(g_b, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if (b(i,j) .ne. a(i,j)) then
               write(6,*) ' copy ', me, i, j, a(i,j), b(i,j)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' copy is OK '
         write(6,*)
      endif
c
c     Check scatter&gather
c
      call ga_sync()
      if (me .eq. 0) then
         write(6,*) '> Checking scatter/gather (might be slow)... '
         call ffflush(6)
         if(me.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
      endif
      call ga_sync()
c
      crap = drand(ga_nodeid()*51 +1) !Different seed for each proc
      do j = 1, 10
       call ga_sync()
#ifndef MIRROR
       itmp = iran(nproc)-1
       if(me.eq.itmp) then
#else
       itmp = iran(lprocs)-1
       if(iproc.eq.itmp) then
#endif
         do loop = 1,m
           ilo = iran(n)
           jlo = iran(n)
           iv(loop) = ilo
           jv(loop) = jlo
         enddo
         call ga_gather(g_a, v, iv, jv, m)
         do loop = 1,m
           ilo= iv(loop)
           jlo= jv(loop)
           call ga_get(g_a,ilo,ilo,jlo,jlo,v(loop),1)
           if(v(loop)  .ne. a(ilo,jlo))then
             write(6,*)me,' gather ', ilo,',',jlo,',', a(ilo,jlo)
     &             ,' ',v(loop)
             call ffflush(6)
             call ga_error('... exiting ',0)
           endif
         enddo
       endif
      enddo
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' gather is  OK'
         write(6,*)
         call ffflush(6)
      endif
c
      do j = 1,10
       call ga_sync()
#ifndef MIRROR
       if(me.eq.iran(ga_nnodes())-1) then
#else
       if(me.eq.iran(lprocs)-1) then
#endif
         do loop = 1,m
           ilo = iran(n)
           jlo = iran(n)
           iv(loop) = ilo
           jv(loop) = jlo
           v(loop) = (1d0,-1d0) *(ilo+jlo)
         enddo
         call ga_scatter(g_a, v, iv, jv, m)
         do loop = 1,m
           ilo= iv(loop)
           jlo= jv(loop)
           call ga_get(g_a,ilo,ilo,jlo,jlo,w(loop),1)
           if(w(loop)  .ne. (1d0,-1d0) *(ilo+jlo) )then
             write(6,*)me,' scatter ', ilo,',',jlo,',',w(loop)
     &             ,' ', (1d0,-1d0) *(ilo+jlo)
             call ffflush(6)
           endif
         enddo
       endif
       call ga_sync()
      enddo
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' scatter is  OK'
         write(6,*)
      endif
c
c     Delete the global arrays
c

      status = ga_destroy(g_b)
      status = ga_destroy(g_a)
c
      istatus = nga_deregister_type(scpl_type)
c
      end
c-----------------------------------------------------------------

      subroutine check_complex()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c
      integer n,m
      parameter (n = 60)
      parameter (m = 2*n)
      double complex a(n,n), b(n,n), v(m),w(m)
#ifdef MIRROR
      integer ndim, dims(2), chunk(2), p_mirror
#else
#  ifdef NEW_API
      integer ndim, dims(2), chunk(2), p_mirror
#  endif
#endif
      integer iv(m), jv(m)
      logical status
      integer g_a, g_b
      integer iran, i,j, loop,nloop,maxloop, ilo, ihi, jlo, jhi, itmp
      integer nproc, me, int, ij, inc, ii, jj, nnodes
      parameter (maxloop = 100)
      integer maxproc
      parameter (maxproc = 4096)
      double precision crap, real
      double precision nwords
      double complex   x, sum1, sum2, factor
      integer lprocs, inode, iproc, lproc
      integer dcpl_type, istatus
#ifdef USE_RESTRICTED
      integer num_rstrctd
      integer rstrctd_list(maxproc/2)
#endif
#ifdef BLOCK_CYCLIC
      integer block_size(2), proc_grid(2)
#endif
      intrinsic int
      iran(i) = int(drand(0)*real(i)) + 1
c
      nproc = ga_nnodes()
      me    = ga_nodeid()
      inode = ga_cluster_nodeid()
      lprocs = ga_cluster_nprocs(inode)
      nnodes = ga_cluster_nnodes()
      iproc = mod(me,lprocs)
      nloop = Min(maxloop,n)
#ifdef USE_RESTRICTED
      num_rstrctd = nproc/2
      if (num_rstrctd.eq.0) num_rstrctd = 1
      do i = 1, num_rstrctd
        rstrctd_list(i) = (num_rstrctd/2) + i-1
      end do
#endif
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
c     a() is a local copy of what the global array should start as
c
      do j = 1, n
         do i = 1, n
#ifndef MIRROR
            a(i,j) = cmplx(dble(i-1), dble((j-1)*n))
#else
            a(i,j) = cmplx(dble(inode),0.0d00)
     +             + cmplx(dble(i-1), dble((j-1)*n))
#endif
            b(i,j) = cmplx(-1d0,1d0)
         enddo
      enddo
c
c     Create  type
c
      dcpl_type = nga_register_type(16)
c
c     Create a global array
c
c     print *,ga_nodeid(), ' creating array'
      call ffflush(6)
c     call setdbg(1)
#ifdef NEW_API
      ndim = 2
      dims(1) = n
      dims(2) = n
      g_a = ga_create_handle()
      call ga_set_data(g_a,ndim,dims,dcpl_type)
      call ga_set_array_name(g_a,'a')
#ifdef USE_RESTRICTED
      call ga_set_restricted(g_a, rstrctd_list, num_rstrctd)
#endif
#ifdef BLOCK_CYCLIC
#ifdef USE_SCALAPACK_DISTR
      call ga_set_block_cyclic_proc_grid(g_a,block_size,proc_grid)
#else
      call ga_set_block_cyclic(g_a,block_size)
#endif
#endif
#  ifdef MIRROR
      p_mirror = ga_pgroup_get_mirror()
      call ga_set_pgroup(g_a,p_mirror)
#  endif
      status = ga_allocate(g_a)
#else
#  ifndef MIRROR
      status = ga_create(dcpl_type, n, n, 'a', 0, 0, g_a)
#  else
      ndim = 2
      dims(1) = n
      dims(2) = n
      chunk(1) = 0
      chunk(2) = 0
      p_mirror = ga_pgroup_get_mirror()
      status = nga_create_config(dcpl_type, ndim, dims, 'a', chunk,
     +                           p_mirror, g_a)
#  endif
#endif
      if (.not. status) then
         write(6,*) ' ga_create failed'
         call ga_error('... exiting ',0)
      endif
#ifdef NEW_API
      g_b = ga_create_handle()
      call ga_set_data(g_b,ndim,dims,dcpl_type)
      call ga_set_array_name(g_b,'b')
#ifdef BLOCK_CYCLIC
#ifdef USE_SCALAPACK_DISTR
      call ga_set_block_cyclic_proc_grid(g_b,block_size,proc_grid)
#else
      call ga_set_block_cyclic(g_b,block_size)
#endif
#endif
#  ifdef MIRROR
      call ga_set_pgroup(g_b,p_mirror)
#  endif
      if (.not.ga_allocate(g_b)) then
#else
#  ifndef MIRROR
      if (.not. ga_create(dcpl_type, n, n, 'b', 0, 0, g_b)) then
#  else
      if (.not. nga_create_config(dcpl_type, ndim, dims, 'b', chunk,
     _                            p_mirror, g_b)) then
#  endif
#endif
         call ga_error('ga_create failed for second array ',0)
      endif

#ifndef MIRROR
      call ga_distribution(g_a,me,ilo, ihi, jlo, jhi)
#else
      lproc = me - ga_cluster_procid(inode,0)
      call ga_distribution(g_a, lproc, ilo, ihi, jlo, jhi)
#endif
      call ga_sync()
c
c     Zero the array
c
      if (me .eq. 0) then
         write(6,21)
 21      format('> Checking zero ... ')
         call ffflush(6)
      endif
      call ga_zero(g_a)
c
c     Check that it is indeed zero
c
      call ga_get(g_a, 1, n, 1, n, b, n)
      call ga_sync()
      do i = 1, n
         do j = 1, n
            if(b(i,j).ne.(0d0,0d0)) then
               write(6,*) me,' zero ', i, j, b(i,j)
               call ffflush(6)
c              call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' ga_zero is OK'
         write(6,*)
      endif
      call ga_sync()
c
c     Each node fills in disjoint sections of the array
c
      if (me .eq. 0) then
         write(6,2)
 2       format('> Checking disjoint put ... ')
         call ffflush(6)
      endif
      call ga_sync()
c
      inc = (n-1)/20 + 1
      ij = 0
      do j = 1, n, inc
         do i = 1, n, inc
#ifndef MIRROR
            if (mod(ij,nproc) .eq. me) then
#else
            if (mod(ij,lprocs) .eq. iproc) then
#endif
               ilo = i
               ihi = min(i+inc, n)
               jlo = j
               jhi = min(j+inc, n)
               call ga_put(g_a, ilo, ihi, jlo, jhi, a(ilo, jlo), n)
            endif
            ij = ij + 1
         enddo
      enddo
      call ga_sync()
c
c     All nodes check all of a
c
      call util_qfill(n*n, (0d0,0d0), b, 1)
      call ga_get(g_a, 1, n, 1, n, b, n)
c
      do i = 1, n
         do j = 1, n
            if (b(i,j) .ne. a(i,j)) then
               write(6,*) ' put ', me, i, j, a(i,j),b(i,j)
               call ffflush(6)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' ga_put is OK'
         write(6,*)
      endif
      call ga_sync()
c
c     Now check nloop random gets from each node
c
      if (me .eq. 0) then
         write(6,5) nloop
 5       format('> Checking random get (',i5,' calls)...')
         call ffflush(6)
      endif
      call ga_sync()
c
      nwords = 0
c
      crap = drand(ga_nodeid()*51 +1) !Different seed for each proc
      do loop = 1, nloop
         ilo = iran(loop)
         ihi = iran(loop)
         if (ihi.lt. ilo) then
            itmp = ihi
            ihi = ilo
            ilo = itmp
         endif
         jlo = iran(loop)
         jhi = iran(loop)
         if (jhi.lt. jlo) then
            itmp = jhi
            jhi = jlo
            jlo = itmp
         endif
c
         nwords = nwords + (ihi-ilo+1)*(jhi-jlo+1)
c
         call util_qfill(n*n, (0.0d0,0d0), b, 1)
         call ga_get(g_a, ilo, ihi, jlo, jhi, b(ilo, jlo), n)
         if (me .eq. 0 .and. mod(loop-1, max(1,nloop/20)).eq.0) then
            write(6,1) loop, me, ilo, ihi, jlo, jhi, nwords
 1          format(' call ',i5, ' node ',i2,' checking get ',4i4,
     $           ' total ',d9.2)
            call ffflush(6)
         endif
         do j = jlo, jhi
            do i = ilo, ihi
               if (b(i,j) .ne. a(i,j)) then
                  write(6,*)'error:', i, j, b(i,j), a(i,j)
                  call ga_error('... exiting ',0)
               endif
            enddo
         enddo
c
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' ga_get is OK'
         write(6,*)
         call ffflush(6)
      endif
      call ga_sync()
c
c     Check the ga_copy function
c
      if (me .eq. 0) then
         write(6,*)
         write(6,*)'> Checking copy'
         write(6,*)
         call ffflush(6)
      endif
      call ga_sync()
#ifndef MIRROR
      if(me.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
#else
      if(iproc.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
#endif
      call ga_copy(g_a, g_b)
      call ga_get(g_b, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if (b(i,j) .ne. a(i,j)) then
               write(6,*) ' copy ', me, i, j, a(i,j), b(i,j)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' copy is OK '
         write(6,*)
      endif
c
c     Check scatter&gather
c
      call ga_sync()
      if (me .eq. 0) then
         write(6,*) '> Checking scatter/gather (might be slow)... '
         call ffflush(6)
         if(me.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
      endif
      call ga_sync()
c
      crap = drand(ga_nodeid()*51 +1) !Different seed for each proc
      do j = 1, 10
       call ga_sync()
#ifndef MIRROR
       itmp = iran(nproc)-1
       if(me.eq.itmp) then
#else
       itmp = iran(lprocs)-1
       if(iproc.eq.itmp) then
#endif
         do loop = 1,m
           ilo = iran(n)
           jlo = iran(n)
           iv(loop) = ilo
           jv(loop) = jlo
         enddo
         call ga_gather(g_a, v, iv, jv, m)
         do loop = 1,m
           ilo= iv(loop)
           jlo= jv(loop)
           call ga_get(g_a,ilo,ilo,jlo,jlo,v(loop),1)
           if(v(loop)  .ne. a(ilo,jlo))then
             write(6,*)me,' gather ', ilo,',',jlo,',', a(ilo,jlo)
     &             ,' ',v(loop)
             call ffflush(6)
             call ga_error('... exiting ',0)
           endif
         enddo
       endif
      enddo
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' gather is  OK'
         write(6,*)
         call ffflush(6)
      endif
c
      do j = 1,10
       call ga_sync()
#ifndef MIRROR
       if(me.eq.iran(ga_nnodes())-1) then
#else
       if(me.eq.iran(lprocs)-1) then
#endif
         do loop = 1,m
           ilo = iran(n)
           jlo = iran(n)
           iv(loop) = ilo
           jv(loop) = jlo
           v(loop) = (1d0,-1d0) *(ilo+jlo)
         enddo
         call ga_scatter(g_a, v, iv, jv, m)
         do loop = 1,m
           ilo= iv(loop)
           jlo= jv(loop)
           call ga_get(g_a,ilo,ilo,jlo,jlo,w(loop),1)
           if(w(loop)  .ne. (1d0,-1d0) *(ilo+jlo) )then
             write(6,*)me,' scatter ', ilo,',',jlo,',',w(loop)
     &             ,' ', (1d0,-1d0) *(ilo+jlo)
             call ffflush(6)
           endif
         enddo
       endif
       call ga_sync()
      enddo
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' scatter is  OK'
         write(6,*)
      endif
c
c     Delete the global arrays
c
      status = ga_destroy(g_b)
      status = ga_destroy(g_a)
c
      istatus = nga_deregister_type(dcpl_type)
c
      end

      subroutine util_qfill(n,val,a,ia)
      implicit none
      double  complex  a(*), val
      integer n, ia, i
c
c     initialise double complex array to scalar value
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

