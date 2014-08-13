#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c $Id: pg2test.F,v 1.6 2004-12-01 22:39:46 manoj Exp $
c vector boxes lack arithmetic precision 
#ifdef CRAY_YMP
# define THRESH 1d-10
# define THRESHF 1e-5
#elif defined(FUJITSU)
# define THRESH 1d-12
# define THRESHF 1e-5
#else
# define THRESH 1d-13
# define THRESHF 1e-5
#endif

#define MISMATCH(x,y) abs(x-y)/max(1d0,abs(x)).gt.THRESH
#define MISMATCHF(x,y) abs(x-y)/max(1.0,abs(x)).gt.THRESHF 

#define NEW_API
c#define GA3
#define NGA_GATSCAT

      program main
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      integer heap, stack, fudge, ma_heap, me, nproc
      integer inode,proclist(100),i,j,nprocs
      integer proc_group(0:100)
      logical status
      parameter (heap=200*200*4, fudge=100, stack=200*200)
      integer n,m,grp_id
      parameter (n = 128)
      integer ndim, dims(2), chunk(2)
      integer g_w, midproc
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
c
c***  Initialize the MA package
c     MA must be initialized before any global array is allocated
c
      ma_heap = heap/nproc + fudge 
      ma_heap = 2*ma_heap
      status = ma_init(MT_DCPL, stack, ma_heap)
      if (.not. status) call ga_error('ma_init failed',-1) 
c
c     we can also use GA_set_memory_limit BEFORE first ga_create call
c
      call GA_set_memory_limit(util_mdtob(ma_heap))
c
c***  Create a global array on the current default (world) group
c
      ndim = 2
      dims(1) = n
      dims(2) = n
      g_w = ga_create_handle()
      call ga_set_data(g_w,ndim,dims,MT_DBL)
      call ga_set_array_name(g_w,'w')
      status = ga_allocate(g_w)
      call ga_zero(g_w)
c
      midproc = nproc/2
      do i = 1, midproc
        proclist(i) = i-1
      end do
      nprocs = midproc
      proc_group(0) = ga_pgroup_create(proclist, nprocs)
      do i = midproc + 1, nproc
        proclist(i-midproc) = i-1
      end do
      nprocs = nproc - midproc
      proc_group(1) = ga_pgroup_create(proclist, nprocs)

c
      write(6,*) 'Groups have been created'
      call ga_sync()
      if (me.ge.midproc) then
          call runtest(proc_group(1), g_w)
      else
      endif
c
c***  Check if memory limits are enforced 
c
c      call check_mem(util_mdtob(ma_heap*ga_nnodes()))
c
c***  Tidy up the GA package
c
      call ga_sync()
      write(0,*) 'Calling ga_terminate'
      call ga_terminate()
c
c***  Tidy up after message-passing library
c
      call MP_FINALIZE()
c
      stop
      end
c
      subroutine runtest(grp_id, g_w)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      integer me, nproc, grp_id, g_w
      logical status

      me = ga_pgroup_nodeid(grp_id)
      nproc = ga_pgroup_nnodes(grp_id)

      if(me.eq.0)then
         print *,' GA initialized '
         call ffflush(6)
      endif
c
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
      call ga_pgroup_sync(grp_id)
c
c***  Check support for double precision arrays
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' CHECKING DOUBLES  '
         write(6,*)
         call ffflush(6)
      endif

      call check_dbl(grp_id, g_w)
c
c***  Check support for double precision complex arrays
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' CHECKING DOUBLE COMPLEX  '
         write(6,*)
         call ffflush(6)
      endif

      call check_complex(grp_id)
c
c***  Check support for integer arrays
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' CHECKING INTEGERS  '
         write(6,*)
         call ffflush(6)
      endif

      call check_int(grp_id)
c
c
c***  Check support for single precision
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' CHECKING SINGLE PRECISION  '
         write(6,*)
         call ffflush(6)
      endif
 
      call check_flt(grp_id) 
c     
      if (me.eq.0) then
         write(6,*)
         write(6,*)'CHECKING Wrappers to Message Passing Collective ops'
         write(6,*)
         call ffflush(6)
      endif

      call check_wrappers(grp_id)
c
      if(me.eq.0) call ga_print_stats()
      if(me.eq.0) print *,' ' 
      if(me.eq.0) print *,'All tests succesful ' 
c
      return
      end


      subroutine check_dbl(grp_id, g_w)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c     
      integer n,m,grp_id
      parameter (n = 128)
      parameter (m = 2*n)
      double precision a(n,n), b(n,n), v(m),w(m)
      integer ndim, dims(2), chunk(2)
      integer iv(m), jv(m)
#ifdef NGA_GATSCAT
      integer ijv(2,m)
#endif
      logical status
      integer g_a, g_b, g_w
      integer iran, i,j, loop,nloop,maxloop, ilo, ihi, jlo, jhi, itmp
      integer nproc, me, int, ij, inc, ii, jj
      parameter (maxloop = 100)
      integer maxproc
      parameter (maxproc = 128)
      double precision crap, sum1, sum2, x
      double precision nwords
      logical found
      integer lprocs,inode,iproc,nnodes
      intrinsic int
      iran(i) = int(drand(0)*real(i)) + 1
c     
      nproc = ga_pgroup_nnodes(grp_id)
      me    = ga_pgroup_nodeid(grp_id)
      lprocs = nproc
      inode = me
      nnodes = ga_cluster_nnodes()
      iproc = me
      nloop = Min(maxloop,n)
c     
c     a() is a local copy of what the global array should start as
c
      do j = 1, n
         do i = 1, n
            a(i,j) = i-1 + (j-1)*n
            b(i,j) =-1. 
         enddo
      enddo
*      write(6,*) ' correct '
*      call output(a, 1, n, 1, n, n, n, 1)
*      call ffflush(6)
c
c     Create a global array
c
c     print *,ga_nodeid(), ' creating array'
*     call ffflush(6)
c     call setdbg(1)
      g_a = ga_create_handle()
      ndim = 2
      dims(1) = n
      dims(2) = n
      call ga_set_data(g_a,ndim,dims,MT_DBL)
      call ga_set_array_name(g_a,'a')
      write(6,*) 'grp_id ',grp_id
      call ga_set_pgroup(g_a,grp_id)
      status = ga_allocate(g_a)
      if (.not. status) then
         write(6,*) ' ga_create failed'
         call ffflush(6)
         call ga_error('... exiting ',0)
      endif
c
c     check if handle is valid. be quiet unless error
c
      if(.not.ga_valid_handle(g_a)) call ga_error("invalid handle",g_a)
c
      call ga_distribution(g_a,me,ilo, ihi, jlo, jhi)
      call ga_pgroup_sync(grp_id)
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
      call ga_pgroup_sync(grp_id)
      do i = 1, n
         do j = 1, n
            if (b(i,j) .ne. 0.0d0) then
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
      call ga_pgroup_sync(grp_id)
c
c     Each node fills in disjoint sections of the array
c
      if (me .eq. 0) then
         write(6,2) 
 2       format(/'> Checking disjoint put ... ')
         call ffflush(6)
      endif
      call ga_pgroup_sync(grp_id)
c
      inc = (n-1)/20 + 1
      inc = n
      ij = 0
      do j = 1, n, inc
         do i = 1, n, inc
            if (mod(ij,nproc) .eq. me) then
               ilo = i
               ihi = min(i+inc, n)
               jlo = j
               jhi = min(j+inc, n)
*               write(6,4) me, ilo, ihi, jlo, jhi
* 4             format(' node ',i2,' checking put ',4i4)
*               call ffflush(6)
               call ga_put(g_a, ilo, ihi, jlo, jhi, a(ilo, jlo), n)
            endif
            ij = ij + 1
         enddo
      enddo
      call ga_pgroup_sync(grp_id)
c
c     All nodes check all of a
c
      call util_dfill(n*n, 0.0d0, b, 1)
*     call ga_print(g_a,1)
      call ga_get(g_a, 1, n, 1, n, b, n)
*      write(6,*) ' after get'
*      call output(b, 1, n, 1, n, n, n, 1)
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
      call ga_pgroup_sync(grp_id)
c
c     Now check nloop random gets from each node
c
      if (me .eq. 0) then
         write(6,5) nloop
 5       format(/'> Checking random get (',i5,' calls)...')
         call ffflush(6)
      endif
      call ga_pgroup_sync(grp_id)
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
         call util_dfill(n*n, 0.0d0, b, 1)
         call ga_get(g_a, ilo, ihi, jlo, jhi, b(ilo, jlo), n)
         if (me .eq. 0 .and. mod(loop-1, max(1,nloop/20)).eq.0) then
            write(6,1) loop, me, ilo, ihi, jlo, jhi, nwords
 1          format(' call ',i5, ' node ',i2,' checking get ',4i4,
     $           ' total ',d9.2)
            call ffflush(6)
         endif
         sum1 = 0.0d0
         do j = jlo, jhi
            do i = ilo, ihi
               sum1 = sum1 + b(i,j)
               if (b(i,j) .ne. a(i,j)) then
                  write(6,*) i, j, b(i,j), a(i,j)
                  call ffflush(6)
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
      call ga_pgroup_sync(grp_id)
c
c     Each node accumulates into disjoint sections of the array
c
      if (me .eq. 0) then
         write(6,9) 
 9       format(/'> Checking accumulate ... ')
         call ffflush(6)
      endif
      call ga_pgroup_sync(grp_id)
c     
      crap = drand(12345)       ! Same seed for each process
      do j = 1, n
         do i = 1, n
c           b(i,j) = drand(0)
            b(i,j) = i+j 
         enddo
      enddo
c
      inc = (n-1)/20 + 1
      ij = 0
      do j = 1, n, inc
         do i = 1, n, inc
c           x = drand(0)
            x = 10.
            ilo = i
            ihi = min(i+inc-1, n)
            if(ihi.eq.n-1)ihi=n
c           ihi = min(i+inc, n)
            jlo = j
            jhi = min(j+inc-1, n)
            if(jhi.eq.n-1)jhi=n
c           jhi = min(j+inc-1, n)
*               call ffflush(6)
            if (mod(ij,nproc) .eq. me) then
c               print *, me, 'checking accumulate ',ilo,ihi,jlo,jhi,x
* 11            format(' node ',i2,' checking accumulate ',4i4)
*               call ffflush(6)
               call ga_acc(g_a, ilo, ihi, jlo, jhi, b(ilo, jlo), n, x)
            endif
            ij = ij + 1
c     
c           Each process applies all updates to its local copy
c
            do jj = jlo, jhi
               do ii = ilo, ihi
                  a(ii,jj) = a(ii,jj) + x * b(ii,jj)
               enddo
            enddo
         enddo
      enddo
      call ga_pgroup_sync(grp_id)
c
c     All nodes check all of a
c
      call ga_get(g_a, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if(MISMATCH(b(i,j),a(i,j)))then
               write(6,*) ' acc ', me, i, j, a(i,j), b(i,j)
               call ffflush(6)
c              call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' disjoint ga_acc is OK'
         write(6,*)
      endif
c
c     overlapping accumulate
      call ga_pgroup_sync(grp_id)
      g_b = ga_create_handle()
      call ga_set_data(g_b,ndim,dims,MT_DBL)
      call ga_set_array_name(g_b,'b')
      call ga_set_pgroup(g_b,grp_id)
      if (.not.ga_allocate(g_b)) then
         call ga_error('ga_create failed for second array ',0)
      endif
c
      call ga_zero(g_b)
      call ga_pgroup_sync(grp_id)
      call ga_acc(g_b, n/2, n/2, n/2, n/2, 1d0, 1, 1d0)
      call ga_pgroup_sync(grp_id)
      if (me.eq.0) then
         call ga_get(g_b, n/2, n/2, n/2, n/2, b(1,1), 1)
         x = abs(b(1,1) -1d0*nproc)
         if(x.gt. 1d-10)then
            write(6,*)'val=',b(1,1),' expected=',nproc, x
            call ga_error('overlapping accumulate failed',0)
         endif
         if (me.eq.0) then
           write(6,*)
           write(6,*) ' overlapping ga_acc is OK'
           write(6,*)
         endif
      endif
c
c     Check the ga_add function
c
      if (me .eq. 0) then
         write(6,91)
 91      format(/'> Checking add ...')
         call ffflush(6)
      endif
c
c
      crap = drand(12345)       ! Everyone has same seed
      do j = 1, n
         do i = 1, n
            b(i,j) = drand(0)
            a(i,j) = 0.1d0*a(i,j) + 0.9d0*b(i,j)
         enddo
      enddo
      
      if (me.eq.0) call ga_put(g_b, 1, n, 1, n, b, n)
      call ga_add(0.1d0, g_a, 0.9d0, g_b, g_b)
      call ga_get(g_b, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if(MISMATCH(b(i,j), a(i,j)))then
               write(6,*) ' add ', me, i, j, a(i,j), b(i,j)
               call ffflush(6)
               call ga_error('... exiting ',0) 
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' add is OK '
         write(6,*)
      endif
      call ga_pgroup_sync(grp_id)
c
c     Check the ddot function
c
      if (me .eq. 0) then
         write(6,19) 
 19      format(/'> Checking ddot ...')
         call ffflush(6)
      endif
      crap = drand(12345)       ! Everyone has same seed
      sum1 = 0.0d0
      do j = 1, n
         do i = 1, n
            b(i,j) = drand(0)
            sum1 = sum1 + a(i,j)*b(i,j)
         enddo
      enddo
      if (me.eq.0) then 
         call ga_put(g_b, 1, n, 1, n, b, n)
         call ga_put(g_a, 1, n, 1, n, a, n)
      endif
      call ga_pgroup_sync(grp_id)
      sum2 = ga_ddot(g_a,g_b)
      if(MISMATCH(sum1, sum2))then
         write(6,*) ' ddot wrong ', sum1, sum2
         call ffflush(6)
         call ga_error('... exiting ',0)
      else if (me.eq.0) then
         write(6,*)
         write(6,*) ' ddot is OK '
         write(6,*)
      endif
c
c     Check the ga_scale function
c
      if (me .eq. 0) then
         write(6,92)
 92      format(/'> Checking scale ...')
         call ffflush(6)
      endif
      call ga_scale(g_a, 0.123d0)
      call ga_get(g_a, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            a(i,j) = a(i,j)*0.123d0
            if (MISMATCH(b(i,j), a(i,j)))then
               write(6,*) ' dscal ', me, i, j, a(i,j), b(i,j)
               call ffflush(6)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' scale is OK '
         write(6,*)
      endif
c
c     Check the ga_copy function
c
      if (me .eq. 0) then
         write(6,*)
         write(6,*)'> Checking copy'
         write(6,*)
         call ffflush(6)
      endif
      if(me.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
      call ga_copy(g_a, g_b)
      call ga_get(g_b, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if (b(i,j) .ne. a(i,j)) then
               write(6,*) ' copy ', me, i, j, a(i,j), b(i,j)
               call ffflush(6)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' copy is OK '
         write(6,*)
      endif
      call ga_pgroup_sync(grp_id)
c
c     Check the ga_copy function between two groups
c
      if (me .eq. 0) then
         write(6,*)
         write(6,*)'> Checking copy between groups'
         write(6,*)
         call ffflush(6)
      endif
      if(me.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
      call ga_copy(g_a, g_w)
      call ga_get(g_w, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if (b(i,j) .ne. a(i,j)) then
               write(6,*) ' copy ', me, i, j, a(i,j), b(i,j)
               call ffflush(6)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' copy is OK from small group to big group'
         write(6,*)
      endif
      call ga_pgroup_sync(grp_id)
      call ga_zero(g_a)
      call ga_copy(g_w, g_a)
      call ga_get(g_a, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if (b(i,j) .ne. a(i,j)) then
               write(6,*) ' copy ', me, i, j, a(i,j), b(i,j)
               call ffflush(6)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' copy is OK from big group to small group'
         write(6,*)
      endif
c
c     Test copy_patch
c
      if (me .eq. 0) then
         write(6,*)
         write(6,*)'> Checking copy_patch between groups'
         write(6,*)
         call ffflush(6)
      endif
c
c   klugy mechanism for zeroing g_w
c
      call ga_zero(g_a)
      call ga_copy(g_a,g_w)
      if(me.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
c
      ilo = n/4+1
      jlo = n/4+1
      ihi = ilo + n/2
      jhi = jlo + n/2
      call ga_copy_patch('n',g_a,ilo,ihi,jlo,jhi,g_w,ilo,ihi,jlo,jhi)
      call ga_get(g_w, 1, n, 1, n, b, n)
      do j = jlo, jhi
         do i = ilo, ihi
            if (b(i,j) .ne. a(i,j)) then
               write(6,*) ' copy ', me, i, j, a(i,j), b(i,j)
               call ffflush(6)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' copy_patch is OK from small group to big group'
         write(6,*)
      endif
c
      call ga_zero(g_a)
      call ga_copy_patch('n',g_w,ilo,ihi,jlo,jhi,g_a,ilo,ihi,jlo,jhi)
      call ga_get(g_a, 1, n, 1, n, b, n)
      do j = jlo, jhi
         do i = ilo, ihi
            if (b(i,j) .ne. a(i,j)) then
               write(6,*) ' copy ', me, i, j, a(i,j), b(i,j)
               call ffflush(6)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' copy_patch is OK from big group to small group'
         write(6,*)
      endif
      if(me.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
c
      call ga_pgroup_sync(grp_id)
      if (me .eq. 0) then
         write(6,*) '> Checking scatter/gather (might be slow)... '
         call ffflush(6)
      endif
      call ga_pgroup_sync(grp_id)
c
      crap = drand(ga_nodeid()*51 +1) !Different seed for each proc
      do j = 1, 10 
       call ga_pgroup_sync(grp_id)
       itmp = iran(nproc)-1
       if(me.eq.itmp) then
#ifndef NGA_GATSCAT
         do loop = 1,m 
           ilo = iran(n) 
           jlo = iran(n)
           iv(loop) = ilo
           jv(loop) = jlo
         enddo 
         write(6,*) 'Calling ga_gather'
         call ga_gather(g_a, v, iv, jv, m)
         write(6,*) 'Completed ga_gather'
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
#else
         do loop = 1,m 
           ilo = iran(n) 
           jlo = iran(n)
           ijv(1,loop) = ilo
           ijv(2,loop) = jlo
         enddo 
         call nga_gather(g_a, v, ijv, m)
         do loop = 1,m
           ilo= ijv(1,loop)
           jlo= ijv(2,loop)
           call ga_get(g_a,ilo,ilo,jlo,jlo,v(loop),1)
           if(v(loop)  .ne. a(ilo,jlo))then
             write(6,*)me,' gather ', ilo,',',jlo,',', a(ilo,jlo)
     &             ,' ',v(loop)
             call ffflush(6)
             call ga_error('... exiting ',0)
           endif
         enddo
#endif
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
       call ga_pgroup_sync(grp_id)
       if(me.eq.iran(ga_pgroup_nnodes(grp_id))-1) then
#ifndef NGA_GATSCAT
         do loop = 1,m
           ilo = iran(n)
           jlo = iran(n)
           iv(loop) = ilo
           jv(loop) = jlo
c          v(loop) = DSIN(a(ilo,jlo)+b(ilo,jlo)) 
           v(loop) = 1d0 *(ilo+jlo)
         enddo 
         call ga_scatter(g_a, v, iv, jv, m)
         do loop = 1,m
           ilo= iv(loop)
           jlo= jv(loop) 
           call ga_get(g_a,ilo,ilo,jlo,jlo,w(loop),1)
c          if(v(loop)  .ne. w(loop))then 
           if(w(loop)  .ne. 1d0 *(ilo+jlo) )then 
             write(6,*)me,' scatter ', ilo,',',jlo,',',w(loop) 
     &             ,' ', 1d0 *(ilo+jlo) 
             call ffflush(6)
             call ga_error('... exiting ',0)
           endif
         enddo 
#else
         do loop = 1,m
           ilo = iran(n)
           jlo = iran(n)
           ijv(1,loop) = ilo
           ijv(2,loop) = jlo
c          v(loop) = DSIN(a(ilo,jlo)+b(ilo,jlo)) 
           v(loop) = 1d0 *(ilo+jlo)
         enddo 
         call nga_scatter(g_a, v, ijv, m)
         do loop = 1,m
           ilo= ijv(1,loop)
           jlo= ijv(2,loop) 
           call ga_get(g_a,ilo,ilo,jlo,jlo,w(loop),1)
c          if(v(loop)  .ne. w(loop))then 
           if(w(loop)  .ne. 1d0 *(ilo+jlo) )then 
             write(6,*)me,' scatter ', ilo,',',jlo,',',w(loop) 
     &             ,' ', 1d0 *(ilo+jlo) 
             call ffflush(6)
             call ga_error('... exiting ',0)
           endif
         enddo 
#endif
       endif
       call ga_pgroup_sync(grp_id)
      enddo
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' scatter is  OK'
         write(6,*)
      endif
c
      call ga_pgroup_sync(grp_id)
c
c     scatter-acc available in GA ver. 3.0
#ifdef GA3
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' checking scatter-accumulate'
         write(6,*)
      endif
c
      crap = drand(1234)
      call ga_zero(g_a)
c
      do j = 1, n
         do i = 1, n
            b(i,j) =0.
         enddo
      enddo
c
      x = .1d0
      ii =n
      do jj = 1,1 
        call ga_pgroup_sync(grp_id)
        do loop = 1, ii 

c         generate unique i,j pairs
10        continue
          i = iran(n)
          j=iran(n)
          if (found(i,j, iv, jv, loop-1) ) goto 10
#ifndef NGA_GATSCAT
          iv(loop) = i
          jv(loop) = j
#else
          ijv(1,loop) = i
          ijv(2,loop) = j
#endif
          v(loop) = 1d0 *(i+j)
          b(i,j) = b(i,j) + nproc*x*v(loop) ! update local ref. array
        enddo


#ifndef NGA_GATSCAT
        call ga_scatter_acc(g_a,v,iv,jv, ii,x)
#else
        call nga_scatter_acc(g_a,v,ijv,ii,x)
#endif

c
        call ga_pgroup_sync(grp_id)
c     
c       check the result
c
        call ga_get(g_a, 1, n, 1,n, a, n)

        do loop = 1,ii
#ifndef NGA_GATSCAT
          i = iv(loop)
          j = jv(loop)
#else
          i = ijv(1,loop)
          j = ijv(2,loop)
#endif
          if(MISMATCH(a(i,j),b(i,j)))then
             print *,'Error',i,j,loop,a(i,j),'expected=',b(i,j)  
*            if(me.eq.0)then
*              do ii=1,loop
*                   print *,'element',ii, iv(ii),jv(ii)
*              enddo
*            endif
             call ga_error('scatter-acc error in trial ',jj) 
          endif
        enddo
        call ga_pgroup_sync(grp_id)

      enddo

      call ga_pgroup_sync(grp_id)
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' scatter-accumulate is  OK'
         write(6,*)
      endif
#endif

c
c     Delete the global arrays
c
      status = ga_destroy(g_b)
      status = ga_destroy(g_a)
c     
      return
      end

c-----------------------------------------------------------------
      subroutine check_complex(grp_id)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c
      integer n,m,grp_id
      parameter (n = 60)
      parameter (m = 2*n)
      double complex a(n,n), b(n,n), v(m),w(m)
      integer ndim, dims(2), chunk(2), p_mirror
      integer iv(m), jv(m)
      logical status
      integer g_a, g_b
      integer iran, i,j, loop,nloop,maxloop, ilo, ihi, jlo, jhi, itmp
      integer nproc, me, int, ij, inc, ii, jj, nnodes
      parameter (maxloop = 100)
      integer maxproc
      parameter (maxproc = 128)
      double precision crap, real
      double precision nwords
      double complex   x, sum1, sum2, factor
      integer lprocs, inode, iproc
      intrinsic int
      iran(i) = int(drand(0)*real(i)) + 1
c
      nproc = ga_pgroup_nnodes(grp_id)
      me    = ga_pgroup_nodeid(grp_id)
      inode = ga_cluster_nodeid()
      lprocs = ga_cluster_nprocs(inode)
      nnodes = ga_cluster_nnodes()
      iproc = mod(me,lprocs)
      nloop = Min(maxloop,n)
c
c     a() is a local copy of what the global array should start as
c
      do j = 1, n
         do i = 1, n
            a(i,j) = cmplx(dble(i-1), dble((j-1)*n))
            b(i,j) = cmplx(-1d0,1d0)
         enddo
      enddo
c
c     Create a global array
c
c     print *,ga_nodeid(), ' creating array'
      call ffflush(6)
c     call setdbg(1)
      ndim = 2
      dims(1) = n
      dims(2) = n
      g_a = ga_create_handle()
      call ga_set_data(g_a,ndim,dims,MT_DCPL)
      call ga_set_array_name(g_a,'a')
      call ga_set_pgroup(g_a,grp_id)
      status = ga_allocate(g_a)
      if (.not. status) then
         write(6,*) ' ga_create failed'
         call ga_error('... exiting ',0)
      endif
      g_b = ga_create_handle()
      call ga_set_data(g_b,ndim,dims,MT_DCPL)
      call ga_set_array_name(g_b,'b')
      call ga_set_pgroup(g_b,grp_id)
      if (.not.ga_allocate(g_b)) then
         call ga_error('ga_create failed for second array ',0)
      endif

      call ga_distribution(g_a,me,ilo, ihi, jlo, jhi)
      call ga_pgroup_sync(grp_id)
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
      call ga_pgroup_sync(grp_id)
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
      call ga_pgroup_sync(grp_id)
c
c     Each node fills in disjoint sections of the array
c
      if (me .eq. 0) then
         write(6,2)
 2       format(/'> Checking disjoint put ... ')
         call ffflush(6)
      endif
      call ga_pgroup_sync(grp_id)
c
      inc = (n-1)/20 + 1
      ij = 0
      do j = 1, n, inc
         do i = 1, n, inc
            if (mod(ij,nproc) .eq. me) then
               ilo = i
               ihi = min(i+inc, n)
               jlo = j
               jhi = min(j+inc, n)
               call ga_put(g_a, ilo, ihi, jlo, jhi, a(ilo, jlo), n)
            endif
            ij = ij + 1
         enddo
      enddo
      call ga_pgroup_sync(grp_id)
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
      call ga_pgroup_sync(grp_id)
c
c     Now check nloop random gets from each node
c
      if (me .eq. 0) then
         write(6,5) nloop
 5       format(/'> Checking random get (',i5,' calls)...')
         call ffflush(6)
      endif
      call ga_pgroup_sync(grp_id)
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
      call ga_pgroup_sync(grp_id)
c
c     Each node accumulates into disjoint sections of the array
c
      if (me .eq. 0) then
         write(6,9)
 9       format(/'> Checking accumulate ... ')
         call ffflush(6)
      endif
      call ga_pgroup_sync(grp_id)
c
      crap = drand(12345)       ! Same seed for each process
      do j = 1, n
         do i = 1, n
            b(i,j) = cmplx(drand(0),drand(1))
            b(i,j) = cmplx(dble(i),dble(j))
         enddo
      enddo
c
      inc = (n-1)/20 + 1
      ij = 0
      do j = 1, n, inc
         do i = 1, n, inc
c           x = cmplx(drand(0),0.333d0)
c           x = cmplx(0.333d0,0)
*           x = cmplx(0d0,0d0)
            x = 0 
            ilo = i
            ihi = min(i+inc-1, n)
            if(ihi.eq.n-1)ihi=n
            jlo = j
            jhi = min(j+inc-1, n)
            if(jhi.eq.n-1)jhi=n
            if (mod(ij,nproc) .eq. me) then
               call ga_acc(g_a, ilo, ihi, jlo, jhi, b(ilo, jlo), n, x)
            endif
            ij = ij + 1
c
c           Each process applies all updates to its local copy
c
            do jj = jlo, jhi
               do ii = ilo, ihi
                  a(ii,jj) = a(ii,jj) + x * b(ii,jj)
               enddo
            enddo
         enddo
      enddo
      call ga_pgroup_sync(grp_id)
c
c     All nodes check all of a
c
      call ga_get(g_a, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if (MISMATCH(b(i,j), a(i,j)))then
               write(6,*) ' acc ', me, i, j, a(i,j), b(i,j)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' disjoint ga_acc is OK'
         write(6,*)
      endif
c
c     overlapping accumulate
c
      call ga_zero(g_b)
      call ga_acc(g_b, n/2, n/2, n/2, n/2, (1d0,-1d0), 1, (1d0,0d0))
      call ga_pgroup_sync(grp_id)
      if (me.eq.0) then
         call ga_get(g_b, n/2, n/2, n/2, n/2, x, 1)
         if (MISMATCH(x, ((1d0,-1d0)*nproc)))then
c        if(error.gt. (1d-8))then
            write(6,*)'val=',x,' expected=(',nproc,',',-nproc,')'
            call ga_error('overlapping accumulate failed',0)
         endif
         write(6,*)
         write(6,*) ' overlapping ga_acc is OK'
         write(6,*)
      endif
c
c     Check the ga_copy function
c
      if (me .eq. 0) then
         write(6,*)
         write(6,*)'> Checking copy'
         write(6,*)
         call ffflush(6)
      endif
      call ga_pgroup_sync(grp_id)
      if(me.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
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
c
c     Check the ga_scale function
c
      if (me .eq. 0) then
         write(6,92)
 92      format(/'> Checking scale ...')
         call ffflush(6)
      endif
      factor = (1d0,-1d0)
      call ga_scale(g_a, factor)
      call ga_get(g_a, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            a(i,j) = a(i,j)*factor
            if (MISMATCH(b(i,j), a(i,j)))then
               write(6,*) ' dscal ', me, i, j, a(i,j), b(i,j)
               call ffflush(6)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' scale is OK '
         write(6,*)
      endif
c
c     Check scatter&gather
c
      call ga_pgroup_sync(grp_id)
      if (me .eq. 0) then
         write(6,*) '> Checking scatter/gather (might be slow)... '
         call ffflush(6)
         if(me.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
      endif
      call ga_pgroup_sync(grp_id)
c
      crap = drand(ga_nodeid()*51 +1) !Different seed for each proc
      do j = 1, 10
       call ga_pgroup_sync(grp_id)
       itmp = iran(nproc)-1
       if(me.eq.itmp) then
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
       call ga_pgroup_sync(grp_id)
       if(me.eq.iran(ga_nnodes())-1) then
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
       call ga_pgroup_sync(grp_id)
      enddo
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' scatter is  OK'
         write(6,*)
      endif
c
c     Check ga_add
c
      if (me .eq. 0) then
         write(6,91)
 91      format(/'> Checking add ...')
         call ffflush(6)
      endif
      call ga_get(g_a, 1, n, 1, n, a, n)
      crap = drand(12345)       ! Everyone has same seed
      do j = 1, n
         do i = 1, n
            b(i,j) = cmplx(drand(0), drand(1))
            a(i,j) = (0.1d0,-.1d0)*a(i,j) + (.9d0,-.9d0)*b(i,j)
         enddo
      enddo

      if (me.eq.0) call ga_put(g_b, 1, n, 1, n, b, n)
      call ga_add((0.1d0,-.1d0), g_a, (0.9d0,-.9d0), g_b, g_b)
      call ga_get(g_b, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if (MISMATCH(b(i,j), a(i,j)))then
               write(6,*) ' add ', me, i, j, a(i,j), b(i,j)
               call ffflush(6)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' add is OK '
         write(6,*)
      endif
      call ga_pgroup_sync(grp_id)
c
c     Check the zdot function
c
      if (me .eq. 0) then
         write(6,19)
 19      format(/'> Checking zdot ...')
         call ffflush(6)
      endif
      crap = drand(12345)       ! Everyone has same seed
      sum1 = (0.0d0,0.d0)
      do j = 1, n
         do i = 1, n
            b(i,j) = cmplx(drand(0), drand(1))
            sum1 = sum1 + a(i,j)*b(i,j)
         enddo
      enddo
      if (me.eq.0) then
         call ga_put(g_b, 1, n, 1, n, b, n)
         call ga_put(g_a, 1, n, 1, n, a, n)
      endif
      call ga_pgroup_sync(grp_id)
      sum2 =  ga_zdot(g_a,g_b)
      if (MISMATCH(sum1, sum2))then
         write(6,*) ' zdot wrong ', sum1, sum2
         call ffflush(6)
         call ga_error('... exiting ',0)
      else if (me.eq.0) then
         write(6,*)
         write(6,*) ' zdot is OK '
         write(6,*)
      endif
c
c     Delete the global arrays
c
  123 continue
      status = ga_destroy(g_b)
      status = ga_destroy(g_a)
c
      return
      end
c-----------------------------------------------------------------




      subroutine check_int(grp_id)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c     
      integer n,grp_id
      parameter (n = 128)
      integer a(n,n), b(n,n)
      logical status
      integer g_a
      integer iran, i, j, loop, nloop, ilo, ihi, jlo, jhi, itmp
      integer nproc, me, int, ij, inc, dimi, dimj, ii, jj
      integer ndim, dims(2), chunk(2), p_mirror, nnodes
      integer lprocs, inode, iproc,lproc, nprocs
      double precision nwords
      parameter (nloop = 100)
      integer maxproc
      parameter (maxproc = 128)
      integer map(5,maxproc), found, np,k
      double precision crap, sum1, real
      integer buf, abs_me
      intrinsic int
      iran(i) = int(drand(0)*real(i)) + 1
c     
      nproc = ga_pgroup_nnodes(grp_id)
      me    = ga_pgroup_nodeid(grp_id)
      abs_me = ga_nodeid()
      inode = ga_cluster_nodeid()
      lprocs = ga_cluster_nprocs(inode)
      nnodes = ga_cluster_nnodes()
      iproc = mod(me,lprocs)
c     
c     a() is a local copy of what the global array should start as
c
      do j = 1, n
         do i = 1, n
            a(i,j) = i-1 + (j-1)*1000
         enddo
      enddo
c     
c     Create a global array
c
      ndim = 2
      dims(1) = n
      dims(2) = n
      chunk(1) = 0
      chunk(2) = 0
      g_a = ga_create_handle()
      call ga_set_data(g_a,ndim,dims,MT_INT)
      call ga_set_array_name(g_a,'a')
      call ga_set_pgroup(g_a,grp_id)
      if (.not.ga_allocate(g_a)) then
         write(6,*) ' ga_create failed'
         call ffflush(6)
         call ga_error('... exiting ',0)
      endif
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
      do i = 1, n
         do j = 1, n
            if (b(i,j) .ne. 0) then
               write(6,*) ' zero ', me, i, j, b(i,j)
               call ffflush(6)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' ga_zero is OK'
         write(6,*)
      endif
      call ga_pgroup_sync(grp_id)
c
c     Each node fills in disjoint sections of the array
c
      if (me .eq. 0) then
         write(6,2) 
 2       format(/'> Checking disjoint put ... ')
         call ffflush(6)
      endif
      call ga_pgroup_sync(grp_id)
c
      inc = (n-1)/20 + 1
      ij = 0
      do j = 1, n, inc
         do i = 1, n, inc
            if (mod(ij,nproc) .eq. me) then
               ilo = i
               ihi = min(i+inc, n)
               jlo = j
               jhi = min(j+inc, n)
c              write(6,4) me, ilo, ihi, jlo, jhi
c4             format(' node ',i2,' checking put ',4i4)
c              call ffflush(6)
               call ga_put(g_a, ilo, ihi, jlo, jhi, a(ilo, jlo), n)
            endif
            ij = ij + 1
         enddo
      enddo
      call ga_pgroup_sync(grp_id)
c
c     All nodes check all of a
c
      if(me.eq.0)then
        call ga_get(g_a, 1, n, 1, n, b, n)
        do i = 1, n
           do j = 1, n
              if (b(i,j) .ne. a(i,j)) then
                 write(6,*) ' put ', me, i, j, a(i,j),b(i,j)
                 call ffflush(6)
                 call ga_error('... exiting ',0)
              endif
           enddo
        enddo
      endif
      call ga_pgroup_sync(grp_id)
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' ga_put is OK'
         write(6,*)
      endif
c
c     Now check nloop random gets from each node
c
      if (me .eq. 0) then
         write(6,5) nloop
 5       format(/'> Checking random get (',i5,' calls)...')
         call ffflush(6)
      endif
      call ga_pgroup_sync(grp_id)
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
         call util_ifill(n*n, 0.0d0, b, 1)
         call ga_get(g_a, ilo, ihi, jlo, jhi, b(ilo, jlo), n)
c
         if (me .eq. 0 .and. mod(loop-1, max(1,nloop/20)).eq.0) then
            write(6,1) loop, me, ilo, ihi, jlo, jhi, nwords
 1          format(' call ',i5, ' node ',i2,' checking get ',4i4,
     $           ' total ',d9.2)
            call ffflush(6)
         endif
c
         sum1 = 0.0d0
         do j = jlo, jhi
            do i = ilo, ihi
               sum1 = sum1 + b(i,j)
               if (b(i,j) .ne. a(i,j)) then
                  write(6,*) 'error ', i, j, b(i,j), a(i,j)
                  call ffflush(6)
                  call ga_error('... exiting ',0)
               endif
            enddo
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' ga_get is OK'
         write(6,*)
      endif
c
      call ga_pgroup_sync(grp_id)
      if (me .eq. 0 .and. n.gt.7) then
         write(6,*)
         write(6,*) '> Checking ga_print_patch --- should see '
         write(6,*)' [2002 3002 4002 5002 6002]'  
         write(6,*)' [2003 3003 4003 5003 6003]'  
         write(6,*)' [2004 3004 4004 5004 6004]'  
         write(6,*)
         call ffflush(6)
      endif
      if(n.gt.5) call ga_print_patch(g_a,3,5,3,7,0)
c
      call ga_pgroup_sync(grp_id)
c
      call ga_pgroup_sync(grp_id)
      if (me .eq. 0) then
         write(6,*)
         write(6,*) '> Checking read_inc ... '
         write(6,*)
         call ffflush(6)
      endif
      call ga_pgroup_sync(grp_id)
c
      crap = drand(ga_nodeid()*51 +1) !Different seed for each proc
      inc =5 
c     every processor will be operating on somebody elses data
c
      lproc = nproc - me - 1
c
      call ga_distribution(g_a,lproc,ilo,ihi,jlo,jhi)
c
      dimi = ihi-ilo
      dimj = jhi-jlo
      write(6,111) me,ilo,ihi,jlo,jhi,dimi,dimj
  111 format(i2,'..',4i8,'.',2i8)
      call ga_pgroup_sync(grp_id)
      if(ilo .gt.0 .and. jhi .gt. 0)then
       do loop = 1,nloop 
         ii= IABS(iran(dimi))
         jj= IABS(iran(dimj))
         i=ilo + Mod(ii,dimi)
         j=jlo + Mod(jj,dimj)
c
c        write(6,*) me,'..',ilo,ihi,jlo,jhi,'.',dimi,dimj,'..',i,j
c        call ffflush(6)
         buf = ga_read_inc(g_a,i,j,inc)
         if(a(i,j).ne. buf)then   
            write(6,*)me,'READ_inc ', i,',',j,',', a(i,j),' ',buf,me 
            call ffflush(6)
         endif
         call ga_get(g_a, i,i,j,j, buf,1) 
         a(i,j) = a(i,j)+inc 
         if(a(i,j).ne.  buf)then   
            write(6,*)me,'read_INC ', i,',',j,',', a(i,j),' ',buf,me
            call ffflush(6)
         endif
       enddo
      endif
      call ga_pgroup_sync(grp_id)
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' read_inc is OK'
         write(6,*)
      endif
c
c
c
#if 000
      if (me.eq.0) then
         write(6,*)
         write(6,*) '> checking ga_fence and ga_lock'
         write(6,*)
         call ffflush(6)
      endif
c
      call ga_zero(g_a)
c
c***  use ga_read_inc and elements g_a(1:2,1) to implement a lock
c***  compute g_a(:,n) = sum (1(:)) for P processors
c     
      status = ga_create_mutexes(1)
      if (.not. status) then
         call ga_error('ga_create_mutexes failed ',0)
      endif

      if(n.lt.2)call ga_error('insufficient n to test ga_fence',n)

      call ga_lock(0)
c     call my_lock(g_a)

c          get original values  g_a(:,n)
           call ga_get(g_a, 1,n, n,n, b,n)
c          add my contribution
           do i =1,n
              b(i,1)= b(i,1)+1
           enddo
c
c          need to use fence to assure that coms complete before leaving
c          Critical Section
c
           call ga_init_fence()
                call ga_put(g_a, 1,n, n,n, b,n)
           call ga_fence()
      call ga_unlock(0)
c     call my_unlock(g_a)
c
333   if(.not.ga_destroy_mutexes())
     $   call ga_error('mutex not destroyed',0)

      call ga_pgroup_sync(grp_id)
      if (me.eq.0) then
         call ga_get(g_a, 1,n, n,n, b,n)
         do i =1,n
              if(b(i,1).ne. nproc)then
                  print *, 'mismatch',b(i,1),nproc
                  call ga_error('fence failed',i)
              endif
         enddo
         write(6,*)
         write(6,*) ' ga_fence and ga_lock are OK'
         write(6,*)
      endif 
#endif
c
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) '> checking ga_locate_region'
         write(6,*)
         call ffflush(6)
      endif

      status = ga_locate_region(g_a, 1, n, 1,n, map,np)
      found = 0
      do j=1,n
         do i=1,n
            b(i,j)=-1
         enddo
      enddo
      if(me.eq.0)call ga_put(g_a,1,n,1,n,b,n)
      call ga_pgroup_sync(grp_id)
      do k = 1, np
         if(map(5,k).eq.me)then
                if(found.eq.1) then
                    write(6,*)'double entry in map for proc ',me
                    call ffflush(6)
                endif
                do j= map(3,k), map(4,k)
                  do i= map(1,k), map(2,k)
                    b(i,j)=1*me
                  enddo
                enddo
                call ga_put(g_a, map(1,k),map(2,k),map(3,k),map(4,k),
     &                 b(map(1,k),map(3,k)),n)
                found = 1
         endif
      enddo
      call ga_pgroup_sync(grp_id)
c
      do k = 1, np
         if(map(5,k).eq.me)then
                call ga_get(g_a, map(1,k),map(2,k),map(3,k),map(4,k),
     &                 a(map(1,k),map(3,k)),n)
                do j= map(3,k), map(4,k)
                      do i= map(1,k), map(2,k)
                        if(b(i,j).ne.a(i,j)) then
                           write(6,*)
     &                     'proc ',me, 'overlap with ',a(i,j)
                           call ffflush(6)
                        endif
                      enddo
                enddo
         endif
      enddo
      call ga_pgroup_sync(grp_id)
c
      if(me.eq.0)then

        call ga_get(g_a,1,n,1,n,a,n)
        do j=1,n
         do i=1,n
            if(a(i,j).eq.-1)then
               write(6,*)'i=',i,' j=',j, ' not assigned '
               call ga_error('... exiting ',0)
            endif
         enddo
        enddo
      endif

      if (me.eq.0) then
         write(6,*)
         write(6,*) ' ga_locate_region is OK'
         write(6,*)
         call ffflush(6)
      endif
c     
c     Delete the global array
c
      status = ga_destroy(g_a)
c     
      return
      end

c---------------------------------------------------------------------

      subroutine check_flt(grp_id)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      integer n, m, grp_id
      parameter (n =10) 
      parameter (m=2*n)
      real a(n,n), b(n,n), v(m), w(m)
      integer iv(m), jv(m)
      logical status
      integer g_a, g_b
      integer iran, i, j, loop, nloop, maxloop, ilo, ihi, jlo, jhi, itmp
      integer nproc, me, int, ij, inc, ii, jj, nnodes
      double precision nwords
      parameter (maxloop = 100)
      integer maxproc
      integer ndim, dims(2), chunk(2), p_mirror
      parameter (maxproc = 128)
      double precision crap
      real x, sum1, sum2
      logical found
      integer lprocs, inode, iproc
      intrinsic int
      iran(i) = int(drand(0)*real(i)) + 1

      nproc = ga_pgroup_nnodes(grp_id)
      me    = ga_pgroup_nodeid(grp_id)
      inode = ga_cluster_nodeid()
      lprocs = ga_cluster_nprocs(inode)
      nnodes = ga_cluster_nnodes()
      iproc = mod(me,lprocs)
      nloop = Min(maxloop,n)
c     
c     a() is a local copy of what the global array should start as
c
      do j = 1, n
         do i = 1, n
            a(i,j) = i-1 + (j-1)*n
            b(i,j) = -1.
         enddo
      enddo

c     
c     Create a global array
c
      ndim = 2
      dims(1) = n
      dims(2) = n
      chunk(1) = 0
      chunk(2) = 0
      g_a = ga_create_handle()
      call ga_set_data(g_a,ndim,dims,MT_REAL)
      call ga_set_array_name(g_a,'a')
      call ga_set_pgroup(g_a,grp_id)
      if (.not.ga_allocate(g_a)) then
         write(6,*) ' ga_create failed'
         call ffflush(6)
         call ga_error('... exiting ',0)
      endif
c
c     check if handle is valid. be quiet unless error
c
      if(.not.ga_valid_handle(g_a)) call ga_error("invalid handle",g_a)
c
      call ga_distribution(g_a,me,ilo, ihi, jlo, jhi)
      call ga_pgroup_sync(grp_id)      
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
      call ga_pgroup_sync(grp_id)
      do i = 1, n
         do j = 1, n
            if (b(i,j) .ne. 0.0) then
               write(6,*) ' zero ', me, i, j, b(i,j)
               call ffflush(6)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' ga_zero is OK'
         write(6,*)
      endif
      call ga_pgroup_sync(grp_id)
c
c     Each node fills in disjoint sections of the array
c
      if (me .eq. 0) then
         write(6,2) 
 2       format(/'> Checking disjoint put ... ')   
         call ffflush(6)
      endif
      call ga_pgroup_sync(grp_id)
c
      inc = (n-1)/20 + 1
      ij = 0
      do j = 1, n, inc
         do i = 1, n, inc
            if (mod(ij,nproc) .eq. me) then
               ilo = i
               ihi = min(i+inc, n)
               jlo = j
               jhi = min(j+inc, n)
c               write(6,4) me, ilo, ihi, jlo, jhi
c 4             format(' node ',i2,' checking put ',4i4)
c               call ffflush(6)
               call ga_put(g_a, ilo, ihi, jlo, jhi, a(ilo, jlo), n)
            endif
            ij = ij + 1
         enddo
      enddo
      call ga_pgroup_sync(grp_id)
c
c     All nodes check all of a
c
      call ga_get(g_a, 1, n, 1, n, b, n)
      do i = 1, n
         do j = 1, n
            if (b(i,j) .ne. a(i,j)) then
               write(6,*) ' put ', me, i, j, a(i,j),b(i,j)
               call ffflush(6)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      call ga_pgroup_sync(grp_id)
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' ga_put is OK'
         write(6,*)
      endif
c
c     Now check nloop random gets from each node
c
      if (me .eq. 0) then
         write(6,5) nloop
 5       format(/'> Checking random get (',i5,' calls)...')
         call ffflush(6)
      endif
      call ga_pgroup_sync(grp_id)
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
         call util_rfill(n*n, 0.0, b, 1)
         call ga_get(g_a, ilo, ihi, jlo, jhi, b(ilo, jlo), n)
c
         if (me .eq. 0 .and. mod(loop-1, max(1,nloop/20)).eq.0) then
            write(6,1) loop, me, ilo, ihi, jlo, jhi, nwords
 1          format(' call ',i5, ' node ',i2,' checking get ',4i4,
     $           ' total ',d9.2)
            call ffflush(6)
         endif
c
         sum1 = 0.0d0
         do j = jlo, jhi
            do i = ilo, ihi
               sum1 = sum1 + b(i,j)
               if (b(i,j) .ne. a(i,j)) then
                  write(6,*) 'error ', i, j, b(i,j), a(i,j)
                  call ffflush(6)
                  call ga_error('... exiting ',0)
               endif
            enddo
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' ga_get is OK'
         write(6,*)
      endif
      call ga_pgroup_sync(grp_id)
c
c     Each node accumulates into disjoint sections of the array
c
      if (me .eq. 0) then
         write(6,9)
 9       format(/'> Checking accumulate ... ')
         call ffflush(6)
      endif
      call ga_pgroup_sync(grp_id)

c
      crap = drand(12345)       ! Same seed for each process
      do j = 1, n
         do i = 1, n
c           b(i,j) = real(drand(0))
            b(i,j) = i+j
         enddo
      enddo
c
      inc = (n-1)/20 + 1
      ij = 0
      do j = 1, n, inc
         do i = 1, n, inc
c           x = real(drand(0))
            x = 10.
            ilo = i
            ihi = min(i+inc-1, n)
            if(ihi.eq.n-1)ihi=n
c           ihi = min(i+inc, n)
            jlo = j
            jhi = min(j+inc-1, n)
            if(jhi.eq.n-1)jhi=n
c           jhi = min(j+inc-1, n)
*               call ffflush(6)
            if (mod(ij,nproc) .eq. me) then
c               print *, me, 'checking accumulate ',ilo,ihi,jlo,jhi,x
* 11            format(' node ',i2,' checking accumulate ',4i4)
*               call ffflush(6)
               call ga_acc(g_a, ilo, ihi, jlo, jhi, b(ilo, jlo), n, x)
            endif
            ij = ij + 1
c
c           Each process applies all updates to its local copy
c
            do jj = jlo, jhi
               do ii = ilo, ihi
                  a(ii,jj) = a(ii,jj) + x * b(ii,jj)
               enddo
            enddo
         enddo
      enddo
      call ga_pgroup_sync(grp_id)                       
c
c     All nodes check all of a
      call ga_get(g_a, 1, n, 1, n, b, n)
c
      do j = 1, n
         do i = 1, n
            if(MISMATCHF(b(i,j),a(i,j)))then
               write(6,*) ' acc ', me, i, j, a(i,j), b(i,j)
               call ffflush(6)
               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' disjoint ga_acc is OK'
         write(6,*)
      endif
c
c     overlapping accumulate
      call ga_pgroup_sync(grp_id)
      g_b = ga_create_handle()
      call ga_set_data(g_b,ndim,dims,MT_REAL)
      call ga_set_array_name(g_b,'b')
      call ga_set_pgroup(g_b,grp_id)
      if (.not.ga_allocate(g_b)) then
         call ga_error('ga_create failed for second array ',0)
      endif
c
      call ga_zero(g_b)
      call ga_acc(g_b, n/2, n/2, n/2, n/2, 1.0, 1, 1.0)
      call ga_pgroup_sync(grp_id)
      if (me.eq.0) then
         call ga_get(g_b, n/2, n/2, n/2, n/2, b(1,1), 1)
         x = abs(b(1,1) -1*nproc)
         if(x.gt. 1e-10)then
            write(6,*)'val=',b(1,1),' expected=',nproc, x
            call ga_error('overlapping accumulate failed',0)
         endif
         write(6,*)
         write(6,*) ' overlapping ga_acc is OK'
         write(6,*)
      endif
c
c     Check the ga_add function
c
      if (me .eq. 0) then
         write(6,91)
 91      format(/'> Checking add ...')
         call ffflush(6)
      endif                           
c
      crap = drand(12345)       ! Everyone has same seed
      do j = 1, n
         do i = 1, n
            b(i,j) = real(drand(0)*real(i)) + 1
            a(i,j) = 0.1*a(i,j) + 0.9*b(i,j)
         enddo
      enddo
      if (me.eq.0) call ga_put(g_b, 1, n, 1, n, b, n)
      call ga_add(0.1, g_a, 0.9, g_b, g_b)
      call ga_get(g_b, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if(MISMATCHF(b(i,j), a(i,j)))then
               write(6,*) ' add ', me, i, j, a(i,j), b(i,j)
               call ffflush(6)
CCC               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' add is OK '
         write(6,*)
      endif
      call ga_pgroup_sync(grp_id)
c
c     Check the sdot function
c
      if (me .eq. 0) then
         write(6,19)
 19      format(/'> Checking sdot ...')
         call ffflush(6)
      endif
      crap = drand(12345)       ! Everyone has same seed
      sum1 = 0.0
      do j = 1, n
         do i = 1, n
            b(i,j) = drand(0)
            sum1 = sum1 + a(i,j)*b(i,j)
         enddo
      enddo
      if (me.eq.0) then
         call ga_put(g_b, 1, n, 1, n, b, n)
         call ga_put(g_a, 1, n, 1, n, a, n)
      endif
      call ga_pgroup_sync(grp_id)
      sum2 = ga_sdot(g_a,g_b)
      if(MISMATCHF(sum1, sum2))then
         write(6,*) ' fdot wrong ', sum1, sum2
         call ffflush(6)
            call ga_error('... exiting ',0)
      else if (me.eq.0) then
         write(6,*)
         write(6,*) ' sdot is OK '
         write(6,*)
      endif
c
c     Check the ga_scale function
c
      if (me .eq. 0) then
         write(6,92)
 92      format(/'> Checking scale ...')
         call ffflush(6)
      endif
      call ga_scale(g_a, 0.123)
      call ga_get(g_a, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            a(i,j) = a(i,j)*0.123
            if (MISMATCHF(b(i,j), a(i,j)))then
               write(6,*) ' dscal ', me, i, j, a(i,j), b(i,j)
               call ffflush(6)
CCC               call ga_error('... exiting ',0)
            endif
         enddo
      enddo
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' scale is OK '
         write(6,*)
      endif
c                                
c     Check the ga_copy function
c
      if (me .eq. 0) then
         write(6,*)
         write(6,*)'> Checking copy'
         write(6,*)
         call ffflush(6)
      endif
      if(me.eq.0) call ga_put(g_a, 1, n, 1, n, a, n)
      call ga_copy(g_a, g_b)
      call ga_get(g_b, 1, n, 1, n, b, n)
      do j = 1, n
         do i = 1, n
            if (b(i,j) .ne. a(i,j)) then
               write(6,*) ' copy ', me, i, j, a(i,j), b(i,j)
               call ffflush(6)
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
      call ga_pgroup_sync(grp_id)
      if (me .eq. 0) then
         write(6,*) '> Checking scatter/gather (might be slow)... '
         call ffflush(6)
      endif
      call ga_pgroup_sync(grp_id)
c
      crap = drand(ga_nodeid()*51 +1) !Different seed for each proc
      do j = 1, 10
       call ga_pgroup_sync(grp_id)
       itmp = iran(nproc)-1
       if(me.eq.itmp) then
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
       call ga_pgroup_sync(grp_id)
       if(me.eq.iran(ga_nnodes())-1) then
         do loop = 1,m
           ilo = iran(n)
           jlo = iran(n)
           iv(loop) = ilo
           jv(loop) = jlo
c          v(loop) = DSIN(a(ilo,jlo)+b(ilo,jlo))
           v(loop) = 1.0 *(ilo+jlo)
         enddo
         call ga_scatter(g_a, v, iv, jv, m)
         do loop = 1,m
           ilo= iv(loop)
           jlo= jv(loop)
           call ga_get(g_a,ilo,ilo,jlo,jlo,w(loop),1)
c          if(v(loop)  .ne. w(loop))then
           if(w(loop)  .ne. 1.0 *(ilo+jlo) )then
             write(6,*)me,' scatter ', ilo,',',jlo,',',w(loop)
     &             ,' ', 1.0 *(ilo+jlo)
             call ffflush(6)
             call ga_error('... exiting ',0)
           endif
         enddo
       endif
       call ga_pgroup_sync(grp_id)
      enddo
c
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' scatter is  OK'
         write(6,*)
      endif
c
      call ga_pgroup_sync(grp_id)
c
c     scatter-acc available in GA ver. 3.0
#ifdef GA3                                         
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' checking scatter-accumulate'
         write(6,*)
      endif
c
      crap = drand(1234)
      call ga_zero(g_a)
c
      do j = 1, n
         do i = 1, n
            b(i,j) =0.
         enddo
      enddo
c
      x = .1d0
      ii =n
      do jj = 1,1
          call ga_pgroup_sync(grp_id)
        do loop = 1, ii
 
c          generate unique i,j pairs
10         continue
              i = iran(n)
              j=iran(n)
           if (found(i,j, iv, jv, loop-1) ) goto 10

           iv(loop) = i
           jv(loop) = j
           v(loop) = 1.0 *(i+j)
           b(i,j) = b(i,j) + nproc*x*v(loop) ! update local ref. array
        enddo
 
        call ga_scatter_acc(g_a,v,iv,jv, ii,x)
 
c
        call ga_pgroup_sync(grp_id)
c
c       check the result
c
        call ga_get(g_a, 1, n, 1,n, a, n)
 
        do loop = 1,ii
          i = iv(loop)
          j = jv(loop)
          if(MISMATCHF(a(i,j),b(i,j)))then
             print *,'Error',i,j,loop,a(i,j),'expected=',b(i,j)
*            if(me.eq.0)then
*              do ii=1,loop
*                   print *,'element',ii, iv(ii),jv(ii)
*              enddo
*            endif
             call ga_error('scatter-acc error in trial ',jj)
          endif
        enddo                                 
          call ga_pgroup_sync(grp_id)
 
      enddo
 
      call ga_pgroup_sync(grp_id)
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' scatter-accumulate is  OK'
         write(6,*)
      endif
#endif                      
c
c     Delete the global array
c
      status = ga_destroy(g_a)
c     
      return
      end
c_____________________________________________________________

      subroutine check_wrappers(grp_id) 
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
      double precision sum
      integer isum, ibuf(2)
      integer me, nproc, grp_id, root
      real fsum
c
      nproc = ga_pgroup_nnodes(grp_id)
      me = ga_pgroup_nodeid(grp_id)
c
      call ga_pgroup_sync(grp_id)
      if (me .eq. 0) then
         write(6,*)
         write(6,*)'> Checking ga_igop'
         write(6,*)
         call ffflush(6)
      endif
      ibuf(1) = 1
      ibuf(2) =  me
      call ga_pgroup_igop(grp_id,10000, ibuf, 2, '+')
      if(ibuf(1).ne.nproc)then
        call ga_error('ga_igop error',isum)
      endif
      if(ibuf(2).ne.((nproc-1)*nproc/2))then
        call ga_error('ga_igop error -2',isum)
      endif
      call ga_pgroup_sync(grp_id)
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' ga_igop is  OK'
         write(6,*)
      endif
      call ga_pgroup_sync(grp_id)
c
      call ga_pgroup_sync(grp_id)
      if (me .eq. 0) then
         write(6,*)
         write(6,*)'> Checking ga_dgop'
         write(6,*)
         call ffflush(6)
      endif
      sum = 1d0 * me
      call ga_pgroup_dgop(grp_id,10000, sum, 1, '+')
      if(Int(sum).ne.((nproc-1)*nproc/2))then
        call ga_error('ga_dgop error',Int(sum))
      endif
      call ga_pgroup_sync(grp_id)
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' ga_dgop is  OK'
         write(6,*)
      endif
c
      call ga_pgroup_sync(grp_id)
      if (me .eq. 0) then
         write(6,*)
         write(6,*)'> Checking ga_sgop'
         write(6,*)
         call ffflush(6)
      endif
      fsum = 1.0 * me
      call ga_pgroup_sgop(grp_id,10000, fsum, 1, '+')
      if(Int(sum).ne.((nproc-1)*nproc/2))then
        call ga_error('ga_fgop error',Int(sum))
      endif
      call ga_pgroup_sync(grp_id)
      if (me.eq.0) then
         write(6,*)
         write(6,*) ' ga_sgop is  OK'
         write(6,*)
      endif
c                                   
      call ga_pgroup_sync(grp_id)
      if (me .eq. 0) then
         write(6,*)
         write(6,*)'> Checking ga_brdcst'
         write(6,*)
         call ffflush(6)
      endif
      if(me.eq.nproc-1)then
        ibuf(1) = me
        ibuf(2) = nproc
      endif
      root = 7
      call ga_pgroup_brdcst(grp_id,1000,ibuf,util_mitob(2),nproc-1)
      if(ibuf(1).ne.nproc-1)call ga_error('ibuf(1) error',ibuf(1))
      if(ibuf(2).ne.nproc)call ga_error('ibuf(2) error',ibuf(2))
      call ga_pgroup_sync(grp_id)
      if (me .eq. 0) then
         write(6,*)
         write(6,*)'> ga_brdcst is OK '
         write(6,*)
         call ffflush(6)
      endif
      call ga_pgroup_sync(grp_id)
      return
      end


      subroutine check_mem(mem_size)
      implicit none
      integer mem_size
#include "mafdecls.fh"
#include "global.fh"
#include "testutil.fh"
c
      integer n,nmax,left,need, me,procs,g_a, g_b 
      logical status
c
      if(.not. ga_memory_limited())return
      me = ga_nodeid()
      procs = ga_nnodes()
      nmax = int(dsqrt(dble(mem_size/util_mitob(1))))
      left =  mem_size/procs - ga_inquire_memory()
      n = nmax/2
      need = util_mdtob(n*n)/procs
c
      if(me.eq.0)then
        write(6,*)' '
        if(ga_uses_ma())then
           write(6,*)' CHECKING GA MEMORY RESTRICTIONS (MA used)'
        else
           write(6,*)' CHECKING GA MEMORY RESTRICTIONS (MA not used)'
        endif
        write(6,*)' '
        write(6,*)' '
        call print_mem_info(n,left, need, mem_size/procs)
      endif
c
      status = ga_create(MT_DBL, n, n, 'a', 0, 0, g_a)
c
      if(me.eq.0) then
        if(status) then
              write(6,*) '  success'
        else
              write(6,*) '  failure'
        endif
        call ffflush(6)
      endif
c
      n = nmax
      left =  mem_size/procs - ga_inquire_memory()
      need = util_mdtob(n*n)/procs
      if(me.eq.0)then
        call print_mem_info(n,left, need, mem_size/procs)
      endif
c
      status = ga_create(MT_DBL, n, n, 'b', 0, 0, g_b)
c
      if(me.eq.0) then
        if(status) then
              write(6,*) '  success'
        else
              write(6,*) '  failure'
        endif
        write(6,*)' '
        write(6,*)' '
        call ffflush(6)
      endif
      status = ga_destroy(g_a)
      end

     

      subroutine print_mem_info(n,left, need, total)
      implicit none
      integer n,left, need, total
c
      write(6,*)' '
      if(left - need .ge. 0) then
        write(6,1)n,n 
1       format('> Creating array ',i4,' by ',i4,' -- should succeed')
      else
        write(6,2)n,n 
2       format('> Creating array ',i4,' by ',i4,' -- SHOULD FAIL')
      endif
      write(6,3) need, left, total
3     format('  (need ',i7,' and  ',i7,' out of ',i7,' bytes are left)')
      write(6,*)' '
      call ffflush(6)
c
      end


      
      subroutine my_lock(g_b)
      implicit none
#include "global.fh"
      integer g_b, val, flag, i
      logical first_time
      double precision  dummy 
      common /lock/ val 
      common /dum/ dummy
      data first_time /.true./
c
c     this awkward initialization is to avoid a weird problem 
C     with block data on SUN
      if(first_time)then
        first_time = .false.
        dummy = .0
      endif
c
      val = ga_read_inc(g_b,1,1, 1)
10    call ga_get(g_b, 2,2,1,1, flag, 1)
           if(flag.eq.val) return 
c
c          to reduce memory stress,  wait a while before retrying
           do i = 1, 100
              dummy = dummy + .1 
           enddo 
      goto 10
      end


      subroutine my_unlock(g_b)
      implicit none
#include "global.fh"
      integer g_b, val
      common /lock/ val
c
      call ga_put(g_b, 2,2,1,1, val+1, 1) 
      end


      logical function found(i,j, iv, jv, n)
      integer n
      integer i,j, iv(n), jv(n)
      integer loop

      found = .false.
      do loop = 1, n
               if(i .eq. iv(loop) .and. j .eq.jv(loop))then
                 found = .true.
                 goto 99
               endif
            enddo
99    continue
      return
      end 
      
            
      subroutine proc_remap()
      implicit none
#include "global.fh"
      integer proc(100),nproc,i
      nproc = ga_nnodes()
      if(nproc.gt.100)
     $   call ga_error("remap requires<=100 processes",nproc)
      do i = 1, nproc
         proc(i) = nproc-i
      enddo
c     call ga_register_proclist(proc,nproc)
      end

      
      subroutine util_rfill(n,val,a,ia)
      implicit none
      real  a(*), val
      integer n, ia, i
c
c     initialise real array to scalar value
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

      
      subroutine util_dfill(n,val,a,ia)
      implicit none
      double precision  a(*), val
      integer n, ia, i
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

      subroutine util_ifill(n,val,a,ia)
      implicit none
      integer n, ia, i, a(*),val
c
c     initialise integer array to scalar value
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


