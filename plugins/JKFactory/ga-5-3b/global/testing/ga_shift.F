#if HAVE_CONFIG_H
#   include "config.fh"
#endif

# define THRESH 1d-13
# define THRESHF 1e-5
#define MISMATCH(x,y) abs(x-y)/max(1d0,abs(x)).gt.THRESH
#define MISMATCHF(x,y) abs(x-y)/max(1.0,abs(x)).gt.THRESHF


      program main
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer DIM, MAXPROC, NITER
      parameter (DIM = 2, MAXPROC = 1000000, NITER=10)
      integer MAXSIZE
      parameter (MAXSIZE = 2**NITER)
      integer heap, stack, fudge, ma_heap
      integer me, nproc, nshift
      integer iter, lsize, ndim, dims(2), lo(2), hi(2), chunk(2), ld
      integer g_a
      integer i, j, ioff
      double precision tbeg, t_init, t_term, t_ma_init, t_msg_init
      double precision t_create, t_destroy, t_put, t_get, t_sync, bw
      double precision a(MAXSIZE,MAXSIZE)
      logical status
c     
c***  Intitialize a message passing library
c
#ifdef MPI
      integer ierr
c     tbeg = ga_wtime()
      call mpi_init(ierr)
#else
c     tbeg = ga_wtime()
      call pbeginf
#endif
c     t_msg_init = ga_wtime() - tbeg
c
c***  Initialize GA
c
c     There are 2 choices: ga_initialize or ga_initialize_ltd.
c     In the first case, there is no explicit limit on memory usage.
c     In the second, user can set limit (per processor) in bytes.
c
      tbeg = ga_wtime()
      call ga_initialize()
      t_init = ga_wtime() - tbeg

      nproc = ga_nnodes()
      me = ga_nodeid()
      nshift = nproc/2
c
c     call ga_dgop(1,t_msg_init,1,'+')
c     t_msg_init = t_msg_init/dble(nproc)
c     if (me.eq.0) then
c       write(6,'(a,f16.6)') 'Time spent in runtime initialization: ',
c    +                       t_msg_init
c     endif
      call ga_dgop(2,t_init,1,'+')
      t_init = t_init/dble(nproc)
      if (me.eq.0) then
        write(6,'(a,f16.6)') 'Time spent in GA initialization:      ',
     +                       t_init
      endif
c     we can also use GA_set_memory_limit BEFORE first ga_create call
c
      heap = 2000*2000*4
      fudge = 100
      stack = 2000*2000
      ma_heap = heap + fudge 
      call GA_set_memory_limit(8*ma_heap)
c
      if(ga_nodeid().eq.0)then
         write(6,'(a)') ' '
         write(6,'(a)') ' GA initialized '
         write(6,'(a)') ' '
         call ffflush(6)
      endif
c
c***  Initialize the MA package
c     MA must be initialized before any global array is allocated
c
      tbeg = ga_wtime()
      status = ma_init(MT_DCPL, stack, ma_heap)
      t_ma_init = ga_wtime() - tbeg
      if (.not. status) call ga_error('ma_init failed',-1) 
      call ga_dgop(3,t_ma_init,1,'+')
      t_ma_init = t_ma_init/dble(nproc)
      if (me.eq.0) then
        write(6,'(a,f16.6)') 'Time spent in MA initialization:      ',
     +                       t_ma_init
      endif
c
      if(me.eq.0)then
        write(6,'(a,i8,a)') 'Using ', nproc, ' process(es)'
        call ffflush(6)
      endif
c
c*** begin looping over tests
c
      lsize = 1
      do iter = 1, NITER
        lsize = 2*lsize
c
c*** create GA that is lsize X lsize*nproc in dimension
c
        if (me.eq.0) then
          write(6,'(a)') ' '
          write(6,'(a,i8)') ' Testing block size of ',lsize*lsize
          write(6,'(a)') ' '
        endif
        ndim = 2
        dims(1) = lsize
        dims(2) = nproc*lsize
        chunk(1) = lsize
        chunk(2) = -1
        call ga_sync
        tbeg = ga_wtime()
        g_a = ga_create_handle()
        call ga_set_data(g_a,ndim,dims,MT_DBL)
        call ga_set_chunk(g_a,chunk)
        status = ga_allocate(g_a)
        if (.not.status) then
          call ga_error('ga_allocate failed for size: ',lsize)
        endif
        t_create = ga_wtime() - tbeg
        call ga_dgop(iter,t_create,1,'+')
        t_create = t_create/dble(nproc)
        if (me.eq.0) then
          write(6,'(a,f16.6)') 'Time spent creating GA:               ',
     +                       t_create
        endif
c
c*** Fill local buffer with values
c
        ioff = mod(me+nshift,nproc)
        ioff = ioff*lsize**2
        do j = 1, lsize
          do i = 1, lsize
            a(i,j) = dble((j-1)*lsize+i+ioff)
          end do
        end do
        ld = MAXSIZE
        lo(1) = 1
        hi(1) = lsize
        lo(2) = mod(me+nshift,nproc)
        lo(2) = lo(2)*lsize + 1
        hi(2) = lo(2) - 1 + lsize
        call ga_zero(g_a)
        tbeg = ga_wtime()
        call nga_put(g_a,lo,hi,a,ld)
        t_put = ga_wtime() - tbeg
        call ga_dgop(iter+1,t_put,1,'+')
        t_put = t_put/dble(nproc)
        if (me.eq.0) then
          write(6,'(a,f16.6)') 'Time spent in put:                    ',
     +                       t_put
        endif
        bw = dble(8*lsize*lsize)/t_put
        bw = bw/1000000.0d00
        if (me.eq.0) then
          write(6,'(a,f16.6)') 'Bandwidth for put (MB/s):             ',
     +                       bw
        endif
        tbeg = ga_wtime()
        call ga_sync
        t_sync = ga_wtime() - tbeg
        call ga_dgop(iter+2,t_sync,1,'+')
        t_sync = t_sync/dble(nproc)
        if (me.eq.0) then
          write(6,'(a,f16.6)') 'Time spent in sync:                   ',
     +                       t_sync
        endif
        bw = dble(8*lsize*lsize)/(t_put+t_sync)
        bw = bw/1000000.0d00
        if (me.eq.0) then
          write(6,'(a,f16.6)') 'Bandwidth for put with sync (MB/s):   ',
     +                       bw
        endif
        do j = 1, lsize
          do i = 1, lsize
            a(i,j) = 0.0d00
          end do
        end do
        call ga_sync
        tbeg = ga_wtime()
        call nga_get(g_a,lo,hi,a,ld)
        t_get = ga_wtime() - tbeg
        ioff = mod(me+nshift,nproc)
        ioff = ioff*lsize**2
        do j = 1, lsize
          do i = 1, lsize
            if (a(i,j).ne.dble((j-1)*lsize+i+ioff)) then
              write(6,'(i4,a,i8,a,i8,a)') me,' mismatch for element (',
     +                                    i,',',j,')' 
            endif
          end do
        end do
        call ga_dgop(iter+3,t_get,1,'+')
        t_get = t_get/dble(nproc)
        if (me.eq.0) then
          write(6,'(a,f16.6)') 'Time spent in get:                    ',
     +                       t_get
        endif
        bw = dble(8*lsize*lsize)/t_get
        bw = bw/1000000.0d00
        if (me.eq.0) then
          write(6,'(a,f16.6)') 'Bandwidth for get (MB/s):             ',
     +                       bw
        endif
        call ga_sync
        tbeg = ga_wtime()
        status = ga_destroy(g_a)
        t_destroy = ga_wtime() - tbeg
        call ga_dgop(iter+4,t_destroy,1,'+')
        t_destroy = t_destroy/dble(nproc)
        if (me.eq.0) then
          write(6,'(a,f16.6)') 'Time spent destroying GA:             ',
     +                       t_destroy
        endif
      end do
c
c***  Tidy up the GA package
c
      call ga_terminate()
c
c***  Tidy up after message-passing library
c
#ifdef MPI
      call mpi_finalize(ierr)
#else
      call pend()
#endif
c
      stop
      end
