#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c #define USE_SUMMA 1
#ifdef  USE_SUMMA
#  define GA_CREATE  ga_create2
#  define GA_DGEMM   ga_summa___
#else
#  define GA_CREATE  ga_create
#  define GA_DGEMM   ga_dgemm
#endif

#define BASE_NAME  'matrixA.da'
#define BASE_NAME1 'matrixB.da'
#define BASE_NAME2 'matrixC.da'
#ifdef  HPIODIR
#  define FNAME   HPIODIR//BASE_NAME
#  define FNAME1  HPIODIR//BASE_NAME1
#  define FNAME2  HPIODIR//BASE_NAME2
#else
#  define FNAME   BASE_NAME
#  define FNAME1  BASE_NAME1
#  define FNAME2  BASE_NAME2
#endif

      program test
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "dra.fh"
      integer max_arrays
      integer status
      double precision max_sz, max_disk, max_mem
      data max_arrays, max_sz, max_disk, max_mem /4,1d8,1d10, 1d6/

#include "mp3.fh"
      if (.not. ma_init(MT_DBL, 800000, 800000))
     $     call ga_error('ma_init failed', 1)
      call ga_initialize()
      if(dra_init(max_arrays, max_sz, max_disk, max_mem).ne.0)
     $           call ga_error('dra_init failed: ',0)
c
      call testit()
c
      status = dra_terminate()
      call ga_terminate()
      call MP_FINALIZE()
      end


      subroutine testit()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "dra.fh"
#include "mp3def.fh"
c
      integer n, nrhs, nb
      integer g_temp
      integer d_a, d_b, d_c
      integer i, j, req
      integer nproc, me
c     
      parameter (n = 500)
      parameter (nrhs = n/2)
      parameter (nb = n/5)
c
      double precision temp(n), mbytes
      double precision tt0, tt1
      logical status
      integer util_mdtob
      external util_mdtob
c     ...................................................
c     
      nproc = ga_nnodes()
      me    = ga_nodeid()
c     
c***  Create global arrays
c
      if (.not. GA_CREATE(MT_DBL, n, n, 'temp', 1, 1, g_temp))
     $     call ga_error(' GA_CREATE failed: temp ',1)
c     
c***  creating disk arrays
c
      if (me.eq.0) print*, 'Creating Disk Array A:',n,' x',n
      if(dra_create(MT_DBL, n, n, 'A', 
     &    FNAME, 
     &    DRA_RW, nb, n, d_a).ne.0)
     $   CALL ga_error('dra_create failed: ',0)
      if (me.eq.0) print*, 'Creating Disk Array B:',n,' x',n
      if(dra_create(MT_DBL, n, nrhs, 'B',
     &    FNAME1, 
     &    DRA_RW, n, nb, d_b).ne.0)
     $   CALL ga_error('dra_create failed: ',0)
      if (me.eq.0) print*, 'Creating Disk Array C:',n,' x',n
      if(dra_create(MT_DBL, n, nrhs, 'C',
     &    FNAME2, 
     &    DRA_RW, nb, nb, d_c).ne.0)
     $   CALL ga_error('dra_create failed: ',0)
c
c*************************  Fill in arrays A & B  ************************
c
      do j = 1+me, n, nproc ! static load balancing
         do i = 1, n
            temp(i) = 1d0*(i +j)
         enddo
         call ga_put(g_temp, 1,n, j,j, temp,n)
      enddo
*     call ga_print(g_temp, 1)
      call ga_sync()
c
c...  writing to disk arrays 
c
      if(me.eq.0) print *, 'writing global array  to disk'
      if (me.eq.0)call ffflush(6)
      tt0 = MP_TIMER()
      if(dra_write(g_temp,d_a,req).ne.0)
     $   CALL ga_error('write failed:',0)
      if(dra_wait(req).ne.0) call ga_error('dra_wait failed: ' ,req)
      if(dra_close(d_a).ne.0)call ga_error('dra_close failed: ',d_a)
      tt1 = MP_TIMER() -tt0
      mbytes = 1e-6*util_mdtob(n*n)
      if (me.eq.0)then
          write(6,100)mbytes,tt1,mbytes/tt1
      endif
c
c     Now, array B
c
      status =  ga_destroy(g_temp)
      if (.not. ga_create(MT_DBL, n, nrhs, 'temp', 1, 1, g_temp))
     $     call ga_error(' GA_CREATE failed: temp ',1)
c
      do j = 1+me, nrhs, nproc ! static load balancing
         do i = 1, n
            temp(i) = dble(j)
         enddo
         call ga_put(g_temp, 1,n, j,j, temp,n)
      enddo
*     call ga_print(g_temp, 1)
      call ga_sync()
c
c...  writing to disk arrays
c
      if(me.eq.0) print *, 'writing global array  to disk'
      if (me.eq.0)call ffflush(6)
      tt0 = MP_TIMER()
      if(dra_write(g_temp,d_b,req).ne.0)
     $   CALL ga_error('write failed:',0)
      if(dra_wait(req).ne.0) call ga_error('dra_wait failed: ' ,req)
      if(dra_close(d_b).ne.0)call ga_error('dra_close failed: ',d_b)
      tt1 = MP_TIMER() -tt0
      mbytes = 1e-6*util_mdtob(n*n)
      if (me.eq.0)then
          write(6,100)mbytes,tt1,mbytes/tt1
      endif
c
c**********************************************************************
c
c***  opening disk arrays
      if(dra_open(FNAME,DRA_R, d_a).ne.0)
     &            call ga_error('dra_open failed',0)
      if(dra_open(FNAME1,DRA_R, d_b).ne.0)
     &            call ga_error('dra_open failed',0)
c
      if (me .eq. 0) then
         print *,' '
         print *,' Checking out-of-core matrix multiplication'
         print *,' '
         call ffflush(6)
      endif
c
      call ga_sync()
      tt0 = MP_TIMER()
      call chunk_mxm(n, nrhs, nb, d_a, d_b, d_c)
      tt1 = MP_TIMER() -tt0
      if(me.eq.0)print *,'time=',tt1,' MFLOPS=', 2d-6*n*n*nrhs/tt1 
c
c     this is a hack to avoid specifying -lglobal twice 
      if(ga_nnodes().le.0) call ga_dgemm()
      call ga_sync()
c
      if(me.eq.0) print *, 'reading global array from disk',g_temp
      if (me.eq.0)call ffflush(6)
      if(dra_read(g_temp,d_c,req).ne.0)call ga_error('read failed:',0)
      if(dra_wait(req).ne.0) call ga_error('dra_wait failed: ' ,req)
c     if(dra_close(d_a).ne.0)call ga_error('dra_close failed: ',d_a)
c
c***  clean-up
c
      if(dra_delete(d_a).ne.0)
     &            call ga_error('dra_delete failed',0)
      if(dra_delete(d_b).ne.0)
     &            call ga_error('dra_delete failed',0)
      if(dra_delete(d_c).ne.0)
     &            call ga_error('dra_delete failed',0)
c
c
*     call ga_print(g_temp, 1)
      do j = 1+me, nrhs, nproc ! check result
         call ga_get(g_temp, 1,n, j,j, temp,n)
         do i = 1, n
            if(ABS(temp(i)-0.5d0*n*j*(n+1+2*i)).gt. 1d-8) then
              print *,i,j, temp(i)
              call ga_error('error :',i)
            endif
         enddo
      enddo
      call ga_sync()
c
c***  delete temp GA array
      status = ga_destroy(g_temp)
      if(me.eq.0) print *, 'all done'
100   format(g11.2,' MB  time=',g11.2,' rate=',g11.3,'MB/s')
      end


c     given current (i,j) compute (ni, nj) - next pair of loop indices
c     i - outermost loop, j- innermost loop
c     return .false. if (i,j) = (ni,nj), otherwise .true.
c     
      logical function Next2D(i, imin, imax, j, jmin, jmax, ni, nj)
      implicit none
      integer i, imin, imax, j, jmin, jmax, ni, nj
c
      if(i .lt. imin .or. j.lt.jmin)then
c          to handle initial out-of-range indices
           ni = imin
           nj = jmin
      else
           ni = i
           nj = j+1
      endif
      if(nj.gt.jmax)then
         nj = jmin
         ni = ni + 1
      endif
      if(ni.gt.imax)ni=imin            ! restart from the beginning 
      Next2D = ni.ne.i .or. nj.ne.j    ! check if we have a new pair
      end
        

c     given current (i,j) block, issues asynchronous read for next block of matrix
c      
      integer function Prefetch_Next(g_a, d_a, i, j, nb, n, m)
      implicit none
#include "global.fh"
#include "dra.fh"
      integer g_a, d_a, i, j, nb, n,m
      integer jlo, jhi, ni, nj, chunks_i, chunks_j, nc
      integer LO, HI
      logical next2d
      external next2d
      LO(i)  = (i-1)*nb +1
      HI(i)  = MIN(LO(i) + nb -1, m)
c
      chunks_i = 1  ! chunking j only - reads full columns
      chunks_j = (m+nb-1)/nb
      if( next2d(i, 1, chunks_i, j, 1, chunks_j, ni, nj)) then
        jlo = LO(nj)
        jhi = HI(nj)
        nc  = jhi -jlo+1
        if(dra_read_section(.false., g_a, 1, n, 1, nc,
     &             d_a, 1, n, jlo, jhi, prefetch_next).ne.0)
     &                       call ga_error('dra_read failed:',0)

      else
        Prefetch_Next = DRA_REQ_INVALID ! there is nothing new to read
      endif
      end
   

c     Out-of-core matrix multiplication C = A*B, with prefetching.
c
c     A is n x n, B and C are n x m. All arrays reside on disk.
c     nb is a block dimension for the in-core matrices, A_incore: nb x n;
c     B_incore: n x nb; C_incore: nb x nb.
c
      subroutine chunk_mxm(n, m, nb, d_a, d_b, d_c)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "dra.fh"
#include "mp3def.fh"

      integer n, m, nb, d_a, d_b, d_c
      integer i, j, chunks_i, chunks_j
      integer ilo, ihi, jlo, jhi, nr, nc
      integer g_a_chunk, g_b_chunk(2), g_c_chunk
      logical status
      integer req, req_b, req_c, Prefetch_Next
      integer nxt_idx, cur_idx
      double precision tot, t0
      logical ga_create2
      integer LO, HI, NEXT
      LO(i,n)  = (i-1)*nb +1
      HI(i,n)  = MIN(LO(i,n) + nb -1, n)
      NEXT(i)  = Mod(i+2,2)+1        ! generates series: 1, 2, 1, 2, 1, ...
c
      chunks_i = (n+nb-1)/nb
      chunks_j = (m+nb-1)/nb
c
c***  create temporary "chunk" GA arrays
c
      if (.not. GA_CREATE(MT_DBL, nb, n, 'a', 1, 1, g_a_chunk))
     $     call ga_error('chunk_mxm: GA_CREATE failed: g_a ',0)

      if (.not. GA_CREATE(MT_DBL, n, nb, 'b(1)', 1, 1, g_b_chunk(1)))
     $     call ga_error('chunk_mxm: GA_CREATE failed: g_b(1) ',0)
      if (.not. GA_CREATE(MT_DBL, n, nb, 'b(2)', 1, 1, g_b_chunk(2)))
     $     call ga_error('chunk_mxm: GA_CREATE failed: g_b(2) ',0)

      if (.not. GA_CREATE(MT_DBL, nb, nb, 'c', 1, 1, g_c_chunk))
     $     call ga_error('chunk_mxm: GA_CREATE failed: g_c ',0)
c
c
      cur_idx = 1
      req_b = Prefetch_Next(g_b_chunk(cur_idx), d_b, 0, 0, nb, n, m)
      if(dra_wait(req_b).ne.0) call ga_error('dra_wait failed:',req)
 
      do i = 1, chunks_i
c
         ilo = LO(i,n)
         ihi = HI(i,n) 
         nr  = ihi-ilo+1
         if(dra_read_section(.false., g_a_chunk, 1, nr, 1, n,
     &                               d_a, ilo, ihi, 1, n, req).ne.0)
     &                          call ga_error('dra_read failed:',0)
         if(dra_wait(req).ne.0)  call ga_error('dra_wait failed:',req)
c
*        call ga_print_patch(g_a_chunk, 1, nr, 1, n, 1)
         do j = 1, chunks_j
c
            jlo = LO(j,m)
            jhi = HI(j,m) 
            nc  = jhi-jlo+1

            if(dra_wait(req_b).ne.0) 
     $                 call ga_error('dra_wait failed:',req_b)
            nxt_idx = NEXT(cur_idx)
            req_b = prefetch_next(g_b_chunk(nxt_idx), d_b, i,j, nb, n,m)
c
c           handle the case when we have the data & there in nothing to prefetch
            if(req_b .eq. DRA_REQ_INVALID)then
                     nxt_idx = cur_idx
                     print *, 'invalid',i,j
            endif
c
*           call ga_print_patch(g_b_chunk(cur_idx), 1, n, 1, nc, 1)
            t0 =MP_TIMER()
            call GA_DGEMM('n','n', nr, nc, n, 1d0, 
     $                     g_a_chunk, g_b_chunk(cur_idx), 0d0,g_c_chunk)
            tot = MP_TIMER() -t0

*           call ga_print_patch(g_c_chunk, 1, nr, 1,nc, 1)
            cur_idx = nxt_idx

            if(dra_write_section(.false., g_c_chunk, 1, nr, 1, nc,
     &                          d_c, ilo, ihi, jlo, jhi, req_c).ne.0)
     &                       call ga_error('dra_write failed for C:',0)
            if(dra_wait(req_c).ne.0)call ga_error('wait failed:',req_c)
            if(ga_nodeid().eq.0)
     &         write(6,'(2i5,5x,2e15.4)')i,j,tot, MP_TIMER()
         enddo
      enddo
      if(dra_wait(req_b).ne.0) call ga_error('dra_wait failed:',req_b)
c
*     call ga_print(g_c, 1)
      status =  ga_destroy(g_c_chunk)
      status =  ga_destroy(g_b_chunk(2)) .and. status
      status =  ga_destroy(g_b_chunk(1)) .and. status
      status =  ga_destroy(g_a_chunk)    .and. status
      if (.not.status) print *,'ga_destroy failed'
      end

