#if HAVE_CONFIG_H
#   include "config.fh"
#endif

#if 0
#define USE_SUMMA 1
#endif

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
      double precision max_sz, max_disk, max_mem
      data max_arrays, max_sz, max_disk, max_mem /4,1d8,1d10, 1d6/
      integer status
      logical ga_create2
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
      integer n, nb
      integer g_temp
      integer d_a, d_b, d_c
      integer i, j, req
      integer nproc, me
      integer util_mdtob
      logical ga_create2
c     
      parameter (n = 10000)
      parameter (nb = n/5)
c
      double precision temp(n), mbytes
      double precision dsin, sum, tt0, tt1
c
      logical status
      external util_mdtob
c     ...................................................
c     
      nproc = ga_nnodes()
      me    = ga_nodeid()
c     
c***  Create global arrays
c
      if (.not. GA_CREATE(MT_DBL, n, n, 'temp', 1, 1, g_temp))
     $     call ga_error(' GA_CREATE failed: temp ',0)
c     
c***  creating disk arrays
c
      if (me.eq.0) print*, 'Creating Disk Array A:',n,' x',n
      if(dra_create(MT_DBL, n, n, 'A', 
     &    FNAME, 
     &    DRA_RW, nb, n, d_a).ne.0)
     $   CALL ga_error('dra_create failed: ',0)
      if (me.eq.0) print*, 'Creating Disk Array B:',n,' x',n
      if(dra_create(MT_DBL, n, n, 'B',
     &    FNAME1, 
     &    DRA_RW, n, nb, d_b).ne.0)
     $   CALL ga_error('dra_create failed: ',0)
      if (me.eq.0) print*, 'Creating Disk Array C:',n,' x',n
      if(dra_create(MT_DBL, n, n, 'C',
     &    FNAME2, 
     &    DRA_RW, nb, nb, d_c).ne.0)
     $   CALL ga_error('dra_create failed: ',0)
c
c*************************  Fill in arrays A & B  ************************
c
      do j = 1+me, n, nproc ! static load balancing
         do i = 1, n
            temp(i) = 1d0/(i +j)
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
c***  delete temp GA array
      status = ga_destroy(g_temp)
c***  opening disk arrays
c
      if(dra_open(FNAME,DRA_R, d_a).ne.0)
     &            call ga_error('dra_open failed',0)
      if(dra_open(FNAME1,DRA_R, d_b).ne.0)
     &            call ga_error('dra_open failed',0)
c
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
      call chunk_mxm(n, nb, d_a, d_b, d_c)
      tt1 = MP_TIMER() -tt0
      if(me.eq.0)print *,'time=',tt1,' MFLOPS=', 2d-6*n*n*n/tt1 
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
c     to help linker when summa is used
      if(ga_nnodes().le.0) call ga_dgemm()
100   format(g11.2,' MB  time=',g11.2,' rate=',g11.3,'MB/s')
      end


c     given current (i,j) compute (ni, nj) - next loop index
c     i - outermost loop, j- innermost loop
c     
      logical function next2d(i, imin, imax, j, jmin, jmax, ni, nj)
      implicit none
      integer i, imin, imax, j, jmin, jmax, ni, nj
c
      if(i .lt. imin .or. j.lt.jmin)then
c          to handle initial (0,0) indices
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
      next2d = ni.le.imax
      return
      end
        

      integer function prefetch_next(g_a, d_a, i, j, nb, n)
      implicit none
#include "global.fh"
#include "dra.fh"
#include "chemio.h"
      integer g_a, d_a, i, j, nb, n
      integer jlo, jhi, ni, nj, chunks, nc
      logical next2d
      integer LO, HI
      integer retval
      LO(i)  = (i-1)*nb +1
      HI(i)  = MIN(LO(i) + nb -1, n)
c
      chunks = (n+nb-1)/nb
      if(next2d(i, 1, chunks, j, 1, chunks, ni, nj)) then
        jlo = LO(nj)
        jhi = HI(nj)
        nc  = jhi -jlo+1
        if(ELIO_OK.eq.dra_read_section(.false., g_a, 1, n, 1, nc,
     &                     d_a, 1, n, jlo, jhi, retval)) then
           prefetch_next = retval
        else
           call ga_error('dra_read failed:',0)
        endif

      else
        prefetch_next = DRA_REQ_INVALID
      endif
      return
      end
   

      subroutine chunk_mxm(n, nb, d_a, d_b, d_c)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "dra.fh"
#include "mp3def.fh"
      integer n, nb, d_a, d_b, d_c
      integer i, j, chunks
      integer ilo, ihi, jlo, jhi, nr, nc
      integer g_a_chunk, g_b_chunk(2), g_c_chunk
      logical status
      integer req, req_b, req_c, prefetch_next
      integer nxt_idx, cur_idx
      double precision tot, t0
      logical ga_create2
      integer LO, HI, NEXT
      LO(i)  = (i-1)*nb +1
      HI(i)  = MIN(LO(i) + nb -1, n)
      NEXT(i)= Mod(i+2,2)+1        ! generates series: 1, 2, 1, 2, 1, ...
c
      chunks = (n+nb-1)/nb
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
      req_b = prefetch_next(g_b_chunk(cur_idx), d_b, 0, 0, nb, n)
      if(dra_wait(req_b).ne.0) call ga_error('dra_wait failed:',req)
 
      do i = 1, chunks
c
         ilo = LO(i)
         ihi = HI(i) 
         nr  = ihi-ilo+1
         if(ELIO_OK.ne.dra_read_section(.false., g_a_chunk, 1, nr, 1, n,
     &                               d_a, ilo, ihi, 1, n, req))
     &                          call ga_error('dra_read failed:',0)
         if(dra_wait(req).ne.0)  call ga_error('dra_wait failed:',req)
c
*        call ga_print_patch(g_a_chunk, 1, nr, 1, n)
         do j = 1, chunks
c
            jlo = LO(j)
            jhi = HI(j) 
            nc  = jhi-jlo+1

            if(dra_wait(req_b).ne.0) 
     $                 call ga_error('dra_wait failed:',req_b)
            nxt_idx = NEXT(cur_idx)
            req_b = prefetch_next(g_b_chunk(nxt_idx), d_b, i, j, nb, n)
c
*           call ga_print_patch(g_b_chunk(cur_idx), 1, n, 1, nc)
            t0 =MP_TIMER()
            call GA_DGEMM('n','n', nr, nc, n, 1d0, 
     $                     g_a_chunk, g_b_chunk(cur_idx), 0d0,g_c_chunk)
            tot = MP_TIMER() -t0

*           call ga_print(g_c_chunk, 1)
            cur_idx = nxt_idx

            if(ELIO_OK.ne.dra_write_section(
     &            .false., g_c_chunk, 1, nr, 1, nc,
     &            d_c, ilo, ihi, jlo, jhi, req_c))
     &         call ga_error('dra_write failed for C:',0)
            if(dra_wait(req_c).ne.0)call ga_error('wait failed:',req_c)
            if(ga_nodeid().eq.0)
     &         write(6,'(2i5,5x,2e15.4)')i,j,tot, MP_TIMER()
         enddo
      enddo
c
*     call ga_print(g_c, 1)
      status =  ga_destroy(g_c_chunk)
      status =  ga_destroy(g_b_chunk(2)) .and. status
      status =  ga_destroy(g_b_chunk(1)) .and. status
      status =  ga_destroy(g_a_chunk)    .and. status
      if(dra_wait(req_b).ne.0) call ga_error('dra_wait failed:',req_b)
      end

