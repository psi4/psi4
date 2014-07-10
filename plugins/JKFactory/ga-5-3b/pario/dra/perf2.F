#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c
c FNAME - filename for test program
c
#define BASE_NAME  'da.try'
#define BASE_NAME1 'da1.try'
#ifdef  HPIODIR
#  define FNAME   HPIODIR//BASE_NAME
#  define FNAME1  HPIODIR//BASE_NAME1
#else
#  define FNAME   BASE_NAME
#  define FNAME1  BASE_NAME1
#endif

      program io
#include "mafdecls.fh"
#include "global.fh"
#include "dra.fh"
      integer status, me
      integer max_arrays
      integer stack, heap
      double precision max_sz, max_disk, max_mem 
      data max_arrays, max_sz, max_disk, max_mem /10,1d8,1d10, 1d6/
#if   defined(IBM)|| defined(CRAY_T3E)
      data stack, heap /9000000, 4000000/
#else
      data stack, heap /1200000, 800000/
#endif
c
#include "mp3.fh"
      call ga_initialize()
      if(.not. ga_uses_ma())then
         stack = 10000 
         heap  = 10000 
      endif
c
      if(ma_init(MT_F_DBL, stack, heap) ) then
        me    = ga_nodeid()
        if(dra_init(max_arrays, max_sz, max_disk, max_mem).ne.0)then
                call ga_error('dra_init failed: ',0)
        endif
        if (me.eq.0) print *,'  '
        if(me.eq.0)print *, 'TESTING PERFORMANCE OF DISK ARRAYS'
        if (me.eq.0) print *,' '
        call test_io_dbl()
        status = dra_terminate()
        call ga_terminate()
      else
        print *,'ma_init failed'
      endif
      if(me.eq.0)print *, 'all done ...'
      call MP_FINALIZE()
      end



      subroutine test_io_dbl
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "dra.fh"
#include "mp3def.fh"
      integer n,m
      parameter (n=5000)
      parameter (m = 2*n)
      double precision   err, tt0, tt1, mbytes
      integer g_a, g_b,  d_a, d_b
      double precision drand
      integer i,j, req, loop
      integer dlo(2),dhi(2)
      integer glo(2),ghi(2)
      integer me, nproc
      integer index, dims(2), reqdims(2), ld(2)
      integer iran
      integer util_mdtob
      logical status
      external  drand
      external util_mdtob
      intrinsic int, dble
      iran(i) = int(drand(0)*dble(i-1)) + 1
c    
      loop  = 30
      req = -1
      nproc = ga_nnodes()
      me    = ga_nodeid()
c    
      if (me.eq.0) print *, 'creating global arrays ',n,' x',n
      if (me.eq.0)call ffflush(6)
      call ga_sync()
      if(.not. ga_create(MT_DBL, n, n, 'a', 1, 1, g_a))
     &         call ga_error('ga_create failed: a', 0)
      if(.not. ga_create(MT_DBL, n, n, 'b', 1, 1, g_b))
     &         call ga_error('ga_create failed: b', 0)
      if (me.eq.0) print *,'done '
      if (me.eq.0)call ffflush(6)
c
c     initialize g_a, g_b with random values
c     ... use ga_access to avoid allocating local buffers for ga_put
c
      call ga_sync()
      call nga_distribution(g_a, me, glo, ghi)
      call nga_access(g_a, glo, ghi, index, ld)
      call fill_random(DBL_MB(index), ghi(1)-glo(1)+1,
     &                 ghi(2)-glo(2)+1, ld(1))
      call ga_sync()
*     if (me.eq.0) print *,'done '
*     if (me.eq.0)call ffflush(6)
c
      call ga_zero(g_b)
c
c
c.......................................................................
      if (me.eq.0) print *, 'creating disk array ',n,' x',n
      if (me.eq.0)call ffflush(6)
      dims(1) = n
      dims(2) = n
      reqdims(1) = n
      reqdims(2) = n
      if(ndra_create(MT_DBL, 2, dims, 'A',
     &      FNAME, 
     &      DRA_RW, reqdims, d_a).ne.0)
     $   CALL ga_error('ndra_create failed: ',0)
c
      if(me.eq.0) print *, 'alligned blocking write'
      if (me.eq.0)call ffflush(6)
      tt0 = MP_TIMER()
      if(ndra_write(g_a, d_a, req).ne.0)
     $   CALL ga_error('ndra_write failed:',0)
      if(dra_wait(req).ne.0) call ga_error('dra_wait failed: ' ,req)
      tt1 = MP_TIMER() -tt0
      mbytes = 1e-6*util_mdtob(n*n)
      if (me.eq.0)then
          write(6,100)mbytes,tt1,mbytes/tt1
      endif
c
      if(dra_close(d_a).ne.0)call ga_error('dra_close failed: ',d_a)
      tt1 = MP_TIMER() -tt0
      if (me.eq.0)then
          write(6,100)mbytes,tt1,mbytes/tt1
      endif
c
      if (me.eq.0) print *,' '
      if (me.eq.0) print *,'disk array closed '
      if (me.eq.0)call ffflush(6)
c.......................................................................
c
c
      if (me.eq.0) print *, 'creating disk array ',m,' x',m
      if (me.eq.0)call ffflush(6)
      dims(1) = m
      dims(2) = m
      reqdims(1) = n
      reqdims(2) = n
      if(ndra_create(MT_DBL, 2, dims, 'B',
     &      FNAME1,
     &      DRA_RW, reqdims, d_b).ne.0)
     $   CALL ga_error('ndra_create failed: ',0)
c
      if(me.eq.0) print *, 'non alligned blocking write'
      if (me.eq.0)call ffflush(6)
c
      glo(1) = 1
      glo(2) = 1
      ghi(1) = n
      ghi(2) = n
      dlo(1) = 2
      dlo(2) = 2
      dhi(1) = n+1
      dhi(2) = n+1
      tt0 = MP_TIMER()
      if(ndra_write_section(.false., g_a, glo, ghi,
     &               d_b, dlo, dhi, req).ne.0)
     &         call  ga_error('ndra_write_section failed:',0)

      if(dra_wait(req).ne.0) call ga_error('dra_wait failed: ' ,req)
      tt1 = MP_TIMER() -tt0
      mbytes = 1e-6*util_mdtob(n*n)
      if (me.eq.0)then
          write(6,100)mbytes,tt1,mbytes/tt1
      endif
c
      if(dra_close(d_b).ne.0)call ga_error('dra_close failed: ',d_b)
      tt1 = MP_TIMER() -tt0
      if (me.eq.0)then
          write(6,100)mbytes,tt1,mbytes/tt1
      endif
c
      if (me.eq.0) print *,' '
      if (me.eq.0) print *,'disk array closed '
      if (me.eq.0)call ffflush(6)
c.......................................................................
c
c
      if (me.eq.0) print *,' '
      if (me.eq.0) print *,'opening disk array'
      if(dra_open(FNAME,DRA_R, d_a).ne.0)
     &            call ga_error('dra_open failed',0)
      if(me.eq.0) print *, 'alligned blocking read'
      if (me.eq.0)call ffflush(6)
      tt0 = MP_TIMER()
      if(ndra_read(g_b, d_a, req).ne.0)
     $   CALL ga_error('ndra_read failed:',0)
      if(dra_wait(req).ne.0) call ga_error('dra_wait failed: ' ,req)
      tt1 = MP_TIMER() -tt0
      if (me.eq.0)then
          write(6,100)mbytes,tt1,mbytes/tt1
      endif
      call ga_dadd(1d0, g_a, -1d0, g_b, g_b)
      err = ga_ddot(g_b, g_b)
      if(err.ne.0) then
          if (me.eq.0) print *,'BTW, we have error =', err
      else
          if (me.eq.0) print *,'OK'
      endif
      if(dra_delete(d_a).ne.0)
     &            call ga_error('dra_delete failed',0)
c.......................................................................
c
      if (me.eq.0) print *,' '
      if (me.eq.0) print *,'opening disk array'
      if(dra_open(FNAME1,DRA_R, d_b).ne.0)
     &            call ga_error('dra_open failed',0)
      if(me.eq.0) print *, 'non alligned blocking read'
      if (me.eq.0)call ffflush(6)
      tt0 = MP_TIMER()
      if(ndra_read_section(.false.,   g_b, glo, ghi,
     &               d_b, dlo, dhi, req).ne.0)
     &        call   ga_error('ndra_read_section failed:',0)
      if(dra_wait(req).ne.0) call ga_error('dra_wait failed: ' ,req)
      tt1 = MP_TIMER() -tt0
      if (me.eq.0)then
          write(6,100)mbytes,tt1,mbytes/tt1
      endif
      call ga_dadd(1d0, g_a, -1d0, g_b, g_b)
      err = ga_ddot(g_b, g_b)
      if(err.ne.0) then
          if (me.eq.0) print *,'BTW, we have error =', err
      else
          if (me.eq.0) print *,'OK'
      endif
      if(dra_delete(d_b).ne.0)
     &            call ga_error('dra_delete failed',0)
c.......................................................................
      status = ga_destroy(g_a)
      status = ga_destroy(g_b)
100   format(g11.2,' MB  time=',g11.2,' rate=',g11.3,'MB/s')
      end



      subroutine swap(a,b)
      integer a, b, temp
         temp = a
         a = b
         b = temp
      end


      subroutine init_char(str, len, char)
      character*(*) str
      character*1 char
      integer i
      do i = 1, len -1 
         str(i:i+1) = char
      enddo
      end
          
      
      subroutine fill_random(a, n,m, ld)
      integer ld, n,m
      double precision a(ld,*), drand, seed
      integer i,j
      external  drand
c
      do j=1,m
         seed = drand(j)
         do i=1,n
            a(i,j)=seed*i
         enddo
      enddo
      end
