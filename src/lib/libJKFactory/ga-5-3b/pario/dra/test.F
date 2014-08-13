#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c
c FNAME - filename for test program
c
#define BASE_NAME 'dra.file'
#ifdef  HPIODIR
#  define FNAME   HPIODIR//BASE_NAME
#else
#  define FNAME   BASE_NAME
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
      data stack, heap /80000, 80000/
c
#include "mp3.fh"
      if(ma_init(MT_F_DBL, stack, heap) ) then
        call ga_initialize()
        me    = ga_nodeid()
        if(dra_init(max_arrays, max_sz, max_disk, max_mem).ne.0)then
                call ga_error('dra_init failed: ',0)
        endif
        if (me.eq.0) print*,' '
        if(me.eq.0)print *, 'TESTING INTEGERS'
        if (me.eq.0) print*,' '
        call test_io_int()
        if (me.eq.0) print*,' '
        if(me.eq.0)print *, 'TESTING DOUBLES'
        if (me.eq.0) print*,' '
        call test_io_dbl()
        status = dra_terminate()
        call ga_terminate()
      else
        print *,'ma_init failed'
      endif
c     if(me.eq.0)print *, 'all done ...'
      call MP_FINALIZE()
      end



      subroutine test_io_int
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "dra.fh"
      integer n,m
      parameter (n = 100)
      parameter (m = 2*n)
      integer a(n,n)
      integer g_a, g_b, d_a
      integer i,j, req, err, dim1, dim2, type
      integer me, nproc
      logical status
      character*1 filename(200)
      character*80 name
c    
      nproc = ga_nnodes()
      me    = ga_nodeid()
      call init_char(name,100, ' ')
      call init_char(filename,200, ' ')
c    
c     a() is a local copy of what the l array should start as
c
      do j = 1, n
         do i = 1, n
            a(i,j) = i-1 + (j-1)*n
         enddo
      enddo
c
      if (me.eq.0) print *, 'Creating global arrays'
      call ga_sync()
      if(.not. ga_create(MT_INT, n, n, 'a', 1, 1, g_a))
     &         call ga_error('ga_create failed: a', 0)
      if(.not. ga_create(MT_INT, n, n, 'b', 1, 1, g_b))
     &         call ga_error('ga_create failed: b', 0)
c
*     call ga_zero(g_a)
*     call ga_zero(g_b)
c
      do j = 1+me, n, nproc
         call ga_put(g_a, 1, n, j, j, a(1, j), n)
      enddo
c
      if(me.eq.0) print *, 'Creating Disk Array ',n,' x',n
      if(dra_create(MT_INT, n, n, 'array A', 
     &   FNAME, 
     &   DRA_RW,    n, 1, d_a).ne.0)
     $   CALL ga_error('dra_create failed: ',0)
      if (me.eq.0) print *,'OK'
      if (me.eq.0) print*,' '
c
      if(me.eq.0) print *, 'Writing Global Array to Disk Array'
      if(dra_write(g_a, d_a,req).ne.0)
     $   CALL ga_error('dra_write failed:',0)
      if(dra_wait(req).ne.0)call ga_error('dra_wait failed: ' ,req)
      if (me.eq.0) print *,'OK'
      if (me.eq.0) print*,' '
      if(me.eq.0) print *, 'Closing Disk Array'
      if(dra_close(d_a).ne.0)call ga_error('dra_close failed: ',d_a)
      if (me.eq.0) print *,'OK'
      if (me.eq.0) print*,' '
c
      if(me.eq.0) print *, 'Opening Existing Disk Array'
      if(dra_open(FNAME ,DRA_R, d_a).ne.0)
     &            call ga_error('dra_open failed',0)
c     
      if(dra_inquire(d_a, type, dim1, dim2, name, filename).ne.0)
     &            call ga_error('dra_inquire failed',0)
      if(dim1.ne.n) call ga_error('dim1 error',dim1)
      if(dim2.ne.n) call ga_error('dim2 error',dim2)
      if(type.ne.MT_INT) call ga_error('type error',type)
      if(me.eq.0) print *, 'array name read from disk is:',name
      call ga_sync()
      if (me.eq.0) print *,'OK'
      if (me.eq.0) print*,' '
c
      if(me.eq.0) print *, 'Checking dra_read'
c     call dra_set_mode(1)
      if(dra_read(g_b, d_a, req).ne.0)
     $   CALL ga_error('dra_read failed:',0)
      if(dra_wait(req).ne.0) call ga_error('dra_wait failed: ' ,req)
c
c     error checking: (g_a - g_b)^2
c
c     call ga_print(g_a)
c     call ga_print(g_b)
      call ga_add(1, g_a, -1, g_b, g_b) 
      err = ga_idot(g_b, g_b) 
c
      if(err.ne.0 )then
          if( me.eq.0) call ga_error('failed', err) 
      else
          if (me.eq.0) print *,'OK'
      endif 
      if (me.eq.0) print*,' '
c
      if(me.eq.0) print *, 'Checking dra_delete'
      if(dra_delete(d_a).ne.0)
     &            call ga_error('dra_delete failed',0)
      if (me.eq.0) print *,'OK'
      if (me.eq.0) print*,' '
      status = ga_destroy(g_a)
      status = ga_destroy(g_b)
      end




      subroutine test_io_dbl
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "dra.fh"
      integer n,m
      parameter (n = 100)
      parameter (m = 2*n)
      double precision a(n,n),  err
      integer g_a, g_b,  d_a
      double precision drand
      integer i,j, req, loop
      integer dilo,dihi,djlo,djhi
      integer gilo,gihi,gjlo,gjhi
      integer ielem, jelem
      integer me, nproc
      integer iran
      logical status
      external  drand
      intrinsic int, dble
      iran(i) = int(drand(0)*dble(i-1)) + 1
c    
      loop  = 30
      nproc = ga_nnodes()
      me    = ga_nodeid()
c    
c     a() is a local copy of what the l array should start as
c
      do j = 1, n
         do i = 1, n
            a(i,j) = dble(i-1 + (j-1)*n)
         enddo
      enddo
c
      call ga_sync()
      if(.not. ga_create(MT_DBL, n, n, 'a', 1, 1, g_a))
     &         call ga_error('ga_create failed: a', 0)
      if(.not. ga_create(MT_DBL, n, n, 'b', 1, 1, g_b))
     &         call ga_error('ga_create failed: b', 0)
c
      call ga_zero(g_a)
      call ga_zero(g_b)
c
      do j = 1+me, n, nproc
         call ga_put(g_a, 1, n, j, j, a(1, j), n)
      enddo
c
      if (me.eq.0) print*, 'Creating Disk Array ',n,' x',n
      if(dra_create(MT_DBL, n, n, 'A', 
     &      FNAME, 
     &      DRA_RW, n, 3, d_a).ne.0)
     $   CALL ga_error('dra_create failed: ',0)
c
      if(me.eq.0) print *, 'Writing Global Array to Disk Array'
      if(dra_write(g_a, d_a,req).ne.0)
     $   CALL ga_error('dra_write failed:',0)
      if (me.eq.0) print*,' '
      if(dra_wait(req).ne.0) call ga_error('dra_wait failed: ' ,req)
c
      if(dra_close(d_a).ne.0)call ga_error('dra_close failed: ',d_a)
c
      if(me.eq.0) print *, 'Checking dra_read'
      if(dra_open(FNAME,DRA_R, d_a).ne.0)
     &            call ga_error('dra_open failed',0)
      if(dra_read(g_b, d_a, req).ne.0)
     $   CALL ga_error('dra_read failed:',0)
      if(dra_wait(req).ne.0) call ga_error('dra_wait failed: ' ,req)
c
c     error checking: (g_a - g_b)^2
c
c     call ga_print(g_a)
c     call ga_print(g_b)
      call ga_add(1d0, g_a, -1d0, g_b, g_b) 
      err = ga_ddot(g_b, g_b) 
      if(err.ne.0) then
          if (me.eq.0) print *,'error =', err 
      else
          if (me.eq.0) print *,'OK'
      endif 
      if (me.eq.0) print*,' '
c
      if(me.eq.0) print *, 'Checking dra_read_section'
c
      call ga_zero(g_b)
      do j = 1, loop
         dilo = iran(n)
         dihi = iran(n) 
         djlo = iran(n)
         djhi = iran(n) 
c
         if(dilo.gt.dihi) call swap(dilo,dihi)
         if(djlo.gt.djhi) call swap(djlo,djhi)
         ielem = dihi - dilo +1
         jelem = djhi - djlo +1
c
         gilo = iran(n-ielem)+1 
         gjlo = iran(n-jelem)+1 
         gihi = gilo+ielem-1
         gjhi = gjlo+jelem-1
c
         if(me.eq.0) then
            write(6,100) gilo,gihi,gjlo,gjhi,dilo,dihi,djlo,djhi 
100         format(1x,'reading  global[',2(i4,':',i4),
     &                ']  from  disk[',2(i4,':',i4),']')
            call ffflush(6)
         endif
c
*        call ga_print(g_b,1)
         if(dra_read_section(.false.,   g_b, gilo, gihi, gjlo, gjhi,
     &                      d_a, dilo, dihi, djlo, djhi, req).ne.0)
     &      call ga_error('dra_read failed:',0)
         if(dra_wait(req).ne.0) call ga_error('dra_wait failed:',req)
*        call ga_print(g_b,1)
*        call ga_print_patch(g_a, dilo, dihi, djlo, djhi,1)
*        call ga_print_patch(g_b, gilo, gihi, gjlo, gjhi,1)
c
         call ga_add_patch(1d0, g_a, dilo, dihi, djlo, djhi,
     &                     -1d0, g_b, gilo, gihi, gjlo, gjhi,
     &                           g_b, gilo, gihi, gjlo, gjhi) 
         err = ga_ddot_patch(g_b,'n', gilo, gihi, gjlo, gjhi,
     &                       g_b,'n', gilo, gihi, gjlo, gjhi) 
         if(err.ne.0 .and. me.eq.0)then
            print *,'error =', err
            call ga_error('failed',0)
         endif
      enddo
      if (me.eq.0) print*,' OK'
      if (me.eq.0) print*,' '
      if(dra_delete(d_a).ne.0)
     &            call ga_error('dra_delete failed',0)
c
c***  now d_a is 4 times larger than g_a
c
      if (me.eq.0) print*, 'Creating New Disk Array ',m,' x',m
      if(dra_create(MT_DBL, m, m, 'A', 
     &      FNAME, 
     &      DRA_RW, n, 2, d_a).ne.0)
     $   CALL ga_error('dra_create failed: ',0)
      if (me.eq.0) print*,' OK'
      if (me.eq.0) print*,' '
c
c   
c     call ga_print(g_a)
c     call ga_zero(g_b)
      if (me.eq.0) print*,'Testing dra_write_section'
      do j = 1, loop
         gilo = iran(n)
         gihi = iran(n)
         gjlo = iran(n)
         gjhi = iran(n)
c
         if(gilo.gt.gihi) call swap(gilo,gihi)
         if(gjlo.gt.gjhi) call swap(gjlo,gjhi)
         ielem = gihi - gilo +1
         jelem = gjhi - gjlo +1
c
         dilo = iran(m-ielem)+1
         djlo = iran(m-jelem)+1
         dihi = dilo+ielem-1
         djhi = djlo+jelem-1
c
         if(me.eq.0) then
            write(6,200) gilo,gihi,gjlo,gjhi,dilo,dihi,djlo,djhi
200         format(1x,'writing  global[',2(i4,':',i4),
     &                ']  to  disk[',2(i4,':',i4),']') 
            call ffflush(6)
         endif
c
c        call dra_set_mode(0)
         if(dra_write_section(.false.,  g_a, gilo, gihi, gjlo, gjhi,
     &                      d_a, dilo, dihi, djlo, djhi, req).ne.0)
     &      call ga_error('dra_write failed:',0)
         if(dra_wait(req).ne.0) call ga_error('dra_wait failed:',req)
c
c***     dra_read was tested already and we use it for testing dra_write_section
c
c        call dra_set_mode(1)
         if(dra_read_section(.false.,   g_b, gilo, gihi, gjlo, gjhi,
     &                      d_a, dilo, dihi, djlo, djhi, req).ne.0)
     &      call ga_error('dra_read failed:',0)
         if(dra_wait(req).ne.0) call ga_error('dra_wait failed:',req)
c
         call ga_add_patch(1d0, g_a, gilo, gihi, gjlo, gjhi,
     &                     -1d0, g_b, gilo, gihi, gjlo, gjhi,
     &                           g_b, gilo, gihi, gjlo, gjhi)
         err = ga_ddot_patch(g_b,'n', gilo, gihi, gjlo, gjhi,
     &                       g_b,'n', gilo, gihi, gjlo, gjhi)
         if(err.ne.0d0 )then
            call ga_print_patch(g_a, gilo, gihi,gjlo,gjhi,1)
            call ga_print_patch(g_b, gilo, gihi,gjlo,gjhi,1)
         endif
         if(err.ne.0d0 .and. me.eq.0)then
            print *,'error =', err
            call ga_error('error in dra_write_section',0)
         endif
      enddo
      if (me.eq.0) print*,' OK'
c
      if(dra_delete(d_a).ne.0)
     &            call ga_error('dra_delete failed',0)
      status = ga_destroy(g_a)
      status = ga_destroy(g_b)
      end



      subroutine swap(a,b)
      integer a, b, temp
         temp = a
         a = b
         b = temp
      end


      subroutine init_char(str, len, char)
      integer len, i
#if defined(CRAY_T3D) || defined(CRAY_T3E)
      character*1 str(len)
#else
      character*(*) str
#endif
      character*1 char
      do i = 1, len -1 
         str(i:i+1) = char
      enddo
      end
          
      
