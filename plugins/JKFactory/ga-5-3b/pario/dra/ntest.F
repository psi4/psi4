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
      parameter (n = 20)
      parameter (m = 2*n)
      integer a(n,n,n)
      integer g_a, g_b, d_a
      integer i,j,k, err, type
      integer me, nproc, dims(3), req(3), ndim
      integer chunk(3), lo(3), hi(3), ld(3)
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
      do k = 1, n
        do j = 1, n
          do i = 1, n
            a(i,j,k) = i-1 + (j-1)*n + (k-1)*n*n
          enddo
        enddo
      enddo
c
      if (me.eq.0) print *, 'Creating global arrays'
      call ga_sync()
      do i = 1, 3
        dims(i) = n
        chunk(i) = 1
        ld(i) = n
      end do
      if(.not. nga_create(MT_INT, 3, dims, 'a', chunk, g_a))
     &         call ga_error('nga_create failed: a', 0)
      if(.not. nga_create(MT_INT, 3, dims, 'b', chunk, g_b))
     &         call ga_error('nga_create failed: b', 0)
c
*     call nga_zero(g_a)
*     call nga_zero(g_b)
c
      write(6,*) 'Imediately before nga_put'
      do j = 1+me, n, nproc
        lo(1) = 1
        lo(2) = 1
        lo(3) = j
        hi(1) = n
        hi(2) = n
        hi(3) = j
        call nga_put(g_a, lo, hi, a(1,1,j), ld)
      enddo
c
      if(me.eq.0) print *, 'Creating Disk Array ',n,' x',n,' x',n
      dims(1) = n
      dims(2) = n
      dims(3) = n
      req(1) = n
      req(2) = n
      req(3) = 1
      ndim = 3
      if(ndra_create(MT_INT, ndim, dims, 'array A', 
     &   FNAME, 
     &   DRA_RW,    req, d_a).ne.0)
     $   CALL ga_error('ndra_create failed: ',0)
      if (me.eq.0) print *,'OK'
      if (me.eq.0) print*,' '
c
      if(me.eq.0) print *, 'Writing Global Array to Disk Array'
      if(ndra_write(g_a, d_a,req).ne.0)
     $   CALL ga_error('ndra_write failed:',0)
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
      if(ndra_inquire(d_a, type, ndim, dims, name, filename).ne.0)
     &            call ga_error('ndra_inquire failed',0)
      if(dims(1).ne.n) call ga_error('dim1 error',dims(1))
      if(dims(2).ne.n) call ga_error('dim2 error',dims(2))
      if(dims(3).ne.n) call ga_error('dim3 error',dims(3))
      if(type.ne.MT_INT) call ga_error('type error',type)
      if(me.eq.0) print *, 'array name read from disk is:',name
      call ga_sync()
      if (me.eq.0) print *,'OK'
      if (me.eq.0) print*,' '
c
      if(me.eq.0) print *, 'Checking ndra_read'
c     call dra_set_mode(1)
      write(6,*) 'About to perform ndra_read'
      if(ndra_read(g_b, d_a, req).ne.0)
     $   CALL ga_error('ndra_read failed:',0)
      write(6,*) 'Finished performing ndra_read'
      call ffflush(6)
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
      parameter (n = 20)
      parameter (m = 2*n)
      double precision a(n,n,n),  err
      integer g_a, g_b,  d_a
      double precision drand
      integer i,j,k, loop
      integer dlo(3),dhi(3)
      integer glo(3),ghi(3)
      integer elem
      integer me, nproc
      integer iran, ndim, dims(3), req(3), chunk(3), ld(3)
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
      do k = 1, n
        do j = 1, n
          do i = 1, n
            a(i,j,k) = dble(i-1 + (j-1)*n + (k-1)*n*n)
          enddo
        enddo
      end do
c
      call ga_sync()
      ndim = 3
      do i = 1, ndim
        dims(i) = n
        chunk(i) = 1
        ld(i) = n
      end do
      write(6,*) 'Creating global arrays'
      write(6,*)
      if(.not. nga_create(MT_DBL, ndim, dims, 'a', chunk, g_a))
     &         call ga_error('ga_create failed: a', 0)
      if(.not. nga_create(MT_DBL, ndim, dims, 'b', chunk, g_b))
     &         call ga_error('ga_create failed: b', 0)
c
      write(6,*) 'Zeroing global arrays'
      write(6,*)
      call ga_zero(g_a)
      call ga_zero(g_b)
c
      do j = 1+me, n, nproc
        dlo(1) = 1
        dlo(2) = 1
        dlo(3) = j
        dhi(1) = n
        dhi(2) = n
        dhi(3) = j
        call nga_put(g_a, dlo, dhi, a(1, 1, j), ld)
      enddo
c
      if (me.eq.0) print*, 'Creating Disk Array ',n,' x',n,' x',n
      req(1) = n
      req(2) = n
      req(3) = 3
      if(ndra_create(MT_DBL, ndim, dims, 'A', 
     &      FNAME, 
     &      DRA_RW, req, d_a).ne.0)
     $   CALL ga_error('ndra_create failed: ',0)
c
      if(me.eq.0) print *, 'Writing Global Array to Disk Array'
      if(ndra_write(g_a, d_a,req).ne.0)
     $   CALL ga_error('ndra_write failed:',0)
      if (me.eq.0) print*,' '
      if(dra_wait(req).ne.0) call ga_error('dra_wait failed: ' ,req)
c
      if(dra_close(d_a).ne.0)call ga_error('dra_close failed: ',d_a)
c
      if(me.eq.0) print *, 'Checking ndra_read'
      if(dra_open(FNAME,DRA_R, d_a).ne.0)
     &            call ga_error('dra_open failed',0)
      if(ndra_read(g_b, d_a, req).ne.0)
     $   CALL ga_error('ndra_read failed:',0)
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
      if(me.eq.0) print *, 'Checking ndra_read_section'
c
      call ga_zero(g_b)
      do j = 1, loop
         do i = 1, ndim
           dlo(i) = iran(n)
           dhi(i) = iran(n)
           if (dlo(i).gt.dhi(i)) call swap(dlo(i),dhi(i))
           elem = dhi(i) - dlo(i) + 1
           glo(i) = iran(n-elem) + 1
           ghi(i) = glo(i) + elem - 1
         end do
c
         if(me.eq.0) then
            write(6,100) (glo(i),ghi(i),i=1,3),(dlo(i),dhi(i),i=1,3)
100         format(1x,'reading  global[',3(i4,':',i4),
     &                ']  from  disk[',3(i4,':',i4),']')
            call ffflush(6)
         endif
c
*        call ga_print(g_b,1)

         if(ndra_read_section(.false.,   g_b, glo, ghi,
     &                      d_a, dlo, dhi, req).ne.0)
     &      call ga_error('ndra_read failed:',0)
         if(dra_wait(req).ne.0) call ga_error('dra_wait failed:',req)
*        call ga_print(g_b,1)
*        call ga_print_patch(g_a, dilo, dihi, djlo, djhi,1)
*        call ga_print_patch(g_b, gilo, gihi, gjlo, gjhi,1)
c
         call nga_add_patch(1d0, g_a, dlo, dhi, -1d0, g_b, glo, ghi,
     &                           g_b, glo, ghi) 
         err = nga_ddot_patch(g_b,'n', glo, ghi, g_b,'n', glo, ghi) 
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
      if (me.eq.0) print*, 'Creating New Disk Array ',m,' x',m,' x',m
      dims(1) = m
      dims(2) = m
      dims(3) = m
      req(1) = n
      req(2) = n
      req(3) = 2
      if(ndra_create(MT_DBL, ndim, dims, 'A', 
     &      FNAME, 
     &      DRA_RW, req, d_a).ne.0)
     $   CALL ga_error('dra_create failed: ',0)
      if (me.eq.0) print*,' OK'
      if (me.eq.0) print*,' '
c
c   
c     call ga_print(g_a)
c     call ga_zero(g_b)
      if (me.eq.0) print*,'Testing ndra_write_section'
      do j = 1, loop
        do i = 1, ndim
          glo(i) = iran(n)
          if(glo(i).gt.ghi(i)) call swap(glo(i),ghi(i))
          elem = ghi(i) - glo(i) +1
          dlo(i) = iran(m-elem)+1
          dhi(i) = dlo(i)+elem-1
        end do
c
         if(me.eq.0) then
            write(6,200) (glo(i),ghi(i),i=1,3),(dlo(i),dhi(i),i=1,3)
200         format(1x,'writing  global[',3(i4,':',i4),
     &                ']  to  disk[',3(i4,':',i4),']') 
            call ffflush(6)
         endif
c
c        call dra_set_mode(0)
c         if(dra_write_section(.false.,  g_a, gilo, gihi, gjlo, gjhi,
c     &                      d_a, dilo, dihi, djlo, djhi, req).ne.0)
         if(ndra_write_section(.false., g_a, glo, ghi,
     &                        d_a, dlo, dhi, req).ne.0)
     &      call ga_error('dra_write failed:',0)
         if(dra_wait(req).ne.0) call ga_error('dra_wait failed:',req)
c
c***     dra_read was tested already and we use it for testing ndra_write_section
c
c        call dra_set_mode(1)
         if(ndra_read_section(.false.,   g_b, glo, ghi,
     &                      d_a, dlo, dhi, req).ne.0)
     &      call ga_error('ndra_read failed:',0)
         if(dra_wait(req).ne.0) call ga_error('dra_wait failed:',req)
c
         call nga_add_patch(1d0, g_a, glo, ghi, -1d0, g_b, glo, ghi,
     &                           g_b, glo, ghi)
         err = nga_ddot_patch(g_b,'n', glo, ghi, g_b,'n', glo, ghi)
cbjp         if(err.ne.0d0 )then
cbjp            call ga_print_patch(g_a, gilo, gihi,gjlo,gjhi,1)
cbjp            call ga_print_patch(g_b, gilo, gihi,gjlo,gjhi,1)
cbjp         endif
         if(err.ne.0d0 .and. me.eq.0)then
            print *,'error =', err
            call ga_error('error in ndra_write_section',0)
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
          
      
