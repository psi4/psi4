#if HAVE_CONFIG_H
#   include "config.fh"
#endif
#define IO_TEST_SZ (1024*1024/8)
#define IO_NUM_FILES 4
#define MAX_ITER 4
#define DMAP(x) dble(IO_TEST_SZ*((x-1)*eaf_util_szint()))
#define MAX_FILE_SZ 48.0d9
#define BASE_NAME 'EAFfile'
#ifdef  HPIODIR
#  define FNAME   HPIODIR//BASE_NAME
#else
#  define FNAME   BASE_NAME
#endif


      program eaftest
c     
c     call test1
c     call test2
c     call test3
      call test4
c     
      print *,' '
      print *,'All tests completed'
      end

      subroutine test2
      implicit none
#include "eaf.fh"
c     
      integer buf(IO_TEST_SZ, MAX_ITER), buf2(IO_TEST_SZ)
      character fname*200
      integer fnum,  iter
      integer i, nap
      real    dummy
      integer fd(IO_NUM_FILES)
      integer buf_sz
      integer aio(IO_NUM_FILES, MAX_ITER)
      integer stat2
c     
      buf_sz = (IO_TEST_SZ*eaf_util_szint())
c
      write(6,*)
      write(6,*) ' TEST 2 ... asynchronous read/write of multiple files'
      write(6,*) ' ----------------------------------------------------'
      write(6,*)
c     
 10   format (A,'.', I1)
c     
      do fnum = 1, IO_NUM_FILES
         write (fname, 10) FNAME, fnum
c     
         call testerr(eaf_open(fname, eaf_rw, fd(fnum)))
c     
         print *,'testing EAF_AWrite', ' file= ',fname
         do iter = 1, MAX_ITER
            do i = 1, IO_TEST_SZ
               buf(i,iter) = i + fnum + (iter*10000)
            enddo

            call testerr(eaf_awrite(fd(fnum),
     $           DMAP(iter),
     $           buf(1,iter), buf_sz, aio(fnum, iter)))

            write(0,*) 'Submitted EAF_AWriteF  fnum=', fnum,
     $           '   iter=',iter, '   aio=',aio(fnum, iter)
         enddo
      enddo

      do fnum = 1, IO_NUM_FILES
         do iter = 1, MAX_ITER
            write(0,*) 'Waiting to close  fnum=', fnum,
     $           '   iter=',iter
            call testerr(eaf_wait(fd(fnum), aio(fnum, iter)))
         enddo
         call eaf_print_stats(fd(fnum))
         call testerr(eaf_close(fd(fnum)))
      enddo

c.................................

      print *,'testing EAF_ARead'
      do fnum = 1, IO_NUM_FILES

         write (fname, 10) FNAME, fnum

         call testerr(eaf_open(fname, eaf_r, fd(fnum)))
         
         do iter = 1, MAX_ITER
            do i = 1, IO_TEST_SZ
               buf(i,iter) = i + fnum + (iter*10000)
            enddo

            call testerr(eaf_aread(fd(fnum), 
     $           DMAP(iter), 
     $           buf2, buf_sz, 
     $           aio(fnum,iter)))

 1000       call testerr(eaf_probe(aio(fnum, iter), stat2))
            write(0,999) stat2
 999        format(i5,i8)
            if(stat2 .ne. 0) then
               do nap=1, 2000000
                  dummy = 2.0**2.0+nap
               enddo
               goto 1000
            endif

#ifndef __crayx1
            call sleep(10)
#endif
            write(0,*) 'Finished reading file=',fnum, 
     $           '  fd=',fd(fnum),
     $           '  iter=',iter

            do i = 1, IO_TEST_SZ
               if(buf(i,iter) .ne. buf2(i)) then 
                  write(6,*) 'Bad read of file=',fnum, 
     $                 '  fd=',fd(fnum),
     $                 '  iter=',iter,'  offset=',i,'   buf1=',
     $                 buf(i,iter),'  buf2=',buf2(i)
                  stop 1
               endif
            enddo
         enddo
      enddo
c     
      do fnum = 1, IO_NUM_FILES
         call eaf_print_stats(fd(fnum))
         call testerr(eaf_close(fd(fnum)))
         write (fname, 10) FNAME, fnum
         call testerr(eaf_delete(fname))
      enddo
c     
      end
      subroutine test1
      implicit none
#include "eaf.fh"
      integer buf(IO_TEST_SZ), buf2(IO_TEST_SZ)
      character fname*200
      integer fnum,  iter
      integer i
      integer fd(IO_NUM_FILES)
      integer buf_sz
c     
      buf_sz = (IO_TEST_SZ*eaf_util_szint())
c
      write(6,*)
      write(6,*) ' TEST 1 ... synchronous read/write of multiple files'
      write(6,*) ' ---------------------------------------------------'
      write(6,*)
c     
      do fnum = 1, IO_NUM_FILES
         write (fname, 10) FNAME, fnum

 10      format (A,'.', I1)
         call testerr(eaf_open(fname, eaf_rw, fd(fnum)))

         print *, 'testing EAF_Write'
         do iter = 1, MAX_ITER
            do i = 1, IO_TEST_SZ
               buf(i) = i + fnum + (iter*10000)
            enddo
            call testerr(eaf_write(fd(fnum), 
     $           DMAP(iter),
     $           buf(1), buf_sz))
         enddo
      enddo
c     
      do fnum = 1, IO_NUM_FILES
         call eaf_print_stats(fd(fnum))
         call testerr(eaf_close(fd(fnum)))
      enddo
c     
c.................................
c     
      do fnum = 1, IO_NUM_FILES
         write (fname, 10) FNAME, fnum
c     
         call testerr(eaf_open(fname, eaf_r, fd(fnum)))
c     
         print *, 'testing EAF_Read'
         do iter = 1, MAX_ITER
            do i = 1, IO_TEST_SZ
               buf(i) = i + fnum + (iter*10000)
            enddo
c     
            call testerr(eaf_read(fd(fnum), 
     $           DMAP(iter), 
     $           buf2(1), buf_sz))
            do i = 1, IO_TEST_SZ
               if(buf(i) .ne. buf2(i)) then 
                  write(6,*) 'Bad read of file=',fnum,
     $                 '  iter=',iter,'  offset=',i,'   buf1=',
     $                 buf(i),'  buf2=',buf2(i)
                  stop 1
               endif
            enddo
         enddo
      enddo
c     
      do fnum = 1, IO_NUM_FILES
         call eaf_print_stats(fd(fnum))
         call testerr(eaf_close(fd(fnum)))
         write (fname, 10) FNAME, fnum
         call testerr(eaf_delete(fname))
      enddo
c     
      end
      subroutine testerr(code)
      implicit none
#include "eaf.fh"
      integer code
      character*80 msg
      integer fstrlen
      external fstrlen
c     
      call eaf_errmsg(code, msg)
      write(6,*) ' code ', code, ' = ', msg(1:fstrlen(msg))
***   stop 1
c
      end
      integer function fstrlen(a)
      implicit none
      character*(*) a
      integer i
      intrinsic len
c     
      do i = len(a),1,-1
         if (a(i:i) .ne. ' ') goto 10
      enddo
c     
 10   fstrlen = i
c     
      end
      subroutine test3
      implicit none
#include "eaf.fh"
      character fname*200
      character*8 fstype
      integer avail, fd
      double precision length
c
      write(6,*)
      write(6,*) ' TEST 3 ... misc. functions'
      write(6,*) ' --------------------------'
      write(6,*)
c
      write (fname, 10) FNAME, 1
 10   format (A,'.', I1)
c
      call testerr(eaf_stat(fname, avail, fstype))
      write(6,*) ' avail = ', avail, ' fstype = ', fstype
c
      call testerr(eaf_open(fname, eaf_rw, fd))
      call testerr(eaf_length(fd,length))
      write(6,*) ' file length = ', length
      call testerr(eaf_truncate(fd,129.0d0))
      call testerr(eaf_length(fd,length))
      write(6,*) ' file length (should be 129) = ', length
      call testerr(eaf_close(fd))
c
      write(6,*) ' This delete should succeed'
      call testerr(eaf_delete(fname))
c
      call testerr(eaf_open(fname, eaf_rw, fd))
      call testerr(eaf_length(fd,length))
      write(6,*) ' should be zero file length = ', length
      call testerr(eaf_close(fd))
      call testerr(eaf_delete(fname))
c
      end
      subroutine test4
      implicit none
#include "eaf.fh"
      character*200 fname
      character*255 fstype
      integer avail, fd, i, nwrite, ierr, iolen4, iolen, loop
      double precision length, offset, ofw
      double precision buf(IO_TEST_SZ)
c
      write(6,*)
      write(6,*) ' TEST 4 ... fill the file system '
      write(6,*) ' ------------------------------- '
      write(6,*)
c
      fname = ' '
      write (fname, 10) FNAME, 1
 10   format (A,'.', I1)
c
      call testerr(eaf_stat(fname, avail, fstype))
      write(6,*) ' avail = ', avail, ' fstype = ', fstype
      avail = 0.95*avail
      length = 1024*1024.0d0*avail
      call testerr(eaf_open(fname, eaf_rw, fd))
      if(length.gt. MAX_FILE_SZ) then
         length = MAX_FILE_SZ
         write(6,*) ' trying max file length ', length,' bytes'
      else
         write(6,*) ' 95% available length ', length,' bytes'
      endif
c
      ofw = 0
      do offset = 0, length, (IO_TEST_SZ)*8
         nwrite = nint(min(length-offset+1.0d0,(IO_TEST_SZ)*8.0d0))
         call util_dfill(nwrite/8, offset, buf, 1)
         if (ofw .le. offset) then
            write(6,*) ' writing at ',offset, nwrite
            ofw = ofw + 100000000.0d0
         endif
         ierr = eaf_write(fd, offset, buf, nwrite)
         if (ierr .ne. 0) then
           call testerr(ierr)
           call ga_error('bad write',ierr)
         endif
      enddo
c
      call eaf_print_stats(fd)
c
      ofw = 0
      do offset = 0, length, (IO_TEST_SZ)*8
         nwrite = nint(min(length-offset+1.0d0,(IO_TEST_SZ)*8.0d0))
         if (ofw .le. offset) then
            write(6,*) ' reading at ', offset, nwrite
            ofw = ofw + 100000000.0d0
         endif
         ierr = eaf_read(fd, offset, buf, nwrite)
         if (ierr .ne. 0) then
           call testerr(ierr)
           call ga_error('bad read',ierr)
         endif
         do i = 1, nwrite/8
            if (buf(i).ne.offset) call ga_error('bad ', i)
         enddo
      enddo
c
      call eaf_print_stats(fd)
c
      call testerr(eaf_close(fd))
c
c     random read and write
c
      buf(1) = eaf_util_random(35321)
      do iolen4 = 0,2
         iolen = 8*1024 * (4**iolen4)   ! 1024,4096,16384 doubles
         write(6,*) ' Random write: len=',iolen
         call testerr(eaf_open(fname, eaf_rw, fd))
         do loop = 1, int(length/iolen)/4
           offset = aint((length-iolen)*eaf_util_random(0))
           ierr = eaf_write(fd, offset, buf, iolen)
           if (ierr .ne. 0) then
             call testerr(ierr)
             call ga_error('bad write',ierr)
           endif
         enddo
         call eaf_print_stats(fd)
         write(6,*) ' Random read : len=',iolen
         call testerr(eaf_open(fname, eaf_rw, fd))
         do loop = 1, int(length/iolen)/4
           offset = aint((length-iolen)*eaf_util_random(0))
           ierr = eaf_read(fd, offset, buf, iolen)
           if (ierr .ne. 0) then
             call testerr(ierr)
             call ga_error('bad write',ierr)
           endif
         enddo
         call eaf_print_stats(fd)
         call testerr(eaf_close(fd))
      enddo
c
      call testerr(eaf_delete(fname))
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
      end

