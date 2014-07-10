#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c
c FNAME - filename for test program
c
#define BASE_NAME 'shared.file'
#ifdef  HPIODIR
#  define FNAME   HPIODIR//BASE_NAME
#else
#  define FNAME   BASE_NAME
#endif
#define SYNCH ga_sync

      program sftest
c
      implicit none
#include "global.fh"
#include "sf.fh"
c
      integer size
      integer dimnsn
      integer maxid
      parameter(maxid=5)           ! max number of outstanding I/O requests
      parameter(dimnsn=8192)
      real*8 buffer(dimnsn,maxid)  ! need buffering for all maxid requests
      double precision tt0,ttw, ttr
c     integer stack, heap
      integer idlist(maxid), curid
      character*80 errmsg
c
      integer me,nproc,handle,i,start,end,j, rc, chunk
c
#include "mp3.fh"
      call ga_initialize()
c
      curid = 0
      me = ga_nodeid()
      nproc = ga_nnodes()
      size  = maxid*dimnsn*nproc
c
      rc = sf_create(FNAME,
     $     dble(16*size),dble(8*size),dble(8*dimnsn),handle)
c
      call SYNCH
      chunk = (size+nproc-1)/nproc
      start = me*chunk+1
      end = min((start+chunk-1),size)
      tt0 = MP_TIMER()
c
      print *,'me=',me,'writing:', start, end
c     everybody writes chunk of data
      if(start.le.end) then
         do i = start, end,dimnsn
            do j = 1, min(dimnsn,(end-i+1))
               buffer(j,curid+1) = dble(i+j-1)
            enddo

            if(curid .eq. maxid)then
               rc=sf_waitall(idlist,maxid)
               curid = 0
            endif
            curid = curid+1
            rc = sf_write(handle,  8*dble(i-1), 
     &           8*dble(min(dimnsn,(end-i+1))), buffer(1,curid), 
     &           idlist(curid))
            if (rc.ne.0)call ga_error('write failed',rc)

         enddo
      endif
c
      call sf_fsync(handle)
      rc=sf_waitall(idlist,curid)
      if(rc.ne.0)call ga_error('waitall failed',rc)
      curid = 0

      ttw = MP_TIMER() -tt0
c
      call ga_dgop(88,ttw,1,'max')
      call SYNCH
c
c
c     everybody reads different chunk of data
      start = (nproc-me-1)*chunk+1
      end = min((start+chunk-1),size)
      print *,'me=',me,'reading:', start, end
      tt0 = MP_TIMER()
      do i = start,end,dimnsn

c           read and test data chunk by chunk
            rc = sf_read(handle, 8*dble(i-1),
     .           8*dble(min(dimnsn,(end-i+1))), buffer, idlist(1))
            if (rc.ne.0)then
               call sf_errmsg(rc,errmsg) 
               print *,'read at offset ',8*dble(i-1),' failed:',errmsg
               call ga_error('read failed',rc)
            endif
            rc=sf_wait(idlist(1))
            if (rc.ne.0)call ga_error('wait failed',rc)
c
            do j = 1,min(dimnsn,(end-i+1))
               if(buffer(j,1).ne.dble(i+j-1)) then
                  print *, me,buffer(j,1), i+j-1,i 
                  stop 'test failed'
               endif
            enddo
      enddo
      ttr = MP_TIMER() -tt0
c
      call ga_dgop(88,ttr,1,'max')
      call SYNCH
c
      rc = sf_destroy(handle)
      if(me.eq.0)then
        write(*,*)' '
        write(*,*)'test passed ', 8*maxid*dimnsn,' bytes'
        write(*,*) 8.0e-6*maxid*dimnsn/ttw,' MB/s write rate'
        write(*,*) 8.0e-6*maxid*dimnsn/ttr,' MB/s read rate'
        write(*,*)' '
      endif
c
      call SYNCH
      call ga_terminate
      call MP_FINALIZE()
      end
