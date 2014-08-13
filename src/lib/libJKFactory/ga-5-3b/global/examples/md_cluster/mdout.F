#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c
c                                   NOTICE
c
c   This software is being made available for internal testing and
c   evaluation purposes only. This software is a pre-release test version
c   which has not yet been authenticated and cleared for publication. Adherence
c   to this notice may be necessary for the author, Battelle Memorial
c   Institute, to successfully assert copyright in and commercialize this
c   software. This software is not intended for duplication or distribution
c   to third parties without the permission of the Manager of Software
c   Products at Pacific Northwest Laboratory, Richland, Washington,  99352.
c
      subroutine mdout
#include "common.fh"
c
      integer pnum,idx,idy,idz,me
c
c   This routine finishes off the calculation by writing
c   out accumulated information, final configurations etc.
c
      pnum = ga_nnodes()
      me = ga_nodeid()
      call factor(pnum,idx,idy,idz)
      if (ga_nodeid().eq.0.and.l_stdio) then
        write(6,100) pnum
        write(6,200) idx,idy,idz
        write(6,300) failcount
      endif
c
      call prtnrg
      call tmout
      call cluster_print_binr
      if (l_rst) call mdrst
c
      return
  100 format('Total number of processors used in simulation      :',i9)
  200 format('Processors configured as ',i2,' X ',i2,' X ',i2,' grid')
  300 format('Total number of collision failures                 :',i9)
      end
c
      subroutine header
#include "common.fh"
c
      integer oldstp
      save oldstp
      data oldstp /0/
c
c   this prints out the current MD time if time has
c   changed since the last call to this subroutine
c
      if (ga_nodeid().ne.0) return
      if (istep.ne.oldstp.and.l_stdio) then
         write(6,1000) float(istep) * tau
         oldstp = istep
      endif
      return
 1000 format(1x,'Statistics for MD time ',f16.4)
      end
c
      subroutine mdrst
#include "common.fh"
c
      double precision x,y,z,rx,ry,rz
      double precision rbuf(3,MAXAT)
      integer i,inode,pnum,me,one,ld2,ld3
      integer rptr,rnode,rtot,ibuf(2,MAXAT)
      character*32 filename
c
c   This routine finishes off the calculation by writing
c   out accumulated information, final configurations etc.
c
c
c   Write out final configuration to file 'md.rst'
c
      pnum = ga_nnodes()
      me = ga_nodeid()
      call sort
c
c   write all coordinates and velocities to buffers
c
      do i = 1, antot
        xcrd(i) = ra(i,1,1)
        ycrd(i) = ra(i,2,1)
        zcrd(i) = ra(i,3,1)
        xfrc(i) = ra(i,1,2)
        yfrc(i) = ra(i,2,2)
        zfrc(i) = ra(i,3,2)
        bidx(i) = aidx(i)
        bat(i) = at(i)
      end do
      btot = antot
c
c  open file from node 0 and write out box size
c
      if (task_id.lt.10) then
        write(filename,100) task_id
      else if (task_id.ge.10.and.task_id.lt.100) then
        write(filename,101) task_id
      else if (task_id.ge.100.and.task_id.lt.1000) then
        write(filename,102) task_id
      else if (task_id.ge.1000.and.task_id.lt.10000) then
        write(filename,103) task_id
      endif
  100 format('md.rst',i1)
  101 format('md.rst',i2)
  102 format('md.rst',i3)
  103 format('md.rst',i4)
      if (me.eq.0) then
        open(unit=10,file=filename,status='unknown')
        write(10,5400) atot
        write(10,5300) xbox,ybox,zbox,r_cluster
      endif
c
c     get data from other nodes and write to md.rst
c
      one = 1
      ld2 = 2
      ld3 = 3
c
c   exchange the size of lists and put coordinate data in
c   global arrays
c
      rptr = gsize_lo(me)
      call nga_put(g_size,rptr,rptr,btot,one)
c
      gcoords_hi(1) = 3
      gcoords_hi(2) = gcoords_lo(2,me) + btot - 1
      giat_hi(1) = 2
      giat_hi(2) = giat_lo(2,me) + btot - 1
      do i = 1, btot
        rbuf(1,i) = xcrd(i)
        rbuf(2,i) = ycrd(i)
        rbuf(3,i) = zcrd(i)
        ibuf(1,i) = bidx(i)
        ibuf(2,i) = bat(i)
      end do
      if (btot.gt.0) call nga_put(g_coords,gcoords_lo(1,me),
     +                            gcoords_hi,rbuf,ld3)
      if (btot.gt.0) call nga_put(g_iat,giat_lo(1,me),
     +                            giat_hi,ibuf,ld2)
c
      gfrc_hi(1) = 3
      gfrc_hi(2) = gfrc_lo(2,me) + btot - 1
      do i = 1, btot
        rbuf(1,i) = xfrc(i)
        rbuf(2,i) = yfrc(i)
        rbuf(3,i) = zfrc(i)
      end do
      if (btot.gt.0) call nga_put(g_frc,gfrc_lo(1,me),gfrc_hi,rbuf,ld3)
      call ga_sync()
c
      do inode = 0, pnum - 1
        if (me.eq.0) then
          rptr = gsize_lo(inode)
          call nga_get(g_size,rptr,rptr,rtot,one)
          gcoords_hi(2) = gcoords_lo(2,inode) + rtot - 1
          giat_hi(2) = giat_lo(2,inode) + rtot - 1
          if (rtot.gt.0) call nga_get(g_coords,gcoords_lo(1,inode),
     +                                gcoords_hi,rbuf,ld3)
          if (rtot.gt.0) call nga_get(g_iat,giat_lo(1,inode),
     +                                giat_hi,ibuf,ld2)
          do i = 1, rtot
            xcrd(i) = rbuf(1,i)
            ycrd(i) = rbuf(2,i)
            zcrd(i) = rbuf(3,i)
            bidx(i) = ibuf(1,i)
            bat(i) = ibuf(2,i)
          end do
          gfrc_hi(2) = gfrc_lo(2,inode) + rtot - 1
          if (rtot.gt.0) call nga_get(g_frc,gfrc_lo(1,inode),
     +                                gfrc_hi,rbuf,ld3)
          do i = 1, rtot
            xfrc(i) = rbuf(1,i)
            yfrc(i) = rbuf(2,i)
            zfrc(i) = rbuf(3,i)
          end do
c
          do i = 1, rtot
            if (istop.eq.1) then
              write(10,5000) bat(i),xcrd(i),ycrd(i),zcrd(i)
            else
              write(10,5100) bat(i),xcrd(i),ycrd(i),zcrd(i),
     +                       xfrc(i),yfrc(i),zfrc(i)
            endif
          end do
        endif
      end do
      if (me.eq.0) close(10)
c
 5000 format (i8,3(1pe13.5))
 5100 format (i8,6(1pe13.5))
 5300 format (4(1pe16.8))
 5400 format (i8)
c
      return
      end
