#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine shuffle
#include "common.fh"
c
      integer snode,rnode,i,inx,iny,inz,pnum,idx,idy,idz,ipx,ipy,ipz
      integer icnt,inode,me
c
c   This subroutine shuffles the coordinates on each of the processors
c   around so that each processor has the coordinates for atoms residing
c   in the physical domain corresponding to that processor. This routine
c   assumes that atoms are completely randomized.
c
      pnum = ga_nnodes()
      me = ga_nodeid()
      call factor(pnum,idx,idy,idz)
      call i_proc_to_xyz(me,ipx,ipy,ipz,idx,idy,idz)
c
c   copy all particles into the buffer
c
      icnt = 0
      do i = 1, antot
        icnt = icnt + 1
        ra(i,1,6) = ra(i,1,1)
        ra(i,2,6) = ra(i,2,1)
        ra(i,3,6) = ra(i,3,1)
        xcrd(icnt) = ra(i,1,6)
        ycrd(icnt) = ra(i,2,6)
        zcrd(icnt) = ra(i,3,6)
        xfrc(icnt) = ra(i,1,2)
        yfrc(icnt) = ra(i,2,2)
        zfrc(icnt) = ra(i,3,2)
        mbuf(icnt) = mass(i)
        bidx(icnt) = aidx(i)
        bat(icnt)  = at(i)
      end do
      btot = antot
      antot = 0
c
c   send buffers to processes controlling neighboring domains,
c   and select out particles that are within that domain
c
c   send along x-axis and gather particles whose x value lies
c   between xmin and xmax
c
      call cull(1)
      do i = 1, idx-1
        inx = ipx - 1
        if (inx.eq.-1) inx = idx - 1
        call i_xyz_to_proc(rnode,inx,ipy,ipz,idx,idy,idz)
        call exchange_buf(rnode)
        call cull(1)
      end do
c
c   copy all particles into the buffer
c
      icnt = 0
      do i = 1, antot
        icnt = icnt + 1
        xcrd(icnt) = ra(i,1,6)
        ycrd(icnt) = ra(i,2,6)
        zcrd(icnt) = ra(i,3,6)
        xfrc(icnt) = ra(i,1,2)
        yfrc(icnt) = ra(i,2,2)
        zfrc(icnt) = ra(i,3,2)
        mbuf(icnt) = mass(i)
        bidx(icnt) = aidx(i)
        bat(icnt)  = at(i)
      end do
      btot = antot
      antot = 0
c
c   send along y-axis and gather particles whose y value lies
c   between ymin and ymax
c
      call cull(2)
      do i = 1, idy-1
        iny = ipy - 1
        if (iny.eq.-1) iny = idy - 1
        call i_xyz_to_proc(rnode,ipx,iny,ipz,idx,idy,idz)
        call exchange_buf(rnode)
        call cull(2)
      end do
c
c   copy all particles into the buffer
c
      icnt = 0
      do i = 1, antot
        icnt = icnt + 1
        xcrd(icnt) = ra(i,1,6)
        ycrd(icnt) = ra(i,2,6)
        zcrd(icnt) = ra(i,3,6)
        xfrc(icnt) = ra(i,1,2)
        yfrc(icnt) = ra(i,2,2)
        zfrc(icnt) = ra(i,3,2)
        mbuf(icnt) = mass(i)
        bidx(icnt) = aidx(i)
        bat(icnt)  = at(i)
      end do
      btot = antot
      antot = 0
c
c   send along z-axis and gather particles whose z value lies
c   between zmin and zmax
c
      call cull(3)
      do i = 1, idz-1
        inz = ipz - 1
        if (inz.eq.-1) inz = idz - 1
        call i_xyz_to_proc(rnode,ipx,ipy,inz,idx,idy,idz)
        call exchange_buf(rnode)
        call cull(3)
      end do
c
c   rearrange data stack so that they are ordered with respect to
c   atom index
c
      call heapsort(0)
      call fixper
c
      return
      end
c
      subroutine exchange_buf(rnode)
#include "common.fh"
c
      double precision buf(3,MAXAT), rbuf(MAXAT)
      integer i, ibuf(MAXAT)
      integer rnode, rptr, rtot, me, one, ld
c
c   exchange the size of lists
c
      me = ga_nodeid()
      one = 1
      ld = 3
c
      do i = 1, btot
        buf(1,i) = xcrd(i)
        buf(2,i) = ycrd(i)
        buf(3,i) = zcrd(i)
      end do
c
      rptr = gsize_lo(me)
      call nga_put(g_size,rptr,rptr,btot,one)
      gcoords_hi(1) = 3
      gcoords_hi(2) = gcoords_lo(2,me) + btot - 1
      if (btot.gt.0) call nga_put(g_coords,gcoords_lo(1,me),
     +                            gcoords_hi,buf,ld)
      call ga_sync()
      rptr = gsize_lo(rnode)
      call nga_get(g_size,rptr,rptr,rtot,one)
      gcoords_hi(2) = gcoords_lo(2,rnode) + rtot - 1
      if (rtot.gt.0) call nga_get(g_coords,gcoords_lo(1,rnode),
     +                            gcoords_hi,buf,ld)
c
      do i = 1, rtot
        xcrd(i) = buf(1,i)
        ycrd(i) = buf(2,i)
        zcrd(i) = buf(3,i)
      end do
      call ga_sync()
c
      do i = 1, btot
        buf(1,i) = xfrc(i)
        buf(2,i) = yfrc(i)
        buf(3,i) = zfrc(i)
      end do
c
      gcoords_hi(1) = 3
      gcoords_hi(2) = gcoords_lo(2,me) + btot - 1
      if (btot.gt.0) call nga_put(g_coords,gcoords_lo(1,me),
     +                            gcoords_hi,buf,ld)
      call ga_sync()
      gcoords_hi(2) = gcoords_lo(2,rnode) + rtot - 1
      if (rtot.gt.0) call nga_get(g_coords,gcoords_lo(1,rnode),
     +                            gcoords_hi,buf,ld)
c
      do i = 1, rtot
        xfrc(i) = buf(1,i)
        yfrc(i) = buf(2,i)
        zfrc(i) = buf(3,i)
      end do
      call ga_sync()
c
      do i = 1, btot
        rbuf(i) = mbuf(i)
      end do
c
      grvec_hi = grvec_lo(me) + btot - 1
      if (btot.gt.0) call nga_put(g_rvec,grvec_lo(me),grvec_hi,rbuf,ld)
      call ga_sync()
      grvec_hi = grvec_lo(rnode) + rtot - 1
      if (rtot.gt.0) call nga_get(g_rvec,grvec_lo(rnode),
     +                            grvec_hi,rbuf,ld)
c
      do i = 1, rtot
        mbuf(i) = rbuf(i)
      end do
      call ga_sync()
c
      do i = 1, btot
        ibuf(i) = bidx(i)
      end do
c
      gindex_hi = gindex_lo(me) + btot - 1
      if (btot.gt.0) call nga_put(g_index,gindex_lo(me),
     +                            gindex_hi,ibuf,ld)
      call ga_sync()
      gindex_hi = gindex_lo(rnode) + rtot - 1
      if (rtot.gt.0) call nga_get(g_index,gindex_lo(rnode),
     +                            gindex_hi,ibuf,ld)
c
      do i = 1, rtot
        bidx(i) = ibuf(i)
      end do
      call ga_sync()
c
      do i = 1, btot
        ibuf(i) = bat(i)
      end do
c
      gindex_hi = gindex_lo(me) + btot - 1
      if (btot.gt.0) call nga_put(g_index,gindex_lo(me),gindex_hi,
     +                            ibuf,ld)
      call ga_sync()
      gindex_hi = gindex_lo(rnode) + rtot - 1
      if (rtot.gt.0) call nga_get(g_index,gindex_lo(rnode),
     +                            gindex_hi,ibuf,ld)
c
      do i = 1, rtot
        bat(i) = ibuf(i)
      end do
      call ga_sync()
c
      btot = rtot
c
      return
      end
c
      subroutine cull(iflg)
#include "common.fh"
c
      double precision xmax,ymax,zmax,xmin,ymin,zmin
      double precision xt,yt,zt
      integer i,pnum,me,idx,idy,idz,ipx,ipy,ipz
      integer icnt,iflg
      logical xflg,yflg,zflg
c
c   This subroutine gathers all the particles in the buffer arrays that
c   have coordinates lying in the domain corresponding to the processor
c   and puts them in the regular particle arrays
c
c   Determine boundaries of the physical domain assigned to the processor
c
      pnum = ga_nnodes()
      me = ga_nodeid()
      call factor(pnum,idx,idy,idz)
      call i_proc_to_xyz(me,ipx,ipy,ipz,idx,idy,idz)
c
c  set logical flags
c
      if (iflg.eq.1) then
        xflg = .true.
        yflg = .false.
        zflg = .false.
      elseif (iflg.eq.2) then
        xflg = .false.
        yflg = .true.
        zflg = .false.
      elseif (iflg.eq.3) then
        xflg = .false.
        yflg = .false.
        zflg = .true.
      else
        call ga_error("Illegal direction in subroutine cull",iflg )
      endif
c
      xmax = xbox*dble(ipx+1)/dble(idx)
      ymax = ybox*dble(ipy+1)/dble(idy)
      zmax = zbox*dble(ipz+1)/dble(idz)
      xmin = xbox*dble(ipx)/dble(idx)
      ymin = ybox*dble(ipy)/dble(idy)
      zmin = zbox*dble(ipz)/dble(idz)
      xmax = xmax - xbox2
      ymax = ymax - ybox2
      zmax = zmax - zbox2
      xmin = xmin - xbox2
      ymin = ymin - ybox2
      zmin = zmin - zbox2
c
c   Locate all particles on the processor that should reside
c   on the processor and move all others to the buffer
c
      icnt = 0
      do i = 1, btot
        xt = xcrd(i) - xbox * anint(xcrd(i)/xbox)
        yt = ycrd(i) - ybox * anint(ycrd(i)/ybox)
        zt = zcrd(i) - zbox * anint(zcrd(i)/zbox)
        if (xt.eq.xbox2) xt = -xbox2
        if (yt.eq.ybox2) yt = -ybox2
        if (zt.eq.zbox2) zt = -zbox2
        if ((xflg.and.xt.lt.xmax.and.xt.ge.xmin).or.
     +      (yflg.and.yt.lt.ymax.and.yt.ge.ymin).or.
     +      (zflg.and.zt.lt.zmax.and.zt.ge.zmin)) then
          antot = antot + 1
          ra(antot,1,6) = xcrd(i)
          ra(antot,2,6) = ycrd(i)
          ra(antot,3,6) = zcrd(i)
          ra(antot,1,2) = xfrc(i)
          ra(antot,2,2) = yfrc(i)
          ra(antot,3,2) = zfrc(i)
          mass(antot) = mbuf(i)
          aidx(antot) = bidx(i)
          at(antot) = bat(i)
        else
          icnt = icnt + 1
          xcrd(icnt) = xcrd(i)
          ycrd(icnt) = ycrd(i)
          zcrd(icnt) = zcrd(i)
          xfrc(icnt) = xfrc(i)
          yfrc(icnt) = yfrc(i)
          zfrc(icnt) = zfrc(i)
          mbuf(icnt) = mbuf(i)
          bidx(icnt) = bidx(i)
          bat(icnt)  = bat(i)
        endif
      end do
      btot = icnt
c      if (btot.gt.MAXAT) then
c        call ga_error("btot greater than MAXAT in cull",btot)
c      endif
c
      return
      end
