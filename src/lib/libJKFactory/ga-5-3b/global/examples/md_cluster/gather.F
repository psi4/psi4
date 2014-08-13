#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine gather
#include "common.fh"
c
      double precision xmax,ymax,zmax,xmin,ymin,zmin
      double precision tbeg,wraptime
      integer snode,rnode,pnum,idx,idy,idz,ipx,ipy,ipz
      integer i,is,ig,inx,iny,inz,stot,inode
      logical nolist
c
c   This subroutine gathers the coordinates of all particles on
c   neighboring processors that lie within an interaction distance
c   of the domain boundaries corresponding to the local processor
c
      tbeg = wraptime()
c
      if (istep.eq.0.or.(mod(istep,ilist).eq.0.and.
     +    t_rmndr.eq.0.0d00)) then
        nolist = .true.
      else
        nolist = .false.
      endif
c
      pnum = ga_nnodes()
      call factor(pnum,idx,idy,idz)
      inode = ga_nodeid()
      call i_proc_to_xyz(inode,ipx,ipy,ipz,idx,idy,idz)
c      if (task_id.gt.3) then
c        write(6,*) ga_pgroup_nodeid(ga_pgroup_get_world()),
c     +     'Nolist is ',nolist
c        write(6,157) ga_pgroup_nodeid(ga_pgroup_get_world()),
c     +     idx,idy,idz
c  157   format(i3,' Proc dimensions are ',3i6)
c      endif
c
c   determine local boundaries
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
      btot = 0
c
c    Send particles along x-axis. Gather all particles that are within
c    an interaction distance of the x-boundaries and send them to the
c    adjacent processors
c
      if (idx.eq.2) then
c
c    both boundaries are adjacent to boundaries on next processor
c
        stot = 0
        if (nolist) then
          ixgru = 0
          do i = 1, antot
            if (xmax-ra(i,1,1).lt.rcmax.or.
     +         (ra(i,1,1)-xmin).lt.rcmax) then
              ixgru = ixgru + 1
              xgru(ixgru) = i
              stot = stot + 1
              is = btot+stot
              xcrd(is) = ra(i,1,1)
              ycrd(is) = ra(i,2,1)
              zcrd(is) = ra(i,3,1)
              bidx(is) = aidx(i)
              bat(is) = at(i)
            endif
          end do
        else
          do i = 1, ixgru
            ig = xgru(i)
            stot = stot + 1
            is = btot+stot
            xcrd(is) = ra(ig,1,1)
            ycrd(is) = ra(ig,2,1)
            zcrd(is) = ra(ig,3,1)
            bidx(is) = aidx(ig)
            bat(is) = at(ig)
          end do
        endif
c
c    send boundary layer particles to next buffer
c
        inx = ipx - 1
        if (inx.eq.-1) inx = idx - 1
        call i_xyz_to_proc(rnode,inx,ipy,ipz,idx,idy,idz)
        call gather_buf(rnode,stot)
      elseif (idx.gt.2) then
c
c   each boundary is adjacent to a different processor
c
        stot = 0
        if (nolist) then
          ixgru = 0
          do i = 1, antot
            if (xmax-ra(i,1,1).lt.rcmax) then
              ixgru = ixgru + 1
              xgru(ixgru) = i
              stot = stot + 1
              is = btot+stot
              xcrd(is) = ra(i,1,1)
              ycrd(is) = ra(i,2,1)
              zcrd(is) = ra(i,3,1)
              bidx(is) = aidx(i)
              bat(is) = at(i)
            endif
          end do
        else
          do i = 1, ixgru
            ig = xgru(i)
            stot = stot + 1
            is = btot+stot
            xcrd(is) = ra(ig,1,1)
            ycrd(is) = ra(ig,2,1)
            zcrd(is) = ra(ig,3,1)
            bidx(is) = aidx(ig)
            bat(is) = at(ig)
          end do
        endif
c
c    send boundary layer particles to next buffer
c
        inx = ipx - 1
        if (inx.eq.-1) inx = idx - 1
        call i_xyz_to_proc(rnode,inx,ipy,ipz,idx,idy,idz)
        call gather_buf(rnode,stot)
c
        stot = 0
        if (nolist) then
          ixgrd = 0
          do i = 1, antot
            if (ra(i,1,1)-xmin.lt.rcmax) then
              ixgrd = ixgrd + 1
              xgrd(ixgrd) = i
              stot = stot + 1
              is = btot+stot
              xcrd(is) = ra(i,1,1)
              ycrd(is) = ra(i,2,1)
              zcrd(is) = ra(i,3,1)
              bidx(is) = aidx(i)
              bat(is) = at(i)
            endif
          end do
        else
          do i = 1, ixgrd
            ig = xgrd(i)
            stot = stot + 1
            is = btot+stot
            xcrd(is) = ra(ig,1,1)
            ycrd(is) = ra(ig,2,1)
            zcrd(is) = ra(ig,3,1)
            bidx(is) = aidx(ig)
            bat(is) = at(ig)
          end do
        endif
c
c    send boundary layer particles to next buffer
c
        inx = ipx + 1
        if (inx.eq.idx) inx = 0
        call i_xyz_to_proc(rnode,inx,ipy,ipz,idx,idy,idz)
        call gather_buf(rnode,stot)
      endif
c
c    repeat for y axis
c
      if (idy.eq.2) then
c
c    both boundaries are adjacent to boundaries on next processor
c
        stot = 0
        if (nolist) then
          iygru = 0
          do i = 1, antot
            if ((ymax-ra(i,2,1).lt.rcmax.or.
     +         (ra(i,2,1)-ymin).lt.rcmax)) then
              iygru = iygru + 1
              ygru(iygru) = i
              stot = stot + 1
              is = btot+stot
              xcrd(is) = ra(i,1,1)
              ycrd(is) = ra(i,2,1)
              zcrd(is) = ra(i,3,1)
              bidx(is) = aidx(i)
              bat(is) = at(i)
            endif
          end do
        else
          do i = 1, iygru
            ig = ygru(i)
            stot = stot + 1
            is = btot+stot
            xcrd(is) = ra(ig,1,1)
            ycrd(is) = ra(ig,2,1)
            zcrd(is) = ra(ig,3,1)
            bidx(is) = aidx(ig)
            bat(is) = at(ig)
          end do
        endif
c
        if (nolist) then
          iygbu = 0
          do i = 1, btot
            if ((ymax-ycrd(i).lt.rcmax.or.
     +         (ycrd(i)-ymin).lt.rcmax)) then
              iygbu = iygbu + 1
              ygbu(iygbu) = i
              stot = stot + 1
              is = btot+stot
              xcrd(is) = xcrd(i)
              ycrd(is) = ycrd(i)
              zcrd(is) = zcrd(i)
              bidx(is) = bidx(i)
              bat(is) = bat(i)
            endif
          end do
        else
          do i = 1, iygbu
            ig = ygbu(i)
            stot = stot + 1
            is = btot+stot
            xcrd(is) = xcrd(ig)
            ycrd(is) = ycrd(ig)
            zcrd(is) = zcrd(ig)
            bidx(is) = bidx(ig)
            bat(is) = bat(ig)
          end do
        endif
c
c    send boundary layer particles to next buffer
c
        iny = ipy - 1
        if (iny.eq.-1) iny = idy - 1
        call i_xyz_to_proc(rnode,ipx,iny,ipz,idx,idy,idz)
        call gather_buf(rnode,stot)
      elseif (idy.gt.2) then
c
c   each boundary is adjacent to a different processor
c
        stot = 0
        if (nolist) then
          iygru = 0
          do i = 1, antot
            if (ymax-ra(i,2,1).lt.rcmax) then
              iygru = iygru + 1
              ygru(iygru) = i
              stot = stot + 1
              is = btot+stot
              xcrd(is) = ra(i,1,1)
              ycrd(is) = ra(i,2,1)
              zcrd(is) = ra(i,3,1)
              bidx(is) = aidx(i)
              bat(is) = at(i)
            endif
          end do
        else
          do i = 1, iygru
            ig = ygru(i)
            stot = stot + 1
            is = btot+stot
            xcrd(is) = ra(ig,1,1)
            ycrd(is) = ra(ig,2,1)
            zcrd(is) = ra(ig,3,1)
            bidx(is) = aidx(ig)
            bat(is) = at(ig)
          end do
        endif
c
        if (nolist) then
          iygbu = 0
          do i = 1, btot
            if (ymax-ycrd(i).lt.rcmax) then
              iygbu = iygbu + 1
              ygbu(iygbu) = i
              stot = stot + 1
              is = btot+stot
              xcrd(is) = xcrd(i)
              ycrd(is) = ycrd(i)
              zcrd(is) = zcrd(i)
              bidx(is) = bidx(i)
              bat(is) = bat(i)
            endif
          end do
        else
          do i = 1, iygbu
            ig = ygbu(i)
            stot = stot + 1
            is = btot+stot
            xcrd(is) = xcrd(ig)
            ycrd(is) = ycrd(ig)
            zcrd(is) = zcrd(ig)
            bidx(is) = bidx(ig)
            bat(is) = bat(ig)
          end do
        endif
c
c    send boundary layer particles to next buffer
c
        iny = ipy - 1
        if (iny.eq.-1) iny = idy - 1
        call i_xyz_to_proc(rnode,ipx,iny,ipz,idx,idy,idz)
        call gather_buf(rnode,stot)
c
        stot = 0
        if (nolist) then
          iygrd = 0
          do i = 1, antot
            if (ra(i,2,1)-ymin.lt.rcmax) then
              iygrd = iygrd + 1
              ygrd(iygrd) = i
              stot = stot + 1    
              is = btot+stot
              xcrd(is) = ra(i,1,1)
              ycrd(is) = ra(i,2,1)
              zcrd(is) = ra(i,3,1)
              bidx(is) = aidx(i) 
              bat(is) = at(i)  
            endif
          end do
        else
          do i = 1, iygrd
            ig = ygrd(i)
            stot = stot + 1    
            is = btot+stot
            xcrd(is) = ra(ig,1,1)
            ycrd(is) = ra(ig,2,1)
            zcrd(is) = ra(ig,3,1)
            bidx(is) = aidx(ig) 
            bat(is) = at(ig)
          end do
        endif
c
        if (nolist) then
          iygbd = 0
          do i = 1, btot
            if (ycrd(i)-ymin.lt.rcmax.and.ycrd(i).gt.ymin) then
              iygbd = iygbd + 1
              ygbd(iygbd) = i
              stot = stot + 1
              is = btot+stot
              xcrd(is) = xcrd(i)
              ycrd(is) = ycrd(i)
              zcrd(is) = zcrd(i)
              bidx(is) = bidx(i)
              bat(is) = bat(i)
            endif
          end do
        else
          do i = 1, iygbd
            ig = ygbd(i)
            stot = stot + 1
            is = btot+stot
            xcrd(is) = xcrd(ig)
            ycrd(is) = ycrd(ig)
            zcrd(is) = zcrd(ig)
            bidx(is) = bidx(ig)
            bat(is) = bat(ig)
          end do
        endif
c
c    send boundary layer particles to next buffer
c
        iny = ipy + 1
        if (iny.eq.idy) iny = 0
        call i_xyz_to_proc(rnode,ipx,iny,ipz,idx,idy,idz)
        call gather_buf(rnode,stot)
      endif
c
c    repeat for z axis
c
      if (idz.eq.2) then
c
c    both boundaries are adjacent to boundaries on next processor
c
        stot = 0
        if (nolist) then
          izgru = 0
          do i = 1, antot
            if ((zmax-ra(i,3,1).lt.rcmax.or.
     +         ra(i,3,1)-zmin.lt.rcmax)) then
              izgru = izgru + 1
              zgru(izgru) = i
              stot = stot + 1
              is = btot+stot
              xcrd(is) = ra(i,1,1)
              ycrd(is) = ra(i,2,1)
              zcrd(is) = ra(i,3,1)
              bidx(is) = aidx(i)
              bat(is) = at(i)
            endif
          end do
        else
          do i = 1, izgru
            ig = zgru(i)
            stot = stot + 1
            is = btot+stot
            xcrd(is) = ra(ig,1,1)
            ycrd(is) = ra(ig,2,1)
            zcrd(is) = ra(ig,3,1)
            bidx(is) = aidx(ig)
            bat(is) = at(ig)
          end do
        endif
c
        if (nolist) then
          izgbu = 0
          do i = 1, btot
            if ((zmax-zcrd(i).lt.rcmax.or.
     +         zcrd(i)-zmin.lt.rcmax)) then
              izgbu = izgbu + 1
              zgbu(izgbu) = i
              stot = stot + 1
              is = btot+stot
              xcrd(is) = xcrd(i)
              ycrd(is) = ycrd(i)
              zcrd(is) = zcrd(i)
              bidx(is) = bidx(i)
              bat(is) = bat(i)
              bidx(is) = bidx(i)
              bat(is) = bat(i)
            endif
          end do
        else
          do i = 1, izgbu
            ig = zgbu(i)
            stot = stot + 1
            is = btot+stot
            xcrd(is) = xcrd(ig)
            ycrd(is) = ycrd(ig)
            zcrd(is) = zcrd(ig)
            bidx(is) = bidx(ig)
            bat(is) = bat(ig)
          end do
        endif
c
c    send boundary layer particles to next buffer
c
        inz = ipz - 1
        if (inz.eq.-1) inz = idz - 1
        call i_xyz_to_proc(rnode,ipx,ipy,inz,idx,idy,idz)
        call gather_buf(rnode,stot)
      elseif (idz.gt.2) then
c
c   each boundary is adjacent to a different processor
c
        stot = 0
        if (nolist) then
          izgru = 0
          do i = 1, antot
            if (zmax-ra(i,3,1).lt.rcmax) then
              izgru = izgru + 1
              zgru(izgru) = i
              stot = stot + 1
              is = btot+stot
              xcrd(is) = ra(i,1,1)
              ycrd(is) = ra(i,2,1)
              zcrd(is) = ra(i,3,1)
              bidx(is) = aidx(i)
              bat(is) = at(i)
            endif
          end do
        else
          do i = 1, izgru
            ig = zgru(i)
            stot = stot + 1
            is = btot+stot
            xcrd(is) = ra(ig,1,1)
            ycrd(is) = ra(ig,2,1)
            zcrd(is) = ra(ig,3,1)
            bidx(is) = aidx(ig)
            bat(is) = at(ig)
          end do
        endif
c
        if (nolist) then
          izgbu = 0
          do i = 1, btot
            if (zmax-zcrd(i).lt.rcmax) then
              izgbu = izgbu + 1
              zgbu(izgbu) = i
              stot = stot + 1
              is = btot+stot
              xcrd(is) = xcrd(i)
              ycrd(is) = ycrd(i)
              zcrd(is) = zcrd(i)
              bidx(is) = bidx(i)
              bat(is) = bat(i)
            endif
          end do
        else
          do i = 1, izgbu
            ig = zgbu(i)
            stot = stot + 1
            is = btot+stot
            xcrd(is) = xcrd(ig)
            ycrd(is) = ycrd(ig)
            zcrd(is) = zcrd(ig)
            bidx(is) = bidx(ig)
            bat(is) = bat(ig)
          end do
        endif
c
c    send boundary layer particles to next buffer
c
        inz = ipz - 1
        if (inz.eq.-1) inz = idz - 1
        call i_xyz_to_proc(rnode,ipx,ipy,inz,idx,idy,idz)
        call gather_buf(rnode,stot)
c
        stot = 0
        if (nolist) then
          izgrd = 0
          do i = 1, antot
            if (ra(i,3,1)-zmin.lt.rcmax) then
              izgrd = izgrd + 1
              zgrd(izgrd) = i
              stot = stot + 1
              is = btot+stot
              xcrd(is) = ra(i,1,1)
              ycrd(is) = ra(i,2,1)
              zcrd(is) = ra(i,3,1)
              bidx(is) = aidx(i)
              bat(is) = at(i)
            endif
          end do
        else
          do i = 1, izgrd
            ig = zgrd(i)
            stot = stot + 1
            is = btot+stot
            xcrd(is) = ra(ig,1,1)
            ycrd(is) = ra(ig,2,1)
            zcrd(is) = ra(ig,3,1)
            bidx(is) = aidx(ig)
            bat(is) = at(ig)
          end do
        endif
c
        if (nolist) then
          izgbd = 0
          do i = 1, btot
            if (zcrd(i)-zmin.lt.rcmax.and.zcrd(i).gt.zmin) then
              izgbd = izgbd + 1
              zgbd(izgbd) = i
              stot = stot + 1
              is = btot+stot
              xcrd(is) = xcrd(i)
              ycrd(is) = ycrd(i)
              zcrd(is) = zcrd(i)
              bidx(is) = bidx(i)
              bat(is) = bat(i)
            endif
          end do
        else
          do i = 1, izgbd
            ig = zgbd(i)
            stot = stot + 1
            is = btot+stot
            xcrd(is) = xcrd(ig)
            ycrd(is) = ycrd(ig)
            zcrd(is) = zcrd(ig)
            bidx(is) = bidx(ig)
            bat(is) = bat(ig)
          end do
        endif
c
c    send boundary layer particles to next buffer
c
        inz = ipz + 1
        if (inz.eq.idz) inz = 0
        call i_xyz_to_proc(rnode,ipx,ipy,inz,idx,idy,idz)
        call gather_buf(rnode,stot)
      endif
c      call heapsort(1)
cdbg      do i = 1, antot
cdbg        write(6,178) inode,i,aidx(i)
cdbg      end do
cdbg  178 format(i3,' Gather aidx(',i3,'): ',i8)
cdbg      do i = 1, btot
cdbg        write(6,179) inode,i,bidx(i)
cdbg      end do
cdbg  179 format(i3,' Gather bidx(',i3,'): ',i8)
c
      tmstat(6) = tmstat(6) + wraptime() - tbeg
c
      return
      end
c
      subroutine gather_buf(rnode,stot)
#include "common.fh"
c
      double precision buf(3,MAXAT)
      integer ibuf(2,MAXAT)
      double precision tbeg,wraptime
      integer rnode,i
      integer me,one,ld2,ld3,is,rtot,stot
c
c   Get data from processor rnode. The amount of data on the
c   current processor that actually needs to be accessed
c   represents stot particles.
c
      me = ga_nodeid()
      one = 1
      ld2 = 2
      ld3 = 3
c
c   exchange the size of lists
c
      tbeg = wraptime()
c
      do i = btot+1, btot+stot
        is = i - btot
        buf(1,is) = xcrd(i)
        buf(2,is) = ycrd(i)
        buf(3,is) = zcrd(i)
        ibuf(1,is) = bidx(i)
        ibuf(2,is) = bat(i)
      end do
c
      gsize_hi = gsize_lo(me)
      call nga_put(g_size,gsize_lo(me),gsize_hi,stot,one)
      gcoords_hi(1) = 3
      gcoords_hi(2) = gcoords_lo(2,me) + stot - 1
      if (stot.gt.0) call nga_put(g_coords,gcoords_lo(1,me),
     +                            gcoords_hi,buf,ld3)
      giat_hi(1) = 2
      giat_hi(2) = giat_lo(2,me) + stot - 1
      if (stot.gt.0) call nga_put(g_iat,giat_lo(1,me),giat_hi,ibuf,ld2)

      call ga_sync()

      gsize_hi = gsize_lo(rnode)
      call nga_get(g_size,gsize_lo(rnode),gsize_hi,rtot,one)
      gcoords_hi(2) = gcoords_lo(2,rnode) + rtot - 1
      if (rtot.gt.0) call nga_get(g_coords,gcoords_lo(1,rnode),
     +                            gcoords_hi,buf,ld3)
      giat_hi(2) = giat_lo(2,rnode) + rtot - 1
      if (rtot.gt.0) call nga_get(g_iat,giat_lo(1,rnode),
     +                            giat_hi,ibuf,ld2)
c
      do i = btot+1, btot+rtot
        is = i - btot
        xcrd(i) = buf(1,is)
        ycrd(i) = buf(2,is)
        zcrd(i) = buf(3,is)
        bidx(i) = ibuf(1,is)
        bat(i) = ibuf(2,is)
      end do
c
      call ga_sync()
      btot = btot+rtot
c
      tmstat(11) = tmstat(11) + wraptime() - tbeg
c
      return
      end
