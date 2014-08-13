#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine scatter
#include "common.fh"
c
      double precision xmax,xmin,ymax,ymin,zmax,zmin
      double precision xcell,ycell,zcell
      double precision tbeg,wraptime
      integer snode,rnode,pnum,inode,idx,idy,idz,ipx,ipy,ipz
      integer stot,i,inx,iny,inz
      logical nolist
c
c   This subroutine scatters the forces of all particles on
c   neighboring processors back to their home processors after the
c   forces have been calculated.
c
      tbeg = wraptime()
      pnum = ga_nnodes()
      call factor(pnum,idx,idy,idz)
      inode = ga_nodeid()
      call i_proc_to_xyz(inode,ipx,ipy,ipz,idx,idy,idz)
c
      if (istep.eq.0.or.(mod(istep,ilist).eq.0.and.
     +    t_rmndr.eq.0.0d00)) then
        nolist = .true.
      else
        nolist = .false.
      endif
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
      xcell = xbox/dble(idx)
      ycell = ybox/dble(idy)
      zcell = zbox/dble(idz)
c
c    Send particles along x-axis.
c
c
      if (idx.gt.1) then
        if (idx.eq.2) then
          stot = 0
          savtot = 0
          if (nolist) then
            do i = 1, btot
              if (xcrd(i).lt.xmin.or.xcrd(i).ge.xmax) then
                xsbu(i) = .true.
                stot = stot + 1
                xcrd(btot+stot) = xcrd(i)
                ycrd(btot+stot) = ycrd(i)
                zcrd(btot+stot) = zcrd(i)
                xfrc(btot+stot) = xfrc(i)
                yfrc(btot+stot) = yfrc(i)
                zfrc(btot+stot) = zfrc(i)
                bidx(btot+stot) = bidx(i)
              else
                xsbu(i) = .false.
                savtot = savtot + 1
                xcrd(savtot) = xcrd(i)
                ycrd(savtot) = ycrd(i)
                zcrd(savtot) = zcrd(i)
                xfrc(savtot) = xfrc(i)
                yfrc(savtot) = yfrc(i)
                zfrc(savtot) = zfrc(i)
                bidx(savtot) = bidx(i)
              endif
            end do
          else
            do i = 1, btot
              if (xsbu(i)) then
                stot = stot + 1
                xfrc(btot+stot) = xfrc(i)
                yfrc(btot+stot) = yfrc(i)
                zfrc(btot+stot) = zfrc(i)
                bidx(btot+stot) = bidx(i)
              else
                savtot = savtot + 1
                xfrc(savtot) = xfrc(i)
                yfrc(savtot) = yfrc(i)
                zfrc(savtot) = zfrc(i)
                bidx(savtot) = bidx(i)
              endif
            end do
          endif
        else
          stot = 0
          savtot = 0
          if (nolist) then
            do i = 1, btot
              if ((xcrd(i).ge.xmax.and.xcrd(i)-xmax.lt.xcell).or.
     +            (xcrd(i).lt.xmin.and.xmin-xcrd(i).gt.xcell)) then
                xsbu(i) = .true.
                stot = stot + 1
                xcrd(btot+stot) = xcrd(i)
                ycrd(btot+stot) = ycrd(i)
                zcrd(btot+stot) = zcrd(i)
                xfrc(btot+stot) = xfrc(i)
                yfrc(btot+stot) = yfrc(i)
                zfrc(btot+stot) = zfrc(i)
                bidx(btot+stot) = bidx(i)
              else
                xsbu(i) = .false.
                savtot = savtot + 1
                xcrd(savtot) = xcrd(i)
                ycrd(savtot) = ycrd(i)
                zcrd(savtot) = zcrd(i)
                xfrc(savtot) = xfrc(i)
                yfrc(savtot) = yfrc(i)
                zfrc(savtot) = zfrc(i)
                bidx(savtot) = bidx(i)
              endif
            end do
          else
            do i = 1, btot
              if (xsbu(i)) then
                stot = stot + 1
                xfrc(btot+stot) = xfrc(i)
                yfrc(btot+stot) = yfrc(i)
                zfrc(btot+stot) = zfrc(i)
                bidx(btot+stot) = bidx(i)
              else
                savtot = savtot + 1
                xfrc(savtot) = xfrc(i)
                yfrc(savtot) = yfrc(i)
                zfrc(savtot) = zfrc(i)
                bidx(savtot) = bidx(i)
              endif
            end do
          endif
        endif
c
        inx = ipx - 1
        if (inx.eq.-1) inx = idx - 1
        call i_xyz_to_proc(rnode,inx,ipy,ipz,idx,idy,idz)
        call scatter_buf(rnode,stot)
c        call heapsort(1)
        call fcull(1)
      endif
c
      if (idx.gt.2) then
c
        stot = 0
        savtot = 0
        if (nolist) then
          do i = 1, btot
            if ((xcrd(i).lt.xmin.and.xmin-xcrd(i).lt.xcell).or.
     +          (xcrd(i).gt.xmax.and.xcrd(i)-xmax.gt.xcell)) then
              xsbd(i) = .true.
              stot = stot + 1
              xcrd(btot+stot) = xcrd(i)
              ycrd(btot+stot) = ycrd(i)
              zcrd(btot+stot) = zcrd(i)
              xfrc(btot+stot) = xfrc(i)
              yfrc(btot+stot) = yfrc(i)
              zfrc(btot+stot) = zfrc(i)
              bidx(btot+stot) = bidx(i)
            else
              xsbd(i) = .false.
              savtot = savtot + 1
              xcrd(savtot) = xcrd(i)
              ycrd(savtot) = ycrd(i)
              zcrd(savtot) = zcrd(i)
              xfrc(savtot) = xfrc(i)
              yfrc(savtot) = yfrc(i)
              zfrc(savtot) = zfrc(i)
              bidx(savtot) = bidx(i)
            endif
          end do
        else
          do i = 1, btot
            if (xsbd(i)) then
              stot = stot + 1
              xfrc(btot+stot) = xfrc(i)
              yfrc(btot+stot) = yfrc(i)
              zfrc(btot+stot) = zfrc(i)
              bidx(btot+stot) = bidx(i)
            else
              savtot = savtot + 1
              xfrc(savtot) = xfrc(i)
              yfrc(savtot) = yfrc(i)
              zfrc(savtot) = zfrc(i)
              bidx(savtot) = bidx(i)
            endif
          end do
        endif
c
        inx = ipx + 1
        if (inx.eq.idx) inx = 0
        call i_xyz_to_proc(rnode,inx,ipy,ipz,idx,idy,idz)
        call scatter_buf(rnode,stot)
c        call heapsort(1)
        call fcull(2)
c
      endif
c
c    Send particles along y-axis.
c
      if (idy.gt.1) then
c
        if (idy.eq.2) then
          stot = 0
          savtot = 0
          if (nolist) then
            do i = 1, btot
              if (ycrd(i).lt.ymin.or.ycrd(i).ge.ymax) then
                ysbu(i) = .true.
                stot = stot + 1
                xcrd(btot+stot) = xcrd(i)
                ycrd(btot+stot) = ycrd(i)
                zcrd(btot+stot) = zcrd(i)
                xfrc(btot+stot) = xfrc(i)
                yfrc(btot+stot) = yfrc(i)
                zfrc(btot+stot) = zfrc(i)
                bidx(btot+stot) = bidx(i)
              else
                ysbu(i) = .false.
                savtot = savtot + 1
                xcrd(savtot) = xcrd(i)
                ycrd(savtot) = ycrd(i)
                zcrd(savtot) = zcrd(i)
                xfrc(savtot) = xfrc(i)
                yfrc(savtot) = yfrc(i)
                zfrc(savtot) = zfrc(i)
                bidx(savtot) = bidx(i)
              endif
            end do
          else
            do i = 1, btot
              if (ysbu(i)) then
                stot = stot + 1
                xfrc(btot+stot) = xfrc(i)
                yfrc(btot+stot) = yfrc(i)
                zfrc(btot+stot) = zfrc(i)
                bidx(btot+stot) = bidx(i)
              else
                savtot = savtot + 1
                xfrc(savtot) = xfrc(i)
                yfrc(savtot) = yfrc(i)
                zfrc(savtot) = zfrc(i)
                bidx(savtot) = bidx(i)
              endif
            end do
          endif
        else
          stot = 0
          savtot = 0
          if (nolist) then
            do i = 1, btot
              if ((ycrd(i).ge.ymax.and.ycrd(i)-ymax.lt.ycell).or.
     +            (ycrd(i).lt.ymin.and.ymin-ycrd(i).gt.ycell)) then
                ysbu(i) = .true.
                stot = stot + 1
                xcrd(btot+stot) = xcrd(i)
                ycrd(btot+stot) = ycrd(i)
                zcrd(btot+stot) = zcrd(i)
                xfrc(btot+stot) = xfrc(i)
                yfrc(btot+stot) = yfrc(i)
                zfrc(btot+stot) = zfrc(i)
                bidx(btot+stot) = bidx(i)
              else
                ysbu(i) = .false.
                savtot = savtot + 1
                xcrd(savtot) = xcrd(i)
                ycrd(savtot) = ycrd(i)
                zcrd(savtot) = zcrd(i)
                xfrc(savtot) = xfrc(i)
                yfrc(savtot) = yfrc(i)
                zfrc(savtot) = zfrc(i)
                bidx(savtot) = bidx(i)
              endif
            end do
          else
            do i = 1, btot
              if (ysbu(i)) then
                stot = stot + 1
                xfrc(btot+stot) = xfrc(i)
                yfrc(btot+stot) = yfrc(i)
                zfrc(btot+stot) = zfrc(i)
                bidx(btot+stot) = bidx(i)
              else
                savtot = savtot + 1
                xfrc(savtot) = xfrc(i)
                yfrc(savtot) = yfrc(i)
                zfrc(savtot) = zfrc(i)
                bidx(savtot) = bidx(i)
              endif
            end do
          endif
        endif
c
        iny = ipy - 1
        if (iny.eq.-1) iny = idy - 1
        call i_xyz_to_proc(rnode,ipx,iny,ipz,idx,idy,idz)
        call scatter_buf(rnode,stot)
c        call heapsort(1)
        call fcull(3)
c
      endif
c
      if (idy.gt.2) then
c
        stot = 0
        savtot = 0
        if (nolist) then
          do i = 1, btot
            if ((ycrd(i).lt.ymin.and.ymin-ycrd(i).lt.ycell).or.
     +          (ycrd(i).gt.ymax.and.ycrd(i)-ymax.gt.ycell)) then
              ysbd(i) = .true.
              stot = stot + 1
              xcrd(btot+stot) = xcrd(i)
              ycrd(btot+stot) = ycrd(i)
              zcrd(btot+stot) = zcrd(i)
              xfrc(btot+stot) = xfrc(i)
              yfrc(btot+stot) = yfrc(i)
              zfrc(btot+stot) = zfrc(i)
              bidx(btot+stot) = bidx(i)
            else
              ysbd(i) = .false.
              savtot = savtot + 1
              xcrd(savtot) = xcrd(i)
              ycrd(savtot) = ycrd(i)
              zcrd(savtot) = zcrd(i)
              xfrc(savtot) = xfrc(i)
              yfrc(savtot) = yfrc(i)
              zfrc(savtot) = zfrc(i)
              bidx(savtot) = bidx(i)
            endif
          end do
        else
          do i = 1, btot
            if (ysbd(i)) then
              stot = stot + 1
              xfrc(btot+stot) = xfrc(i)
              yfrc(btot+stot) = yfrc(i)
              zfrc(btot+stot) = zfrc(i)
              bidx(btot+stot) = bidx(i)
            else
              savtot = savtot + 1
              xfrc(savtot) = xfrc(i)
              yfrc(savtot) = yfrc(i)
              zfrc(savtot) = zfrc(i)
              bidx(savtot) = bidx(i)
            endif
          end do
        endif
c
        iny = ipy + 1
        if (iny.eq.idy) iny = 0
        call i_xyz_to_proc(rnode,ipx,iny,ipz,idx,idy,idz)
        call scatter_buf(rnode,stot)
c        call heapsort(1)
        call fcull(4)
c
      endif
c
c   send particles along z-axis
c
      if (idz.gt.1) then
c
        if (idz.eq.2) then
          stot = 0
          savtot = 0
          if (nolist) then
            do i = 1, btot
              if (zcrd(i).lt.zmin.or.zcrd(i).ge.zmax) then
                zsbu(i) = .true.
                stot = stot + 1
                xcrd(btot+stot) = xcrd(i)
                ycrd(btot+stot) = ycrd(i)
                zcrd(btot+stot) = zcrd(i)
                xfrc(btot+stot) = xfrc(i)
                yfrc(btot+stot) = yfrc(i)
                zfrc(btot+stot) = zfrc(i)
                bidx(btot+stot) = bidx(i)
              else
                zsbu(i) = .false.
                savtot = savtot + 1
                xcrd(savtot) = xcrd(i)
                ycrd(savtot) = ycrd(i)
                zcrd(savtot) = zcrd(i)
                xfrc(savtot) = xfrc(i)
                yfrc(savtot) = yfrc(i)
                zfrc(savtot) = zfrc(i)
                bidx(savtot) = bidx(i)
              endif
            end do
          else
            do i = 1, btot
              if (zsbu(i)) then
                stot = stot + 1
                xfrc(btot+stot) = xfrc(i)
                yfrc(btot+stot) = yfrc(i)
                zfrc(btot+stot) = zfrc(i)
                bidx(btot+stot) = bidx(i)
              else
                savtot = savtot + 1
                xfrc(savtot) = xfrc(i)
                yfrc(savtot) = yfrc(i)
                zfrc(savtot) = zfrc(i)
                bidx(savtot) = bidx(i)
              endif
            end do
          endif
        else
          stot = 0
          savtot = 0
          if (nolist) then
            do i = 1, btot
              if ((zcrd(i).ge.zmax.and.zcrd(i)-zmax.lt.zcell).or.
     +            (zcrd(i).lt.zmin.and.zmin-zcrd(i).gt.zcell)) then
                zsbu(i) = .true.
                stot = stot + 1
                xcrd(btot+stot) = xcrd(i)
                ycrd(btot+stot) = ycrd(i)
                zcrd(btot+stot) = zcrd(i)
                xfrc(btot+stot) = xfrc(i)
                yfrc(btot+stot) = yfrc(i)
                zfrc(btot+stot) = zfrc(i)
                bidx(btot+stot) = bidx(i)
              else
                zsbu(i) = .false.
                savtot = savtot + 1
                xcrd(savtot) = xcrd(i)
                ycrd(savtot) = ycrd(i)
                zcrd(savtot) = zcrd(i)
                xfrc(savtot) = xfrc(i)
                yfrc(savtot) = yfrc(i)
                zfrc(savtot) = zfrc(i)
                bidx(savtot) = bidx(i)
              endif
            end do
          else
            do i = 1, btot
              if (zsbu(i)) then
                stot = stot + 1
                xfrc(btot+stot) = xfrc(i)
                yfrc(btot+stot) = yfrc(i)
                zfrc(btot+stot) = zfrc(i)
                bidx(btot+stot) = bidx(i)
              else
                savtot = savtot + 1
                xfrc(savtot) = xfrc(i)
                yfrc(savtot) = yfrc(i)
                zfrc(savtot) = zfrc(i)
                bidx(savtot) = bidx(i)
              endif
            end do
          endif
        endif
c
        inz = ipz - 1
        if (inz.eq.-1) inz = idz - 1
        call i_xyz_to_proc(rnode,ipx,ipy,inz,idx,idy,idz)
        call scatter_buf(rnode,stot)
c        call heapsort(1)
        call fcull(5)
c
      endif
c
      if (idz.gt.2) then 
c
        stot = 0
        savtot = 0
        if (nolist) then
          do i = 1, btot
            if ((zcrd(i).lt.zmin.and.zmin-zcrd(i).lt.zcell).or.
     +          (zcrd(i).gt.zmax.and.zcrd(i)-zmax.gt.zcell)) then
              zsbd(i) = .true.
              stot = stot + 1
              xcrd(btot+stot) = xcrd(i)
              ycrd(btot+stot) = ycrd(i)
              zcrd(btot+stot) = zcrd(i)
              xfrc(btot+stot) = xfrc(i)
              yfrc(btot+stot) = yfrc(i)
              zfrc(btot+stot) = zfrc(i)
              bidx(btot+stot) = bidx(i)
            else
              zsbd(i) = .false.
              savtot = savtot + 1
              xcrd(savtot) = xcrd(i)
              ycrd(savtot) = ycrd(i)
              zcrd(savtot) = zcrd(i)
              xfrc(savtot) = xfrc(i)
              yfrc(savtot) = yfrc(i)
              zfrc(savtot) = zfrc(i)
              bidx(savtot) = bidx(i)
            endif
          end do
        else
          do i = 1, btot
            if (zsbd(i)) then
              stot = stot + 1
              xfrc(btot+stot) = xfrc(i)
              yfrc(btot+stot) = yfrc(i)
              zfrc(btot+stot) = zfrc(i)
              bidx(btot+stot) = bidx(i)
            else
              savtot = savtot + 1
              xfrc(savtot) = xfrc(i)
              yfrc(savtot) = yfrc(i)
              zfrc(savtot) = zfrc(i)
              bidx(savtot) = bidx(i)
            endif
          end do
        endif
c
        inz = ipz + 1 
        if (inz.eq.idz) inz = 0
        call i_xyz_to_proc(rnode,ipx,ipy,inz,idx,idy,idz) 
        call scatter_buf(rnode,stot) 
c        call heapsort(1)
        call fcull(6)
c 
      endif
c
      tmstat(7) = tmstat(7) + wraptime() - tbeg
c
      return
      end
c
      subroutine scatter_buf(rnode,stot)
#include "common.fh"
c
      double precision cbuf(3,MAXAT),fbuf(3,MAXAT)
      integer ibuf(MAXAT)
      double precision tbeg, tbeg1, wraptime
      integer i,rnode
      integer me,one,ld,is,rtot,stot
      logical crdflag
c
c   Get the forces from processor rnode. The amount of data
c   on the current processor that actually needs to be
c   accessed represents stot particles
c
      me = ga_nodeid()
      one = 1
      ld = 3
c
c   exchange the size of lists
c
      tbeg = wraptime()
c
      if (istep.eq.0.or.(mod(istep,ilist).eq.0.and.
     +    t_rmndr.eq.0.0d00)) then
        crdflag = .true.
      else
        crdflag = .false.
      endif
c
      if (crdflag) then
        do i = btot + 1, btot + stot
          is = i - btot
          cbuf(1,is) = xcrd(i)
          cbuf(2,is) = ycrd(i)
          cbuf(3,is) = zcrd(i)
        end do
      endif
      do i = btot + 1, btot + stot
        is = i - btot
        fbuf(1,is) = xfrc(i)
        fbuf(2,is) = yfrc(i)
        fbuf(3,is) = zfrc(i)
        ibuf(is) = bidx(i)
      end do
c
      gsize_hi = gsize_lo(me)
      call nga_put(g_size,gsize_lo(me),gsize_hi,stot,one)
      if (crdflag) then
        gcoords_hi(1) = 3
        gcoords_hi(2) = gcoords_lo(2,me) + stot - 1
        if (stot.gt.0) call nga_put(g_coords,gcoords_lo(1,me),
     +                              gcoords_hi,cbuf,ld)
      endif
      gfrc_hi(1) = 3
      gfrc_hi(2) = gfrc_lo(2,me) + stot - 1
      if (stot.gt.0) call nga_put(g_frc,gfrc_lo(1,me),gfrc_hi,fbuf,ld)
      gindex_hi = gindex_lo(me) + stot - 1
      if (stot.gt.0) call nga_put(g_index,gindex_lo(me),
     +                            gindex_hi,ibuf,ld)

      tbeg1 = wraptime()
      call ga_sync()
      tmstat(16) = tmstat(16) + wraptime() - tbeg1
      btot = savtot

      gsize_hi = gsize_lo(rnode)
      call nga_get(g_size,gsize_lo(rnode),gsize_hi,rtot,one)
      if (crdflag) then
        gcoords_hi(2) = gcoords_lo(2,rnode) + rtot - 1
        if (rtot.gt.0) call nga_get(g_coords,gcoords_lo(1,rnode),
     +                              gcoords_hi,cbuf,ld)
      endif
      gfrc_hi(2) = gfrc_lo(2,rnode) + rtot - 1
      if (rtot.gt.0) call nga_get(g_frc,gfrc_lo(1,rnode),
     +                            gfrc_hi,fbuf,ld)
      gindex_hi = gindex_lo(rnode) + rtot - 1
      if (rtot.gt.0) call nga_get(g_index,gindex_lo(rnode),
     +                            gindex_hi,ibuf,ld)
c
      if (crdflag) then
        do i = btot+1, btot+rtot
          is = i - btot
          xcrd(i) = cbuf(1,is)
          ycrd(i) = cbuf(2,is)
          zcrd(i) = cbuf(3,is)
        end do
      endif
      do i = btot+1, btot+rtot
        is = i - btot
        xfrc(i) = fbuf(1,is)
        yfrc(i) = fbuf(2,is)
        zfrc(i) = fbuf(3,is)
        bidx(i) = ibuf(is)
      end do
c
      call ga_sync()
      btot = btot+rtot
      if (btot.gt.MAXAT) then
        call ga_error("Array overflow in scatter_buf",btot)
      endif
c
      tmstat(12) = tmstat(12) + wraptime() - tbeg
c
      return
      end
c
      subroutine fcull(iflg)
#include "common.fh"
c
      double precision xmax,xmin,ymax,ymin,zmax,zmin
      integer pnum,inode,idx,idy,idz,ipx,ipy,ipz
      integer i,is,icnt,iflg,get_hash_a
      logical nolist
c
c   This subroutine culls the forces of particles in the buffers
c   back on to the particles on the home processor
c
      pnum = ga_nnodes()
      call factor(pnum,idx,idy,idz)
      inode = ga_nodeid()
      call i_proc_to_xyz(inode,ipx,ipy,ipz,idx,idy,idz)
c
      if (istep.eq.0.or.(mod(istep,ilist).eq.0.and.
     +    t_rmndr.eq.0.0d00)) then
        nolist = .true.
      else
        nolist = .false.
      endif
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
      icnt = 0
      is = 1
      if (nolist) then
        do i = 1, btot
          if ((xcrd(i).lt.xmax.and.xcrd(i).ge.xmin).and.
     +        (ycrd(i).lt.ymax.and.ycrd(i).ge.ymin).and.
     +        (zcrd(i).lt.zmax.and.zcrd(i).ge.zmin)) then
            clist(iflg,i) = .true.
            is = get_hash_a(bidx(i))
            if (is.gt.0) then
              ra(is,1,4) = ra(is,1,4) + xfrc(i)
              ra(is,2,4) = ra(is,2,4) + yfrc(i)
              ra(is,3,4) = ra(is,3,4) + zfrc(i)
            else
              call ga_error("Unknown atom index: ",i)
            endif
          else
            clist(iflg,i) = .false.
            icnt = icnt + 1
            xcrd(icnt) = xcrd(i)
            ycrd(icnt) = ycrd(i)
            zcrd(icnt) = zcrd(i)
            xfrc(icnt) = xfrc(i)
            yfrc(icnt) = yfrc(i)
            zfrc(icnt) = zfrc(i)
            bidx(icnt) = bidx(i)
          endif
        end do
      else
        do i = 1, btot
          if (clist(iflg,i)) then
            is = get_hash_a(bidx(i))
            if (is.gt.0) then
              ra(is,1,4) = ra(is,1,4) + xfrc(i)
              ra(is,2,4) = ra(is,2,4) + yfrc(i)
              ra(is,3,4) = ra(is,3,4) + zfrc(i)
            else
              call ga_error("Unknown atom index: ",i)
            endif
          else
            icnt = icnt + 1
            xfrc(icnt) = xfrc(i)
            yfrc(icnt) = yfrc(i)
            zfrc(icnt) = zfrc(i)
            bidx(icnt) = bidx(i)
          endif
        end do
      endif
      btot = icnt
      if (btot.gt.MAXAT) then
        call ga_error("Array overflow in fcull",btot)
      endif
c
      return
      end
