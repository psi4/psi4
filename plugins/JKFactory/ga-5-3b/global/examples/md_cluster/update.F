#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine update
#include "common.fh"
c
      double precision xmax,ymax,zmax,xmin,ymin,zmin
      double precision xcell,ycell,zcell
      double precision tbeg,wraptime
      integer snode,rnode,pnum,inode,idx,idy,idz,ipx,ipy,ipz
      integer icnt,i,j,k,inx,iny,inz
c
c   This subroutine gathers coordinates of particles that have moved off
c   the processor and sends then to neighboring processors. It also recieves
c   the coordinates of particles that have moved onto the processor.
c
      tbeg = wraptime()
      call clear_hash
      pnum = ga_nnodes()
      call factor(pnum,idx,idy,idz)
      inode = ga_nodeid()
      call i_proc_to_xyz(inode,ipx,ipy,ipz,idx,idy,idz)
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
      btot = 0
c
c    Send particles along x-axis. Gather all particles that lie within
c    processor domain
c
      if (idx.ge.2) then
c
c    collect all particles that have moved off of processor
c
        if (idx.eq.2) then
c
c    both boundaries are adjacent to boundaries on next processor
c
          icnt = 0
          do i = 1, antot
            if (ra(i,1,1).ge.xmax.or.ra(i,1,1).lt.xmin) then
              btot = btot + 1
              xcrd(btot) = ra(i,1,6)
              ycrd(btot) = ra(i,2,6)
              zcrd(btot) = ra(i,3,6)
              xfrc(btot) = ra(i,1,2)
              yfrc(btot) = ra(i,2,2)
              zfrc(btot) = ra(i,3,2)
              xacc(btot) = ra(i,1,3)
              yacc(btot) = ra(i,2,3)
              zacc(btot) = ra(i,3,3)
              bat(btot) = at(i)    
              mbuf(btot) = mass(i)
              bidx(btot) = aidx(i)
            else
              icnt = icnt + 1
              do j = 1, 3   
                do k = 1, 6
                  ra(icnt,j,k) = ra(i,j,k)
                end do
              end do 
              at(icnt) = at(i)
              mass(icnt) = mass(i)
              aidx(icnt) = aidx(i)
            endif    
          end do
          antot = icnt
        else
c
c    processor bounded by two separate processors
c
          icnt = 0 
          do i = 1, antot 
            if ((ra(i,1,1).ge.xmax.and.ra(i,1,1)-xmax.lt.xcell).or.
     +          (ra(i,1,1).lt.xmin.and.xmin-ra(i,1,1).gt.xcell)) then
              btot = btot + 1
              xcrd(btot) = ra(i,1,6)
              ycrd(btot) = ra(i,2,6)
              zcrd(btot) = ra(i,3,6)
              xfrc(btot) = ra(i,1,2)
              yfrc(btot) = ra(i,2,2)
              zfrc(btot) = ra(i,3,2)
              xacc(btot) = ra(i,1,3)
              yacc(btot) = ra(i,2,3)
              zacc(btot) = ra(i,3,3)
              bat(btot) = at(i)
              mbuf(btot) = mass(i)
              bidx(btot) = aidx(i)
            else
              icnt = icnt + 1
              do j = 1, 3
                do k = 1, 6
                  ra(icnt,j,k) = ra(i,j,k)
                end do
              end do
              at(icnt) = at(i)
              mass(icnt) = mass(i)
              aidx(icnt) = aidx(i)
            endif
          end do
          antot = icnt
        endif
c
c    send particles in positive x direction
c
        inx = ipx - 1
        if (inx.eq.-1) inx = idx - 1
        call i_xyz_to_proc(rnode,inx,ipy,ipz,idx,idy,idz)
        call update_buf(rnode)
        call ucull
      endif
c
      if (idx.gt.2) then
c
c    each boundary is adjacent to a different processor
c
        icnt = 0  
        do i = 1, antot  
          if ((ra(i,1,1).lt.xmin.and.xmin-ra(i,1,1).lt.xcell).or.
     +        (ra(i,1,1).gt.xmax.and.ra(i,1,1)-xmax.gt.xcell)) then
            btot = btot + 1 
            xcrd(btot) = ra(i,1,6) 
            ycrd(btot) = ra(i,2,6) 
            zcrd(btot) = ra(i,3,6) 
            xfrc(btot) = ra(i,1,2) 
            yfrc(btot) = ra(i,2,2) 
            zfrc(btot) = ra(i,3,2) 
            xacc(btot) = ra(i,1,3) 
            yacc(btot) = ra(i,2,3) 
            zacc(btot) = ra(i,3,3)
            bat(btot) = at(i)
            mbuf(btot) = mass(i)
            bidx(btot) = aidx(i)
          else
            icnt = icnt + 1
            do j = 1, 3
              do k = 1, 6
                ra(icnt,j,k) = ra(i,j,k)
              end do
            end do
            at(icnt) = at(i)
            mass(icnt) = mass(i)
            aidx(icnt) = aidx(i)
          endif
        end do
        antot = icnt
c
c    send particles in negative x direction
c
        inx = ipx + 1
        if (inx.eq.idx) inx = 0
        call i_xyz_to_proc(rnode,inx,ipy,ipz,idx,idy,idz)
        call update_buf(rnode)
        call ucull
      endif
c
c    repeat for y axis
c
      if (idy.ge.2) then
c
c    both boundaries are adjacent to boundaries on next processor
c
        if (idy.eq.2) then
          icnt = 0
          do i = 1, antot
            if (ra(i,2,1).ge.ymax.or.ra(i,2,1).lt.ymin) then
              btot = btot + 1
              xcrd(btot) = ra(i,1,6)
              ycrd(btot) = ra(i,2,6)
              zcrd(btot) = ra(i,3,6)
              xfrc(btot) = ra(i,1,2)
              yfrc(btot) = ra(i,2,2)
              zfrc(btot) = ra(i,3,2)
              xacc(btot) = ra(i,1,3)
              yacc(btot) = ra(i,2,3)
              zacc(btot) = ra(i,3,3)
              bat(btot) = at(i)
              mbuf(btot) = mass(i)
              bidx(btot) = aidx(i)
            else
              icnt = icnt + 1
              do j = 1, 3
                do k = 1, 6
                  ra(icnt,j,k) = ra(i,j,k)
                end do
              end do
              at(icnt) = at(i)
              mass(icnt) = mass(i)
              aidx(icnt) = aidx(i)
            endif
          end do
          antot = icnt
        else
c
c    processor bounded by two separate processors
c
          icnt = 0 
          do i = 1, antot 
            if ((ra(i,2,1).ge.ymax.and.ra(i,2,1)-ymax.lt.ycell).or.
     +          (ra(i,2,1).lt.ymin.and.ymin-ra(i,2,1).gt.ycell)) then
              btot = btot + 1
              xcrd(btot) = ra(i,1,6)
              ycrd(btot) = ra(i,2,6)
              zcrd(btot) = ra(i,3,6)
              xfrc(btot) = ra(i,1,2)
              yfrc(btot) = ra(i,2,2)
              zfrc(btot) = ra(i,3,2)
              xacc(btot) = ra(i,1,3)
              yacc(btot) = ra(i,2,3)
              zacc(btot) = ra(i,3,3)
              bat(btot) = at(i)
              mbuf(btot) = mass(i)
              bidx(btot) = aidx(i)
            else
              icnt = icnt + 1
              do j = 1, 3
                do k = 1, 6
                  ra(icnt,j,k) = ra(i,j,k)
                end do
              end do
              at(icnt) = at(i)
              mass(icnt) = mass(i)
              aidx(icnt) = aidx(i)
            endif
          end do
          antot = icnt
        endif
c
c    send particles in positive y direction
c
        iny = ipy - 1
        if (iny.eq.-1) iny = idy - 1
        call i_xyz_to_proc(rnode,ipx,iny,ipz,idx,idy,idz)
        call update_buf(rnode)
        call ucull
      endif
      icnt = 0
      if (idy.gt.2) then
c
c    each boundary is adjacent to a different processor
c
        icnt = 0  
        do i = 1, antot  
          if ((ra(i,2,1).lt.ymin.and.ymin-ra(i,2,1).lt.ycell).or.
     +        (ra(i,2,1).gt.ymax.and.ra(i,2,1)-ymax.gt.ycell)) then
            btot = btot + 1 
            xcrd(btot) = ra(i,1,6)
            ycrd(btot) = ra(i,2,6)
            zcrd(btot) = ra(i,3,6)
            xfrc(btot) = ra(i,1,2)
            yfrc(btot) = ra(i,2,2)
            zfrc(btot) = ra(i,3,2)
            xacc(btot) = ra(i,1,3)
            yacc(btot) = ra(i,2,3)
            zacc(btot) = ra(i,3,3)
            bat(btot) = at(i)
            mbuf(btot) = mass(i)
            bidx(btot) = aidx(i)
          else
            icnt = icnt + 1
            do j = 1, 3
              do k = 1, 6
                ra(icnt,j,k) = ra(i,j,k)
              end do
            end do
            at(icnt) = at(i)
            mass(icnt) = mass(i)
            aidx(icnt) = aidx(i)
          endif
        end do
        antot = icnt
c
c    send particles in negative y direction
c
        iny = ipy + 1
        if (iny.eq.idy) iny = 0
        call i_xyz_to_proc(rnode,ipx,iny,ipz,idx,idy,idz)
        call update_buf(rnode)
        call ucull
      endif
c
c    repeat for z axis
c
      if (idz.ge.2) then
c
c    both boundaries are adjacent to boundaries on next processor
c
        if (idz.eq.2) then
          icnt = 0
          do i = 1, antot
            if (ra(i,3,1).ge.zmax.or.ra(i,3,1).lt.zmin) then
              btot = btot + 1
              xcrd(btot) = ra(i,1,6)
              ycrd(btot) = ra(i,2,6)
              zcrd(btot) = ra(i,3,6)
              xfrc(btot) = ra(i,1,2)
              yfrc(btot) = ra(i,2,2)
              zfrc(btot) = ra(i,3,2)
              xacc(btot) = ra(i,1,3)
              yacc(btot) = ra(i,2,3)
              zacc(btot) = ra(i,3,3)
              bat(btot) = at(i)
              mbuf(btot) = mass(i)
              bidx(btot) = aidx(i)
            else
              icnt = icnt + 1
              do j = 1, 3
                do k = 1, 6
                  ra(icnt,j,k) = ra(i,j,k)
                end do
              end do
              at(icnt) = at(i)
              mass(icnt) = mass(i)
              aidx(icnt) = aidx(i)
            endif
          end do
          antot = icnt
        else
c
c    processor bounded by two separate processors
c
          icnt = 0
          do i = 1, antot
            if ((ra(i,3,1).ge.zmax.and.ra(i,3,1)-zmax.lt.zcell).or.
     +          (ra(i,3,1).lt.zmin.and.zmin-ra(i,3,1).gt.zcell)) then
              btot = btot + 1
              xcrd(btot) = ra(i,1,6)
              ycrd(btot) = ra(i,2,6)
              zcrd(btot) = ra(i,3,6)
              xfrc(btot) = ra(i,1,2)
              yfrc(btot) = ra(i,2,2)
              zfrc(btot) = ra(i,3,2)
              xacc(btot) = ra(i,1,3)
              yacc(btot) = ra(i,2,3)
              zacc(btot) = ra(i,3,3)
              bat(btot) = at(i)
              mbuf(btot) = mass(i)
              bidx(btot) = aidx(i)
            else
              icnt = icnt + 1
              do j = 1, 3
                do k = 1, 6
                  ra(icnt,j,k) = ra(i,j,k)
                end do
              end do
              at(icnt) = at(i)
              mass(icnt) = mass(i)
              aidx(icnt) = aidx(i)
            endif
          end do
          antot = icnt
        endif
c
c    send particles in positive z direction
c
        inz = ipz - 1
        if (inz.eq.-1) inz = idz - 1
        call i_xyz_to_proc(rnode,ipx,ipy,inz,idx,idy,idz)
        call update_buf(rnode)
        call ucull
      endif
      if (idz.gt.2) then
c
c    each boundary is adjacent to a different processor
c
        icnt = 0  
        do i = 1, antot  
          if ((ra(i,3,1).lt.zmin.and.zmin-ra(i,3,1).lt.zcell).or.
     +        (ra(i,3,1).gt.zmax.and.ra(i,3,1)-zmax.gt.zcell)) then
            btot = btot + 1 
            xcrd(btot) = ra(i,1,6)
            ycrd(btot) = ra(i,2,6)
            zcrd(btot) = ra(i,3,6)
            xfrc(btot) = ra(i,1,2)
            yfrc(btot) = ra(i,2,2)
            zfrc(btot) = ra(i,3,2)
            xacc(btot) = ra(i,1,3)
            yacc(btot) = ra(i,2,3)
            zacc(btot) = ra(i,3,3)
            bat(btot) = at(i)
            mbuf(btot) = mass(i)
            bidx(btot) = aidx(i)
          else
            icnt = icnt + 1
            do j = 1, 3
              do k = 1, 6
                ra(icnt,j,k) = ra(i,j,k)
              end do
            end do
            at(icnt) = at(i)
            mass(icnt) = mass(i)
            aidx(icnt) = aidx(i)
          endif
        end do
        antot = icnt
c
c    send particles in negative z direction
c
        inz = ipz + 1
        if (inz.eq.idz) inz = 0
        call i_xyz_to_proc(rnode,ipx,ipy,inz,idx,idy,idz)
        call update_buf(rnode)
        call ucull
      endif
c
c      call heapsort(0)
      do i = 1, antot
        call add_hash_a(aidx(i),i)
      end do
c
      tmstat(8) = tmstat(8) + wraptime() - tbeg
c
      return
      end
c
      subroutine update_buf(rnode)
#include "common.fh"
c
      integer ibuf(2,MAXAT)
      double precision cbuf(3,MAXAT),fbuf(3,MAXAT),abuf(3,MAXAT)
      double precision vbuf(MAXAT)
      double precision tbeg,wraptime
      integer i, stot, rtot
      integer me,rnode,one,is,ld2,ld3

c
c   Reassign all particles to the appropriate node.
c
      me = ga_nodeid()
      one = 1
      ld2 = 2
      ld3 = 3
c
      tbeg = wraptime()
c
      stot = btot
      do i = 1, stot
        cbuf(1,i) = xcrd(i)
        cbuf(2,i) = ycrd(i)
        cbuf(3,i) = zcrd(i)
        fbuf(1,i) = xfrc(i)
        fbuf(2,i) = yfrc(i)
        fbuf(3,i) = zfrc(i)
        abuf(1,i) = xacc(i)
        abuf(2,i) = yacc(i)
        abuf(3,i) = zacc(i)
        vbuf(i) = mbuf(i)
        ibuf(1,i) = bidx(i)
        ibuf(2,i) = bat(i)
      end do
c
      gsize_hi = gsize_lo(me)
      call nga_put(g_size,gsize_lo(me),gsize_hi,stot,one)
      gcoords_hi(1) = 3
      gcoords_hi(2) = gcoords_lo(2,me) + stot - 1
      if (stot.gt.0) call nga_put(g_coords,gcoords_lo(1,me),
     +                            gcoords_hi,cbuf,ld3)
      gfrc_hi(1) = 3
      gfrc_hi(2) = gfrc_lo(2,me) + stot - 1
      if (stot.gt.0) call nga_put(g_frc,gfrc_lo(1,me),gfrc_hi,fbuf,ld3)
      gacc_hi(1) = 3
      gacc_hi(2) = gacc_lo(2,me) + stot - 1
      if (stot.gt.0) call nga_put(g_acc,gacc_lo(1,me),gacc_hi,abuf,ld3)
      grvec_hi = grvec_lo(me) + stot - 1
      if (stot.gt.0) call nga_put(g_rvec,grvec_lo(me),grvec_hi,vbuf,one)
      giat_hi(1) = 2
      giat_hi(2) = giat_lo(2,me) + stot - 1
      if (stot.gt.0) call nga_put(g_iat,giat_lo(1,me),giat_hi,ibuf,ld2)
c
      call ga_sync()
c
      gsize_hi = gsize_lo(rnode)
      call nga_get(g_size,gsize_lo(rnode),gsize_hi,rtot,one)
      gcoords_hi(2) = gcoords_lo(2,rnode) + rtot - 1
      if (rtot.gt.0) call nga_get(g_coords,gcoords_lo(1,rnode),
     +                            gcoords_hi,cbuf,ld3)
      gfrc_hi(2) = gfrc_lo(2,rnode) + rtot - 1
      if (rtot.gt.0) call nga_get(g_frc,gfrc_lo(1,rnode),
     +                            gfrc_hi,fbuf,ld3)
      gacc_hi(2) = gacc_lo(2,rnode) + rtot - 1
      if (rtot.gt.0) call nga_get(g_acc,gacc_lo(1,rnode),
     +                            gacc_hi,abuf,ld3)
      grvec_hi = grvec_lo(rnode) + rtot - 1
      if (rtot.gt.0) call nga_get(g_rvec,grvec_lo(rnode),
     +                            grvec_hi,vbuf,one)
      giat_hi(2) = giat_lo(2,rnode) + rtot - 1
      if (rtot.gt.0) call nga_get(g_iat,giat_lo(1,rnode),
     +                            giat_hi,ibuf,ld2)
c
      do i = 1, rtot
        xcrd(i) = cbuf(1,i)
        ycrd(i) = cbuf(2,i)
        zcrd(i) = cbuf(3,i)
        xfrc(i) = fbuf(1,i)
        yfrc(i) = fbuf(2,i)
        zfrc(i) = fbuf(3,i)
        xacc(i) = abuf(1,i)
        yacc(i) = abuf(2,i)
        zacc(i) = abuf(3,i)
        mbuf(i) = vbuf(i)
        bidx(i) = ibuf(1,i)
        bat(i) = ibuf(2,i)
      end do
c
      btot = rtot
      call ga_sync()
c
      tmstat(13) = tmstat(13) + wraptime() - tbeg
c
      return
      end
c
      subroutine ucull
#include "common.fh"
c
      double precision xmax,ymax,zmax,xmin,ymin,zmin
      double precision xt,yt,zt
      integer pnum,inode,idx,idy,idz,ipx,ipy,ipz
      integer i,icnt
c
      pnum = ga_nnodes()
      call factor(pnum,idx,idy,idz)
      inode = ga_nodeid()
      call i_proc_to_xyz(inode,ipx,ipy,ipz,idx,idy,idz)
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
c   select out particles from buffers that belong on the home array
c
      icnt = 0
      do i = 1, btot
        xt = xcrd(i) - xbox*anint(xcrd(i)/xbox)
        yt = ycrd(i) - ybox*anint(ycrd(i)/ybox)
        zt = zcrd(i) - zbox*anint(zcrd(i)/zbox)
        antot = antot + 1
        ra(antot,1,6) = xcrd(i)
        ra(antot,2,6) = ycrd(i)
        ra(antot,3,6) = zcrd(i)
        ra(antot,1,2) = xfrc(i)
        ra(antot,2,2) = yfrc(i)
        ra(antot,3,2) = zfrc(i)
        ra(antot,1,3) = xacc(i)
        ra(antot,2,3) = yacc(i)
        ra(antot,3,3) = zacc(i)
        ra(antot,1,1) = xt
        ra(antot,2,1) = yt
        ra(antot,3,1) = zt
        at(antot) = bat(i)
        mass(antot) = mbuf(i)
        aidx(antot) = bidx(i)
      end do
      btot = icnt
c
      return
      end
