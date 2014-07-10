#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine cell_list
#include "common.fh"
c
      integer MAXCL,MSKIN
      parameter (MAXCL=50,MSKIN=5)
      double precision xu,yu,zu,xl,yl,zl
      double precision rct2sd,rbcut,rbcut2
      double precision r2(MAXAT),dx(MAXAT),dy(MAXAT),dz(MAXAT)
      double precision br2(MAXAT),bdx(MAXAT),bdy(MAXAT),bdz(MAXAT)
      double precision xcl,ycl,zcl,xmax,ymax,zmax,xmin,ymin,zmin
      double precision rix,riy,riz
      double precision tbeg,wraptime
      integer link(MAXAT)
      integer ltop(-MSKIN:MAXCL+MSKIN,-MSKIN:MAXCL+MSKIN,
     +             -MSKIN:MAXCL+MSKIN)
      integer btop(-MSKIN:MAXCL+MSKIN,-MSKIN:MAXCL+MSKIN,
     +             -MSKIN:MAXCL+MSKIN),blink(MAXAT)
      integer acml(MAXAT),bacml(MAXAT)
      integer kxmax,kymax,kzmax,kc,kx,ky,kz,jx,jy,jz,ltot,bltot
      integer i,j,k,ipx,ipy,ipz,nlx,nly,nlz,ix,iy,iz
      integer idx,idy,idz,pnum,inode,ll,kk,l,lb,jj
      double precision rchk
c
c    This subroutine locates the nearest neighbors of each atom
c    using a link cell algorithm and constructs a neighbor list
c
      tbeg = wraptime()
c
      pnum = ga_nnodes()
      call factor(pnum,idx,idy,idz)
      inode = ga_nodeid()
      call i_proc_to_xyz(inode,ipx,ipy,ipz,idx,idy,idz)
c
c   determine the size of each spatial decomposition zone
c   and then use this to divide each zone into cells
c
      rcell = 1.1d00
      rbcut = 1.5d00
      rchk = 0.0d00
c
      rbcut2 = rbcut**2
      xcl = xbox / dble(idx)
      ycl = ybox / dble(idy)
      zcl = zbox / dble(idz)
      if (int(xcl/rcell).ge.MAXCL) then
        rcell = xcl/(dble(MAXCL)) + 0.1
      endif
      if (int(ycl/rcell).ge.MAXCL) then
        rcell = ycl/(dble(MAXCL)) + 0.1
      endif
      if (int(zcl/rcell).ge.MAXCL) then
        rcell = zcl/(dble(MAXCL)) + 0.1
      endif
      nlx = int(xcl/rcell)
      nly = int(ycl/rcell)
      nlz = int(zcl/rcell)
      xcl = xcl / dble(nlx)
      ycl = ycl / dble(nly)
      zcl = zcl / dble(nlz)
c
c  If system is small, just use brute force neighbor list construction
c  and skip the link cell algorithm
c
      ix = int(rcmax/xcl)
      if (rcmax-xcl*dble(ix).gt.0.0d00) ix = ix+1
      iy = int(rcmax/ycl)
      if (rcmax-ycl*dble(iy).gt.0.0d00) iy = iy+1
      iz = int(rcmax/zcl)
      if (rcmax-zcl*dble(iz).gt.0.0d00) iz = iz+1
      if ((idx.eq.2.and.xcl*dble(ix).ge.0.5d00*xbox/dble(idx)).or.
     +    (idy.eq.2.and.ycl*dble(iy).ge.0.5d00*ybox/dble(idy)).or.
     +    (idz.eq.2.and.zcl*dble(iz).ge.0.5d00*zbox/dble(idz))) then
        call small_list
        return
      endif

c
c   find link cells
c
      call getcell
c
c   find location of cell edges
c
      xmax = xbox * dble(ipx+1)/dble(idx)
      ymax = ybox * dble(ipy+1)/dble(idy)
      zmax = zbox * dble(ipz+1)/dble(idz)
      xmax = xmax - xbox2
      ymax = ymax - ybox2
      zmax = zmax - zbox2
      xmin = xbox * dble(ipx)/dble(idx)
      ymin = ybox * dble(ipy)/dble(idy)
      zmin = zbox * dble(ipz)/dble(idz)
      xmin = xmin - xbox2
      ymin = ymin - ybox2
      zmin = zmin - zbox2
c
c   calculate link cell indices. Start by setting all
c   elements of ltop() equal to 0.
c
      if (nlx.ge.MAXCL) then
        write(6,105) ga_nodeid(),istep
        write(6,*) 'NLX = ',nlx,' ON NODE',ga_nodeid()
        call ga_error("Out of bound cell index",0)
      elseif (nly.ge.MAXCL) then
        write(6,105) ga_nodeid(),istep
        write(6,*) 'NLY = ',nly,' ON NODE',ga_nodeid()
        call ga_error("Out of bound cell index",0)
      elseif (nlz.ge.MAXCL) then
        write(6,105) ga_nodeid(),istep
        write(6,*) 'NLZ = ',nlz,' ON NODE',ga_nodeid()
        call ga_error("Out of bound cell index",0)
      endif
      do i = -MSKIN, nlx+MSKIN
        do j = -MSKIN, nly+MSKIN
          do k = -MSKIN, nlz+MSKIN
            ltop(i,j,k) = 0
            btop(i,j,k) = 0
          end do
        end do
      end do
  105 format(i5,' Out of bound failure at step ',i8)
c
c  find leading index for each cell and store it in
c  ltop(). Store link list in list(). Store atoms
c  in both local and buffer arrays.
c
      do i = 1, antot
        ix = int((ra(i,1,1)-xmin)/xcl)+1
        iy = int((ra(i,2,1)-ymin)/ycl)+1
        iz = int((ra(i,3,1)-zmin)/zcl)+1
        j = ltop(ix,iy,iz)
        ltop(ix,iy,iz) = i
        link(i) = j
#if 0
        if (ix.lt.1.or.ix.gt.nlx.or.iy.lt.1.or.iy.gt.nly.or.
     +     iz.lt.1.or.iz.gt.nlz) then
          write(6,*) 'i = ',i
          write(6,*) 'antot = ',antot
          write(6,*) 'atot = ',atot
          write(6,*) 'ra(i,1,1) = ',ra(i,1,1)
          write(6,*) 'ra(i,2,1) = ',ra(i,2,1)
          write(6,*) 'ra(i,3,1) = ',ra(i,3,1)
          write(6,*) 'xmin = ',xmin
          write(6,*) 'ymin = ',ymin
          write(6,*) 'zmin = ',zmin
          write(6,*) 'xcl = ',xcl
          write(6,*) 'ycl = ',ycl
          write(6,*) 'zcl = ',zcl
          write(6,*) 'ix ',ix,nlx
          write(6,*) 'iy ',iy,nly
          write(6,*) 'iz ',iz,nlz
          write(6,105) ga_nodeid(),istep
          write(6,*) 'OUT OF BOUNDS AT 1 ON NODE ',ga_nodeid()
          call ga_error("Out of bounds error",0)
        endif
#endif
      end do
      do i = 1, btot
        if ((xcrd(i).lt.xmin.and.xmin-xcrd(i).lt.rcmax).or.
     +      (xcrd(i).gt.xmax.and.xcrd(i)-xmax.gt.rcmax)) then
          xl = xmin - xcrd(i)
          xl = xl - xbox * anint(xl/xbox)
          ix = -int(xl/xcl)
        elseif ((xcrd(i).ge.xmax.and.xcrd(i)-xmax.lt.rcmax).or.
     +      (xcrd(i).lt.xmin.and.xmin-xcrd(i).gt.rcmax)) then
          xu = xcrd(i) - xmax
          xu = xu - xbox * anint(xu/xbox)
          ix = nlx + int(xu/xcl) + 1
        else
          ix = int((xcrd(i)-xmin)/xcl)+1
        endif
        if ((ycrd(i).lt.ymin.and.ymin-ycrd(i).lt.rcmax).or.
     +      (ycrd(i).gt.ymax.and.ycrd(i)-ymax.gt.rcmax)) then
          yl = ymin - ycrd(i)
          yl = yl - ybox * anint(yl/ybox)
          iy = -int(yl/ycl)
        elseif ((ycrd(i).ge.ymax.and.ycrd(i)-ymax.lt.rcmax).or.
     +      (ycrd(i).lt.ymin.and.ymin-ycrd(i).gt.rcmax)) then
          yu = ycrd(i) - ymax
          yu = yu - ybox * anint(yu/ybox)
          iy = nly + int(yu/ycl) + 1
        else
          iy = int((ycrd(i)-ymin)/ycl)+1
        endif
        if ((zcrd(i).lt.zmin.and.zmin-zcrd(i).lt.rcmax).or.
     +      (zcrd(i).gt.zmax.and.zcrd(i)-zmax.gt.rcmax)) then
          zl = zmin - zcrd(i)
          zl = zl - zbox * anint(zl/zbox)
          iz = -int(zl/zcl)
        elseif ((zcrd(i).ge.zmax.and.zcrd(i)-zmax.lt.rcmax).or.
     +      (zcrd(i).lt.zmin.and.zmin-zcrd(i).gt.rcmax)) then
          zu = zcrd(i) - zmax
          zu = zu - zbox * anint(zu/zbox)
          iz = nlz + int(zu/zcl) + 1
        else
          iz = int((zcrd(i)-zmin)/zcl)+1
        endif
        j = btop(ix,iy,iz)
        btop(ix,iy,iz) = i
        blink(i) = j
#if 0
        if (ix.ge.1.and.ix.le.nlx.and.iy.ge.1.and.iy.le.nly.and.
     +     iz.ge.1.and.iz.le.nlz) then
          write(6,105) ga_nodeid(),istep
          write(6,*) 'OUT OF BOUNDS AT 2 ON NODE ',ga_nodeid()
          call ga_error("Out of bounds error",0)
        endif
#endif
      end do
c
c  Link cell lists are complete. Now find neighbors and construct
c  neighbor lists.
c
c  Primary loop over all local atoms
c
      ll = 1
      kk = 1
      do i = 1, antot
c
c    Locate all atoms in both the local and buffer lists that
c    interact with atom i. Start by collecting all atoms in the
c    same cell
c
        rix = ra(i,1,1)
        riy = ra(i,2,1)
        riz = ra(i,3,1)
c
        l = 0
        lb = 0
        j = link(i)
        if (j.gt.0) then
   99     l = l + 1
          acml(l) = j
          dx(l) = ra(j,1,1)
          dy(l) = ra(j,2,1)
          dz(l) = ra(j,3,1)
          j = link(j)
          if (j.gt.0) go to 99
        endif
c
c    Collect atoms in half of the adjacent cells.
c    Find indices of cell containing atom i
c
        ix = int((ra(i,1,1)-xmin)/xcl)+1
        iy = int((ra(i,2,1)-ymin)/ycl)+1
        iz = int((ra(i,3,1)-zmin)/zcl)+1
c
c   loop over neighboring cells
c
        do kc = 1, kcmax
          jx = ix + nix(kc)
          if (jx.gt.nlx.and.idx.eq.1) jx = jx - nlx
          if (jx.lt.1.and.idx.eq.1) jx = nlx + jx
          jy = iy + niy(kc)
          if (jy.gt.nly.and.idy.eq.1) jy = jy - nly
          if (jy.lt.1.and.idy.eq.1) jy = nly + jy
          jz = iz + niz(kc)
          if (jz.gt.nlz.and.idz.eq.1) jz = jz - nlz
          if (jz.lt.1.and.idz.eq.1) jz = nlz + jz
          if ((jx.ge.1.and.jx.le.nlx).and.
     +        (jy.ge.1.and.jy.le.nly).and.
     +        (jz.ge.1.and.jz.le.nlz)) then
c
c  neighboring cell is local
c
            j = ltop(jx,jy,jz)
c
c    gather atoms in neighboring cell
c
            if (j.gt.0) then
  199         l = l + 1
              acml(l) = j
              dx(l) = ra(j,1,1)
              dy(l) = ra(j,2,1)
              dz(l) = ra(j,3,1)
              j = link(j)
              if (j.gt.0) go to 199
            endif
          else
c
c  neighboring cell is in buffer
c
            j = btop(jx,jy,jz)
c
c    gather atoms in neighboring cell
c
            if (j.gt.0) then
  299         lb = lb + 1
              bacml(lb) = j
              bdx(lb) = xcrd(j)
              bdy(lb) = ycrd(j)
              bdz(lb) = zcrd(j)
              j = blink(j)
              if (j.gt.0) go to 299
            endif
          endif
        end do
        ltot = l
        bltot = lb
c
        do l = 1, ltot
          dx(l) = rix - dx(l)
          dx(l) = dx(l) - xbox * anint(dx(l)/xbox)
          dy(l) = riy - dy(l)
          dy(l) = dy(l) - ybox * anint(dy(l)/ybox)
          dz(l) = riz - dz(l)
          dz(l) = dz(l) - zbox * anint(dz(l)/zbox)
          r2(l) = dx(l)**2 + dy(l)**2 + dz(l)**2
        end do
        do lb = 1, bltot
          bdx(lb) = rix - bdx(lb)
          bdx(lb) = bdx(lb) - xbox * anint(bdx(lb)/xbox)
          bdy(lb) = riy - bdy(lb)
          bdy(lb) = bdy(lb) - ybox * anint(bdy(lb)/ybox)
          bdz(lb) = riz - bdz(lb)
          bdz(lb) = bdz(lb) - zbox * anint(bdz(lb)/zbox)
          br2(lb) = bdx(lb)**2 + bdy(lb)**2 + bdz(lb)**2
        end do
        if (ltot.gt.antot.or.bltot.gt.btot) then
          write(6,*) 'NODE ',ga_nodeid(),ltot,bltot,antot,btot
          stop
        endif
c
c   accumulate all atoms within an interaction distance of
c   atom i
c
        nafirst(i) = ll
        bafirst(i) = kk
        do 300 j = 1, ltot
          jj = acml(j)
          rct2sd = (acut(at(i),at(jj))+1.0d00)**2
          if (rct2sd.gt.rchk) rchk = rct2sd
c
          if (r2(j).gt.rct2sd) go to 300
c          if (atml(i).eq.atml(jj)) then
c            if (r2(j).le.rbcut2) then
c              blist(kk) = jj
c              kk = kk + 1
c            endif
c            go to 300
c          endif
c
          nblist(ll) = jj
          ll = ll + 1
  300   continue
        nalast(i) = ll - 1
        nbfirst(i) = ll
        balast(i) = kk - 1
        bbfirst(i) = kk
        do 400 j = 1, bltot
          jj = bacml(j)
          rct2sd = (acut(at(i),bat(jj))+1.0d00)**2
          if (rct2sd.gt.rchk) rchk = rct2sd
          if (br2(j).gt.rct2sd) go to 400
c          if (atml(i).eq.batml(jj)) then
c            if (br2(j).le.rbcut2) then
c              blist(kk) = jj
c              kk = kk + 1
c            endif
c            go to 400
c          endif
c
          nblist(ll) = jj
          ll = ll + 1
  400   continue
        nblast(i) = ll - 1
        bblast(i) = kk - 1
      end do
      tmstat(23) = tmstat(23) + wraptime() - tbeg
c
      return
      end
c
      subroutine getcell
#include "common.fh"
c
c  This subroutine finds the total number of nearest of cells
c  required by the link-cell routine
c
      double precision xcell,ycell,zcell,cx,cy,cz,rc,rcmax2
      integer kxmin,kymin,kzmin,kxmax,kymax,kzmax,kc,kx,ky,kz
      integer pnum,idx,idy,idz,nlx,nly,nlz
c
      pnum = ga_nnodes()
      call factor(pnum,idx,idy,idz)
c
      xcell = xbox / dble(idx)
      ycell = ybox / dble(idy)
      zcell = zbox / dble(idz)
      nlx = int(xcell/rcell)
      nly = int(ycell/rcell)
      nlz = int(zcell/rcell)
      xcell = xcell / dble(nlx)
      ycell = ycell / dble(nly)
      zcell = zcell / dble(nlz)
c
      kxmax = int(rcmax/xcell)+1
      kymax = int(rcmax/ycell)+1
      kzmax = int(rcmax/zcell)+1
c
      rcmax2 = rcmax**2
c
      kc = 0
      do kx = 0, kxmax
        kymin = -kymax
        if (kx.eq.0) kymin = 0
        do ky = kymin, kymax
          kzmin = -kzmax
          if (kx.eq.0.and.ky.eq.0) kzmin = 1
          do kz = kzmin, kzmax
            if (kx.eq.0) then
              cx = 0.0d00
            else
              cx = xcell*dble(kx-1)
            endif
            if (ky.gt.0) then
              cy = ycell*dble(ky-1)
            elseif (ky.lt.0) then
              cy = ycell*dble(-ky-1)
            else
              cy = 0.0d00
            endif
            if (kz.gt.0) then
              cz = zcell*dble(kz-1)
            elseif (kz.lt.0) then
              cz = zcell*dble(-kz-1)
            else
              cz = 0.0d00
            endif
            rc = cx**2+cy**2+cz**2
            if (rc.lt.rcmax2) then
              kc = kc + 1
              nix(kc) = kx
              niy(kc) = ky
              niz(kc) = kz
            endif
          end do
        end do
      end do
      kcmax = kc
      return
      end
c
      subroutine small_list
#include "common.fh"
c
      integer MAXCL,MSKIN
      parameter (MAXCL=50,MSKIN=5)
      double precision rcmax2
      double precision r2,dx,dy,dz
      double precision rix,riy,riz,rjx,rjy,rjz
      double precision tbeg,wraptime
      integer acml(MAXAT),bacml(MAXAT)
      integer i,j,ll
c
c    This subroutine locates the nearest neighbors of each atom
c    using a direct search over all pairs. It is used to construct
c    the neighbor list for small systems where the link cell
c    algorithm may fail.
c
      tbeg = wraptime()
c
c  Find neighbors and construct neighbor lists.
c  Loop over all local atoms
c
      rcmax2 = rcmax**2
      ll = 1
      do i = 1, antot
c
c    Locate all atoms in both the local and buffer lists that
c    interact with atom i.
c
        rix = ra(i,1,1)
        riy = ra(i,2,1)
        riz = ra(i,3,1)
        nafirst(i) = ll
        do j = i+1, antot
          rjx = ra(j,1,1)
          rjy = ra(j,2,1)
          rjz = ra(j,3,1)
          dx = rix - rjx
          dx = dx - xbox * anint(dx/xbox)
          dy = riy - rjy
          dy = dy - ybox * anint(dy/ybox)
          dz = riz - rjz
          dz = dz - zbox * anint(dz/zbox)
          r2 = dx**2 + dy**2 + dz**2
          if (r2.le.rcmax2) then
            nblist(ll) = j
            ll = ll + 1
          endif
        end do
        nalast(i) = ll - 1
        nbfirst(i) = ll
        do j = 1, btot
          if (bidx(j).gt.aidx(i)) then
            rjx = xcrd(j)
            rjy = ycrd(j)
            rjz = zcrd(j)
            dx = rix - rjx
            dx = dx - xbox * anint(dx/xbox)
            dy = riy - rjy
            dy = dy - ybox * anint(dy/ybox)
            dz = riz - rjz
            dz = dz - zbox * anint(dz/zbox)
            r2 = dx**2 + dy**2 + dz**2
            if (r2.le.rcmax2) then
              nblist(ll) = j
              ll = ll + 1
            endif
          endif
        end do
        nblast(i) = ll - 1
      end do
c
      tmstat(23) = tmstat(23) + wraptime() - tbeg
c
      return
      end
