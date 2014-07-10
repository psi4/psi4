#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine pairs
#include "common.fh"
c
      integer MAXCL,MSKIN
      parameter (MAXCL=50,MSKIN=5)
      double precision savd(MAXAT),savr(MAXAT),ft(MAXAT,3)
      double precision rsqi,r6i,dedr,dedrd,dedrr,phid,phir
      double precision xu,yu,zu,xl,yl,zl,rcut2
      double precision tprssr(MAXAT,3),a12(MAXAT),a6(MAXAT)
      double precision rct2sd,rsqcrt,rsq,b12(MAXAT),b6(MAXAT)
      double precision r2(MAXAT),dx(MAXAT),dy(MAXAT),dz(MAXAT)
      double precision br2(MAXAT),bdx(MAXAT),bdy(MAXAT),bdz(MAXAT)
      double precision xcl,ycl,zcl,xmax,ymax,zmax,xmin,ymin,zmin
      double precision rcp(MAXAT),dcp(MAXAT),brcp(MAXAT),bdcp(MAXAT)
      double precision one,six,twelve
      double precision rix,riy,riz
      double precision tbeg,wraptime
      integer acml(MAXAT),bacml(MAXAT)
      integer i,j,l,fat,lat,jp,jtmp,jat,ltot,bltot
c
      if (istep.eq.0.or.(mod(istep,ilist).eq.0.and.
     +    t_rmndr.eq.0.0d00)) call cell_list
      tbeg = wraptime()
c
c   this subroutine calculates the forces on individual particles
c   due to the Van der Waals interactions
c
c   calculate some parameters
c
      one = 1.0d00
      six = 6.0d00
      twelve = 12.0d00
c
c  Primary loop over all local atoms
c
      do i = 1, antot
c
c   evaluate forces for both local and buffer atoms
c
        fat = nafirst(i)
        lat = nalast(i)
c
c   initialize temporary storage arrays
c
        l = 0
        do j = fat, lat
          l = l + 1
          savr(l) = 0.0d00
          savd(l) = 0.0d00
          ft(l,1) = 0.0d00
          ft(l,2) = 0.0d00
          ft(l,3) = 0.0d00
          tprssr(l,1) = 0.0d00
          tprssr(l,2) = 0.0d00
          tprssr(l,3) = 0.0d00
        end do
c
        jp = 0
        do jat = fat, lat
          j = nblist(jat)
          jp = jp + 1
c
          acml(jp) = j
c
c    calculate interatomic distances
c
          dx(jp) = ra(i,1,6) - ra(j,1,6)
          dx(jp) = dx(jp) - xbox * anint(dx(jp) / xbox)
          dy(jp) = ra(i,2,6) - ra(j,2,6)
          dy(jp) = dy(jp) - ybox * anint(dy(jp) / ybox)
          dz(jp) = ra(i,3,6) - ra(j,3,6)
          dz(jp) = dz(jp) - zbox * anint(dz(jp) / zbox)
c
          r2(jp) = dx(jp)**2 + dy(jp)**2 + dz(jp)**2
        end do
        ltot = jp
c
c  check displacements against cutoff
c
        jp = 0
        do 1200 j = 1, ltot
          rcut2 = acut2(at(i),at(acml(j)))
          if (r2(j).ge.rcut2) go to 1200
          jp = jp + 1
          jtmp = acml(j)
          acml(jp) = jtmp
c
          r2(jp) = r2(j)
c
          a12(jp) = e12(at(jtmp),at(i))
          a6(jp) = e6(at(jtmp),at(i))
c
          dcp(jp) = -dcorr(at(jtmp),at(i))
          rcp(jp) = -rcorr(at(jtmp),at(i))
c
          dx(jp) = dx(j)
          dy(jp) = dy(j)
          dz(jp) = dz(j)
 1200   continue
        ltot = jp
c
c   evaluate interaction between all local pairs
c
c
        do j = 1, ltot
          rsq = r2(j)
          rsqi = one / rsq
          r6i = rsqi * rsqi * rsqi
c
          phir =  a12(j) * r6i * r6i
          phid = -a6(j) * r6i
          savr(j) = savr(j) + phir + rcp(j)
          savd(j) = savd(j) + phid + dcp(j)
          dedrr = twelve * phir
          dedrd = six * phid
          dedr = (dedrr + dedrd) * rsqi
          ft(j,1) = dedr * dx(j)
          ft(j,2) = dedr * dy(j)
          ft(j,3) = dedr * dz(j)
          tprssr(j,1) = dedr*dx(j)**2
          tprssr(j,2) = dedr*dy(j)**2
          tprssr(j,3) = dedr*dz(j)**2
        end do
c
c   sum forces on atom i and contributions to the energy and
c   pressure
c
        do j = 1, ltot
          ra(i,1,4) = ra(i,1,4) + ft(j,1)
          ra(i,2,4) = ra(i,2,4) + ft(j,2)
          ra(i,3,4) = ra(i,3,4) + ft(j,3)
          esvr = esvr + savr(j)
          esvd = esvd + savd(j)
          presf = presf + tprssr(j,1) + tprssr(j,2) + tprssr(j,3)
          presfx = presfx + tprssr(j,1)
          presfy = presfy + tprssr(j,2)
          presfz = presfz + tprssr(j,3)
        end do
c
c   scatter forces to atoms j
c
        do jp = 1, ltot
          j = acml(jp)
          ra(j,1,4) = ra(j,1,4) - ft(jp,1)
          ra(j,2,4) = ra(j,2,4) - ft(jp,2)
          ra(j,3,4) = ra(j,3,4) - ft(jp,3)
        end do
c
c   Calculate interactions between locally held atoms and remote
c   atoms. Start by re-initializing temporary storage arrays.
c
        fat = nbfirst(i)
        lat = nblast(i)
c
        l = 0
        do j = fat, lat
          l = l + 1
          savr(l) = 0.0d00
          savd(l) = 0.0d00
          ft(l,1) = 0.0d00
          ft(l,2) = 0.0d00
          ft(l,3) = 0.0d00
          tprssr(l,1) = 0.0d00
          tprssr(l,2) = 0.0d00
          tprssr(l,3) = 0.0d00
        end do
c
        jp = 0
        do jat = fat, lat
          j = nblist(jat)
          jp = jp + 1
c
          bacml(jp) = j
c
c    calculate interatomic distances
c
          bdx(jp) = ra(i,1,6) - xcrd(j)
          bdx(jp) = bdx(jp) - xbox * anint(bdx(jp) / xbox)
          bdy(jp) = ra(i,2,6) - ycrd(j)
          bdy(jp) = bdy(jp) - ybox * anint(bdy(jp) / ybox)
          bdz(jp) = ra(i,3,6) - zcrd(j)
          bdz(jp) = bdz(jp) - zbox * anint(bdz(jp) / zbox)
c
          br2(jp) = bdx(jp)**2 + bdy(jp)**2 + bdz(jp)**2
        end do
        bltot = jp
c
c  check displacements against cutoff
c
        jp = 0
        do 1250 j = 1, bltot
          rcut2 = acut2(at(i),bat(bacml(j)))
          if (br2(j).ge.rcut2) go to 1250
          jp = jp + 1
          jtmp = bacml(j)
          bacml(jp) = jtmp
c
          br2(jp) = br2(j)
c
          b12(jp) = e12(bat(jtmp),at(i))
          b6(jp) = e6(bat(jtmp),at(i))
c
          bdcp(jp) = -dcorr(bat(jtmp),at(i))
          brcp(jp) = -rcorr(bat(jtmp),at(i))
c
          bdx(jp) = bdx(j)
          bdy(jp) = bdy(j)
          bdz(jp) = bdz(j)
 1250   continue
        bltot = jp
c
c   evaluate interaction between all local pairs
c
        do j = 1, bltot
          rsq = br2(j)
          rsqi = one / rsq
          r6i = rsqi * rsqi * rsqi
c
          phir =  b12(j) * r6i * r6i
          phid = -b6(j) * r6i
          savr(j) = savr(j) + phir + brcp(j)
          savd(j) = savd(j) + phid + bdcp(j)
          dedrr = twelve * phir
          dedrd = six * phid
          dedr = (dedrr + dedrd) * rsqi
          ft(j,1) = dedr * bdx(j)
          ft(j,2) = dedr * bdy(j)
          ft(j,3) = dedr * bdz(j)
          tprssr(j,1) = dedr*bdx(j)**2
          tprssr(j,2) = dedr*bdy(j)**2
          tprssr(j,3) = dedr*bdz(j)**2
        end do
c
c   sum forces on atom i and contributions to the energy and
c   pressure
c
        do j = 1, bltot
          ra(i,1,4) = ra(i,1,4) + ft(j,1)
          ra(i,2,4) = ra(i,2,4) + ft(j,2)
          ra(i,3,4) = ra(i,3,4) + ft(j,3)
          esvr = esvr + savr(j)
          esvd = esvd + savd(j)
          presf = presf + tprssr(j,1) + tprssr(j,2) + tprssr(j,3)
          presfx = presfx + tprssr(j,1)
          presfy = presfy + tprssr(j,2)
          presfz = presfz + tprssr(j,3)
        end do
c
c   scatter forces to atoms j
c
        do jp = 1, bltot
          j = bacml(jp)
          xfrc(j) = xfrc(j) - ft(jp,1)
          yfrc(j) = yfrc(j) - ft(jp,2)
          zfrc(j) = zfrc(j) - ft(jp,3)
        end do
      end do
c
      tmstat(3) = tmstat(3) + wraptime() - tbeg
c
      return
      end
