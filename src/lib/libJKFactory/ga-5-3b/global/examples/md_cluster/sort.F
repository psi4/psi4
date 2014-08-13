#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine sort
#include "common.fh"
c
      integer snode,rnode,pnum,idx,idy,idz,ipx,ipy,ipz,ilast,ifirst
      integer i,icnt,inode,mynode
c
c   This subroutine sorts the particle data onto each of the processors
c   so that each processor has a list of consecutively numbered particles.
c
      pnum = ga_nnodes()
      call factor(pnum,idx,idy,idz)
      inode = ga_nodeid()
      call i_proc_to_xyz(inode,ipx,ipy,ipz,idx,idy,idz)
c
c   Find indices of first and last particle that should be on this
c   processor
c
      ilast = nint(dble((ga_nodeid()+1)*atot)/dble(pnum))
      ifirst = nint(dble(ga_nodeid()*atot)/dble(pnum))
      ifirst = ifirst + 1
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
c   send buffers to all nodes in systolic loop
c   and select out particles that fall within the appropriate
c   range of indices
c
      mynode = ga_nodeid()
      call icull(ifirst,ilast)
      do inode = 1, pnum-1
        snode = ga_nodeid() + 1
        if (snode.eq.pnum) snode = 0
        rnode = ga_nodeid() - 1
        if (rnode.eq.-1) rnode = pnum - 1
        call exchange_buf(rnode)
        call icull(ifirst,ilast)
      end do
c
c   rearrange data stack so that they are ordered with respect to
c   atom index
c
      call heapsort(0)
      call fixper
      icnt = 0
c
      return
      end
c
      subroutine icull(ifirst,ilast)
#include "common.fh"
c
      integer ifirst,ilast,i,icnt
c
c   This subroutine gathers all the particles in the buffer arrays that
c   have indices lying in the domain [ifirst,ilast]
c   and puts them in the regular particle arrays
c
      icnt = 0
      do i = 1, btot
        if (bidx(i).ge.ifirst.and.bidx(i).le.ilast) then
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
c
      return
      end
