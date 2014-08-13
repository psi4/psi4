#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      program main

c
#include "common"
c
      integer i, j, gdim, type
      integer gdims(NDIM+1), gwidth(NDIM+1)
      integer pdims(NDIM+1), mcnt, mapc(5000)
      integer g_fg, g_fld, g_bc, ptr_fg, ptr_fld, ptr_bc
      integer ld_fg(NDIM+1), ld_fld(NDIM+1), ld_bc(NDIM)
      integer heap, stack, me, nproc
c
#include "mafdecls.fh"
#include "global.fh"
c
c  Initialize a message passing library
c
#include "mp3.fh"
      size(1) = NSIZE
      size(2) = NSIZE
c
c   Initialize global arrays
c
      call ga_initialize()
c
      nproc = ga_nnodes()
      me = ga_nodeid()
c
      if (ga_uses_ma()) then
        heap = (size(1)+2)*(size(2)+2)*34/nproc
      else
        heap = 100000
      endif
      stack = 50000
c
      if (.not.ma_init(MT_DBL, stack, heap))
     +  call ga_error("ma init failed", -1)

c
c initialize global arrays
c
      do i = 1, NDIM
        dims(i) = size(i)
        gdims(i) = size(i)
        width(i) = 1
        gwidth(i) = 1
      end do
c
c  evaluate distribution of processors
c
      gdim = NDIM
      call factor(nproc,gdim,pdims)
      mcnt = 1
      do i = 1, NDIM
        do j = 0, pdims(i)-1
          mapc(mcnt) = ((dble(j)/dble(pdims(i)))*dble(gdims(i)))+1
          mcnt = mcnt + 1
        end do
      end do
      do i = 0, pdims(1)-1
        mapc(mcnt) = ((dble(i)/dble(pdims(1)))*dble(NDIM))+1
        mcnt = mcnt + 1
      end do
      do i = 0, pdims(2)-1
        mapc(mcnt) = ((dble(i)/dble(pdims(2)))*dble(NDIM))+1
        mcnt = mcnt + 1
      end do
c
c   Create global arrays. Start by creating array for LB distribution
c   functions. The last dimension runs over the distribution function
c   indices. The first 9 elements are the actual distribution elements,
c   the next 9 indices are the equilibrium distribution elements,
c   and the last 9 elements are temporary storage space used for doing
c   the streaming updates.
c
      type = MT_DBL
      gdim = 3
      gdims(3) = 27
      gwidth(3) = 0
      pdims(3) = 1
      mapc(mcnt) = 1
      if (.not.nga_create_ghosts_irreg(type, gdim, gdims, gwidth,
     +    "lb_dist", mapc, pdims, g_fg))
     +    call ga_error("g_fg init failed",me)
c
c   Create global array to hold density, momentum, pressure,
c   and relaxation parameters. These are stored at each point
c   and indexed by the last indice as density, p_x, p_y,
c   pressure, t_rho.
c
      type = MT_DBL
      gdim = 3
      gdims(3) = 6
      gwidth(3) = 0
      pdims(3) = 1
      mapc(mcnt) = 1
      if (.not.nga_create_ghosts_irreg(type, gdim, gdims, gwidth,
     +    "fields", mapc, pdims, g_fld))
     +    call ga_error("g_fld init failed",me)
c
c   Create global array to hold boundary condition data.
c
      type = MT_INT
      gdim = 2
      if (.not.nga_create_ghosts_irreg(type, gdim, gdims, gwidth,
     +    "bc_mask", mapc, pdims, g_bc))
     +    call ga_error("g_bc init failed",me)
c
c   Find pointers to global array data
c
      call nga_access_ghosts(g_fg,dims_fg,ptr_fg,ld_fg)
      call nga_access_ghosts(g_fld,dims_fld,ptr_fld,ld_fld)
      call nga_access_ghosts(g_bc,dims_bc,ptr_bc,ld_bc)
      do i = 1, NDIM
        dims(i) = dims_fg(i)
      end do
c
c   Call routine to run main simulation
c
      call boltzmann(g_fg,  dbl_mb(ptr_fg),  ld_fg(1),  ld_fg(2),
     +               g_fld, dbl_mb(ptr_fld), ld_fld(1), ld_fld(2),
     +               g_bc,  int_mb(ptr_bc),  ld_bc(1))
c
c   Close out calculation
c
      call ga_terminate()
      call MP_FINALIZE()
      stop
      end
c
      subroutine factor(p,ndim,dims)
      implicit none
      integer i,j,p,ndim,dims(*),imin,mdim
      integer ip,ifac,pmax,prime(1000)
      integer fac(1000)
c
      i = 1
      ip = p
      do i = 1, ndim
        dims(i) = 1
      end do
c
c    factor p completely
c    first, find all prime numbers less than or equal to p
c
      pmax = 0
      do i = 2, p
        do j = 1, pmax
          if (mod(i,prime(j)).eq.0) go to 100
        end do
        pmax = pmax + 1
        prime(pmax) = i
  100   continue
      end do
c
c    find all prime factors of p
c
      ifac = 0
      do i = 1, pmax
  200   if (mod(ip,prime(i)).eq.0) then
          ifac = ifac + 1
          fac(ifac) = prime(i)
          ip = ip/prime(i)
          go to 200
        endif
      end do
c
c    determine dimensions of processor grid
c
      do i = ifac, 1, -1
c
c    find dimension with minimum value
c
        imin = dims(1)
        mdim = 1
        do j = 2, ndim
          if (dims(j).lt.imin) then
            imin = dims(j)
            mdim = j
          endif
        end do
        dims(mdim) = dims(mdim)*fac(i)
      end do
c
      return
      end
