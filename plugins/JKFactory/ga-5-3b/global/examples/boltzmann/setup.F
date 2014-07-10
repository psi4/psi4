#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine setup(g_fg,  fg,  ld_fg1,  ld_fg2,
     +                 g_fld, fld, ld_fld1, ld_fld2,
     +                 g_bc,  bc,  ld_bc1)
#include "common"
c
      integer ld_fg1, ld_fg2, ld_fld1, ld_fld2, ld_bc1
      double precision fg(ld_fg1,ld_fg2, *)
      double precision fld(ld_fld1, ld_fld2, *)
      integer bc(ld_bc1, *)
      integer g_fg, g_fld, g_bc 
c
#include "mafdecls.fh"
#include "global.fh"
c
      double precision rho0, ux0, uy0, t_rho, pi
      integer i, j, ii, jj, gtmp, me
      integer ld(NDIM)
c
c   Set simulation parameters
c
      nsteps = 5000
      viscosity = 1.0d00
      delta_t = 1.0d00
      xmax = 256.0
c
      rho0 = 2.7
      tmprtr0 = 1.0
      ux0 = 0.0
      uy0 = 0.0
      t_rho = 1.0d00
      uxbc = 0.5
      rhobc = rho0
c
      rgas = 1.0
      a_vdw = 0.0d00
      b_vdw = 0.0d00
c
      delta_x = xmax/dble(size(1)-1)
      cspd = sqrt(2.0d00)*delta_x/delta_t
      if (ga_nodeid().eq.0) write(6,*) 'Value of tau_rho is ',
     +      6.0d00*viscosity/(cspd**2*delta_t) + 0.5d00
c
c   Find low and high indices of locally held data
c
      me = ga_nodeid()
      call nga_distribution(g_fg, me, lo, hi)
      call nga_access(g_fg, lo, hi, gtmp, ld)
c
c   Initialize boundary array
c
      call ga_zero(g_bc)
      do jj = width(2) + 1, dims(2) - width(2)
        j = jj - width(2) - 1 + lo(2)
        do ii = width(1) + 1, dims(1) - width(1)
          i = ii - width(1) - 1 + lo(1)
          if (i.eq.1) then
            bc(ii,jj) = 1
          else if (i.eq.size(1)) then
            bc(ii,jj) = 1
          else if (j.eq.1) then
            bc(ii,jj) = 1
          else if (j.eq.size(2)) then
            bc(ii,jj) = 2
          else
            bc(ii,jj) = 0
          endif
c          bc(ii,jj) = 0
        end do
      end do
      call ga_update_ghosts(g_bc)
c
c     Create initial distribution of density and velocities
c
      pi = 4.0d00*atan(1.0d00)
      rtot = 0.0d00
      do jj = width(2) + 1, dims(2) - width(2)
        j = jj - width(2) - 1 + lo(2)
        do ii = width(1) + 1, dims(1) - width(1)
          i = ii - width(1) - 1 + lo(1)
          fld(ii,jj,1) = rho0 + 0.0*cos(2.0d00*pi*dble(j-1)
     +                 / dble(size(2)-1))
          fld(ii,jj,2) = ux0
          fld(ii,jj,3) = uy0
          fld(ii,jj,4) = rho0*rgas*tmprtr0/(1.0d00-b_vdw*rho0)
     +                 - a_vdw*rho0**2
          fld(ii,jj,5) = 6.0d00*viscosity/(cspd**2*delta_t) + 0.5d00
          rtot = rtot + fld(ii,jj,1)
        end do
      end do
c
c     initialize lattice parameters
c
      call initpar
c
c     evaluate equilibrium distribution
c
      call equil(g_fg,  fg,  ld_fg1,  ld_fg2,
     +           g_fld, fld, ld_fld1, ld_fld2,
     +           g_bc,  bc,  ld_bc1)
      do jj = width(2)+1, dims(2) - width(2)
        do j = 1, 9
          do ii = width(1)+1, dims(1) - width(1)
            fg(ii,jj,j) = fg(ii,jj,j+9)
          end do
        end do
      end do
      call properties(g_fg,  fg,  ld_fg1,  ld_fg2,
     +                g_fld, fld, ld_fld1, ld_fld2,
     +                g_bc,  bc,  ld_bc1)
c        
      return
      end
