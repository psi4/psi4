#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine timestep(g_fg,  fg,  ld_fg1,  ld_fg2,
     +                    g_fld, fld, ld_fld1, ld_fld2,
     +                    g_bc,  bc,  ld_bc1)
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
#define MP_DEFINES_ONLY
#include "mp3.fh"
#undef MP_DEFINES_ONLY
c
c   Advance simulation one timestep
c
      integer ii, jj, i, j, k, ix, iy
      double precision start1, start2
c
c   Make backup copies of distribution and update ghost cells
c
      start1 = MP_TIMER()
      do jj = width(2) + 1, dims(2) - width(2)
        do i = 1, 9
          j = i + 18 
          do ii = width(1) + 1, dims(1) - width(1)
            fg(ii,jj,j) = fg(ii,jj,i)
          end do
        end do
      end do
      start2 = MP_TIMER()
      call ga_mask_sync(.false.,.false.)
      call ga_update_ghosts(g_fg)
      tstats(2) = tstats(2) + MP_TIMER() - start2
c
c   Perform streaming operation
c
      do jj = width(2) + 1, dims(2) - width(2)
        do ii = width(1) + 1, dims(1) - width(1)
          if (bc(ii,jj).eq.0) then
            do i = 2, 9
              j = i + 18 
              k = i - 1
              ix = nint(ei(k,1))
              iy = nint(ei(k,2))
              fg(ii,jj,i) = fg(ii-ix,jj-iy,j)
            end do
          else
            call get_patch(g_fg,  fg,  ld_fg1,  ld_fg2,
     +                     g_fld, fld, ld_fld1, ld_fld2,
     +                     g_bc,  bc,  ld_bc1, ii, jj)
            do i = 2, 9
              k = i - 1
              ix = nint(ei(k,1))
              iy = nint(ei(k,2))
              fg(ii,jj,i) = fgp(-ix,-iy,i)
            end do
            
          endif
        end do
      end do
c
c   Update properties
c
      call properties(g_fg,  fg,  ld_fg1,  ld_fg2,
     +                g_fld, fld, ld_fld1, ld_fld2,
     +                g_bc,  bc,  ld_bc1)
c
c   Perform relaxation
c
      call equil(g_fg,  fg,  ld_fg1,  ld_fg2,
     +           g_fld, fld, ld_fld1, ld_fld2,
     +           g_bc,  bc,  ld_bc1)
      do jj = width(2) + 1, dims(2) - width(2)
        do i = 1, 9
          do ii = width(1) + 1, dims(1) - width(1)
            if (fld(ii,jj,5).gt.0.0d00) then
              fg(ii,jj,i) = fg(ii,jj,i+18)
     +                    - (fg(ii,jj,i+18)-fg(ii,jj,i+9))
     +                    / fld(ii,jj,5)
            endif
c            if (bc(ii,jj).eq.2) then
c              fg(ii,jj,i) = fg(ii,jj,i+9)
c            endif
          end do
        end do
      end do
c
      tstats(1) = tstats(1) + MP_TIMER() - start1
c
      return
      end
