#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      subroutine boltzmann(g_fg,  fg,  ld_fg1,  ld_fg2,
     +                     g_fld, fld, ld_fld1, ld_fld2,
     +                     g_bc,  bc,  ld_bc1)
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
      integer i, me, nprocs
      me = ga_nodeid()
      nprocs = ga_nnodes()
c
c     initialize problem
c
      call setup(g_fg,  fg,  ld_fg1,  ld_fg2,
     +           g_fld, fld, ld_fld1, ld_fld2,
     +           g_bc,  bc,  ld_bc1)
c
      if (me.eq.0) then
        write(6,*) 'Total mass is ',rtot
        write(6,*) 'Total x-momentum is ',uxtot
        write(6,*) 'Total y-momentum is ',uytot
      endif
      do i = 1, nsteps
        call timestep(g_fg,  fg,  ld_fg1,  ld_fg2,
     +                g_fld, fld, ld_fld1, ld_fld2,
     +                g_bc,  bc,  ld_bc1)
        if (mod(i,200).eq.0) then
          if (me.eq.0) then
            write(6,*) 'Completed step ',i
            write(6,*) 'Total density is ',rtot
            write(6,*) 'Total x-momentum is ',uxtot
            write(6,*) 'Total y-momentum is ',uytot
          endif
          call vorticity(g_fg,  fg,  ld_fg1,  ld_fg2,
     +                   g_fld, fld, ld_fld1, ld_fld2,
     +                   g_bc,  bc,  ld_bc1)
          call printdat(g_fg,  fg,  ld_fg1,  ld_fg2,
     +                  g_fld, fld, ld_fld1, ld_fld2,
     +                  g_bc,  bc,  ld_bc1)
        endif
      end do 
      call vorticity(g_fg,  fg,  ld_fg1,  ld_fg2,
     +               g_fld, fld, ld_fld1, ld_fld2,
     +               g_bc,  bc,  ld_bc1)
      call printdat(g_fg,  fg,  ld_fg1,  ld_fg2,
     +              g_fld, fld, ld_fld1, ld_fld2,
     +              g_bc,  bc,  ld_bc1)
c
      call ga_dgop(MT_DBL,tstats,50,'+')
      if (me.eq.0) write(6,100) tstats(1)/dble(nprocs)
      if (me.eq.0) write(6,200) tstats(2)/dble(nprocs)
c
  100 format('Average time in lattice Boltzmann updates :',f12.4)
  200 format('Average time in ghost cell updates        :',f12.4)
      return
      end
