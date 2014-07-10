c     
c     md example program due to dieter heermann
c     restructured and pdf added by rjh december 1990
c
c     the message passing version distributes the computation
c     of the forces (order npart**2) over the processes, assuming
c     that each process has all of the co-ordinates.
c     A global add then gives each process all the information
c     needed to compute the next update.
c     None of the order npart work has been parallelized so that
c     will begin to dominate on many processors.
c     
      implicit double precision (a-h, o-z)
c      parameter (mm = 13, lenpdf = 256)
c      parameter (mm = 8, lenpdf = 256)
c      parameter (mm = 6, lenpdf = 256)
C$Id: md.f,v 1.2 1995-02-02 23:24:18 d3g681 Exp $
      parameter (mm = 4, lenpdf = 256)
c      parameter (mm = 3, lenpdf = 256)
      parameter (npart = 4*mm*mm*mm)
      parameter (maxint = npart*150)
      include 'msgtypesf.h'
c     
      dimension  x(1:npart,1:3), vh(1:npart,1:3),f(1:npart,1:3),
     &     pdf(1:lenpdf+1), times(8), inter(2,maxint)
      data times/8*0.0d0/
c
c     initalize message passing environment
c
      call pbeginf
      me = nodeid()
c     
c     parameter definition (density, volume, temperature ...)
c     
      den = 0.83134d0
      side = (dble(npart) / den)**0.3333333d0
      tref = 0.722d0
      rcoff = min(3.5d0, side/2.0d0)
c     islow = 1 to match the rather large original timestep
      islow = 4
      h = 0.064d0 / islow
      irep = 50
      istop = 400
      iprint = 10
      ineigh = 10 * (1 + islow/4)
      movemx = 800
      delpdf = 0.5d0*side/lenpdf
      rdelp = 1.0d0 / delpdf
c
      if (me.eq.0)write(6,1) npart,side,rcoff,tref,h,delpdf,irep,istop,
     &     iprint,ineigh,movemx
 1    format(' molecular dynamics simulation example program'/
     &     ' ---------------------------------------------'//
     &     ' number of particles is ............ ',i6/
     &     ' side length of the box is ......... ',f13.6/
     &     ' cut off is ........................ ',f13.6/
     &     ' reduced temperature is ............ ',f13.6/
     &     ' basic time step is ................ ',f13.6/
     &     ' pdf sampling interval ............. ',f13.6/
     &     ' temperature scale interval ........ ',i6/
     &     ' stop scaling at move .............. ',i6/
     &     ' print interval .................... ',i6/
     &     ' update neighbor list every ........ ',i6, ' steps'/
     &     ' total no. of steps ................ ',i6)
c     call flush(6)
c     
      a = side / dble(mm)
      sideh = side * 0.5d0
      hsq = h * h
      hsq2 = hsq * 0.5d0
      npartm = npart - 1
      rcoffs = rcoff * rcoff
      tscale = 16.0d0 / (1.0d0 * npart - 1.0d0)
      vaver = 1.13d0 * sqrt(tref / 24.0d0)
      ekinavg = 0.0d0
c     
c     generate fcc lattice for atoms inside box
c     
      rjunk = timer()
      call fcc(x, npart, mm, a)
      times(1) = times(1) + timer()
c     
c     initialise velocites and forces (which are zero in fcc positions)
c     
      call mxwell(vh,3*npart,h,tref)
      call dfill(3*npart, 0.0d0, f, 1)
      times(2) = times(2) + timer()
c     
c     start of md. 
c     
      if (me.eq.0) write(6,3)
 3    format(//1x,'   i  ','     ke    ','      pe     ','      e     ',
     &     '    temp   ','   pres   ','   vel    ','  rp'/
     &     1x,' -----','  ----------','  ----------','  ----------',
     &        '  --------','  --------','  --------','  ----')
c     call flush(6)
c     
      do 200 move = 1,movemx
         if (move.eq.1 .or. move.eq.(istop+1))
     $        call dfill(lenpdf, 0.0d0, pdf, 1)
c     
c     move the particles and partially update velocities
c     
         call domove(3*npart,x,vh,f,side)
         times(3) = times(3) + timer()
c     
c     compute forces in the new positions and accumulate the pdf
c     virial and potential energy. Have to get the full forces
c     on each node, hence the dgop (global add).
c     
         if (mod(move-1,ineigh).eq.0) then
            call neigh(npart, x, side, rcoff, ninter, inter, maxint,
     $           pdf, rdelp, lenpdf)
            times(8) = times(8) + timer()
         endif
         call forces(npart, x, f, vir, epot, side, rcoff, ninter, inter)
         call dgop(1+MSGDBL, f, 3*npart, '+')
         times(4) = times(4) + timer()
c
c     scale forces, complete update of velocites and compute k.e.
c
         call mkekin(npart,f,vh,hsq2,hsq,ekin)
         ekinavg = ekinavg + ekin
         times(5) = times(5) + timer()
c
c     average the velocity and temperature scale if desired
c
         if ((move.le.istop) .and. (mod(move, irep).eq.0)) then
            call velavg(npart, vh, vaver, count, vel, h)
            sc = sqrt( tref / (tscale * ekinavg / irep) )
            call dscal(3*npart, sc, vh, 1)
            ekinavg = 0.0d0
         endif
         times(6) = times(6) + timer()
c
c     printout information if desired ... have to do global
c     sum to get full potential energy and virial
c
         if (mod(move, iprint) .eq. 0) then
            call velavg(npart, vh, vaver, count, vel, h)
            call dgop(2+MSGDBL, epot, 1, '+')
            call dgop(2+MSGDBL, vir, 1, '+')
            if (me.eq.0)
     $           call prnout(move, ekin, epot, tscale, vir, vel, count,
     $           npart, den)
            times(7) = times(7) + timer()
         endif
c     
 200  continue
c     
c     print out the pdf at the end of the calculation
c     have first to get contribution from everyone with global add
c     
      call dgop(2+MSGDBL, pdf, lenpdf, '+')
      if (me.eq.0)
     $     call prnpdf(lenpdf, pdf, side, delpdf, npart, movemx-istop,
     $     ineigh)
      times(7) = times(7) + timer()
      if (me.eq.0) write(6,431) nnodes(),(times(i),i=1,8)
 431  format('  nproc    geom  mxwell  domove  forces    ekin  velscl',
     &     '   print   neigh'
     &     /1x,i6,8f8.2)
c
      if (me.eq.0) call stats
      call pend
      call fexit
      end
      subroutine mxwell(vh,n3,h,tref)
      implicit double precision (a-h, o-z)
      dimension vh(1:n3)
c
c sample maxwell distribution at temperature tref
c
c alliant 3
      iseed = 4711
      ujunk = drand48(iseed)
      iseed = 0
      npart = n3/3
      iof1 = npart
      iof2 = npart*2
      tscale = 16.0d0 / (1.0d0 * npart - 1.0d0)
      do 10 i = 1,n3,2
c
c cray 2
c1         u1 = ranf()
c          u2 = ranf()
c alliant 2
1         u1 = drand48(iseed)
          u2 = drand48(iseed)
          v1 = 2.0d0 * u1 - 1.0d0
          v2 = 2.0d0 * u2 - 1.0d0
          s = v1*v1 + v2*v2
          if (s.ge.1.0) goto 1
c
          r = sqrt(-2.0d0*dlog(s)/s)
          vh(i) = v1 * r
          vh(i+1) = v2 * r
10    continue
c
      ekin = 0.0d0
      sp = 0.0d0
      do 20 i = 1,npart
          sp = sp + vh(i)
20    continue
      sp = sp / npart
      do 21 i = 1,npart
          vh(i) = vh(i) - sp
          ekin = ekin + vh(i)*vh(i)
21    continue
      sp = 0.0d0
      do 22 i = iof1 + 1,iof2
          sp = sp + vh(i)
22    continue
      sp = sp / npart
      do 23 i = iof1 + 1,iof2
          vh(i) = vh(i) - sp
          ekin = ekin + vh(i)*vh(i)
23    continue
      sp = 0.0d0
      do 24 i = iof2 + 1,n3
          sp = sp + vh(i)
24    continue
      sp = sp / npart
      do 25 i = iof2 + 1,n3
          vh(i) = vh(i) - sp
          ekin = ekin + vh(i)*vh(i)
25    continue
      ts = tscale * ekin
      sc = h * sqrt(tref/ts)
      do 30 i = 1,n3
          vh(i) = vh(i) * sc
30    continue
c
      end
      subroutine domove(n3,x,vh,f,side)
      implicit double precision (a-h, o-z)
      dimension x(n3),vh(n3),f(n3)
c
c     move particles
c
      do 10 i = 1,n3
         x(i) = x(i) + vh(i) + f(i)
c     periodic boundary conditions
         if (x(i).lt.0.0d0) x(i) = x(i) + side
         if (x(i).gt.side) x(i) = x(i) - side
c     partial velocity updates
         vh(i) = vh(i) + f(i)
c     initialise forces for next iteration
         f(i) = 0.0d0
 10   continue
c
      end    
      subroutine mkekin(npart,f,vh,hsq2,hsq,ekin)
      implicit double precision (a-h, o-z)
      dimension f(1:npart,3),vh(1:npart,3)
c
c     scale forces, update velocites and compute k.e.
c
      sum = 0.0d0
      do 10 ix = 1,3
         do 20 i = 1,npart
            f(i,ix) = f(i,ix) * hsq2
            vold = vh(i,ix)
            vh(i,ix) = vh(i,ix) + f(i,ix)
            sum = sum + vh(i,ix) * vh(i,ix)
 20      continue
 10   continue
      ekin = sum / hsq
c
      end
      subroutine fcc(x, npart, mm, a)
      implicit double precision (a-h, o-z)
      dimension x(1:npart, 3)
c     
c     generate fcc lattice for atoms inside box
c     
      ijk = 0
      do 10 lg = 0,1
         do 11 i = 0,mm-1
            do 12 j = 0,mm-1
               do 13 k = 0,mm-1
                  ijk = ijk + 1
                  x(ijk,1) = i * a + lg * a * 0.5d0
                  x(ijk,2) = j * a + lg * a * 0.5d0
                  x(ijk,3) = k * a
 13            continue
 12         continue
 11      continue
 10   continue
      do 20 lg = 1,2
         do 21 i = 0,mm-1
            do 22 j = 0,mm-1
               do 23 k = 0,mm-1
                  ijk = ijk + 1
                  x(ijk,1) = i * a + (2-lg) * a * 0.5d0
                  x(ijk,2) = j * a + (lg-1) * a * 0.5d0
                  x(ijk,3) = k * a + a * 0.5d0
 23            continue
 22         continue
 21      continue
20    continue
c
      end
      subroutine dfill(n,val,a,ia)
      implicit double precision (a-h, o-z)
      dimension a(*)
c
c     initialise double precision array to scalar value
c
      do 10 i = 1,(n-1)*ia+1,ia
         a(i) = val
 10   continue
c
      end
      subroutine prnout(move, ekin, epot, tscale, vir, vel, count,
     $     npart, den)
      implicit double precision (a-h, o-z)
c
c     printout interesting (?) information at current timestep
c
      ek = 24.0d0 * ekin
      epot = 4.0d0 * epot
      etot = ek + epot
      temp = tscale * ekin
      pres = den * 16.0d0 * (ekin - vir) / npart
      vel = vel / npart
      rp = (count / dble(npart)) * 100.0d0
      write(6,2) move,ek,epot,etot,temp,pres,vel,rp
 2    format(1x,i6,3f12.4,f10.4,f10.4,f10.4,f6.1)
c     call flush(6)
c
      end
      subroutine velavg(npart, vh, vaver, count, vel, h)
      implicit double precision (a-h, o-z)
      dimension vh(npart, 3)
c     
c     compute average velocity
c     
      vaverh = vaver*h
      vel = 0.0d0
      count = 0.0d0
      do 10 i = 1,npart
         sq = sqrt(vh(i,1)**2 + vh(i,2)**2 + vh(i,3)**2)
         if (sq.gt.vaverh) count = count + 1.0d0
         vel = vel + sq
 10   continue
      vel = vel / h
c     
      end
      subroutine prnpdf(lenpdf, pdf, side, delpdf, npart, nmove, ineigh)
      implicit double precision (a-h, o-z)
      dimension pdf(lenpdf)
c     
c     final scaling and printout of the pdf
c     
      write(6,1)
 1    format(/' pair distribution function written to file pdf.dat'/)
      open(1, file='pdf.dat', form='formatted', status='unknown',
     $     err=999)
c     
      coord = 0.0d0
      volfac = side*side*side / (4.0d0*delpdf*delpdf*delpdf*3.141593d0)
      facnn = 2.0d0 / dble(npart*nmove/ineigh)
      facn = 1.0d0 / dble(npart)
      do 10 i = 1,lenpdf
         ri = dble(i)
         grfac = volfac / (ri*ri)
         func = pdf(i) * facnn
         coord = coord + func
         pdf(i) = func * grfac * facn
         write(1,2) dble(i)*delpdf,pdf(i),coord
 10   continue
 2    format(1x,f7.3,f13.6,4x,f9.2)
      close(1)
      return
c
 999  write(6,*) ' error opening pdf.dat'
      call parerr(999)
c     
      end
      subroutine neigh(npart, x, side, rcoff, ninter, inter, maxint,
     $     pdf, rdelp, lenpdf)
      implicit double precision (a-h, o-z)
      dimension x(npart, 3), inter(2,maxint), pdf(lenpdf)
c     
c     Form my part of the neighbour list and also compute the pair
c     distribution function.
c     
c     npart = no. of particles
c     x(,)  = coords
c     side  = side of box
c     rcoff = cutoff for force
c     ninter= returns no. of interactions
c     inter(,) = returns interactions
c     maxint   = size of inter
c     
      me = nodeid()
      nproc = nnodes()
c     
      sideh = 0.5d0*side
      rcoffs = (rcoff*1.2d0)**2
      ninter = 0
c
c     Get better work distribution by having the same
c     processor handle particles (i) and (npart-i)
c     Note that assume that npart is even.
c
      do 270 ii = me+1, npart/2, nproc
         do 275 icase = 1, 2
            if (icase .eq. 1) then
               i = ii
            else
               i = npart - ii
            endif
            xi = x(i,1)
            yi = x(i,2)
            zi = x(i,3)
            do 280 j = i+1,npart
               ij = ij + 1
               xx = xi - x(j,1)
               yy = yi - x(j,2)
               zz = zi - x(j,3)
               if (xx .lt. -sideh) xx = xx + side
               if (xx .gt.  sideh) xx = xx - side
               if (yy .lt. -sideh) yy = yy + side
               if (yy .gt.  sideh) yy = yy - side
               if (zz .lt. -sideh) zz = zz + side
               if (zz .gt.  sideh) zz = zz - side
               rd = xx*xx + yy*yy + zz*zz
               ipdf = min(sqrt(rd),sideh) * rdelp + 1
               pdf(ipdf) = pdf(ipdf) + 1.0d0
               if (rd .le. rcoffs) then
                  ninter = ninter + 1
                  if (ninter.gt.maxint) then
                     write(6,*) ' too many interactions ', ninter
                     call parerr(1)
                  endif
                  inter(1,ninter) = i
                  inter(2,ninter) = j
               endif
 280        continue
 275     continue
 270  continue
c     
c$$$      ij = ninter
c$$$      call igop(99, ij, 1, '+')
c$$$      if (me .eq. 0) then
c$$$      write(6,*) ' No. of interactions per particle = ',
c$$$     $        ij/npart
c$$$      endif
c     
      end
      subroutine forces(npart, x, f, vir, epot, side, rcoff, 
     $     ninter, inter)
      implicit double precision (a-h, o-z)
      logical oshift
      parameter (oshift = .true.)
      dimension x(npart, 3), f(npart, 3), inter(2, ninter)
c     
c     compute forces driven by the neighbour list
c     
      vir = 0.0d0
      epot = 0.0d0
      sideh = 0.5d0*side
      rcoffs = rcoff*rcoff
c
c     for shifted potential ... set oshift true to enable
c
      if (oshift) then
         rc6  = 1.0d0 / rcoff**6
         rc12 = 1.0d0 / rcoff**12
         ecut = rc12 - rc6
         fcut = (rc12 - 0.5d0*rc6)/rcoff
         efcut = 12.0d0 * fcut
      endif
c
      do 10 ij = 1, ninter
         i = inter(1,ij)
         j = inter(2,ij)
         xx = x(i,1) - x(j,1)
         yy = x(i,2) - x(j,2)
         zz = x(i,3) - x(j,3)
         if (xx .lt. -sideh) xx = xx + side
         if (xx .gt.  sideh) xx = xx - side
         if (yy .lt. -sideh) yy = yy + side
         if (yy .gt.  sideh) yy = yy - side
         if (zz .lt. -sideh) zz = zz + side
         if (zz .gt.  sideh) zz = zz - side
         rd = xx*xx + yy*yy + zz*zz
         if (rd .le. rcoffs) then
            rrd = 1.0d0/rd
            rrd2 = rrd*rrd
            rrd3 = rrd2*rrd
            rrd4 = rrd2*rrd2
            rrd6 = rrd2*rrd4
            rrd7 = rrd6*rrd
            epot = epot + (rrd6 - rrd3)
            r148 = rrd7 - 0.5d0*rrd4
            if (oshift) then
               r = sqrt(rd)
               epot = epot - ecut + efcut*(r-rcoff)
               r148 = r148 - fcut / r
            endif
            vir = vir - rd*r148
            forcex = xx * r148
            f(i,1) = f(i,1) + forcex
            f(j,1) = f(j,1) - forcex
            forcey = yy * r148
            f(i,2) = f(i,2) + forcey
            f(j,2) = f(j,2) - forcey
            forcez = zz * r148
            f(i,3) = f(i,3) + forcez
            f(j,3) = f(j,3) - forcez
         endif
 10   continue
c     
      end
