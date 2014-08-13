      program he
C$Id: mc.f,v 1.2 1995-02-02 23:24:16 d3g681 Exp $
      implicit double precision (a-h, o-z)
      include 'msgtypesf.h'
      parameter (maxpt=512)
      dimension x(6,maxpt), psix(maxpt), ex(maxpt)
      data esq/0.0d0/
c
c     Monte Carlo Test Program ... evaluates the expectation 
c     value of the energy for a simple wavefunction for helium.
c
c     Initalize the message passing environment
c
      call pbeginf
      call evon
c
c     process zero reads in input and then broadcasts to others
c
      if (nodeid().eq.0) then
        write(6,1)
1       format(' He atom variational monte carlo '//
     &         ' Input neq and nstep '$)
c       call flush(6)
        read(5,*) neq, nstep
        nstep = ( (nstep+99)/100 ) * 100
        neq = ( (neq+99)/100 ) * 100
      endif
      call brdcst(1+MSGINT, neq, mitob(1), 0)
      call brdcst(1+MSGINT, nstep, mitob(1), 0)
c
c     Divide up the work ... each process does a subset of the points
c     or configurations ... make sure that total is indeed maxpt
c
      nproc = nnodes()
      npoint = maxpt / nproc
      nremain = maxpt - npoint*nproc
      if (nodeid().lt.nremain) npoint = npoint + 1
c
c     Actually do the work. Routine init generates the intial points
c     in a 2x2x2 cube. Equilibriate for neq moves then compute averages
c     for nstep moves.
c
      rjunk = timer()
      call init(npoint, x, psix)
      call mover(x, psix, ex, e, esq, npoint, neq)
      call mover(x, psix, ex, e, esq, npoint, nstep)
c
c     Sum the results from all the processors
c
      call dgop(2+MSGDBL, e, 1, '+')
      call dgop(3+MSGDBL, esq, 1, '+')
c
c     Write out the results before terminating
c
      if (nodeid().eq.0) then
         e = e / dble(nproc)
         esq = esq / dble(nproc)
         err = sqrt((esq-e*e) / dble(nproc*nstep/100))
         used = timer()
         write(6,2) e, err, used, nproc
 2       format(' energy =',f10.6,' +/-',f9.6,': used =',f7.2,' secs',
     $        ': nproc =',i3)
      endif
c
      call pend
      call fexit
c
      end
      subroutine mover(x, psix, ex, e, esq, npoint, nstep)
      implicit double precision (a-h, o-z)
      dimension x(6,npoint), psix(npoint), ex(npoint)
c
c     move the set of points nstep times accumulating averages
c     for the energy and square of the energy
c     
      e = 0.0d0
      esq = 0.0d0
      eb = 0.0d0
      do 10 istep = 1, nstep
c
c     sample a new set of points
c
         do 20 ipoint = 1, npoint
            call sample(x(1, ipoint), psix(ipoint), ex(ipoint))
 20      continue
c
c     accumulate average(s)
c
         do 30 ipoint = 1, npoint
            eb = eb + ex(ipoint)
 30      continue
c
c     block averages every 100 moves to reduce statistical correlation
c
         if (mod(istep,100).eq.0) then
            eb = eb / dble(npoint*100)
            e = e + eb
            esq = esq + eb*eb
            eb = 0.0d0
         endif
 10   continue
c
c     compute final averages
c
      e = e / dble(nstep/100)
      esq = esq / dble(nstep/100)
c
      end
      subroutine sample(x, psix, ex)
      implicit double precision (a-h, o-z)
      dimension x(6), xnew(6)
c
c     sample a new point ... i.e. move current point by a
c     random amount and accept the move according to the
c     ratio of the square of the wavefunction at the new
c     point and the old point.
c
c     generate trial point and info at the new point
c
      do 10 i = 1,6
         xnew(i) = x(i) + (drand48(0)-0.5d0)*0.3d0
 10   continue
      call rvals(xnew, r1, r2, r12, r1dr2)
      psinew = psi(r1, r2, r12)
c
c     accept or reject the move
c
      prob = min((psinew / psix)**2, 1.0d0)
      if (prob .gt. drand48(0)) then
         do 20 i = 1,6
            x(i) = xnew(i)
 20      continue
         psix = psinew
      else
         call rvals(x, r1, r2, r12, r1dr2)
      endif
      ex = elocal(r1, r2, r12, r1dr2)
c     
      end
      subroutine rvals(x, r1, r2, r12, r1dr2)
      implicit double precision (a-h, o-z)
      dimension x(6)
c
c     compute required distances etc.
c
      r1 = dsqrt(x(1)**2 + x(2)**2 + x(3)**2)
      r2 = dsqrt(x(4)**2 + x(5)**2 + x(6)**2)
      r12 = dsqrt((x(1)-x(4))**2 + (x(2)-x(5))**2 + (x(3)-x(6))**2)
      r1dr2 = x(1)*x(4) + x(2)*x(5) + x(3)*x(6)
c
      end
      double precision function psi(r1, r2, r12)
      implicit double precision (a-h, o-z)
c
c     compute value of the trial wavefunction
c
      psi = dexp(-2.0d0*(r1+r2)) * (1.0d0 + 0.5d0*r12)
c
      end
      double precision function elocal(r1, r2, r12, r1dr2)
      implicit double precision (a-h, o-z)
c
c     compute local energy = (H psi) / psi 
c
      f = r12*(1.0d0 + 0.5d0*r12)
      g = 0.5d0*r12 + r1 +r2 - r1dr2*(1.0d0/r1 + 1.0d0/r2)
      elocal = -4.0d0 + g / f
c
      end
      subroutine init(npoint, x, psix)
      implicit double precision (a-h, o-z)
      dimension x(6,npoint), psix(npoint)
c
c     distribute points in a 2x2x2 cube.
c
c     in parallel version make random no. seed depend on the
c     process number ... for a production code should be more
c     sophisticated than this.
c
      call srand48(2*nodeid()+1)
c      
      do 10 ipoint = 1,npoint
         do 20 i = 1,6
            x(i,ipoint) = (drand48(0) - 0.5d0) * 2.0d0
 20      continue
         call rvals(x(1,ipoint), r1, r2, r12, r1dr2)
         psix(ipoint) = psi(r1, r2, r12)
 10   continue
c
      end
