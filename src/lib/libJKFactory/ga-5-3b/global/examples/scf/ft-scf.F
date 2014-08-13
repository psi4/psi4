#if HAVE_CONFIG_H
#   include "config.fh"
#endif
      program scf
C$Id: ft-scf.F,v 1.1 2005-07-07 17:41:07 vinod Exp $
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"

#define USE_TRANSFORM 1
      integer heap, stack
      data tinit, tonel, ttwoel, tdiag, tdens, tprint /6*0.0d0/
      data eone, etwo, energy, deltad /4*0.0d0/
c     
c     initalize the parallel message passing environment
c     
#ifdef MPI
      integer ierr
      call mpi_init(ierr)
#else
      call pbeginf
#endif
      call ga_set_spare_procs(2);
      call ga_initialize()
c
c   Allocate memory
c
      heap = 2000000
      stack = 2000000
      if (.not.ma_init(MT_DBL, stack, heap))
     +   call ga_error("ma_init failed",-1)
      me = ga_nodeid()
      nproc = ga_nnodes()
c     
c     initialize a bunch of stuff and initial density matrix
c     
      rjunk = timer()
c
c     get input from file be.inpt
c
      call input
c
c     create and allocate global arrays
c
      call setarrays
      write(6,*) 'done with setarrays'
      call ininrm
      write(6,*) '1-done with setarrays'
c
c     create initial guess for density matrix by using single atom
c     densities
c
      call denges
      tinit = timer()
#if USE_TRANSFORM
c
c     make initial orthogonal orbital set for solution method using
c     similarity transform
c
      call makeob
#endif
c     
c     make info for sparsity test
c     
         call makesz(schwmax)
c     
c     iterate
      call iterate(10)
c     
      call ga_igop(6, icut1, 1, '+')
      call ga_igop(7, icut2, 1, '+')
      call ga_igop(8, icut3, 1, '+')
      if (me.eq.0) then
c     
c     print out timing information
c     
         call prnfin(energy)
         write(6,1) tinit, tonel, ttwoel, tdiag, tdens, tprint,
     $        nproc
 1       format(/'   init   onel  twoel   diag   dens   print  ncpu'/
     $        '  ------ ------ ------ ------ ------ ------ ------'/
     $        1x, 6f7.2, i7/)
c     
c     print out information on # integrals evaulated each iteration
c     
         nints = nnbfn*(nnbfn+1)/2
         nints = nbfn**4
         frac  = dble(icut3)/dble(nints)
         write(6,2) icut1, icut2, icut3, nints, frac
 2       format(/'       No. of integrals screened or computed '
     $        /'       -------------------------------------'/
     $        /1x,'   #ij test   #kl test   #compute     #total',
     $        '  fraction',
     $        /1x,'  ---------  ---------  ---------  ---------',
     $        '  --------',
     $        /1x,4(2x,i9),f9.3)
         call stats
      endif
c     
      call closearrays
      call ga_terminate
#ifdef MPI
      call mpi_finalize
#else
      call pend
#endif
c     
      end
c
      subroutine iterate(niter)
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      do 10 iter = 1, niter
#ifdef CHECKPOINT_SCF
         call ga_checkpoint_arrays(arraylist,numarrays)
#endif
c
c     make the one particle contribution to the fock matrix
c     and get the partial contribution to the energy
c     
         call oneel(schwmax, eone)
         tonel = tonel + timer()
c     
c     compute the two particle contributions to the fock matrix and
c     get the total energy.
c     
         call twoel(schwmax, etwo)
         ttwoel = ttwoel + timer()
c     
c     Diagonalize the fock matrix. The diagonalizers used in this
c     subroutine are actually sequential, not parallel.
c
         call diagon(tester,iter)
         tdiag = tdiag + timer()
c     
c     make the new density matrix in g_work from orbitals in g_orbs,
c     compute the norm of the change in the density matrix and
c     then update the density matrix in g_dens with damping.
c     
         call makden
         deltad = dendif()
         if (iter.eq.1) then
            scale = 0.0d0
         else if (iter .le. 5) then
            if (nbfn .gt. 60) then
               scale = 0.5d0
            else
               scale = 0.0d0
            endif
         else
            scale = 0.0d0
         endif
         call damp(scale)
         tdens = tdens + timer()
c     
c     add up energy and print out convergence information
c     
         if (me.eq.0) then
           energy = enrep + eone + etwo
           call prnout(iter, energy, deltad, tester)
           tprint = tprint + timer()
         endif
c     
c     if converged then exit iteration loop
c     
         if (deltad .lt. tol) goto 20
 10   continue
      if(me.eq.0)
     $     write(6,*) ' SCF failed to converge in ', niter, ' iters'
c     
c     finished ... print out eigenvalues and occupied orbitals
c     
 20   continue
      end

      subroutine makesz(schwmax)
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      dimension work(ichunk,ichunk)
      integer lo(2),hi(2),i,j,iloc,jloc,ld
      logical dotask
c
c     schwarz(ij) = (ij|ij) for sparsity test
c
      icut1 = 0
      icut2 = 0
      icut3 = 0
c
      call ga_zero(g_schwarz)
      call ga_zero(g_counter)
      schwmax = 0.0d0
      dotask = next_chunk(lo,hi)
      ld = ichunk
      do while (dotask)
        do i = lo(1), hi(1)
          iloc = i - lo(1) + 1
          do j = lo(2), hi(2)
            jloc = j - lo(2) + 1
            call g(gg,i,j,i,j)
            work(iloc,jloc) = sqrt(gg)
            schwmax = max(schwmax, work(iloc,jloc))
          end do
        end do
        call nga_put(g_schwarz,lo,hi,work,ld)
        dotask = next_chunk(lo,hi)
      end do
      call ga_dgop(11,schwmax,1,'max')
c
      return
      end
c
      subroutine ininrm
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c
c     write a little welcome message
c
      if (ga_nodeid().eq.0) write(6,1) natom, nocc, nbfn, tol
1     format(/' Example Direct Self Consistent Field Program '/
     $        ' -------------------------------------------- '//
     $        ' no. of atoms ............... ',i3/
     $        ' no. of occupied orbitals ... ',i3/
     $        ' no. of basis functions ..... ',i3/
     $        ' convergence threshold ...... ',d9.2//)
c
c     generate normalisation coefficients for the basis functions
c     and the index array iky
c
      do 10 i = 1, nbfn
         iky(i) = i*(i-1)/2
 10   continue
c
      do 20 i = 1, nbfn
         rnorm(i) = (expnt(i)*2.0d0/pi)**0.75d0
 20   continue
c
c     initialize common for computing f0
c
      call setfm
c
      end
      double precision function h(i,j)
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
cvd$r novector
cvd$r noconcur
c
c     generate the one particle hamiltonian matrix element
c     over the normalized primitive 1s functions i and j
c
      f0val = 0.0d0
      sum = 0.0d0
      rab2 = (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2
      facij = expnt(i)*expnt(j)/(expnt(i)+expnt(j))
      expij = exprjh(-facij*rab2)
      repij = (2.0d0*pi/(expnt(i)+expnt(j))) * expij
c
c     first do the nuclear attraction integrals
c
      do 10 iat = 1, natom
         xp = (x(i)*expnt(i) + x(j)*expnt(j))/(expnt(i)+expnt(j))
         yp = (y(i)*expnt(i) + y(j)*expnt(j))/(expnt(i)+expnt(j))
         zp = (z(i)*expnt(i) + z(j)*expnt(j))/(expnt(i)+expnt(j))
         rpc2 = (xp-ax(iat))**2 + (yp-ay(iat))**2 + (zp-az(iat))**2
c
         call f0(f0val, (expnt(i)+expnt(j))*rpc2)
         sum = sum - repij * q(iat) * f0val
 10   continue
c
c     add on the kinetic energy term
c
      sum = sum + facij*(3.0d0-2.0d0*facij*rab2) *
     $     (pi/(expnt(i)+expnt(j)))**1.5d0 * expij
c
c     finally multiply by the normalization constants
c
      h = sum * rnorm(i) * rnorm(j)
c
      end
      double precision function s(i,j)
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c
c     generate the overlap matrix element between the normalized
c     primitve gaussian 1s functions i and j
c
      rab2 = (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2
      facij = expnt(i)*expnt(j)/(expnt(i)+expnt(j))
      s = (pi/(expnt(i)+expnt(j)))**1.5d0 * exprjh(-facij*rab2) *
     $     rnorm(i)*rnorm(j)
c
      end
      subroutine makden
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c      dimension work(maxnbfn,maxnbfn), torbs(maxnbfn,maxnbfn)
      dimension work(ichunk,ichunk), orbsi(ichunk,maxnbfn)
      dimension orbsj(ichunk,maxnbfn)
      integer lo(2), hi(2), tlo(2), thi(2), i, j, iloc, jloc, ld
      logical dotask
c
c     generate density matrix from orbitals in g_orbs. the first
c     nocc orbitals are doubly occupied.
c
      call ga_zero(g_counter)
      dotask = next_chunk(lo,hi)
      ld = ichunk
      do while (dotask)
        tlo(1) = lo(1)
        thi(1) = hi(1)
        tlo(2) = 1
        thi(2) = nocc
        call nga_get(g_orbs,tlo,thi,orbsi,ld)
        tlo(1) = lo(2)
        thi(1) = hi(2)
        call nga_get(g_orbs,tlo,thi,orbsj,ld)
        do i = lo(1), hi(1)
          iloc = i - lo(1) + 1
          do j = lo(2), hi(2)
            jloc = j - lo(2) + 1
            p = 0.0d00
            do k = 1, nocc
              p = p + orbsi(iloc,k)*orbsj(jloc,k)
            end do
            work(iloc,jloc) = 2.0d00*p
          end do
        end do
        call nga_put(g_work,lo,hi,work,ld)
        dotask = next_chunk(lo,hi)
      end do
      return
      end
c
      subroutine oneel(schwmax, eone)
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer lo(2), hi(2), i, j, iloc, jloc, ld
      dimension dens(nnbfn), fock(nnbfn), schwarz(nnbfn)
      dimension work(ichunk,ichunk),tfock(ichunk,ichunk)
      logical dotask
c
c     fill in the one-electron part of the fock matrix and
c     compute the one-electron energy contribution
c
      me = ga_nodeid()
      nproc = ga_nnodes()
c
      call ga_zero(g_counter)
      dotask = next_chunk(lo,hi)
      ld = ichunk
      do while (dotask)
        call nga_get(g_schwarz,lo,hi,work,ld)
        do j = lo(2), hi(2)
          jloc = j - lo(2) + 1
          do i = lo(1), hi(1)
            iloc = i - lo(1) + 1
            tfock(iloc,jloc) = 0.0d00
            if (work(iloc,jloc)*schwmax.gt.tol2e)
     +        tfock(iloc,jloc) = h(i,j)
          end do
        end do
        call nga_put(g_fock,lo,hi,tfock,ld)
        dotask = next_chunk(lo,hi)
      end do
      eone = 0.5d00*contract_matrices(g_fock,g_dens)
c
      end
      integer function nxtask(nproc)
      parameter (ichunk = 10)
      save icount, nleft
      data nleft, icount /0, 0/
c
c     wrapper round nxtval() to increase granularity
c     and thus reduce no. of requests to shared counter
c
      if(nproc.gt.0) then
         if(nleft.eq.0) then
            icount = nxtval(nproc) * ichunk
            nleft = ichunk
         endif
         nxtask = icount
         icount = icount + 1
         nleft = nleft -1
      else
          nleft = 0
          nxtask = 0
          junk = nxtval(nproc)
      endif
c
c     following does dumb static load balancing
c
c$$$      if(nproc.gt.0) then
c$$$         if (nleft .eq. 0) then
c$$$            icount = ga_nodeid()
c$$$            nleft = 1
c$$$         endif
c$$$         nxtask = icount
c$$$         icount = icount + ga_nnodes()
c$$$      else
c$$$          nleft = 0
c$$$          nxtask = 0
c$$$      endif
      end
c
      logical function next_chunk(lo,hi)
#include "cscf.h"
      integer one
      parameter (one = 1)
      integer imax, lo(2), hi(2), ilo, jlo
      itask = nga_read_inc(g_counter,one,one)
      imax = nbfn/ichunk
      if (itask.lt.imax*imax) then
        if (nbfn - ichunk*imax.gt.0) imax = imax + 1
        ilo = mod(itask,imax)
        jlo = (itask-ilo)/imax
        lo(1) = ilo*ichunk + 1
        lo(2) = jlo*ichunk + 1
        hi(1) = min((ilo+1)*ichunk,nbfn)
        hi(2) = min((jlo+1)*ichunk,nbfn)
        next_chunk = .true.
      else
        next_chunk = .false.
      endif
      return
      end
c
      logical function next_4chunk(lo,hi,ilo,jlo,klo,llo)
#include "cscf.h"
      integer one
      parameter (one = 1)
      integer imax, lo(4), hi(4), ilo, jlo, klo, llo, itmp
      itask = nga_read_inc(g_counter,one,one)
      imax = nbfn/ichunk
      if (nbfn - ichunk*imax.gt.0) imax = imax + 1
      if (itask.lt.imax**4) then
        ilo = mod(itask,imax)
        itmp = (itask - ilo)/imax
        jlo = mod(itmp,imax)
        itmp = (itmp - jlo)/imax
        klo = mod(itmp,imax)
        llo = (itmp - klo)/imax
        lo(1) = ilo*ichunk + 1
        lo(2) = jlo*ichunk + 1
        lo(3) = klo*ichunk + 1
        lo(4) = llo*ichunk + 1
        hi(1) = min((ilo+1)*ichunk,nbfn)
        hi(2) = min((jlo+1)*ichunk,nbfn)
        hi(3) = min((klo+1)*ichunk,nbfn)
        hi(4) = min((llo+1)*ichunk,nbfn)
        next_4chunk = .true.
      else
        next_4chunk = .false.
      endif
      return
      end
c
      subroutine clean_chunk(chunk)
#include "cscf.h"
      double precision chunk(ichunk,ichunk)
      integer i,j
      do j = 1, ichunk
        do i = 1, ichunk
          chunk(i,j) = 0.0d00
        end do
      end do
      return
      end
c
      subroutine twoel(schwmax, etwo)
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      double precision f_ij(ichunk,ichunk),d_kl(ichunk,ichunk)
      double precision f_ik(ichunk,ichunk),d_jl(ichunk,ichunk)
      double precision s_ij(ichunk,ichunk),s_kl(ichunk,ichunk)
      double precision schwmax, one
      integer nproc,lo(4),hi(4),ld,ich
      integer lo_ik(2),hi_ik(2),lo_jl(2),hi_jl(2)
      integer i,j,k,l,iloc,jloc,kloc,lloc,it,jt,kt,lt
      logical dotask, next_4chunk
c     
c     add in the two-electron contribution to the fock matrix
c     
      nproc = ga_nnodes()
      one = 1.0d00
c
      call ga_zero(g_counter)
      ld = maxnbfn
      ich = ichunk
      dotask = next_4chunk(lo,hi,it,jt,kt,lt)
      itask = 0
      do while (dotask)
        lo_ik(1) = lo(1)
        lo_ik(2) = lo(3)
        hi_ik(1) = hi(1)
        hi_ik(2) = hi(3)
        lo_jl(1) = lo(2)
        lo_jl(2) = lo(4)
        hi_jl(1) = hi(2)
        hi_jl(2) = hi(4)
        call nga_get(g_schwarz,lo,hi,s_ij,ich)
        call nga_get(g_schwarz,lo(3),hi(3),s_kl,ich)
        call nga_get(g_dens,lo(3),hi(3),d_kl,ich)
        call nga_get(g_dens,lo_jl,hi_jl,d_jl,ich)
        itask = itask + 1
        call clean_chunk(f_ij)
        call clean_chunk(f_ik)
        do i = lo(1), hi(1)
          iloc = i-lo(1) + 1
          do j = lo(2), hi(2)
            jloc = j-lo(2) + 1
            if (s_ij(iloc,jloc)*schwmax .lt. tol2e) then
              icut1 = icut1 + (hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)
            else
              do k = lo(3), hi(3)
                kloc = k-lo(3) + 1
                do l = lo(4), hi(4)
                  lloc = l-lo(4) + 1
                  if (s_ij(iloc,jloc)*s_kl(kloc,lloc).lt.tol2e) then
                    icut2 = icut2 + 1
                  else
                    call g(gg, i, j, k, l)
                    f_ij(iloc,jloc) = f_ij(iloc,jloc)
     +                              + gg*d_kl(kloc,lloc)
                    f_ik(iloc,kloc) = f_ik(iloc,kloc)
     +                              - 0.5d00*gg*d_jl(jloc,lloc)
                    icut3 = icut3 + 1
                  endif
                end do
              end do
            endif
          end do
        end do
        call nga_acc(g_fock,lo,hi,f_ij,ich,one)
        call nga_acc(g_fock,lo_ik,hi_ik,f_ik,ich,one)
        dotask = next_4chunk(lo,hi,it,jt,kt,lt)
      end do
      etwo = 0.5d00*contract_matrices(g_fock,g_dens)
      return
      end
c
      subroutine damp(fac)
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c
c    create damped density matrix as a linear combination of
c    old density matrix and density matrix formed from new orbitals
c
      ofac = 1.0d0 - fac
      call ga_add(fac,g_dens,ofac,g_work,g_dens)
      return
      end
c
      subroutine prnout(iter, energy, deltad, tester)
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c
c     printout results of each iteration
c
      if (ga_nodeid().ne.0) return
      write(6,1) iter, energy, deltad, tester
      call flush(6)
1     format(' iter=',i3,', energy=',f13.8,', deltad=',d9.2,
     $     ', deltaf=',d9.2)
      return
      end
c
      double precision function dendif()
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      double precision xdiff
      dimension dens_c(ichunk,ichunk),work_c(ichunk,ichunk)
      integer lo(2), hi(2), i, j, ld
      logical dotask
c
c     compute largest change in density matrix elements
c
      denmax = 0.0d0
      call ga_zero(g_counter)
      dotask = next_chunk(lo,hi)
      ld = ichunk
      do while(dotask)
        call nga_get(g_dens,lo,hi,dens_c,ld)
        call nga_get(g_work,lo,hi,work_c,ld)
        do j = 1, hi(2)-lo(2)+1
          do i = 1, hi(1)-lo(1)+1
            xdiff = abs(dens_c(i,j)-work_c(i,j))
            if (xdiff.gt.denmax) denmax = xdiff
          end do
        end do
        dotask = next_chunk(lo,hi)
      end do
      call ga_dgop(1,denmax,1,'max')
      dendif = denmax
      return
      end
c
      double precision function testfock()
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      double precision xmax, xtmp
      dimension work(ichunk,ichunk)
      integer lo(2), hi(2), i, j, iloc, jloc, ld
      logical dotask
c
c     compute largest change in density matrix elements
c
      xmax = 0.0d0
      call ga_zero(g_counter)
      dotask = next_chunk(lo,hi)
      ld = ichunk
      do while(dotask)
        call nga_get(g_fock,lo,hi,work,ld)
        do j = lo(2), hi(2)
          jloc = j - lo(2) + 1
          do i = lo(1), hi(1)
            iloc = i - lo(1) + 1
            if (i.ne.j) then
              xtmp = abs(work(iloc,jloc))
              if (xtmp.gt.xmax) xmax = xtmp
            endif
          end do
        end do
        dotask = next_chunk(lo,hi)
      end do
      call ga_dgop(1,xmax,1,'max')
      testfock = xmax
      return
      end
c
      subroutine shiftfock(shift)
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      double precision shift
      dimension work(ichunk,ichunk)
      integer lo(2), hi(2), i, j, iloc, jloc, ld, icnt
      logical dotask
c
c     compute largest change in density matrix elements
c
      call ga_zero(g_counter)
      dotask = next_chunk(lo,hi)
      ld = ichunk
      do while(dotask)
        call nga_get(g_fock,lo,hi,work,ld)
        icnt = 0
        do j = lo(2), hi(2)
          jloc = j - lo(2) + 1
          do i = lo(1), hi(1)
            iloc = i - lo(1) + 1
            if (i.eq.j.and.i.gt.nocc) then
              work(iloc,jloc) = work(iloc,jloc) + shift
              icnt = icnt + 1
            endif
          end do
        end do
        if (icnt.gt.0) call nga_put(g_fock,lo,hi,work,ld)
        dotask = next_chunk(lo,hi)
      end do
      return
      end
c
      subroutine prnfin(energy)
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      dimension orbs(maxnbfn, maxnbfn)
      integer lo(2),hi(2),ld
c
c     printout final results
c
      if (ga_nodeid().ne.0) return
      write(6,1) energy
 1    format(//' final energy = ',f16.11//' eigenvalues')
      call output(eigv, 1, min(nbfn,nocc+5), 1, 1, nbfn, 1, 1)
      write(6,2)
2     format(//' eigenvectors ')
      lo(1) = 1
      lo(2) = 1
      hi(1) = nbfn
      hi(2) = nbfn
      ld = nbfn
      call nga_get(g_orbs,lo,hi,orbs,ld) 
      call output(orbs, 1, nbfn, 1, nocc, nbfn, nbfn, 1)
c
      return
      end
      subroutine g(value,i,j,k,l)
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c
c     compute the two electon integral (ij|kl) over normalized
c     primitive 1s gaussians
c
      f0val = 0.0d0
      rab2 = (x(i)-x(j))**2 + (y(i)-y(j))**2 + (z(i)-z(j))**2
      rcd2 = (x(k)-x(l))**2 + (y(k)-y(l))**2 + (z(k)-z(l))**2
      facij = expnt(i)*expnt(j)/(expnt(i)+expnt(j))
      fackl = expnt(k)*expnt(l)/(expnt(k)+expnt(l))
      exijkl = exprjh(- facij*rab2 - fackl*rcd2)
      denom = (expnt(i)+expnt(j))*(expnt(k)+expnt(l)) *
     $        sqrt(expnt(i)+expnt(j)+expnt(k)+expnt(l))
      fac = (expnt(i)+expnt(j))*(expnt(k)+expnt(l)) /
     $        (expnt(i)+expnt(j)+expnt(k)+expnt(l))
c
      xp = (x(i)*expnt(i) + x(j)*expnt(j))/(expnt(i)+expnt(j))
      yp = (y(i)*expnt(i) + y(j)*expnt(j))/(expnt(i)+expnt(j))
      zp = (z(i)*expnt(i) + z(j)*expnt(j))/(expnt(i)+expnt(j))
      xq = (x(k)*expnt(k) + x(l)*expnt(l))/(expnt(k)+expnt(l))
      yq = (y(k)*expnt(k) + y(l)*expnt(l))/(expnt(k)+expnt(l))
      zq = (z(k)*expnt(k) + z(l)*expnt(l))/(expnt(k)+expnt(l))
      rpq2 = (xp-xq)**2 + (yp-yq)**2 + (zp-zq)**2
c
      call f0(f0val, fac*rpq2)
      value = (2.0d0 * pi**2.5d0 / denom) * exijkl * f0val *
     $    rnorm(i)*rnorm(j)*rnorm(k)*rnorm(l)
      return
      end
c
      subroutine diagon(tester, iter)
c      subroutine diagon(fock, orbs, evals, work, tester, iter)
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      double precision r_zero, r_one, shift, tester
c
#if USE_TRANSFORM
c
c     use similarity transform to solve standard eigenvalue problem
c     (overlap matrix has been transformed out of the problem)
c
      r_one = 1.0d00
      r_zero = 0.0d00
      call ga_dgemm('n','n',nbfn,nbfn,nbfn,r_one,g_fock,g_orbs,
     +               r_zero,g_tfock)
      call ga_dgemm('t','n',nbfn,nbfn,nbfn,r_one,g_orbs,g_tfock,
     +               r_zero,g_fock)
      tester = testfock()
      shift = 0.0d00
      if (tester.gt.0.3d0) then
        shift = 0.3d0
      endif
      if (iter.ge.2.and.shift.ne.0.0d00) then
        call shiftfock(shift)
      endif
      call ga_copy(g_orbs,g_tfock)
      call ga_diag_std_seq(g_fock, g_work, eigv)
c
c     Back transform eigenvectors
c
      call ga_dgemm('n','n',nbfn,nbfn,nbfn,r_one,g_tfock,g_work,
     +               r_zero,g_orbs)
      if (iter.ge.2.and.shift.ne.0.0d00) then
         do 50 i = nocc+1, nbfn
            eigv(i) = eigv(i) - shift
 50      continue
      endif
#else
c
c     Keep remaking overlap matrix since ga_diag_seq does not
c     guarantee that g_ident is preserved.
c
         call makoverlap
         call ga_diag_seq(g_fock, g_ident, g_orbs, eigv)
         tester = 0.0d00
#endif
      return
      end
c
      subroutine makeob
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      double precision work(ichunk,ichunk),orbs(ichunk,ichunk)
      double precision eval(maxnbfn)
      integer lo(2),hi(2),ld,me,i,j,iloc,jloc
      logical dotask
c
c     generate set of orthonormal vectors by creating a random
c     symmetric matrix and solving associated generalized eigenvalue
c     problem using the correct overlap matrix.
c
      me = ga_nodeid()
      call ga_zero(g_counter)
      dotask = next_chunk(lo,hi)
      ld = ichunk
      do while (dotask)
        do j = lo(2), hi(2)
          jloc = j - lo(2) + 1
          do i = lo(1), hi(1)
            iloc = i - lo(1) + 1
            work(iloc,jloc) = s(i,j)
            orbs(iloc,jloc) = drand48(0)
          end do
        end do
        call nga_put(g_ident,lo,hi,work,ld)
        call nga_put(g_fock,lo,hi,orbs,ld)
        dotask = next_chunk(lo,hi)
      end do
      call ga_symmetrize(g_fock)
      call ga_diag_seq(g_fock, g_ident, g_orbs, eval)
c
      return
      end
c
      subroutine denges
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
c
c     Form guess density from superposition of atomic densities in the AO
c     basis set ... instead of doing the atomic SCF hardwire for this
c     small basis set for the Be atom.
c
      integer one, itask, lo(2), hi(2), ld
      dimension atdens(15,15)
      data atdens/
     $     0.000002,0.000027,0.000129,0.000428,0.000950,0.001180,
     $     0.000457,-0.000270,-0.000271,0.000004,0.000004,0.000004,
     $     0.000004,0.000004,0.000004,0.000027,0.000102,0.000987,
     $     0.003269,0.007254,0.009007,0.003492,-0.002099,-0.002108,
     $     0.000035,0.000035,0.000035,0.000035,0.000035,0.000035,
     $     0.000129,0.000987,0.002381,0.015766,0.034988,0.043433,
     $     0.016835,-0.010038,-0.010082,0.000166,0.000166,0.000166,
     $     0.000166,0.000166,0.000166,0.000428,0.003269,0.015766,
     $     0.026100,0.115858,0.144064,0.055967,-0.035878,-0.035990,
     $     0.000584,0.000584,0.000584,0.000584,0.000584,0.000584,
     $     0.000950,0.007254,0.034988,0.115858,0.128586,0.320120,
     $     0.124539,-0.083334,-0.083536,0.001346,0.001346,0.001346,
     $     0.001346,0.001346,0.001346,0.001180,0.009007,0.043433,
     $     0.144064,0.320120,0.201952,0.159935,-0.162762,-0.162267,
     $     0.002471,0.002471,0.002471,0.002471,0.002471,0.002471,
     $     0.000457,0.003492,0.016835,0.055967,0.124539,0.159935,
     $     0.032378,-0.093780,-0.093202,0.001372,0.001372,0.001372,
     $     0.001372,0.001372,0.001372,-0.000270,-0.002099,-0.010038,
     $     -0.035878,-0.083334,-0.162762,-0.093780,0.334488,0.660918,
     $     -0.009090,-0.009090,-0.009090,-0.009090,-0.009090,-0.009090,
     $     -0.000271,-0.002108,-0.010082,-0.035990,-0.083536,-0.162267,
     $     -0.093202,0.660918,0.326482,-0.008982,-0.008982,-0.008981,
     $     -0.008981,-0.008981,-0.008982,0.000004,0.000035,0.000166,
     $     0.000584,0.001346,0.002471,0.001372,-0.009090,-0.008982,
     $     0.000062,0.000124,0.000124,0.000124,0.000124,0.000124,
     $     0.000004,0.000035,0.000166,0.000584,0.001346,0.002471,
     $     0.001372,-0.009090,-0.008982,0.000124,0.000062,0.000124,
     $     0.000124,0.000124,0.000124,0.000004,0.000035,0.000166,
     $     0.000584,0.001346,0.002471,0.001372,-0.009090,-0.008981,
     $     0.000124,0.000124,0.000062,0.000124,0.000124,0.000124,
     $     0.000004,0.000035,0.000166,0.000584,0.001346,0.002471,
     $     0.001372,-0.009090,-0.008981,0.000124,0.000124,0.000124,
     $     0.000062,0.000124,0.000124,0.000004,0.000035,0.000166,
     $     0.000584,0.001346,0.002471,0.001372,-0.009090,-0.008981,
     $     0.000124,0.000124,0.000124,0.000124,0.000062,0.000124,
     $     0.000004,0.000035,0.000166,0.000584,0.001346,0.002471,
     $     0.001372,-0.009090,-0.008982,0.000124,0.000124,0.000124,
     $     0.000124,0.000124,0.000062/
c
c   Create initial guess for density matrix in global array
c
      call ga_zero(g_dens)
      call ga_zero(g_counter)
      one = 1
      ld = 15
c
c   Correct for a factor of two along the diagonal
c
      do i = 1, ld
        atdens(i,i) = 2.0d00*atdens(i,i)
      end do
      itask = nga_read_inc(g_counter,one,one)
      do while(itask.lt.natom)
        ioff = itask*15
        lo(1) = ioff+1
        lo(2) = ioff+1
        hi(1) = ioff+15
        hi(2) = ioff+15
        call nga_put(g_dens,lo,hi,atdens,ld)
        itask = nga_read_inc(g_counter,one,one)
      end do
      call ga_sync
      return
      end
c
      subroutine setarrays
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer one, two, dims(2)
      logical status
      one = 1
      two = 2
      g_counter = ga_create_handle()
      call ga_set_data(g_counter,one,one,MT_INT)
      status = ga_allocate(g_counter)
      call ga_zero(g_counter)

      dims(1) = nbfn
      dims(2) = nbfn
      g_dens = ga_create_handle()
      call ga_set_data(g_dens, two, dims, MT_DBL)
      status = ga_allocate(g_dens)
      call ga_zero(g_dens)

      g_schwarz = ga_create_handle()
      call ga_set_data(g_schwarz, two, dims, MT_DBL)
      status = ga_allocate(g_schwarz)
      call ga_zero(g_schwarz)

      g_fock = ga_create_handle()
      call ga_set_data(g_fock, two, dims, MT_DBL)
      status = ga_allocate(g_fock)
      call ga_zero(g_fock)

      g_tfock = ga_create_handle()
      call ga_set_data(g_tfock, two, dims, MT_DBL)
      status = ga_allocate(g_tfock)
      call ga_zero(g_tfock)

      g_work = ga_create_handle()
      call ga_set_data(g_work, two, dims, MT_DBL)
      status = ga_allocate(g_work)
      call ga_zero(g_work)
#ifdef CHECKPOINT_SCF
      arraylist(1) = g_counter
      arraylist(2) = g_dens
      arraylist(3) = g_schwarz
      arraylist(4) = g_fock
      arraylist(5) = g_work
      write (6,*) arraylist(1),arraylist(2),arraylist(3)
      call ga_checkpoint_arrays(arraylist,numarrays)
#endif

      g_ident = ga_create_handle()
      call ga_set_data(g_ident, two, dims, MT_DBL)
      status = ga_allocate(g_ident)
      call ga_zero(g_ident)

      g_orbs = ga_create_handle()
      call ga_set_data(g_orbs, two, dims, MT_DBL)
      status = ga_allocate(g_orbs)
      call ga_zero(g_orbs)

      return
      end
      subroutine closearrays
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      logical status
c
      status = ga_destroy(g_counter)
      status = ga_destroy(g_dens)
      status = ga_destroy(g_schwarz)
      status = ga_destroy(g_fock)
      status = ga_destroy(g_tfock)
      status = ga_destroy(g_work)
      status = ga_destroy(g_ident)
      status = ga_destroy(g_orbs)
c
      return
      end
c
      subroutine makoverlap
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer me, lo(2), hi(2), ptr, ld(2)
      integer ld1, ld2
      me = ga_nodeid()
      call nga_distribution(g_ident, me, lo, hi)
      call nga_access(g_ident, lo, hi, ptr, ld)
      ld1 = hi(1) - lo(1) + 1
      ld2 = hi(2) - lo(2) + 1
      call setoverlap(dbl_mb(ptr),lo,hi,ld1,ld2)
      call nga_release(g_ident)
      return
      end
c
      subroutine setoverlap(a,lo,hi,ld1,ld2)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer lo(2), hi(2)
      integer ld1, ld2, ii, jj
      double precision a(ld1,ld2)
      do i = 1, ld1
        ii = i + lo(1) - 1
        do j = 1, ld2
          jj = j + lo(2) - 1
#if USE_TRANSFORM
          if (ii.eq.jj) then
            a(i,j) = 1.0d00
          else
            a(i,j) = 0.0d00
          endif
#else
          a(i,j) = s(ii,jj)
#endif
        end do
      end do
      return
      end
c
      subroutine print_ga_block(g_a)
      implicit double precision(a-h,o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer lo(2), hi(2), ptr, ld1, ld2
c
      me = ga_nodeid()
      call nga_distribution(g_a, me, lo, hi)
      ld1 = hi(1) - lo(1) + 1
      ld2 = hi(2) - lo(2) + 1
      call nga_access(g_a, lo, hi, ptr, ld)
      call dump_chunk(dbl_mb(ptr),ld1,ld2)
      call nga_release(g_a)
c
      return
      end
c
      subroutine print_ga_block_ij(g_a,tlo)
      implicit double precision(a-h,o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer lo(2), hi(2), ptr, ld1, ld2
c
      me = ga_nodeid()
      call nga_distribution(g_a, me, lo, hi)
      ld1 = hi(1) - lo(1) + 1
      ld2 = hi(2) - lo(2) + 1
      call nga_access(g_a, tlo, hi, ptr, ld)
      call dump_chunk(dbl_mb(ptr),ld1,ld2)
      call nga_release(g_a)
c
      return
      end
c
      subroutine dump_chunk(a,ld1,ld2)
      implicit double precision (a-h, o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer ld1, ld2
      double precision a(ld1, ld2)
      do i = 1, min(10,ld1)
        write(6,100) (a(i,j), j = 1, min(10,ld2))
      end do
      write(6,*)
  100 format(10f10.4)
      return
      end
c
      double precision function contract_matrices(g_a,g_b)
      implicit double precision(a-h,o-z)
#include "cscf.h"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer lo(2), hi(2), ptr_a, ptr_b, ld, ld1, ld2
      double precision a(ichunk,ichunk),b(ichunk,ichunk)
      double precision value
      logical dotask
c
c   evalute sum_ij a_ij*b_ij
c
      value = 0.0d00
      call ga_zero(g_counter)
      dotask = next_chunk(lo,hi)
      ld = ichunk
      do while (dotask)
        call nga_get(g_a,lo,hi,a,ld)
        call nga_get(g_b,lo,hi,b,ld)
        do j = 1, hi(2)-lo(2)+1
          do i = 1, hi(1)-lo(1)+1
            value = value + a(i,j)*b(i,j)
          end do
        end do
        dotask = next_chunk(lo,hi)
      end do
      call ga_dgop(3,value,1,'+')
      contract_matrices=value
c
      return
      end
