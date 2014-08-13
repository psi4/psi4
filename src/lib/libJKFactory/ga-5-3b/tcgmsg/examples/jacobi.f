      program main
C$Id: jacobi.f,v 1.2 1995-02-02 23:24:15 d3g681 Exp $
      implicit double precision (a-h, o-z)
      dimension q(1)
      include 'msgtypesf.h'
      parameter (niter = 30, maxp = 2047)
      dimension jlo(0:maxp), nj(0:maxp)
      data iqoff, iqaddr/0, 0/
c
c     Simple iterative (Jacobi) linear equation solver (Ax = b)
c
c     Initialize the message passing environment
c
      call pbeginf
c     call evon
c
c     Read in the dimension  of the matrix on process 0 and
c     broadcast its value to everyone else
c
      if (nodeid() .eq. 0) then
         write(6,1)
 1       format(' Input dimension of matrix: '$)
         read(5,*) n
      endif
      call brdcst(1+MSGINT, n, mitob(1), 0)
c
c     makjlo assigns a range of columns to each process.
c     ncolp is the number that this process is to do.
c
      call makjlo(n, jlo, nj)
      ncolp = nj(nodeid())
c
c     allocate memory
c
      need = n*ncolp + 2*n + 2*ncolp
      call getmem(need, q, iqaddr, iqoff)
      if (iqaddr.eq.0) call parerr(999)
      ia = iqoff + 1
      ib = ia + n*ncolp
      ix = ib + ncolp
      is = ix + ncolp
      iw = is + n
c
c     Make matrix (a), rhs vector (b) and initial guess (x)
c
      call makabx(n, jlo(nodeid()), ncolp, q(ia), q(ib), q(ix))
c
c     Do niter iterations of the Jacobi algorithm.
c     Synchronize first for accurate timings.
c     
      call synch(13)
      rjunk = timer()
      call jacobi(n, jlo, nj, q(ia), q(ib), q(ix), q(is), q(iw), niter)
      used = timer()
c
c     Print out results
c
      rmflop = dble(niter*2*n)*dble(n) / (used*1.0d6)
      if (nodeid() .eq. 0) write(6,2) n, used, nnodes(), rmflop
 2    format(' N=',i4,' used ',f6.2,' secs with ',i3,' processes',
     $     ', mflop=', f8.2)
c
      call pend
      call fexit
      end
      subroutine jacobi(n, jlo, nj, a, b, x, s, work, niter)
      implicit double precision (a-h, o-z)
      dimension a(n, *), b(*), x(*), s(n), work(n),
     $     jlo(0:*), nj(0:*)
c
c     Apply niter iterations of the Jacobi algorithm.
c
      me = nodeid()
      do 10 iter = 1, niter
c
c     Compute matrix vector product Ax ... this is the real work
c
c     Do the part that we have
c
         call mxv(a, n, x, nj(me), s)
c
c     Now we have to add up the result over all the processors
c     Call dgop for simple but inefficient version. Call mxvadd
c     for a much more efficient version
c
c         call dgop(2, s, n, '+')
c     
         call mxvadd(s, work, jlo, nj)
c
c     Compute our part of the update vector and compute
c     the residual error (the error requires a global sum)
c
         err = 0.0d0
         do 20 j = jlo(me), jlo(me) + nj(me) - 1
            jj = j - jlo(me) + 1
            x(jj) = x(jj) + (b(jj)-s(j))/a(j,jj)
            err = err + abs((b(jj)-s(j)))
 20      continue
c
c     Write out results every now and again
c
         if (mod(iter,10).eq.0) then
            call dgop(3, err, 1, '+')
            if (nodeid().eq.0) write(6,1) iter, err
 1          format(' Iteration',i3,', Error',d9.2)
         endif
 10   continue
c
      end
      subroutine makabx(n, jlo, nj, a, b, x)
      implicit double precision (a-h, o-z)
      dimension a(n, nj), b(nj), x(nj)
c
      jhi = jlo + nj - 1
      do 10 j = jlo, jhi
         jj = j - jlo + 1
cvd$  novector
         do 20 i = 1, n
            a(i,jj) = dble(i+j) / dble(abs(i-j)*n+n/50+1)
 20      continue
         b(jj) = dble(mod(j,3))
         x(jj) = b(jj) / a(j,jj)
 10   continue
c
      end
      subroutine makjlo(n, jlo, nj)
      dimension jlo(0:*), nj(0:*)
c
      ncolp = n / nnodes()
      next  = n - (ncolp*nnodes())
      jjlo = 1
      do 10 iproc = 0, nnodes()-1
         jlo(iproc) = jjlo
         jjlo = jjlo + ncolp
         if (iproc.lt.next) jjlo = jjlo + 1
         nj(iproc) = jjlo - jlo(iproc)
 10   continue
c
      end
      subroutine mxvadd(s, work, jlo, nj)
      implicit real*8 (a-h, o-z)
      include 'msgtypesf.h'
      dimension s(*), work(*), jlo(0:*), nj(0:*)
      logical synch
      parameter (synch=.true.)
c
c     We have in s(1:n) this process's contribution to the
c     matrix vector product A*x where we had nj(me) elements
c     of x starting at element jlo(me). Each process needs
c     to end up with the same elements of the full result
c     vector s.
c
c     Thus we need to send to each process (ip) the elements
c     s(jlo(ip)+k-1), k=1,nj(ip). And we need to receive from
c     each process our piece of s which we add onto our result
c     vector.
c
c     If communication is synchronous then we must explictly pair up
c     send/receive requests on this process with the matching
c     receive/send operations on other processes.
c
c     There is potential for much more asynch stuff but the damned
c     iPSC-i860 hangs (irreproducibly) if we send off too many
c     unresolved asynchronous messages (how many is too much?).
c
      me = nodeid()
      nproc = nnodes()
      nn = nproc + mod(nproc,2)
c
      if (synch) then
        do 10 iter = 1, nn-1
          call pairup(nn, me, iter, ip)
          if (ip.lt.nproc) then
            if (me. lt. ip) then
              call snd(3+MSGDBL, s(jlo(ip)), mdtob(nj(ip)), ip, 1)
              call rcv(3+MSGDBL, work, mdtob(nj(me)), lenmes, ip,
     &                 node, 1)
            else if (me.gt.ip) then
               call rcv(3+MSGDBL, work, mdtob(nj(me)), lenmes, ip,
     &                  node, 1)
               call snd(3+MSGDBL, s(jlo(ip)), mdtob(nj(ip)), ip, 1)
            endif
            call vvadd(nj(me), s(jlo(me)), work)
          endif
 10     continue
      else
        do 20 iter = 1, nn-1
          call pairup(nn, me, iter, ip)
          if (ip.lt.nproc) then
            call snd(3+MSGDBL, s(jlo(ip)), mdtob(nj(ip)), ip, 0)
            call rcv(3+MSGDBL, work, mdtob(nj(me)), lenmes, ip, node, 1)
            call vvadd(nj(me), s(jlo(me)), work)
          endif
 20     continue
        call waitcom(-1)
      endif
c
      end
      subroutine pairup(n, me, iter, ipair)
c
c     one of many ways of generating maximally overlapped pairs
c     (not all that good on a hypercube though!)
c
      if (iter.eq.1) then
        ipair = mod(n+1-me,n)
      else if (me.eq.0) then
        ipair = iter
      else if (me.eq.iter) then
        ipair = 0
      else
        if (ipair.eq.0) ipair = me
        ipair = ipair + 2
        if (ipair.ge.n) ipair = ipair + 1 - n
      endif
      end
      subroutine vvadd(n, a, b)
      implicit real*8 (a-h, o-z)
      dimension a(*), b(*)
c
      do 10 i = 1,n
        a(i) = a(i) + b(i)
 10   continue
c
      end
