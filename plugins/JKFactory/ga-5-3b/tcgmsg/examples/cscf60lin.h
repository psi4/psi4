/*$Id: cscf60lin.h,v 1.2 1995-02-02 23:24:04 d3g681 Exp $*/
c
c     include file defining common /cscf/
c
c     natom    = no. of atoms             (constant parameter)
c     nbfn     = no. of basis functions   (constant parameter)
c     nnbfn    = nbfn*(nbfn+1)/2          (constant parameter)
c     nocc     = no. of occupied orbitals (constant parameter)
c     mxiter   = maximim no. of iterations(constant parameter)
c     tol      = convergence criterion    (constant parameter)
c     pi       = a familiar constant      (constant parameter)
c     tol2e    = 2-e integral screening   (constant parameter)
c
c     the remainder is initialized in block data or in the
c     routine ininrm (rnorm and iky)
c
c     enrep        = nuclear repulsion energy
c     q(1:natom)   = nuclear charge of atom
c     ax(1:natom)  = x co-ordinate of atom
c     ay(1:natom)  = y ...
c     az(1:natom)  = z ...
c     x(1:nbfn)    = x co-ordinate of basis function
c     y(1:nbfn)    = y ...
c     z(1:nbfn)    = z ...
c     expnt(1:nbfn)= exponent of gaussian
c     rnorm(1:nbfn)= normalization constant of gaussian
c     iky(1:nbfn)  = iky(i) = i*(i-1)/2 to speed up fock build
c     icut1        = no. of successful ij   2-e screening tests
c     icut2        = no. of successful ijkl 2-e screening tests
c     icut3        = no. of 2-e integrals computed
c
      parameter (natom = 4, nbfn = 60, nocc = 8, mxiter = 20)
      parameter (nnbfn = nbfn*(nbfn+1)/2, pi = 3.141592653589793d0)
      parameter (tol= 1.0d-4)
      parameter (tol2e=1.0d-6)
c
      common /cscf/
     $     enrep, q(natom), ax(natom), ay(natom), az(natom),
     $     x(nbfn), y(nbfn), z(nbfn), expnt(nbfn), rnorm(nbfn),
     $     iky(nbfn), icut1, icut2, icut3
c

