c
c     include file defining common /cscf/
c
c     constant parameters are set in cscf, parameters are set
c     after reading input.
c
c     maxatom  = maximum no. of atoms     (constant parameter)
c     maxnbfn  = maximum no. of bas. fnct.(constant parameter)
c     maxnnbfn = maxnbfn*(maxnbfn+1)/2    (constant parameter)
c     natom    = no. of atoms             (parameter)
c     nbfn     = no. of basis functions   (parameter)
c     nnbfn    = nbfn*(nbfn+1)/2          (parameter)
c     nocc     = no. of occupied orbitals (parameter)
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
      parameter (maxatom = 286)        !cste original value 50
      parameter (maxnbfn =15*maxatom, mxiter = 30)
      parameter (maxnnbfn = maxnbfn*(maxnbfn+1)/2)
      parameter (pi = 3.141592653589793d0)
      parameter (tol= 0.5d-3)
      parameter (tol2e=1.0d-6)
c
      common /cscf/
     $     enrep, q(maxatom), ax(maxatom), ay(maxatom), az(maxatom),
     $     x(maxnbfn), y(maxnbfn), z(maxnbfn), expnt(maxnbfn),
     $     rnorm(maxnbfn),iky(maxnbfn), icut1, icut2, icut3, icut4,
     $     natom, nocc, nbfn, nnbfn
      double precision enrep, q, ax, ay, az, x, y, z, expnt, rnorm
      integer*8 iky, icut1, icut2, icut3, icut4, natom, nocc, nbfn,
     $     nnbfn
c
c    Global array parameters used in calculations:
c
c    ichunk:    chunk size for distributing workload
c
c    g_counter: global array used to assign next task
c    g_dens:    global array used to store density matrix
c    g_fock:    global array used to store fock matrix
c    g_tfock:   global array used to store transformed fock matrix
c    g_schwarz: global array used to store schwarz matrix
c    g_work:    global array used to store work matrix
c    g_ident:   global array used to store identity matrix
c    g_orbs:    global array used to store orbital vectors
c
      parameter (ichunk = 20)        !cste original value 10
      common /g_arrays/ eigv(maxnbfn),
     $     g_counter, g_dens, g_fock, g_tfock, g_schwarz, g_work,
     $     g_ident, g_orbs
      double precision eigv
      integer g_counter, g_dens, g_fock, g_tfock, g_schwarz, g_work,
     $     g_ident, g_orbs
