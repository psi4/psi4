      block data
C$Id: blkdat15.f,v 1.2 1995-02-02 23:23:47 d3g681 Exp $
      implicit double precision (a-h, o-z)
      include 'cscf.h'
c
c     initalize data in common ... clumsy but avoids code to read in data
c
c     one be atom with 15 orbitals
c
c     have 9s functions on each center and simulate p's by having
c     s function at +- 1 in each of x, y, z
c
      data ax /0.0d0/, ay /0.0d0/, az /0.0d0/, q /4.0d0/, enrep/0.0d0/
      data x /9*0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/
      data y /9*0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0, 0.0d0, 0.0d0/
      data z /9*0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.6d0, -1.6d0/
      data expnt /1741.0d0, 262.1d0, 60.33d0, 17.62d0, 5.933d0, 2.185d0,
     $     0.859, 0.1806d0, 0.05835d0, 6*0.3d0/
c
      end
