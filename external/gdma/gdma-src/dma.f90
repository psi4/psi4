MODULE DMA

!  Distributed Multipole Analysis
!
!  Copyright (C) 2005-08  Anthony J. Stone
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to
!  the Free Software Foundation, Inc., 51 Franklin Street,
!  Fifth Floor, Boston, MA 02110-1301, USA.

use iso_c_binding
USE atom_grids, ONLY : grid, ng, make_grid, Lebedev,                   &
    n_a, n_r, k_mu, start
IMPLICIT NONE

PRIVATE
PUBLIC dma_main, c, ex, cs, cp, iax, kmin, kmax, kstart, katom, ktype, &
    kng, kloc, maxcen, maxbfn, maxs, nat, name, nshell, title, &
    zan, Qfactor, echarge, bohr, rfact, punchfile

INTEGER, PARAMETER :: dp=kind(1d0)

!  MAXS:   Maximum number of sites for a normal DMA.
!  MAXCEN: Maximum number of atoms.
!  MAXBFN: Maximum number of basis functions (after converting spherical
!          basis functions to cartesian form).
!  MAXEXP: Maximum number of primitives.

INTEGER :: maxs, maxcen, maxbfn
! INTEGER, PARAMETER :: MAXS=100, MAXCEN=100, MAXBFN=2000
! INTEGER, PARAMETER :: MAXEXP=5000, MAXSHL=1000

! bohr radius x 1e10
REAL(dp), PARAMETER :: bohr=0.529177211d0
REAL(dp), PARAMETER :: rtpi=1.7724538509055160272D0
REAL(dp) :: Qfactor(0:20)=1d0
!  elementary charge x 1e20
REAL(dp), PARAMETER :: echarge=16.0217653d0

INTEGER :: nat, nshell
REAL(dp), ALLOCATABLE :: zan(:), c(:,:)

INTEGER, ALLOCATABLE :: kstart(:), katom(:), ktype(:), kng(:),&
    kloc(:), kmin(:), kmax(:)
REAL(dp), ALLOCATABLE :: ex(:), cs(:), cp(:)

! REAL(dp) :: ex(MAXEXP)=0d0, cs(MAXEXP)=0d0, cp(MAXEXP)=0d0!,            &
!    cd(MAXEXP)=0d0

!  IX(m) is the power of x in the mth basis function in the list s, px,
!  py, ...
INTEGER :: IX(56),IY(56),IZ(56)
DATA IX(1:20) /0,  1,0,0,  2,0,0,1,1,0,  3,0,0,2,2,1,0,1,0,1/
DATA IY(1:20) /0,  0,1,0,  0,2,0,1,0,1,  0,3,0,1,0,2,2,0,1,1/
DATA IZ(1:20) /0,  0,0,1,  0,0,2,0,1,1,  0,0,3,0,1,0,1,2,2,1/
DATA IX(21:35) /4,0,0,3,3,1,0,1,0,2,2,0,2,1,1/
DATA IY(21:35) /0,4,0,1,0,3,3,0,1,2,0,2,1,2,1/
DATA IZ(21:35) /0,0,4,0,1,0,1,3,3,0,2,2,1,1,2/
DATA IX(36:56) /5,0,0,4,4,1,0,1,0,3,3,2,0,2,0,3,1,1,2,2,1/
DATA IY(36:56) /0,5,0,1,0,4,4,0,1,2,0,3,3,0,2,1,3,1,2,1,2/
DATA IZ(36:56) /0,0,5,0,1,0,1,4,4,0,2,0,2,3,3,1,1,3,1,2,2/

INTEGER, ALLOCATABLE :: limit(:)
REAL(dp), ALLOCATABLE :: xs(:,:), radius(:), q(:,:)
REAL(dp) :: rt(0:20), binom(0:20,0:20), rtbinom(0:20,0:20),         &
    d(56,56), qt(0:121)

LOGICAL:: slice, linear, planar, general
INTEGER :: ns, lmax, perp, mindc, maxdc
REAL(dp) :: tol, spread

!  Charge density at grid points
REAL(dp), ALLOCATABLE :: rho(:)

!  Crossover value for exponents
REAL(dp) :: bigexp
REAL(dp),PARAMETER :: bigexp_default=4.0

!  Neglect electron density terms involving exp(-e) with e>etol
REAL(dp) :: etol=36

CHARACTER(LEN=8), ALLOCATABLE :: name(:)
CHARACTER(LEN=80) :: title(2)
LOGICAL :: nuclei=.true.

INTEGER, ALLOCATABLE :: iax(:)

REAL(dp) :: rfact=1d0

CHARACTER(LEN=40) :: punchfile="gdma.punch", tempfile=""
CHARACTER(LEN=8) :: runit

!  Addresses of Gauss-Hermite points and weights.
INTEGER :: mink(20), maxk(20)
DATA mink /1,2,4, 7,11,16,22,29,37,46,56,67,79, 92,106,121,137,154,172,191/
DATA maxk /1,3,6,10,15,21,28,36,45,55,66,78,91,105,120,136,153,171,190,210/

!  Gauss-Hermite points and weights.
REAL(dp) :: H(210), W(210)

DATA H(1) /                                                       &
    0.00000000000000D+00/
DATA W(1) /                                                       &
    1.77245385090552D+00/

DATA H(2:3) /                                          &
    -7.07106781186548D-01, 7.07106781186547D-01/
DATA W(2:3) /                                          &
    8.86226925452758D-01, 8.86226925452758D-01/

DATA H(4:6) /                                          &
    -1.22474487139159D+00, 0.00000000000000D+00, 1.22474487139159D+00/
DATA W(4:6) /                                          &
    2.95408975150920D-01, 1.18163590060368D+00, 2.95408975150920D-01/

DATA H(7:10) /                                          &
    -1.65068012388578D+00,-5.24647623275290D-01, 5.24647623275290D-01,&
    1.65068012388578D+00/
DATA W(7:10) /                                          &
    8.13128354472457D-02, 8.04914090005513D-01, 8.04914090005513D-01,&
    8.13128354472455D-02/

DATA H(11:15) /                                          &
    -2.02018287045608D+00,-9.58572464613819D-01, 0.00000000000000D+00,&
    9.58572464613817D-01, 2.02018287045608D+00/
DATA W(11:15) /                                          &
    1.99532420590462D-02, 3.93619323152241D-01, 9.45308720482942D-01,&
    3.93619323152243D-01, 1.99532420590461D-02/

DATA H(16:21) /                                          &
    -2.35060497367449D+00,-1.33584907401370D+00,-4.36077411927616D-01,&
    4.36077411927616D-01, 1.33584907401369D+00, 2.35060497367449D+00/
DATA W(16:21) /                                          &
    4.53000990550896D-03, 1.57067320322857D-01, 7.24629595224393D-01,&
    7.24629595224393D-01, 1.57067320322858D-01, 4.53000990550891D-03/

DATA H(22:28) /                                          &
    -2.65196135683523D+00,-1.67355162876747D+00,-8.16287882858965D-01,&
    0.00000000000000D+00, 8.16287882858964D-01, 1.67355162876747D+00,&
    2.65196135683523D+00/
DATA W(22:28) /                                          &
    9.71781245099546D-04, 5.45155828191274D-02, 4.25607252610128D-01,&
    8.10264617556808D-01, 4.25607252610129D-01, 5.45155828191276D-02,&
    9.71781245099541D-04/

DATA H(29:36) /                                          &
    -2.93063742025724D+00,-1.98165675669584D+00,-1.15719371244678D+00,&
    -3.81186990207322D-01, 3.81186990207321D-01, 1.15719371244678D+00,&
    1.98165675669584D+00, 2.93063742025724D+00/
DATA W(29:36) /                                          &
    1.99604072211375D-04, 1.70779830074137D-02, 2.07802325814892D-01,&
    6.61147012558242D-01, 6.61147012558242D-01, 2.07802325814893D-01,&
    1.70779830074137D-02, 1.99604072211373D-04/

DATA H(37:45) /                                          &
    -3.19099320178152D+00,-2.26658058453184D+00,-1.46855328921666D+00,&
    -7.23551018752838D-01, 0.00000000000000D+00, 7.23551018752836D-01,&
    1.46855328921666D+00, 2.26658058453184D+00, 3.19099320178152D+00/
DATA W(37:45) /                                          &
    3.96069772632663D-05, 4.94362427553707D-03, 8.84745273943779D-02,&
    4.32651559002556D-01, 7.20235215606051D-01, 4.32651559002557D-01,&
    8.84745273943776D-02, 4.94362427553706D-03, 3.96069772632660D-05/

DATA H(46:55) /                                          &
    -3.43615911883773D+00,-2.53273167423278D+00,-1.75668364929988D+00,&
    -1.03661082978951D+00,-3.42901327223704D-01, 3.42901327223704D-01,&
    1.03661082978951D+00, 1.75668364929988D+00, 2.53273167423278D+00,&
    3.43615911883773D+00/
DATA W(46:55) /                                          &
    7.64043285523306D-06, 1.34364574678130D-03, 3.38743944554813D-02,&
    2.40138611082315D-01, 6.10862633735326D-01, 6.10862633735326D-01,&
    2.40138611082317D-01, 3.38743944554818D-02, 1.34364574678127D-03,&
    7.64043285523290D-06/

DATA H(56:66) /                                          &
    -3.66847084655957D+00,-2.78329009978165D+00,-2.02594801582575D+00,&
    -1.32655708449493D+00,-6.56809566882100D-01, 0.00000000000000D+00,&
    6.56809566882098D-01, 1.32655708449493D+00, 2.02594801582575D+00,&
    2.78329009978165D+00, 3.66847084655958D+00/
DATA W(56:66) /                                          &
    1.43956039371434D-06, 3.46819466323358D-04, 1.19113954449119D-02,&
    1.17227875167709D-01, 4.29359752356125D-01, 6.54759286914592D-01,&
    4.29359752356126D-01, 1.17227875167710D-01, 1.19113954449119D-02,&
    3.46819466323356D-04, 1.43956039371433D-06/

DATA H(67:78) /                                          &
    -3.88972489786977D+00,-3.02063702512088D+00,-2.27950708050106D+00,&
    -1.59768263515260D+00,-9.47788391240164D-01,-3.14240376254359D-01,&
    3.14240376254358D-01, 9.47788391240162D-01, 1.59768263515260D+00,&
    2.27950708050105D+00, 3.02063702512088D+00, 3.88972489786978D+00/
DATA W(67:78) /                                          &
    2.65855168435650D-07, 8.57368704358828D-05, 3.90539058462913D-03,&
    5.16079856158845D-02, 2.60492310264161D-01, 5.70135236262480D-01,&
    5.70135236262480D-01, 2.60492310264162D-01, 5.16079856158849D-02,&
    3.90539058462919D-03, 8.57368704358820D-05, 2.65855168435644D-07/

DATA H(79:91) /                                          &
    -4.10133759617864D+00,-3.24660897837240D+00,-2.51973568567823D+00,&
    -1.85310765160151D+00,-1.22005503659075D+00,-6.05763879171060D-01,&
    0.00000000000000D+00, 6.05763879171059D-01, 1.22005503659074D+00,&
    1.85310765160151D+00, 2.51973568567823D+00, 3.24660897837240D+00,&
    4.10133759617863D+00/
DATA W(79:91) /                                          &
    4.82573185007318D-08, 2.04303604027083D-05, 1.20745999271943D-03,&
    2.08627752961704D-02, 1.40323320687024D-01, 4.21616296898543D-01,&
    6.04393187921162D-01, 4.21616296898544D-01, 1.40323320687025D-01,&
    2.08627752961704D-02, 1.20745999271942D-03, 2.04303604027081D-05,&
    4.82573185007342D-08/

DATA H(92:105) /                                          &
    -4.30444857047362D+00,-3.46265693360226D+00,-2.74847072498540D+00,&
    -2.09518325850771D+00,-1.47668273114114D+00,-8.78713787329399D-01,&
    -2.91745510672562D-01, 2.91745510672561D-01, 8.78713787329397D-01,&
    1.47668273114114D+00, 2.09518325850771D+00, 2.74847072498539D+00,&
    3.46265693360226D+00, 4.30444857047362D+00/
DATA W(92:105) /                                          &
    8.62859116812580D-09, 4.71648435501913D-06, 3.55092613551935D-04,&
    7.85005472645823D-03, 6.85055342234660D-02, 2.73105609064247D-01,&
    5.36405909712091D-01, 5.36405909712091D-01, 2.73105609064248D-01,&
    6.85055342234660D-02, 7.85005472645827D-03, 3.55092613551940D-04,&
    4.71648435501930D-06, 8.62859116812587D-09/

DATA H(106:120) /                                          &
    -4.49999070730939D+00,-3.66995037340445D+00,-2.96716692790559D+00,&
    -2.32573248617385D+00,-1.71999257518649D+00,-1.13611558521092D+00,&
    -5.65069583255576D-01, 0.00000000000000D+00, 5.65069583255574D-01,&
    1.13611558521092D+00, 1.71999257518648D+00, 2.32573248617385D+00,&
    2.96716692790560D+00, 3.66995037340444D+00, 4.49999070730938D+00/
DATA W(106:120) /                                          &
    1.52247580425354D-09, 1.05911554771112D-06, 1.00004441232505D-04,&
    2.77806884291283D-03, 3.07800338725465D-02, 1.58488915795936D-01,&
    4.12028687498899D-01, 5.64100308726418D-01, 4.12028687498900D-01,&
    1.58488915795938D-01, 3.07800338725467D-02, 2.77806884291286D-03,&
    1.00004441232503D-04, 1.05911554771115D-06, 1.52247580425362D-09/

DATA H(121:136) /                                          &
    -4.68873893930580D+00,-3.86944790486012D+00,-3.17699916197995D+00,&
    -2.54620215784747D+00,-1.95178799091625D+00,-1.38025853919888D+00,&
    -8.22951449144656D-01,-2.73481046138152D-01, 2.73481046138152D-01,&
    8.22951449144654D-01, 1.38025853919888D+00, 1.95178799091625D+00,&
    2.54620215784747D+00, 3.17699916197995D+00, 3.86944790486011D+00,&
    4.68873893930581D+00/
DATA W(121:136) /                                          &
    2.65480747401156D-10, 2.32098084486533D-07, 2.71186009253801D-05,&
    9.32284008624214D-04, 1.28803115355103D-02, 8.38100413989866D-02,&
    2.80647458528534D-01, 5.07929479016614D-01, 5.07929479016614D-01,&
    2.80647458528535D-01, 8.38100413989867D-02, 1.28803115355103D-02,&
    9.32284008624218D-04, 2.71186009253805D-05, 2.32098084486537D-07,&
    2.65480747401144D-10/

DATA H(137:153) /                                          &
    -4.87134519367440D+00,-4.06194667587546D+00,-3.37893209114149D+00,&
    -2.75776291570388D+00,-2.17350282666661D+00,-1.61292431422123D+00,&
    -1.06764872574345D+00,-5.31633001342655D-01, 0.00000000000000D+00,&
    5.31633001342653D-01, 1.06764872574345D+00, 1.61292431422123D+00,&
    2.17350282666661D+00, 2.75776291570388D+00, 3.37893209114148D+00,&
    4.06194667587546D+00, 4.87134519367439D+00/
DATA W(137:153) /                                          &
    4.58057893079898D-11, 4.97707898163125D-08, 7.11228914002172D-06,&
    2.98643286697763D-04, 5.06734995762769D-03, 4.09200341497568D-02,&
    1.72648297670098D-01, 4.01826469470412D-01, 5.30917937624864D-01,&
    4.01826469470413D-01, 1.72648297670098D-01, 4.09200341497570D-02,&
    5.06734995762773D-03, 2.98643286697762D-04, 7.11228914002204D-06,&
    4.97707898163130D-08, 4.58057893079903D-11/

DATA H(154:171) /                                          &
    -5.04836400887446D+00,-4.24811787356812D+00,-3.57376906848625D+00,&
    -2.96137750553160D+00,-2.38629908916668D+00,-1.83553160426162D+00,&
    -1.30092085838962D+00,-7.76682919267412D-01,-2.58267750519097D-01,&
    2.58267750519096D-01, 7.76682919267409D-01, 1.30092085838961D+00,&
    1.83553160426162D+00, 2.38629908916668D+00, 2.96137750553159D+00,&
    3.57376906848626D+00, 4.24811787356811D+00, 5.04836400887446D+00/
DATA W(154:171) /                                          &
    7.82819977211642D-12, 1.04672057957931D-08, 1.81065448109359D-06,&
    9.18112686793001D-05, 1.88852263026852D-03, 1.86400423875452D-02,&
    9.73017476413160D-02, 2.84807285669980D-01, 4.83495694725456D-01,&
    4.83495694725456D-01, 2.84807285669981D-01, 9.73017476413168D-02,&
    1.86400423875452D-02, 1.88852263026848D-03, 9.18112686793015D-05,&
    1.81065448109355D-06, 1.04672057957932D-08, 7.82819977211684D-12/

DATA H(172:190) /                                          &
    -5.22027169053746D+00,-4.42853280660377D+00,-3.76218735196401D+00,&
    -3.15784881834759D+00,-2.59113378979453D+00,-2.04923170985061D+00,&
    -1.52417061939353D+00,-1.01036838713431D+00,-5.03520163423888D-01,&
    0.00000000000000D+00, 5.03520163423886D-01, 1.01036838713431D+00,&
    1.52417061939353D+00, 2.04923170985061D+00, 2.59113378979453D+00,&
    3.15784881834760D+00, 3.76218735196401D+00, 4.42853280660376D+00,&
    5.22027169053747D+00/
DATA W(172:190) /                                          &
    1.32629709449876D-12, 2.16305100986376D-09, 4.48824314722341D-07,&
    2.72091977631637D-05, 6.70877521407221D-04, 7.98886677772321D-03,&
    5.08103869090527D-02, 1.83632701306997D-01, 3.91608988613031D-01,&
    5.02974888276187D-01, 3.91608988613031D-01, 1.83632701306999D-01,&
    5.08103869090531D-02, 7.98886677772327D-03, 6.70877521407217D-04,&
    2.72091977631626D-05, 4.48824314722354D-07, 2.16305100986386D-09,&
    1.32629709449867D-12/

DATA H(191:210) /                                          &
    -5.38748089001122D+00,-4.60368244955073D+00,-3.94476404011561D+00,&
    -3.34785456738321D+00,-2.78880605842812D+00,-2.25497400208927D+00,&
    -1.73853771211658D+00,-1.23407621539532D+00,-7.37473728545394D-01,&
    -2.45340708300902D-01, 2.45340708300900D-01, 7.37473728545390D-01,&
    1.23407621539532D+00, 1.73853771211658D+00, 2.25497400208927D+00,&
    2.78880605842812D+00, 3.34785456738320D+00, 3.94476404011561D+00,&
    4.60368244955074D+00, 5.38748089001122D+00/
DATA W(191:210) /                                          &
    2.22939364553440D-13, 4.39934099227386D-10, 1.08606937076940D-07,&
    7.80255647853246D-06, 2.28338636016365D-04, 3.24377334223799D-03,&
    2.48105208874644D-02, 1.09017206020024D-01, 2.86675505362834D-01,&
    4.62243669600611D-01, 4.62243669600611D-01, 2.86675505362836D-01,&
    1.09017206020025D-01, 2.48105208874642D-02, 3.24377334223800D-03,&
    2.28338636016370D-04, 7.80255647853270D-06, 1.08606937076942D-07,&
    4.39934099227335D-10, 2.22939364553444D-13/

CONTAINS
!-------------------------------------------------------------
! Some getter functions to access the DMA information from C
INTEGER(C_INT) FUNCTION get_nsites() BIND(c, name='get_nsites')
    get_nsites = ns
END FUNCTION get_nsites

INTEGER(C_INT) FUNCTION get_order(site) BIND(c, name='get_order')
    INTEGER(C_INT), VALUE, INTENT(IN) :: site
    get_order = limit(site)
END FUNCTION get_order

REAL(C_DOUBLE) FUNCTION get_dma_value(site, addr) BIND(c, name='get_dma_value')
    INTEGER(C_INT), VALUE, INTENT(IN) :: site, addr
    get_dma_value = q(addr,site)
END FUNCTION get_dma_value

REAL(C_DOUBLE) FUNCTION get_tot_value(addr) BIND(c, name='get_tot_value')
    INTEGER(C_INT), VALUE, INTENT(IN) :: addr
    get_tot_value = qt(addr)
END FUNCTION get_tot_value
!-----------------------------------------------------------------   DMA

SUBROUTINE dma_main(w,kp,infile,outfile)
USE input
IMPLICIT NONE
!-----------------------------------------------------
!     Copyright A J Stone University of Cambridge 1983
!     Modifications and interface to Cadpac, R D Amos
!     Version for Cadpac5 , R D Amos, June 1990
!-----------------------------------------------------
REAL(dp) :: w(*)
INTEGER :: kp, prank=5, infile, outfile

!                      Distributed Multipole Analysis

!  The charge distribution is analysed by treating the overlap density
!  of each pair of primitives as a multipole expansion about the centre
!  of the overlap distribution. This expansion (which terminates) is
!  then moved to the nearest of a list of expansion sites, giving a
!  non-terminating but rapidly convergent expansion. By default, the
!  expansion sites comprise the nuclei, but sites can be deleted or
!  additional ones added. An additional site at the centre of the
!  molecule is recommended for example if there is no atom there. It is
!  possible to specify the maximum rank of multipole which can be
!  generated at any site; multipoles of higher rank which would have
!  been moved to such limited sites are moved to other other sites
!  instead.
!  Normally the SCF density matrix is analysed, but it is possible to
!  specify other density matrices, generated for example by
!  perturbations arising from correlation or external fields, or by
!  CASSCF or MCSCF calculation.

LOGICAL eof

CHARACTER(LEN=8) :: aa, posn
CHARACTER(LEN=16) :: wa, wb, wc

LOGICAL :: check

INTEGER :: iw=6, kr=0, kw=0, nerror=0, itol, zs(MAXS)
INTEGER :: i, j, k, l, m, ok
REAL(dp) :: r, ox=0d0, oy=0d0, oz=0d0

!  Directive syntax:
!   MULTIPOLES
!   [options]
!   START

!  Options specifying the choice of sites should be given first, followed
!  by other options as required. Sites are specified by:

!  ATOMS   (default) Move all contributions to nearest atom.

!  DELETE name
!          Delete all sites with the name given.  DELETE ALL deletes
!          all sites. DELETE CHARGE imol  deletes nuclear charges on
!          atoms in molecule imol=1,2. If only one molecule is present
!          imol can be omitted.

!  ADD name x y z [LIMIT lmax] [RADIUS radius]
!          Add a new site at (x,y,z) with the name specified.  The
!          multipole rank is limited to lmax if a value is specified,
!          and a relative radius can be specified also (see below).

!  The program automatically takes advantage of the geometry of the molecule
!  and the sites chosen to use the fastest algorithm. For calculation of
!  distributed multipoles there is a faster procedure that is suitable
!  when the molecule is linear and arranged parallel to the z axis.
!  However it is possible to override the automatic
!  selection. In particular, it is necessary to specify GENERAL for a linear
!  molecule that is subject to an external field, or is not in a singlet
!  sigma state. If such a molecule is treated as linear, there may be
!  non-vanishing multipole moments which will not be calculated; however
!  the Ql0 will still be correct.

!  GENERAL Ignore any special features of the geometry and use the most
!          general version of the calculation.

!  LIMIT lmax [name name ...]
!          Limit the rank of multipoles on sites with any of the names
!          given to lmax at most.  (Contributions with higher ranks are
!          moved to other sites.) If no name is given the limit applies
!          all sites. Default (and maximum) is 20 for the linear
!          version, 10 otherwise.

!  RADIUS radius name name ...
!          Specify a radius for sites with any of the names given.
!          The actual distances from an overlap centre to the sites are
!          scaled by dividing by the radii of the sites, and
!          the contributions are moved to the site which is closest, in
!          terms of scaled distances, to the overlap centre. The default
!          is that all sites have radius 0.65 angstrom, except for 
!          hydrogen, which has radius 0.325 angstrom.

!  BIGEXP e
!          If the sum of the exponents of a pair of primitives is greater
!          than e, and they are on the same atom, use conventional DMA to
!          handle the multipoles. If it is less than e, evaluate the
!          electron density associated with this pair of primitives and
!          obtain the multipole contributions by integration over a grid.
!          Default value 1.0.

!  PUNCH [file] [RANK prank] [APPEND]
!          Write the multipole moments to the specified file.
!          The moments are given in a form suitable for re-input to the
!          ORIENT program, together with the site names and positions.
!          Multipoles above rank prank (default 5) are not output. APPEND
!          causes the output to be appended to the specified file. 

!  REPORT  Print the multipole contributions of each pair of primitives
!          as the calculation proceeds.

!  NONUCLEAR
!          The nuclear contribution is not to be evaluated.

!  SLICE
!          (For LINEAR molecules only.) Evaluate integrals over slices
!          of space each containing one DMA site. The positions of the
!          dividing planes between the slices are determined by reference
!          to the site radii. In the calculation of DMA matrices for
!          linear molecules this procedure is always used.

lmax=20
iw=outfile
kr=0
kw=0
linear=.false.
planar=.false.
general=.false.
perp=0
itol=18
nuclei=.true.
slice=.false.
spread=1.0d0
bigexp=bigexp_default
lebedev=.true.
nerror=0
iax(1)=0
do i=1,maxbfn
  iax(i+1)=iax(i)+i
end do
if (rfact==1d0) then
  runit="bohr"
else
  runit="angstrom"
end if

if (allocated(limit)) deallocate(limit,xs,radius,q)
allocate(limit(maxs),xs(3,maxs),radius(maxs),q(0:121,maxs),stat=ok)
if (ok>0) call die("Can't allocate site arrays")

!  ATOMS (default choice of sites)
call atom_sites

!  Read input keywords
do
  if (nerror .eq. 1) then
    write(outfile, '(a)') 'Syntax checking only from this point.'
    nerror=nerror+1
  endif
  call read_line(eof, infile)
  call readu(wa)
  select case(wa)
  case('ATOMS')  !  Default choice of sites
    call atom_sites
  case('MODIFY','END','NOTE',"","!")
    !  Ignore
  case('START')
    exit
  case('LINEAR')
    write(outfile, '(/A)') 'The LINEAR option is no longer needed'
  case('PLANAR')
    call readu(wb)
    select case(wb)
    case('X','YZ')
      i=1
      wc='yz'
    case('Y','XZ')
      i=2
      wc='xz'
    case("",'Z','XY')
      i=4
      wc='xy'
    end select
    if (linear .and. (i .eq. 1 .or. i .eq. 2)) then
      linear=.false.
      planar=.true.
      perp=i
    else if (.not. planar .or. perp .ne. i) then
      write (iw,'(3a)')                                          &
          'The sites are not arranged parallel to the ',             &
          WC(1:2),' plane.'
      call die('Calculation abandoned')
    endif

  case('GENERAL')
    !  Override automatic switches
    general=.true.
    linear=.false.
    planar=.false.
    perp=0

  case('SPREAD')
    call readf(spread,rfact)
    if (spread .eq. 0.0d0) spread=1.0d0

  case("GRID")
    do while (item < nitems)
      call readu(wa)
      select case(wa)
      case("RADIAL")
        call readi(n_r)
      case("LEBEDEV")
        call readi(n_a)
        Lebedev=.true.
      case("GAUSS-LEGENDRE", "LEGENDRE")
        Lebedev=.false.
        call readi(n_a)
      case("BECKE","SMOOTHING")
        call readi(k_mu) ! Default is 3
!     case("RADII")
!       call readu(wa)
!       select case(wa)
!       case("SLATER","ON")
!         Slater=.true.
!       case("EQUAL","OFF")
!         Slater=.false.
!       case default
!         call die ("Grid radii option "//trim(wa)//" not recognised",.true.)
!       end select
      case default
        call die ("Grid option "//trim(wa)//" not recognised",.true.)
      end select
    end do

  case('REPORT')
    if (nitems .gt. 1) then
      call readi(kr)
    else
      kr=7
    endif
    if (kr .gt. 1) kw=iw

  case('PUNCH')
    if (kp .eq. 0) then
      posn="rewind"
    else
      posn="append"
    endif
    do while (item .lt. nitems)
      call readu(wa)
      select case(wa)
      case("APPEND")
        posn="append"
      case("RANK")
        call readi(prank)
      case default
        call reread(-1)
        call reada(tempfile)
        if (tempfile .ne. punchfile .and. kp .ne. 0) then
          close(kp)
          kp=0
          posn="rewind"
        endif
        punchfile=tempfile
      end select
    end do
    if (kp .eq. 0) then
      inquire(file=punchfile,exist=check)
      if (check) then
        open(unit=7,file=punchfile,position=posn)
      else
        open(unit=7,file=punchfile)
      endif
      kp=7
    endif

  case('NONUCLEAR')
    nuclei=.false.

  case('DELETE')  !  Delete sites
    call readu(wb)
    if (wb .eq. 'ALL') then
      !  Delete all sites
      ns=0
    else
      !  Delete sites matching name given
      k=0
      call reread(-1)
      call reada(aa)
      do i=1,ns
        if (aa .ne. name(i)) then
          k=k+1
          if (k .ne. i) then
            name(k)=name(i)
            do j=1,3
              xs(j,k)=xs(j,i)
            end do
            limit(k)=limit(i)
            radius(k)=radius(i)
          endif
        endif
      end do
      ns=k
    endif
    !  Recalculate linear/planar switches
    call planes

  case('SLICE')
    slice=.true.

  !  LIMIT limit [site site ...]
  case("LIMIT")
    call readi(l)
    l=min(l,20)
    if (item < nitems) then
      do while (item < nitems)
        call reada(aa)
        do i=1,ns
          if (name(i) .eq. aa) limit(i)=l
        end do
      end do
    else
      lmax=l
      do i=1,ns
        limit(i)=min(lmax,limit(i))
      end do
    endif

  !  RADIUS name radius name radius ...
  case('RADIUS')
    ! check=.false.
    do while (item<nitems)
      call reada(aa)
      call readf(r,rfact)
      do i=1,ns
        if (name(i) .eq. aa) then
          radius(i)=r
          check=.true.
        endif
      end do
      ! if (.not. check) then
      !   print '(3a)', 'Name ', trim(aa), ' not found.'
      !   nerror=nerror+1
      ! endif
    end do

  case('ORIGIN')
    call readf(ox,rfact)
    call readf(oy,rfact)
    call readf(oz,rfact)

  case('ADD')
    !  ADD site x y z [LIMIT limit] [RADIUS radius]
    if (ns .ge. maxs) then
      write (iw,'(a,i4)') 'Too many sites -- maximum is', MAXS
      nerror=nerror+1
    else
      ns=ns+1
    endif
    call reada(aa)
    name(ns)=aa
    do i=1,3
      call readf(xs(i,ns),rfact)
    end do
    limit(ns)=lmax
    radius(ns)=0.65d0/bohr
    zs(ns)=0
    do while (item < nitems)
      call readu(aa)
      select case(aa)
      case("LIMIT")
        call readi(limit(ns))
      case("RADIUS")
        call readf(radius(ns),rfact)
      end select
    end do
!  Recalculate linear/planar switches
    call planes

! case("RADII")
!   call readu(aa)
!   select case(aa)
!   case("EQUAL")
!     do i=1,ns
!       radius(i)=0.5d0/bohr
!     end do
!   case("DEFAULT","SLATER")
!     do i=1,ns
!       radius(i)=slater_radius(zs(i))
!     end do
!   end select

  case("BIGEXP","SWITCH")
    call readf(bigexp)
    if (bigexp < 0d0) then
      write(outfile,"(a)") "Switch value must not be negative"
      nerror=nerror+1
    end if
    if (bigexp > 0d0) general=.true.

  case default
    write(outfile,'(a,a)') 'Unrecognised DMA keyword ', wa
    nerror=nerror+1
  end select
end do

!  START

if (nerror .gt. 0) then
  call die ('Multipoles directive not executed.',.false.)
endif

binom(0,0)=1.0d0
rtbinom(0,0)=1d0
rt(0)=0.0d0
do k=1,20
  rt(k)=sqrt(dble(k))
  binom(k,0)=1.0d0
  rtbinom(k,0)=1.0d0
  binom(k-1,k)=0.0d0
  do m=1,k
    binom(k,m)=binom(k-1,m-1)+binom(k-1,m)
    rtbinom(k,m)=sqrt(binom(k,m))
  end do
end do

if (.not. linear) lmax=min0(lmax,10)

l=0
do i=1,ns
  limit(i)=min0(lmax,limit(i))
  l=max0(l,limit(i))
end do
lmax=l

tol=2.30258d0*itol

!  Standard DMA

write (iw,"(/25x,a/)") "Distributed Multipole Analysis"

Q=0.0d0

if (bigexp > 0d0) then
  write(outfile,"(a,a,f0.5)") "Standard DMA for products of primitives",  &
      " with exponent greater than ", bigexp
  general=.true.
  call make_grid(ns, zs, xs, radius, outfile)
  if (allocated(rho)) then
    deallocate(rho)
  end if
  allocate(rho(ng), stat=ok)
  if (ok>0) then
    write(outfile,"(a)") "Allocation of density grid failed"
    stop
  end if
else
  write(outfile, "(a)") "Standard DMA"
end if

if (general) then
  !  Override automatic switches
  linear=.false.
  planar=.false.
  perp=0
end if

write(outfile, "(/2a)") "Positions and radii in ", trim(runit)

if (Qfactor(0) .ne. 1d0) then
  write(outfile, "(/a)") "Multipole moments are in SI units, multiplied by 10^(10k+20) for rank k"
else
  write(outfile, "(a)") "Multipole moments in atomic units, ea_0^k for rank k"
end if

if (linear) then
  call dmaql0(w,kw)
else
  call dmaqlm(w,kw)
endif

qt=0.0d0

if (kp .gt. 0) then
  !  Heading for 'punch' results
  write (kp,'(2a)') ("! ", trim(title(i)), i=1,2)
  write (kp,"(2a)") "Units ", trim(runit)
endif

do i=1,ns
!  Print results for site I
  write (iw,"(/a8,3(a,f10.6),1x,a/8x,a,i3,a,f7.3,1x,a)")               &
      name(i),                                                         &
      '   x =', xs(1,i)*rfact, '  y =', xs(2,i)*rfact,                 &
      '  z =', xs(3,i)*rfact, trim(runit), &
      '   Maximum rank =', limit(i),                                   &
      '   Radius =', radius(i)*rfact, trim(runit)
  call printq(q(0:,i),limit(i), linear, iw)
  if (kp .gt. 0) then
    !  'Punch' results
    if (linear) then
      write (kp,'(/a8,f16.10,4x,a,i3)')                                &
          name(i), xs(3,i)*rfact, 'Rank', min(prank,limit(i))
      write (kp,'(f15.10)') (q(l,i), l=0,min(prank,limit(i)))
    else
      write (kp,'(/a8,3f16.10,4x,a,i3)')                               &
          name(i), xs(:,i)*rfact, 'Rank', min(prank,limit(i))
      do l=0,min(prank,limit(i))
        write (kp,'(5(f15.10))') (q(k,i), k=l**2+1,(l+1)**2)
      end do
    endif
  endif
!  Shift total for site I to origin
  if (linear) then
    call shiftz(q(0:,i),0,limit(i), qt,lmax, xs(3,i)-oz)
  else
    call shiftq(q(1:,i),0,limit(i), qt(1:),lmax,                       &
        xs(1,i)-ox, xs(2,i)-oy, xs(3,i)-oz)
  endif
end do

!  Print total multipoles
if (linear) then
  write (iw,"(/a,f11.6,1x,a)")                                         &
      'Total multipoles referred to origin at z =', oz*rfact, trim(runit)
else
  write (iw,'(/a/10x,3(a,f11.6),1x,a)')                                &
      "Total multipoles referred to origin at", " x =", ox*rfact,      &
      ',  y = ', oy*rfact, ',  z = ', oz*rfact, trim(runit)
endif
call printq(qt,lmax, linear, iw)

if (bigexp > 0d0) deallocate(rho)

CONTAINS

SUBROUTINE atom_sites

mindc=1
maxdc=maxcen
do i=1,nat
  xs(:,i)=c(:,i)
  limit(i)=lmax
  zs(i)=nint(zan(i))
  if (zs(i)==1) then
    !  Hydrogen atoms have default radius 0.325 Angstrom
    radius(i)=0.325d0/bohr
  else
    !  Otherwise the default radius is 0.65 Angstrom
    radius(i)=0.65d0/bohr
  end if
end do
ns=nat

call planes

END SUBROUTINE atom_sites

END SUBROUTINE dma_main

!---------------------------------------------------------------  DMAQL0

SUBROUTINE dmaql0(densty,iw)
IMPLICIT NONE
!-----------------------------------------------------
!     Copyright A J Stone University of Cambridge 1983
!     Modifications and interface to Cadpac, R D Amos
!     Version for Cadpac5 , R D Amos, June 1990
!-----------------------------------------------------
REAL(dp) :: densty(*)
INTEGER :: iw

!  Calculate multipole moments, and shift them to the nearest site. In
!  this routine, appropriate for linear molecules, only the moments Qlm
!  with m=0 are calculated and stored, and they are held in the order
!  Q0, Q1, Q2, ...
!  If IW>0, print details of the progress of the calculation.
!  If NUCLEI, include the nuclear charges in the calculation.
!  BINOM(K,M) contains the binomial coefficient (K).
!                                               (M)
!  RT(K) contains sqrt(K).

LOGICAL :: ieqj, iieqjj, skip

REAL(dp) :: gx(0:20), gz(0:20,0:5,0:5), sep(0:maxs),                   &
    qt(0:20)
INTEGER :: sort(maxs)

INTEGER :: i, i1, i2, ia, ib, ii, ii1, ii2, ig, iq, iqmin, is,         &
    j, j1, j2, jb, jj, jj1, jj2, jg, jgmax, js, k, k1, k2,             &
    la, lb, loci, locj, m, mx, my, mini, maxi, minj, maxj, nq
REAL(dp) :: aa, ai, arri, aj, ci, cj, dum, fac, g,                     &
    p, paz, pq, pz, ps, rr, s, t,                                      &
    z1, z2, za, zas, zb, zbs, zp, zi, zj, zji, zs

if (slice) then
  lmax=min(lmax,11)
  do is=1,ns
    limit(is)=min(limit(is),lmax)
  end do
  !  Sort DMA sites into order of increasing z.
  do is=1,ns
    sort(is)=is
  end do
  do is=2,ns
    js=is
    do while (js>1)
      if (xs(3,sort(js)) .ge. xs(3,sort(js-1))) exit
      m=sort(js)
      sort(js)=sort(js-1)
      sort(js-1)=m
      js=js-1
    end do
  end do
!  Coordinates of separating planes
  sep(0)=-1.0d6
  do is=1,ns-1
    i1=sort(is)
    i2=sort(is+1)
    sep(is)=(radius(i1)*xs(3,i2)+radius(i2)*xs(3,i1))                  &
        /(radius(i1)+radius(i2))
  end do
  sep(ns)=1.0d6
  if (iw .gt. 0) write (iw,"(a//(28x, f14.4 / 1x, i4, i8, f10.4))")    &
      ' slice        site               separator',                    &
      sep(0), (is, sort(is), xs(3,sort(is)), sep(is), is=1,ns)
  if (iw .gt. 0) write (iw,'(a/a/)')                                   &
      '   Atoms   Shells Primitives  Position',                        &
      'Slice     Multipole contributions ...'
else
  if (iw .gt. 0) write (iw,'(/a,a/)')                                  &
      '    Atoms   Shells Primitives  Position',                       &
      '  Multipole contributions ...'
endif

!  Clear temporary multipole array. Note that subsequently it is cleared
!  whenever multipoles are moved from it to an expansion site.
qt(0:lmax)=0.0d0

!  Loop over pairs of atoms
katom(nshell+1)=0
do i=1,nat
  zi=c(3,i)
  if (nuclei .and. i .ge. mindc .and. i .le. maxdc) then
    qt(0)=zan(i)
    call movez(qt, zi, iw)
  endif
  !  Find shells for atom i
  ii1=0
  ishells: do ii=1,nshell
    if (katom(ii) .eq. i) then
      ii1=ii
      do k=ii1,nshell+1
        if (katom(k) .ne. i) then
          ii2=k-1
          exit ishells
        end if
      end do
    end if
  end do ishells
  if (ii1 .eq. 0) cycle  !  No basis functions on atom I

  !  Loop over shells for atom i
  do ii=ii1,ii2
    i1=kstart(ii)
    i2=i1+kng(ii)-1
    la=ktype(ii)-1
    mini=kmin(ii)
    maxi=kmax(ii)
    loci=kloc(ii)-mini

    !  Loop over atoms j
    do j=1,i
      ieqj=i .eq. j
      zj=c(3,j)
      zji=zi-zj
      rr=zji**2

      !  Find shells for atom j
      jj1=0
      jshells: do jj=1,nshell
        if (katom(jj) .eq. j) then
          jj1=jj
          do k=jj1,nshell+1
            if (katom(k) .ne. j) then
              jj2=k-1
              exit jshells
            end if
          end do
        end if
      end do jshells
      if (jj1 .eq. 0) cycle  !  No basis functions on atom j

      !  Loop over shells for atom j
      if (ieqj) jj2=ii
      do jj=jj1,jj2
        j1=kstart(jj)
        j2=j1+kng(jj)-1
        lb=ktype(jj)-1
        minj=kmin(jj)
        maxj=kmax(jj)
        locj=kloc(jj)-minj
        iieqjj=ii.eq.jj
        !  Set up temporary density matrix for this pair of shells
        do ib=mini,maxi
          m=iax(loci+ib)
          if (iieqjj) then
            do jb=minj,ib
              d(ib,jb)=densty(m+locj+jb)
              d(jb,ib)=densty(m+locj+jb)
            end do
          else
            do jb=minj,maxj
              d(ib,jb)=densty(m+locj+jb)
            end do
          end if
        end do
        !  Insert factors of sqrt(3) for xy, xz and yz functions, if present,
        !  and factors of sqrt(5) and sqrt(15) for f functions. This is to
        !  compensate for the use of the same normalisation factor in all the
        !  d integrals and in all the f integrals.
        select case(la)
        case(2)
          d(8:10,minj:maxj)=rt(3)*d(8:10,minj:maxj)
        case(3)
          d(14:19,minj:maxj)=rt(5)*d(14:19,minj:maxj)
          d(20,minj:maxj)=rt(15)*d(20,minj:maxj)
        case(4)
          d(24:29,minj:maxj)=rt(7)*d(24:29,minj:maxj)
          d(30:32,minj:maxj)=(rt(5)*rt(7)/rt(3))*d(30:32,minj:maxj)
          d(33:35,minj:maxj)=rt(5)*rt(7)*d(33:35,minj:maxj)
        end select
        select case(lb)
        case(2)
          d(mini:maxi,8:10)=rt(3)*d(mini:maxi,8:10)
        case(3)
          d(mini:maxi,14:19)=rt(5)*d(mini:maxi,14:19)
          d(mini:maxi,20)=rt(15)*d(mini:maxi,20)
        case(4)
          d(mini:maxi,24:29)=rt(7)*d(mini:maxi,24:29)
          d(mini:maxi,30:32)=(rt(5)*rt(7)/rt(3))*d(mini:maxi,30:32)
          d(mini:maxi,33:35)=rt(5)*rt(7)*d(mini:maxi,33:35)
        end select

        !  I primitive
        jgmax=j2
        do ig=i1,i2
          ai=ex(ig)
          arri=ai*rr
          !  J primitive
          if (iieqjj) jgmax=ig
          do jg=j1,jgmax
            aj=ex(jg)
            aa=ai+aj
            dum=aj*arri/aa
            if (dum.gt.tol) cycle
            fac=dexp(-dum)
!  Introduce factor of 2 if (a) shells are different (II.NE.JJ) because
!  loops over atoms and shells count each pair only once, and (b) if
!  II.EQ.JJ but IG.NE.JG, because loops over primitives count each pair
!  only once.  In fact IG.NE.JG covers both cases.  However, different
!  atoms may use the same shells and primitives if there is symmetry, an
!  there is always a factor of 2 if the atoms are different.
            if (ig .ne. jg .or. .not. ieqj) fac=2.0d0*fac
            !  ZP is the position of the overlap centre.
            p=aj/aa
            zp=zi-p*zji
            za=zi-zp
            zb=zj-zp
            if (iw .gt. 0) write (iw,"(3(i5,i4), f11.4)")              &
                i,j, ii,jj, ig,jg, zp
!  Use numerical integration to evaluate the multipole integrals
!  over x and y (they are the same).
            t=dsqrt(1.0d0/aa)
!  The x integrals involve polynomials up to order LA+LB+LMAX,
!  for which NQ+1 integration points are required.
            if (slice) then
              nq=(la+lb+lmax)/2
            else
              nq=la+lb
            endif
            k1=mink(nq+1)
            k2=maxk(nq+1)
            gx(0:20)=0.0d0
!  In GX it is only necessary to accumulate the sums of even powers of
!  xk:
!           sum(k) gk * xk**(IQ),  IQ even,
!  since xA=xB=0. The same values serve for GY.
            do k=k1,k2
              s=h(k)*t
              g=w(k)*t
              ps=g
              iqmin=min((la+lb+lmax),20)
              do iq=0,iqmin,2
                gx(iq)=gx(iq)+ps
                ps=ps*(s**2)
              end do
            end do

            if (slice) then

!  Use the analytical formula, and recursion formulae, to evaluate
!  integrals over slices of space each containing one DMA site.
!  Loop over slices
              do is=1,ns
                z1=sep(is-1)-zp
                z2=sep(is)-zp
                call dmaerf(aa,la,lb, za,zb, z1,z2, gz, skip)
!  If SKIP, all the integrals in GZ are negligible.
                if (skip) cycle

!  Now these basic integrals are used to construct the multipole moments
!  for the overlap density corresponding to each pair of basis functions
!  in the pair of shells.
                qt(0:lmax)=0.0d0
                ci=cs(ig)
                do ia=mini,maxi
                  if (ia .gt. 1 .and. la .eq. 1) ci=cp(ig)
                  ! if (ia .gt. 4) ci=cd(ig)
                  cj=cs(jg)
                  do jb=minj,maxj
                    if (jb .gt. 1 .and. lb .eq. 1) cj=cp(jg)
                    ! if (jb .gt. 4) cj=cd(jg)
!  The integral of (x**i)*(y**j)*(z**k) over the current pair of atomic
!  primitive orbitals is F*GX(MX+i)*GY(MY+j)*GZ(k,IZ(IA),IZ(IB)), and is
!  zero unless both MX and MY are even.
                    mx=ix(ia)+ix(jb)
                    my=iy(ia)+iy(jb)
                    if (mod(mx,2) .eq. 0 .and. mod(my,2) .eq. 0)       &
                        call addql0(qt, lmax, -fac*ci*cj*d(ia,jb),     &
                        gx(mx),gx(my),gz(0,iz(ia),iz(jb)))
!  End of loop over basis functions
                  end do
                end do
                if (iw .gt. 0) write (iw,"(a,i3, 11f11.7)")            &
                    "slice", is, (qt(iq), iq=0,min(lmax,10))
!  Move multipoles to expansion site contained in this slice.
!  Note that they are currently referred to the overlap centre P.
                zs=xs(3,sort(is))
                call shiftz(qt,0,lmax, q(0:,sort(is)),lmax, zp-zs)
!  End of loop over slices
              end do


            else  !  if not slice


!  Clear integral arrays
              gz(0:20,0:5,0:5)=0.0d0

!  The following loop runs through the integration points, accumulating
!  in GZ(IQ,IA,IB) the quantity
!       sum(k) gk * (zk-zA)**IA * (zk-zB)**IB * zk**IQ
!  where gk=w(k)/sqrt(AA) and zk=h(k)/sqrt(AA).
              do k=k1,k2
                s=h(k)*t
                g=w(k)*t
                zas=s-za
                zbs=s-zb
                paz=g
                do ia=0,la
                  pz=paz
                  do ib=0,lb
                    pq=pz
                    do iq=0,la+lb
                      gz(iq,ia,ib)=gz(iq,ia,ib)+pq
                      pq=pq*s
                    end do
                    pz=pz*zbs
                  end do
                  paz=paz*zas
                end do
              end do

!  Now these basic integrals are used to construct the multipole moments
!  for the overlap density corresponding to each pair of basis functions
!  in the pair of shells.
              qt(0:lmax)=0.0d0
              ci=cs(ig)
              do ia=mini,maxi
                if (ia .gt. 1 .and. la .eq. 1) ci=cp(ig)
                ! if (ia .gt. 4) ci=cd(ig)
                cj=cs(jg)
                do jb=minj,maxj
                  if (jb .gt. 1 .and. lb .eq. 1) cj=cp(jg)
                  ! if (jb .gt. 4) cj=cd(jg)
!  The integral of (x**i)*(y**j)*(z**k) over the current pair of atomic
!  primitive orbitals is F*GX(MX+i)*GY(MY+j)*GZ(k,IZ(IA),IZ(JB)), and
!  is zero unless both MX and MY are even.
                  mx=ix(ia)+ix(jb)
                  my=iy(ia)+iy(jb)
                  if (mod(mx,2) .eq. 0 .and. mod(my,2) .eq. 0)         &
                      call addql0 (qt, min(nq,lmax), -fac*ci*cj*d(ia,jb), &
                      gx(mx),gx(my),gz(0,iz(ia),iz(jb)))
                  if (iw .gt. 0) write (iw,"(1p,3e10.2,2i3,1p,4e10.2)") &
                      fac, ci, cj, ia, jb, d(ia,jb),                   &
                      gx(mx), gx(my), gz(0,iz(ia),iz(jb))
!  End of loop over basis functions
                end do
              end do
              if (iw .gt. 0) write (iw,"(3(i5,i4), f11.4, 3x, 11f12.8)") &
                  i,j, ii,jj, ig,jg, zp, (qt(iq), iq=0,10)
!  Move multipoles to nearest site.
              call movez(qt, zp, iw)
            endif
!  End of loop over primitives
          end do
        end do
        if (iw .gt. 0) write (iw,"(/(5f15.8))")                        &
            (q(0:4,ia), ia=1,ns)
        if (iw .gt. 0) write (iw,"(1x)")
!  End of loop over shells
      end do
    end do
  end do
!  End of loop over atoms
end do

END SUBROUTINE dmaql0

!-----------------------------------------------------------------ADDQL0

SUBROUTINE ADDQL0(Q, L, F, GX,GY,GZ)
IMPLICIT NONE
!-----------------------------------------------------
!     Copyright A J Stone University of Cambridge 1983
!     Modifications and interface to Cadpac, R D Amos
!     Version for Cadpac5 , R D Amos, June 1990
!-----------------------------------------------------
REAL(dp),INTENT(IN) :: GX(0:*), GY(0:*), GZ(0:*), F
INTEGER,INTENT(IN) :: L
REAL(dp),INTENT(INOUT) :: Q(0:*)

REAL(dp) :: xy0, xy2, xy4, xy6, xy8, xy10

GO TO (100,101,102,103,104,105,106,107,108,109,110,111), L+1
111   CONTINUE
110   XY10=GX(10)*GY(0)+5.0D0*GX(8)*GY(2)+10.0D0*GX(6)*GY(4)+           &
                  10.0D0*GX(4)*GY(6)+5.0D0*GX(2)*GY(8)+GX(0)*GY(10)
109   CONTINUE
108   XY8=GX(8)*GY(0)+4.0D0*GX(6)*GY(2)+6.0D0*GX(4)*GY(4)+              &
                                4.0D0*GX(2)*GY(6)+GX(0)*GY(8)
107   CONTINUE
106   XY6=GX(6)*GY(0)+3.0D0*GX(4)*GY(2)+3.0D0*GX(2)*GY(4)+GX(0)*GY(6)
105   CONTINUE
104   XY4=GX(4)*GY(0)+2.0D0*GX(2)*GY(2)+GX(0)*GY(4)
103   CONTINUE
102   XY2=GX(2)*GY(0)+GX(0)*GY(2)
101   CONTINUE
100   XY0=GX(0)*GY(0)
GO TO (200,201,202,203,204,205,206,207,208,209,210,211), L+1
211   Q(11)=Q(11)+0.00390625D0*F*                                       &
          (256.0D0*XY0*GZ(11)                                           &
          -7040.0D0*XY2*GZ(9)                                           &
          +31680.0D0*XY4*GZ(7)                                          &
          -36960.0D0*XY6*GZ(5)                                          &
          +11550.0D0*XY8*GZ(3)                                          &
          -693.0D0*XY10*GZ(1))
210   Q(10)=Q(10)+0.00390625D0*F*                                       &
          (256.0D0*XY0*GZ(10)                                           &
          -5760.0D0*XY2*GZ(8)                                           &
          +20160.0D0*XY4*GZ(6)                                          &
          -16800.0D0*XY6*GZ(4)                                          &
          +3150.0D0*XY8*GZ(2)                                           &
          -63.0D0*XY10*GZ(0))
209   Q(9)=Q(9)+0.0078125D0*F*                                          &
          (128.0D0*XY0*GZ(9)                                            &
          -2304.0D0*XY2*GZ(7)                                           &
          +6048.0D0*XY4*GZ(5)                                           &
          -3360.0D0*XY6*GZ(3)                                           &
          +315.0D0*XY8*GZ(1))
208   Q(8)=Q(8)+0.0078125D0*F*                                          &
          (128.0D0*XY0*GZ(8)                                            &
          -1792.0D0*XY2*GZ(6)                                           &
          +3360.0D0*XY4*GZ(4)                                           &
          -1120.0D0*XY6*GZ(2)                                           &
          +35.0D0*XY8*GZ(0))
207   Q(7)=Q(7)+0.0625D0*F*                                             &
          (16.0D0*XY0*GZ(7)                                             &
          -168.0D0*XY2*GZ(5)                                            &
          +210.0D0*XY4*GZ(3)                                            &
          -35.0D0*XY6*GZ(1))
206   Q(6)=Q(6)+0.0625D0*F                                              &
          *(16.0D0*XY0*GZ(6)                                            &
           -120.0D0*XY2*GZ(4)                                           &
           +90.0D0*XY4*GZ(2)                                            &
           -5.0D0*XY6*GZ(0))
205   Q(5)=Q(5)+0.125D0*F                                               &
          *(8.0D0*GZ(5)*XY0                                             &
           -40.0D0*GZ(3)*XY2                                            &
           +15.0D0*GZ(1)*XY4)
204   Q(4)=Q(4)+0.125D0*F*(8.0D0*XY0*GZ(4)                              &
                       -24.0D0*XY2*GZ(2)                                &
                        +3.0D0*XY4*GZ(0))
203   Q(3)=Q(3)+F*(XY0*GZ(3)-1.5D0*XY2*GZ(1))
202   Q(2)=Q(2)+F*(XY0*GZ(2)-0.5D0*XY2*GZ(0))
201   Q(1)=Q(1)+F*XY0*GZ(1)
200   Q(0)=Q(0)+F*XY0*GZ(0)

END SUBROUTINE ADDQL0

!-----------------------------------------------------------------SHIFTZ

SUBROUTINE shiftz(q1,l1,m1, q2,m2, z)
IMPLICIT NONE

!-----------------------------------------------------
!     Copyright A J Stone University of Cambridge 1983
!     Modifications and interface to Cadpac, R D Amos
!     Version for Cadpac5 , R D Amos, June 1990
!-----------------------------------------------------
!  Shift those multipoles in Q1 with L1 <= n <= M1 to the origin to
!  which the multipoles Q2 are referred, and add the shifted values
!  to Q2.
!  Values are required in Q2 only up to n=M2.
!  Z gives the position of Q1 relative to an origin at Q2.
!  Qn(0) = sum(s) [ (n!/s!(n-s)!) * Qs(z) * z**(n-s) ]

INTEGER, INTENT(IN) :: l1, m1, m2
REAL(dp), INTENT(IN) :: q1(0:20), z
REAL(dp), INTENT(INOUT) :: q2(0:20)

REAL(dp) :: zs(0:20)
INTEGER :: i, n, s, s1, s2

if (l1 .gt. m1) return
zs(0)=1.0d0
do i=1,m2
  zs(i)=z*zs(i-1)
end do
!  Charge on Q1
if (l1 .eq. 0) then
  q2(0)=q2(0)+q1(0)
  do n=1,m2
    q2(n)=q2(n)+q1(0)*zs(n)
  end do
end if
if (m1 .le. 0) return
s1=max0(1,l1)
s2=min0(m1,m2)
do s=s1,s2
  q2(s)=q2(s)+q1(s)
  do n=s+1,m2
    q2(n)=q2(n)+binom(n,s)*q1(s)*zs(n-s)
  end do
end do

END SUBROUTINE shiftz

!----------------------------------------------------------------- MOVEZ

SUBROUTINE movez (qp, p, iw)
IMPLICIT NONE
!-----------------------------------------------------
!     Copyright A J Stone University of Cambridge 1983
!     Modifications and interface to Cadpac, R D Amos
!     Version for Cadpac5 , R D Amos, June 1990
!-----------------------------------------------------
!  Move the multipole contributions in QP, which are referred to z=P,
!  to the nearest site.
REAL(dp) :: qp(0:), p
INTEGER, INTENT(IN) :: iw

REAL(dp) :: r(maxs)
REAL(dp), PARAMETER :: eps=1d-8
INTEGER :: i, j, k, l, n, M(2), low

!  J is one of the sites with the largest value of LIMIT.
J=1
do i=1,ns
  r(i)=dabs(xs(3,i)-p)/radius(i)
  if (limit(i) .gt. limit(j)) j=i
end do
!  LOW is the lowest rank of multipole to be moved at this stage
low=0
do
  k=j
  do i=1,ns
    if (r(i) .lt. r(k) .and. limit(i) .ge. low) k=i
  end do
  n=1
  m(1)=k
  do i=1,ns
    if (r(i) .gt. r(k)+eps .or. i .eq. k .or.                         &
        limit(i) .lt. low .or. limit(i) .ne. limit(k)) cycle
    n=2
    m(2)=i
  end do
!  Multipoles of ranks LOW to LIMIT(K) are to be moved at this stage
  IF (IW .GT. 0) WRITE (6,1001) P, LOW, LIMIT(K),                   &
      (M(I), XS(3,M(I)), I=1,N)
1001  FORMAT (' From', F7.3, ': ranks', I3, ' to', I3,                  &
          ' to be moved to site ', I1, ' at', F7.3:                     &
          ' and site', I2, ' at', F7.3)
  if (n .eq. 2) then
    do i=low,limit(k)
      qp(i)=0.5d0*qp(i)
    end do
  end if
  do i=1,n
    k=m(i)
    call shiftz(qp,low,limit(k), q(0:,k),lmax, p-xs(3,k))
  end do
  do i=low,limit(k)
    qp(i)=0.0d0
  end do
!  At this point the new sites carry multipoles of all ranks up to
!  LMAX.  Shift any of ranks higher than LIMIT(K) back to the original
!  origin.  Note that following this process QP carries multipoles of
!  all ranks from LIMIT(K)+1 to LMAX.
  if (limit(k) .ge. lmax) return
  do i=1,n
    k=m(i)
    call shiftz(q(0:,k),limit(k)+1,lmax, qp,lmax, xs(3,k)-p)
    do l=limit(k)+1,lmax
      q(l,k)=0.0d0
    end do
  end do
  low=limit(k)+1
end do
END SUBROUTINE movez

!------------------------------------------------------------------------

REAL(dp) FUNCTION derf(x)
IMPLICIT NONE

REAL(dp), INTENT(IN) :: x
!  Error function
REAL(dp) :: c(18) =                                                    &
    (/1.9449071068178803d0,4.20186582324414d-2,-1.86866103976769d-2,   &
    5.1281061839107d-3,-1.0683107461726d-3,1.744737872522d-4,          &
    -2.15642065714d-5,1.7282657974d-6,-2.00479241d-8,                  &
    -1.64782105d-8,2.0008475d-9,2.57716d-11,-3.06343d-11,              &
    1.9158d-12,3.703d-13,-5.43d-14,-4.0d-15,1.2d-15/)
REAL(dp) :: d(17) =                                                    &
    (/1.4831105640848036d0,-3.010710733865950d-1,6.89948306898316d-2,  &
    -1.39162712647222d-2,2.4207995224335d-3,-3.658639685849d-4,        &
    4.86209844323d-5,-5.7492565580d-6,6.113243578d-7,                  &
    -5.89910153d-8,5.2070091d-9,-4.232976d-10,3.18811d-11,             &
    -2.2361d-12,1.467d-13,-9.0d-15,5.0d-16/)
INTEGER ::  ncfc=18, ncfd=17
REAL(dp) :: xup=6.25d0,sqrtpi=1.7724538509055160d0
REAL(dp), PARAMETER :: zero=0.0d0, one=1.0d0, two=2.0d0,               &
    three=3.0d0, twenty=20.0d0, half = 0.5d0

REAL(dp) :: bj, bjp1, bjp2, x2, xv
INTEGER :: j

xv = dabs(x)
if (xv.ge.xup) then
  derf = dsign(one,x)
else if (xv.le.two) then
  x2 = x*x - two
  bjp2 = zero
  bjp1 = d(ncfd)
  j = ncfd - 1
  do
    bj = x2*bjp1 - bjp2 + d(j)
    if (j.eq.1) exit
    bjp2 = bjp1
    bjp1 = bj
    j = j - 1
  end do
  derf = half*(bj-bjp2)*x
else
  x2 = two - twenty/(xv+three)
  bjp2 = zero
  bjp1 = c(ncfc)
  j = ncfc - 1
  do
    bj = x2*bjp1 - bjp2 + c(j)
    if (j.eq.1) exit
    bjp2 = bjp1
    bjp1 = bj
    j = j - 1
  end do
  x2 = half*(bj-bjp2)/xv*dexp(-x*x)/sqrtpi
  derf = (one-x2)*dsign(one,x)
end IF
END FUNCTION derf

!-----------------------------------------------------------------DMAERF

SUBROUTINE dmaerf(aa, la,lb, za,zb, z1,z2, gz, skip)
!-----------------------------------------------------
!     Copyright A J Stone University of Cambridge 1983
!     Modifications and interface to Cadpac, R D Amos
!     Version for Cadpac5 , R D Amos, June 1990
!-----------------------------------------------------
!  C1  modified to remove the calculation of integrals of powers of x,
!      which can be done more generally and efficiently using the usual
!      numerical method.

!  Returns in GZ(n,k,l) the integral (Z1,Z2) of

!                 z**n (z-zA)**k (z-zB)**l exp(-AA*z**2).
IMPLICIT NONE

REAL(dp), INTENT(IN) :: aa, za, zb, z1, z2
INTEGER, INTENT(IN) :: la, lb
REAL(dp), INTENT(OUT) :: gz(0:20,0:5,0:5)
LOGICAL :: skip

REAL(dp) :: t, v, e1, e2, p1, p2
INTEGER :: k, l, n

v=dsqrt(aa)
t=rtpi/v

gz(0,0,0)=0.5d0*t*(derf(v*z2)-derf(v*z1))
skip=(dabs(gz(0,0,0)) .lt. 1.0d-8)
e1=dexp(-aa*z1**2)/(2d0*aa)
e2=dexp(-aa*z2**2)/(2d0*aa)
gz(1,0,0)=e1-e2
skip=skip .and. (dabs(gz(1,0,0)) .lt. 1.0d-8)
p1=e1
p2=e2
do n=2,16
  p1=p1*z1
  p2=p2*z2
  gz(n,0,0)=p1-p2+(n-1)*gz(n-2,0,0)/(2d0*aa)
  skip=skip .and. (dabs(gz(n,0,0)) .lt. 1.0d-8)
end do
if (skip) return

do n=15,0,-1
  do l=1,min(lb,16-n)
    gz(n,0,l)=gz(n+1,0,l-1)-zb*gz(n,0,l-1)
  end do
  do k=1,min(la,16-n)
    gz(n,k,0)=gz(n+1,k-1,0)-za*gz(n,k-1,0)
    do l=1,min(lb,16-n-la)
      gz(n,k,l)=gz(n+1,k,l-1)-zb*gz(n,k,l-1)
    end do
  end do
end do

END SUBROUTINE dmaerf

!-----------------------------------------------------------------DMAQLM

SUBROUTINE dmaqlm(densty,iw)
IMPLICIT NONE
REAL(dp), INTENT(IN) :: densty(*)
INTEGER, INTENT(IN) :: iw

!-----------------------------------------------------
!     Copyright A J Stone University of Cambridge 1983
!     Modifications and interface to Cadpac, R D Amos
!     Version for Cadpac5 , R D Amos, June 1990
!-----------------------------------------------------
!  Calculate multipole moments, and shift them to the nearest site.
!  Normal routine, for use with non-linear molecules.
!  If IW>0, print details of the progress of the calculation.
!  If NUCLEI, include the nuclear charges in the calculation.
!  The multipoles Qlmc and Qlms are sqrt(2) times the real and
!  imaginary parts of the complex multipole [Ql,-m]*=(-1)**m Qlm,
!  except for m=0.
!  Qlm=<Rlm*p>, where p is the charge density operator and the Rlm
!  are
!           R00  = 1
!           R10  = z
!           R11c = x
!           R11s = y
!           R20  = (1/2)*(3*z**2-r**2)
!           R21c = sqrt(3)*x*z
!           R21s = sqrt(3)*y*z
!           R22c = sqrt(3/4)*(x**2-y**2)
!           R22s = sqrt(3)*x*y
!           etc.
!  The multipoles are stored in the order Q00, Q10, Q11c, Q11s, Q20, ...
!  The program calculates multipoles up to rank 6 for each pair of
!  primitives. This is enough to cope with s, p, d and f functions.
!  Multipoles up to rank LMAX are retained when shifting to a new
!  origin. The maximum value allowed for LMAX is currently 10.
!  The basis functions are:
!               1   2   3   4   5   6   7   8   9  10
!               s   px  py  pz dxx dyy dzz dxy dxz dyz
!               11  12  13  14  15  16  17  18  19  20
!               xxx yyy zzz xxy xxz xyy yyz xzz yzz xyz
!  Note that integrals involving dxy, dxz and dyz must be multiplied by
!  sqrt(3), because of the difference in normalization from dxx etc.
!  Similarly those involving the f functions xxy,xxz,xyy,yyz,xzz and yzz
!  must be multiplied by sqrt(5) and those involving xyz by sqrt(15).


REAL(dp) :: qt
COMMON/BIG/qt(121)

!  L1=L+1, LL1=2L+1, where L is the maximum angular momentum of basis
!  function that can be handled, i.e. 4 at present.
INTEGER, PARAMETER :: L1=5, LL1=9
!  Dimension is (2L+1)*(L+1)**2, where L is the highest angular momentum
!  basis function
REAL(dp) :: gx(0:l1*l1*ll1), gy(0:l1*l1*ll1), gz(0:l1*l1*ll1)
REAL(dp) :: xai(0:5)=1d0, yai(0:5)=1d0, zai(0:5)=1d0,                  &
    xbj(0:5)=1d0, ybj(0:5)=1d0, zbj(0:5)=1d0

LOGICAL :: ieqj, iieqjj, do_quadrature

INTEGER :: i, i1, i2, ia, ib, ii, ii1, ii2, ig, iq,                    &
    j, j1, j2, jb, jj, jj1, jj2, jg, jgmax, k, k1, k2,                 &
    l, la, lb, loci, locj, lq,                                  &
    m, ma, mb, mx, my, mz, mini, maxi, minj, maxj, nq
REAL(dp) :: aa, ai, arri, aj, ci, cj, dum, e, fac, f, g,               &
    p, pax, pay, paz, pq, px, py, pz, rr, s, t, xi, xj, xa, xb,        &
    xji, xk, xp, xas, xbs, yi, yj, ya, yb, yji, yk, yp, yas, ybs,      &
    za, zas, zb, zbs, zk, zp, zi, zj, zji

do_quadrature=.false.

if (iw .gt. 0) write (iw,"(a/a)")                                      &
    '    Atoms   Shells Primitives            Position',               &
    '     Multipole contributions ...'
katom(nshell+1)=0
!  Loop over pairs of atoms
do i=1,nat
  xi=c(1,i)
  yi=c(2,i)
  zi=c(3,i)
  qt=0.0d0
  if (nuclei .and. i .ge. mindc .and. i .le. maxdc) then
    qt(1)=zan(i)
    call moveq (xi,yi,zi)
  endif

  !  Find shells for atom i
  ii1=0
  ishells: do ii=1,nshell
    if (katom(ii) .eq. i) then
      ii1=ii
      do k=ii1,nshell+1
        if (katom(k) .ne. i) then
          ii2=k-1
          exit ishells
        end if
      end do
    end if
  end do ishells
  if (ii1 .eq. 0) cycle  !  No basis functions on atom I

  !  Loop over shells for atom i
  do ii=ii1,ii2
    i1=kstart(ii)
    i2=i1+kng(ii)-1
    la=ktype(ii)-1
    mini=kmin(ii)
    maxi=kmax(ii)
    loci=kloc(ii)-mini

    !  Loop over atoms j
    do j=1,i
      ieqj=i .eq. j
      xj=c(1,j)
      yj=c(2,j)
      zj=c(3,j)
      xji=xi-xj
      yji=yi-yj
      zji=zi-zj
      rr=xji**2+yji**2+zji**2

      !  Find shells for atom j
      jj1=0
      jshells: do jj=1,nshell
        if (katom(jj) .eq. j) then
          jj1=jj
          do k=jj1,nshell+1
            if (katom(k) .ne. j) then
              jj2=k-1
              exit jshells
            end if
          end do
        end if
      end do jshells
      if (jj1 .eq. 0) cycle  !  No basis functions on atom j

      !  Loop over shells for atom j
      if (ieqj) jj2=ii
      do jj=jj1,jj2
        j1=kstart(jj)
        j2=j1+kng(jj)-1
        lb=ktype(jj)-1
        minj=kmin(jj)
        maxj=kmax(jj)
        locj=kloc(jj)-minj
        iieqjj=ii.eq.jj

        !  Set up temporary density matrix for this pair of shells
        do ib=mini,maxi
          m=iax(loci+ib)
          if (iieqjj) then
            do jb=minj,ib
              d(ib,jb)=densty(m+locj+jb)
              d(jb,ib)=densty(m+locj+jb)
            end do
          else
            do jb=minj,maxj
              d(ib,jb)=densty(m+locj+jb)
            end do
          end if
        end do
        !  Insert factors of sqrt(3) for xy, xz and yz functions, if present,
        !  and factors of sqrt(5) and sqrt(15) for f functions. This is to
        !  compensate for the use of the same normalisation factor in all the
        !  d integrals and in all the f integrals.
        select case(la)
        case(2)
          d(8:10,minj:maxj)=rt(3)*d(8:10,minj:maxj)
        case(3)
          d(14:19,minj:maxj)=rt(5)*d(14:19,minj:maxj)
          d(20,minj:maxj)=rt(15)*d(20,minj:maxj)
        case(4)
          d(24:29,minj:maxj)=rt(7)*d(24:29,minj:maxj)
          d(30:32,minj:maxj)=(rt(5)*rt(7)/rt(3))*d(30:32,minj:maxj)
          d(33:35,minj:maxj)=rt(5)*rt(7)*d(33:35,minj:maxj)
        end select
        select case(lb)
        case(2)
          d(mini:maxi,8:10)=rt(3)*d(mini:maxi,8:10)
        case(3)
          d(mini:maxi,14:19)=rt(5)*d(mini:maxi,14:19)
          d(mini:maxi,20)=rt(15)*d(mini:maxi,20)
        case(4)
          d(mini:maxi,24:29)=rt(7)*d(mini:maxi,24:29)
          d(mini:maxi,30:32)=(rt(5)*rt(7)/rt(3))*d(mini:maxi,30:32)
          d(mini:maxi,33:35)=rt(5)*rt(7)*d(mini:maxi,33:35)
        end select

        !  I primitive
        jgmax=j2
        do ig=i1,i2
          ai=ex(ig)
          arri=ai*rr
          !  J primitive
          if (iieqjj) jgmax=ig
          do jg=j1,jgmax
            aj=ex(jg)
            aa=ai+aj
            dum=aj*arri/aa
            if (dum.gt.tol) cycle
            fac=dexp(-dum)
!  Introduce factor of 2 if (a) shells are different (II.NE.JJ) because
!  loops over atoms and shells count each pair only once, and (b) if
!  II.EQ.JJ but IG.NE.JG, because loops over primitives count each pair
!  only once.  In fact IG.NE.JG covers both cases.  However, different
!  atoms may use the same shells and primitives if there is symmetry, and
!  there is always a factor of 2 if the atoms are different.
            if (ig .ne. jg .or. .not. ieqj) fac=2.0d0*fac

            if (ex(ig)+ex(jg) > bigexp) then
!  At least one primitive with large exponent.
!           if (ieqj .and. ex(ig)+ex(jg) > bigexp) then
!  Both primitives on the same atom, at least one with large exponent.
!  Use original DMA algorithm.

!  Vectors A and B are the positions of atoms I and J relative to the
!  centre P of the overlap distribution.
            p=aj/aa
            xa=p*xji
            ya=p*yji
            za=p*zji
            xb=xa-xji
            yb=ya-yji
            zb=za-zji
            xp=xi-xa
            yp=yi-ya
            zp=zi-za
            t=dsqrt(1.0d0/aa)
            if (iw .gt. 0) write (iw,"(3(i5,i4), 3x, 3f10.5)")      &
                i,j, ii,jj, ig,jg, xp,yp,zp
!  LQ is the maximum rank of multipole to which these functions
!  contribute. The integrals involve polynomials up to order 2(NQ-1),
!  for which NQ integration points are required.
            lq=la+lb
            nq=lq+1
            k1=mink(nq)
            k2=maxk(nq)
            !  Clear integral arrays
            gx=0.0d0
            gy=0.0d0
            gz=0.0d0

!  The following loop runs through the integration points, accumulating
!  in GX(LL1*L1*(IA-1)+L1*(IB-1)+IQ) the quantity
!       sum(k) gk * (xk-xA)**(IA-1) * (xk-xB)**(IB-1) * xk**(IQ-1)
!  where gk=w(k)/sqrt(AA) and xk=h(k)/sqrt(AA), and L1=L+1 and
!  LL1=2L+1, where L is the maximum angular momentum to be handled,
!  i.e. 3 at present.
!  Similar expressions for y and z integrals are formed in GY and GZ.
            do k=k1,k2
              s=h(k)*t
              g=w(k)*t
              xas=s-xa
              yas=s-ya
              zas=s-za
              xbs=s-xb
              ybs=s-yb
              zbs=s-zb
              ma=0
              pax=g
              pay=g
              paz=g
              do ia=0,la
                mb=ma
                px=pax
                py=pay
                pz=paz
                do ib=0,lb
                  pq=1.0d0
                  do iq=1,nq
                    gx(mb+iq)=gx(mb+iq)+px*pq
                    gy(mb+iq)=gy(mb+iq)+py*pq
                    gz(mb+iq)=gz(mb+iq)+pz*pq
                    pq=pq*s
                  end do
                  mb=mb+ll1
                  px=px*xbs
                  py=py*ybs
                  pz=pz*zbs
                end do
                ma=ma+l1*ll1
                pax=pax*xas
                pay=pay*yas
                paz=paz*zas
              end do
            end do
!  Now these basic integrals are used to construct the multipole moments
!  for the overlap density corresponding to each pair of basis functions
!  in the pair of shells.
            qt=0.0d0
            ci=cs(ig)
            do ia=mini,maxi
              if (ia .gt. 1 .and. la .eq. 1) ci=cp(ig)
              ! if (ia .gt. 4) ci=cd(ig)
              cj=cs(jg)
              do jb=minj,maxj
                if (jb .gt. 1 .and. lb .eq. 1) cj=cp(jg)
                ! if (jb .gt. 4) cj=cd(jg)
                f=-fac*ci*cj*d(ia,jb)
                mx=(ix(ia)*l1+ix(jb))*ll1+1
                my=(iy(ia)*l1+iy(jb))*ll1+1
                mz=(iz(ia)*l1+iz(jb))*ll1+1
!  Now the integral of (x**i)*(y**j)*(z**k) over the current pair of
!  atomic primitive orbitals is F*GX(MX+i)*GY(MY+j)*GZ(MZ+k).
                call addqlm(lq, f, gx(mx:),gy(my:),gz(mz:))
!  End of loop over basis functions
              end do
            end do

            if (iw .gt. 0) write (iw,1003) qt(1:nq**2)
1003        format (f10.6: / 3f10.6: / 5f10.6: / 7f10.6: / 9f10.6: /   &
                11f10.6: / 13f10.6: / 15f10.6: / 17f10.6: / 19f10.6: / &
                21f10.6)
!  Move multipoles to expansion centre nearest to overlap centre P.
            call moveq(xp,yp,zp)

            else
!  General case. Use integration over grid to assign multipole moments
!  to atoms.
              if (.not. do_quadrature) then
                !  Clear grid of density values
                rho(:)=0d0
                do_quadrature=.true.
              end if

!  P is the overlap centre
              p=aj/aa
              xp=xi-p*xji
              yp=yi-p*yji
              zp=zi-p*zji
              if (iw > 0) write(iw, "(2i5, 3f20.15)") ig, jg, xp, yp, zp
              do k=1,ng
                xk=grid(1,k)
                yk=grid(2,k)
                zk=grid(3,k)
                xa=xk-xi
                ya=yk-yi
                za=zk-zi
                xb=xk-xj
                yb=yk-yj
                zb=zk-zj
                e=aa*((xk-xp)**2+(yk-yp)**2+(zk-zp)**2)
                if (e > etol) cycle
                e=exp(-e)
                ci=cs(ig)
                do l=1,la
                  xai(l)=xa*xai(l-1)
                  yai(l)=ya*yai(l-1)
                  zai(l)=za*zai(l-1)
                end do
                do ia=mini,maxi
                  if (ia .gt. 1 .and. la .eq. 1) ci=cp(ig)
                  ! if (ia .gt. 4) ci=cd(ig)
                  cj=cs(jg)
                  do l=1,lb
                    xbj(l)=xb*xbj(l-1)
                    ybj(l)=yb*ybj(l-1)
                    zbj(l)=zb*zbj(l-1)
                  end do
                  do jb=minj,maxj
                    if (jb .gt. 1 .and. lb .eq. 1) cj=cp(jg)
                    ! if (jb .gt. 4) cj=cd(jg)
                    f=-fac*ci*cj*d(ia,jb)
                    ! print "(4i5)", k, ia, jb
                    rho(k)=rho(k)+e*f*                                 &
                        xai(ix(ia))*yai(iy(ia))*zai(iz(ia))*           &
                        xbj(ix(jb))*ybj(iy(jb))*zbj(iz(jb))
                  end do
                end do
              end do

!  Now we have added in the overlap charge density for this pair of
!  primitives

            end if
            !  End of loop over primitives
          end do
        end do

        !  End of loop over shells for atom j
      end do
      !  End of loop over atoms j
    end do
    !  End of loop over shells for atom i
  end do
  !  End of loop over atoms i
end do

if (do_quadrature) then
  !  Loop over sites to evaluate multipole contributions from charge density
  gx(0)=1d0
  gy(0)=1d0
  gz(0)=1d0
  do i=1,ns
    xi=xs(1,i)
    yi=xs(2,i)
    zi=xs(3,i)
    qt=0d0
    !  Work backwards (most distant points first) to avoid loss of
    !  significance.
    do k=start(i+1)-1,start(i),-1
      ! print "(i6,f12.7)", k, rho(k)
      xk=grid(1,k)-xi
      yk=grid(2,k)-yi
      zk=grid(3,k)-zi
      do l=1,lmax
        gx(l)=gx(l-1)*xk
        gy(l)=gy(l-1)*yk
        gz(l)=gz(l-1)*zk
      end do
      call addqlm(lmax, grid(4,k)*rho(k), gx, gy, gz)
    end do
    call moveq(xi,yi,zi)
  end do
end if

if (allocated(grid)) deallocate(grid)
if (allocated(start)) deallocate(start)

END SUBROUTINE dmaqlm

!-----------------------------------------------------------------ADDQLM

SUBROUTINE addqlm(l, f, gx,gy,gz)
IMPLICIT NONE
REAL(dp), INTENT(IN) :: gx(0:), gy(0:), gz(0:), f
INTEGER, INTENT(IN) :: l
!-----------------------------------------------------
!     Copyright A J Stone University of Cambridge 1983
!     Modifications and interface to Cadpac, R D Amos
!     Version for Cadpac5 , R D Amos, June 1990
!-----------------------------------------------------

REAL(dp) :: Q00, Q10,Q11c,Q11s, Q20,Q21c,Q21s,Q22c,Q22s,               &
    Q30,Q31c,Q31s,Q32c,Q32s,Q33c,Q33s,                                 &
    Q40,Q41c,Q41s,Q42c,Q42s,Q43c,Q43s,Q44c,Q44s,                       &
    Q50,Q51c,Q51s,Q52c,Q52s,Q53c,Q53s,Q54c,Q54s,Q55c,Q55s,             &
    Q60,Q61c,Q61s,Q62c,Q62s,Q63c,Q63s,Q64c,Q64s,Q65c,Q65s,Q66c,Q66s,   &
    Q70,Q71c,Q71s,Q72c,Q72s,Q73c,Q73s,Q74c,Q74s,Q75c,Q75s,Q76c,Q76s,   &
    Q77c,Q77s,                                                         &
    Q80,Q81c,Q81s,Q82c,Q82s,Q83c,Q83s,Q84c,Q84s,Q85c,Q85s,Q86c,Q86s,   &
    Q87c,Q87s,Q88c,Q88s,                                               &
    Q90,Q91c,Q91s,Q92c,Q92s,Q93c,Q93s,Q94c,Q94s,Q95c,Q95s,Q96c,Q96s,   &
    Q97c,Q97s,Q98c,Q98s,Q99c,Q99s,                                     &
    QA0,QA1c,QA1s,QA2c,QA2s,QA3c,QA3s,QA4c,QA4s,QA5c,QA5s,QA6c,QA6s,   &
    QA7c,QA7s,QA8c,QA8s,QA9c,QA9s,QAAc,QAAs

COMMON/BIG/Q00, Q10,Q11c,Q11s, Q20,Q21c,Q21s,Q22c,Q22s,                &
    Q30,Q31c,Q31s,Q32c,Q32s,Q33c,Q33s,                                 &
    Q40,Q41c,Q41s,Q42c,Q42s,Q43c,Q43s,Q44c,Q44s,                       &
    Q50,Q51c,Q51s,Q52c,Q52s,Q53c,Q53s,Q54c,Q54s,Q55c,Q55s,             &
    Q60,Q61c,Q61s,Q62c,Q62s,Q63c,Q63s,Q64c,Q64s,Q65c,Q65s,Q66c,Q66s,   &
    Q70,Q71c,Q71s,Q72c,Q72s,Q73c,Q73s,Q74c,Q74s,Q75c,Q75s,Q76c,Q76s,   &
    Q77c,Q77s,                                                         &
    Q80,Q81c,Q81s,Q82c,Q82s,Q83c,Q83s,Q84c,Q84s,Q85c,Q85s,Q86c,Q86s,   &
    Q87c,Q87s,Q88c,Q88s,                                               &
    Q90,Q91c,Q91s,Q92c,Q92s,Q93c,Q93s,Q94c,Q94s,Q95c,Q95s,Q96c,Q96s,   &
    Q97c,Q97s,Q98c,Q98s,Q99c,Q99s,                                     &
    QA0,QA1c,QA1s,QA2c,QA2s,QA3c,QA3s,QA4c,QA4s,QA5c,QA5s,QA6c,QA6s,   &
    QA7c,QA7s,QA8c,QA8s,QA9c,QA9s,QAAc,QAAs
    

go to (100,101,102,103,104,105,106,107,108,109,110), min(l,lmax)+1
110                                                                     &
QA0=QA0+f*(256d0*gx(0)*gy(0)*gz(10)-5760d0*gx(0)*gy(2)*gz(8)            &
    -5760d0*gx(2)*gy(0)*gz(8)+20160d0*gx(0)*gy(4)*gz(6)                 &
    +40320d0*gx(2)*gy(2)*gz(6)+20160d0*gx(4)*gy(0)*gz(6)                &
    -16800d0*gx(0)*gy(6)*gz(4)-50400d0*gx(2)*gy(4)*gz(4)                &
    -50400d0*gx(4)*gy(2)*gz(4)-16800d0*gx(6)*gy(0)*gz(4)                &
    +3150d0*gx(0)*gy(8)*gz(2)+12600d0*gx(2)*gy(6)*gz(2)                 &
    +18900d0*gx(4)*gy(4)*gz(2)+12600d0*gx(6)*gy(2)*gz(2)                &
    +3150d0*gx(8)*gy(0)*gz(2)-63d0*gx(0)*gy(10)*gz(0)                   &
    -315d0*gx(2)*gy(8)*gz(0)-630d0*gx(4)*gy(6)*gz(0)                    &
    -630d0*gx(6)*gy(4)*gz(0)-315d0*gx(8)*gy(2)*gz(0)                    &
    -63d0*gx(10)*gy(0)*gz(0))/256d0
QA1s=QA1s+f*rt(5)*rt(11)*(128d0*gx(0)*gy(1)*gz(9)                       &
    -1152d0*gx(0)*gy(3)*gz(7)-1152d0*gx(2)*gy(1)*gz(7)                  &
    +2016d0*gx(0)*gy(5)*gz(5)+4032d0*gx(2)*gy(3)*gz(5)                  &
    +2016d0*gx(4)*gy(1)*gz(5)-840d0*gx(0)*gy(7)*gz(3)                   &
    -2520d0*gx(2)*gy(5)*gz(3)-2520d0*gx(4)*gy(3)*gz(3)                  &
    -840d0*gx(6)*gy(1)*gz(3)+63d0*gx(0)*gy(9)*gz(1)                     &
    +252d0*gx(2)*gy(7)*gz(1)+378d0*gx(4)*gy(5)*gz(1)                    &
    +252d0*gx(6)*gy(3)*gz(1)+63d0*gx(8)*gy(1)*gz(1))/128d0
QA1c=QA1c+f*rt(5)*rt(11)*(128d0*gx(1)*gy(0)*gz(9)                       &
    -1152d0*gx(1)*gy(2)*gz(7)-1152d0*gx(3)*gy(0)*gz(7)                  &
    +2016d0*gx(1)*gy(4)*gz(5)+4032d0*gx(3)*gy(2)*gz(5)                  &
    +2016d0*gx(5)*gy(0)*gz(5)-840d0*gx(1)*gy(6)*gz(3)                   &
    -2520d0*gx(3)*gy(4)*gz(3)-2520d0*gx(5)*gy(2)*gz(3)                  &
    -840d0*gx(7)*gy(0)*gz(3)+63d0*gx(1)*gy(8)*gz(1)                     &
    +252d0*gx(3)*gy(6)*gz(1)+378d0*gx(5)*gy(4)*gz(1)                    &
    +252d0*gx(7)*gy(2)*gz(1)+63d0*gx(9)*gy(0)*gz(1))/128d0
QA2s=QA2s+f*rt(3)*rt(5)*rt(11)*(768d0*gx(1)*gy(1)*gz(8)                 &
    -3584d0*gx(1)*gy(3)*gz(6)-3584d0*gx(3)*gy(1)*gz(6)                  &
    +3360d0*gx(1)*gy(5)*gz(4)+6720d0*gx(3)*gy(3)*gz(4)                  &
    +3360d0*gx(5)*gy(1)*gz(4)-672d0*gx(1)*gy(7)*gz(2)                   &
    -2016d0*gx(3)*gy(5)*gz(2)-2016d0*gx(5)*gy(3)*gz(2)                  &
    -672d0*gx(7)*gy(1)*gz(2)+14d0*gx(1)*gy(9)*gz(0)                     &
    +56d0*gx(3)*gy(7)*gz(0)+84d0*gx(5)*gy(5)*gz(0)                      &
    +56d0*gx(7)*gy(3)*gz(0)+14d0*gx(9)*gy(1)*gz(0))/256d0
QA2c=QA2c+f*rt(3)*rt(5)*rt(11)*(-384d0*gx(0)*gy(2)*gz(8)                &
    +384d0*gx(2)*gy(0)*gz(8)+1792d0*gx(0)*gy(4)*gz(6)                   &
    -1792d0*gx(4)*gy(0)*gz(6)-1680d0*gx(0)*gy(6)*gz(4)                  &
    -1680d0*gx(2)*gy(4)*gz(4)+1680d0*gx(4)*gy(2)*gz(4)                  &
    +1680d0*gx(6)*gy(0)*gz(4)+336d0*gx(0)*gy(8)*gz(2)                   &
    +672d0*gx(2)*gy(6)*gz(2)-672d0*gx(6)*gy(2)*gz(2)                    &
    -336d0*gx(8)*gy(0)*gz(2)-7d0*gx(0)*gy(10)*gz(0)                     &
    -21d0*gx(2)*gy(8)*gz(0)-14d0*gx(4)*gy(6)*gz(0)                      &
    +14d0*gx(6)*gy(4)*gz(0)+21d0*gx(8)*gy(2)*gz(0)                      &
    +7d0*gx(10)*gy(0)*gz(0))/256d0
QA3s=QA3s+f*rt(2)*rt(3)*rt(5)*rt(11)*rt(13)*(-64d0*gx(0)*gy(3)*gz(7)    &
    +192d0*gx(2)*gy(1)*gz(7)+168d0*gx(0)*gy(5)*gz(5)                    &
    -336d0*gx(2)*gy(3)*gz(5)-504d0*gx(4)*gy(1)*gz(5)                    &
    -84d0*gx(0)*gy(7)*gz(3)+84d0*gx(2)*gy(5)*gz(3)                      &
    +420d0*gx(4)*gy(3)*gz(3)+252d0*gx(6)*gy(1)*gz(3)                    &
    +7d0*gx(0)*gy(9)*gz(1)-42d0*gx(4)*gy(5)*gz(1)                       &
    -56d0*gx(6)*gy(3)*gz(1)-21d0*gx(8)*gy(1)*gz(1))/128d0
QA3c=QA3c+f*rt(2)*rt(3)*rt(5)*rt(11)*rt(13)*(-192d0*gx(1)*gy(2)*gz(7)   &
    +64d0*gx(3)*gy(0)*gz(7)+504d0*gx(1)*gy(4)*gz(5)                     &
    +336d0*gx(3)*gy(2)*gz(5)-168d0*gx(5)*gy(0)*gz(5)                    &
    -252d0*gx(1)*gy(6)*gz(3)-420d0*gx(3)*gy(4)*gz(3)                    &
    -84d0*gx(5)*gy(2)*gz(3)+84d0*gx(7)*gy(0)*gz(3)                      &
    +21d0*gx(1)*gy(8)*gz(1)+56d0*gx(3)*gy(6)*gz(1)                      &
    +42d0*gx(5)*gy(4)*gz(1)-7d0*gx(9)*gy(0)*gz(1))/128d0
QA4s=QA4s+f*rt(3)*rt(5)*rt(11)*rt(13)*(-448d0*gx(1)*gy(3)*gz(6)         &
    +448d0*gx(3)*gy(1)*gz(6)+672d0*gx(1)*gy(5)*gz(4)                    &
    -672d0*gx(5)*gy(1)*gz(4)-168d0*gx(1)*gy(7)*gz(2)                    &
    -168d0*gx(3)*gy(5)*gz(2)+168d0*gx(5)*gy(3)*gz(2)                    &
    +168d0*gx(7)*gy(1)*gz(2)+4d0*gx(1)*gy(9)*gz(0)                      &
    +8d0*gx(3)*gy(7)*gz(0)-8d0*gx(7)*gy(3)*gz(0)-4d0*gx(9)*gy(1)*gz(0)) &
    /128d0
QA4c=QA4c+f*rt(3)*rt(5)*rt(11)*rt(13)*(112d0*gx(0)*gy(4)*gz(6)          &
    -672d0*gx(2)*gy(2)*gz(6)+112d0*gx(4)*gy(0)*gz(6)                    &
    -168d0*gx(0)*gy(6)*gz(4)+840d0*gx(2)*gy(4)*gz(4)                    &
    +840d0*gx(4)*gy(2)*gz(4)-168d0*gx(6)*gy(0)*gz(4)                    &
    +42d0*gx(0)*gy(8)*gz(2)-168d0*gx(2)*gy(6)*gz(2)                     &
    -420d0*gx(4)*gy(4)*gz(2)-168d0*gx(6)*gy(2)*gz(2)                    &
    +42d0*gx(8)*gy(0)*gz(2)-gx(0)*gy(10)*gz(0)+3d0*gx(2)*gy(8)*gz(0)    &
    +14d0*gx(4)*gy(6)*gz(0)+14d0*gx(6)*gy(4)*gz(0)                      &
    +3d0*gx(8)*gy(2)*gz(0)-gx(10)*gy(0)*gz(0))/128d0
QA5s=QA5s+f*rt(2)*rt(3)*rt(11)*rt(13)*(168d0*gx(0)*gy(5)*gz(5)          &
    -1680d0*gx(2)*gy(3)*gz(5)+840d0*gx(4)*gy(1)*gz(5)                   &
    -140d0*gx(0)*gy(7)*gz(3)+1260d0*gx(2)*gy(5)*gz(3)                   &
    +700d0*gx(4)*gy(3)*gz(3)-700d0*gx(6)*gy(1)*gz(3)                    &
    +15d0*gx(0)*gy(9)*gz(1)-120d0*gx(2)*gy(7)*gz(1)                     &
    -210d0*gx(4)*gy(5)*gz(1)+75d0*gx(8)*gy(1)*gz(1))/128d0
QA5c=QA5c+f*rt(2)*rt(3)*rt(11)*rt(13)*(840d0*gx(1)*gy(4)*gz(5)          &
    -1680d0*gx(3)*gy(2)*gz(5)+168d0*gx(5)*gy(0)*gz(5)                   &
    -700d0*gx(1)*gy(6)*gz(3)+700d0*gx(3)*gy(4)*gz(3)                    &
    +1260d0*gx(5)*gy(2)*gz(3)-140d0*gx(7)*gy(0)*gz(3)                   &
    +75d0*gx(1)*gy(8)*gz(1)-210d0*gx(5)*gy(4)*gz(1)                     &
    -120d0*gx(7)*gy(2)*gz(1)+15d0*gx(9)*gy(0)*gz(1))/128d0
QA6s=QA6s+f*rt(2)*rt(3)*rt(5)*rt(11)*rt(13)*(1344d0*gx(1)*gy(5)*gz(4)   &
    -4480d0*gx(3)*gy(3)*gz(4)+1344d0*gx(5)*gy(1)*gz(4)                  &
    -576d0*gx(1)*gy(7)*gz(2)+1344d0*gx(3)*gy(5)*gz(2)                   &
    +1344d0*gx(5)*gy(3)*gz(2)-576d0*gx(7)*gy(1)*gz(2)                   &
    +18d0*gx(1)*gy(9)*gz(0)-24d0*gx(3)*gy(7)*gz(0)                      &
    -84d0*gx(5)*gy(5)*gz(0)-24d0*gx(7)*gy(3)*gz(0)                      &
    +18d0*gx(9)*gy(1)*gz(0))/512d0
QA6c=QA6c+f*rt(2)*rt(3)*rt(5)*rt(11)*rt(13)*(-224d0*gx(0)*gy(6)*gz(4)   &
    +3360d0*gx(2)*gy(4)*gz(4)-3360d0*gx(4)*gy(2)*gz(4)                  &
    +224d0*gx(6)*gy(0)*gz(4)+96d0*gx(0)*gy(8)*gz(2)                     &
    -1344d0*gx(2)*gy(6)*gz(2)+1344d0*gx(6)*gy(2)*gz(2)                  &
    -96d0*gx(8)*gy(0)*gz(2)-3d0*gx(0)*gy(10)*gz(0)                      &
    +39d0*gx(2)*gy(8)*gz(0)+42d0*gx(4)*gy(6)*gz(0)                      &
    -42d0*gx(6)*gy(4)*gz(0)-39d0*gx(8)*gy(2)*gz(0)                      &
    +3d0*gx(10)*gy(0)*gz(0))/512d0
QA7s=QA7s+f*rt(2)*rt(3)*rt(5)*rt(11)*rt(13)*rt(17)*(                    &
    -16d0*gx(0)*gy(7)*gz(3)+336d0*gx(2)*gy(5)*gz(3)                     &
    -560d0*gx(4)*gy(3)*gz(3)+112d0*gx(6)*gy(1)*gz(3)                    &
    +3d0*gx(0)*gy(9)*gz(1)-60d0*gx(2)*gy(7)*gz(1)                       &
    +42d0*gx(4)*gy(5)*gz(1)+84d0*gx(6)*gy(3)*gz(1)                      &
    -21d0*gx(8)*gy(1)*gz(1))/256d0
QA7c=QA7c+f*rt(2)*rt(3)*rt(5)*rt(11)*rt(13)*rt(17)*(                    &
    -112d0*gx(1)*gy(6)*gz(3)+560d0*gx(3)*gy(4)*gz(3)                    &
    -336d0*gx(5)*gy(2)*gz(3)+16d0*gx(7)*gy(0)*gz(3)                     &
    +21d0*gx(1)*gy(8)*gz(1)-84d0*gx(3)*gy(6)*gz(1)                      &
    -42d0*gx(5)*gy(4)*gz(1)+60d0*gx(7)*gy(2)*gz(1)                      &
    -3d0*gx(9)*gy(0)*gz(1))/256d0
QA8s=QA8s+f*rt(5)*rt(11)*rt(13)*rt(17)*(-144d0*gx(1)*gy(7)*gz(2)        &
    +1008d0*gx(3)*gy(5)*gz(2)-1008d0*gx(5)*gy(3)*gz(2)                  &
    +144d0*gx(7)*gy(1)*gz(2)+8d0*gx(1)*gy(9)*gz(0)                      &
    -48d0*gx(3)*gy(7)*gz(0)+48d0*gx(7)*gy(3)*gz(0)                      &
    -8d0*gx(9)*gy(1)*gz(0))/256d0
QA8c=QA8c+f*rt(5)*rt(11)*rt(13)*rt(17)*(18d0*gx(0)*gy(8)*gz(2)          &
    -504d0*gx(2)*gy(6)*gz(2)+1260d0*gx(4)*gy(4)*gz(2)                   &
    -504d0*gx(6)*gy(2)*gz(2)+18d0*gx(8)*gy(0)*gz(2)-gx(0)*gy(10)*gz(0)  &
    +27d0*gx(2)*gy(8)*gz(0)-42d0*gx(4)*gy(6)*gz(0)                      &
    -42d0*gx(6)*gy(4)*gz(0)+27d0*gx(8)*gy(2)*gz(0)-gx(10)*gy(0)*gz(0))  &
    /256d0
QA9s=QA9s+f*rt(2)*rt(5)*rt(11)*rt(13)*rt(17)*rt(19)*(gx(0)*gy(9)*gz(1)  &
    -36d0*gx(2)*gy(7)*gz(1)+126d0*gx(4)*gy(5)*gz(1)                     &
    -84d0*gx(6)*gy(3)*gz(1)+9d0*gx(8)*gy(1)*gz(1))/256d0
QA9c=QA9c+f*rt(2)*rt(5)*rt(11)*rt(13)*rt(17)*rt(19)*(                   &
    9d0*gx(1)*gy(8)*gz(1)-84d0*gx(3)*gy(6)*gz(1)                        &
    +126d0*gx(5)*gy(4)*gz(1)-36d0*gx(7)*gy(2)*gz(1)+gx(9)*gy(0)*gz(1))  &
    /256d0
QAAs=QAAs+f*rt(2)*rt(11)*rt(13)*rt(17)*rt(19)*(10d0*gx(1)*gy(9)*gz(0)   &
    -120d0*gx(3)*gy(7)*gz(0)+252d0*gx(5)*gy(5)*gz(0)                    &
    -120d0*gx(7)*gy(3)*gz(0)+10d0*gx(9)*gy(1)*gz(0))/512d0
QAAc=QAAc+f*rt(2)*rt(11)*rt(13)*rt(17)*rt(19)*(-gx(0)*gy(10)*gz(0)      &
    +45d0*gx(2)*gy(8)*gz(0)-210d0*gx(4)*gy(6)*gz(0)                     &
    +210d0*gx(6)*gy(4)*gz(0)-45d0*gx(8)*gy(2)*gz(0)+gx(10)*gy(0)*gz(0)) &
    /512d0

109                                                                     &
Q90=Q90+f*(128d0*gx(0)*gy(0)*gz(9)-2304d0*gx(0)*gy(2)*gz(7)             &
    -2304d0*gx(2)*gy(0)*gz(7)+6048d0*gx(0)*gy(4)*gz(5)                  &
    +12096d0*gx(2)*gy(2)*gz(5)+6048d0*gx(4)*gy(0)*gz(5)                 &
    -3360d0*gx(0)*gy(6)*gz(3)-10080d0*gx(2)*gy(4)*gz(3)                 &
    -10080d0*gx(4)*gy(2)*gz(3)-3360d0*gx(6)*gy(0)*gz(3)                 &
    +315d0*gx(0)*gy(8)*gz(1)+1260d0*gx(2)*gy(6)*gz(1)                   &
    +1890d0*gx(4)*gy(4)*gz(1)+1260d0*gx(6)*gy(2)*gz(1)                  &
    +315d0*gx(8)*gy(0)*gz(1))/128d0
Q91s=Q91s+f*rt(5)*(384d0*gx(0)*gy(1)*gz(8)-2688d0*gx(0)*gy(3)*gz(6)     &
    -2688d0*gx(2)*gy(1)*gz(6)+3360d0*gx(0)*gy(5)*gz(4)                  &
    +6720d0*gx(2)*gy(3)*gz(4)+3360d0*gx(4)*gy(1)*gz(4)                  &
    -840d0*gx(0)*gy(7)*gz(2)-2520d0*gx(2)*gy(5)*gz(2)                   &
    -2520d0*gx(4)*gy(3)*gz(2)-840d0*gx(6)*gy(1)*gz(2)                   &
    +21d0*gx(0)*gy(9)*gz(0)+84d0*gx(2)*gy(7)*gz(0)                      &
    +126d0*gx(4)*gy(5)*gz(0)+84d0*gx(6)*gy(3)*gz(0)                     &
    +21d0*gx(8)*gy(1)*gz(0))/128d0
Q91c=Q91c+f*rt(5)*(384d0*gx(1)*gy(0)*gz(8)-2688d0*gx(1)*gy(2)*gz(6)     &
    -2688d0*gx(3)*gy(0)*gz(6)+3360d0*gx(1)*gy(4)*gz(4)                  &
    +6720d0*gx(3)*gy(2)*gz(4)+3360d0*gx(5)*gy(0)*gz(4)                  &
    -840d0*gx(1)*gy(6)*gz(2)-2520d0*gx(3)*gy(4)*gz(2)                   &
    -2520d0*gx(5)*gy(2)*gz(2)-840d0*gx(7)*gy(0)*gz(2)                   &
    +21d0*gx(1)*gy(8)*gz(0)+84d0*gx(3)*gy(6)*gz(0)                      &
    +126d0*gx(5)*gy(4)*gz(0)+84d0*gx(7)*gy(2)*gz(0)                     &
    +21d0*gx(9)*gy(0)*gz(0))/128d0
Q92s=Q92s+f*rt(2)*rt(5)*rt(11)*(192d0*gx(1)*gy(1)*gz(7)                 &
    -672d0*gx(1)*gy(3)*gz(5)-672d0*gx(3)*gy(1)*gz(5)                    &
    +420d0*gx(1)*gy(5)*gz(3)+840d0*gx(3)*gy(3)*gz(3)                    &
    +420d0*gx(5)*gy(1)*gz(3)-42d0*gx(1)*gy(7)*gz(1)                     &
    -126d0*gx(3)*gy(5)*gz(1)-126d0*gx(5)*gy(3)*gz(1)                    &
    -42d0*gx(7)*gy(1)*gz(1))/64d0
Q92c=Q92c+f*rt(2)*rt(5)*rt(11)*(-96d0*gx(0)*gy(2)*gz(7)                 &
    +96d0*gx(2)*gy(0)*gz(7)+336d0*gx(0)*gy(4)*gz(5)                     &
    -336d0*gx(4)*gy(0)*gz(5)-210d0*gx(0)*gy(6)*gz(3)                    &
    -210d0*gx(2)*gy(4)*gz(3)+210d0*gx(4)*gy(2)*gz(3)                    &
    +210d0*gx(6)*gy(0)*gz(3)+21d0*gx(0)*gy(8)*gz(1)                     &
    +42d0*gx(2)*gy(6)*gz(1)-42d0*gx(6)*gy(2)*gz(1)                      &
    -21d0*gx(8)*gy(0)*gz(1))/64d0
Q93s=Q93s+f*rt(2)*rt(3)*rt(5)*rt(7)*rt(11)*(-64d0*gx(0)*gy(3)*gz(6)     &
    +192d0*gx(2)*gy(1)*gz(6)+120d0*gx(0)*gy(5)*gz(4)                    &
    -240d0*gx(2)*gy(3)*gz(4)-360d0*gx(4)*gy(1)*gz(4)                    &
    -36d0*gx(0)*gy(7)*gz(2)+36d0*gx(2)*gy(5)*gz(2)                      &
    +180d0*gx(4)*gy(3)*gz(2)+108d0*gx(6)*gy(1)*gz(2)+gx(0)*gy(9)*gz(0)  &
    -6d0*gx(4)*gy(5)*gz(0)-8d0*gx(6)*gy(3)*gz(0)-3d0*gx(8)*gy(1)*gz(0)) &
    /128d0
Q93c=Q93c+f*rt(2)*rt(3)*rt(5)*rt(7)*rt(11)*(-192d0*gx(1)*gy(2)*gz(6)    &
    +64d0*gx(3)*gy(0)*gz(6)+360d0*gx(1)*gy(4)*gz(4)                     &
    +240d0*gx(3)*gy(2)*gz(4)-120d0*gx(5)*gy(0)*gz(4)                    &
    -108d0*gx(1)*gy(6)*gz(2)-180d0*gx(3)*gy(4)*gz(2)                    &
    -36d0*gx(5)*gy(2)*gz(2)+36d0*gx(7)*gy(0)*gz(2)                      &
    +3d0*gx(1)*gy(8)*gz(0)+8d0*gx(3)*gy(6)*gz(0)+6d0*gx(5)*gy(4)*gz(0)  &
    -gx(9)*gy(0)*gz(0))/128d0
Q94s=Q94s+f*rt(5)*rt(7)*rt(11)*rt(13)*(-96d0*gx(1)*gy(3)*gz(5)          &
    +96d0*gx(3)*gy(1)*gz(5)+96d0*gx(1)*gy(5)*gz(3)                      &
    -96d0*gx(5)*gy(1)*gz(3)-12d0*gx(1)*gy(7)*gz(1)                      &
    -12d0*gx(3)*gy(5)*gz(1)+12d0*gx(5)*gy(3)*gz(1)                      &
    +12d0*gx(7)*gy(1)*gz(1))/64d0
Q94c=Q94c+f*rt(5)*rt(7)*rt(11)*rt(13)*(24d0*gx(0)*gy(4)*gz(5)           &
    -144d0*gx(2)*gy(2)*gz(5)+24d0*gx(4)*gy(0)*gz(5)                     &
    -24d0*gx(0)*gy(6)*gz(3)+120d0*gx(2)*gy(4)*gz(3)                     &
    +120d0*gx(4)*gy(2)*gz(3)-24d0*gx(6)*gy(0)*gz(3)                     &
    +3d0*gx(0)*gy(8)*gz(1)-12d0*gx(2)*gy(6)*gz(1)                       &
    -30d0*gx(4)*gy(4)*gz(1)-12d0*gx(6)*gy(2)*gz(1)                      &
    +3d0*gx(8)*gy(0)*gz(1))/64d0
Q95s=Q95s+f*rt(2)*rt(11)*rt(13)*(168d0*gx(0)*gy(5)*gz(4)                &
    -1680d0*gx(2)*gy(3)*gz(4)+840d0*gx(4)*gy(1)*gz(4)                   &
    -84d0*gx(0)*gy(7)*gz(2)+756d0*gx(2)*gy(5)*gz(2)                     &
    +420d0*gx(4)*gy(3)*gz(2)-420d0*gx(6)*gy(1)*gz(2)                    &
    +3d0*gx(0)*gy(9)*gz(0)-24d0*gx(2)*gy(7)*gz(0)                       &
    -42d0*gx(4)*gy(5)*gz(0)+15d0*gx(8)*gy(1)*gz(0))/128d0
Q95c=Q95c+f*rt(2)*rt(11)*rt(13)*(840d0*gx(1)*gy(4)*gz(4)                &
    -1680d0*gx(3)*gy(2)*gz(4)+168d0*gx(5)*gy(0)*gz(4)                   &
    -420d0*gx(1)*gy(6)*gz(2)+420d0*gx(3)*gy(4)*gz(2)                    &
    +756d0*gx(5)*gy(2)*gz(2)-84d0*gx(7)*gy(0)*gz(2)                     &
    +15d0*gx(1)*gy(8)*gz(0)-42d0*gx(5)*gy(4)*gz(0)                      &
    -24d0*gx(7)*gy(2)*gz(0)+3d0*gx(9)*gy(0)*gz(0))/128d0
Q96s=Q96s+f*rt(2)*rt(3)*rt(5)*rt(11)*rt(13)*(84d0*gx(1)*gy(5)*gz(3)     &
    -280d0*gx(3)*gy(3)*gz(3)+84d0*gx(5)*gy(1)*gz(3)                     &
    -18d0*gx(1)*gy(7)*gz(1)+42d0*gx(3)*gy(5)*gz(1)                      &
    +42d0*gx(5)*gy(3)*gz(1)-18d0*gx(7)*gy(1)*gz(1))/64d0
Q96c=Q96c+f*rt(2)*rt(3)*rt(5)*rt(11)*rt(13)*(-14d0*gx(0)*gy(6)*gz(3)    &
    +210d0*gx(2)*gy(4)*gz(3)-210d0*gx(4)*gy(2)*gz(3)                    &
    +14d0*gx(6)*gy(0)*gz(3)+3d0*gx(0)*gy(8)*gz(1)                       &
    -42d0*gx(2)*gy(6)*gz(1)+42d0*gx(6)*gy(2)*gz(1)                      &
    -3d0*gx(8)*gy(0)*gz(1))/64d0
Q97s=Q97s+f*rt(2)*rt(5)*rt(11)*rt(13)*(-48d0*gx(0)*gy(7)*gz(2)          &
    +1008d0*gx(2)*gy(5)*gz(2)-1680d0*gx(4)*gy(3)*gz(2)                  &
    +336d0*gx(6)*gy(1)*gz(2)+3d0*gx(0)*gy(9)*gz(0)                      &
    -60d0*gx(2)*gy(7)*gz(0)+42d0*gx(4)*gy(5)*gz(0)                      &
    +84d0*gx(6)*gy(3)*gz(0)-21d0*gx(8)*gy(1)*gz(0))/256d0
Q97c=Q97c+f*rt(2)*rt(5)*rt(11)*rt(13)*(-336d0*gx(1)*gy(6)*gz(2)         &
    +1680d0*gx(3)*gy(4)*gz(2)-1008d0*gx(5)*gy(2)*gz(2)                  &
    +48d0*gx(7)*gy(0)*gz(2)+21d0*gx(1)*gy(8)*gz(0)                      &
    -84d0*gx(3)*gy(6)*gz(0)-42d0*gx(5)*gy(4)*gz(0)                      &
    +60d0*gx(7)*gy(2)*gz(0)-3d0*gx(9)*gy(0)*gz(0))/256d0
Q98s=Q98s+f*rt(5)*rt(11)*rt(13)*rt(17)*(-24d0*gx(1)*gy(7)*gz(1)         &
    +168d0*gx(3)*gy(5)*gz(1)-168d0*gx(5)*gy(3)*gz(1)                    &
    +24d0*gx(7)*gy(1)*gz(1))/128d0
Q98c=Q98c+f*rt(5)*rt(11)*rt(13)*rt(17)*(3d0*gx(0)*gy(8)*gz(1)           &
    -84d0*gx(2)*gy(6)*gz(1)+210d0*gx(4)*gy(4)*gz(1)                     &
    -84d0*gx(6)*gy(2)*gz(1)+3d0*gx(8)*gy(0)*gz(1))/128d0
Q99s=Q99s+f*rt(2)*rt(5)*rt(11)*rt(13)*rt(17)*(gx(0)*gy(9)*gz(0)         &
    -36d0*gx(2)*gy(7)*gz(0)+126d0*gx(4)*gy(5)*gz(0)                     &
    -84d0*gx(6)*gy(3)*gz(0)+9d0*gx(8)*gy(1)*gz(0))/256d0
Q99c=Q99c+f*rt(2)*rt(5)*rt(11)*rt(13)*rt(17)*(9d0*gx(1)*gy(8)*gz(0)     &
    -84d0*gx(3)*gy(6)*gz(0)+126d0*gx(5)*gy(4)*gz(0)                     &
    -36d0*gx(7)*gy(2)*gz(0)+gx(9)*gy(0)*gz(0))/256d0

108                                                                     &
Q80=Q80+f*(128d0*gx(0)*gy(0)*gz(8)-1792d0*gx(0)*gy(2)*gz(6)             &
    -1792d0*gx(2)*gy(0)*gz(6)+3360d0*gx(0)*gy(4)*gz(4)                  &
    +6720d0*gx(2)*gy(2)*gz(4)+3360d0*gx(4)*gy(0)*gz(4)                  &
    -1120d0*gx(0)*gy(6)*gz(2)-3360d0*gx(2)*gy(4)*gz(2)                  &
    -3360d0*gx(4)*gy(2)*gz(2)-1120d0*gx(6)*gy(0)*gz(2)                  &
    +35d0*gx(0)*gy(8)*gz(0)+140d0*gx(2)*gy(6)*gz(0)                     &
    +210d0*gx(4)*gy(4)*gz(0)+140d0*gx(6)*gy(2)*gz(0)                    &
    +35d0*gx(8)*gy(0)*gz(0))/128d0
Q81s=Q81s+f*(192d0*gx(0)*gy(1)*gz(7)-1008d0*gx(0)*gy(3)*gz(5)           &
    -1008d0*gx(2)*gy(1)*gz(5)+840d0*gx(0)*gy(5)*gz(3)                   &
    +1680d0*gx(2)*gy(3)*gz(3)+840d0*gx(4)*gy(1)*gz(3)                   &
    -105d0*gx(0)*gy(7)*gz(1)-315d0*gx(2)*gy(5)*gz(1)                    &
    -315d0*gx(4)*gy(3)*gz(1)-105d0*gx(6)*gy(1)*gz(1))/32d0
Q81c=Q81c+f*(192d0*gx(1)*gy(0)*gz(7)-1008d0*gx(1)*gy(2)*gz(5)           &
    -1008d0*gx(3)*gy(0)*gz(5)+840d0*gx(1)*gy(4)*gz(3)                   &
    +1680d0*gx(3)*gy(2)*gz(3)+840d0*gx(5)*gy(0)*gz(3)                   &
    -105d0*gx(1)*gy(6)*gz(1)-315d0*gx(3)*gy(4)*gz(1)                    &
    -315d0*gx(5)*gy(2)*gz(1)-105d0*gx(7)*gy(0)*gz(1))/32d0
Q82s=Q82s+f*rt(2)*rt(5)*rt(7)*(192d0*gx(1)*gy(1)*gz(6)                  &
    -480d0*gx(1)*gy(3)*gz(4)-480d0*gx(3)*gy(1)*gz(4)                    &
    +180d0*gx(1)*gy(5)*gz(2)+360d0*gx(3)*gy(3)*gz(2)                    &
    +180d0*gx(5)*gy(1)*gz(2)-6d0*gx(1)*gy(7)*gz(0)                      &
    -18d0*gx(3)*gy(5)*gz(0)-18d0*gx(5)*gy(3)*gz(0)                      &
    -6d0*gx(7)*gy(1)*gz(0))/64d0
Q82c=Q82c+f*rt(2)*rt(5)*rt(7)*(-96d0*gx(0)*gy(2)*gz(6)                  &
    +96d0*gx(2)*gy(0)*gz(6)+240d0*gx(0)*gy(4)*gz(4)                     &
    -240d0*gx(4)*gy(0)*gz(4)-90d0*gx(0)*gy(6)*gz(2)                     &
    -90d0*gx(2)*gy(4)*gz(2)+90d0*gx(4)*gy(2)*gz(2)                      &
    +90d0*gx(6)*gy(0)*gz(2)+3d0*gx(0)*gy(8)*gz(0)                       &
    +6d0*gx(2)*gy(6)*gz(0)-6d0*gx(6)*gy(2)*gz(0)-3d0*gx(8)*gy(0)*gz(0)) &
    /64d0
Q83s=Q83s+f*rt(3)*rt(5)*rt(7)*rt(11)*(-16d0*gx(0)*gy(3)*gz(5)           &
    +48d0*gx(2)*gy(1)*gz(5)+20d0*gx(0)*gy(5)*gz(3)                      &
    -40d0*gx(2)*gy(3)*gz(3)-60d0*gx(4)*gy(1)*gz(3)                      &
    -3d0*gx(0)*gy(7)*gz(1)+3d0*gx(2)*gy(5)*gz(1)                        &
    +15d0*gx(4)*gy(3)*gz(1)+9d0*gx(6)*gy(1)*gz(1))/32d0
Q83c=Q83c+f*rt(3)*rt(5)*rt(7)*rt(11)*(-48d0*gx(1)*gy(2)*gz(5)           &
    +16d0*gx(3)*gy(0)*gz(5)+60d0*gx(1)*gy(4)*gz(3)                      &
    +40d0*gx(3)*gy(2)*gz(3)-20d0*gx(5)*gy(0)*gz(3)                      &
    -9d0*gx(1)*gy(6)*gz(1)-15d0*gx(3)*gy(4)*gz(1)                       &
    -3d0*gx(5)*gy(2)*gz(1)+3d0*gx(7)*gy(0)*gz(1))/32d0
Q84s=Q84s+f*rt(7)*rt(11)*(-480d0*gx(1)*gy(3)*gz(4)                      &
    +480d0*gx(3)*gy(1)*gz(4)+288d0*gx(1)*gy(5)*gz(2)                    &
    -288d0*gx(5)*gy(1)*gz(2)-12d0*gx(1)*gy(7)*gz(0)                     &
    -12d0*gx(3)*gy(5)*gz(0)+12d0*gx(5)*gy(3)*gz(0)                      &
    +12d0*gx(7)*gy(1)*gz(0))/64d0
Q84c=Q84c+f*rt(7)*rt(11)*(120d0*gx(0)*gy(4)*gz(4)                       &
    -720d0*gx(2)*gy(2)*gz(4)+120d0*gx(4)*gy(0)*gz(4)                    &
    -72d0*gx(0)*gy(6)*gz(2)+360d0*gx(2)*gy(4)*gz(2)                     &
    +360d0*gx(4)*gy(2)*gz(2)-72d0*gx(6)*gy(0)*gz(2)                     &
    +3d0*gx(0)*gy(8)*gz(0)-12d0*gx(2)*gy(6)*gz(0)                       &
    -30d0*gx(4)*gy(4)*gz(0)-12d0*gx(6)*gy(2)*gz(0)                      &
    +3d0*gx(8)*gy(0)*gz(0))/64d0
Q85s=Q85s+f*rt(7)*rt(11)*rt(13)*(12d0*gx(0)*gy(5)*gz(3)                 &
    -120d0*gx(2)*gy(3)*gz(3)+60d0*gx(4)*gy(1)*gz(3)                     &
    -3d0*gx(0)*gy(7)*gz(1)+27d0*gx(2)*gy(5)*gz(1)                       &
    +15d0*gx(4)*gy(3)*gz(1)-15d0*gx(6)*gy(1)*gz(1))/32d0
Q85c=Q85c+f*rt(7)*rt(11)*rt(13)*(60d0*gx(1)*gy(4)*gz(3)                 &
    -120d0*gx(3)*gy(2)*gz(3)+12d0*gx(5)*gy(0)*gz(3)                     &
    -15d0*gx(1)*gy(6)*gz(1)+15d0*gx(3)*gy(4)*gz(1)                      &
    +27d0*gx(5)*gy(2)*gz(1)-3d0*gx(7)*gy(0)*gz(1))/32d0
Q86s=Q86s+f*rt(2)*rt(3)*rt(11)*rt(13)*(84d0*gx(1)*gy(5)*gz(2)           &
    -280d0*gx(3)*gy(3)*gz(2)+84d0*gx(5)*gy(1)*gz(2)                     &
    -6d0*gx(1)*gy(7)*gz(0)+14d0*gx(3)*gy(5)*gz(0)                       &
    +14d0*gx(5)*gy(3)*gz(0)-6d0*gx(7)*gy(1)*gz(0))/64d0
Q86c=Q86c+f*rt(2)*rt(3)*rt(11)*rt(13)*(-14d0*gx(0)*gy(6)*gz(2)          &
    +210d0*gx(2)*gy(4)*gz(2)-210d0*gx(4)*gy(2)*gz(2)                    &
    +14d0*gx(6)*gy(0)*gz(2)+gx(0)*gy(8)*gz(0)-14d0*gx(2)*gy(6)*gz(0)    &
    +14d0*gx(6)*gy(2)*gz(0)-gx(8)*gy(0)*gz(0))/64d0
Q87s=Q87s+f*rt(5)*rt(11)*rt(13)*(-3d0*gx(0)*gy(7)*gz(1)                 &
    +63d0*gx(2)*gy(5)*gz(1)-105d0*gx(4)*gy(3)*gz(1)                     &
    +21d0*gx(6)*gy(1)*gz(1))/32d0
Q87c=Q87c+f*rt(5)*rt(11)*rt(13)*(-21d0*gx(1)*gy(6)*gz(1)                &
    +105d0*gx(3)*gy(4)*gz(1)-63d0*gx(5)*gy(2)*gz(1)                     &
    +3d0*gx(7)*gy(0)*gz(1))/32d0
Q88s=Q88s+f*rt(5)*rt(11)*rt(13)*(-24d0*gx(1)*gy(7)*gz(0)                &
    +168d0*gx(3)*gy(5)*gz(0)-168d0*gx(5)*gy(3)*gz(0)                    &
    +24d0*gx(7)*gy(1)*gz(0))/128d0
Q88c=Q88c+f*rt(5)*rt(11)*rt(13)*(3d0*gx(0)*gy(8)*gz(0)                  &
    -84d0*gx(2)*gy(6)*gz(0)+210d0*gx(4)*gy(4)*gz(0)                     &
    -84d0*gx(6)*gy(2)*gz(0)+3d0*gx(8)*gy(0)*gz(0))/128d0

107                                                                     &
Q70=Q70+f*(16d0*gx(0)*gy(0)*gz(7)-168d0*gx(0)*gy(2)*gz(5)               &
    -168d0*gx(2)*gy(0)*gz(5)+210d0*gx(0)*gy(4)*gz(3)                    &
    +420d0*gx(2)*gy(2)*gz(3)+210d0*gx(4)*gy(0)*gz(3)                    &
    -35d0*gx(0)*gy(6)*gz(1)-105d0*gx(2)*gy(4)*gz(1)                     &
    -105d0*gx(4)*gy(2)*gz(1)-35d0*gx(6)*gy(0)*gz(1))/16d0
Q71s=Q71s+f*rt(7)*(64d0*gx(0)*gy(1)*gz(6)-240d0*gx(0)*gy(3)*gz(4)       &
    -240d0*gx(2)*gy(1)*gz(4)+120d0*gx(0)*gy(5)*gz(2)                    &
    +240d0*gx(2)*gy(3)*gz(2)+120d0*gx(4)*gy(1)*gz(2)                    &
    -5d0*gx(0)*gy(7)*gz(0)-15d0*gx(2)*gy(5)*gz(0)                       &
    -15d0*gx(4)*gy(3)*gz(0)-5d0*gx(6)*gy(1)*gz(0))/32d0
Q71c=Q71c+f*rt(7)*(64d0*gx(1)*gy(0)*gz(6)-240d0*gx(1)*gy(2)*gz(4)       &
    -240d0*gx(3)*gy(0)*gz(4)+120d0*gx(1)*gy(4)*gz(2)                    &
    +240d0*gx(3)*gy(2)*gz(2)+120d0*gx(5)*gy(0)*gz(2)                    &
    -5d0*gx(1)*gy(6)*gz(0)-15d0*gx(3)*gy(4)*gz(0)                       &
    -15d0*gx(5)*gy(2)*gz(0)-5d0*gx(7)*gy(0)*gz(0))/32d0
Q72s=Q72s+f*rt(2)*rt(3)*rt(7)*(96d0*gx(1)*gy(1)*gz(5)                   &
    -160d0*gx(1)*gy(3)*gz(3)-160d0*gx(3)*gy(1)*gz(3)                    &
    +30d0*gx(1)*gy(5)*gz(1)+60d0*gx(3)*gy(3)*gz(1)                      &
    +30d0*gx(5)*gy(1)*gz(1))/32d0
Q72c=Q72c+f*rt(2)*rt(3)*rt(7)*(-48d0*gx(0)*gy(2)*gz(5)                  &
    +48d0*gx(2)*gy(0)*gz(5)+80d0*gx(0)*gy(4)*gz(3)                      &
    -80d0*gx(4)*gy(0)*gz(3)-15d0*gx(0)*gy(6)*gz(1)                      &
    -15d0*gx(2)*gy(4)*gz(1)+15d0*gx(4)*gy(2)*gz(1)                      &
    +15d0*gx(6)*gy(0)*gz(1))/32d0
Q73s=Q73s+f*rt(3)*rt(7)*(-80d0*gx(0)*gy(3)*gz(4)                        &
    +240d0*gx(2)*gy(1)*gz(4)+60d0*gx(0)*gy(5)*gz(2)                     &
    -120d0*gx(2)*gy(3)*gz(2)-180d0*gx(4)*gy(1)*gz(2)                    &
    -3d0*gx(0)*gy(7)*gz(0)+3d0*gx(2)*gy(5)*gz(0)                        &
    +15d0*gx(4)*gy(3)*gz(0)+9d0*gx(6)*gy(1)*gz(0))/32d0
Q73c=Q73c+f*rt(3)*rt(7)*(-240d0*gx(1)*gy(2)*gz(4)                       &
    +80d0*gx(3)*gy(0)*gz(4)+180d0*gx(1)*gy(4)*gz(2)                     &
    +120d0*gx(3)*gy(2)*gz(2)-60d0*gx(5)*gy(0)*gz(2)                     &
    -9d0*gx(1)*gy(6)*gz(0)-15d0*gx(3)*gy(4)*gz(0)                       &
    -3d0*gx(5)*gy(2)*gz(0)+3d0*gx(7)*gy(0)*gz(0))/32d0
Q74s=Q74s+f*rt(3)*rt(7)*rt(11)*(-40d0*gx(1)*gy(3)*gz(3)                 &
    +40d0*gx(3)*gy(1)*gz(3)+12d0*gx(1)*gy(5)*gz(1)                      &
    -12d0*gx(5)*gy(1)*gz(1))/16d0
Q74c=Q74c+f*rt(3)*rt(7)*rt(11)*(10d0*gx(0)*gy(4)*gz(3)                  &
    -60d0*gx(2)*gy(2)*gz(3)+10d0*gx(4)*gy(0)*gz(3)                      &
    -3d0*gx(0)*gy(6)*gz(1)+15d0*gx(2)*gy(4)*gz(1)                       &
    +15d0*gx(4)*gy(2)*gz(1)-3d0*gx(6)*gy(0)*gz(1))/16d0
Q75s=Q75s+f*rt(3)*rt(7)*rt(11)*(12d0*gx(0)*gy(5)*gz(2)                  &
    -120d0*gx(2)*gy(3)*gz(2)+60d0*gx(4)*gy(1)*gz(2)-gx(0)*gy(7)*gz(0)   &
    +9d0*gx(2)*gy(5)*gz(0)+5d0*gx(4)*gy(3)*gz(0)-5d0*gx(6)*gy(1)*gz(0)) &
    /32d0
Q75c=Q75c+f*rt(3)*rt(7)*rt(11)*(60d0*gx(1)*gy(4)*gz(2)                  &
    -120d0*gx(3)*gy(2)*gz(2)+12d0*gx(5)*gy(0)*gz(2)                     &
    -5d0*gx(1)*gy(6)*gz(0)+5d0*gx(3)*gy(4)*gz(0)+9d0*gx(5)*gy(2)*gz(0)  &
    -gx(7)*gy(0)*gz(0))/32d0
Q76s=Q76s+f*rt(2)*rt(3)*rt(7)*rt(11)*rt(13)*(6d0*gx(1)*gy(5)*gz(1)      &
    -20d0*gx(3)*gy(3)*gz(1)+6d0*gx(5)*gy(1)*gz(1))/32d0
Q76c=Q76c+f*rt(2)*rt(3)*rt(7)*rt(11)*rt(13)*(-gx(0)*gy(6)*gz(1)         &
    +15d0*gx(2)*gy(4)*gz(1)-15d0*gx(4)*gy(2)*gz(1)+gx(6)*gy(0)*gz(1))   &
    /32d0
Q77s=Q77s+f*rt(3)*rt(11)*rt(13)*(-gx(0)*gy(7)*gz(0)                     &
    +21d0*gx(2)*gy(5)*gz(0)-35d0*gx(4)*gy(3)*gz(0)                      &
    +7d0*gx(6)*gy(1)*gz(0))/32d0
Q77c=Q77c+f*rt(3)*rt(11)*rt(13)*(-7d0*gx(1)*gy(6)*gz(0)                 &
    +35d0*gx(3)*gy(4)*gz(0)-21d0*gx(5)*gy(2)*gz(0)+gx(7)*gy(0)*gz(0))   &
    /32d0

106                                                                     &
Q60=Q60+f*(16d0*gx(0)*gy(0)*gz(6)-120d0*gx(0)*gy(2)*gz(4)               &
    -120d0*gx(2)*gy(0)*gz(4)+90d0*gx(0)*gy(4)*gz(2)                     &
    +180d0*gx(2)*gy(2)*gz(2)+90d0*gx(4)*gy(0)*gz(2)                     &
    -5d0*gx(0)*gy(6)*gz(0)-15d0*gx(2)*gy(4)*gz(0)                       &
    -15d0*gx(4)*gy(2)*gz(0)-5d0*gx(6)*gy(0)*gz(0))/16d0
Q61s=Q61s+f*rt(3)*rt(7)*(8d0*gx(0)*gy(1)*gz(5)-20d0*gx(0)*gy(3)*gz(3)   &
    -20d0*gx(2)*gy(1)*gz(3)+5d0*gx(0)*gy(5)*gz(1)                       &
    +10d0*gx(2)*gy(3)*gz(1)+5d0*gx(4)*gy(1)*gz(1))/8d0
Q61c=Q61c+f*rt(3)*rt(7)*(8d0*gx(1)*gy(0)*gz(5)-20d0*gx(1)*gy(2)*gz(3)   &
    -20d0*gx(3)*gy(0)*gz(3)+5d0*gx(1)*gy(4)*gz(1)                       &
    +10d0*gx(3)*gy(2)*gz(1)+5d0*gx(5)*gy(0)*gz(1))/8d0
Q62s=Q62s+f*rt(2)*rt(3)*rt(5)*rt(7)*(32d0*gx(1)*gy(1)*gz(4)             &
    -32d0*gx(1)*gy(3)*gz(2)-32d0*gx(3)*gy(1)*gz(2)                      &
    +2d0*gx(1)*gy(5)*gz(0)+4d0*gx(3)*gy(3)*gz(0)+2d0*gx(5)*gy(1)*gz(0)) &
    /32d0
Q62c=Q62c+f*rt(2)*rt(3)*rt(5)*rt(7)*(-16d0*gx(0)*gy(2)*gz(4)            &
    +16d0*gx(2)*gy(0)*gz(4)+16d0*gx(0)*gy(4)*gz(2)                      &
    -16d0*gx(4)*gy(0)*gz(2)-gx(0)*gy(6)*gz(0)-gx(2)*gy(4)*gz(0)         &
    +gx(4)*gy(2)*gz(0)+gx(6)*gy(0)*gz(0))/32d0
Q63s=Q63s+f*rt(2)*rt(3)*rt(5)*rt(7)*(-8d0*gx(0)*gy(3)*gz(3)             &
    +24d0*gx(2)*gy(1)*gz(3)+3d0*gx(0)*gy(5)*gz(1)                       &
    -6d0*gx(2)*gy(3)*gz(1)-9d0*gx(4)*gy(1)*gz(1))/16d0
Q63c=Q63c+f*rt(2)*rt(3)*rt(5)*rt(7)*(-24d0*gx(1)*gy(2)*gz(3)            &
    +8d0*gx(3)*gy(0)*gz(3)+9d0*gx(1)*gy(4)*gz(1)+6d0*gx(3)*gy(2)*gz(1)  &
    -3d0*gx(5)*gy(0)*gz(1))/16d0
Q64s=Q64s+f*rt(7)*(-120d0*gx(1)*gy(3)*gz(2)+120d0*gx(3)*gy(1)*gz(2)     &
    +12d0*gx(1)*gy(5)*gz(0)-12d0*gx(5)*gy(1)*gz(0))/16d0
Q64c=Q64c+f*rt(7)*(30d0*gx(0)*gy(4)*gz(2)-180d0*gx(2)*gy(2)*gz(2)       &
    +30d0*gx(4)*gy(0)*gz(2)-3d0*gx(0)*gy(6)*gz(0)                       &
    +15d0*gx(2)*gy(4)*gz(0)+15d0*gx(4)*gy(2)*gz(0)                      &
    -3d0*gx(6)*gy(0)*gz(0))/16d0
Q65s=Q65s+f*rt(2)*rt(7)*rt(11)*(3d0*gx(0)*gy(5)*gz(1)                   &
    -30d0*gx(2)*gy(3)*gz(1)+15d0*gx(4)*gy(1)*gz(1))/16d0
Q65c=Q65c+f*rt(2)*rt(7)*rt(11)*(15d0*gx(1)*gy(4)*gz(1)                  &
    -30d0*gx(3)*gy(2)*gz(1)+3d0*gx(5)*gy(0)*gz(1))/16d0
Q66s=Q66s+f*rt(2)*rt(3)*rt(7)*rt(11)*(6d0*gx(1)*gy(5)*gz(0)             &
    -20d0*gx(3)*gy(3)*gz(0)+6d0*gx(5)*gy(1)*gz(0))/32d0
Q66c=Q66c+f*rt(2)*rt(3)*rt(7)*rt(11)*(-gx(0)*gy(6)*gz(0)                &
    +15d0*gx(2)*gy(4)*gz(0)-15d0*gx(4)*gy(2)*gz(0)+gx(6)*gy(0)*gz(0))   &
    /32d0

105                                                                     &
Q50=Q50+f*(8d0*gx(0)*gy(0)*gz(5)-40d0*gx(0)*gy(2)*gz(3)                 &
    -40d0*gx(2)*gy(0)*gz(3)+15d0*gx(0)*gy(4)*gz(1)                      &
    +30d0*gx(2)*gy(2)*gz(1)+15d0*gx(4)*gy(0)*gz(1))/8d0
Q51s=Q51s+f*rt(3)*rt(5)*(8d0*gx(0)*gy(1)*gz(4)-12d0*gx(0)*gy(3)*gz(2)   &
    -12d0*gx(2)*gy(1)*gz(2)+gx(0)*gy(5)*gz(0)+2d0*gx(2)*gy(3)*gz(0)     &
    +gx(4)*gy(1)*gz(0))/8d0
Q51c=Q51c+f*rt(3)*rt(5)*(8d0*gx(1)*gy(0)*gz(4)-12d0*gx(1)*gy(2)*gz(2)   &
    -12d0*gx(3)*gy(0)*gz(2)+gx(1)*gy(4)*gz(0)+2d0*gx(3)*gy(2)*gz(0)     &
    +gx(5)*gy(0)*gz(0))/8d0
Q52s=Q52s+f*rt(3)*rt(5)*rt(7)*(4d0*gx(1)*gy(1)*gz(3)                    &
    -2d0*gx(1)*gy(3)*gz(1)-2d0*gx(3)*gy(1)*gz(1))/4d0
Q52c=Q52c+f*rt(3)*rt(5)*rt(7)*(-2d0*gx(0)*gy(2)*gz(3)                   &
    +2d0*gx(2)*gy(0)*gz(3)+gx(0)*gy(4)*gz(1)-gx(4)*gy(0)*gz(1))/4d0
Q53s=Q53s+f*rt(2)*rt(5)*rt(7)*(-8d0*gx(0)*gy(3)*gz(2)                   &
    +24d0*gx(2)*gy(1)*gz(2)+gx(0)*gy(5)*gz(0)-2d0*gx(2)*gy(3)*gz(0)     &
    -3d0*gx(4)*gy(1)*gz(0))/16d0
Q53c=Q53c+f*rt(2)*rt(5)*rt(7)*(-24d0*gx(1)*gy(2)*gz(2)                  &
    +8d0*gx(3)*gy(0)*gz(2)+3d0*gx(1)*gy(4)*gz(0)+2d0*gx(3)*gy(2)*gz(0)  &
    -gx(5)*gy(0)*gz(0))/16d0
Q54s=Q54s+f*rt(5)*rt(7)*(-12d0*gx(1)*gy(3)*gz(1)                        &
    +12d0*gx(3)*gy(1)*gz(1))/8d0
Q54c=Q54c+f*rt(5)*rt(7)*(3d0*gx(0)*gy(4)*gz(1)-18d0*gx(2)*gy(2)*gz(1)   &
    +3d0*gx(4)*gy(0)*gz(1))/8d0
Q55s=Q55s+f*rt(2)*rt(7)*(3d0*gx(0)*gy(5)*gz(0)-30d0*gx(2)*gy(3)*gz(0)   &
    +15d0*gx(4)*gy(1)*gz(0))/16d0
Q55c=Q55c+f*rt(2)*rt(7)*(15d0*gx(1)*gy(4)*gz(0)-30d0*gx(3)*gy(2)*gz(0)  &
    +3d0*gx(5)*gy(0)*gz(0))/16d0

!  Hexadecapole terms
104                                                                     &
Q40=Q40+f*(8d0*gx(0)*gy(0)*gz(4)-24d0*gx(0)*gy(2)*gz(2)                 &
    -24d0*gx(2)*gy(0)*gz(2)+3d0*gx(0)*gy(4)*gz(0)                       &
    +6d0*gx(2)*gy(2)*gz(0)+3d0*gx(4)*gy(0)*gz(0))/8d0
Q41s=Q41s+f*rt(2)*rt(5)*(4d0*gx(0)*gy(1)*gz(3)-3d0*gx(0)*gy(3)*gz(1)    &
    -3d0*gx(2)*gy(1)*gz(1))/4d0
Q41c=Q41c+f*rt(2)*rt(5)*(4d0*gx(1)*gy(0)*gz(3)-3d0*gx(1)*gy(2)*gz(1)    &
    -3d0*gx(3)*gy(0)*gz(1))/4d0
Q42s=Q42s+f*rt(5)*(12d0*gx(1)*gy(1)*gz(2)-2d0*gx(1)*gy(3)*gz(0)         &
    -2d0*gx(3)*gy(1)*gz(0))/4d0
Q42c=Q42c+f*rt(5)*(-6d0*gx(0)*gy(2)*gz(2)+6d0*gx(2)*gy(0)*gz(2)         &
    +gx(0)*gy(4)*gz(0)-gx(4)*gy(0)*gz(0))/4d0
Q43s=Q43s+f*rt(2)*rt(5)*rt(7)*(-gx(0)*gy(3)*gz(1)                       &
    +3d0*gx(2)*gy(1)*gz(1))/4d0
Q43c=Q43c+f*rt(2)*rt(5)*rt(7)*(-3d0*gx(1)*gy(2)*gz(1)                   &
    +gx(3)*gy(0)*gz(1))/4d0
Q44s=Q44s+f*rt(5)*rt(7)*(-4d0*gx(1)*gy(3)*gz(0)+4d0*gx(3)*gy(1)*gz(0))  &
    /8d0
Q44c=Q44c+f*rt(5)*rt(7)*(gx(0)*gy(4)*gz(0)-6d0*gx(2)*gy(2)*gz(0)        &
    +gx(4)*gy(0)*gz(0))/8d0

!  Octopole terms
103                                                                     &
Q30=Q30+f*(2d0*gx(0)*gy(0)*gz(3)-3d0*gx(0)*gy(2)*gz(1)                  &
    -3d0*gx(2)*gy(0)*gz(1))/2d0
Q31s=Q31s+f*rt(2)*rt(3)*(4d0*gx(0)*gy(1)*gz(2)-gx(0)*gy(3)*gz(0)        &
    -gx(2)*gy(1)*gz(0))/4d0
Q31c=Q31c+f*rt(2)*rt(3)*(4d0*gx(1)*gy(0)*gz(2)-gx(1)*gy(2)*gz(0)        &
    -gx(3)*gy(0)*gz(0))/4d0
Q32s=Q32s+f*rt(3)*rt(5)*(2d0*gx(1)*gy(1)*gz(1))/2d0
Q32c=Q32c+f*rt(3)*rt(5)*(-gx(0)*gy(2)*gz(1)+gx(2)*gy(0)*gz(1))/2d0
Q33s=Q33s+f*rt(2)*rt(5)*(-gx(0)*gy(3)*gz(0)+3d0*gx(2)*gy(1)*gz(0))/4d0
Q33c=Q33c+f*rt(2)*rt(5)*(-3d0*gx(1)*gy(2)*gz(0)+gx(3)*gy(0)*gz(0))/4d0

!  Quadrupole terms
102                                                                     &
Q20=Q20+f*(2d0*gx(0)*gy(0)*gz(2)-gx(0)*gy(2)*gz(0)-gx(2)*gy(0)*gz(0))/2d0
Q21s=Q21s+f*rt(3)*(gx(0)*gy(1)*gz(1))
Q21c=Q21c+f*rt(3)*(gx(1)*gy(0)*gz(1))
Q22s=Q22s+f*rt(3)*(2d0*gx(1)*gy(1)*gz(0))/2d0
Q22c=Q22c+f*rt(3)*(-gx(0)*gy(2)*gz(0)+gx(2)*gy(0)*gz(0))/2d0

!  Dipole terms
101   Q10=Q10+f*(gx(0)*gy(0)*gz(1))
Q11s=Q11s+f*(gx(0)*gy(1)*gz(0))
Q11c=Q11c+f*(gx(1)*gy(0)*gz(0))
!  Monopole term
100   Q00=Q00+F*GX(0)*GY(0)*GZ(0)

      END SUBROUTINE addqlm

!-----------------------------------------------------------------SHIFTQ

SUBROUTINE shiftq (q1,l1,m1, q2,m2, x,y,z)
IMPLICIT NONE
REAL(dp), INTENT(IN) :: q1(121), x, y, z
INTEGER, INTENT(IN) :: l1, m1, m2
REAL(dp), INTENT(INOUT) :: q2(121)
!-----------------------------------------------------
!     Copyright A J Stone University of Cambridge 1983
!     Modifications and interface to Cadpac, R D Amos
!     Version for Cadpac5 , R D Amos, June 1990
!-----------------------------------------------------
!  Shift the multipoles Q1 relative to the point (X,Y,Z) to the point
!  (0,0,0) and add them to the multipole expansion Q2
!  Multipoles of ranks L1 through M1 are to be tranferred from Q1
!  and added to Q2, keeping ranks up to M2 in the transferred expansion.
!  Q2(L,M) = sum(K,Q) sqrt | (L+M) (L-M) | * Q1(K,Q) * R(L-K,M-Q)
!                          | (K+Q) (K-Q) |
!  where the quantities in the square root are binomial coefficients
!  and the R() are regular solid harmonics of (x,y,z).

!  Real multipoles and harmonics are in the order
!  Q00, Q10, Q11c, Q11s, Q20, ...
!  Complex multipoles and harmonics are in the order
!  Q00, Q1-1, Q10, Q11, Q2-2, ...
REAL(dp), PARAMETER :: RTHALF=0.7071067811865475244D0, EPS=0.0D0

REAL(dp) :: R(121), r2, a, s
INTEGER :: i, jb, k, k1, kb, km, kmax, l, lb, lm, m, n2,         &
    qmin, qmax, qq, t1, t2
REAL(dp) :: RC(121,2), QC(121,2), QZ(121,2)

IF (L1 .GT. M1 .OR. L1 .GT. M2) RETURN
!  Estimate largest significant transferred multipole.  The magnitude of
!  the contribution of the Q(k,q) to Q2(l,m) is of order
!  |Q(k)|*R**(l-k), multiplied by binomial coefficient factors which we
!  estimate as 2**(l-k).  If the total of such estimates for rank l is
!  greater than sqrt(EPS), transfers up to rank l are calculated
!  explicitly.
!  If EPS is zero, this procedure is bypassed.
n2=m2
k1=max0(1,l1)
if (eps .ne. 0.0d0 .and. m2 .gt. 0) then
  r2=4.0d0*(x**2+y**2+z**2)
  n2=0
  a=0.0d0
  if (l1 .eq. 0) a=q1(1)**2
  do k=k1,m2
    a=a*r2
    if (k .le. m1) then
      t1=k**2+1
      t2=(k+1)**2
      do i=t1,t2
        a=a+q1(i)**2
      end do
    end if
    if (a .gt. eps) n2=k
  end do
end if
!  Evaluate solid harmonics in real form
call solidh (x,y,z, n2, r,121)
!  Construct complex solid harmonics RC
rc(1,1)=r(1)
rc(1,2)=0.0d0
do k=1,n2
  kb=k**2+k+1
  km=k**2+1
  rc(kb,1)=r(km)
  rc(kb,2)=0.0d0
  km=km+1
  s=rthalf
  do m=1,k
    s=-s
    rc(kb-m,1)=rthalf*r(km)
    rc(kb-m,2)=-rthalf*r(km+1)
    rc(kb+m,1)=s*r(km)
    rc(kb+m,2)=s*r(km+1)
    km=km+2
  end do
end do
!  Construct complex multipoles QC corresponding to original
!  real multipoles Q1
if (l1 .eq. 0) then
  qc(1,1)=q1(1)
  qc(1,2)=0.0d0
endif
if (m1 .gt. 0) then
  do k=k1,m1
    kb=k**2+k+1
    km=k**2+1
    qc(kb,1)=q1(km)
    qc(kb,2)=0.0d0
    km=km+1
    s=rthalf
    do m=1,k
      s=-s
      qc(kb-m,1)=rthalf*q1(km)
      qc(kb-m,2)=-rthalf*q1(km+1)
      qc(kb+m,1)=s*q1(km)
      qc(kb+m,2)=s*q1(km+1)
      km=km+2
    end do
  end do
end if
!  Construct shifted complex multipoles QZ (only for non-negative M)
if (l1 .eq. 0) then
  qz(1,1)=qc(1,1)
  qz(1,2)=qc(1,2)
endif
do l=k1,n2
  kmax=min0(l,m1)
  lb=l**2+l+1
  lm=lb
  m=0
  do while (m .le. l)
    qz(lm,1)=0.0d0
    qz(lm,2)=0.0d0
    if (l1 .eq. 0) then
      qz(lm,1)=qc(1,1)*rc(lm,1)-qc(1,2)*rc(lm,2)
      qz(lm,2)=qc(1,1)*rc(lm,2)+qc(1,2)*rc(lm,1)
    endif
      do k=k1,kmax
        qmin=max0(-k,k-l+m)
        qmax=min0(k,l-k+m)
        kb=k**2+k+1
        jb=(l-k)**2+(l-k)+1
        do qq=qmin,qmax
          qz(lm,1)=qz(lm,1)+rtbinom(l+m,k+qq)*rtbinom(l-m,k-qq)            &
              *(qc(kb+qq,1)*rc(jb+m-qq,1)-qc(kb+qq,2)*rc(jb+m-qq,2))
          qz(lm,2)=qz(lm,2)+rtbinom(l+m,k+qq)*rtbinom(l-m,k-qq)            &
              *(qc(kb+qq,1)*rc(jb+m-qq,2)+qc(kb+qq,2)*rc(jb+m-qq,1))
        end do
      end do
      m=m+1
    lm=lm+1
  end do
end do
!  Construct real multipoles and add to Q2
if (l1 .eq. 0) q2(1)=q2(1)+qz(1,1)
do k=k1,n2
  kb=k**2+k+1
  km=k**2+1
  q2(km)=q2(km)+qz(kb,1)
  s=1.0d0/rthalf
  km=km+1
  do m=1,k
    s=-s
    q2(km)=q2(km)+s*qz(kb+m,1)
    q2(km+1)=q2(km+1)+s*qz(kb+m,2)
    km=km+2
  end do
end do

END SUBROUTINE shiftq

!----------------------------------------------------------------- MOVEQ

SUBROUTINE moveq (x,y,z)
IMPLICIT NONE
REAL(dp), INTENT(IN) :: x, y, z
!-----------------------------------------------------
!     Copyright A J Stone University of Cambridge 1983
!     Modifications and interface to Cadpac, R D Amos
!     Version for Cadpac5 , R D Amos, June 1990
!-----------------------------------------------------
!  Move the set of multipoles in QT to the nearest multipole site
!  (but if two or more sites are almost equidistant, move a fraction
!  to each).

REAL(dp) :: qt
COMMON/BIG/qt(121)

REAL(dp) :: an, rr(maxs)
REAL(dp), PARAMETER :: eps=1d-6
INTEGER :: i, j, k, l, low, lp1sq, m(6), n, t1, t2

j=1
do i=1,ns
  rr(i)=((x-xs(1,i))**2 + (y-xs(2,i))**2 + (z-xs(3,i))**2)         &
      /radius(i)**2
  if (limit(i) .gt. limit(j)) j=i
end do
lp1sq=(lmax+1)**2
low=0
do
  k=j
  do i=1,ns
    if (rr(i) .lt. rr(k) .and. limit(i) .ge. low) k=i
  end do
  t1=low**2+1
  t2=(limit(k)+1)**2
  if (rr(k) .le. 1d-6) then
    do i=t1,t2
      q(i,k)=q(i,k)+qt(i)
      qt(i)=0.0d0
    end do
    if (limit(k) .ge. lmax) return
  else
    n=1
    m(1)=k
    do i=1,ns
      if (rr(i) .gt. rr(k)+eps .or. i .eq. k .or.                      &
          limit(i) .ne. limit(k) .or. limit(i) .lt. low) cycle
      n=n+1
      m(n)=i
    end do
    if (n .gt. 1) then
      an=1.0d0/dble(n)
      do k=t1,t2
        qt(k)=an*qt(k)
      end do
    end if
    do i=1,n
      k=m(i)
      call shiftq (qt,low,limit(k), q(1:,k),lmax,                      &
          x-xs(1,k),y-xs(2,k),z-xs(3,k))
    end do
    do i=t1,t2
      qt(i)=0.0d0
    end do
    if (limit(k) .ge. lmax) return
    t1=t2+1
    do i=1,n
      k=m(i)
      call shiftq(q(1:,k),limit(k)+1,lmax, qt,lmax,                    &
          xs(1,k)-x,xs(2,k)-y,xs(3,k)-z)
      do l=t1,lp1sq
        q(l,k)=0.0d0
      end do
    end do
  end if
  low=limit(k)+1
end do

END SUBROUTINE moveq

!-----------------------------------------------------------------SOLIDH

SUBROUTINE solidh (X,Y,Z, J, R,MAX)
USE input, ONLY : die
IMPLICIT NONE
REAL(dp), INTENT(IN) :: x, y, z
INTEGER, INTENT(IN) :: j, max
REAL(dp), INTENT(OUT) :: r(max)
!-----------------------------------------------------
!     Copyright A J Stone University of Cambridge 1983
!     Modifications and interface to Cadpac, R D Amos
!     Version for Cadpac5 , R D Amos, June 1990
!-----------------------------------------------------
!  Computes regular solid harmonics r**k Ckq(theta,phi) for ranks k up
!  to J, if J >= 0;
!  or irregular solid harmonics r**(-k-1) Ckq(theta,phi) for ranks k up
!  to |J|, if J < 0.

!  Locations in R are used as follows:
!        1    2    3    4    5    6    7    8    9   10   11  ...
!  kq = 00   10   11c  11s  20   21c  21s  22c  22s  30   31c ...
!  R(k,0) is real and is left in location k**2 + 1.
!  R(k,mc) and r(k,ms) are sqrt(2) times the real and imaginary parts
!  respectively of the complex solid harmonic R(k,-m)* = (-1)**m R(k,m),
!  and are left in locations K**2 + 2m and K**2 + 2m + 1 respectively.

INTEGER :: k, l, lk, ln, lp, m, n
REAL(dp) :: a2kp1, rr, rfx, rfy, rfz, s

l=iabs(j)
if ((l+1)**2 .gt. max) then
  write (6,'(a,i3)')                                           &
      'Insufficient array space for harmonics up to rank',L
  call die('Consult authors')
endif
rr=x**2+y**2+z**2
if (j .ge. 0) then
!  Regular
  r(1)=1.0d0
  r(2)=z
  r(3)=x
  r(4)=y
  rfz=z
  rfx=x
  rfy=y
else
!  Irregular
  rr=1.0d0/rr
  rfx=x*rr
  rfy=y*rr
  rfz=z*rr
  r(1)=dsqrt(rr)
  r(2)=rfz*r(1)
  r(3)=rfx*r(1)
  r(4)=rfy*r(1)
endif
!  Remaining values are found using recursion formulae, relating
!  the new set N to the current set K and the previous set P.
k=1
do while (k<l)
  n=k+1
  ln=n*n+1
  lk=k*k+1
  lp=(k-1)**2+1
  a2kp1=k+k+1
!  Obtain R(k+1,0) from R(k,0)*R(1,0) and R(k-1,0)
  r(ln)=(a2kp1*r(lk)*rfz-k*rr*r(lp))/(k+1)
  m=1
  ln=ln+1
  lk=lk+1
  lp=lp+1
  if (k .gt. 1) then
    do while (m<k)
!  Obtain R(k+1,m) from R(k,m)*R(1,0) and R(k-1,m)
      r(ln)=(a2kp1*r(lk)*rfz-rt(k+m)*rt(k-m)*rr*r(lp))                  &
                                           /(rt(n+m)*rt(n-m))
      r(ln+1)=(a2kp1*r(lk+1)*rfz-rt(k+m)*rt(k-m)*rr*r(lp+1))            &
                                           /(rt(n+m)*rt(n-m))
      m=m+1
      ln=ln+2
      lk=lk+2
      lp=lp+2
    end do
  end if
!  Obtain R(k+1,k) from R(k,k)*R(1,0)
  r(ln)=rt(n+k)*r(lk)*rfz
  r(ln+1)=rt(n+k)*r(lk+1)*rfz
  ln=ln+2
!  Obtain R(k+1,k+1) from R(k,k)*R(1,1)
  s=rt(n+k)/rt(n+n)
  r(ln)=s*(rfx*r(lk)-rfy*r(lk+1))
  r(ln+1)=s*(rfx*r(lk+1)+rfy*r(lk))
  k=k+1
end do

END SUBROUTINE solidh

!-----------------------------------------------------------------PRINTQ

SUBROUTINE printq(q,lm, linear, iw)
IMPLICIT NONE
REAL(dp), INTENT(IN) :: q(0:121)
INTEGER, INTENT(IN) :: lm, iw
LOGICAL, INTENT(IN) :: linear

!-----------------------------------------------------
!     Copyright A J Stone University of Cambridge 1983
!     Modifications and interface to Cadpac, R D Amos
!     Version for Cadpac5 , R D Amos, June 1990
!-----------------------------------------------------

INTEGER :: p(31)
LOGICAL :: big

CHARACTER(LEN=2) :: ql(15)=(/'Q1','Q2','Q3','Q4','Q5','Q6','Q7',       &
    'Q8','Q9','QA','QB','QC','QD','QE','QF'/)
CHARACTER(LEN=2) :: qm(31)=(/'0 ','1c','1s','2c','2s','3c','3s',       &
    '4c','4s','5c','5s','6c','6s','7c','7s','8c','8s','9c','9s',       &
    'Ac','As', 'Bc','Bs','Cc','Cs','Dc','Ds','Ec','Es','Fc','Fs'/)
CHARACTER :: LQ='Q'
INTEGER :: i, k, l, ll1, n
REAL(dp) :: qs, qsq, qf

if (linear) then
  write (iw,1006) (lq, k, q(k)*Qfactor(k), k=0,lm)
1006 FORMAT (3(5x, a1, i1, ' =', f14.8)/                               &
         3(5x, a1, i1, ' =', f14.8)/                                   &
         3(5x, a1, i1, ' =', f14.8)/                                   &
         5x, a1, i1, ' =', f14.8, 2(4x, a1, i2, ' =', f14.8)/          &
         (3(4x, a1, i2, ' =', f14.8)))
ELSE
  write (iw,'(19x, a, f11.6)') 'Q00  =', q(1)*Qfactor(0)
  k=1
  do l=1,lm
    qf=Qfactor(l)
    ll1=l+l+1
    n=0
    qsq=0d0
    do i=1,ll1
      qsq=qsq+q(k+i)**2
      if (dabs(q(k+i)) .ge. 5d-7) then
        n=n+1
        p(n)=i
      end if
    end do
    qs=dsqrt(qsq)
    big=(qs .ge. 1d3)
    if (n .gt. 0 .and. big) write (iw,1001) ql(l), qs*qf,              &
        (ql(l), qm(p(i)), q(p(i)+k)*qf, i=1,n)
    if (n .eq. 0 .and. .not. big) write (iw,1002) ql(l), qs*qf
    if (n .gt. 0 .and. .not. big) write (iw,1002) ql(l), qs*qf,        &
        (ql(l), qm(p(i)), q(p(i)+k)*qf, i=1,n)
1001 FORMAT ('|', a2, '| =', 1p,e11.3:, 3(2x, 2a2, ' =', e11.3:)       &
         / (17x, 3(2x, 2a2, ' =', e11.3:)))
1002 FORMAT ('|', a2, '| =', f11.6:, 3(2x, 2a2, ' =', f11.6:)          &
         / (17x, 3(2x, 2a2, ' =', f11.6:)))
    k=k+ll1
  end do
endif

END SUBROUTINE printq

!-----------------------------------------------------------------PLANES

SUBROUTINE planes

!  Look for special geometries (planar or linear). On exit, LINEAR is
!  true if the sites are all on a line parallel to the z axis. PLANAR
!  is true if the sites are all in a plane parallel to one of the
!  coordinate planes. PERP=1 if all sites have the same x, 2 if all sites
!  have the same y, 4 if all sites have the same z. For linear molecules,
!  PERP=3, i.e. 1+2. Note that this categorization applies to the sites,
!  not to the atoms, so it is only relevant for calculating matrices, not
!  for distributed multipoles.

IMPLICIT NONE

REAL(dp) :: xsq, ysq, zsq, xsum, ysum, zsum
INTEGER :: is

xsq=0d0
xsum=0d0
ysq=0d0
ysum=0d0
zsq=0d0
zsum=0d0
do is=1,ns
  xsum=xsum+xs(1,is)
  xsq=xsq+xs(1,is)**2
  ysum=ysum+xs(2,is)
  ysq=ysq+xs(2,is)**2
  zsum=zsum+xs(3,is)
  zsq=zsq+xs(3,is)**2
end do
perp=0
if (ns*xsq-xsum**2 .lt. 1D-8) perp=perp+1
if (ns*ysq-ysum**2 .lt. 1D-8) perp=perp+2
if (ns*zsq-zsum**2 .lt. 1D-8) perp=perp+4
if (perp == 7) then
  !  Atom
  linear=.false.
  general=.true.
else if (perp == 3) then
  linear=.true.
else if (perp .ne. 0) then
  planar=.true.
end if

END SUBROUTINE planes

END MODULE dma
