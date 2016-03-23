module gdma

!  Distributed Multipole Analysis for Gaussian Wavefunctions
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
!
!  This version has been modified by Andy Simmonett (03/16) to link
!  into Psi4, rather than serve as a standalone executable.

USE input
use iso_c_binding

USE dma
USE atom_grids, ONLY: debug_grid => debug
USE timing, ONLY: start_timer, timer, time_and_date
IMPLICIT NONE

INTEGER, PARAMETER :: dp=kind(1d0)

CHARACTER(LEN=100) :: file
CHARACTER(LEN=80) :: buffer
CHARACTER(LEN=20) :: key
CHARACTER(LEN=8) :: whichd="SCF"
CHARACTER(LEN=24) :: datestring

!  Maximum number of sites is number of atoms + nextra
INTEGER :: nextra=16
INTEGER :: ncoorb, maxl, cmax, nprim, nx, num, ich, mul
INTEGER, ALLOCATABLE :: shell_type(:)
INTEGER :: i, j, k, kp=0
LOGICAL :: eof, fchk, first, ok=.false.
INTEGER open_status, infile

REAL(dp), ALLOCATABLE :: densty(:,:), dtri(:)
INTEGER :: ir=5 ! Input stream

LOGICAL :: verbose=.false., debug(0:2)=.false.

CONTAINS


subroutine run_gdma(c_outfilename, c_datfilename) bind(c, name='run_gdma')
!character(kind=c_char,len=1), intent(in) :: c_outfilename
CHARACTER(kind=C_CHAR) :: c_outfilename(*), c_datfilename(*)
character(len=:), allocatable :: outfilename, datfilename
integer i, nchars
integer outfile

i = 1
do
   if (c_outfilename(i) == c_null_char) exit
   i = i + 1
end do
nchars = i - 1  ! Exclude null character from Fortran string
allocate(character(len=nchars) :: outfilename)
outfilename = transfer(c_outfilename(1:nchars), outfilename)
i = 1
do
   if (c_datfilename(i) == c_null_char) exit
   i = i + 1
end do
nchars = i - 1  ! Exclude null character from Fortran string
allocate(character(len=nchars) :: datfilename)
datfilename = transfer(c_datfilename(1:nchars), datfilename)

!
! Added file IO (ACS 03/16)
!
infile = 51
outfile = 52
open (unit=infile, file=datfilename, status='old', &
    iostat=open_status, action='read', position='rewind')
if ( open_status /= 0 ) then
    write(outfile, *) 'Could not open GDMA input for reading.', &
    'unit = ', infile
    stop
endif
open (unit=outfile, file=outfilename, status='old', &
    iostat=open_status, action='write', position='append')
if ( open_status /= 0 ) then
    write(outfile, *) 'Could not open psi4 output for writing.', &
    'unit = ', infile
    stop
endif
!!!
write(outfile, "(15x,a/)")                                             &
    "                      G D M A",                                   &
    "                  by Anthony Stone",                              &
    "            version 2.2.06, 22 June 2011",                        &
    "Distributed Multipoles from Gaussian wavefunctions"

call time_and_date(datestring)
write(outfile, "(/2A)") "Starting at ", datestring

call start_timer

punchfile="dma.punch"
nat=0
fchk=.false.
first=.true.
do
  call read_line(eof, infile)
  if (eof) exit
  call readu(key)
  select case(key)
  case("","NOTE","!")
    cycle
  case("VERBOSE")
    verbose=.true.
  case("QUIET")
    debug=.false.
    verbose=.false.
  case("DEBUG")
    debug(0)=.true.
    do while (item<nitems)
      call readi(k)
      if (k>0) then
        debug(k)=.true.
      else
        debug(-k)=.false.
      end if
    end do
    debug_grid=.true.
    verbose=.true.
  case("ANGSTROM")
    rfact=bohr
  case("BOHR")
    rfact=1d0
  case("SI")
    Qfactor(0)=echarge
    do k=1,20
      Qfactor(k)=Qfactor(k-1)*bohr
    end do
  case("AU")
    Qfactor=1d0
  case("COMMENT","TITLE")
    call reada(buffer)
    write(outfile, "(/a/)") trim(buffer)
  case("DENSITY")
    if (fchk) call die                                         &
        ("Specify density to use before reading data file",.true.)
    call readu(whichd)
  case("FILE","READ")
    nat=0
    fchk=.false.
    ok=.false.
    first=.true.
    if (allocated(dtri)) deallocate(dtri)
    do while (item<nitems)
      call readu(key)
      select case(key)
      case("DENSITY")
        call readu(whichd)
      case default
        call reread(-1)
        call reada(file)
        open(unit=9,file=file,status="old",iostat=k)
        if (k .ne. 0) then
          call die("Can't open file "//file,.true.)
        endif
      end select
    end do
    ir=9
    call get_data(whichd,ok,outfile)
    close(9)
    ir=5
    fchk=.true.
  case ("HERE")
    call get_data(whichd,ok,outfile)
    fchk=.true.
  case("NAMES")
    if (.not. fchk) call die                                   &
        ("Read data file before specifying atom names",.false.)
    call read_line(eof, infile)
    do i=1,nat
      call geta(name(i))
    end do
  case("GO","START","MULTIPOLES")
    if (.not. ok) then
      call die (trim(whichd)//" density not found",.false.)
    endif
    if (first) then
      ! convert density matrix to triangular form
      allocate(dtri(nx))
      k=0
      do i=1,num
        do j=1,i
          k=k+1
          dtri(k)=densty(i,j)
        end do
      end do
      deallocate(densty)
      first=.false.
    endif
    write(outfile, "(//2A/)") "Using "//trim(whichd)//" density matrix",  &
        " from file "//trim(file)
    call dma_main(dtri,kp,infile, outfile)
    !call timer
  case("RESET")
    nat=0
    fchk=.false.
    ok=.false.
    first=.true.
    deallocate(dtri)
    whichd="SCF"
  case("FINISH")
    exit
  case default
    call die("Keyword "//trim(key)//" not recognized",.true.)
  end select
end do

call time_and_date(datestring)
write(outfile, "(/2A)") "Finished at ", datestring
close(outfile)
close(infile)
end subroutine run_gdma


!-----------------------------------------------------------------------

SUBROUTINE get_data(whichd,ok,outfile)

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: whichd
INTEGER, INTENT(IN) :: outfile
LOGICAL, INTENT(OUT) :: ok

INTEGER :: atom, i, j, k, n, nn, aok
REAL(dp) :: e, rt3v2, td(5,6), tf(7,10), tg(9,15)
REAL(dp), ALLOCATABLE :: temp(:,:)
LOGICAL eof
CHARACTER :: text*40, buffer*80, ww*2, density_header*24, type*1
CHARACTER(LEN=8) :: dummy
REAL(dp), PARAMETER :: PI=3.14159265358979d0
CHARACTER(LEN=2), DIMENSION(54) :: element=(/"H ", "He",               &
    "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne",                    &
    "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar",                    &
    "K ", "Ca", "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni",        &
    "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",                    &
    "Rb", "Sr", "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",        &
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I ", "Xe"/)
REAL(dp), PARAMETER :: rt2=1.4142135623731d0,                          &
    rt3=1.73205080756888d0, rt5=2.23606797749979d0, rt7=2.64575131106459d0
REAL(dp), PARAMETER :: rt10=rt2*rt5, rt14=rt2*rt7, rt21=rt3*rt7,       &
    rt35=rt5*rt7
INTEGER, PARAMETER :: v400=1, v040=2, v004=3, v310=4, v301=5,          &
    v130=6, v031=7, v103=8, v013=9, v220=10, v202=11, v022=12,         &
    v211=13, v121=14, v112=15
CHARACTER(LEN=5) :: label(-5:5)=(/"h(s)","g(s)","f(s)","d(s)","sp  ",  &
    "s   ","p   ","d   ","f   ","g   ","h   "/)

!  Conversion from normalised spherical form to normalised Cartesian
!  Schlegel & Frisch, IJQC (1995) 54, 83-87.
rt3v2=rt3/2d0
!  d functions
!   1   2   3   4   5   6
!   xx  yy  zz  xy  xz  yz
!  200 020 002 110 101 011
td(1,:)=(/-0.5d0, -0.5d0, 1d0, 0d0, 0d0, 0d0/)
td(2,:)=(/0d0,    0d0,    0d0, 0d0, 1d0, 0d0/)
td(3,:)=(/0d0,    0d0,    0d0, 0d0, 0d0, 1d0/)
td(4,:)=(/rt3v2,  -rt3v2, 0d0, 0d0, 0d0, 0d0/)
td(5,:)=(/0d0,    0d0,    0d0, 1d0, 0d0, 0d0/)

!  f functions
!   1   2   3   4   5   6   7   8   9   10
!  xxx yyy zzz xxy xxz xyy yyz xzz yzz xyz
!  300 030 003 210 201 120 021 102 012 111
tf(:,:)=0d0
! 30
tf(1,3)=1d0; tf(1,5)=-1.5d0/sqrt(5d0); tf(1,7)=-1.5d0/sqrt(5d0)
! 31c ( F+1 in Gaussian notation )
tf(2,1)=-sqrt(3d0/8d0); tf(2,6)=-sqrt(1.2d0)/4d0; tf(2,8)=sqrt(1.2d0)
! 31s ( F-1 )
tf(3,2)=-sqrt(3d0/8d0); tf(3,4)=-sqrt(1.2d0)/4d0; tf(3,9)=sqrt(1.2d0)
! 32c ( F+2 )
tf(4,5)=sqrt(0.75d0); tf(4,7)=-sqrt(0.75d0)
! 32s
tf(5,10)=1d0
! 33c
tf(6,1)=sqrt(10d0)/4d0; tf(6,6)=-0.75d0*sqrt(2d0)
! 33s
tf(7,2)=-sqrt(10d0)/4d0; tf(7,4)=0.75d0*sqrt(2d0)

!  g functions
!   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
! xxxx yyyy zzzz xxxy xxxz xyyy yyyz xzzz yzzz xxyy xxzz yyzz xxyz xyyz xyzz
!  400  040  004  310  301  130  031  103  013  220  202  022  211  121  112
!     12   23  D     23 23,12  13  23,12 12  D      23  12
tg=0d0
!  40
tg(1,v400)=0.375d0; tg(1,v040)=0.375d0; tg(1,v004)=1d0
tg(1,v220)=0.75d0*rt3/rt35; tg(1,v202)=-3d0*rt3/rt35; tg(1,v022)=-3d0*rt3/rt35
!  41c
tg(2,v103)=rt10/rt7; tg(2,v301)=-0.75d0*rt10/rt7; tg(2,v121)=-0.75d0*rt2/rt7
!  41s
tg(3,v013)=rt10/rt7; tg(3,v031)=-0.75d0*rt10/rt7; tg(3,v211)=-0.75d0*rt2/rt7
!  42c
tg(4,v202)=1.5d0*rt3/rt7; tg(4,v022)=-1.5d0*rt3/rt7; tg(4,v400)=-rt5/4d0; tg(4,v040)=rt5/4d0
!  42s
tg(5,v112)=3d0/rt7; tg(5,v310)=-rt5/(2d0*rt7); tg(5,v130)=-rt5/(2d0*rt7)
!  43c
tg(6,v301)=rt10/4d0; tg(6,v121)=-0.75d0*rt2
!  43s
tg(7,v031)=-rt10/4d0; tg(7,v211)=0.75d0*rt2
!  44c
tg(8,v400)=rt35/8d0; tg(8,v040)=rt35/8d0; tg(8,v220)=-0.75*rt3
!  44s
tg(9,v310)=rt5/2d0; tg(9,v130)=-rt5/2d0

! select case(whichg)
! case("94","G94")
!   rfact=1d0
! case("98","G98","03","G03")
!   rfact=1d0
! end select

ok=.false.
call stream(ir)
read (ir,"(A80/A80)") title(1), title(2)
density_header="Total "//trim(whichd)//" Density"
if (verbose) then
  write(outfile,"(a/a/a/)") "Gaussian header:", trim(title(1)), trim(title(2))
  write(outfile,"(3a)") 'Looking for "', trim(density_header), '"'
end if
do
  read (ir,"(A)",iostat=k) buffer
  if (k .ne. 0) then
    call stream(5)
    exit
  end if
  if (debug(0)) write(outfile, "(a)") buffer
  text=buffer(1:40)
  type=buffer(44:44)
  ww=buffer(48:49)
  if (ww .eq. "N=") then
    read(buffer,"(55X,I6)") nn
  else
    nn=0
  endif
  if (debug(0)) write(outfile, "(a)") text
  select case(text)
  case("")
    cycle
  case("END","End","end")
    call stream(5)
    exit
  case("Number of atoms")
    read(buffer,"(55X,I6)") nat
    if (verbose) write(outfile, "(i0,a)") nat, " atoms"
    if (allocated(zan)) deallocate(zan,c)
    allocate(zan(nat),c(3,nat),stat=aok)
    if (aok>0) call die("Can't allocate atom arrays")
    maxcen=nat
    maxs=maxcen+nextra  !  Arbitrary limit on number of sites
    if (allocated(name)) deallocate(name)
    allocate(name(maxs),stat=aok)
    if (aok>0) call die("Can't allocate site-name array")
  case("Charge")
    read(buffer,"(55X,I6)") ich
    if (verbose) write(outfile, "(a,i0)") "Charge ", ich
  case("Multiplicity")
    read(buffer,"(55X,I6)") mul
    if (verbose) write(outfile, "(a,i0)") "Multiplicity ", mul
  case("Number of basis functions")
    read(buffer,"(55X,I6)") ncoorb
    !  This number may be increased following conversion from
    !  spherical to cartesian
    if (verbose) write(outfile, "(i0,a)") ncoorb, " basis functions"
  case("Number of contracted shells")
    read(buffer,"(55X,I6)") nshell
    if (verbose) write(outfile, "(i0,a)") nshell, " shells"
    if (allocated(kstart))                                             &
        deallocate(kstart,katom,ktype,kng,kloc,kmin,kmax,shell_type)
    allocate (kstart(nshell), katom(nshell+1), ktype(nshell),          &
        kng(nshell), kloc(nshell), kmin(nshell), kmax(nshell),         &
        shell_type(nshell),stat=aok)
    if (aok>0) call die("Can't allocate shell arrays")
    shell_type=0
  case("Highest angular momentum")
    read(buffer,"(55X,I6)") maxl
    if (verbose) write(outfile, "(a,i0)") "Highest angular momentum ", maxl
    if (maxl .gt. 4) call die                                 &
        ("Sorry -- GDMA can only handle s, p, d, f and g basis functions",.false.)
  case("Largest degree of contraction")
    read(buffer,"(55X,I6)") cmax
    if (verbose) write(outfile, "(a,i0)") "Largest contraction depth ", cmax
    if (cmax .gt. 16) call die                                &
        ("Sorry -- maximum contraction depth is 16",.false.)
  case("Number of primitive shells")
    read(buffer,"(55X,I6)") nprim
    if (verbose) write(outfile, "(i0,a)") nprim, " primitive shells"
    if (allocated(ex)) deallocate(ex,cs,cp)
    allocate(ex(nprim), cs(nprim), cp(nprim), stat=aok)
    if (aok>0) then
      call die("Can't allocate arrays for primitives.")
    end if
    ex=0d0; cs=0d0; cp=0d0
  case("Atomic numbers")
    call read_line(eof)
    do i=1,nat
      call geti(k)
      name(i)=element(k)
    end do
    if (verbose) write(outfile, "(a,20a3/(17x,20a3))")                    &
        "Atoms:            ", name(1:nat)
  case("Nuclear charges")
    call read_line(eof)
    do i=1,nat
      call getf(zan(i))
    end do
    if (verbose) write(outfile, "(a,20i3/(16x,20i3))")                    &
        "Nuclear charges:", nint(zan(1:nat))
  case("Current cartesian coordinates")
    call read_line(eof)
    if (verbose) write(outfile, "(a)") "Atom  Z   Position (a.u.)"
    do i=1,nat
      do j=1,3
        call getf(c(j,i),rfact)
      end do
      if (verbose) write(outfile, "(a3,i4,3f10.5)") name(i), nint(zan(i)), c(:,i)
    end do
  case("Shell types")
    call read_line(eof)
    ! num is the original number of basis functions. n is the
    ! increased number when conversion from spherical to
    ! cartesian is taken into account.
    num=0
    n=0
    do i=1,nshell
      call geti(shell_type(i))
      select case(shell_type(i))
      case(0)
        num=num+1; n=n+1
      case(1)
        num=num+3; n=n+3
      case(-1)
        num=num+4; n=n+4
      case(2)
        num=num+6; n=n+6
      case(-2)
        num=num+5; n=n+6
      case(3)
        num=num+10; n=n+10
      case(-3)
        num=num+7; n=n+10
      case(4)
        num=num+15; n=n+15
      case(-4)
        num=num+9; n=n+15
      end select
    end do
    if (verbose) write(outfile, "(i0,a)") num, " basis functions"
    if ((verbose) .and. n>num)                           &
        write(outfile, "(a,i0,a)") "(", n, " after conversion to cartesian)"
    maxbfn=n
    if (allocated(iax)) deallocate(iax)
    allocate(iax(n+1), stat=aok)
    if (aok>0) then
      call die("Can't allocate IAX array")
    end if
  case("Number of primitives per shell")
    call read_line(eof)
    k=1
    do i=1,nshell
      call geti(j)
      kstart(i)=k
      k=k+j
      kng(i)=j
    end do
    if (verbose) then
      write(outfile,"(a,20i3/(19x,20i3))") "Contraction depths:", kng(1:nshell)
      write(outfile,"(a,i0)") "Total number of primitives required: ", k-1
    end if
    if (k .ne. nprim+1) call die                              &
        ("Shell contractions do not match number of primitives",.false.)
  case("Shell to atom map")
    call read_line(eof)
    do i=1,nshell
      call geti(katom(i))
    end do
    if (verbose) then
      write(outfile,"(a,120i3)") "shell to atom", katom(1:nshell)
    endif
  case("Primitive exponents")
    call read_line(eof)
    do i=1,nprim
      call getf(ex(i))
    end do
!    print "(a,20F10.6)", "primitive exps", ex(1:nprim)
  case("Contraction coefficients")
    call read_line(eof)
    do i=1,nshell
      do j=kstart(i),kstart(i)+kng(i)-1
        call getf(e)
        if (shell_type(i) .eq. 1) then
          cp(j)=e
        else
          cs(j)=e
        end if
        ! select case(shell_type(i))
        ! case(-1,0)
        !   cs(j)=e
        ! case(1)
        !   cp(j)=e
        ! case(2,-2,3,-3,4,-4,5,-5)
        !   cd(j)=e
        ! end select
      end do
    end do
  case("P(S=P) Contraction coefficients")
    call read_line(eof)
    do i=1,nshell
      do j=kstart(i),kstart(i)+kng(i)-1
        call getf(e)
        if (shell_type(i) .eq. -1) then
          cp(j)=e
        end if
      end do
    end do
!   case("Alpha MO coefficients","Beta MO coefficients")
!     if (debug(2)) then
!       print "(/a)", text
!       allocate(temp(ncoorb,ncoorb))
!       do i=1,ncoorb
!         do j=1,ncoorb
!           call getf(temp(j,i))
!         end do
!       end do
!       call matwrtt(temp,1,ncoorb,1,ncoorb,format="6f12.5")
!       deallocate(temp)
!     end if
  case default
    if (text .eq. density_header) then
      ncoorb=n
      nx=n*(n+1)/2
      allocate(densty(n,n),temp(n,n))
      call read_line(eof)
      do i=1,num
        do j=1,i
          call getf(densty(i,j)); densty(j,i)=densty(i,j)
        end do
      end do
      ok=.true.
    else
      !  Ignore this section
      if (nn .gt. 0) then
        if (type .eq. "I") then
          do i=1,(nn+5)/6
            read(ir,"(A8)") dummy
          end do
        else if (type .eq. "R") then
          do i=1,(nn+4)/5
            read(ir,"(A8)") dummy
          end do
        end if
      endif
    endif
  end select
end do

if (verbose) then
  atom=0
  do i=1,nshell
    if (katom(i) .ne. atom) then
      atom=katom(i)
      write(outfile, "(a)") name(atom)
    end if
    write(outfile, "(a,i0,3x,a)") "Shell ", i, label(shell_type(i))
    do j=kstart(i),kstart(i)+kng(i)-1
      select case (shell_type(i))
      case (-1)
        write(outfile, "(i10,f16.8, 2f14.8)") j, ex(j), cs(j), cp(j)
      case(0)
        write(outfile, "(i10,f16.8, 2f14.8)") j, ex(j), cs(j)
      case(1)
        write(outfile, "(i10,f16.8, 2f14.8)") j, ex(j), cp(j)
      case(2,-2,3,-3,4,-4)
        write(outfile, "(i10,f16.8, 2f14.8)") j, ex(j), cs(j)
      end select
    end do
  end do
end if
!  We use unnormalized primitive functions, so we transfer the
!  normalising factor to the contraction coefficients. This is
!  the factor for z^n exp(-e*r^2). General formula is
!  (4e)^(n/2).(2e/pi)^{3/4}/sqrt{(2n-1)!!}
do i=1,nshell
  do j=kstart(i),kstart(i)+kng(i)-1
    e=ex(j)
    select case(abs(shell_type(i)))
    case(0,1)
      cs(j)=cs(j)*sqrt(sqrt((2d0*e/pi)**3))
      cp(j)=cp(j)*sqrt(4d0*e*sqrt((2d0*e/pi)**3))
    case(5) ! h shell
      cs(j)=cs(j)*(4d0*e)**2*sqrt(4d0*e*sqrt((2d0*e/pi)**3)/945d0)
    case(4) ! g shell
      cs(j)=cs(j)*(4d0*e)**2*sqrt(sqrt((2d0*e/pi)**3)/105d0)
    case(3) ! f shell
      cs(j)=cs(j)*4d0*e*sqrt(4d0*e*sqrt((2d0*e/pi)**3)/15d0)
    case(2) ! d shell
      cs(j)=cs(j)*4d0*e*sqrt(sqrt((2d0*e/pi)**3)/3d0)
    end select
  end do
end do

if (.not. ok) return
!     call matwrtt(densty,1,num,1,num,format='5F10.5', cols=5)

!  Deal with shell types, transforming from spherical to cartesian
!  basis if necessary
k=0
do i=1,nshell
  kloc(i)=k+1 ! First basis function for shell i
  select case(shell_type(i))
  case(-1) ! sp shell
    kmin(i)=1
    kmax(i)=4
    ktype(i)=2
  case(0) ! s shell
    kmin(i)=1
    kmax(i)=1
    ktype(i)=1
  case(1) ! p shell
    kmin(i)=2
    kmax(i)=4
    ktype(i)=2
  case(2,-2) ! d shell
    kmin(i)=5
    kmax(i)=10
    ktype(i)=3
    if (shell_type(i) .lt. 0) then ! Spherical d shell
      temp(1:num,1:k)=densty(1:num,1:k)
      temp(1:num,k+1:k+6)=matmul(densty(1:num,k+1:k+5),td)
      temp(1:num,k+7:num+1)=densty(1:num,k+6:num)
      num=num+1
      densty(1:k,1:num)=temp(1:k,1:num)
      densty(k+1:k+6,1:num)=matmul(transpose(td),temp(k+1:k+5,1:num))
      densty(k+7:num,1:num)=temp(k+6:num-1,1:num)
    endif
  case(3,-3) ! f shell
    kmin(i)=11
    kmax(i)=20
    ktype(i)=4
    if (shell_type(i) .lt. 0) then ! Spherical f shell
      temp(1:num,1:k)=densty(1:num,1:k)
      temp(1:num,k+1:k+10)=matmul(densty(1:num,k+1:k+7),tf)
      if (i<nshell) temp(1:num,k+11:num+3)=densty(1:num,k+8:num)
      num=num+3
      densty(1:k,1:num)=temp(1:k,1:num)
      densty(k+1:k+10,1:num)=matmul(transpose(tf),temp(k+1:k+7,1:num))
      if (i<nshell) densty(k+11:num,1:num)=temp(k+8:num-3,1:num)
    endif
  case(4,-4) ! g shell
    kmin(i)=21
    kmax(i)=35
    ktype(i)=5
    ! print "(a,i0,a,i0)", "num = ", num, "  k = ", k
    if (shell_type(i) .lt. 0) then ! Spherical g shell
      temp(1:num,1:k)=densty(1:num,1:k)
      temp(1:num,k+1:k+15)=matmul(densty(1:num,k+1:k+9),tg)
      if (i<nshell) temp(1:num,k+16:num+6)=densty(1:num,k+10:num)
      num=num+6
      densty(1:k,1:num)=temp(1:k,1:num)
      densty(k+1:k+15,1:num)=matmul(transpose(tg),temp(k+1:k+9,1:num))
      if (i<nshell) densty(k+16:num,1:num)=temp(k+10:num-6,1:num)
    endif
  case default
    write (buffer,"(a,i0)") "Unrecognized or unimplemented shell type ", i
    call die(trim(buffer),.false.)
  end select
  k=k+kmax(i)-kmin(i)+1
end do
if (k .ne. n .or. num .ne. n) call die                             &
    ("Mismatch in number of basis functions",.false.)
if (debug(1)) call matwrtt(densty,1,n,1,n,format='5F10.5', iformat="5i10", cols=5)

deallocate(temp)

END SUBROUTINE get_data

!----------------------------------------------------------------

SUBROUTINE matwrtt(c,i1,i2,j1,j2,format,cols,iformat)

!  Print rows I1 to I2, columns J1 to J2, of the lower triangle of
!  the symmetric matrix C

IMPLICIT NONE
REAL(dp) :: c(:,:)
INTEGER i1, i2, j1, j2
CHARACTER(LEN=*), OPTIONAL :: format, iformat
INTEGER, OPTIONAL :: cols

INTEGER i, j, jstart, jfinis, ncols
CHARACTER(LEN=20) :: fmt, ifmt

fmt="1p6g12.4"
ifmt="12i12"
if (present(format)) fmt=format
if (present(iformat)) ifmt=iformat
ncols=6
if (present(cols)) ncols=cols

jfinis=j1-1
do while (jfinis .lt. j2)
  jstart=jfinis+1
  jfinis=min(j2,jfinis+ncols)
  write (6,"(/1x,"//ifmt//")") (j,j=jstart,jfinis)
  write (6,'(1x)')
  do i=max(i1,jstart),i2
    write (6,"(1x,i3,1x,"//fmt//")")                            &
        i,(c(i,j),j=jstart,min(i,jfinis))
  end do
end do

END SUBROUTINE matwrtt

END module gdma
