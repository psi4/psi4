MODULE timing

!  Timing routines
!
!  Copyright (C) 2005 Anthony J. Stone
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
!  Fifth Floor, Boston, MA 02110-1301, USA, or see
!  http://www.gnu.org/copyleft/gpl.html

REAL, SAVE ::  start_time, last_time, now

PRIVATE
PUBLIC start_timer, timer, time_and_date, start_time, now, last_time
#if defined(RS6000) || defined(ETIME) || defined(DTIME)
PUBLIC cpu_time
#endif

CONTAINS

SUBROUTINE start_timer
IMPLICIT NONE

call cpu_time(start_time)
last_time=start_time

END SUBROUTINE start_timer

!=======================================================================

SUBROUTINE timer
IMPLICIT NONE
!     CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: string

call cpu_time(now)
!     if (present(string)) then
!       print "(A,I6,A,F6.3,A)", string,                                &
!           floor((now-last_time)/60d0), "m", mod(now-last_time,60.0), "s"
!     else
print "(/A,I6,A,F6.3,2A,I6,A,F6.3,A)",                          &
    "CPU time used: ",                                          &
    floor((now-last_time)/60d0), "m", mod(now-last_time,60.0), "s", &
    "          Total: ",                                        &
    floor((now-start_time)/60d0), "m", mod(now-start_time,60.0), "s"
!     endif
last_time=now

END SUBROUTINE timer

!=======================================================================

!  cpu_time is an intrinsic routine in F95
#ifdef RS6000
SUBROUTINE cpu_time(seconds)
IMPLICIT NONE
REAL :: seconds
seconds=0.01*mclock()
END SUBROUTINE cpu_time
#elif defined(ETIME)
SUBROUTINE cpu_time(seconds)
IMPLICIT NONE
REAL :: time(2), total, seconds
REAL, EXTERNAL :: etime

total=etime(time)
seconds=time(1)

END SUBROUTINE cpu_time
#elif defined(DTIME)
SUBROUTINE cpu_time(seconds)
IMPLICIT NONE
REAL :: time(2), total, seconds
REAL, EXTERNAL :: dtime

total=dtime(time)
seconds=time(1)

END SUBROUTINE cpu_time
#endif

!=======================================================================

SUBROUTINE time_and_date(string)

IMPLICIT NONE
CHARACTER(LEN=*), INTENT(OUT) :: string
CHARACTER(LEN=3) :: month(12) = (/ "Jan", "Feb", "Mar", "Apr",    &
    "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" /)

INTEGER :: v(8)

call date_and_time(values=v)
string = trim(stri(v(5))) // ":" // trim(stri(v(6))) // ":" //    &
    trim(stri(v(7))) // " on " // trim(stri(v(3))) // " " //      &
    month(v(2)) // " " // trim(stri(v(1)))

CONTAINS

FUNCTION stri(i)

CHARACTER(LEN=8) :: stri
INTEGER, INTENT(IN) :: i

write (stri,"(i8)") i
stri=adjustl(stri)

END FUNCTION stri

END SUBROUTINE time_and_date

END MODULE timing
