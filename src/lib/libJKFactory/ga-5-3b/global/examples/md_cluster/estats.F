#if HAVE_CONFIG_H
#   include "config.fh"
#endif
c
c                                   NOTICE
c
c   This software is being made available for internal testing and
c   evaluation purposes only. This software is a pre-release test version
c   which has not yet been authenticated and cleared for publication. Adherence
c   to this notice may be necessary for the author, Battelle Memorial
c   Institute, to successfully assert copyright in and commercialize this
c   software. This software is not intended for duplication or distribution
c   to third parties without the permission of the Manager of Software
c   Products at Pacific Northwest Laboratory, Richland, Washington,  99352.
c
      subroutine estat
#include "common.fh"
c
c    This subroutine accumulates energy statistics for the simulation.
c    The averages are stored in the array anrg(). The contents of nrg()
c    and anrg() are listed below. Temporary values for the averages are
c    stored in tmpnrg() and the uncertainties are stored anrg2().
c    Corresponding entries in anrg(), tmpnrg(), and anrg2() are all indexed
c    by the same variable.
c
c                   Contents of nrg()
c
c    nrg(1): the total number of degrees of freedom in the system
c
c    nrg(2): the total number of degrees of freedom in the system
c            times Rgas/2
c
c    nrg(3): the total energy of the system
c
c    nrg(4): the total kinetic energy of the system
c
c    nrg(5): temperature (in Kelvin) of the system
c
c    nrg(6): potential energy of the system
c
c    nrg(7): current pressure (in atmospheres)
c
c    nrg(8): the total translational kinetic energy of the system
c
c    nrg(9): the instantaneous value of the Hamiltonian
c
c    nrg(10): the instantaneous value of the virial
c
c    nrg(11): the instantaneous value of the Enthalpy
c
c    nrg(12): the instantaneous coulomb energy of the system
c
c    nrg(13): the instantaneous repulsion energy of the system
c
c    nrg(14): the instantaneous dispersion energy of the system
c
c    nrg(15): the virial contribution to the pressure (in atmospheres)
c
c    nrg(16): the instantaneous torsional energy of the system
c
c    nrg(17): instantaneous energy from harmonic oscillators
c
c    nrg(18): x component of virial
c
c    nrg(19): y component of virial
c
c    nrg(20): z component of virial
c
c    nrg(21): instantaneous energy from bond angle bends
c
c    nrg(22): x component of total momentum
c
c    nrg(23): y component of total momentum
c
c    nrg(24): z component of total momentum
c
c    nrg(25): total mass
c
c                   Contents of anrg(), tmpnrg(), and anrg2()
c
c    anrg(1): the average total energy of the system
c
c    anrg(2): the average total kinetic energy of the system
c
c    anrg(3): the average temperature (in Kelvin) of the system
c
c    anrg(4): the average potential energy of the system
c
c    anrg(5): the average pressure (in atmospheres)
c
c    anrg(6): the average total translational kinetic energy of the system
c
c    anrg(7): the average value of the Hamiltonian
c
c    anrg(8): the average value of the Enthalpy
c
c    anrg(9): the average coulomb energy of the system
c
c    anrg(10): the average repulsion energy of the system
c
c    anrg(11): the average dispersion energy of the system
c
c    anrg(12): the average x dimension of the system
c
c    anrg(13): the average y dimension of the system
c
c    anrg(14): the average z dimension of the system
c
c    anrg(15): the average value of the time scaling parameter
c
c    anrg(16): the average volume of the system
c
c    anrg(17): mean square total energy
c
c    anrg(18): mean square Hamiltonian
c
c    anrg(19): mean square kinetic energy
c
c    anrg(20): mean square potential energy
c
c    anrg(21): mean square pressure (in atmospheres**2)
c
c    anrg(22): mean square x dimension
c
c    anrg(23): mean square y dimension
c
c    anrg(24): mean square z dimension
c
c    anrg(25): mean square scaling parameter
c
c    anrg(26): mean square temperature
c
c    anrg(27): mean square Enthalpy
c
c    arng(28): the average torsional energy
c
c    arng(29): the average harmonic oscillator energy
c
c    arng(30): the average bond angle bend energy
c
      nrgcnt = nrgcnt + 1
c
      anrg(1) = anrg(1) + nrg(3)
      anrg(2) = anrg(2) + nrg(4)
      anrg(3) = anrg(3) + nrg(5)
      anrg(4) = anrg(4) + nrg(6)
      anrg(5) = anrg(5) + nrg(7)
      anrg(6) = anrg(6) + nrg(8)
      anrg(7) = anrg(7) + nrg(9)
      anrg(8) = anrg(8) + nrg(11)
      anrg(9) = anrg(9) + nrg(12)
      anrg(10) = anrg(10) + nrg(13)
      anrg(11) = anrg(11) + nrg(14)
      anrg(12) = anrg(12) + xbox
      anrg(13) = anrg(13) + ybox
      anrg(14) = anrg(14) + zbox
      anrg(15) = anrg(15) + scal1
      anrg(16) = anrg(16) + vol1
      anrg(17) = anrg(17) + nrg(3)**2
      anrg(18) = anrg(18) + nrg(9)**2
      anrg(19) = anrg(19) + nrg(4)**2
      anrg(20) = anrg(20) + nrg(6)**2
      anrg(21) = anrg(21) + nrg(7)**2
      anrg(22) = anrg(22) + xbox**2
      anrg(23) = anrg(23) + ybox**2
      anrg(24) = anrg(24) + zbox**2
      anrg(25) = anrg(25) + scal1**2
      anrg(26) = anrg(26) + nrg(5)**2
      anrg(27) = anrg(27) + nrg(11)**2
      anrg(28) = anrg(28) + nrg(16)
      anrg(29) = anrg(29) + nrg(17)
      anrg(30) = anrg(30) + nrg(21)
c
      return
      end
c
      subroutine prtnrg
#include "common.fh"
c
      double precision bstat,b2stat
      integer i,ibstat
c
      if (nrgcnt.eq.0) return
c
      do 25 i = 1, MAXNRG
        anrg(i) = anrg(i) / dble(nrgcnt)
   25 continue
c
c      call ga_dgop(3,bstat,1,'+')
c      call ga_dgop(3,b2stat,1,'+')
c      call ga_igop(3,ibstat,1,'+')
c
      if (ga_nodeid().eq.0.and.l_stdio) then
c        if (ibstat.gt.0) bstat = bstat/dble(ibstat)
c        if (ibstat.gt.0) b2stat = b2stat/dble(ibstat)
c        b2stat = sqrt(abs(b2stat-bstat**2))
c
        write(6,50)
        write(6,100) anrg(1)
        write(6,200) anrg(2)
        write(6,300) anrg(3)
c        write(6,600) anrg(6)
        write(6,400) anrg(4)
c        write(6,900) anrg(9)
        write(6,1000) anrg(10)
        write(6,1100) anrg(11)
c        write(6,2800) anrg(28)
        write(6,2900) anrg(29)
        write(6,3000) anrg(30)
        write(6,500) anrg(5)
        write(6,700) anrg(7)
        write(6,800) anrg(8)
c        write(6,850) bstat
        write(6,1150)
        write(6,1200) anrg(12)
        write(6,1300) anrg(13)
        write(6,1400) anrg(14)
        write(6,1500) anrg(15)
        write(6,1600) anrg(16)
        write(6,1650)
        write(6,1700) sqrt(abs(anrg(17)-anrg(1)**2))
        write(6,1800) sqrt(abs(anrg(18)-anrg(7)**2))
        write(6,1900) sqrt(abs(anrg(19)-anrg(2)**2))
        write(6,2000) sqrt(abs(anrg(20)-anrg(4)**2))
        write(6,2600) sqrt(abs(anrg(26)-anrg(3)**2))
        write(6,2700) sqrt(abs(anrg(27)-anrg(8)**2))
        write(6,2100) sqrt(abs(anrg(21)-anrg(5)**2))
        write(6,2200) sqrt(abs(anrg(22)-anrg(12)**2))
        write(6,2300) sqrt(abs(anrg(23)-anrg(13)**2))
        write(6,2400)
     +                   sqrt(abs(anrg(24)-anrg(14)**2))
        write(6,2500) sqrt(abs(anrg(25)-anrg(15)**2))
        write(6,2550) b2stat
      endif
c
   50 format('Average Energy statistics for simulation:')
  100 format('                   <Energy>                :',f16.4)
  200 format('                   <Kinetic Energy>        :',f16.4)
  300 format('                   <Temperature>           :',f16.4)
  400 format('                   <Potential Energy>      :',f16.4)
  500 format('                   <Pressure>              :',f16.4)
  600 format('                       <T_trans>           :',f16.4)
  700 format('                   <Hamiltonian>           :',f16.4)
  800 format('                   <Enthalpy>              :',f16.4)
  850 format('                   <Bond Length>           :',f16.4)
  900 format('                       <Coulomb Energy>    :',f16.4)
 1000 format('                       <Repulsion Energy>  :',f16.4)
 1100 format('                       <Dispersion Energy> :',f16.4)
 1150 format('Average system parameters for simulation :')
 1200 format('                   <X>                     :',f16.4)
 1300 format('                   <Y>                     :',f16.4)
 1400 format('                   <Z>                     :',f16.4)
 1500 format('                   <S>                     :',f16.4)
 1600 format('                   <V>                     :',f16.4)
 1650 format('Root mean square fluctuations:')
 1700 format('                   Energy                  :',f16.4)
 1800 format('                   Hamiltonian             :',f16.4)
 1900 format('                   Kinetic Energy          :',f16.4)
 2000 format('                   Potential Energy        :',f16.4)
 2100 format('                   Pressure                :',f16.4)
 2200 format('                   X                       :',f16.4)
 2300 format('                   Y                       :',f16.4)
 2400 format('                   Z                       :',f16.4)
 2500 format('                   S                       :',f16.4)
 2550 format('                   Bond Length             :',f16.4)
 2600 format('                   Temperature             :',f16.4)
 2700 format('                   Enthalpy                :',f16.4)
 2800 format('                       <Torsional Energy>  :',f16.4)
 2900 format('                       <Vibrational Energy>:',f16.4)
 3000 format('                       <Angle Bend Energy> :',f16.4)
      return
      end
c
      subroutine nrgsum
#include "common.fh"
c
      double precision rbuf(50)
c
c   this subroutine does a global add on some of the quantities stored
c   in the array nrg()
c
      rbuf(1) = nrg(3)
      rbuf(2) = nrg(4)
      rbuf(3) = nrg(6)
      rbuf(4) = nrg(8)
      rbuf(5) = nrg(9)
      rbuf(6) = nrg(10)
      rbuf(7) = nrg(11)
      rbuf(8) = nrg(12)
      rbuf(9) = nrg(13)
      rbuf(10) = nrg(14)
      rbuf(11) = nrg(17)
      rbuf(12) = nrg(21)
      rbuf(13) = nrg(22)
      rbuf(14) = nrg(23)
      rbuf(15) = nrg(24)
      rbuf(16) = nrg(25)
      rbuf(17) = nrg(15)
c
      call ga_dgop(3,rbuf,20,'+')
c
      nrg(3) = rbuf(1)
      nrg(4) = rbuf(2)
      nrg(6) = rbuf(3)
      nrg(8) = rbuf(4)
      nrg(9) = rbuf(5)
      nrg(10) = rbuf(6)
      nrg(11) = rbuf(7)
      nrg(12) = rbuf(8)
      nrg(13) = rbuf(9)
      nrg(14) = rbuf(10)
      nrg(17) = rbuf(11)
      nrg(21) = rbuf(12)
      nrg(22) = rbuf(13)
      nrg(23) = rbuf(14)
      nrg(24) = rbuf(15)
      nrg(25) = rbuf(16)
      nrg(15) = rbuf(17)
c
      return
      end
