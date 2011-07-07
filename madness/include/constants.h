/*
  This file is part of MADNESS.
  
  Copyright (C) 2007,2010 Oak Ridge National Laboratory
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
  
  For more information please contact:
  
  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367
  
  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
  
  $Id$
*/
#ifndef MADNESS_CONSTANTS_H
#define MADNESS_CONSTANTS_H

/// \file constants.h
/// \brief Defines common mathematical and physical constants


namespace madness {
    namespace constants {

        // Mathematical constants
        const double pi = 3.14159265358979323846264338328; ///< Mathematical constant pi


        // Misc physical constants and units
        const double calorie_joule_relationship = 4.184000000 ; // J  ... i.e.,  1 kcal/mol = 4.184 kJ/mol

        const double Debye = 3.335640035e-30 ; // Cm  (coulomb metre)


        // The following physical constants and units were obtained from NIST http://physics.nist.gov/constants
        // on 2/2/2008.  The comments contain the uncertainty and units.

        const double atomic_mass_constant = 1.660538782e-27 ; //    0.000000083e-27    kg

        const double atomic_unit_of_1st_hyperpolarizablity = 3.206361533e-53 ; //   0.000000081e-53    C^3 m^3 J^-2

        const double atomic_unit_of_2nd_hyperpolarizablity = 6.23538095e-65  ; //    0.00000031e-65     C^4 m^4 J^-3

        const double atomic_unit_of_action = 1.054571628e-34 ; //   0.000000053e-34    J s

        const double atomic_unit_of_charge = 1.602176487e-19 ; //   0.000000040e-19    C

        const double atomic_unit_of_charge_density = 1.081202300e12 ; //     0.000000027 e12     C m^-3

        const double atomic_unit_of_current = 6.62361763e-3 ; //     0.00000017e-3      A

        const double atomic_unit_of_electric_dipole_moment = 8.47835281e-30 ; //    0.00000021e-30     C m

        const double atomic_unit_of_electric_quadrupole_moment = 4.48655107e-40 ; //    0.00000011e-40     C m^2

        const double atomic_unit_of_electric_field = 5.14220632e11 ; //      0.00000013 e11      V m^-1

        const double atomic_unit_of_electric_field_gradient = 9.71736166e21 ; //  0.00000024 e21      V m^-2

        const double atomic_unit_of_electric_polarizablity = 1.6487772536e-41 ; //   0.0000000034e-41   C^2 m^2 J^-1

        const double atomic_unit_of_electric_potential = 27.21138386 ; //        0.00000068          V

        const double atomic_unit_of_energy = 4.35974394e-18 ; //     0.00000022e-18     J

        const double atomic_unit_of_force = 8.23872206e-8 ; //     0.00000041e-8      N

        const double atomic_unit_of_length = 0.52917720859e-10 ; // 0.00000000036e-10 m

        const double atomic_unit_of_magnetic_dipole_moment = 1.854801830e-23 ; //   0.000000046e-23    J T^-1

        const double atomic_unit_of_magnetic_flux_density = 2.350517382e5 ; //      0.000000059 e5      T

        const double atomic_unit_of_magnetizability = 7.891036433e-29 ; //    0.000000027e-29    J T^-2

        const double atomic_unit_of_mass = 9.10938215e-31 ; //    0.00000045e-31     kg

        const double atomic_unit_of_momentum = 1.992851565e-24 ; //   0.000000099e-24    kg m s^-1

        const double atomic_unit_of_permittivity = 1.112650056e-10 ; // (exact)               F m^-1

        const double atomic_unit_of_time = 2.418884326505e-17 ; // 0.000000000016e-17 s

        const double atomic_unit_of_velocity = 2.1876912541e6 ; //  0.0000000015 e6     m s^-1

        const double Avogadro_constant =  6.02214179e23 ; //     0.00000030 e23      mol^-1

        const double Bohr_magneton = 927.400915e-26  ; //    0.000023e-26        J T^-1

        const double Bohr_radius = 0.52917720859e-10 ; // 0.00000000036e-10 m

        const double Boltzmann_constant = 1.3806504e-23 ; //       0.0000024e-23       J K^-1

        const double Compton_wavelength = 2.4263102175e-12 ; //  0.0000000033e-12   m

        const double conductance_quantum = 7.7480917004e-5 ; //  0.0000000053e-5    S

        const double electron_g_factor = -2.0023193043622 ; //  0.0000000000015    

        const double electron_gyromagnetic_ratio = 1.760859770e11 ; //   0.000000044 e11     s^-1 T^-1

        const double electron_magnetic_moment = -928.476377e-26 ; //    0.000023e-26        J T^-1

        const double electron_magnetic_moment_anomaly =  1.15965218111e-3 ; // 0.00000000074e-3  

        const double electron_magnetic_moment_to_Bohr_magneton_ratio =  -1.00115965218111 ; // 0.000000000000 74  

        const double electron_magnetic_moment_to_nuclear_magneton_ratio = -1838.28197092 ; //     0.00000080          

        const double electron_mass = 9.10938215e-31 ; //    0.00000045e-31     kg

        const double electron_proton_mass_ratio = 5.4461702177e-4 ; //   0.0000000024e-4    

        const double electron_volt = 1.602176487e-19 ; //   0.000000040e-19    J

        const double electron_volt_hartree_relationship = 3.674932540e-2 ; //     0.000000092e-2     E_h

        const double electron_volt_hertz_relationship = 2.417989454e14 ; //     0.000000060 e14     Hz

        const double electron_volt_joule_relationship = 1.602176487e-19 ; //   0.000000040e-19    J

        const double elementary_charge = 1.602176487e-19 ; //    0.000000040e-19    C

        const double Faraday_constant = 96485.3399 ; //           0.0024                C mol^-1

        const double Fermi_coupling_constant = 1.16637e-5 ; // 0.00001e-5          GeV^-2

        const double fine_structure_constant = 7.2973525376e-3 ;  //   0.0000000050e-3    

        const double hartree_electron_volt_relationship = 27.21138386 ; //        0.00000068          eV

        const double Hartree_energy = 4.35974394e-18 ; //     0.00000022e-18     J

        const double hartree_hertz_relationship = 6.579683920722e15 ; // 0.000000000044 e15 Hz

        const double hartree_inverse_meter_relationship = 2.194746313705e7 ; //  0.000000000015 e7  m^-1

        const double hartree_joule_relationship = 4.35974394e-18 ; //    0.00000022e-18     J

        const double hartree_kelvin_relationship = 3.1577465e5 ; //        0.0000055 e5         K

        const double hertz_electron_volt_relationship = 4.13566733e-15 ; //     0.00000010e-15     eV

        const double hertz_hartree_relationship = 1.519829846006e-16 ; // 0.000000000010e-16 E_h

        const double hertz_joule_relationship = 6.62606896e-34 ; //  0.00000033e-34     J

        const double hertz_kelvin_relationship = 4.7992374e-11 ; //      0.0000084e-11       K

        const double nuclear_magneton = 5.05078324e-27 ; //     0.00000013e-27     J T^-1

        const double nuclear_magneton_in_eV_per_T = 3.1524512326e-8 ; //   0.0000000045e-8    eV T^-1

        const double Planck_constant = 6.62606896e-34 ; //    0.00000033e-34     J s

        const double Planck_constant_over_2_pi = 1.054571628e-34 ; //   0.000000053e-34    J s

        const double proton_electron_mass_ratio = 1836.15267247 ; //      0.00000080          

        const double proton_g_factor = 5.585694713 ; //        0.000000046         

        const double proton_gyromagnetic_ratio = 2.675222099e8 ; //     0.000000070 e8      s^-1 T^-1

        const double proton_magnetic_moment = 1.410606662e-26 ; //   0.000000037e-26    J T^-1

        const double proton_mass = 1.672621637e-27 ; //   0.000000083e-27    kg

        const double Rydberg_constant = 10973731.568527 ; //    0.000073             m^-1

        const double speed_of_light_in_vacuum = 299792458 ; //           (exact)               m s^-1

        const double Stefan_Boltzmann_constant = 5.670400e-8 ; //      0.000040e-8         W m^-2 K^-4

        const double unified_atomic_mass_unit = 1.660538782e-27 ; //   0.000000083e-27    kg

    }

}







#endif
