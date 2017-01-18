#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

"""
| Database (Truhlar) of several classes of noncovalent interactions.
| Geometries from Truhlar and coworkers at site http://comp.chem.umn.edu/database_noncov/noncovalent.htm
| Reference energies from Truhlar and coworkers at site http://comp.chem.umn.edu/database_noncov/noncovalent.htm
| First comprehensive citation JPCA 109 5656 (2005). 

- **cp**  ``'off'`` || ``'on'``

- **rlxd** ``'off'`` || ``'on'``

- **benchmark**

  - ``'<benchmark_name>'`` <Reference>.
  - |dl| ``'<default_benchmark_name>'`` |dr| <Reference>.

- **subset**

  - ``'small'`` 3: HF-HF, He-Ne, HCCH-HCCH
  - ``'large'`` 1: BzBz_PD
  - ``'HB6'`` hydrogen-bonded
  - ``'CT7'`` charge-transfer
  - ``'DI6'`` dipole-interacting
  - ``'WI7'`` weakly interacting
  - ``'PPS5'`` pi-pi stacking

"""
import qcdb

# <<< NCB31 Database Module >>>
dbse = 'NCB31'

# <<< Database Members >>>
HRXN_SM = ['HB6-2', 'WI7-1', 'PPS5-1']
HRXN_LG = ['PPS5-5']
HB6 = ['HB6-1', 'HB6-2', 'HB6-3', 'HB6-4', 'HB6-5', 'HB6-6']
CT7 = ['CT7-1', 'CT7-2', 'CT7-3', 'CT7-4', 'CT7-5', 'CT7-6', 'CT7-7']
DI6 = ['DI6-1', 'DI6-2', 'DI6-3', 'DI6-4', 'DI6-5', 'DI6-6']
WI7 = ['WI7-1', 'WI7-2', 'WI7-3', 'WI7-4', 'WI7-5', 'WI7-6', 'WI7-7']
PPS5 = ['PPS5-1', 'PPS5-2', 'PPS5-3', 'PPS5-4', 'PPS5-5']
HRXN = sum([HB6, CT7, DI6, WI7, PPS5], [])

# <<< Chemical Systems Involved >>>
RXNM = {}        # reaction matrix of reagent contributions per reaction
RXNM_CPRLX = {}  # reaction matrix of reagent contributions per reaction for counterpoise- and deformation-corrected
ACTV = {}        # order of active reagents per reaction
ACTV_CP = {}     # order of active reagents per counterpoise-corrected reaction
ACTV_SA = {}     # order of active reagents for non-supramolecular calculations
ACTV_RLX = {}    # order of active reagents for deformation-corrected reaction
ACTV_CPRLX = {}  # order of active reagents for counterpoise- and deformation-corrected reaction

hold = {}
hold['CT7-1'] = ['C2H4', 'F2']
hold['CT7-2'] = ['NH3', 'F2']
hold['CT7-3'] = ['HCCH', 'ClF']
hold['CT7-4'] = ['HCN', 'ClF']
hold['CT7-5'] = ['NH3', 'Cl2']
hold['CT7-6'] = ['H2O', 'ClF']
hold['CT7-7'] = ['NH3', 'ClF']
hold['DI6-1'] = ['H2S', 'H2S']
hold['DI6-2'] = ['HCl', 'HCl']
hold['DI6-3'] = ['HCl', 'H2S']
hold['DI6-4'] = ['CH3Cl', 'HCl']
hold['DI6-5'] = ['HCN', 'CH3SH']
hold['DI6-6'] = ['CH3SH', 'HCl']
hold['HB6-1'] = ['NH3', 'NH3']
hold['HB6-2'] = ['HF', 'HF']
hold['HB6-3'] = ['H2O', 'H2O']
hold['HB6-4'] = ['NH3', 'H2O']
hold['HB6-5'] = ['HCONH2', 'HCONH2']
hold['HB6-6'] = ['HCOOH', 'HCOOH']
hold['PPS5-1'] = ['HCCH', 'HCCH']
hold['PPS5-2'] = ['C2H4', 'C2H4']
hold['PPS5-3'] = ['Bz', 'Bz']
hold['PPS5-4'] = ['Bz', 'Bz']
hold['PPS5-5'] = ['Bz', 'Bz']
hold['WI7-1'] = ['He', 'Ne']
hold['WI7-2'] = ['He', 'Ar']
hold['WI7-3'] = ['Ne', 'Ne']
hold['WI7-4'] = ['Ne', 'Ar']
hold['WI7-5'] = ['CH4', 'Ne']
hold['WI7-6'] = ['Bz', 'Ne']
hold['WI7-7'] = ['CH4', 'CH4']

for rxn in HRXN:
    RXNM[      '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                         '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                         '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                         '%s-%s-monoA-unCP' % (dbse, rxn) : -1,
                                         '%s-%s-monoB-unCP' % (dbse, rxn) : -1,
                                         '%s-%s-mono-RLX'   % (dbse, hold[rxn][0]) : -1,
                                         '%s-%s-mono-RLX'   % (dbse, hold[rxn][1]) : -1 }

    RXNM_CPRLX['%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                         '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                         '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                         '%s-%s-monoA-unCP' % (dbse, rxn) : +1,
                                         '%s-%s-monoB-unCP' % (dbse, rxn) : +1,
                                         '%s-%s-mono-RLX'   % (dbse, hold[rxn][0]) : -1,
                                         '%s-%s-mono-RLX'   % (dbse, hold[rxn][1]) : -1 }

    ACTV_SA[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

    ACTV[      '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                         '%s-%s-monoA-unCP' % (dbse, rxn),
                                         '%s-%s-monoB-unCP' % (dbse, rxn) ]

    ACTV_CP[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                         '%s-%s-monoA-CP'   % (dbse, rxn),
                                         '%s-%s-monoB-CP'   % (dbse, rxn) ]

    ACTV_RLX[  '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                         '%s-%s-mono-RLX'   % (dbse, hold[rxn][0]),
                                         '%s-%s-mono-RLX'   % (dbse, hold[rxn][1]) ]

    ACTV_CPRLX['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                         '%s-%s-monoA-CP'   % (dbse, rxn),
                                         '%s-%s-monoB-CP'   % (dbse, rxn),
                                         '%s-%s-monoA-unCP' % (dbse, rxn),
                                         '%s-%s-monoB-unCP' % (dbse, rxn),
                                         '%s-%s-mono-RLX'   % (dbse, hold[rxn][0]),
                                         '%s-%s-mono-RLX'   % (dbse, hold[rxn][1]) ]

# <<< Reference Values [kcal/mol] >>>
BIND = {}
nan = float('NaN')
BIND['%s-%s'            % (dbse, 'CT7-1'                 )] =      -1.06
BIND['%s-%s'            % (dbse, 'CT7-2'                 )] =      -1.81
BIND['%s-%s'            % (dbse, 'CT7-3'                 )] =      -3.81
BIND['%s-%s'            % (dbse, 'CT7-4'                 )] =      -4.86
BIND['%s-%s'            % (dbse, 'CT7-5'                 )] =      -4.88
BIND['%s-%s'            % (dbse, 'CT7-6'                 )] =      -5.36
BIND['%s-%s'            % (dbse, 'CT7-7'                 )] =     -10.62
BIND['%s-%s'            % (dbse, 'DI6-1'                 )] =      -1.66
BIND['%s-%s'            % (dbse, 'DI6-2'                 )] =      -2.01
BIND['%s-%s'            % (dbse, 'DI6-3'                 )] =      -3.35
BIND['%s-%s'            % (dbse, 'DI6-4'                 )] =      -3.55
BIND['%s-%s'            % (dbse, 'DI6-5'                 )] =      -3.59
BIND['%s-%s'            % (dbse, 'DI6-6'                 )] =      -4.16
BIND['%s-%s'            % (dbse, 'HB6-1'                 )] =      -3.15
BIND['%s-%s'            % (dbse, 'HB6-2'                 )] =      -4.57
BIND['%s-%s'            % (dbse, 'HB6-3'                 )] =      -4.97
BIND['%s-%s'            % (dbse, 'HB6-4'                 )] =      -6.41
BIND['%s-%s'            % (dbse, 'HB6-5'                 )] =     -14.94
BIND['%s-%s'            % (dbse, 'HB6-6'                 )] =     -16.15
BIND['%s-%s'            % (dbse, 'PPS5-1'                )] =      -1.34
BIND['%s-%s'            % (dbse, 'PPS5-2'                )] =      -1.42
BIND['%s-%s'            % (dbse, 'PPS5-3'                )] =      -1.81
BIND['%s-%s'            % (dbse, 'PPS5-4'                )] =      -2.74
BIND['%s-%s'            % (dbse, 'PPS5-5'                )] =      -2.78
BIND['%s-%s'            % (dbse, 'WI7-1'                 )] =      -0.04
BIND['%s-%s'            % (dbse, 'WI7-2'                 )] =      -0.06
BIND['%s-%s'            % (dbse, 'WI7-3'                 )] =      -0.08
BIND['%s-%s'            % (dbse, 'WI7-4'                 )] =      -0.13
BIND['%s-%s'            % (dbse, 'WI7-5'                 )] =      -0.22
BIND['%s-%s'            % (dbse, 'WI7-6'                 )] =      -0.47
BIND['%s-%s'            % (dbse, 'WI7-7'                 )] =      -0.51

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, 'CT7-1'                 )] = """Ethene-Fluorine Molecule Complex (C2H4-F2) """
TAGL['%s-%s-dimer'      % (dbse, 'CT7-1'                 )] = """Dimer from Ethene-Fluorine Molecule Complex (C2H4-F2) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'CT7-1'                 )] = """Monomer A from Ethene-Fluorine Molecule Complex (C2H4-F2) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'CT7-1'                 )] = """Monomer B from Ethene-Fluorine Molecule Complex (C2H4-F2) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'CT7-1'                 )] = """Monomer A from Ethene-Fluorine Molecule Complex (C2H4-F2) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'CT7-1'                 )] = """Monomer B from Ethene-Fluorine Molecule Complex (C2H4-F2) """
TAGL['%s-%s'            % (dbse, 'CT7-2'                 )] = """Ammonia-Fluorine Molecule Complex (NH3-F2) """
TAGL['%s-%s-dimer'      % (dbse, 'CT7-2'                 )] = """Dimer from Ammonia-Fluorine Molecule Complex (NH3-F2) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'CT7-2'                 )] = """Monomer A from Ammonia-Fluorine Molecule Complex (NH3-F2) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'CT7-2'                 )] = """Monomer B from Ammonia-Fluorine Molecule Complex (NH3-F2) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'CT7-2'                 )] = """Monomer A from Ammonia-Fluorine Molecule Complex (NH3-F2) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'CT7-2'                 )] = """Monomer B from Ammonia-Fluorine Molecule Complex (NH3-F2) """
TAGL['%s-%s'            % (dbse, 'CT7-3'                 )] = """Ethine-Chlorine Monofluoride Complex (HCCH-ClF) """
TAGL['%s-%s-dimer'      % (dbse, 'CT7-3'                 )] = """Dimer from Ethine-Chlorine Monofluoride Complex (HCCH-ClF) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'CT7-3'                 )] = """Monomer A from Ethine-Chlorine Monofluoride Complex (HCCH-ClF) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'CT7-3'                 )] = """Monomer B from Ethine-Chlorine Monofluoride Complex (HCCH-ClF) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'CT7-3'                 )] = """Monomer A from Ethine-Chlorine Monofluoride Complex (HCCH-ClF) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'CT7-3'                 )] = """Monomer B from Ethine-Chlorine Monofluoride Complex (HCCH-ClF) """
TAGL['%s-%s'            % (dbse, 'CT7-4'                 )] = """Hydrogen Cyanide-Chlorine Monofluoride Complex (HCN-ClF) """
TAGL['%s-%s-dimer'      % (dbse, 'CT7-4'                 )] = """Dimer from Hydrogen Cyanide-Chlorine Monofluoride Complex (HCN-ClF) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'CT7-4'                 )] = """Monomer A from Hydrogen Cyanide-Chlorine Monofluoride Complex (HCN-ClF) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'CT7-4'                 )] = """Monomer B from Hydrogen Cyanide-Chlorine Monofluoride Complex (HCN-ClF) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'CT7-4'                 )] = """Monomer A from Hydrogen Cyanide-Chlorine Monofluoride Complex (HCN-ClF) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'CT7-4'                 )] = """Monomer B from Hydrogen Cyanide-Chlorine Monofluoride Complex (HCN-ClF) """
TAGL['%s-%s'            % (dbse, 'CT7-5'                 )] = """Ammonia-Chlorine Molecule (NH3-Cl2) """
TAGL['%s-%s-dimer'      % (dbse, 'CT7-5'                 )] = """Dimer from Ammonia-Chlorine Molecule (NH3-Cl2) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'CT7-5'                 )] = """Monomer A from Ammonia-Chlorine Molecule (NH3-Cl2) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'CT7-5'                 )] = """Monomer B from Ammonia-Chlorine Molecule (NH3-Cl2) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'CT7-5'                 )] = """Monomer A from Ammonia-Chlorine Molecule (NH3-Cl2) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'CT7-5'                 )] = """Monomer B from Ammonia-Chlorine Molecule (NH3-Cl2) """
TAGL['%s-%s'            % (dbse, 'CT7-6'                 )] = """Water-Chlorine Monofluoride Complex (H2O-ClF) """
TAGL['%s-%s-dimer'      % (dbse, 'CT7-6'                 )] = """Dimer from Water-Chlorine Monofluoride Complex (H2O-ClF) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'CT7-6'                 )] = """Monomer A from Water-Chlorine Monofluoride Complex (H2O-ClF) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'CT7-6'                 )] = """Monomer B from Water-Chlorine Monofluoride Complex (H2O-ClF) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'CT7-6'                 )] = """Monomer A from Water-Chlorine Monofluoride Complex (H2O-ClF) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'CT7-6'                 )] = """Monomer B from Water-Chlorine Monofluoride Complex (H2O-ClF) """
TAGL['%s-%s'            % (dbse, 'CT7-7'                 )] = """Ammonia-Chlorine Monofluoride Complex (NH3-ClF) """
TAGL['%s-%s-dimer'      % (dbse, 'CT7-7'                 )] = """Dimer from Ammonia-Chlorine Monofluoride Complex (NH3-ClF) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'CT7-7'                 )] = """Monomer A from Ammonia-Chlorine Monofluoride Complex (NH3-ClF) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'CT7-7'                 )] = """Monomer B from Ammonia-Chlorine Monofluoride Complex (NH3-ClF) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'CT7-7'                 )] = """Monomer A from Ammonia-Chlorine Monofluoride Complex (NH3-ClF) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'CT7-7'                 )] = """Monomer B from Ammonia-Chlorine Monofluoride Complex (NH3-ClF) """
TAGL['%s-%s'            % (dbse, 'DI6-1'                 )] = """Hydrogen Sulfide Dimer (H2S-H2S) """
TAGL['%s-%s-dimer'      % (dbse, 'DI6-1'                 )] = """Dimer from Hydrogen Sulfide Dimer (H2S-H2S) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'DI6-1'                 )] = """Monomer A from Hydrogen Sulfide Dimer (H2S-H2S) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'DI6-1'                 )] = """Monomer B from Hydrogen Sulfide Dimer (H2S-H2S) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'DI6-1'                 )] = """Monomer A from Hydrogen Sulfide Dimer (H2S-H2S) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'DI6-1'                 )] = """Monomer B from Hydrogen Sulfide Dimer (H2S-H2S) """
TAGL['%s-%s'            % (dbse, 'DI6-2'                 )] = """Hydrogen Chloride Dimer (HCl-HCl) """
TAGL['%s-%s-dimer'      % (dbse, 'DI6-2'                 )] = """Dimer from Hydrogen Chloride Dimer (HCl-HCl) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'DI6-2'                 )] = """Monomer A from Hydrogen Chloride Dimer (HCl-HCl) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'DI6-2'                 )] = """Monomer B from Hydrogen Chloride Dimer (HCl-HCl) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'DI6-2'                 )] = """Monomer A from Hydrogen Chloride Dimer (HCl-HCl) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'DI6-2'                 )] = """Monomer B from Hydrogen Chloride Dimer (HCl-HCl) """
TAGL['%s-%s'            % (dbse, 'DI6-3'                 )] = """Hydrogen Chloride-Hydrogen Sulfide Complex (HCl-H2S) """
TAGL['%s-%s-dimer'      % (dbse, 'DI6-3'                 )] = """Dimer from Hydrogen Chloride-Hydrogen Sulfide Complex (HCl-H2S) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'DI6-3'                 )] = """Monomer A from Hydrogen Chloride-Hydrogen Sulfide Complex (HCl-H2S) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'DI6-3'                 )] = """Monomer B from Hydrogen Chloride-Hydrogen Sulfide Complex (HCl-H2S) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'DI6-3'                 )] = """Monomer A from Hydrogen Chloride-Hydrogen Sulfide Complex (HCl-H2S) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'DI6-3'                 )] = """Monomer B from Hydrogen Chloride-Hydrogen Sulfide Complex (HCl-H2S) """
TAGL['%s-%s'            % (dbse, 'DI6-4'                 )] = """Methyl Chloride-Hydrogen Chloride (CH3Cl-HCl) """
TAGL['%s-%s-dimer'      % (dbse, 'DI6-4'                 )] = """Dimer from Methyl Chloride-Hydrogen Chloride (CH3Cl-HCl) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'DI6-4'                 )] = """Monomer A from Methyl Chloride-Hydrogen Chloride (CH3Cl-HCl) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'DI6-4'                 )] = """Monomer B from Methyl Chloride-Hydrogen Chloride (CH3Cl-HCl) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'DI6-4'                 )] = """Monomer A from Methyl Chloride-Hydrogen Chloride (CH3Cl-HCl) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'DI6-4'                 )] = """Monomer B from Methyl Chloride-Hydrogen Chloride (CH3Cl-HCl) """
TAGL['%s-%s'            % (dbse, 'DI6-5'                 )] = """Hydrogen Cyanide-Methanethiol (HCN-CH3SH) """
TAGL['%s-%s-dimer'      % (dbse, 'DI6-5'                 )] = """Dimer from Hydrogen Cyanide-Methanethiol (HCN-CH3SH) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'DI6-5'                 )] = """Monomer A from Hydrogen Cyanide-Methanethiol (HCN-CH3SH) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'DI6-5'                 )] = """Monomer B from Hydrogen Cyanide-Methanethiol (HCN-CH3SH) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'DI6-5'                 )] = """Monomer A from Hydrogen Cyanide-Methanethiol (HCN-CH3SH) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'DI6-5'                 )] = """Monomer B from Hydrogen Cyanide-Methanethiol (HCN-CH3SH) """
TAGL['%s-%s'            % (dbse, 'DI6-6'                 )] = """Methanethiol-Hydrogen Chloride Complex (CH3SH-HCl) """
TAGL['%s-%s-dimer'      % (dbse, 'DI6-6'                 )] = """Dimer from Methanethiol-Hydrogen Chloride Complex (CH3SH-HCl) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'DI6-6'                 )] = """Monomer A from Methanethiol-Hydrogen Chloride Complex (CH3SH-HCl) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'DI6-6'                 )] = """Monomer B from Methanethiol-Hydrogen Chloride Complex (CH3SH-HCl) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'DI6-6'                 )] = """Monomer A from Methanethiol-Hydrogen Chloride Complex (CH3SH-HCl) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'DI6-6'                 )] = """Monomer B from Methanethiol-Hydrogen Chloride Complex (CH3SH-HCl) """
TAGL['%s-%s'            % (dbse, 'HB6-1'                 )] = """Ammonia Dimer (NH3-NH3) """
TAGL['%s-%s-dimer'      % (dbse, 'HB6-1'                 )] = """Dimer from Ammonia Dimer (NH3-NH3) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'HB6-1'                 )] = """Monomer A from Ammonia Dimer (NH3-NH3) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'HB6-1'                 )] = """Monomer B from Ammonia Dimer (NH3-NH3) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'HB6-1'                 )] = """Monomer A from Ammonia Dimer (NH3-NH3) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'HB6-1'                 )] = """Monomer B from Ammonia Dimer (NH3-NH3) """
TAGL['%s-%s'            % (dbse, 'HB6-2'                 )] = """Hydrogen Fluoride Dimer (HF-HF) """
TAGL['%s-%s-dimer'      % (dbse, 'HB6-2'                 )] = """Dimer from Hydrogen Fluoride Dimer (HF-HF) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'HB6-2'                 )] = """Monomer A from Hydrogen Fluoride Dimer (HF-HF) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'HB6-2'                 )] = """Monomer B from Hydrogen Fluoride Dimer (HF-HF) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'HB6-2'                 )] = """Monomer A from Hydrogen Fluoride Dimer (HF-HF) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'HB6-2'                 )] = """Monomer B from Hydrogen Fluoride Dimer (HF-HF) """
TAGL['%s-%s'            % (dbse, 'HB6-3'                 )] = """Water Dimer (H2O-H2O) """
TAGL['%s-%s-dimer'      % (dbse, 'HB6-3'                 )] = """Dimer from Water Dimer (H2O-H2O) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'HB6-3'                 )] = """Monomer A from Water Dimer (H2O-H2O) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'HB6-3'                 )] = """Monomer B from Water Dimer (H2O-H2O) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'HB6-3'                 )] = """Monomer A from Water Dimer (H2O-H2O) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'HB6-3'                 )] = """Monomer B from Water Dimer (H2O-H2O) """
TAGL['%s-%s'            % (dbse, 'HB6-4'                 )] = """Ammonia-Water Complex (NH3-H2O) """
TAGL['%s-%s-dimer'      % (dbse, 'HB6-4'                 )] = """Dimer from Ammonia-Water Complex (NH3-H2O) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'HB6-4'                 )] = """Monomer A from Ammonia-Water Complex (NH3-H2O) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'HB6-4'                 )] = """Monomer B from Ammonia-Water Complex (NH3-H2O) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'HB6-4'                 )] = """Monomer A from Ammonia-Water Complex (NH3-H2O) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'HB6-4'                 )] = """Monomer B from Ammonia-Water Complex (NH3-H2O) """
TAGL['%s-%s'            % (dbse, 'HB6-5'                 )] = """Formamide Dimer (HCONH2-HCONH2) """
TAGL['%s-%s-dimer'      % (dbse, 'HB6-5'                 )] = """Dimer from Formamide Dimer (HCONH2-HCONH2) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'HB6-5'                 )] = """Monomer A from Formamide Dimer (HCONH2-HCONH2) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'HB6-5'                 )] = """Monomer B from Formamide Dimer (HCONH2-HCONH2) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'HB6-5'                 )] = """Monomer A from Formamide Dimer (HCONH2-HCONH2) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'HB6-5'                 )] = """Monomer B from Formamide Dimer (HCONH2-HCONH2) """
TAGL['%s-%s'            % (dbse, 'HB6-6'                 )] = """Formic Acid Dimer (HCOOH-HCOOH) """
TAGL['%s-%s-dimer'      % (dbse, 'HB6-6'                 )] = """Dimer from Formic Acid Dimer (HCOOH-HCOOH) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'HB6-6'                 )] = """Monomer A from Formic Acid Dimer (HCOOH-HCOOH) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'HB6-6'                 )] = """Monomer B from Formic Acid Dimer (HCOOH-HCOOH) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'HB6-6'                 )] = """Monomer A from Formic Acid Dimer (HCOOH-HCOOH) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'HB6-6'                 )] = """Monomer B from Formic Acid Dimer (HCOOH-HCOOH) """
TAGL['%s-%s'            % (dbse, 'PPS5-1'                )] = """Ethine Dimer (HCCH-HCCH) """
TAGL['%s-%s-dimer'      % (dbse, 'PPS5-1'                )] = """Dimer from Ethine Dimer (HCCH-HCCH) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'PPS5-1'                )] = """Monomer A from Ethine Dimer (HCCH-HCCH) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'PPS5-1'                )] = """Monomer B from Ethine Dimer (HCCH-HCCH) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'PPS5-1'                )] = """Monomer A from Ethine Dimer (HCCH-HCCH) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'PPS5-1'                )] = """Monomer B from Ethine Dimer (HCCH-HCCH) """
TAGL['%s-%s'            % (dbse, 'PPS5-2'                )] = """Ethene Dimer (C2H4-C2H4) """
TAGL['%s-%s-dimer'      % (dbse, 'PPS5-2'                )] = """Dimer from Ethene Dimer (C2H4-C2H4) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'PPS5-2'                )] = """Monomer A from Ethene Dimer (C2H4-C2H4) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'PPS5-2'                )] = """Monomer B from Ethene Dimer (C2H4-C2H4) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'PPS5-2'                )] = """Monomer A from Ethene Dimer (C2H4-C2H4) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'PPS5-2'                )] = """Monomer B from Ethene Dimer (C2H4-C2H4) """
TAGL['%s-%s'            % (dbse, 'PPS5-3'                )] = """Sandwich Benzene Dimer (BzBz_S) """
TAGL['%s-%s-dimer'      % (dbse, 'PPS5-3'                )] = """Dimer from Sandwich Benzene Dimer (BzBz_S) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'PPS5-3'                )] = """Monomer A from Sandwich Benzene Dimer (BzBz_S) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'PPS5-3'                )] = """Monomer B from Sandwich Benzene Dimer (BzBz_S) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'PPS5-3'                )] = """Monomer A from Sandwich Benzene Dimer (BzBz_S) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'PPS5-3'                )] = """Monomer B from Sandwich Benzene Dimer (BzBz_S) """
TAGL['%s-%s'            % (dbse, 'PPS5-4'                )] = """T-Shaped Benzene Dimer (BzBz_T) """
TAGL['%s-%s-dimer'      % (dbse, 'PPS5-4'                )] = """Dimer from T-Shaped Benzene Dimer (BzBz_T) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'PPS5-4'                )] = """Monomer A from T-Shaped Benzene Dimer (BzBz_T) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'PPS5-4'                )] = """Monomer B from T-Shaped Benzene Dimer (BzBz_T) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'PPS5-4'                )] = """Monomer A from T-Shaped Benzene Dimer (BzBz_T) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'PPS5-4'                )] = """Monomer B from T-Shaped Benzene Dimer (BzBz_T) """
TAGL['%s-%s'            % (dbse, 'PPS5-5'                )] = """Parallel-Displaced Benzene Dimer (BzBz_PD) """
TAGL['%s-%s-dimer'      % (dbse, 'PPS5-5'                )] = """Dimer from Parallel-Displaced Benzene Dimer (BzBz_PD) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'PPS5-5'                )] = """Monomer A from Parallel-Displaced Benzene Dimer (BzBz_PD) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'PPS5-5'                )] = """Monomer B from Parallel-Displaced Benzene Dimer (BzBz_PD) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'PPS5-5'                )] = """Monomer A from Parallel-Displaced Benzene Dimer (BzBz_PD) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'PPS5-5'                )] = """Monomer B from Parallel-Displaced Benzene Dimer (BzBz_PD) """
TAGL['%s-%s'            % (dbse, 'WI7-1'                 )] = """Helium-Neon Complex (He-Ne) """
TAGL['%s-%s-dimer'      % (dbse, 'WI7-1'                 )] = """Dimer from Helium-Neon Complex (He-Ne) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'WI7-1'                 )] = """Monomer A from Helium-Neon Complex (He-Ne) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'WI7-1'                 )] = """Monomer B from Helium-Neon Complex (He-Ne) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'WI7-1'                 )] = """Monomer A from Helium-Neon Complex (He-Ne) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'WI7-1'                 )] = """Monomer B from Helium-Neon Complex (He-Ne) """
TAGL['%s-%s'            % (dbse, 'WI7-2'                 )] = """Helium-Argon Complex (He-Ar) """
TAGL['%s-%s-dimer'      % (dbse, 'WI7-2'                 )] = """Dimer from Helium-Argon Complex (He-Ar) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'WI7-2'                 )] = """Monomer A from Helium-Argon Complex (He-Ar) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'WI7-2'                 )] = """Monomer B from Helium-Argon Complex (He-Ar) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'WI7-2'                 )] = """Monomer A from Helium-Argon Complex (He-Ar) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'WI7-2'                 )] = """Monomer B from Helium-Argon Complex (He-Ar) """
TAGL['%s-%s'            % (dbse, 'WI7-3'                 )] = """Neon Dimer (Ne-Ne) """
TAGL['%s-%s-dimer'      % (dbse, 'WI7-3'                 )] = """Dimer from Neon Dimer (Ne-Ne) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'WI7-3'                 )] = """Monomer A from Neon Dimer (Ne-Ne) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'WI7-3'                 )] = """Monomer B from Neon Dimer (Ne-Ne) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'WI7-3'                 )] = """Monomer A from Neon Dimer (Ne-Ne) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'WI7-3'                 )] = """Monomer B from Neon Dimer (Ne-Ne) """
TAGL['%s-%s'            % (dbse, 'WI7-4'                 )] = """Neon-Argon Complex (Ne-Ar) """
TAGL['%s-%s-dimer'      % (dbse, 'WI7-4'                 )] = """Dimer from Neon-Argon Complex (Ne-Ar) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'WI7-4'                 )] = """Monomer A from Neon-Argon Complex (Ne-Ar) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'WI7-4'                 )] = """Monomer B from Neon-Argon Complex (Ne-Ar) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'WI7-4'                 )] = """Monomer A from Neon-Argon Complex (Ne-Ar) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'WI7-4'                 )] = """Monomer B from Neon-Argon Complex (Ne-Ar) """
TAGL['%s-%s'            % (dbse, 'WI7-5'                 )] = """Methane-Neon Complex (CH4-Ne) """
TAGL['%s-%s-dimer'      % (dbse, 'WI7-5'                 )] = """Dimer from Methane-Neon Complex (CH4-Ne) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'WI7-5'                 )] = """Monomer A from Methane-Neon Complex (CH4-Ne) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'WI7-5'                 )] = """Monomer B from Methane-Neon Complex (CH4-Ne) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'WI7-5'                 )] = """Monomer A from Methane-Neon Complex (CH4-Ne) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'WI7-5'                 )] = """Monomer B from Methane-Neon Complex (CH4-Ne) """
TAGL['%s-%s'            % (dbse, 'WI7-6'                 )] = """Benzene-Neon Complex (Bz-Ne) """
TAGL['%s-%s-dimer'      % (dbse, 'WI7-6'                 )] = """Dimer from Benzene-Neon Complex (Bz-Ne) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'WI7-6'                 )] = """Monomer A from Benzene-Neon Complex (Bz-Ne) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'WI7-6'                 )] = """Monomer B from Benzene-Neon Complex (Bz-Ne) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'WI7-6'                 )] = """Monomer A from Benzene-Neon Complex (Bz-Ne) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'WI7-6'                 )] = """Monomer B from Benzene-Neon Complex (Bz-Ne) """
TAGL['%s-%s'            % (dbse, 'WI7-7'                 )] = """Methane Dimer (CH4-CH4) """
TAGL['%s-%s-dimer'      % (dbse, 'WI7-7'                 )] = """Dimer from Methane Dimer (CH4-CH4) """
TAGL['%s-%s-monoA-CP'   % (dbse, 'WI7-7'                 )] = """Monomer A from Methane Dimer (CH4-CH4) """
TAGL['%s-%s-monoB-CP'   % (dbse, 'WI7-7'                 )] = """Monomer B from Methane Dimer (CH4-CH4) """
TAGL['%s-%s-monoA-unCP' % (dbse, 'WI7-7'                 )] = """Monomer A from Methane Dimer (CH4-CH4) """
TAGL['%s-%s-monoB-unCP' % (dbse, 'WI7-7'                 )] = """Monomer B from Methane Dimer (CH4-CH4) """
TAGL['%s-%s-mono-RLX'   % (dbse, 'HCCH'                  )] = """Ethine Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'C2H4'                  )] = """Ethene Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'Bz'                    )] = """Benzene Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'CH3Cl'                 )] = """Methyl Chloride Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'CH3SH'                 )] = """Methanethiol Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'CH4'                   )] = """Methane Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'F2'                    )] = """Fluorine Molecule Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'H2O'                   )] = """Water Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'H2S'                   )] = """Hydrogen Sulfide Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'HCl'                   )] = """Hydrogen Chloride Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'HCN'                   )] = """Hydrogen Cyanide Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'HCONH2'                )] = """Formamide Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'HCOOH'                 )] = """Formic Acid Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'He'                    )] = """Helium Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'Ne'                    )] = """Neon Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'Ar'                    )] = """Argon Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'HF'                    )] = """Hydrogen Fluoride Relaxed Monomer """
TAGL['%s-%s-mono-RLX'   % (dbse, 'NH3'                   )] = """Ammonia Relaxed Monomer """

# <<< Geometry Specification Strings >>>
GEOS = {}

GEOS['%s-%s-%s' % (dbse, 'CT7-1', 'dimer')] = qcdb.Molecule("""
0 1
C        0.00000000    -2.19285000    -0.66839500
C       -0.00000000    -2.19286000     0.66839500
H       -0.92518700    -2.19231600    -1.23398200
H        0.92518700    -2.19232500    -1.23398300
H       -0.92518700    -2.19232000     1.23398200
H        0.92518700    -2.19231100     1.23398200
--
0 1
F        0.00000000     0.78568800     0.00000000
F        0.00000000     2.20564800     0.00000100
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'CT7-2', 'dimer')] = qcdb.Molecule("""
0 1
N        0.00000000     0.00000000    -2.14998500
H        0.00000000     0.93965200    -2.53440100
H        0.81376200    -0.46982600    -2.53440100
H       -0.81376200    -0.46982600    -2.53440100
--
0 1
F        0.00000000     0.00000000     0.54577100
F        0.00000000     0.00000000     1.97124000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'CT7-3', 'dimer')] = qcdb.Molecule("""
0 1
H        0.00000000     1.67189100    -2.21255500
C        0.00000000     0.60529300    -2.19955900
C        0.00000000    -0.60529300    -2.19955900
H        0.00000000    -1.67189100    -2.21255500
--
0 1
Cl       0.00000000    -0.00000000     0.61188000
F        0.00000000    -0.00000000     2.26865100

units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'CT7-4', 'dimer')] = qcdb.Molecule("""
0 1
N        0.00000000     0.00000000    -1.83951900
C        0.00000000     0.00000000    -2.99573100
H        0.00000000     0.00000000    -4.06502600
--
0 1
F       -0.00000000     0.00000000     2.42592000
Cl      -0.00000000     0.00000000     0.76957400
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'CT7-5', 'dimer')] = qcdb.Molecule("""
0 1
N        0.00000000     0.00000000    -2.83845100
H        0.00000000     0.94268700    -3.21538300
H        0.81639100    -0.47134300    -3.21538300
H       -0.81639100    -0.47134300    -3.21538300
--
0 1
Cl       0.00000000     0.00000000    -0.15004400
Cl       0.00000000     0.00000000     1.88623900
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'CT7-6', 'dimer')] = qcdb.Molecule("""
0 1
O        2.23981900     0.00002700    -0.08823100
H        2.60088700     0.76196300     0.37705500
H        2.60108700    -0.76172700     0.37719400
--
0 1
Cl      -0.31586800    -0.00006600    -0.01691400
F       -1.97230800     0.00007400     0.02657000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'CT7-7', 'dimer')] = qcdb.Molecule("""
0 1
N        0.00000000     0.00000000    -2.05789900
H        0.00000000     0.94960500    -2.41448800
H        0.82238200    -0.47480300    -2.41448800
H       -0.82238200    -0.47480300    -2.41448800
--
0 1
Cl       0.00000000     0.00000000     0.24385500
F        0.00000000     0.00000000     1.94480300
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'DI6-1', 'dimer')] = qcdb.Molecule("""
0 1
S       -2.03099600     0.10323300    -0.00078200
H       -1.93402000    -0.81846200     0.96967600
H       -1.94045000    -0.83661600    -0.95429900
--
0 1
S        2.07983800    -0.08511200     0.00018100
H        2.33915400     1.23101900    -0.00221400
H        0.75384800     0.13412100    -0.00353700
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'DI6-2', 'dimer')] = qcdb.Molecule("""
0 1
Cl       1.86082400    -0.06541100    -0.00006800
H        1.75394100     1.21098100     0.00034100
--
0 1
Cl      -1.92526600     0.00557100    -0.00009700
H       -0.65842700    -0.19370300     0.00247600
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'DI6-3', 'dimer')] = qcdb.Molecule("""
0 1
Cl      -1.91163600    -0.00001100     0.00349800
H       -0.62731700    -0.00005800    -0.10405100
--
0 1
S        1.84252900     0.00001300    -0.10154300
H        1.82277900    -0.96181000     0.83465000
H        1.82187700     0.96186000     0.83462200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'DI6-4', 'dimer')] = qcdb.Molecule("""
0 1
C       -1.49512800     1.12579900    -0.00000200
Cl      -1.40247600    -0.66254400     0.00013900
H       -0.48106900     1.51836100    -0.00121600
H       -2.02718100     1.43516300     0.89531200
H       -2.02924000     1.43492300    -0.89417200
--
0 1
Cl       2.13960800     0.03729800    -0.00013800
H        0.97700200    -0.51405400     0.00007200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'DI6-5', 'dimer')] = qcdb.Molecule("""
0 1
C        1.99644300     0.05718500    -0.00648300
N        2.98021800     0.65834500     0.10945000
H        1.07234100    -0.48518900    -0.10641600
--
0 1
S       -1.51439900    -0.79999400    -0.11697900
C       -1.57014400     1.01297400     0.01160700
H       -1.55457900    -1.05260000     1.20049200
H       -1.54556000     1.39238100    -1.01019600
H       -0.70866100     1.40255300     0.55309700
H       -2.49314500     1.33992300     0.48665400
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'DI6-6', 'dimer')] = qcdb.Molecule("""
0 1
C       -1.44764800     1.15564900     0.01851300
S       -1.41459500    -0.65984600    -0.08354400
H       -1.46628400     1.51681600    -1.00988000
H       -0.55297100     1.53526500     0.51001200
H       -2.34423900     1.49773300     0.53186300
H       -1.37736100    -0.89092100     1.23821400
--
0 1
Cl       2.12576600     0.02408100     0.00315600
H        0.92223800    -0.44463500    -0.09824700
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'HB6-1', 'dimer')] = qcdb.Molecule("""
0 1
N        1.57522500     0.00008500    -0.04260700
H        2.13110800     0.81394900    -0.28661400
H        1.49645000    -0.00293600     0.97025700
H        2.13172100    -0.81189200    -0.29145300
--
0 1
N       -1.68824500     0.00008300     0.10484800
H       -2.12640300    -0.81268000    -0.31731000
H       -2.12744200     0.81184200    -0.31815800
H       -0.71429700     0.00054300    -0.19240700
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'HB6-2', 'dimer')] = qcdb.Molecule("""
0 1
F        1.32373600    -0.09022600    -0.00000700
H        1.74043700     0.73339000     0.00001300
--
0 1
F       -1.45719500     0.01925700    -0.00001100
H       -0.53931000    -0.09466400     0.00014500
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'HB6-3', 'dimer')] = qcdb.Molecule("""
0 1
O        1.53175000     0.00592200    -0.12088000
H        0.57596800    -0.00524900     0.02496600
H        1.90624900    -0.03756100     0.76321800
--
0 1
O       -1.39622600    -0.00499000     0.10676600
H       -1.78937200    -0.74228300    -0.37100900
H       -1.77703700     0.77763800    -0.30426400
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'HB6-4', 'dimer')] = qcdb.Molecule("""
0 1
N       -1.39559100    -0.02156400     0.00003700
H       -1.62981100     0.96109600    -0.10622400
H       -1.86276700    -0.51254400    -0.75597400
H       -1.83354700    -0.33077000     0.86230700
--
0 1
O        1.56850100     0.10589200     0.00000500
H        0.60673600    -0.03396200    -0.00062800
H        1.94051900    -0.78000500     0.00022200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'HB6-5', 'dimer')] = qcdb.Molecule("""
0 1
O       -1.14108700     1.44521200     0.00000000
C       -0.06175400     2.03094700     0.00000000
H       -0.01368700     3.13016900     0.00000000
N        1.14108700     1.43587700     0.00000000
H        1.21768600     0.41652700     0.00000000
H        1.97144600     2.00209500     0.00000000
--
0 1
O        1.14108700    -1.44521200     0.00000000
C        0.06175400    -2.03094700     0.00000000
H        0.01368700    -3.13016900     0.00000000
N       -1.14108700    -1.43587700     0.00000000
H       -1.21768600    -0.41652700     0.00000000
H       -1.97144600    -2.00209500     0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'HB6-6', 'dimer')] = qcdb.Molecule("""
0 1
C       -0.12023400     1.91407000     0.00000000
H       -0.16729500     3.00701800     0.00000000
O       -1.12185700     1.22098200     0.00000000
O        1.12185700     1.48048900     0.00000000
H        1.12758200     0.48902400     0.00000000
--
0 1
O        1.12185700    -1.22098200     0.00000000
C        0.12023400    -1.91407000     0.00000000
O       -1.12185700    -1.48048900     0.00000000
H       -1.12758200    -0.48902400     0.00000000
H        0.16729500    -3.00701800     0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'PPS5-1', 'dimer')] = qcdb.Molecule("""
0 1
C       -0.41254600     1.67817500     0.00000000
C        0.41254600     2.56162700     0.00000000
H       -1.13202600     0.89080900     0.00000000
H        1.13465100     3.34577000     0.00000000
--
0 1
C        0.41254600    -1.67817500     0.00000000
C       -0.41254600    -2.56162700     0.00000000
H        1.13202600    -0.89080900     0.00000000
H       -1.13465100    -3.34577000     0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'PPS5-2', 'dimer')] = qcdb.Molecule("""
0 1
C        1.85776800     0.47280300     0.47242500
C        1.85776800    -0.47280300    -0.47242500
H        0.93377200     0.87468800     0.87406300
H        2.78381800     0.87170900     0.87155600
H        2.78381800    -0.87170900    -0.87155600
H        0.93377200    -0.87468800    -0.87406300
--
0 1
C       -1.85776800     0.47280300    -0.47242500
C       -1.85776800    -0.47280300     0.47242500
H       -2.78381800     0.87170900    -0.87155600
H       -0.93377200     0.87468800    -0.87406300
H       -0.93377200    -0.87468800     0.87406300
H       -2.78381800    -0.87170900     0.87155600
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'PPS5-3', 'dimer')] = qcdb.Molecule("""
0 1
C        0.00000000     1.95000000     1.39150000
H        0.00000000     1.95000000     2.47150000
C        1.20507435     1.95000000     0.69575000
H        2.14038179     1.95000000     1.23575000
C        1.20507435     1.95000000    -0.69575000
H        2.14038179     1.95000000    -1.23575000
C       -0.00000000     1.95000000    -1.39150000
H       -0.00000000     1.95000000    -2.47150000
C       -1.20507435     1.95000000    -0.69575000
H       -2.14038179     1.95000000    -1.23575000
C       -1.20507435     1.95000000     0.69575000
H       -2.14038179     1.95000000     1.23575000
--
0 1
C       -1.20507435    -1.95000000    -0.69575000
H       -2.14038179    -1.95000000    -1.23575000
C       -0.00000000    -1.95000000    -1.39150000
H       -0.00000000    -1.95000000    -2.47150000
C        1.20507435    -1.95000000    -0.69575000
H        2.14038179    -1.95000000    -1.23575000
C        1.20507435    -1.95000000     0.69575000
H        2.14038179    -1.95000000     1.23575000
C       -0.00000000    -1.95000000     1.39150000
H       -0.00000000    -1.95000000     2.47150000
C       -1.20507435    -1.95000000     0.69575000
H       -2.14038179    -1.95000000     1.23575000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'PPS5-4', 'dimer')] = qcdb.Molecule("""
0 1
C        1.39150000    -0.00000000     2.49575000
H        2.47150000    -0.00000000     2.49575000
C        0.69575000     1.20507435     2.49575000
H        1.23575000     2.14038179     2.49575000
C        0.69575000    -1.20507435     2.49575000
H        1.23575000    -2.14038179     2.49575000
C       -0.69575000     1.20507435     2.49575000
H       -1.23575000     2.14038179     2.49575000
C       -0.69575000    -1.20507435     2.49575000
H       -1.23575000    -2.14038179     2.49575000
C       -1.39150000    -0.00000000     2.49575000
H       -2.47150000    -0.00000000     2.49575000
--
0 1
C        0.00000000     0.00000000    -1.10425000
C       -0.00000000    -1.20507435    -1.80000000
H       -0.00000000    -2.14038179    -1.26000000
H        0.00000000     0.00000000    -0.02425000
C       -0.00000000    -1.20507435    -3.19150000
H       -0.00000000    -2.14038179    -3.73150000
C       -0.00000000     0.00000000    -3.88725000
H       -0.00000000     0.00000000    -4.96725000
C       -0.00000000     1.20507435    -3.19150000
H        0.00000000     2.14038179    -3.73150000
C        0.00000000     1.20507435    -1.80000000
H        0.00000000     2.14038179    -1.26000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'PPS5-5', 'dimer')] = qcdb.Molecule("""
0 1
C       -0.80000000     1.80000000     1.39150000
H       -0.80000000     1.80000000     2.47150000
C        0.40507435     1.80000000     0.69575000
H        1.34038179     1.80000000     1.23575000
C       -2.00507435     1.80000000     0.69575000
H       -2.94038179     1.80000000     1.23575000
C        0.40507435     1.80000000    -0.69575000
H        1.34038179     1.80000000    -1.23575000
C       -2.00507435     1.80000000    -0.69575000
H       -2.94038179     1.80000000    -1.23575000
C       -0.80000000     1.80000000    -1.39150000
H       -0.80000000     1.80000000    -2.47150000
--
0 1
C        0.80000000    -1.80000000    -1.39150000
C        2.00507435    -1.80000000    -0.69575000
H        2.94038179    -1.80000000    -1.23575000
H        0.80000000    -1.80000000    -2.47150000
C        2.00507435    -1.80000000     0.69575000
H        2.94038179    -1.80000000     1.23575000
C        0.80000000    -1.80000000     1.39150000
H        0.80000000    -1.80000000     2.47150000
C       -0.40507435    -1.80000000     0.69575000
H       -1.34038179    -1.80000000     1.23575000
C       -0.40507435    -1.80000000    -0.69575000
H       -1.34038179    -1.80000000    -1.23575000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'WI7-1', 'dimer')] = qcdb.Molecule("""
0 1
He       0.00000000     0.00000000     0.00000000
--
0 1
Ne       3.03100000     0.00000000     0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'WI7-2', 'dimer')] = qcdb.Molecule("""
0 1
He       0.00000000     0.00000000     0.00000000
--
0 1
Ar       3.48000000     0.00000000     0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'WI7-3', 'dimer')] = qcdb.Molecule("""
0 1
Ne       0.00000000     0.00000000     0.00000000
--
0 1
Ne       3.09100000     0.00000000     0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'WI7-4', 'dimer')] = qcdb.Molecule("""
0 1
Ne       0.00000000     0.00000000     0.00000000
--
0 1
Ar       3.48900000     0.00000000     0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'WI7-5', 'dimer')] = qcdb.Molecule("""
0 1
Ne       0.00070500    -0.03504900    -1.74260200
--
0 1
C       -0.00070500     0.03504800     1.74257700
H       -0.00115700     0.05752400     2.83186300
H       -0.02121400     1.05430800     1.35836800
H       -0.87960700    -0.50371400     1.39016200
H        0.89915700    -0.46792400     1.39016200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'WI7-6', 'dimer')] = qcdb.Molecule("""
0 1
C        0.00000000     1.39566300    -0.61935100
C       -1.20868000     0.69783100    -0.61935100
C       -1.20868000    -0.69783100    -0.61935100
C       -0.00000000    -1.39566300    -0.61935100
C        1.20868000    -0.69783100    -0.61935100
C        1.20868000     0.69783100    -0.61935100
H        0.00000000     2.48003700    -0.61754900
H       -2.14777500     1.24001800    -0.61754900
H       -2.14777500    -1.24001800    -0.61754900
H       -0.00000000    -2.48003700    -0.61754900
H        2.14777500    -1.24001800    -0.61754900
H        2.14777500     1.24001800    -0.61754900
--
0 1
Ne       0.00000000     0.00000000     2.60019400
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'WI7-7', 'dimer')] = qcdb.Molecule("""
0 1
C       -0.00000000     0.00000000     1.80727900
H       -0.00000000     1.02664300     1.44240000
H       -0.88909900    -0.51332200     1.44240000
H       -0.00000000     0.00000000     2.89684300
H        0.88909900    -0.51332200     1.44240000
--
0 1
C       -0.00000000    -0.00000000    -1.80727900
H        0.88909900     0.51332200    -1.44240000
H       -0.00000000    -0.00000000    -2.89684300
H       -0.88909900     0.51332200    -1.44240000
H       -0.00000000    -1.02664300    -1.44240000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'HCCH', 'mono-RLX')] = qcdb.Molecule("""
0 1
C        0.00000400    -0.60420400     0.00000000
C        0.00000400     0.60419800     0.00000000
H        0.00679500    -1.67012800     0.00000000
H       -0.00683900     1.67016300     0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'C2H4', 'mono-RLX')] = qcdb.Molecule("""
0 1
C        0.00000000     0.00000000     0.66807800
C        0.00000000     0.00000000    -0.66807800
H        0.00000000     0.92453300     1.23491900
H        0.00000000    -0.92453300     1.23491900
H        0.00000000     0.92453300    -1.23491900
H        0.00000000    -0.92453300    -1.23491900
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'Bz', 'mono-RLX')] = qcdb.Molecule("""
0 1
C        0.00000000     1.39567100    -0.61715800
C       -1.20868600     0.69783500    -0.61715800
C       -1.20868600    -0.69783500    -0.61715800
C        0.00000000    -1.39567100    -0.61715800
C        1.20868600    -0.69783500    -0.61715800
C        1.20868600     0.69783500    -0.61715800
H        0.00000000     2.47987600    -0.61699800
H       -2.14763600     1.23993800    -0.61699800
H       -2.14763600    -1.23993800    -0.61699800
H        0.00000000    -2.47987600    -0.61699800
H        2.14763600    -1.23993800    -0.61699800
H        2.14763600     1.23993800    -0.61699800
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'CH3Cl', 'mono-RLX')] = qcdb.Molecule("""
0 1
C        0.00000000     0.00000000    -1.12626800
Cl       0.00000000     0.00000000     0.65820600
H        0.00000000     1.03097000    -1.47059600
H        0.89284600    -0.51548500    -1.47059600
H       -0.89284600    -0.51548500    -1.47059600
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'CH3SH', 'mono-RLX')] = qcdb.Molecule("""
0 1
C       -0.04788200     1.15150600     0.00000000
S       -0.04788200    -0.66495900     0.00000000
H        1.28433700    -0.82104700     0.00000000
H       -1.09471300     1.45662100     0.00000000
H        0.43188500     1.54736900     0.89371000
H        0.43188500     1.54736900    -0.89371000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'CH4', 'mono-RLX')] = qcdb.Molecule("""
0 1
C        0.00000000     0.00000000     0.00000000
H        0.00000000    -1.08947061     0.00000000
H       -1.02716274     0.36315688     0.00000000
H        0.34238759     0.36315688     0.96841832
H        0.34238759     0.36315688    -0.96841832
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'F2', 'mono-RLX')] = qcdb.Molecule("""
0 1
F        0.00000000     0.00000000     1.41423000
F        0.00000000     0.00000000     0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'H2O', 'mono-RLX')] = qcdb.Molecule("""
0 1
O        0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000     0.96183119
H        0.00000000     0.93357861    -0.23140921 
O
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'H2S', 'mono-RLX')] = qcdb.Molecule("""
0 1
S        0.00000000     0.00000000     0.10389400
H        0.00000000     0.96116200    -0.83115300
H        0.00000000    -0.96116200    -0.83115300
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'HCl', 'mono-RLX')] = qcdb.Molecule("""
0 1
Cl       0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000     1.27907275
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'HCN', 'mono-RLX')] = qcdb.Molecule("""
0 1
C        0.00000000     0.00000000    -0.50103200
N        0.00000000     0.00000000     0.65706900
H        0.00000000     0.00000000    -1.57005300
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'HCONH2', 'mono-RLX')] = qcdb.Molecule("""
0 1
C       -0.16068500     0.38839900    -0.00053800
O       -1.19570500    -0.24639200     0.00018900
N        1.08330000    -0.15841900    -0.00029100
H       -0.13991800     1.49035000     0.00139300
H        1.18225800    -1.16041500     0.00111600
H        1.90431600     0.41973500     0.00124500
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'HCOOH', 'mono-RLX')] = qcdb.Molecule("""
0 1
C       -0.13470200     0.40125100    -0.00024900
O       -1.13426200    -0.26458200     0.00006900
O        1.11868000    -0.09107500     0.00005600
H       -0.10761700     1.49546500     0.00051300
H        1.04048400    -1.05771400    -0.00002000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'He', 'mono-RLX')] = qcdb.Molecule("""
0 1
He       0.00000000     0.00000000     0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'Ne', 'mono-RLX')] = qcdb.Molecule("""
0 1
Ne       0.00000000     0.00000000     0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'Ar', 'mono-RLX')] = qcdb.Molecule("""
0 1
Ar       0.00000000     0.00000000     0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'HF', 'mono-RLX')] = qcdb.Molecule("""
0 1
F        0.00000000     0.00000000     0.00000000
H        0.00000000     0.00000000     0.92073754
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'NH3', 'mono-RLX')] = qcdb.Molecule("""
0 1
N        0.00000000     0.00000000     0.11501300
H        0.00000000     0.93975200    -0.26836400
H        0.81385000    -0.46987600    -0.26836400
H       -0.81385000    -0.46987600    -0.26836400
units angstrom
""")

# <<< Derived Geometry Strings >>>
for rxn in HRXN:
    GEOS['%s-%s-monoA-unCP' % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(1)
    GEOS['%s-%s-monoB-unCP' % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(2)
    GEOS['%s-%s-monoA-CP'   % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(1, 2)
    GEOS['%s-%s-monoB-CP'   % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(2, 1)
