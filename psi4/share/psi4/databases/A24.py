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
| Database (Hobza) of interaction energies for bimolecular complexes.
| Geometries from <Reference>.
| Reference interaction energies from Rezac and Hobza, JCTC (in press).


- **cp**  ``'off'`` <erase this comment and after unless on is a valid option> || ``'on'``

- **rlxd** ``'off'`` <erase this comment and after unless on is valid option> || ``'on'``


- **benchmark**

  - ``'<benchmark_name>'`` <Reference>.
  - |dl| ``'<default_benchmark_name>'`` |dr| <Reference>.

- **subset**

  - ``'small'`` <members_description>
  - ``'large'`` <members_description>
  - ``'<subset>'`` <members_description>

"""
import re
import qcdb

# <<< A24 Database Module >>>
dbse = 'A24'

# <<< Database Members >>>
HRXN = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
HRXN_SM = []
HRXN_LG = []

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV_CP = {}  # order of active reagents per counterpoise-corrected reaction
ACTV_SA = {}  # order of active reagents for non-supermolecular calculations
for rxn in HRXN:

    RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,
                                      '%s-%s-monoA-CP'   % (dbse, rxn) : -1,
                                      '%s-%s-monoB-CP'   % (dbse, rxn) : -1,
                                      '%s-%s-monoA-unCP' % (dbse, rxn) : -1,
                                      '%s-%s-monoB-unCP' % (dbse, rxn) : -1 }

    ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]

    ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                      '%s-%s-monoA-CP'   % (dbse, rxn),
                                      '%s-%s-monoB-CP'   % (dbse, rxn) ]

    ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),
                                      '%s-%s-monoA-unCP' % (dbse, rxn),
                                      '%s-%s-monoB-unCP' % (dbse, rxn) ]

# <<< Reference Values [kcal/mol] from Rezac and Hobza dx.doi.org/10.1021/ct400057w >>>
BIND = {}
BIND['%s-%s'            % (dbse, 1  )] =    -6.524
BIND['%s-%s'            % (dbse, 2  )] =    -5.014
BIND['%s-%s'            % (dbse, 3  )] =    -4.749
BIND['%s-%s'            % (dbse, 4  )] =    -4.572
BIND['%s-%s'            % (dbse, 5  )] =    -3.157
BIND['%s-%s'            % (dbse, 6  )] =    -1.679
BIND['%s-%s'            % (dbse, 7  )] =    -0.779
BIND['%s-%s'            % (dbse, 8  )] =    -0.672
BIND['%s-%s'            % (dbse, 9  )] =    -4.474
BIND['%s-%s'            % (dbse, 10 )] =    -2.578
BIND['%s-%s'            % (dbse, 11 )] =    -1.629
BIND['%s-%s'            % (dbse, 12 )] =    -1.537
BIND['%s-%s'            % (dbse, 13 )] =    -1.389
BIND['%s-%s'            % (dbse, 14 )] =    -1.110
BIND['%s-%s'            % (dbse, 15 )] =    -0.514
BIND['%s-%s'            % (dbse, 16 )] =    -1.518
BIND['%s-%s'            % (dbse, 17 )] =    -0.837
BIND['%s-%s'            % (dbse, 18 )] =    -0.615
BIND['%s-%s'            % (dbse, 19 )] =    -0.538
BIND['%s-%s'            % (dbse, 20 )] =    -0.408
BIND['%s-%s'            % (dbse, 21 )] =    -0.370
BIND['%s-%s'            % (dbse, 22 )] =    0.784
BIND['%s-%s'            % (dbse, 23 )] =    0.897
BIND['%s-%s'            % (dbse, 24 )] =    1.075


# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, 1)] = """ water_ammonia_Cs """
TAGL['%s-%s-dimer'      % (dbse, 1)] = """Dimer from water_ammonia_Cs """
TAGL['%s-%s-monoA-CP'   % (dbse, 1)] = """Monomer A water_ammonia_Cs  """
TAGL['%s-%s-monoB-CP'   % (dbse, 1)] = """Monomer B water_ammonia_Cs  """
TAGL['%s-%s-monoA-unCP' % (dbse, 1)] = """Monomer A water_ammonia_Cs  """
TAGL['%s-%s-monoB-unCP' % (dbse, 1)] = """Monomer B water_ammonia_Cs  """
TAGL['%s-%s'            % (dbse, 2)] = """ water_water_Cs """
TAGL['%s-%s-dimer'      % (dbse, 2)] = """Dimer from water_water_Cs """
TAGL['%s-%s-monoA-CP'   % (dbse, 2)] = """Monomer A from water_water_Cs """
TAGL['%s-%s-monoB-CP'   % (dbse, 2)] = """Monomer B from water_water_Cs """
TAGL['%s-%s-monoA-unCP' % (dbse, 2)] = """Monomer A from water_water_Cs """
TAGL['%s-%s-monoB-unCP' % (dbse, 2)] = """Monomer B from water_water_Cs """
TAGL['%s-%s'            % (dbse, 3)] = """ HCN_HCN_Cxv """
TAGL['%s-%s-dimer'      % (dbse, 3)] = """Dimer from HCN_HCN_Cxv """
TAGL['%s-%s-monoA-CP'   % (dbse, 3)] = """Monomer A from HCN_HCN_Cxv """
TAGL['%s-%s-monoB-CP'   % (dbse, 3)] = """Monomer B from HCN_HCN_Cxv """
TAGL['%s-%s-monoA-unCP' % (dbse, 3)] = """Monomer A from HCN_HCN_Cxv """
TAGL['%s-%s-monoB-unCP' % (dbse, 3)] = """Monomer B from HCN_HCN_Cxv """
TAGL['%s-%s'            % (dbse, 4)] = """ HF_HF_Cs """
TAGL['%s-%s-dimer'      % (dbse, 4)] = """Dimer from HF_HF_Cs  """
TAGL['%s-%s-monoA-CP'   % (dbse, 4)] = """Monomer A from HF_HF_Cs """
TAGL['%s-%s-monoB-CP'   % (dbse, 4)] = """Monomer B from HF_HF_Cs """
TAGL['%s-%s-monoA-unCP' % (dbse, 4)] = """Monomer A from HF_HF_Cs """
TAGL['%s-%s-monoB-unCP' % (dbse, 4)] = """Monomer B from HF_HF_Cs """
TAGL['%s-%s'            % (dbse, 5)] = """ ammonia_ammonia_C2h """
TAGL['%s-%s-dimer'      % (dbse, 5)] = """Dimer from ammonia_ammonia_C2h """
TAGL['%s-%s-monoA-CP'   % (dbse, 5)] = """Monomer A from ammonia_ammonia_C2h """
TAGL['%s-%s-monoB-CP'   % (dbse, 5)] = """Monomer B from ammonia_ammonia_C2h """
TAGL['%s-%s-monoA-unCP' % (dbse, 5)] = """Monomer A from ammonia_ammonia_C2h """
TAGL['%s-%s-monoB-unCP' % (dbse, 5)] = """Monomer B from ammonia_ammonia_C2h """
TAGL['%s-%s'            % (dbse, 6)] = """ methane_HF_C3v """
TAGL['%s-%s-dimer'      % (dbse, 6)] = """Dimer from methane_HF_C3v """
TAGL['%s-%s-monoA-CP'   % (dbse, 6)] = """Monomer A from methane_HF_C3v """
TAGL['%s-%s-monoB-CP'   % (dbse, 6)] = """Monomer B from methane_HF_C3v """
TAGL['%s-%s-monoA-unCP' % (dbse, 6)] = """Monomer A from methane_HF_C3v """
TAGL['%s-%s-monoB-unCP' % (dbse, 6)] = """Monomer B from methane_HF_C3v """
TAGL['%s-%s'            % (dbse, 7)] = """ ammmonia_methane_C3v """
TAGL['%s-%s-dimer'      % (dbse, 7)] = """Dimer from ammmonia_methane_C3v """
TAGL['%s-%s-monoA-CP'   % (dbse, 7)] = """Monomer A from ammmonia_methane_C3v """
TAGL['%s-%s-monoB-CP'   % (dbse, 7)] = """Monomer B from ammmonia_methane_C3v """
TAGL['%s-%s-monoA-unCP' % (dbse, 7)] = """Monomer A from ammmonia_methane_C3v """
TAGL['%s-%s-monoB-unCP' % (dbse, 7)] = """Monomer B from ammmonia_methane_C3v """
TAGL['%s-%s'            % (dbse, 8)] = """ methane_water_Cs """
TAGL['%s-%s-dimer'      % (dbse, 8)] = """Dimer from methane_water_Cs """
TAGL['%s-%s-monoA-CP'   % (dbse, 8)] = """Monomer A from methane_water_Cs """
TAGL['%s-%s-monoB-CP'   % (dbse, 8)] = """Monomer B from methane_water_Cs """
TAGL['%s-%s-monoA-unCP' % (dbse, 8)] = """Monomer A from methane_water_Cs """
TAGL['%s-%s-monoB-unCP' % (dbse, 8)] = """Monomer B from methane_water_Cs """
TAGL['%s-%s'            % (dbse, 9)] = """ formaldehyde_formaldehyde_Cs """
TAGL['%s-%s-dimer'      % (dbse, 9)] = """Dimer from formaldehyde_formaldehyde_Cs """
TAGL['%s-%s-monoA-CP'   % (dbse, 9)] = """Monomer A from formaldehyde_formaldehyde_Cs """
TAGL['%s-%s-monoB-CP'   % (dbse, 9)] = """Monomer B from formaldehyde_formaldehyde_Cs """
TAGL['%s-%s-monoA-unCP' % (dbse, 9)] = """Monomer A from formaldehyde_formaldehyde_Cs """
TAGL['%s-%s-monoB-unCP' % (dbse, 9)] = """Monomer B from formaldehyde_formaldehyde_Cs """
TAGL['%s-%s'            % (dbse, 10)] = """ ethene_wat_Cs """
TAGL['%s-%s-dimer'      % (dbse, 10)] = """Dimer from ethene_wat_Cs """
TAGL['%s-%s-monoA-CP'   % (dbse, 10)] = """Monomer A from ethene_wat_Cs """
TAGL['%s-%s-monoB-CP'   % (dbse, 10)] = """Monomer B from ethene_wat_Cs """
TAGL['%s-%s-monoA-unCP' % (dbse, 10)] = """Monomer A from ethene_wat_Cs """
TAGL['%s-%s-monoB-unCP' % (dbse, 10)] = """Monomer B from ethene_wat_Cs """
TAGL['%s-%s'            % (dbse, 11)] = """ ethene_formaldehyde_Cs """
TAGL['%s-%s-dimer'      % (dbse, 11)] = """Dimer from ethene_formaldehyde_Cs """
TAGL['%s-%s-monoA-CP'   % (dbse, 11)] = """Monomer A from ethene_formaldehyde_Cs """
TAGL['%s-%s-monoB-CP'   % (dbse, 11)] = """Monomer B from ethene_formaldehyde_Cs """
TAGL['%s-%s-monoA-unCP' % (dbse, 11)] = """Monomer A from ethene_formaldehyde_Cs """
TAGL['%s-%s-monoB-unCP' % (dbse, 11)] = """Monomer B from ethene_formaldehyde_Cs """
TAGL['%s-%s'            % (dbse, 12)] = """ ethyne_ethyne_C2v """
TAGL['%s-%s-dimer'      % (dbse, 12)] = """Dimer from ethyne_ethyne_C2v """
TAGL['%s-%s-monoA-CP'   % (dbse, 12)] = """Monomer A from ethyne_ethyne_C2v """
TAGL['%s-%s-monoB-CP'   % (dbse, 12)] = """Monomer B from ethyne_ethyne_C2v """
TAGL['%s-%s-monoA-unCP' % (dbse, 12)] = """Monomer A from ethyne_ethyne_C2v """
TAGL['%s-%s-monoB-unCP' % (dbse, 12)] = """Monomer B from ethyne_ethyne_C2v """
TAGL['%s-%s'            % (dbse, 13)] = """ ethene_ammonia_Cs """
TAGL['%s-%s-dimer'      % (dbse, 13)] = """Dimer from ethene_ammonia_Cs """
TAGL['%s-%s-monoA-CP'   % (dbse, 13)] = """Monomer A from ethene_ammonia_Cs """
TAGL['%s-%s-monoB-CP'   % (dbse, 13)] = """Monomer B from ethene_ammonia_Cs """
TAGL['%s-%s-monoA-unCP' % (dbse, 13)] = """Monomer A from ethene_ammonia_Cs """
TAGL['%s-%s-monoB-unCP' % (dbse, 13)] = """Monomer B from ethene_ammonia_Cs """
TAGL['%s-%s'            % (dbse, 14)] = """ ethene_ethene_C2v """
TAGL['%s-%s-dimer'      % (dbse, 14)] = """Dimer from ethene_ethene_C2v """
TAGL['%s-%s-monoA-CP'   % (dbse, 14)] = """Monomer A from ethene_ethene_C2v """
TAGL['%s-%s-monoB-CP'   % (dbse, 14)] = """Monomer B from ethene_ethene_C2v """
TAGL['%s-%s-monoA-unCP' % (dbse, 14)] = """Monomer A from ethene_ethene_C2v """
TAGL['%s-%s-monoB-unCP' % (dbse, 14)] = """Monomer B from ethene_ethene_C2v """
TAGL['%s-%s'            % (dbse, 15)] = """ methane_ethene_Cs """
TAGL['%s-%s-dimer'      % (dbse, 15)] = """Dimer from methane_ethene_Cs """
TAGL['%s-%s-monoA-CP'   % (dbse, 15)] = """Monomer A from methane_ethene_Cs """
TAGL['%s-%s-monoB-CP'   % (dbse, 15)] = """Monomer B from methane_ethene_Cs """
TAGL['%s-%s-monoA-unCP' % (dbse, 15)] = """Monomer A from methane_ethene_Cs """
TAGL['%s-%s-monoB-unCP' % (dbse, 15)] = """Monomer B from methane_ethene_Cs """
TAGL['%s-%s'            % (dbse, 16)] = """ borane_methane_Cs """
TAGL['%s-%s-dimer'      % (dbse, 16)] = """Dimer from borane_methane_Cs """
TAGL['%s-%s-monoA-CP'   % (dbse, 16)] = """Monomer A from borane_methane_Cs """
TAGL['%s-%s-monoB-CP'   % (dbse, 16)] = """Monomer B from borane_methane_Cs """
TAGL['%s-%s-monoA-unCP' % (dbse, 16)] = """Monomer A from borane_methane_Cs """
TAGL['%s-%s-monoB-unCP' % (dbse, 16)] = """Monomer B from borane_methane_Cs """
TAGL['%s-%s'            % (dbse, 17)] = """ methane_ethane_Cs """
TAGL['%s-%s-dimer'      % (dbse, 17)] = """Dimer from methane_ethane_Cs """
TAGL['%s-%s-monoA-CP'   % (dbse, 17)] = """Monomer A from methane_ethane_Cs """
TAGL['%s-%s-monoB-CP'   % (dbse, 17)] = """Monomer B from methane_ethane_Cs """
TAGL['%s-%s-monoA-unCP' % (dbse, 17)] = """Monomer A from methane_ethane_Cs """
TAGL['%s-%s-monoB-unCP' % (dbse, 17)] = """Monomer B from methane_ethane_Cs """
TAGL['%s-%s'            % (dbse, 18)] = """ methane_ethane_C3 """
TAGL['%s-%s-dimer'      % (dbse, 18)] = """Dimer from methane_ethane_C3 """
TAGL['%s-%s-monoA-CP'   % (dbse, 18)] = """Monomer A from methane_ethane_C3 """
TAGL['%s-%s-monoB-CP'   % (dbse, 18)] = """Monomer B from methane_ethane_C3 """
TAGL['%s-%s-monoA-unCP' % (dbse, 18)] = """Monomer A from methane_ethane_C3 """
TAGL['%s-%s-monoB-unCP' % (dbse, 18)] = """Monomer B from methane_ethane_C3 """
TAGL['%s-%s'            % (dbse, 19)] = """ methane_methane_D3d """
TAGL['%s-%s-dimer'      % (dbse, 19)] = """Dimer from methane_methane_D3d """
TAGL['%s-%s-monoA-CP'   % (dbse, 19)] = """Monomer A from methane_methane_D3d """
TAGL['%s-%s-monoB-CP'   % (dbse, 19)] = """Monomer B from methane_methane_D3d """
TAGL['%s-%s-monoA-unCP' % (dbse, 19)] = """Monomer A from methane_methane_D3d """
TAGL['%s-%s-monoB-unCP' % (dbse, 19)] = """Monomer B from methane_methane_D3d """
TAGL['%s-%s'            % (dbse, 20)] = """ methane_Ar_C3v """
TAGL['%s-%s-dimer'      % (dbse, 20)] = """Dimer from methane_Ar_C3v """
TAGL['%s-%s-monoA-CP'   % (dbse, 20)] = """Monomer A from methane_Ar_C3v """
TAGL['%s-%s-monoB-CP'   % (dbse, 20)] = """Monomer B from methane_Ar_C3v """
TAGL['%s-%s-monoA-unCP' % (dbse, 20)] = """Monomer A from methane_Ar_C3v """
TAGL['%s-%s-monoB-unCP' % (dbse, 20)] = """Monomer B from methane_Ar_C3v """
TAGL['%s-%s'            % (dbse, 21)] = """ ethene_Ar_C2v """
TAGL['%s-%s-dimer'      % (dbse, 21)] = """Dimer from ethene_Ar_C2v """
TAGL['%s-%s-monoA-CP'   % (dbse, 21)] = """Monomer A from ethene_Ar_C2v """
TAGL['%s-%s-monoB-CP'   % (dbse, 21)] = """Monomer B from ethene_Ar_C2v """
TAGL['%s-%s-monoA-unCP' % (dbse, 21)] = """Monomer A from ethene_Ar_C2v """
TAGL['%s-%s-monoB-unCP' % (dbse, 21)] = """Monomer B from ethene_Ar_C2v """
TAGL['%s-%s'            % (dbse, 22)] = """ ethene_ethyne_C2v """
TAGL['%s-%s-dimer'      % (dbse, 22)] = """Dimer from ethene_ethyne_C2v """
TAGL['%s-%s-monoA-CP'   % (dbse, 22)] = """Monomer A from ethene_ethyne_C2v """
TAGL['%s-%s-monoB-CP'   % (dbse, 22)] = """Monomer B from ethene_ethyne_C2v """
TAGL['%s-%s-monoA-unCP' % (dbse, 22)] = """Monomer A from ethene_ethyne_C2v """
TAGL['%s-%s-monoB-unCP' % (dbse, 22)] = """Monomer B from ethene_ethyne_C2v """
TAGL['%s-%s'            % (dbse, 23)] = """ ethene_ethene_D2h """
TAGL['%s-%s-dimer'      % (dbse, 23)] = """Dimer from ethene_ethene_D2h """
TAGL['%s-%s-monoA-CP'   % (dbse, 23)] = """Monomer A from ethene_ethene_D2h """
TAGL['%s-%s-monoB-CP'   % (dbse, 23)] = """Monomer B from ethene_ethene_D2h """
TAGL['%s-%s-monoA-unCP' % (dbse, 23)] = """Monomer A from ethene_ethene_D2h """
TAGL['%s-%s-monoB-unCP' % (dbse, 23)] = """Monomer B from ethene_ethene_D2h """
TAGL['%s-%s'            % (dbse, 24)] = """ ethyne_ethyne_D2h """
TAGL['%s-%s-dimer'      % (dbse, 24)] = """Dimer from ethyne_ethyne_D2h """
TAGL['%s-%s-monoA-CP'   % (dbse, 24)] = """Monomer A from ethyne_ethyne_D2h """
TAGL['%s-%s-monoB-CP'   % (dbse, 24)] = """Monomer B from ethyne_ethyne_D2h """
TAGL['%s-%s-monoA-unCP' % (dbse, 24)] = """Monomer A from ethyne_ethyne_D2h """
TAGL['%s-%s-monoB-unCP' % (dbse, 24)] = """Monomer B from ethyne_ethyne_D2h """

# <<< Geometry Specification Strings >>>
GEOS = {}

GEOS['%s-%s-dimer' % (dbse, '1')] = qcdb.Molecule("""
0 1
O          0.00000000      -0.05786571      -1.47979303
H          0.00000000       0.82293384      -1.85541474
H          0.00000000       0.07949567      -0.51934253
--
0 1
N          0.00000000       0.01436394       1.46454628
H          0.00000000      -0.98104857       1.65344779
H         -0.81348351       0.39876776       1.92934049
H          0.81348351       0.39876776       1.92934049
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '2')] = qcdb.Molecule("""
0 1
O         -0.06699914       0.00000000       1.49435474
H          0.81573427       0.00000000       1.86586639
H          0.06885510       0.00000000       0.53914277
--
0 1
O          0.06254775       0.00000000      -1.42263208
H         -0.40696540      -0.76017841      -1.77174450
H         -0.40696540       0.76017841      -1.77174450
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '3')] = qcdb.Molecule("""
0 1
H          0.00000000       0.00000000       3.85521306
C          0.00000000       0.00000000       2.78649976
N          0.00000000       0.00000000       1.63150791
--
0 1
H          0.00000000       0.00000000      -0.59377492
C          0.00000000       0.00000000      -1.66809824
N          0.00000000       0.00000000      -2.82525056
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '4')] = qcdb.Molecule("""
0 1
H          0.00000000       0.80267982       1.69529329
F          0.00000000      -0.04596666       1.34034818
--
0 1
H          0.00000000      -0.12040787      -0.49082840
F          0.00000000       0.00976945      -1.40424978
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '5')] = qcdb.Molecule("""
0 1
N         -0.04998129      -1.58709323       0.00000000
H          0.12296265      -2.16846018       0.81105976
H          0.12296265      -2.16846018      -0.81105976
H          0.65988580      -0.86235298       0.00000000
--
0 1
N          0.04998129       1.58709323       0.00000000
H         -0.12296265       2.16846018       0.81105976
H         -0.65988580       0.86235298       0.00000000
H         -0.12296265       2.16846018      -0.81105976
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '6')] = qcdb.Molecule("""
0 1
C          0.00000000      -0.00000000       1.77071609
H          0.51593378      -0.89362352       1.42025061
H         -0.00000000       0.00000000       2.85805859
H          0.51593378       0.89362352       1.42025061
H         -1.03186756       0.00000000       1.42025061
--
0 1
H         -0.00000000       0.00000000      -0.54877328
F         -0.00000000       0.00000000      -1.46803256
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '7')] = qcdb.Molecule("""
0 1
N         -0.00000000       0.00000000       1.84833659
H          0.93730979      -0.00000000       2.23206741
H         -0.46865489      -0.81173409       2.23206741
H         -0.46865489       0.81173409       2.23206741
--
0 1
H          0.00000000      -0.00000000      -0.94497174
C          0.00000000      -0.00000000      -2.03363752
H          0.51251439       0.88770096      -2.40095125
H          0.51251439      -0.88770096      -2.40095125
H         -1.02502878       0.00000000      -2.40095125
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '8')] = qcdb.Molecule("""
0 1
C          0.00069016       0.00000000      -1.99985520
H         -0.50741740       0.88759452      -2.37290605
H          1.03052749       0.00000000      -2.35282982
H         -0.01314396       0.00000000      -0.91190852
H         -0.50741740      -0.88759452      -2.37290605
--
0 1
O         -0.00472553       0.00000000       1.71597466
H          0.03211863       0.75755459       2.30172044
H          0.03211863      -0.75755459       2.30172044
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '9')] = qcdb.Molecule("""
0 1
C          0.00000000       0.60123980      -1.35383976
O          0.00000000      -0.59301814      -1.55209021
H          0.93542250       1.17427624      -1.26515132
H         -0.93542250       1.17427624      -1.26515132
--
0 1
C          0.00000000      -0.60200476       1.55228866
O          0.00000000       0.59238638       1.35511328
H          0.00000000      -1.00937982       2.57524635
H          0.00000000      -1.32002906       0.71694997
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '10')] = qcdb.Molecule("""
0 1
C          0.01058825      -0.66806246       1.29820809
C          0.01058825       0.66806246       1.29820809
H          0.86863216       1.23267933       0.95426815
H         -0.84608285       1.23258495       1.64525385
H         -0.84608285      -1.23258495       1.64525385
H          0.86863216      -1.23267933       0.95426815
--
0 1
H         -0.79685627       0.00000000      -2.50911038
O          0.04347445       0.00000000      -2.04834054
H         -0.19067546       0.00000000      -1.11576944
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '11')] = qcdb.Molecule("""
0 1
C          0.00000000      -0.59797089       1.47742864
C          0.00000000       0.42131196       2.33957848
H          0.92113351      -1.02957102       1.10653516
H         -0.92113351      -1.02957102       1.10653516
H         -0.92393815       0.85124826       2.70694633
H          0.92393815       0.85124826       2.70694633
--
0 1
O          0.00000000      -0.51877334      -1.82845679
C          0.00000000       0.68616220      -1.73709412
H          0.00000000       1.33077474      -2.63186355
H          0.00000000       1.18902807      -0.75645498
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '12')] = qcdb.Molecule("""
0 1
C          0.00000000       0.60356400      -2.18173438
H          0.00000000       1.66847581      -2.18429610
C          0.00000000      -0.60356400      -2.18173438
H          0.00000000      -1.66847581      -2.18429610
--
0 1
C         -0.00000000       0.00000000       1.57829513
H         -0.00000000       0.00000000       0.51136193
C         -0.00000000       0.00000000       2.78576543
H         -0.00000000       0.00000000       3.85017859
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '13')] = qcdb.Molecule("""
0 1
C          0.00000000      -0.59662248       1.58722206
C          0.00000000       0.68258238       1.20494642
H          0.92312147       1.22423658       1.04062463
H         -0.92312147       1.22423658       1.04062463
H         -0.92388993      -1.13738548       1.75121281
H          0.92388993      -1.13738548       1.75121281
--
0 1
N          0.00000000      -0.00401379      -2.31096701
H         -0.81122549      -0.45983060      -2.71043881
H          0.00000000      -0.22249432      -1.32128161
H          0.81122549      -0.45983060      -2.71043881
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '14')] = qcdb.Molecule("""
0 1
H          0.92444510      -1.23172221      -1.90619313
H         -0.92444510      -1.23172221      -1.90619313
H         -0.92444510       1.23172221      -1.90619313
H          0.92444510       1.23172221      -1.90619313
C          0.00000000       0.66728778      -1.90556520
C          0.00000000      -0.66728778      -1.90556520
--
0 1
H         -0.00000000       1.23344948       2.82931792
H          0.00000000       1.22547148       0.97776199
H         -0.00000000      -1.22547148       0.97776199
H         -0.00000000      -1.23344948       2.82931792
C         -0.00000000      -0.66711698       1.90601042
C         -0.00000000       0.66711698       1.90601042
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '15')] = qcdb.Molecule("""
0 1
C          0.00000000       0.64634385      -1.60849815
C          0.00000000      -0.67914355      -1.45381675
H         -0.92399961      -1.24016223      -1.38784883
H          0.92399961      -1.24016223      -1.38784883
H          0.92403607       1.20737602      -1.67357285
H         -0.92403607       1.20737602      -1.67357285
--
0 1
H          0.00000000       0.08295411       1.59016711
C          0.00000000       0.02871509       2.67711785
H          0.88825459       0.52261990       3.06664029
H         -0.88825459       0.52261990       3.06664029
H          0.00000000      -1.01394800       2.98955227
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '16')] = qcdb.Molecule("""
0 1
C          0.00346000       0.00000000       1.38045208
H          0.84849635       0.00000000       0.68958651
H          0.39513333       0.00000000       2.39584935
H         -0.60268447      -0.88994299       1.22482674
H         -0.60268447       0.88994299       1.22482674
--
0 1
B         -0.00555317       0.00000000      -1.59887976
H          0.58455128      -1.03051800      -1.67949525
H          0.58455128       1.03051800      -1.67949525
H         -1.18903148       0.00000000      -1.47677217
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '17')] = qcdb.Molecule("""
0 1
C          0.00000000      -0.06374421       2.42054090
H          0.00000000       1.02169396       2.34238038
H          0.88828307      -0.46131911       1.93307194
H         -0.88828307      -0.46131911       1.93307194
H          0.00000000      -0.35363606       3.46945195
--
0 1
C          0.00000000       0.78133572      -1.13543912
H          0.00000000       1.37465349      -2.05114442
H         -0.88043002       1.06310554      -0.55580918
C          0.00000000      -0.71332890      -1.44723686
H          0.88043002       1.06310554      -0.55580918
H          0.00000000      -1.30641812      -0.53140693
H         -0.88100343      -0.99533072      -2.02587154
H          0.88100343      -0.99533072      -2.02587154
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '18')] = qcdb.Molecule("""
0 1
C         -0.00000000       0.00000000      -2.85810471
H          0.39304720      -0.94712229      -2.49369739
H          0.62370837       0.81395000      -2.49369739
H         -1.01675556       0.13317229      -2.49369739
H          0.00000000      -0.00000000      -3.94634214
--
0 1
C          0.00000000      -0.00000000       0.76143405
C         -0.00000000      -0.00000000       2.28821715
H         -0.61711193      -0.80824397       0.36571527
H         -0.39140385       0.93855659       0.36571527
H          1.00851577      -0.13031262       0.36571527
H         -1.00891703       0.13031295       2.68258296
H          0.39160418      -0.93890425       2.68258296
H          0.61731284       0.80859130       2.68258296
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '19')] = qcdb.Molecule("""
0 1
C         -0.00000000       0.00000000       1.81901457
H          0.51274115       0.88809373       1.45476743
H          0.51274115      -0.88809373       1.45476743
H         -1.02548230       0.00000000       1.45476743
H          0.00000000      -0.00000000       2.90722072
--
0 1
C          0.00000000      -0.00000000      -1.81901457
H         -0.00000000       0.00000000      -2.90722072
H         -0.51274115       0.88809373      -1.45476743
H         -0.51274115      -0.88809373      -1.45476743
H          1.02548230      -0.00000000      -1.45476743
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '20')] = qcdb.Molecule("""
0 1
C         -0.00000000       0.00000000      -2.62458428
H          0.51286762       0.88831278      -2.26110195
H          0.51286762      -0.88831278      -2.26110195
H         -0.00000000       0.00000000      -3.71273928
H         -1.02573525       0.00000000      -2.26110195
--
0 1
AR        -0.00000000       0.00000000       1.05395172
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '21')] = qcdb.Molecule("""
0 1
C          0.00000000       0.66718073      -2.29024825
C          0.00000000      -0.66718073      -2.29024825
H         -0.92400768       1.23202333      -2.28975239
H          0.92400768       1.23202333      -2.28975239
H         -0.92400768      -1.23202333      -2.28975239
H          0.92400768      -1.23202333      -2.28975239
--
0 1
AR        -0.00000000       0.00000000       1.60829261
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '22')] = qcdb.Molecule("""
0 1
H         -0.92396100       1.23195600      -1.68478123
H          0.92396100       1.23195600      -1.68478123
H          0.92396100      -1.23195600      -1.68478123
H         -0.92396100      -1.23195600      -1.68478123
C          0.00000000       0.66717600      -1.68478123
C          0.00000000      -0.66717600      -1.68478123
--
0 1
H         -0.00000000      -1.66786500       1.81521877
H         -0.00000000       1.66786500       1.81521877
C         -0.00000000      -0.60339700       1.81521877
C         -0.00000000       0.60339700       1.81521877
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '23')] = qcdb.Molecule("""
0 1
H         -0.92396100       1.23195600      -1.75000000
H          0.92396100       1.23195600      -1.75000000
H          0.92396100      -1.23195600      -1.75000000
H         -0.92396100      -1.23195600      -1.75000000
C          0.00000000       0.66717600      -1.75000000
C         -0.00000000      -0.66717600      -1.75000000
--
0 1
H         -0.92396100       1.23195600       1.75000000
H          0.92396100       1.23195600       1.75000000
H          0.92396100      -1.23195600       1.75000000
H         -0.92396100      -1.23195600       1.75000000
C          0.00000000       0.66717600       1.75000000
C         -0.00000000      -0.66717600       1.75000000
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '24')] = qcdb.Molecule("""
0 1
H         -0.00000000      -1.66786500      -1.75000000
H          0.00000000       1.66786500      -1.75000000
C         -0.00000000      -0.60339700      -1.75000000
C          0.00000000       0.60339700      -1.75000000
--
0 1
H         -0.00000000      -1.66786500       1.75000000
H          0.00000000       1.66786500       1.75000000
C         -0.00000000      -0.60339700       1.75000000
C          0.00000000       0.60339700       1.75000000
units angstrom
""")


# <<< Derived Geometry Strings >>>
for rxn in HRXN:
    GEOS['%s-%s-monoA-unCP' % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(1)
    GEOS['%s-%s-monoB-unCP' % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(2)
    GEOS['%s-%s-monoA-CP'   % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(1, 2)
    GEOS['%s-%s-monoB-CP'   % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(2, 1)
