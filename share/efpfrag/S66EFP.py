#
#@BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
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
#@END LICENSE
#

"""
| Database (Hobza) of interaction energies for bimolecular complexes.
| Geometries and reference energies from Rezac et al. JCTC 7 2427 (2011).

- **cp**  ``'off'`` || ``'on'``

- **rlxd** ``'off'``

- **subset**

  - ``'small'``
  - ``'large'``
  - ``'HB'`` hydrogen-bonded systems
  - ``'MX'`` mixed-influence systems
  - ``'DD'`` dispersion-dominated systems

"""
import re
import qcdb

# <<< S66 Database Module >>>
dbse = 'S66'

# <<< Database Members >>>
HRXN = range(1, 67)
HRXN_SM = [1, 12, 59]
HRXN_LG = [26, 34]
HB = range(1, 24)
MX = range(47, 67)
DD = range(24, 47)

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

# <<< Reference Values [kcal/mol] >>>
BIND = {}
BIND['%s-%s'            % (dbse, '1'                     )] =   -4.918
BIND['%s-%s'            % (dbse, '2'                     )] =   -5.592
BIND['%s-%s'            % (dbse, '3'                     )] =   -6.908
BIND['%s-%s'            % (dbse, '4'                     )] =   -8.103
BIND['%s-%s'            % (dbse, '5'                     )] =   -5.757
BIND['%s-%s'            % (dbse, '6'                     )] =   -7.554
BIND['%s-%s'            % (dbse, '7'                     )] =   -8.230
BIND['%s-%s'            % (dbse, '8'                     )] =   -5.009
BIND['%s-%s'            % (dbse, '9'                     )] =   -3.059
BIND['%s-%s'            % (dbse, '10'                    )] =   -4.160
BIND['%s-%s'            % (dbse, '11'                    )] =   -5.419
BIND['%s-%s'            % (dbse, '12'                    )] =   -7.266
BIND['%s-%s'            % (dbse, '13'                    )] =   -6.187
BIND['%s-%s'            % (dbse, '14'                    )] =   -7.454
BIND['%s-%s'            % (dbse, '15'                    )] =   -8.630
BIND['%s-%s'            % (dbse, '16'                    )] =   -5.124
BIND['%s-%s'            % (dbse, '17'                    )] =  -17.182
BIND['%s-%s'            % (dbse, '18'                    )] =   -6.857
BIND['%s-%s'            % (dbse, '19'                    )] =   -7.410
BIND['%s-%s'            % (dbse, '20'                    )] =  -19.093
BIND['%s-%s'            % (dbse, '21'                    )] =  -16.265
BIND['%s-%s'            % (dbse, '22'                    )] =  -19.491
BIND['%s-%s'            % (dbse, '23'                    )] =  -19.189
BIND['%s-%s'            % (dbse, '24'                    )] =   -2.822
BIND['%s-%s'            % (dbse, '25'                    )] =   -3.895
BIND['%s-%s'            % (dbse, '26'                    )] =   -9.829
BIND['%s-%s'            % (dbse, '27'                    )] =   -3.439
BIND['%s-%s'            % (dbse, '28'                    )] =   -5.713
BIND['%s-%s'            % (dbse, '29'                    )] =   -6.819
BIND['%s-%s'            % (dbse, '30'                    )] =   -1.432
BIND['%s-%s'            % (dbse, '31'                    )] =   -3.380
BIND['%s-%s'            % (dbse, '32'                    )] =   -3.738
BIND['%s-%s'            % (dbse, '33'                    )] =   -1.872
BIND['%s-%s'            % (dbse, '34'                    )] =   -3.776
BIND['%s-%s'            % (dbse, '35'                    )] =   -2.613
BIND['%s-%s'            % (dbse, '36'                    )] =   -1.777
BIND['%s-%s'            % (dbse, '37'                    )] =   -2.404
BIND['%s-%s'            % (dbse, '38'                    )] =   -2.997
BIND['%s-%s'            % (dbse, '39'                    )] =   -3.575
BIND['%s-%s'            % (dbse, '40'                    )] =   -2.895
BIND['%s-%s'            % (dbse, '41'                    )] =   -4.848
BIND['%s-%s'            % (dbse, '42'                    )] =   -4.138
BIND['%s-%s'            % (dbse, '43'                    )] =   -3.712
BIND['%s-%s'            % (dbse, '44'                    )] =   -2.005
BIND['%s-%s'            % (dbse, '45'                    )] =   -1.748
BIND['%s-%s'            % (dbse, '46'                    )] =   -4.264
BIND['%s-%s'            % (dbse, '47'                    )] =   -2.876
BIND['%s-%s'            % (dbse, '48'                    )] =   -3.535
BIND['%s-%s'            % (dbse, '49'                    )] =   -3.331
BIND['%s-%s'            % (dbse, '50'                    )] =   -2.867
BIND['%s-%s'            % (dbse, '51'                    )] =   -1.524
BIND['%s-%s'            % (dbse, '52'                    )] =   -4.707
BIND['%s-%s'            % (dbse, '53'                    )] =   -4.361
BIND['%s-%s'            % (dbse, '54'                    )] =   -3.277
BIND['%s-%s'            % (dbse, '55'                    )] =   -4.188
BIND['%s-%s'            % (dbse, '56'                    )] =   -3.231
BIND['%s-%s'            % (dbse, '57'                    )] =   -5.282
BIND['%s-%s'            % (dbse, '58'                    )] =   -4.146
BIND['%s-%s'            % (dbse, '59'                    )] =   -2.850
BIND['%s-%s'            % (dbse, '60'                    )] =   -4.868
BIND['%s-%s'            % (dbse, '61'                    )] =   -2.912
BIND['%s-%s'            % (dbse, '62'                    )] =   -3.534
BIND['%s-%s'            % (dbse, '63'                    )] =   -3.801
BIND['%s-%s'            % (dbse, '64'                    )] =   -2.999
BIND['%s-%s'            % (dbse, '65'                    )] =   -3.991
BIND['%s-%s'            % (dbse, '66'                    )] =   -3.968

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, '1'                     )] = """Water Dimer """
TAGL['%s-%s-dimer'      % (dbse, '1'                     )] = """Dimer from Water Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '1'                     )] = """Monomer A from Water Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '1'                     )] = """Monomer B from Water Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '1'                     )] = """Monomer A from Water Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '1'                     )] = """Monomer B from Water Dimer """
TAGL['%s-%s'            % (dbse, '2'                     )] = """Water-Methanol """
TAGL['%s-%s-dimer'      % (dbse, '2'                     )] = """Dimer from Water-Methanol """
TAGL['%s-%s-monoA-CP'   % (dbse, '2'                     )] = """Monomer A from Water-Methanol """
TAGL['%s-%s-monoB-CP'   % (dbse, '2'                     )] = """Monomer B from Water-Methanol """
TAGL['%s-%s-monoA-unCP' % (dbse, '2'                     )] = """Monomer A from Water-Methanol """
TAGL['%s-%s-monoB-unCP' % (dbse, '2'                     )] = """Monomer B from Water-Methanol """
TAGL['%s-%s'            % (dbse, '3'                     )] = """Water-Methylamine """
TAGL['%s-%s-dimer'      % (dbse, '3'                     )] = """Dimer from Water-Methylamine """
TAGL['%s-%s-monoA-CP'   % (dbse, '3'                     )] = """Monomer A from Water-Methylamine """
TAGL['%s-%s-monoB-CP'   % (dbse, '3'                     )] = """Monomer B from Water-Methylamine """
TAGL['%s-%s-monoA-unCP' % (dbse, '3'                     )] = """Monomer A from Water-Methylamine """
TAGL['%s-%s-monoB-unCP' % (dbse, '3'                     )] = """Monomer B from Water-Methylamine """
TAGL['%s-%s'            % (dbse, '4'                     )] = """Water-N-methylacetamide """
TAGL['%s-%s-dimer'      % (dbse, '4'                     )] = """Dimer from Water-N-methylacetamide """
TAGL['%s-%s-monoA-CP'   % (dbse, '4'                     )] = """Monomer A from Water-N-methylacetamide """
TAGL['%s-%s-monoB-CP'   % (dbse, '4'                     )] = """Monomer B from Water-N-methylacetamide """
TAGL['%s-%s-monoA-unCP' % (dbse, '4'                     )] = """Monomer A from Water-N-methylacetamide """
TAGL['%s-%s-monoB-unCP' % (dbse, '4'                     )] = """Monomer B from Water-N-methylacetamide """
TAGL['%s-%s'            % (dbse, '5'                     )] = """Methanol Dimer """
TAGL['%s-%s-dimer'      % (dbse, '5'                     )] = """Dimer from Methanol Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '5'                     )] = """Monomer A from Methanol Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '5'                     )] = """Monomer B from Methanol Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '5'                     )] = """Monomer A from Methanol Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '5'                     )] = """Monomer B from Methanol Dimer """
TAGL['%s-%s'            % (dbse, '6'                     )] = """Methanol-Methylamine """
TAGL['%s-%s-dimer'      % (dbse, '6'                     )] = """Dimer from Methanol-Methylamine """
TAGL['%s-%s-monoA-CP'   % (dbse, '6'                     )] = """Monomer A from Methanol-Methylamine """
TAGL['%s-%s-monoB-CP'   % (dbse, '6'                     )] = """Monomer B from Methanol-Methylamine """
TAGL['%s-%s-monoA-unCP' % (dbse, '6'                     )] = """Monomer A from Methanol-Methylamine """
TAGL['%s-%s-monoB-unCP' % (dbse, '6'                     )] = """Monomer B from Methanol-Methylamine """
TAGL['%s-%s'            % (dbse, '7'                     )] = """Methanol-N-methylacetamide """
TAGL['%s-%s-dimer'      % (dbse, '7'                     )] = """Dimer from Methanol-N-methylacetamide """
TAGL['%s-%s-monoA-CP'   % (dbse, '7'                     )] = """Monomer A from Methanol-N-methylacetamide """
TAGL['%s-%s-monoB-CP'   % (dbse, '7'                     )] = """Monomer B from Methanol-N-methylacetamide """
TAGL['%s-%s-monoA-unCP' % (dbse, '7'                     )] = """Monomer A from Methanol-N-methylacetamide """
TAGL['%s-%s-monoB-unCP' % (dbse, '7'                     )] = """Monomer B from Methanol-N-methylacetamide """
TAGL['%s-%s'            % (dbse, '8'                     )] = """Methanol-Water """
TAGL['%s-%s-dimer'      % (dbse, '8'                     )] = """Dimer from Methanol-Water """
TAGL['%s-%s-monoA-CP'   % (dbse, '8'                     )] = """Monomer A from Methanol-Water """
TAGL['%s-%s-monoB-CP'   % (dbse, '8'                     )] = """Monomer B from Methanol-Water """
TAGL['%s-%s-monoA-unCP' % (dbse, '8'                     )] = """Monomer A from Methanol-Water """
TAGL['%s-%s-monoB-unCP' % (dbse, '8'                     )] = """Monomer B from Methanol-Water """
TAGL['%s-%s'            % (dbse, '9'                     )] = """Methylamine-Methanol """
TAGL['%s-%s-dimer'      % (dbse, '9'                     )] = """Dimer from Methylamine-Methanol """
TAGL['%s-%s-monoA-CP'   % (dbse, '9'                     )] = """Monomer A from Methylamine-Methanol """
TAGL['%s-%s-monoB-CP'   % (dbse, '9'                     )] = """Monomer B from Methylamine-Methanol """
TAGL['%s-%s-monoA-unCP' % (dbse, '9'                     )] = """Monomer A from Methylamine-Methanol """
TAGL['%s-%s-monoB-unCP' % (dbse, '9'                     )] = """Monomer B from Methylamine-Methanol """
TAGL['%s-%s'            % (dbse, '10'                    )] = """Methylamine Dimer """
TAGL['%s-%s-dimer'      % (dbse, '10'                    )] = """Dimer from Methylamine Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '10'                    )] = """Monomer A from Methylamine Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '10'                    )] = """Monomer B from Methylamine Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '10'                    )] = """Monomer A from Methylamine Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '10'                    )] = """Monomer B from Methylamine Dimer """
TAGL['%s-%s'            % (dbse, '11'                    )] = """Methylamine-N-methylacetamide """
TAGL['%s-%s-dimer'      % (dbse, '11'                    )] = """Dimer from Methylamine-N-methylacetamide """
TAGL['%s-%s-monoA-CP'   % (dbse, '11'                    )] = """Monomer A from Methylamine-N-methylacetamide """
TAGL['%s-%s-monoB-CP'   % (dbse, '11'                    )] = """Monomer B from Methylamine-N-methylacetamide """
TAGL['%s-%s-monoA-unCP' % (dbse, '11'                    )] = """Monomer A from Methylamine-N-methylacetamide """
TAGL['%s-%s-monoB-unCP' % (dbse, '11'                    )] = """Monomer B from Methylamine-N-methylacetamide """
TAGL['%s-%s'            % (dbse, '12'                    )] = """Methylamine-Water """
TAGL['%s-%s-dimer'      % (dbse, '12'                    )] = """Dimer from Methylamine-Water """
TAGL['%s-%s-monoA-CP'   % (dbse, '12'                    )] = """Monomer A from Methylamine-Water """
TAGL['%s-%s-monoB-CP'   % (dbse, '12'                    )] = """Monomer B from Methylamine-Water """
TAGL['%s-%s-monoA-unCP' % (dbse, '12'                    )] = """Monomer A from Methylamine-Water """
TAGL['%s-%s-monoB-unCP' % (dbse, '12'                    )] = """Monomer B from Methylamine-Water """
TAGL['%s-%s'            % (dbse, '13'                    )] = """N-methylacetamide-Methanol """
TAGL['%s-%s-dimer'      % (dbse, '13'                    )] = """Dimer from N-methylacetamide-Methanol """
TAGL['%s-%s-monoA-CP'   % (dbse, '13'                    )] = """Monomer A from N-methylacetamide-Methanol """
TAGL['%s-%s-monoB-CP'   % (dbse, '13'                    )] = """Monomer B from N-methylacetamide-Methanol """
TAGL['%s-%s-monoA-unCP' % (dbse, '13'                    )] = """Monomer A from N-methylacetamide-Methanol """
TAGL['%s-%s-monoB-unCP' % (dbse, '13'                    )] = """Monomer B from N-methylacetamide-Methanol """
TAGL['%s-%s'            % (dbse, '14'                    )] = """N-methylacetamide-Methylamine """
TAGL['%s-%s-dimer'      % (dbse, '14'                    )] = """Dimer from N-methylacetamide-Methylamine """
TAGL['%s-%s-monoA-CP'   % (dbse, '14'                    )] = """Monomer A from N-methylacetamide-Methylamine """
TAGL['%s-%s-monoB-CP'   % (dbse, '14'                    )] = """Monomer B from N-methylacetamide-Methylamine """
TAGL['%s-%s-monoA-unCP' % (dbse, '14'                    )] = """Monomer A from N-methylacetamide-Methylamine """
TAGL['%s-%s-monoB-unCP' % (dbse, '14'                    )] = """Monomer B from N-methylacetamide-Methylamine """
TAGL['%s-%s'            % (dbse, '15'                    )] = """N-methylacetamide Dimer """
TAGL['%s-%s-dimer'      % (dbse, '15'                    )] = """Dimer from N-methylacetamide Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '15'                    )] = """Monomer A from N-methylacetamide Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '15'                    )] = """Monomer B from N-methylacetamide Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '15'                    )] = """Monomer A from N-methylacetamide Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '15'                    )] = """Monomer B from N-methylacetamide Dimer """
TAGL['%s-%s'            % (dbse, '16'                    )] = """N-methylacetamide-Water """
TAGL['%s-%s-dimer'      % (dbse, '16'                    )] = """Dimer from N-methylacetamide-Water """
TAGL['%s-%s-monoA-CP'   % (dbse, '16'                    )] = """Monomer A from N-methylacetamide-Water """
TAGL['%s-%s-monoB-CP'   % (dbse, '16'                    )] = """Monomer B from N-methylacetamide-Water """
TAGL['%s-%s-monoA-unCP' % (dbse, '16'                    )] = """Monomer A from N-methylacetamide-Water """
TAGL['%s-%s-monoB-unCP' % (dbse, '16'                    )] = """Monomer B from N-methylacetamide-Water """
TAGL['%s-%s'            % (dbse, '17'                    )] = """Uracil Dimer, HB """
TAGL['%s-%s-dimer'      % (dbse, '17'                    )] = """Dimer from Uracil Dimer, HB """
TAGL['%s-%s-monoA-CP'   % (dbse, '17'                    )] = """Monomer A from Uracil Dimer, HB """
TAGL['%s-%s-monoB-CP'   % (dbse, '17'                    )] = """Monomer B from Uracil Dimer, HB """
TAGL['%s-%s-monoA-unCP' % (dbse, '17'                    )] = """Monomer A from Uracil Dimer, HB """
TAGL['%s-%s-monoB-unCP' % (dbse, '17'                    )] = """Monomer B from Uracil Dimer, HB """
TAGL['%s-%s'            % (dbse, '18'                    )] = """Water-Pyridine """
TAGL['%s-%s-dimer'      % (dbse, '18'                    )] = """Dimer from Water-Pyridine """
TAGL['%s-%s-monoA-CP'   % (dbse, '18'                    )] = """Monomer A from Water-Pyridine """
TAGL['%s-%s-monoB-CP'   % (dbse, '18'                    )] = """Monomer B from Water-Pyridine """
TAGL['%s-%s-monoA-unCP' % (dbse, '18'                    )] = """Monomer A from Water-Pyridine """
TAGL['%s-%s-monoB-unCP' % (dbse, '18'                    )] = """Monomer B from Water-Pyridine """
TAGL['%s-%s'            % (dbse, '19'                    )] = """Methanol-Pyridine """
TAGL['%s-%s-dimer'      % (dbse, '19'                    )] = """Dimer from Methanol-Pyridine """
TAGL['%s-%s-monoA-CP'   % (dbse, '19'                    )] = """Monomer A from Methanol-Pyridine """
TAGL['%s-%s-monoB-CP'   % (dbse, '19'                    )] = """Monomer B from Methanol-Pyridine """
TAGL['%s-%s-monoA-unCP' % (dbse, '19'                    )] = """Monomer A from Methanol-Pyridine """
TAGL['%s-%s-monoB-unCP' % (dbse, '19'                    )] = """Monomer B from Methanol-Pyridine """
TAGL['%s-%s'            % (dbse, '20'                    )] = """Acetic Acid Dimer """
TAGL['%s-%s-dimer'      % (dbse, '20'                    )] = """Dimer from Acetic Acid Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '20'                    )] = """Monomer A from Acetic Acid Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '20'                    )] = """Monomer B from Acetic Acid Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '20'                    )] = """Monomer A from Acetic Acid Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '20'                    )] = """Monomer B from Acetic Acid Dimer """
TAGL['%s-%s'            % (dbse, '21'                    )] = """Acetamide Dimer """
TAGL['%s-%s-dimer'      % (dbse, '21'                    )] = """Dimer from Acetamide Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '21'                    )] = """Monomer A from Acetamide Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '21'                    )] = """Monomer B from Acetamide Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '21'                    )] = """Monomer A from Acetamide Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '21'                    )] = """Monomer B from Acetamide Dimer """
TAGL['%s-%s'            % (dbse, '22'                    )] = """Acetic Acid-Uracil """
TAGL['%s-%s-dimer'      % (dbse, '22'                    )] = """Dimer from Acetic Acid-Uracil """
TAGL['%s-%s-monoA-CP'   % (dbse, '22'                    )] = """Monomer A from Acetic Acid-Uracil """
TAGL['%s-%s-monoB-CP'   % (dbse, '22'                    )] = """Monomer B from Acetic Acid-Uracil """
TAGL['%s-%s-monoA-unCP' % (dbse, '22'                    )] = """Monomer A from Acetic Acid-Uracil """
TAGL['%s-%s-monoB-unCP' % (dbse, '22'                    )] = """Monomer B from Acetic Acid-Uracil """
TAGL['%s-%s'            % (dbse, '23'                    )] = """Acetamide-Uracil """
TAGL['%s-%s-dimer'      % (dbse, '23'                    )] = """Dimer from Acetamide-Uracil """
TAGL['%s-%s-monoA-CP'   % (dbse, '23'                    )] = """Monomer A from Acetamide-Uracil """
TAGL['%s-%s-monoB-CP'   % (dbse, '23'                    )] = """Monomer B from Acetamide-Uracil """
TAGL['%s-%s-monoA-unCP' % (dbse, '23'                    )] = """Monomer A from Acetamide-Uracil """
TAGL['%s-%s-monoB-unCP' % (dbse, '23'                    )] = """Monomer B from Acetamide-Uracil """
TAGL['%s-%s'            % (dbse, '24'                    )] = """Benzene Dimer, pi-pi """
TAGL['%s-%s-dimer'      % (dbse, '24'                    )] = """Dimer from Benzene Dimer, pi-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '24'                    )] = """Monomer A from Benzene Dimer, pi-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '24'                    )] = """Monomer B from Benzene Dimer, pi-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '24'                    )] = """Monomer A from Benzene Dimer, pi-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '24'                    )] = """Monomer B from Benzene Dimer, pi-pi """
TAGL['%s-%s'            % (dbse, '25'                    )] = """Pyridine Dimer, pi-pi """
TAGL['%s-%s-dimer'      % (dbse, '25'                    )] = """Dimer from Pyridine Dimer, pi-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '25'                    )] = """Monomer A from Pyridine Dimer, pi-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '25'                    )] = """Monomer B from Pyridine Dimer, pi-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '25'                    )] = """Monomer A from Pyridine Dimer, pi-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '25'                    )] = """Monomer B from Pyridine Dimer, pi-pi """
TAGL['%s-%s'            % (dbse, '26'                    )] = """Uracil Dimer, pi-pi """
TAGL['%s-%s-dimer'      % (dbse, '26'                    )] = """Dimer from Uracil Dimer, pi-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '26'                    )] = """Monomer A from Uracil Dimer, pi-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '26'                    )] = """Monomer B from Uracil Dimer, pi-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '26'                    )] = """Monomer A from Uracil Dimer, pi-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '26'                    )] = """Monomer B from Uracil Dimer, pi-pi """
TAGL['%s-%s'            % (dbse, '27'                    )] = """Benzene-Pyridine, pi-pi """
TAGL['%s-%s-dimer'      % (dbse, '27'                    )] = """Dimer from Benzene-Pyridine, pi-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '27'                    )] = """Monomer A from Benzene-Pyridine, pi-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '27'                    )] = """Monomer B from Benzene-Pyridine, pi-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '27'                    )] = """Monomer A from Benzene-Pyridine, pi-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '27'                    )] = """Monomer B from Benzene-Pyridine, pi-pi """
TAGL['%s-%s'            % (dbse, '28'                    )] = """Benzene-Uracil, pi-pi """
TAGL['%s-%s-dimer'      % (dbse, '28'                    )] = """Dimer from Benzene-Uracil, pi-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '28'                    )] = """Monomer A from Benzene-Uracil, pi-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '28'                    )] = """Monomer B from Benzene-Uracil, pi-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '28'                    )] = """Monomer A from Benzene-Uracil, pi-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '28'                    )] = """Monomer B from Benzene-Uracil, pi-pi """
TAGL['%s-%s'            % (dbse, '29'                    )] = """Pyridine-Uracil, pi-pi """
TAGL['%s-%s-dimer'      % (dbse, '29'                    )] = """Dimer from Pyridine-Uracil, pi-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '29'                    )] = """Monomer A from Pyridine-Uracil, pi-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '29'                    )] = """Monomer B from Pyridine-Uracil, pi-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '29'                    )] = """Monomer A from Pyridine-Uracil, pi-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '29'                    )] = """Monomer B from Pyridine-Uracil, pi-pi """
TAGL['%s-%s'            % (dbse, '30'                    )] = """Benzene-Ethene """
TAGL['%s-%s-dimer'      % (dbse, '30'                    )] = """Dimer from Benzene-Ethene """
TAGL['%s-%s-monoA-CP'   % (dbse, '30'                    )] = """Monomer A from Benzene-Ethene """
TAGL['%s-%s-monoB-CP'   % (dbse, '30'                    )] = """Monomer B from Benzene-Ethene """
TAGL['%s-%s-monoA-unCP' % (dbse, '30'                    )] = """Monomer A from Benzene-Ethene """
TAGL['%s-%s-monoB-unCP' % (dbse, '30'                    )] = """Monomer B from Benzene-Ethene """
TAGL['%s-%s'            % (dbse, '31'                    )] = """Uracil-Ethene """
TAGL['%s-%s-dimer'      % (dbse, '31'                    )] = """Dimer from Uracil-Ethene """
TAGL['%s-%s-monoA-CP'   % (dbse, '31'                    )] = """Monomer A from Uracil-Ethene """
TAGL['%s-%s-monoB-CP'   % (dbse, '31'                    )] = """Monomer B from Uracil-Ethene """
TAGL['%s-%s-monoA-unCP' % (dbse, '31'                    )] = """Monomer A from Uracil-Ethene """
TAGL['%s-%s-monoB-unCP' % (dbse, '31'                    )] = """Monomer B from Uracil-Ethene """
TAGL['%s-%s'            % (dbse, '32'                    )] = """Uracil-Ethyne """
TAGL['%s-%s-dimer'      % (dbse, '32'                    )] = """Dimer from Uracil-Ethyne """
TAGL['%s-%s-monoA-CP'   % (dbse, '32'                    )] = """Monomer A from Uracil-Ethyne """
TAGL['%s-%s-monoB-CP'   % (dbse, '32'                    )] = """Monomer B from Uracil-Ethyne """
TAGL['%s-%s-monoA-unCP' % (dbse, '32'                    )] = """Monomer A from Uracil-Ethyne """
TAGL['%s-%s-monoB-unCP' % (dbse, '32'                    )] = """Monomer B from Uracil-Ethyne """
TAGL['%s-%s'            % (dbse, '33'                    )] = """Pyridine-Ethene """
TAGL['%s-%s-dimer'      % (dbse, '33'                    )] = """Dimer from Pyridine-Ethene """
TAGL['%s-%s-monoA-CP'   % (dbse, '33'                    )] = """Monomer A from Pyridine-Ethene """
TAGL['%s-%s-monoB-CP'   % (dbse, '33'                    )] = """Monomer B from Pyridine-Ethene """
TAGL['%s-%s-monoA-unCP' % (dbse, '33'                    )] = """Monomer A from Pyridine-Ethene """
TAGL['%s-%s-monoB-unCP' % (dbse, '33'                    )] = """Monomer B from Pyridine-Ethene """
TAGL['%s-%s'            % (dbse, '34'                    )] = """Pentane Dimer """
TAGL['%s-%s-dimer'      % (dbse, '34'                    )] = """Dimer from Pentane Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '34'                    )] = """Monomer A from Pentane Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '34'                    )] = """Monomer B from Pentane Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '34'                    )] = """Monomer A from Pentane Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '34'                    )] = """Monomer B from Pentane Dimer """
TAGL['%s-%s'            % (dbse, '35'                    )] = """Neopentane-Pentane """
TAGL['%s-%s-dimer'      % (dbse, '35'                    )] = """Dimer from Neopentane-Pentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '35'                    )] = """Monomer A from Neopentane-Pentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '35'                    )] = """Monomer B from Neopentane-Pentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '35'                    )] = """Monomer A from Neopentane-Pentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '35'                    )] = """Monomer B from Neopentane-Pentane """
TAGL['%s-%s'            % (dbse, '36'                    )] = """Neopentane Dimer """
TAGL['%s-%s-dimer'      % (dbse, '36'                    )] = """Dimer from Neopentane Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '36'                    )] = """Monomer A from Neopentane Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '36'                    )] = """Monomer B from Neopentane Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '36'                    )] = """Monomer A from Neopentane Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '36'                    )] = """Monomer B from Neopentane Dimer """
TAGL['%s-%s'            % (dbse, '37'                    )] = """Cyclopentane-Neopentane """
TAGL['%s-%s-dimer'      % (dbse, '37'                    )] = """Dimer from Cyclopentane-Neopentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '37'                    )] = """Monomer A from Cyclopentane-Neopentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '37'                    )] = """Monomer B from Cyclopentane-Neopentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '37'                    )] = """Monomer A from Cyclopentane-Neopentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '37'                    )] = """Monomer B from Cyclopentane-Neopentane """
TAGL['%s-%s'            % (dbse, '38'                    )] = """Cyclopentane Dimer """
TAGL['%s-%s-dimer'      % (dbse, '38'                    )] = """Dimer from Cyclopentane Dimer """
TAGL['%s-%s-monoA-CP'   % (dbse, '38'                    )] = """Monomer A from Cyclopentane Dimer """
TAGL['%s-%s-monoB-CP'   % (dbse, '38'                    )] = """Monomer B from Cyclopentane Dimer """
TAGL['%s-%s-monoA-unCP' % (dbse, '38'                    )] = """Monomer A from Cyclopentane Dimer """
TAGL['%s-%s-monoB-unCP' % (dbse, '38'                    )] = """Monomer B from Cyclopentane Dimer """
TAGL['%s-%s'            % (dbse, '39'                    )] = """Benzene-Cyclopentane """
TAGL['%s-%s-dimer'      % (dbse, '39'                    )] = """Dimer from Benzene-Cyclopentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '39'                    )] = """Monomer A from Benzene-Cyclopentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '39'                    )] = """Monomer B from Benzene-Cyclopentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '39'                    )] = """Monomer A from Benzene-Cyclopentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '39'                    )] = """Monomer B from Benzene-Cyclopentane """
TAGL['%s-%s'            % (dbse, '40'                    )] = """Benzene-Neopentane """
TAGL['%s-%s-dimer'      % (dbse, '40'                    )] = """Dimer from Benzene-Neopentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '40'                    )] = """Monomer A from Benzene-Neopentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '40'                    )] = """Monomer B from Benzene-Neopentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '40'                    )] = """Monomer A from Benzene-Neopentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '40'                    )] = """Monomer B from Benzene-Neopentane """
TAGL['%s-%s'            % (dbse, '41'                    )] = """Uracil-Pentane """
TAGL['%s-%s-dimer'      % (dbse, '41'                    )] = """Dimer from Uracil-Pentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '41'                    )] = """Monomer A from Uracil-Pentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '41'                    )] = """Monomer B from Uracil-Pentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '41'                    )] = """Monomer A from Uracil-Pentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '41'                    )] = """Monomer B from Uracil-Pentane """
TAGL['%s-%s'            % (dbse, '42'                    )] = """Uracil-Cyclopentane """
TAGL['%s-%s-dimer'      % (dbse, '42'                    )] = """Dimer from Uracil-Cyclopentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '42'                    )] = """Monomer A from Uracil-Cyclopentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '42'                    )] = """Monomer B from Uracil-Cyclopentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '42'                    )] = """Monomer A from Uracil-Cyclopentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '42'                    )] = """Monomer B from Uracil-Cyclopentane """
TAGL['%s-%s'            % (dbse, '43'                    )] = """Uracil-Neopentane """
TAGL['%s-%s-dimer'      % (dbse, '43'                    )] = """Dimer from Uracil-Neopentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '43'                    )] = """Monomer A from Uracil-Neopentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '43'                    )] = """Monomer B from Uracil-Neopentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '43'                    )] = """Monomer A from Uracil-Neopentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '43'                    )] = """Monomer B from Uracil-Neopentane """
TAGL['%s-%s'            % (dbse, '44'                    )] = """Ethene-Pentane """
TAGL['%s-%s-dimer'      % (dbse, '44'                    )] = """Dimer from Ethene-Pentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '44'                    )] = """Monomer A from Ethene-Pentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '44'                    )] = """Monomer B from Ethene-Pentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '44'                    )] = """Monomer A from Ethene-Pentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '44'                    )] = """Monomer B from Ethene-Pentane """
TAGL['%s-%s'            % (dbse, '45'                    )] = """Ethyne-Pentane """
TAGL['%s-%s-dimer'      % (dbse, '45'                    )] = """Dimer from Ethyne-Pentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '45'                    )] = """Monomer A from Ethyne-Pentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '45'                    )] = """Monomer B from Ethyne-Pentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '45'                    )] = """Monomer A from Ethyne-Pentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '45'                    )] = """Monomer B from Ethyne-Pentane """
TAGL['%s-%s'            % (dbse, '46'                    )] = """N-methylacetamide-Pentane """
TAGL['%s-%s-dimer'      % (dbse, '46'                    )] = """Dimer from N-methylacetamide-Pentane """
TAGL['%s-%s-monoA-CP'   % (dbse, '46'                    )] = """Monomer A from N-methylacetamide-Pentane """
TAGL['%s-%s-monoB-CP'   % (dbse, '46'                    )] = """Monomer B from N-methylacetamide-Pentane """
TAGL['%s-%s-monoA-unCP' % (dbse, '46'                    )] = """Monomer A from N-methylacetamide-Pentane """
TAGL['%s-%s-monoB-unCP' % (dbse, '46'                    )] = """Monomer B from N-methylacetamide-Pentane """
TAGL['%s-%s'            % (dbse, '47'                    )] = """Benzene Dimer, CH-pi """
TAGL['%s-%s-dimer'      % (dbse, '47'                    )] = """Dimer from Benzene Dimer, CH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '47'                    )] = """Monomer A from Benzene Dimer, CH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '47'                    )] = """Monomer B from Benzene Dimer, CH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '47'                    )] = """Monomer A from Benzene Dimer, CH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '47'                    )] = """Monomer B from Benzene Dimer, CH-pi """
TAGL['%s-%s'            % (dbse, '48'                    )] = """Pyridine Dimer, CH-pi """
TAGL['%s-%s-dimer'      % (dbse, '48'                    )] = """Dimer from Pyridine Dimer, CH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '48'                    )] = """Monomer A from Pyridine Dimer, CH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '48'                    )] = """Monomer B from Pyridine Dimer, CH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '48'                    )] = """Monomer A from Pyridine Dimer, CH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '48'                    )] = """Monomer B from Pyridine Dimer, CH-pi """
TAGL['%s-%s'            % (dbse, '49'                    )] = """Benzene-Pyridine, CH-pi """
TAGL['%s-%s-dimer'      % (dbse, '49'                    )] = """Dimer from Benzene-Pyridine, CH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '49'                    )] = """Monomer A from Benzene-Pyridine, CH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '49'                    )] = """Monomer B from Benzene-Pyridine, CH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '49'                    )] = """Monomer A from Benzene-Pyridine, CH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '49'                    )] = """Monomer B from Benzene-Pyridine, CH-pi """
TAGL['%s-%s'            % (dbse, '50'                    )] = """Benzene-Ethyne, CH-pi """
TAGL['%s-%s-dimer'      % (dbse, '50'                    )] = """Dimer from Benzene-Ethyne, CH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '50'                    )] = """Monomer A from Benzene-Ethyne, CH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '50'                    )] = """Monomer B from Benzene-Ethyne, CH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '50'                    )] = """Monomer A from Benzene-Ethyne, CH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '50'                    )] = """Monomer B from Benzene-Ethyne, CH-pi """
TAGL['%s-%s'            % (dbse, '51'                    )] = """Ethyne Dimer, CH-pi """
TAGL['%s-%s-dimer'      % (dbse, '51'                    )] = """Dimer from Ethyne Dimer, CH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '51'                    )] = """Monomer A from Ethyne Dimer, CH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '51'                    )] = """Monomer B from Ethyne Dimer, CH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '51'                    )] = """Monomer A from Ethyne Dimer, CH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '51'                    )] = """Monomer B from Ethyne Dimer, CH-pi """
TAGL['%s-%s'            % (dbse, '52'                    )] = """Benzene-Acetic Acid, OH-pi """
TAGL['%s-%s-dimer'      % (dbse, '52'                    )] = """Dimer from Benzene-Acetic Acid, OH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '52'                    )] = """Monomer A from Benzene-Acetic Acid, OH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '52'                    )] = """Monomer B from Benzene-Acetic Acid, OH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '52'                    )] = """Monomer A from Benzene-Acetic Acid, OH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '52'                    )] = """Monomer B from Benzene-Acetic Acid, OH-pi """
TAGL['%s-%s'            % (dbse, '53'                    )] = """Benzene-Acetamide, NH-pi """
TAGL['%s-%s-dimer'      % (dbse, '53'                    )] = """Dimer from Benzene-Acetamide, NH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '53'                    )] = """Monomer A from Benzene-Acetamide, NH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '53'                    )] = """Monomer B from Benzene-Acetamide, NH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '53'                    )] = """Monomer A from Benzene-Acetamide, NH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '53'                    )] = """Monomer B from Benzene-Acetamide, NH-pi """
TAGL['%s-%s'            % (dbse, '54'                    )] = """Benzene-Water, OH-pi """
TAGL['%s-%s-dimer'      % (dbse, '54'                    )] = """Dimer from Benzene-Water, OH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '54'                    )] = """Monomer A from Benzene-Water, OH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '54'                    )] = """Monomer B from Benzene-Water, OH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '54'                    )] = """Monomer A from Benzene-Water, OH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '54'                    )] = """Monomer B from Benzene-Water, OH-pi """
TAGL['%s-%s'            % (dbse, '55'                    )] = """Benzene-Methanol, OH-pi """
TAGL['%s-%s-dimer'      % (dbse, '55'                    )] = """Dimer from Benzene-Methanol, OH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '55'                    )] = """Monomer A from Benzene-Methanol, OH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '55'                    )] = """Monomer B from Benzene-Methanol, OH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '55'                    )] = """Monomer A from Benzene-Methanol, OH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '55'                    )] = """Monomer B from Benzene-Methanol, OH-pi """
TAGL['%s-%s'            % (dbse, '56'                    )] = """Benzene-Methylamine, NH-pi """
TAGL['%s-%s-dimer'      % (dbse, '56'                    )] = """Dimer from Benzene-Methylamine, NH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '56'                    )] = """Monomer A from Benzene-Methylamine, NH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '56'                    )] = """Monomer B from Benzene-Methylamine, NH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '56'                    )] = """Monomer A from Benzene-Methylamine, NH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '56'                    )] = """Monomer B from Benzene-Methylamine, NH-pi """
TAGL['%s-%s'            % (dbse, '57'                    )] = """Benzene-N-methylacetamide, NH-pi """
TAGL['%s-%s-dimer'      % (dbse, '57'                    )] = """Dimer from Benzene-N-methylacetamide, NH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '57'                    )] = """Monomer A from Benzene-N-methylacetamide, NH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '57'                    )] = """Monomer B from Benzene-N-methylacetamide, NH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '57'                    )] = """Monomer A from Benzene-N-methylacetamide, NH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '57'                    )] = """Monomer B from Benzene-N-methylacetamide, NH-pi """
TAGL['%s-%s'            % (dbse, '58'                    )] = """Pyridine Dimer, CH-N """
TAGL['%s-%s-dimer'      % (dbse, '58'                    )] = """Dimer from Pyridine Dimer, CH-N """
TAGL['%s-%s-monoA-CP'   % (dbse, '58'                    )] = """Monomer A from Pyridine Dimer, CH-N """
TAGL['%s-%s-monoB-CP'   % (dbse, '58'                    )] = """Monomer B from Pyridine Dimer, CH-N """
TAGL['%s-%s-monoA-unCP' % (dbse, '58'                    )] = """Monomer A from Pyridine Dimer, CH-N """
TAGL['%s-%s-monoB-unCP' % (dbse, '58'                    )] = """Monomer B from Pyridine Dimer, CH-N """
TAGL['%s-%s'            % (dbse, '59'                    )] = """Ethyne-Water, CH-O """
TAGL['%s-%s-dimer'      % (dbse, '59'                    )] = """Dimer from Ethyne-Water, CH-O """
TAGL['%s-%s-monoA-CP'   % (dbse, '59'                    )] = """Monomer A from Ethyne-Water, CH-O """
TAGL['%s-%s-monoB-CP'   % (dbse, '59'                    )] = """Monomer B from Ethyne-Water, CH-O """
TAGL['%s-%s-monoA-unCP' % (dbse, '59'                    )] = """Monomer A from Ethyne-Water, CH-O """
TAGL['%s-%s-monoB-unCP' % (dbse, '59'                    )] = """Monomer B from Ethyne-Water, CH-O """
TAGL['%s-%s'            % (dbse, '60'                    )] = """Ethyne-Acetic Acid, OH-pi """
TAGL['%s-%s-dimer'      % (dbse, '60'                    )] = """Dimer from Ethyne-Acetic Acid, OH-pi """
TAGL['%s-%s-monoA-CP'   % (dbse, '60'                    )] = """Monomer A from Ethyne-Acetic Acid, OH-pi """
TAGL['%s-%s-monoB-CP'   % (dbse, '60'                    )] = """Monomer B from Ethyne-Acetic Acid, OH-pi """
TAGL['%s-%s-monoA-unCP' % (dbse, '60'                    )] = """Monomer A from Ethyne-Acetic Acid, OH-pi """
TAGL['%s-%s-monoB-unCP' % (dbse, '60'                    )] = """Monomer B from Ethyne-Acetic Acid, OH-pi """
TAGL['%s-%s'            % (dbse, '61'                    )] = """Pentane-Acetic Acid """
TAGL['%s-%s-dimer'      % (dbse, '61'                    )] = """Dimer from Pentane-Acetic Acid """
TAGL['%s-%s-monoA-CP'   % (dbse, '61'                    )] = """Monomer A from Pentane-Acetic Acid """
TAGL['%s-%s-monoB-CP'   % (dbse, '61'                    )] = """Monomer B from Pentane-Acetic Acid """
TAGL['%s-%s-monoA-unCP' % (dbse, '61'                    )] = """Monomer A from Pentane-Acetic Acid """
TAGL['%s-%s-monoB-unCP' % (dbse, '61'                    )] = """Monomer B from Pentane-Acetic Acid """
TAGL['%s-%s'            % (dbse, '62'                    )] = """Pentane-Acetamide """
TAGL['%s-%s-dimer'      % (dbse, '62'                    )] = """Dimer from Pentane-Acetamide """
TAGL['%s-%s-monoA-CP'   % (dbse, '62'                    )] = """Monomer A from Pentane-Acetamide """
TAGL['%s-%s-monoB-CP'   % (dbse, '62'                    )] = """Monomer B from Pentane-Acetamide """
TAGL['%s-%s-monoA-unCP' % (dbse, '62'                    )] = """Monomer A from Pentane-Acetamide """
TAGL['%s-%s-monoB-unCP' % (dbse, '62'                    )] = """Monomer B from Pentane-Acetamide """
TAGL['%s-%s'            % (dbse, '63'                    )] = """Benzene-Acetic Acid """
TAGL['%s-%s-dimer'      % (dbse, '63'                    )] = """Dimer from Benzene-Acetic Acid """
TAGL['%s-%s-monoA-CP'   % (dbse, '63'                    )] = """Monomer A from Benzene-Acetic Acid """
TAGL['%s-%s-monoB-CP'   % (dbse, '63'                    )] = """Monomer B from Benzene-Acetic Acid """
TAGL['%s-%s-monoA-unCP' % (dbse, '63'                    )] = """Monomer A from Benzene-Acetic Acid """
TAGL['%s-%s-monoB-unCP' % (dbse, '63'                    )] = """Monomer B from Benzene-Acetic Acid """
TAGL['%s-%s'            % (dbse, '64'                    )] = """N-methylacetamide-Ethene """
TAGL['%s-%s-dimer'      % (dbse, '64'                    )] = """Dimer from N-methylacetamide-Ethene """
TAGL['%s-%s-monoA-CP'   % (dbse, '64'                    )] = """Monomer A from N-methylacetamide-Ethene """
TAGL['%s-%s-monoB-CP'   % (dbse, '64'                    )] = """Monomer B from N-methylacetamide-Ethene """
TAGL['%s-%s-monoA-unCP' % (dbse, '64'                    )] = """Monomer A from N-methylacetamide-Ethene """
TAGL['%s-%s-monoB-unCP' % (dbse, '64'                    )] = """Monomer B from N-methylacetamide-Ethene """
TAGL['%s-%s'            % (dbse, '65'                    )] = """Pyridine-Ethyne """
TAGL['%s-%s-dimer'      % (dbse, '65'                    )] = """Dimer from Pyridine-Ethyne """
TAGL['%s-%s-monoA-CP'   % (dbse, '65'                    )] = """Monomer A from Pyridine-Ethyne """
TAGL['%s-%s-monoB-CP'   % (dbse, '65'                    )] = """Monomer B from Pyridine-Ethyne """
TAGL['%s-%s-monoA-unCP' % (dbse, '65'                    )] = """Monomer A from Pyridine-Ethyne """
TAGL['%s-%s-monoB-unCP' % (dbse, '65'                    )] = """Monomer B from Pyridine-Ethyne """
TAGL['%s-%s'            % (dbse, '66'                    )] = """Methylamine-Pyridine """
TAGL['%s-%s-dimer'      % (dbse, '66'                    )] = """Dimer from Methylamine-Pyridine """
TAGL['%s-%s-monoA-CP'   % (dbse, '66'                    )] = """Monomer A from Methylamine-Pyridine """
TAGL['%s-%s-monoB-CP'   % (dbse, '66'                    )] = """Monomer B from Methylamine-Pyridine """
TAGL['%s-%s-monoA-unCP' % (dbse, '66'                    )] = """Monomer A from Methylamine-Pyridine """
TAGL['%s-%s-monoB-unCP' % (dbse, '66'                    )] = """Monomer B from Methylamine-Pyridine """

# <<< Geometry Specification Strings >>>
GEOS = {}

GEOS['%s-%s-dimer' % (dbse, '1')] = qcdb.Molecule("""
efp water
--
efp water
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '2')] = qcdb.Molecule("""
efp water
--
efp methanol
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '3')] = qcdb.Molecule("""
efp water
--
efp methylamine
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '4')] = qcdb.Molecule("""
efp water
--
efp nmethylactamide
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '5')] = qcdb.Molecule("""
efp methanol
--
efp methanol
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '6')] = qcdb.Molecule("""
efp methanol
--
efp methylamine
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '7')] = qcdb.Molecule("""
efp methanol
--
efp nmethylacetamide
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '8')] = qcdb.Molecule("""
efp methanol
--
efp water
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '9')] = qcdb.Molecule("""
efp methylamine
--
efp methanol
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '10')] = qcdb.Molecule("""
efp methylamine
--
efp methylamine
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '11')] = qcdb.Molecule("""
efp methylamine
--
efp nmethylacetamide
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '12')] = qcdb.Molecule("""
efp methylamine
--
efp water
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '13')] = qcdb.Molecule("""
efp nmethylacetamide
--
efp methanol
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '14')] = qcdb.Molecule("""
efp nmethylacetamide
--
efp methylamine
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '15')] = qcdb.Molecule("""
efp nmethylacetamide
--
efp nmethylacetamide
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '16')] = qcdb.Molecule("""
efp nmethylacetamide
--
efp water
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '17')] = qcdb.Molecule("""
efp uracil-gp
--
efp uracil-gp
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '18')] = qcdb.Molecule("""
efp water
--
efp pyridine
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '19')] = qcdb.Molecule("""
efp methanol
--
efp pyridine
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '20')] = qcdb.Molecule("""
efp aceticacid-hb
--
efp aceticacid-hb
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '21')] = qcdb.Molecule("""
efp acetamide-hb
--
efp acetamide-hb
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '22')] = qcdb.Molecule("""
efp aceticacid-hb
--
efp uracil-gp
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '23')] = qcdb.Molecule("""
efp acetamide-hb
--
efp uracil-gp
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '24')] = qcdb.Molecule("""
efp benzene
--
efp benzene
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '25')] = qcdb.Molecule("""
efp pyridine
--
efp pyridine
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '26')] = qcdb.Molecule("""
efp uracil-gp
--
efp uracil-gp
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '27')] = qcdb.Molecule("""
efp benzene
--
efp pyridine
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '28')] = qcdb.Molecule("""
efp benzene
--
efp uracil-gp
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '29')] = qcdb.Molecule("""
efp pyridine
--
efp uracil-gp
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '30')] = qcdb.Molecule("""
efp benzene
--
efp ethene
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '31')] = qcdb.Molecule("""
efp uracil-gp
--
efp ethene
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '32')] = qcdb.Molecule("""
efp uracil-gp
--
efp ethyne
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '33')] = qcdb.Molecule("""
efp pyridine
--
efp ethene
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '34')] = qcdb.Molecule("""
efp pentane
--
efp pentane
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '35')] = qcdb.Molecule("""
efp neopentane
--
efp pentane
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '36')] = qcdb.Molecule("""
efp neopentane
--
efp neopentane
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '37')] = qcdb.Molecule("""
efp cyclopentane
--
efp neopentane
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '38')] = qcdb.Molecule("""
efp cyclopentane
--
efp cyclopentane
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '39')] = qcdb.Molecule("""
efp benzene
--
efp cyclopentane
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '40')] = qcdb.Molecule("""
efp benzene
--
efp neopentane
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '41')] = qcdb.Molecule("""
efp uracil-gp
--
efp pentane
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '42')] = qcdb.Molecule("""
efp uracil-gp
--
efp cyclopentane
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '43')] = qcdb.Molecule("""
efp uracil-gp
--
efp neopentane
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '44')] = qcdb.Molecule("""
efp ethene
--
efp pentane
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '45')] = qcdb.Molecule("""
efp ethyne
--
efp pentane
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '46')] = qcdb.Molecule("""
efp nmethylacetamide
--
efp pentane
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '47')] = qcdb.Molecule("""
efp benzene
--
efp benzene
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '48')] = qcdb.Molecule("""
efp pyridine
--
efp pyridine
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '49')] = qcdb.Molecule("""
efp benzene
--
efp pyridine
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '50')] = qcdb.Molecule("""
efp benzene
--
efp ethyne
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '51')] = qcdb.Molecule("""
efp ethyne
--
efp ethyne
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '52')] = qcdb.Molecule("""
efp benzene
--
efp aceticacid-gp
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '53')] = qcdb.Molecule("""
efp benzene
--
efp acetamide-gp
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '54')] = qcdb.Molecule("""
efp benzene
--
efp water
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '55')] = qcdb.Molecule("""
efp benzene
--
efp methanol
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '56')] = qcdb.Molecule("""
efp benzene
--
efp methylamine
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '57')] = qcdb.Molecule("""
efp benzene
--
efp nmethylacetamide
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '58')] = qcdb.Molecule("""
efp pyridine
--
efp pyridine
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '59')] = qcdb.Molecule("""
efp ethyne
--
efp water
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '60')] = qcdb.Molecule("""
efp ethyne
--
efp aceticacid-gp
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '61')] = qcdb.Molecule("""
efp pentane
--
efp aceticacid-gp
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '62')] = qcdb.Molecule("""
efp pentane
--
efp acetamide-gp
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '63')] = qcdb.Molecule("""
efp benzene
--
efp aceticacid-gp
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '64')] = qcdb.Molecule("""
efp nmethylacetamide
--
efp ethene
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '65')] = qcdb.Molecule("""
efp pyridine
--
efp ethyne
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '66')] = qcdb.Molecule("""
efp methylamine
--
efp pyridine
units angstrom
""")

# <<< Derived Geometry Strings >>>
for rxn in HRXN:
    GEOS['%s-%s-monoA-unCP' % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(1)
    GEOS['%s-%s-monoB-unCP' % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(2)
    GEOS['%s-%s-monoA-CP'   % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(1, 2)
    GEOS['%s-%s-monoB-CP'   % (dbse, rxn)] = GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(2, 1)

#########################################################################

# <<< Supplementary Quantum Chemical Results >>>
DATA = {}

DATA['NUCLEAR REPULSION ENERGY'] = {}
DATA['NUCLEAR REPULSION ENERGY']['S66-1-dimer'                    ] =      36.51369349
DATA['NUCLEAR REPULSION ENERGY']['S66-1-monoA-unCP'               ] =       9.15671411
DATA['NUCLEAR REPULSION ENERGY']['S66-1-monoB-unCP'               ] =       9.17259114
DATA['NUCLEAR REPULSION ENERGY']['S66-2-dimer'                    ] =      79.98338083
DATA['NUCLEAR REPULSION ENERGY']['S66-2-monoA-unCP'               ] =       9.14996836
DATA['NUCLEAR REPULSION ENERGY']['S66-2-monoB-unCP'               ] =      40.29463192
DATA['NUCLEAR REPULSION ENERGY']['S66-3-dimer'                    ] =      79.77996002
DATA['NUCLEAR REPULSION ENERGY']['S66-3-monoA-unCP'               ] =       9.12565570
DATA['NUCLEAR REPULSION ENERGY']['S66-3-monoB-unCP'               ] =      42.06267577
DATA['NUCLEAR REPULSION ENERGY']['S66-4-dimer'                    ] =     246.86074225
DATA['NUCLEAR REPULSION ENERGY']['S66-4-monoA-unCP'               ] =       9.13184124
DATA['NUCLEAR REPULSION ENERGY']['S66-4-monoB-unCP'               ] =     180.56084030
DATA['NUCLEAR REPULSION ENERGY']['S66-5-dimer'                    ] =     129.52156842
DATA['NUCLEAR REPULSION ENERGY']['S66-5-monoA-unCP'               ] =      40.41731272
DATA['NUCLEAR REPULSION ENERGY']['S66-5-monoB-unCP'               ] =      40.29806380
DATA['NUCLEAR REPULSION ENERGY']['S66-6-dimer'                    ] =     131.81617640
DATA['NUCLEAR REPULSION ENERGY']['S66-6-monoA-unCP'               ] =      40.42467073
DATA['NUCLEAR REPULSION ENERGY']['S66-6-monoB-unCP'               ] =      42.05202847
DATA['NUCLEAR REPULSION ENERGY']['S66-7-dimer'                    ] =     313.95975412
DATA['NUCLEAR REPULSION ENERGY']['S66-7-monoA-unCP'               ] =      40.41876218
DATA['NUCLEAR REPULSION ENERGY']['S66-7-monoB-unCP'               ] =     180.73873695
DATA['NUCLEAR REPULSION ENERGY']['S66-8-dimer'                    ] =      78.74537406
DATA['NUCLEAR REPULSION ENERGY']['S66-8-monoA-unCP'               ] =      40.42326344
DATA['NUCLEAR REPULSION ENERGY']['S66-8-monoB-unCP'               ] =       9.17236900
DATA['NUCLEAR REPULSION ENERGY']['S66-9-dimer'                    ] =     129.31867271
DATA['NUCLEAR REPULSION ENERGY']['S66-9-monoA-unCP'               ] =      42.10593235
DATA['NUCLEAR REPULSION ENERGY']['S66-9-monoB-unCP'               ] =      40.34710761
DATA['NUCLEAR REPULSION ENERGY']['S66-10-dimer'                   ] =     131.71717765
DATA['NUCLEAR REPULSION ENERGY']['S66-10-monoA-unCP'              ] =      42.09217552
DATA['NUCLEAR REPULSION ENERGY']['S66-10-monoB-unCP'              ] =      42.05982938
DATA['NUCLEAR REPULSION ENERGY']['S66-11-dimer'                   ] =     320.50976921
DATA['NUCLEAR REPULSION ENERGY']['S66-11-monoA-unCP'              ] =      42.09328618
DATA['NUCLEAR REPULSION ENERGY']['S66-11-monoB-unCP'              ] =     180.72211450
DATA['NUCLEAR REPULSION ENERGY']['S66-12-dimer'                   ] =      81.87844165
DATA['NUCLEAR REPULSION ENERGY']['S66-12-monoA-unCP'              ] =      42.04336531
DATA['NUCLEAR REPULSION ENERGY']['S66-12-monoB-unCP'              ] =       9.12312499
DATA['NUCLEAR REPULSION ENERGY']['S66-13-dimer'                   ] =     314.84789007
DATA['NUCLEAR REPULSION ENERGY']['S66-13-monoA-unCP'              ] =     180.80545988
DATA['NUCLEAR REPULSION ENERGY']['S66-13-monoB-unCP'              ] =      40.30378877
DATA['NUCLEAR REPULSION ENERGY']['S66-14-dimer'                   ] =     315.64348724
DATA['NUCLEAR REPULSION ENERGY']['S66-14-monoA-unCP'              ] =     180.81499576
DATA['NUCLEAR REPULSION ENERGY']['S66-14-monoB-unCP'              ] =      42.03791353
DATA['NUCLEAR REPULSION ENERGY']['S66-15-dimer'                   ] =     540.42243680
DATA['NUCLEAR REPULSION ENERGY']['S66-15-monoA-unCP'              ] =     180.53794513
DATA['NUCLEAR REPULSION ENERGY']['S66-15-monoB-unCP'              ] =     180.54327910
DATA['NUCLEAR REPULSION ENERGY']['S66-16-dimer'                   ] =     243.51194018
DATA['NUCLEAR REPULSION ENERGY']['S66-16-monoA-unCP'              ] =     180.57089645
DATA['NUCLEAR REPULSION ENERGY']['S66-16-monoB-unCP'              ] =       9.17374713
DATA['NUCLEAR REPULSION ENERGY']['S66-17-dimer'                   ] =    1040.55250335
DATA['NUCLEAR REPULSION ENERGY']['S66-17-monoA-unCP'              ] =     357.25263911
DATA['NUCLEAR REPULSION ENERGY']['S66-17-monoB-unCP'              ] =     357.22824169
DATA['NUCLEAR REPULSION ENERGY']['S66-18-dimer'                   ] =     269.39653929
DATA['NUCLEAR REPULSION ENERGY']['S66-18-monoA-unCP'              ] =       9.12915636
DATA['NUCLEAR REPULSION ENERGY']['S66-18-monoB-unCP'              ] =     206.28546361
DATA['NUCLEAR REPULSION ENERGY']['S66-19-dimer'                   ] =     337.49486033
DATA['NUCLEAR REPULSION ENERGY']['S66-19-monoA-unCP'              ] =      40.42190801
DATA['NUCLEAR REPULSION ENERGY']['S66-19-monoB-unCP'              ] =     206.28426737
DATA['NUCLEAR REPULSION ENERGY']['S66-20-dimer'                   ] =     381.47467603
DATA['NUCLEAR REPULSION ENERGY']['S66-20-monoA-unCP'              ] =     121.35354216
DATA['NUCLEAR REPULSION ENERGY']['S66-20-monoB-unCP'              ] =     121.35037507
DATA['NUCLEAR REPULSION ENERGY']['S66-21-dimer'                   ] =     373.66110820
DATA['NUCLEAR REPULSION ENERGY']['S66-21-monoA-unCP'              ] =     121.85534909
DATA['NUCLEAR REPULSION ENERGY']['S66-21-monoB-unCP'              ] =     121.85562743
DATA['NUCLEAR REPULSION ENERGY']['S66-22-dimer'                   ] =     685.96293615
DATA['NUCLEAR REPULSION ENERGY']['S66-22-monoA-unCP'              ] =     121.30606379
DATA['NUCLEAR REPULSION ENERGY']['S66-22-monoB-unCP'              ] =     357.30242624
DATA['NUCLEAR REPULSION ENERGY']['S66-23-dimer'                   ] =     682.46450694
DATA['NUCLEAR REPULSION ENERGY']['S66-23-monoA-unCP'              ] =     121.91206440
DATA['NUCLEAR REPULSION ENERGY']['S66-23-monoB-unCP'              ] =     357.16987646
DATA['NUCLEAR REPULSION ENERGY']['S66-24-dimer'                   ] =     623.71187998
DATA['NUCLEAR REPULSION ENERGY']['S66-24-monoA-unCP'              ] =     203.71200257
DATA['NUCLEAR REPULSION ENERGY']['S66-24-monoB-unCP'              ] =     203.71172379
DATA['NUCLEAR REPULSION ENERGY']['S66-25-dimer'                   ] =     637.14156863
DATA['NUCLEAR REPULSION ENERGY']['S66-25-monoA-unCP'              ] =     206.22564193
DATA['NUCLEAR REPULSION ENERGY']['S66-25-monoB-unCP'              ] =     206.22748415
DATA['NUCLEAR REPULSION ENERGY']['S66-26-dimer'                   ] =    1163.54572871
DATA['NUCLEAR REPULSION ENERGY']['S66-26-monoA-unCP'              ] =     357.16027337
DATA['NUCLEAR REPULSION ENERGY']['S66-26-monoB-unCP'              ] =     357.16027370
DATA['NUCLEAR REPULSION ENERGY']['S66-27-dimer'                   ] =     630.67443466
DATA['NUCLEAR REPULSION ENERGY']['S66-27-monoA-unCP'              ] =     203.68422363
DATA['NUCLEAR REPULSION ENERGY']['S66-27-monoB-unCP'              ] =     206.25955744
DATA['NUCLEAR REPULSION ENERGY']['S66-28-dimer'                   ] =     878.32907732
DATA['NUCLEAR REPULSION ENERGY']['S66-28-monoA-unCP'              ] =     203.65134501
DATA['NUCLEAR REPULSION ENERGY']['S66-28-monoB-unCP'              ] =     357.16948119
DATA['NUCLEAR REPULSION ENERGY']['S66-29-dimer'                   ] =     885.28192562
DATA['NUCLEAR REPULSION ENERGY']['S66-29-monoA-unCP'              ] =     206.16040036
DATA['NUCLEAR REPULSION ENERGY']['S66-29-monoB-unCP'              ] =     357.23565563
DATA['NUCLEAR REPULSION ENERGY']['S66-30-dimer'                   ] =     327.62509332
DATA['NUCLEAR REPULSION ENERGY']['S66-30-monoA-unCP'              ] =     203.74228045
DATA['NUCLEAR REPULSION ENERGY']['S66-30-monoB-unCP'              ] =      33.43000301
DATA['NUCLEAR REPULSION ENERGY']['S66-31-dimer'                   ] =     518.26358403
DATA['NUCLEAR REPULSION ENERGY']['S66-31-monoA-unCP'              ] =     357.18726739
DATA['NUCLEAR REPULSION ENERGY']['S66-31-monoB-unCP'              ] =      33.40409180
DATA['NUCLEAR REPULSION ENERGY']['S66-32-dimer'                   ] =     495.33117294
DATA['NUCLEAR REPULSION ENERGY']['S66-32-monoA-unCP'              ] =     357.24995067
DATA['NUCLEAR REPULSION ENERGY']['S66-32-monoB-unCP'              ] =      24.63459975
DATA['NUCLEAR REPULSION ENERGY']['S66-33-dimer'                   ] =     332.11307535
DATA['NUCLEAR REPULSION ENERGY']['S66-33-monoA-unCP'              ] =     206.29228895
DATA['NUCLEAR REPULSION ENERGY']['S66-33-monoB-unCP'              ] =      33.42391806
DATA['NUCLEAR REPULSION ENERGY']['S66-34-dimer'                   ] =     577.94330068
DATA['NUCLEAR REPULSION ENERGY']['S66-34-monoA-unCP'              ] =     185.63664994
DATA['NUCLEAR REPULSION ENERGY']['S66-34-monoB-unCP'              ] =     185.63558546
DATA['NUCLEAR REPULSION ENERGY']['S66-35-dimer'                   ] =     574.13141612
DATA['NUCLEAR REPULSION ENERGY']['S66-35-monoA-unCP'              ] =     185.63471242
DATA['NUCLEAR REPULSION ENERGY']['S66-35-monoB-unCP'              ] =     199.36895747
DATA['NUCLEAR REPULSION ENERGY']['S66-36-dimer'                   ] =     573.01241887
DATA['NUCLEAR REPULSION ENERGY']['S66-36-monoA-unCP'              ] =     199.35493735
DATA['NUCLEAR REPULSION ENERGY']['S66-36-monoB-unCP'              ] =     199.35496470
DATA['NUCLEAR REPULSION ENERGY']['S66-37-dimer'                   ] =     569.42803611
DATA['NUCLEAR REPULSION ENERGY']['S66-37-monoA-unCP'              ] =     188.28929834
DATA['NUCLEAR REPULSION ENERGY']['S66-37-monoB-unCP'              ] =     199.34481507
DATA['NUCLEAR REPULSION ENERGY']['S66-38-dimer'                   ] =     562.36494675
DATA['NUCLEAR REPULSION ENERGY']['S66-38-monoA-unCP'              ] =     188.38358820
DATA['NUCLEAR REPULSION ENERGY']['S66-38-monoB-unCP'              ] =     188.37865241
DATA['NUCLEAR REPULSION ENERGY']['S66-39-dimer'                   ] =     594.82529945
DATA['NUCLEAR REPULSION ENERGY']['S66-39-monoA-unCP'              ] =     203.67735882
DATA['NUCLEAR REPULSION ENERGY']['S66-39-monoB-unCP'              ] =     188.40454306
DATA['NUCLEAR REPULSION ENERGY']['S66-40-dimer'                   ] =     598.08168004
DATA['NUCLEAR REPULSION ENERGY']['S66-40-monoA-unCP'              ] =     203.68538784
DATA['NUCLEAR REPULSION ENERGY']['S66-40-monoB-unCP'              ] =     199.37329650
DATA['NUCLEAR REPULSION ENERGY']['S66-41-dimer'                   ] =     843.32242800
DATA['NUCLEAR REPULSION ENERGY']['S66-41-monoA-unCP'              ] =     357.06617642
DATA['NUCLEAR REPULSION ENERGY']['S66-41-monoB-unCP'              ] =     185.61673585
DATA['NUCLEAR REPULSION ENERGY']['S66-42-dimer'                   ] =     830.51659591
DATA['NUCLEAR REPULSION ENERGY']['S66-42-monoA-unCP'              ] =     357.04169352
DATA['NUCLEAR REPULSION ENERGY']['S66-42-monoB-unCP'              ] =     188.33728572
DATA['NUCLEAR REPULSION ENERGY']['S66-43-dimer'                   ] =     830.36688604
DATA['NUCLEAR REPULSION ENERGY']['S66-43-monoA-unCP'              ] =     357.12713115
DATA['NUCLEAR REPULSION ENERGY']['S66-43-monoB-unCP'              ] =     199.36153551
DATA['NUCLEAR REPULSION ENERGY']['S66-44-dimer'                   ] =     303.64951312
DATA['NUCLEAR REPULSION ENERGY']['S66-44-monoA-unCP'              ] =      33.42556566
DATA['NUCLEAR REPULSION ENERGY']['S66-44-monoB-unCP'              ] =     185.65594848
DATA['NUCLEAR REPULSION ENERGY']['S66-45-dimer'                   ] =     285.69697355
DATA['NUCLEAR REPULSION ENERGY']['S66-45-monoA-unCP'              ] =      24.64923587
DATA['NUCLEAR REPULSION ENERGY']['S66-45-monoB-unCP'              ] =     185.73197134
DATA['NUCLEAR REPULSION ENERGY']['S66-46-dimer'                   ] =     576.36980953
DATA['NUCLEAR REPULSION ENERGY']['S66-46-monoA-unCP'              ] =     180.49044991
DATA['NUCLEAR REPULSION ENERGY']['S66-46-monoB-unCP'              ] =     185.67687994
DATA['NUCLEAR REPULSION ENERGY']['S66-47-dimer'                   ] =     592.90348525
DATA['NUCLEAR REPULSION ENERGY']['S66-47-monoA-unCP'              ] =     203.66921988
DATA['NUCLEAR REPULSION ENERGY']['S66-47-monoB-unCP'              ] =     203.67694204
DATA['NUCLEAR REPULSION ENERGY']['S66-48-dimer'                   ] =     601.34387795
DATA['NUCLEAR REPULSION ENERGY']['S66-48-monoA-unCP'              ] =     206.19608668
DATA['NUCLEAR REPULSION ENERGY']['S66-48-monoB-unCP'              ] =     206.19869697
DATA['NUCLEAR REPULSION ENERGY']['S66-49-dimer'                   ] =     596.54644729
DATA['NUCLEAR REPULSION ENERGY']['S66-49-monoA-unCP'              ] =     203.65045916
DATA['NUCLEAR REPULSION ENERGY']['S66-49-monoB-unCP'              ] =     206.22459403
DATA['NUCLEAR REPULSION ENERGY']['S66-50-dimer'                   ] =     300.96547874
DATA['NUCLEAR REPULSION ENERGY']['S66-50-monoA-unCP'              ] =     203.65156163
DATA['NUCLEAR REPULSION ENERGY']['S66-50-monoB-unCP'              ] =      24.63554547
DATA['NUCLEAR REPULSION ENERGY']['S66-51-dimer'                   ] =      73.51391626
DATA['NUCLEAR REPULSION ENERGY']['S66-51-monoA-unCP'              ] =      24.65072244
DATA['NUCLEAR REPULSION ENERGY']['S66-51-monoB-unCP'              ] =      24.64312912
DATA['NUCLEAR REPULSION ENERGY']['S66-52-dimer'                   ] =     488.72204285
DATA['NUCLEAR REPULSION ENERGY']['S66-52-monoA-unCP'              ] =     203.60587521
DATA['NUCLEAR REPULSION ENERGY']['S66-52-monoB-unCP'              ] =     121.22680816
DATA['NUCLEAR REPULSION ENERGY']['S66-53-dimer'                   ] =     475.54833273
DATA['NUCLEAR REPULSION ENERGY']['S66-53-monoA-unCP'              ] =     203.61290966
DATA['NUCLEAR REPULSION ENERGY']['S66-53-monoB-unCP'              ] =     121.83743933
DATA['NUCLEAR REPULSION ENERGY']['S66-54-dimer'                   ] =     274.02041197
DATA['NUCLEAR REPULSION ENERGY']['S66-54-monoA-unCP'              ] =     203.63390042
DATA['NUCLEAR REPULSION ENERGY']['S66-54-monoB-unCP'              ] =       9.16766818
DATA['NUCLEAR REPULSION ENERGY']['S66-55-dimer'                   ] =     349.34385129
DATA['NUCLEAR REPULSION ENERGY']['S66-55-monoA-unCP'              ] =     203.62143957
DATA['NUCLEAR REPULSION ENERGY']['S66-55-monoB-unCP'              ] =      40.41522246
DATA['NUCLEAR REPULSION ENERGY']['S66-56-dimer'                   ] =     347.25412940
DATA['NUCLEAR REPULSION ENERGY']['S66-56-monoA-unCP'              ] =     203.65859480
DATA['NUCLEAR REPULSION ENERGY']['S66-56-monoB-unCP'              ] =      42.10725315
DATA['NUCLEAR REPULSION ENERGY']['S66-57-dimer'                   ] =     584.88796485
DATA['NUCLEAR REPULSION ENERGY']['S66-57-monoA-unCP'              ] =     203.60060155
DATA['NUCLEAR REPULSION ENERGY']['S66-57-monoB-unCP'              ] =     180.55180987
DATA['NUCLEAR REPULSION ENERGY']['S66-58-dimer'                   ] =     577.23538658
DATA['NUCLEAR REPULSION ENERGY']['S66-58-monoA-unCP'              ] =     206.16864626
DATA['NUCLEAR REPULSION ENERGY']['S66-58-monoB-unCP'              ] =     206.16860003
DATA['NUCLEAR REPULSION ENERGY']['S66-59-dimer'                   ] =      53.29797952
DATA['NUCLEAR REPULSION ENERGY']['S66-59-monoA-unCP'              ] =      24.62604423
DATA['NUCLEAR REPULSION ENERGY']['S66-59-monoB-unCP'              ] =       9.17684034
DATA['NUCLEAR REPULSION ENERGY']['S66-60-dimer'                   ] =     206.60195669
DATA['NUCLEAR REPULSION ENERGY']['S66-60-monoA-unCP'              ] =      24.62574637
DATA['NUCLEAR REPULSION ENERGY']['S66-60-monoB-unCP'              ] =     121.22795347
DATA['NUCLEAR REPULSION ENERGY']['S66-61-dimer'                   ] =     475.00612950
DATA['NUCLEAR REPULSION ENERGY']['S66-61-monoA-unCP'              ] =     185.62492607
DATA['NUCLEAR REPULSION ENERGY']['S66-61-monoB-unCP'              ] =     121.23972648
DATA['NUCLEAR REPULSION ENERGY']['S66-62-dimer'                   ] =     478.48168724
DATA['NUCLEAR REPULSION ENERGY']['S66-62-monoA-unCP'              ] =     185.65184859
DATA['NUCLEAR REPULSION ENERGY']['S66-62-monoB-unCP'              ] =     121.86597939
DATA['NUCLEAR REPULSION ENERGY']['S66-63-dimer'                   ] =     496.78090588
DATA['NUCLEAR REPULSION ENERGY']['S66-63-monoA-unCP'              ] =     203.66095658
DATA['NUCLEAR REPULSION ENERGY']['S66-63-monoB-unCP'              ] =     121.23566219
DATA['NUCLEAR REPULSION ENERGY']['S66-64-dimer'                   ] =     300.38789564
DATA['NUCLEAR REPULSION ENERGY']['S66-64-monoA-unCP'              ] =     180.56185111
DATA['NUCLEAR REPULSION ENERGY']['S66-64-monoB-unCP'              ] =      33.41895147
DATA['NUCLEAR REPULSION ENERGY']['S66-65-dimer'                   ] =     292.14525417
DATA['NUCLEAR REPULSION ENERGY']['S66-65-monoA-unCP'              ] =     206.26607138
DATA['NUCLEAR REPULSION ENERGY']['S66-65-monoB-unCP'              ] =      24.59915901
DATA['NUCLEAR REPULSION ENERGY']['S66-66-dimer'                   ] =     349.09867633
DATA['NUCLEAR REPULSION ENERGY']['S66-66-monoA-unCP'              ] =      42.09376472
DATA['NUCLEAR REPULSION ENERGY']['S66-66-monoB-unCP'              ] =     206.23491680
DATA['NUCLEAR REPULSION ENERGY']['S66-1-monoA-CP'                 ] =       9.15671411
DATA['NUCLEAR REPULSION ENERGY']['S66-1-monoB-CP'                 ] =       9.17259114
DATA['NUCLEAR REPULSION ENERGY']['S66-2-monoA-CP'                 ] =       9.14996836
DATA['NUCLEAR REPULSION ENERGY']['S66-2-monoB-CP'                 ] =      40.29463192
DATA['NUCLEAR REPULSION ENERGY']['S66-3-monoA-CP'                 ] =       9.12565570
DATA['NUCLEAR REPULSION ENERGY']['S66-3-monoB-CP'                 ] =      42.06267577
DATA['NUCLEAR REPULSION ENERGY']['S66-4-monoA-CP'                 ] =       9.13184124
DATA['NUCLEAR REPULSION ENERGY']['S66-4-monoB-CP'                 ] =     180.56084030
DATA['NUCLEAR REPULSION ENERGY']['S66-5-monoA-CP'                 ] =      40.41731272
DATA['NUCLEAR REPULSION ENERGY']['S66-5-monoB-CP'                 ] =      40.29806380
DATA['NUCLEAR REPULSION ENERGY']['S66-6-monoA-CP'                 ] =      40.42467073
DATA['NUCLEAR REPULSION ENERGY']['S66-6-monoB-CP'                 ] =      42.05202847
DATA['NUCLEAR REPULSION ENERGY']['S66-7-monoA-CP'                 ] =      40.41876218
DATA['NUCLEAR REPULSION ENERGY']['S66-7-monoB-CP'                 ] =     180.73873695
DATA['NUCLEAR REPULSION ENERGY']['S66-8-monoA-CP'                 ] =      40.42326344
DATA['NUCLEAR REPULSION ENERGY']['S66-8-monoB-CP'                 ] =       9.17236900
DATA['NUCLEAR REPULSION ENERGY']['S66-9-monoA-CP'                 ] =      42.10593235
DATA['NUCLEAR REPULSION ENERGY']['S66-9-monoB-CP'                 ] =      40.34710761
DATA['NUCLEAR REPULSION ENERGY']['S66-10-monoA-CP'                ] =      42.09217552
DATA['NUCLEAR REPULSION ENERGY']['S66-10-monoB-CP'                ] =      42.05982938
DATA['NUCLEAR REPULSION ENERGY']['S66-11-monoA-CP'                ] =      42.09328618
DATA['NUCLEAR REPULSION ENERGY']['S66-11-monoB-CP'                ] =     180.72211450
DATA['NUCLEAR REPULSION ENERGY']['S66-12-monoA-CP'                ] =      42.04336531
DATA['NUCLEAR REPULSION ENERGY']['S66-12-monoB-CP'                ] =       9.12312499
DATA['NUCLEAR REPULSION ENERGY']['S66-13-monoA-CP'                ] =     180.80545988
DATA['NUCLEAR REPULSION ENERGY']['S66-13-monoB-CP'                ] =      40.30378877
DATA['NUCLEAR REPULSION ENERGY']['S66-14-monoA-CP'                ] =     180.81499576
DATA['NUCLEAR REPULSION ENERGY']['S66-14-monoB-CP'                ] =      42.03791353
DATA['NUCLEAR REPULSION ENERGY']['S66-15-monoA-CP'                ] =     180.53794513
DATA['NUCLEAR REPULSION ENERGY']['S66-15-monoB-CP'                ] =     180.54327910
DATA['NUCLEAR REPULSION ENERGY']['S66-16-monoA-CP'                ] =     180.57089645
DATA['NUCLEAR REPULSION ENERGY']['S66-16-monoB-CP'                ] =       9.17374713
DATA['NUCLEAR REPULSION ENERGY']['S66-17-monoA-CP'                ] =     357.25263911
DATA['NUCLEAR REPULSION ENERGY']['S66-17-monoB-CP'                ] =     357.22824169
DATA['NUCLEAR REPULSION ENERGY']['S66-18-monoA-CP'                ] =       9.12915636
DATA['NUCLEAR REPULSION ENERGY']['S66-18-monoB-CP'                ] =     206.28546361
DATA['NUCLEAR REPULSION ENERGY']['S66-19-monoA-CP'                ] =      40.42190801
DATA['NUCLEAR REPULSION ENERGY']['S66-19-monoB-CP'                ] =     206.28426737
DATA['NUCLEAR REPULSION ENERGY']['S66-20-monoA-CP'                ] =     121.35354216
DATA['NUCLEAR REPULSION ENERGY']['S66-20-monoB-CP'                ] =     121.35037507
DATA['NUCLEAR REPULSION ENERGY']['S66-21-monoA-CP'                ] =     121.85534909
DATA['NUCLEAR REPULSION ENERGY']['S66-21-monoB-CP'                ] =     121.85562743
DATA['NUCLEAR REPULSION ENERGY']['S66-22-monoA-CP'                ] =     121.30606379
DATA['NUCLEAR REPULSION ENERGY']['S66-22-monoB-CP'                ] =     357.30242624
DATA['NUCLEAR REPULSION ENERGY']['S66-23-monoA-CP'                ] =     121.91206440
DATA['NUCLEAR REPULSION ENERGY']['S66-23-monoB-CP'                ] =     357.16987646
DATA['NUCLEAR REPULSION ENERGY']['S66-24-monoA-CP'                ] =     203.71200257
DATA['NUCLEAR REPULSION ENERGY']['S66-24-monoB-CP'                ] =     203.71172379
DATA['NUCLEAR REPULSION ENERGY']['S66-25-monoA-CP'                ] =     206.22564193
DATA['NUCLEAR REPULSION ENERGY']['S66-25-monoB-CP'                ] =     206.22748415
DATA['NUCLEAR REPULSION ENERGY']['S66-26-monoA-CP'                ] =     357.16027337
DATA['NUCLEAR REPULSION ENERGY']['S66-26-monoB-CP'                ] =     357.16027370
DATA['NUCLEAR REPULSION ENERGY']['S66-27-monoA-CP'                ] =     203.68422363
DATA['NUCLEAR REPULSION ENERGY']['S66-27-monoB-CP'                ] =     206.25955744
DATA['NUCLEAR REPULSION ENERGY']['S66-28-monoA-CP'                ] =     203.65134501
DATA['NUCLEAR REPULSION ENERGY']['S66-28-monoB-CP'                ] =     357.16948119
DATA['NUCLEAR REPULSION ENERGY']['S66-29-monoA-CP'                ] =     206.16040036
DATA['NUCLEAR REPULSION ENERGY']['S66-29-monoB-CP'                ] =     357.23565563
DATA['NUCLEAR REPULSION ENERGY']['S66-30-monoA-CP'                ] =     203.74228045
DATA['NUCLEAR REPULSION ENERGY']['S66-30-monoB-CP'                ] =      33.43000301
DATA['NUCLEAR REPULSION ENERGY']['S66-31-monoA-CP'                ] =     357.18726739
DATA['NUCLEAR REPULSION ENERGY']['S66-31-monoB-CP'                ] =      33.40409180
DATA['NUCLEAR REPULSION ENERGY']['S66-32-monoA-CP'                ] =     357.24995067
DATA['NUCLEAR REPULSION ENERGY']['S66-32-monoB-CP'                ] =      24.63459975
DATA['NUCLEAR REPULSION ENERGY']['S66-33-monoA-CP'                ] =     206.29228895
DATA['NUCLEAR REPULSION ENERGY']['S66-33-monoB-CP'                ] =      33.42391806
DATA['NUCLEAR REPULSION ENERGY']['S66-34-monoA-CP'                ] =     185.63664994
DATA['NUCLEAR REPULSION ENERGY']['S66-34-monoB-CP'                ] =     185.63558546
DATA['NUCLEAR REPULSION ENERGY']['S66-35-monoA-CP'                ] =     185.63471242
DATA['NUCLEAR REPULSION ENERGY']['S66-35-monoB-CP'                ] =     199.36895747
DATA['NUCLEAR REPULSION ENERGY']['S66-36-monoA-CP'                ] =     199.35493735
DATA['NUCLEAR REPULSION ENERGY']['S66-36-monoB-CP'                ] =     199.35496470
DATA['NUCLEAR REPULSION ENERGY']['S66-37-monoA-CP'                ] =     188.28929834
DATA['NUCLEAR REPULSION ENERGY']['S66-37-monoB-CP'                ] =     199.34481507
DATA['NUCLEAR REPULSION ENERGY']['S66-38-monoA-CP'                ] =     188.38358820
DATA['NUCLEAR REPULSION ENERGY']['S66-38-monoB-CP'                ] =     188.37865241
DATA['NUCLEAR REPULSION ENERGY']['S66-39-monoA-CP'                ] =     203.67735882
DATA['NUCLEAR REPULSION ENERGY']['S66-39-monoB-CP'                ] =     188.40454306
DATA['NUCLEAR REPULSION ENERGY']['S66-40-monoA-CP'                ] =     203.68538784
DATA['NUCLEAR REPULSION ENERGY']['S66-40-monoB-CP'                ] =     199.37329650
DATA['NUCLEAR REPULSION ENERGY']['S66-41-monoA-CP'                ] =     357.06617642
DATA['NUCLEAR REPULSION ENERGY']['S66-41-monoB-CP'                ] =     185.61673585
DATA['NUCLEAR REPULSION ENERGY']['S66-42-monoA-CP'                ] =     357.04169352
DATA['NUCLEAR REPULSION ENERGY']['S66-42-monoB-CP'                ] =     188.33728572
DATA['NUCLEAR REPULSION ENERGY']['S66-43-monoA-CP'                ] =     357.12713115
DATA['NUCLEAR REPULSION ENERGY']['S66-43-monoB-CP'                ] =     199.36153551
DATA['NUCLEAR REPULSION ENERGY']['S66-44-monoA-CP'                ] =      33.42556566
DATA['NUCLEAR REPULSION ENERGY']['S66-44-monoB-CP'                ] =     185.65594848
DATA['NUCLEAR REPULSION ENERGY']['S66-45-monoA-CP'                ] =      24.64923587
DATA['NUCLEAR REPULSION ENERGY']['S66-45-monoB-CP'                ] =     185.73197134
DATA['NUCLEAR REPULSION ENERGY']['S66-46-monoA-CP'                ] =     180.49044991
DATA['NUCLEAR REPULSION ENERGY']['S66-46-monoB-CP'                ] =     185.67687994
DATA['NUCLEAR REPULSION ENERGY']['S66-47-monoA-CP'                ] =     203.66921988
DATA['NUCLEAR REPULSION ENERGY']['S66-47-monoB-CP'                ] =     203.67694204
DATA['NUCLEAR REPULSION ENERGY']['S66-48-monoA-CP'                ] =     206.19608668
DATA['NUCLEAR REPULSION ENERGY']['S66-48-monoB-CP'                ] =     206.19869697
DATA['NUCLEAR REPULSION ENERGY']['S66-49-monoA-CP'                ] =     203.65045916
DATA['NUCLEAR REPULSION ENERGY']['S66-49-monoB-CP'                ] =     206.22459403
DATA['NUCLEAR REPULSION ENERGY']['S66-50-monoA-CP'                ] =     203.65156163
DATA['NUCLEAR REPULSION ENERGY']['S66-50-monoB-CP'                ] =      24.63554547
DATA['NUCLEAR REPULSION ENERGY']['S66-51-monoA-CP'                ] =      24.65072244
DATA['NUCLEAR REPULSION ENERGY']['S66-51-monoB-CP'                ] =      24.64312912
DATA['NUCLEAR REPULSION ENERGY']['S66-52-monoA-CP'                ] =     203.60587521
DATA['NUCLEAR REPULSION ENERGY']['S66-52-monoB-CP'                ] =     121.22680816
DATA['NUCLEAR REPULSION ENERGY']['S66-53-monoA-CP'                ] =     203.61290966
DATA['NUCLEAR REPULSION ENERGY']['S66-53-monoB-CP'                ] =     121.83743933
DATA['NUCLEAR REPULSION ENERGY']['S66-54-monoA-CP'                ] =     203.63390042
DATA['NUCLEAR REPULSION ENERGY']['S66-54-monoB-CP'                ] =       9.16766818
DATA['NUCLEAR REPULSION ENERGY']['S66-55-monoA-CP'                ] =     203.62143957
DATA['NUCLEAR REPULSION ENERGY']['S66-55-monoB-CP'                ] =      40.41522246
DATA['NUCLEAR REPULSION ENERGY']['S66-56-monoA-CP'                ] =     203.65859480
DATA['NUCLEAR REPULSION ENERGY']['S66-56-monoB-CP'                ] =      42.10725315
DATA['NUCLEAR REPULSION ENERGY']['S66-57-monoA-CP'                ] =     203.60060155
DATA['NUCLEAR REPULSION ENERGY']['S66-57-monoB-CP'                ] =     180.55180987
DATA['NUCLEAR REPULSION ENERGY']['S66-58-monoA-CP'                ] =     206.16864626
DATA['NUCLEAR REPULSION ENERGY']['S66-58-monoB-CP'                ] =     206.16860003
DATA['NUCLEAR REPULSION ENERGY']['S66-59-monoA-CP'                ] =      24.62604423
DATA['NUCLEAR REPULSION ENERGY']['S66-59-monoB-CP'                ] =       9.17684034
DATA['NUCLEAR REPULSION ENERGY']['S66-60-monoA-CP'                ] =      24.62574637
DATA['NUCLEAR REPULSION ENERGY']['S66-60-monoB-CP'                ] =     121.22795347
DATA['NUCLEAR REPULSION ENERGY']['S66-61-monoA-CP'                ] =     185.62492607
DATA['NUCLEAR REPULSION ENERGY']['S66-61-monoB-CP'                ] =     121.23972648
DATA['NUCLEAR REPULSION ENERGY']['S66-62-monoA-CP'                ] =     185.65184859
DATA['NUCLEAR REPULSION ENERGY']['S66-62-monoB-CP'                ] =     121.86597939
DATA['NUCLEAR REPULSION ENERGY']['S66-63-monoA-CP'                ] =     203.66095658
DATA['NUCLEAR REPULSION ENERGY']['S66-63-monoB-CP'                ] =     121.23566219
DATA['NUCLEAR REPULSION ENERGY']['S66-64-monoA-CP'                ] =     180.56185111
DATA['NUCLEAR REPULSION ENERGY']['S66-64-monoB-CP'                ] =      33.41895147
DATA['NUCLEAR REPULSION ENERGY']['S66-65-monoA-CP'                ] =     206.26607138
DATA['NUCLEAR REPULSION ENERGY']['S66-65-monoB-CP'                ] =      24.59915901
DATA['NUCLEAR REPULSION ENERGY']['S66-66-monoA-CP'                ] =      42.09376472
DATA['NUCLEAR REPULSION ENERGY']['S66-66-monoB-CP'                ] =     206.23491680
