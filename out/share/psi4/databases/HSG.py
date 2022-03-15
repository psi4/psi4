#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2021 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

"""
| Database (Merz) of interaction energies for bimolecular complexes from protein-indinavir reaction site.
| Geometries from and original reference energies from Faver et al. JCTC 7 790 (2011).
| Revised reference interaction energies (HSGA) from Marshall et al. JCP 135 194102 (2011).

- **cp**  ``'off'`` || ``'on'``

- **rlxd** ``'off'``

- **benchmark**

  - ``'HSG0'`` Faver et al. JCTC 7 790 (2011).
  - |dl| ``'HSGA'`` |dr| Marshall et al. JCP 135 194102 (2011).

- **subset**

  - ``'small'``
  - ``'large'``

"""
import qcdb

# <<< HSG Database Module >>>
dbse = 'HSG'

# <<< Database Memobers >>>
HRXN = range(1, 22)
HRXN_SM = [6, 15]
HRXN_LG = [14]

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV_CP = {}  # order of active reagents per counterpoise-corrected reaction
ACTV_SA = {}  # order of active reagents for non-supramolecular calculations
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

# <<< Reference Values >>>
BIND = {}
# Original publication
BIND_HSG0 = {}
BIND_HSG0['%s-%s' % (dbse,  1)] =   -0.519
BIND_HSG0['%s-%s' % (dbse,  2)] =   -2.181
BIND_HSG0['%s-%s' % (dbse,  3)] =   -2.451
BIND_HSG0['%s-%s' % (dbse,  4)] =  -16.445
BIND_HSG0['%s-%s' % (dbse,  5)] =  -18.984
BIND_HSG0['%s-%s' % (dbse,  6)] =   -6.009
BIND_HSG0['%s-%s' % (dbse,  7)] =   -3.301
BIND_HSG0['%s-%s' % (dbse,  8)] =   -0.554
BIND_HSG0['%s-%s' % (dbse,  9)] =   -5.038
BIND_HSG0['%s-%s' % (dbse, 10)] =   -7.532
BIND_HSG0['%s-%s' % (dbse, 11)] =   -6.279
BIND_HSG0['%s-%s' % (dbse, 12)] =    0.305
BIND_HSG0['%s-%s' % (dbse, 13)] =   -2.087
BIND_HSG0['%s-%s' % (dbse, 14)] =   -1.376
BIND_HSG0['%s-%s' % (dbse, 15)] =   -0.853
BIND_HSG0['%s-%s' % (dbse, 16)] =   -1.097
BIND_HSG0['%s-%s' % (dbse, 17)] =   -1.504
BIND_HSG0['%s-%s' % (dbse, 18)] =   -0.473
BIND_HSG0['%s-%s' % (dbse, 19)] =   -1.569
BIND_HSG0['%s-%s' % (dbse, 20)] =    0.391
BIND_HSG0['%s-%s' % (dbse, 21)] =   -9.486
# Current revision
BIND_HSGA = {}
BIND_HSGA['%s-%s' % (dbse,  1)] =   -0.518
BIND_HSGA['%s-%s' % (dbse,  2)] =   -2.283
BIND_HSGA['%s-%s' % (dbse,  3)] =   -2.478
BIND_HSGA['%s-%s' % (dbse,  4)] =  -16.526
BIND_HSGA['%s-%s' % (dbse,  5)] =  -19.076
BIND_HSGA['%s-%s' % (dbse,  6)] =   -5.998
BIND_HSGA['%s-%s' % (dbse,  7)] =   -3.308
BIND_HSGA['%s-%s' % (dbse,  8)] =   -0.581
BIND_HSGA['%s-%s' % (dbse,  9)] =   -5.066
BIND_HSGA['%s-%s' % (dbse, 10)] =   -7.509
BIND_HSGA['%s-%s' % (dbse, 11)] =   -6.274
BIND_HSGA['%s-%s' % (dbse, 12)] =    0.302
BIND_HSGA['%s-%s' % (dbse, 13)] =   -2.103
BIND_HSGA['%s-%s' % (dbse, 14)] =   -1.378
BIND_HSGA['%s-%s' % (dbse, 15)] =   -0.856
BIND_HSGA['%s-%s' % (dbse, 16)] =   -1.100
BIND_HSGA['%s-%s' % (dbse, 17)] =   -1.534
BIND_HSGA['%s-%s' % (dbse, 18)] =   -0.472
BIND_HSGA['%s-%s' % (dbse, 19)] =   -1.598
BIND_HSGA['%s-%s' % (dbse, 20)] =    0.378
BIND_HSGA['%s-%s' % (dbse, 21)] =   -9.538
# Set default
BIND = BIND_HSGA

# <<< Coment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse,  1)] = 'ala29-big'
TAGL['%s-%s-dimer'      % (dbse,  1)] = 'ala29-big'
TAGL['%s-%s-monoA-CP'   % (dbse,  1)] = 'indinavir from ala29-big'
TAGL['%s-%s-monoB-CP'   % (dbse,  1)] = 'alanine from ala29-big'
TAGL['%s-%s-monoA-unCP' % (dbse,  1)] = 'indinavir from ala29-big'
TAGL['%s-%s-monoB-unCP' % (dbse,  1)] = 'alanine from ala29-big'
TAGL['%s-%s'            % (dbse,  2)] = 'ala128-small'
TAGL['%s-%s-dimer'      % (dbse,  2)] = 'ala128-small'
TAGL['%s-%s-monoA-CP'   % (dbse,  2)] = 'alanine from ala128-small'
TAGL['%s-%s-monoB-CP'   % (dbse,  2)] = 'indinavir from ala128-small'
TAGL['%s-%s-monoA-unCP' % (dbse,  2)] = 'alanine from ala128-small'
TAGL['%s-%s-monoB-unCP' % (dbse,  2)] = 'indinavir from ala128-small'
TAGL['%s-%s'            % (dbse,  3)] = 'arg8'
TAGL['%s-%s-dimer'      % (dbse,  3)] = 'arg8'
TAGL['%s-%s-monoA-CP'   % (dbse,  3)] = 'arginine from arg8'
TAGL['%s-%s-monoB-CP'   % (dbse,  3)] = 'indinavir from arg8'
TAGL['%s-%s-monoA-unCP' % (dbse,  3)] = 'arginine from arg8'
TAGL['%s-%s-monoB-unCP' % (dbse,  3)] = 'indinavir from arg8'
TAGL['%s-%s'            % (dbse,  4)] = 'ash26-asp125'
TAGL['%s-%s-dimer'      % (dbse,  4)] = 'ash26-asp125'
TAGL['%s-%s-monoA-CP'   % (dbse,  4)] = 'aspartic acid from ash26-asp125'
TAGL['%s-%s-monoB-CP'   % (dbse,  4)] = 'indinavir from ash26-asp125'
TAGL['%s-%s-monoA-unCP' % (dbse,  4)] = 'aspartic acid from ash26-asp125'
TAGL['%s-%s-monoB-unCP' % (dbse,  4)] = 'indinavir from ash26-asp125'
TAGL['%s-%s'            % (dbse,  5)] = 'asp129-big'
TAGL['%s-%s-dimer'      % (dbse,  5)] = 'asp129-big'
TAGL['%s-%s-monoA-CP'   % (dbse,  5)] = 'aspartic acid from asp129-big'
TAGL['%s-%s-monoB-CP'   % (dbse,  5)] = 'indinavir from asp129-big'
TAGL['%s-%s-monoA-unCP' % (dbse,  5)] = 'aspartic acid from asp129-big'
TAGL['%s-%s-monoB-unCP' % (dbse,  5)] = 'indinavir from asp129-big'
TAGL['%s-%s'            % (dbse,  6)] = 'asp130'
TAGL['%s-%s-dimer'      % (dbse,  6)] = 'asp130'
TAGL['%s-%s-monoA-CP'   % (dbse,  6)] = 'aspartic acid from asp130'
TAGL['%s-%s-monoB-CP'   % (dbse,  6)] = 'indinavir from asp130'
TAGL['%s-%s-monoA-unCP' % (dbse,  6)] = 'aspartic acid from asp130'
TAGL['%s-%s-monoB-unCP' % (dbse,  6)] = 'indinavir from asp130'
TAGL['%s-%s'            % (dbse,  7)] = 'gly28-big'
TAGL['%s-%s-dimer'      % (dbse,  7)] = 'gly28-big'
TAGL['%s-%s-monoA-CP'   % (dbse,  7)] = 'glycine from gly28-big'
TAGL['%s-%s-monoB-CP'   % (dbse,  7)] = 'indinavir from gly28-big'
TAGL['%s-%s-monoA-unCP' % (dbse,  7)] = 'glycine from gly28-big'
TAGL['%s-%s-monoB-unCP' % (dbse,  7)] = 'indinavir from gly28-big'
TAGL['%s-%s'            % (dbse,  8)] = 'gly50-ring-big'
TAGL['%s-%s-dimer'      % (dbse,  8)] = 'gly50-ring-big'
TAGL['%s-%s-monoA-CP'   % (dbse,  8)] = 'glycine from gly50-ring-big'
TAGL['%s-%s-monoB-CP'   % (dbse,  8)] = 'indinavir from gly50-ring-big'
TAGL['%s-%s-monoA-unCP' % (dbse,  8)] = 'glycine from gly50-ring-big'
TAGL['%s-%s-monoB-unCP' % (dbse,  8)] = 'indinavir from gly50-ring-big'
TAGL['%s-%s'            % (dbse,  9)] = 'gly50-v1'
TAGL['%s-%s-dimer'      % (dbse,  9)] = 'gly50-v1'
TAGL['%s-%s-monoA-CP'   % (dbse,  9)] = 'glycine from gly50-v1'
TAGL['%s-%s-monoB-CP'   % (dbse,  9)] = 'indinavir from gly50-v1'
TAGL['%s-%s-monoA-unCP' % (dbse,  9)] = 'glycine from gly50-v1'
TAGL['%s-%s-monoB-unCP' % (dbse,  9)] = 'indinavir from gly50-v1'
TAGL['%s-%s'            % (dbse, 10)] = 'gly127'
TAGL['%s-%s-dimer'      % (dbse, 10)] = 'gly127'
TAGL['%s-%s-monoA-CP'   % (dbse, 10)] = 'indinavir from gly127'
TAGL['%s-%s-monoB-CP'   % (dbse, 10)] = 'glycine from gly127'
TAGL['%s-%s-monoA-unCP' % (dbse, 10)] = 'indinavir from gly127'
TAGL['%s-%s-monoB-unCP' % (dbse, 10)] = 'glycine from gly127'
TAGL['%s-%s'            % (dbse, 11)] = 'gly148'
TAGL['%s-%s-dimer'      % (dbse, 11)] = 'gly148'
TAGL['%s-%s-monoA-CP'   % (dbse, 11)] = 'glycine from gly148'
TAGL['%s-%s-monoB-CP'   % (dbse, 11)] = 'indinavir from gly148'
TAGL['%s-%s-monoA-unCP' % (dbse, 11)] = 'glycine from gly148'
TAGL['%s-%s-monoB-unCP' % (dbse, 11)] = 'indinavir from gly148'
TAGL['%s-%s'            % (dbse, 12)] = 'ile48-big'
TAGL['%s-%s-dimer'      % (dbse, 12)] = 'ile48-big'
TAGL['%s-%s-monoA-CP'   % (dbse, 12)] = 'isoleucine from ile48-big'
TAGL['%s-%s-monoB-CP'   % (dbse, 12)] = 'indinavir from ile48-big'
TAGL['%s-%s-monoA-unCP' % (dbse, 12)] = 'isoleucine from ile48-big'
TAGL['%s-%s-monoB-unCP' % (dbse, 12)] = 'indinavir from ile48-big'
TAGL['%s-%s'            % (dbse, 13)] = 'ile147'
TAGL['%s-%s-dimer'      % (dbse, 13)] = 'ile147'
TAGL['%s-%s-monoA-CP'   % (dbse, 13)] = 'isoleucine from ile147'
TAGL['%s-%s-monoB-CP'   % (dbse, 13)] = 'indinavir from ile147'
TAGL['%s-%s-monoA-unCP' % (dbse, 13)] = 'isoleucine from ile147'
TAGL['%s-%s-monoB-unCP' % (dbse, 13)] = 'indinavir from ile147'
TAGL['%s-%s'            % (dbse, 14)] = 'ile150-big'
TAGL['%s-%s-dimer'      % (dbse, 14)] = 'ile150-big'
TAGL['%s-%s-monoA-CP'   % (dbse, 14)] = 'isoleucine from ile150-big'
TAGL['%s-%s-monoB-CP'   % (dbse, 14)] = 'indinavir from ile150-big'
TAGL['%s-%s-monoA-unCP' % (dbse, 14)] = 'isoleucine from ile150-big'
TAGL['%s-%s-monoB-unCP' % (dbse, 14)] = 'indinavir from ile150-big'
TAGL['%s-%s'            % (dbse, 15)] = 'ile184'
TAGL['%s-%s-dimer'      % (dbse, 15)] = 'ile184'
TAGL['%s-%s-monoA-CP'   % (dbse, 15)] = 'isoleucine from ile184'
TAGL['%s-%s-monoB-CP'   % (dbse, 15)] = 'indinavir from ile184'
TAGL['%s-%s-monoA-unCP' % (dbse, 15)] = 'isoleucine from ile184'
TAGL['%s-%s-monoB-unCP' % (dbse, 15)] = 'indinavir from ile184'
TAGL['%s-%s'            % (dbse, 16)] = 'leu23-big'
TAGL['%s-%s-dimer'      % (dbse, 16)] = 'leu23-big'
TAGL['%s-%s-monoA-CP'   % (dbse, 16)] = 'leucine from leu23-big'
TAGL['%s-%s-monoB-CP'   % (dbse, 16)] = 'indinavir from leu23-big'
TAGL['%s-%s-monoA-unCP' % (dbse, 16)] = 'leucine from leu23-big'
TAGL['%s-%s-monoB-unCP' % (dbse, 16)] = 'indinavir from leu23-big'
TAGL['%s-%s'            % (dbse, 17)] = 'pro181'
TAGL['%s-%s-dimer'      % (dbse, 17)] = 'pro181'
TAGL['%s-%s-monoA-CP'   % (dbse, 17)] = 'proline from pro181'
TAGL['%s-%s-monoB-CP'   % (dbse, 17)] = 'indinavir from pro181'
TAGL['%s-%s-monoA-unCP' % (dbse, 17)] = 'proline from pro181'
TAGL['%s-%s-monoB-unCP' % (dbse, 17)] = 'indinavir from pro181'
TAGL['%s-%s'            % (dbse, 18)] = 'val33-big'
TAGL['%s-%s-dimer'      % (dbse, 18)] = 'val33-big'
TAGL['%s-%s-monoA-CP'   % (dbse, 18)] = 'valine from val33-big'
TAGL['%s-%s-monoB-CP'   % (dbse, 18)] = 'indinavir from val33-big'
TAGL['%s-%s-monoA-unCP' % (dbse, 18)] = 'valine from val33-big'
TAGL['%s-%s-monoB-unCP' % (dbse, 18)] = 'indinavir from val33-big'
TAGL['%s-%s'            % (dbse, 19)] = 'val83'
TAGL['%s-%s-dimer'      % (dbse, 19)] = 'val83'
TAGL['%s-%s-monoA-CP'   % (dbse, 19)] = 'valine from val83'
TAGL['%s-%s-monoB-CP'   % (dbse, 19)] = 'indinavir from val83'
TAGL['%s-%s-monoA-unCP' % (dbse, 19)] = 'valine from val83'
TAGL['%s-%s-monoB-unCP' % (dbse, 19)] = 'indinavir from val83'
TAGL['%s-%s'            % (dbse, 20)] = 'val132'
TAGL['%s-%s-dimer'      % (dbse, 20)] = 'val132'
TAGL['%s-%s-monoA-CP'   % (dbse, 20)] = 'valine from val132'
TAGL['%s-%s-monoB-CP'   % (dbse, 20)] = 'indinavir from val132'
TAGL['%s-%s-monoA-unCP' % (dbse, 20)] = 'valine from val132'
TAGL['%s-%s-monoB-unCP' % (dbse, 20)] = 'indinavir from val132'
TAGL['%s-%s'            % (dbse, 21)] = 'wat200'
TAGL['%s-%s-dimer'      % (dbse, 21)] = 'wat200'
TAGL['%s-%s-monoA-CP'   % (dbse, 21)] = 'water from wat200'
TAGL['%s-%s-monoB-CP'   % (dbse, 21)] = 'indinavir from wat200'
TAGL['%s-%s-monoA-unCP' % (dbse, 21)] = 'water from wat200'
TAGL['%s-%s-monoB-unCP' % (dbse, 21)] = 'indinavir from wat200'

# <<< Geometry Specification Strings >>>
GEOS = {}

GEOS['%s-%s-dimer' % (dbse, '1')] = qcdb.Molecule("""
0 1
C   13.03200       29.07900       6.986000
H   12.30800       29.25100       7.790000
H   13.47200       28.08100       7.080000
H   13.82700       29.84100       7.035000
H   12.50772       29.16746       6.023030
--
0 1
C   10.60200       24.81800       6.466000
O   10.95600       23.84000       7.103000
N   10.17800       25.94300       7.070000
C   10.09100       26.25600       8.476000
C   9.372000       27.59000       8.640000
C   11.44600       26.35600       9.091000
C   9.333000       25.25000       9.282000
H   9.874000       26.68900       6.497000
H   9.908000       28.37100       8.093000
H   8.364000       27.46400       8.233000
H   9.317000       27.84600       9.706000
H   9.807000       24.28200       9.160000
H   9.371000       25.57400       10.32900
H   8.328000       25.26700       8.900000
H   11.28800       26.57600       10.14400
H   11.97000       27.14900       8.585000
H   11.93200       25.39300       8.957000
H   10.61998       24.85900       5.366911
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '2')] = qcdb.Molecule("""
0 1
C   18.71400       22.19500       2.742000
H   18.37900       21.58700       3.577000
C   17.68800       22.11500       1.586000
H   17.61600       21.07600       1.227000
H   16.69600       22.44400       1.940000
H   18.00000       22.76400       0.747000
H   18.77673       23.23495       3.094948
H   19.70954       21.82087       2.461043
--
0 1
C   16.65000       19.51600       5.550000
C   18.05900       18.95300       5.945000
O   19.08300       19.90900       5.610000
C   18.24000       17.70600       5.049000
C   17.41000       17.97200       3.810000
C   17.48600       17.37200       2.527000
C   16.60600       17.86100       1.547000
C   15.70300       18.89200       1.854000
C   15.65400       19.47000       3.123000
C   16.51300       18.99700       4.107000
H   15.85600       19.11800       6.209000
H   18.17800       18.65700       6.992000
H   17.85200       16.83600       5.584000
H   19.29200       17.55300       4.850000
H   18.19800       16.55200       2.242000
H   16.61400       17.45600       0.523000
H   15.03200       19.27500       1.092000
H   14.96900       20.29200       3.327000
H   19.96600       19.54000       5.876000
H   16.54908       20.59917       5.712962
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '3')] = qcdb.Molecule("""
1 1
C   23.73500       21.90400       8.645000
H   24.33900       21.64000       9.515000
H   23.04400       22.70000       8.945000
N   22.96900       20.71700       8.247000
H   22.85100       20.56200       7.249000
C   22.40300       19.85100       9.070000
N   22.40300       20.05000       10.36600
H   22.82700       20.88300       10.74700
H   21.93200       19.41900       10.98100
N   21.82000       18.77600       8.615000
H   21.76800       18.61500       7.604000
H   21.33700       18.14400       9.214000
H   24.38340       22.25566       7.828969
--
0 1
C   16.73700       21.75300       8.985000
C   18.06300       21.93100       8.570000
C   19.04400       21.01900       8.966000
C   18.68400       19.93600       9.775000
C   17.35900       19.76800       10.19300
C   16.38200       20.67900       9.796000
H   15.33000       20.56400       10.09500
H   17.07400       18.92500       10.82100
H   19.43700       19.21300       10.07200
H   20.08100       21.14800       8.627000
H   18.32800       22.76900       7.913000
H   15.93631       22.42849       8.649437
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '4')] = qcdb.Molecule("""
-1 1
C   17.05600       28.65300       6.834000
H   17.72900       28.22900       7.569000
H   16.32100       29.27500       7.342000
C   16.35100       27.45400       6.256000
O   16.17800       26.43900       6.902000
O   15.98200       27.55700       4.965000
H   15.73800       26.67800       4.650000
C   16.27300       25.57900       0.088000
H   16.75700       24.66100      -0.278000
H   15.39700       25.75100      -0.577000
C   15.87600       25.38300       1.569000
O   16.42900       26.07300       2.466000
O   14.98200       24.56700       1.861000
H   17.61665       29.26091       6.108662
H   16.97158       26.42544       0.013713047
--
0 1
C   14.25800       24.02900       5.093000
O   15.51000       24.53800       4.641000
H   15.42000       24.70300       3.667000
H   14.02700       23.02800       4.754000
H   13.45976       24.69373       4.731124
H   14.36576       23.85731       6.174161
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '5')] = qcdb.Molecule("""
-1 1
C   18.71400       22.19500       2.742000
H   18.37900       21.58700       3.577000
C   20.10300       21.67300       2.350000
O   20.78600       22.26500       1.513000
N   20.55800       20.55500       2.927000
H   20.07200       20.11800       3.686000
C   21.79700       19.92300       2.527000
H   22.55800       20.68900       2.504000
C   22.17900       18.80900       3.507000
H   21.42700       18.01600       3.405000
H   23.13400       18.39700       3.140000
C   22.26300       19.27700       4.986000
O   23.05600       20.18800       5.322000
O   21.52500       18.75500       5.855000
H   21.73715       19.49687       1.514662
H   18.77673       23.23495       3.094948
H   17.98479       22.13814       1.920400
--
0 1
C   18.05900       18.95300       5.945000
H   18.17800       18.65700       6.992000
O   19.08300       19.90900       5.610000
H   19.96600       19.54000       5.876000
H   17.07047       19.34799       5.667876
H   18.18777       18.06583       5.307547
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '6')] = qcdb.Molecule("""
-1 1
C   20.08800       16.66100      -0.398000
H   19.04500       16.60200      -0.082000
H   20.18300       16.24400      -1.387000
C   20.90200       15.70700       0.505000
O   21.91000       15.13700       0.031000
O   20.37700       15.22700       1.542000
H   20.29988       17.73389      -0.5163442
--
0 1
C   17.41000       17.97200       3.810000
C   17.48600       17.37200       2.527000
C   16.60600       17.86100       1.547000
C   15.70300       18.89200       1.854000
C   15.65400       19.47000       3.123000
C   16.51300       18.99700       4.107000
H   18.19800       16.55200       2.242000
H   16.61400       17.45600       0.523000
H   15.03200       19.27500       1.092000
H   14.96900       20.29200       3.327000
H   16.61088       19.36781       5.137980
H   18.01270       17.77884       4.709692
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '7')] = qcdb.Molecule("""
0 1
C   15.27800       30.38900       2.305000
O   14.37900       30.99200       1.712000
N   15.20400       29.08900       2.640000
H   15.90100       28.63400       3.249000
C   14.02400       28.31500       2.332000
H   13.65700       28.66700       1.375000
H   14.32500       27.27400       2.236000
C   12.93200       28.48000       3.398000
O   11.76000       28.20200       3.138000
N   13.27200       28.98600       4.593000
H   14.23600       29.23200       4.806000
H   12.53585       29.15009       5.393723
H   16.18758       30.94244       2.581356
--
0 1
C   14.25800       24.02900       5.093000
O   15.51000       24.53800       4.641000
C   13.12200       24.97500       4.578000
N   11.86200       24.24500       4.359000
C   11.78200       23.78100       2.948000
C   10.62700       24.87500       4.938000
H   15.42000       24.70300       3.667000
H   14.02700       23.02800       4.754000
H   12.93300       25.79400       5.257000
H   13.40100       25.40500       3.596000
H   10.49500       25.96800       4.694000
H   12.00900       24.60600       2.286000
H   12.53400       23.01100       2.841000
H   10.60902       24.83400       6.037089
H   9.792144       24.23767       4.611163
H   10.82891       23.34622       2.612456
H   14.36576       23.85731       6.174161
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '8')] = qcdb.Molecule("""
0 1
C   4.290000       24.10300       10.08600
O   4.354000       23.13400       10.84900
N   4.507000       24.01000       8.777000
H   4.471000       24.81400       8.166000
C   4.858000       22.75500       8.154000
H   4.799000       21.96700       8.885000
H   4.136000       22.58200       7.364000
C   6.276000       22.76600       7.611000
O   6.647000       23.71300       6.921000
N   7.066000       21.74000       7.938000
H   6.754000       21.01200       8.573000
H   4.054781       25.09671       10.49492
H   8.028304       21.54906       7.440493
--
0 1
C   6.218000       24.82200       2.171000
C   6.715000       23.76400       2.930000
C   5.918000       23.22600       3.934000
C   4.663000       23.78000       4.166000
C   4.226000       24.84700       3.389000
N   4.990000       25.38300       2.391000
H   6.798000       25.19700       1.318000
H   3.211000       25.24700       3.575000
H   4.018000       23.38100       4.943000
H   6.274000       22.37100       4.504000
H   7.673767       23.27955       2.693205
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '9')] = qcdb.Molecule("""
0 1
N   9.334000       19.62700       6.121000
H   10.08000       19.74000       6.781000
C   8.252000       20.39100       6.228000
O   7.244000       20.21700       5.545000
C   8.331000       21.48900       7.284000
H   8.671000       22.40300       6.814000
H   9.046000       21.19200       8.047000
N   7.066000       21.74000       7.938000
H   6.754000       21.01200       8.573000
O   6.647000       23.71300       6.921000
C   6.276000       22.76600       7.611000
H   9.495216       18.80266       5.410733
H   5.248769       22.75803       8.004361
--
0 1
C   10.62700       24.87500       4.938000
C   10.60200       24.81800       6.466000
O   10.95600       23.84000       7.103000
N   10.17800       25.94300       7.070000
C   10.09100       26.25600       8.476000
H   9.874000       26.68900       6.497000
H   10.49500       25.96800       4.694000
H   11.53119       24.41376       4.514093
H   9.792144       24.23767       4.611163
H   9.533425       25.51600       9.068883
H   9.572130       27.21869       8.594352
H   11.09040       26.32976       8.929603
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '10')] = qcdb.Molecule("""
0 1
C   18.71400       22.19500       2.742000
H   18.37900       21.58700       3.577000
N   18.79700       23.57100       3.209000
H   18.86300       24.28300       2.488000
C   18.84100       23.91000       4.508000
O   18.83700       23.07300       5.410000
C   18.95700       25.38900       4.894000
H   18.06200       25.68500       5.433000
H   19.81900       25.52500       5.542000
H   19.09001       26.05938       4.032084
H   19.70954       21.82087       2.461043
H   17.98479       22.13814       1.920400
--
0 1
C   16.65000       19.51600       5.550000
H   15.85600       19.11800       6.209000
H   17.38100       21.50400       5.804000
N   16.51500       20.96500       5.768000
C   15.41800       21.48800       6.353000
O   14.37000       20.87200       6.467000
C   15.57100       22.84700       6.986000
H   16.50000       23.29500       6.652000
H   14.74390       23.51764       6.710063
H   15.61789       22.72907       8.078654
H   17.63853       19.12101       5.827124
H   16.55212       19.14519       4.519021
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '11')] = qcdb.Molecule("""
0 1
C   14.75600       16.28200       8.071000
O   15.85100       16.84300       8.024000
N   13.59000       16.92100       8.167000
H   12.70700       16.41600       8.193000
C   13.49400       18.33500       8.443000
H   14.48600       18.75400       8.594000
H   13.03600       18.80900       7.582000
C   12.63300       18.57700       9.678000
O   12.60000       17.78900       10.62400
N   11.87400       19.66100       9.642000
H   11.86900       20.24000       8.807000
H   14.71385       15.18329       8.038301
H   11.09386       19.92029       10.37286
--
0 1
C   16.65000       19.51600       5.550000
N   16.51500       20.96500       5.768000
H   17.38100       21.50400       5.804000
H   15.85600       19.11800       6.209000
C   15.41800       21.48800       6.353000
O   14.37000       20.87200       6.467000
C   15.57100       22.84700       6.986000
H   16.50000       23.29500       6.652000
H   14.74390       23.51764       6.710063
H   15.61789       22.72907       8.078654
H   17.63853       19.12101       5.827124
H   16.55212       19.14519       4.519021
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '12')] = qcdb.Molecule("""
0 1
C   5.250000       26.25800       10.94600
H   5.928000       26.22800       10.10000
C   5.965000       25.62100       12.14400
H   5.388000       25.81100       13.04900
H   6.001000       24.53800       11.98600
C   7.381000       26.17300       12.29600
H   7.851000       25.68400       13.15200
H   7.950000       25.96300       11.38000
H   7.328000       27.25600       12.46000
H   4.327890       25.69684       10.73431
H   4.991234       27.31620       11.09851
--
0 1
C   10.60200       24.81800       6.466000
O   10.95600       23.84000       7.103000
N   10.17800       25.94300       7.070000
C   10.09100       26.25600       8.476000
C   9.372000       27.59000       8.640000
C   11.44600       26.35600       9.091000
C   9.333000       25.25000       9.282000
H   9.874000       26.68900       6.497000
H   9.908000       28.37100       8.093000
H   8.364000       27.46400       8.233000
H   9.317000       27.84600       9.706000
H   9.807000       24.28200       9.160000
H   9.371000       25.57400       10.32900
H   8.328000       25.26700       8.900000
H   11.28800       26.57600       10.14400
H   11.97000       27.14900       8.585000
H   11.93200       25.39300       8.957000
H   10.61998       24.85900       5.366911
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '13')] = qcdb.Molecule("""
0 1
C   16.05600       13.91300       3.701000
H   16.15900       14.95600       4.033000
C   17.23700       13.52900       2.786000
H   17.24200       14.19600       1.903000
H   18.17200       13.64500       3.338000
H   17.11900       12.49000       2.453000
C   14.73500       13.74900       2.932000
H   14.59700       12.70200       2.661000
H   13.89800       14.05900       3.553000
C   14.73600       14.57900       1.670000
H   13.78900       14.41800       1.168000
H   14.86200       15.62400       1.955000
H   15.57300       14.24500       1.050000
H   16.03595       13.25839       4.584791
--
0 1
C   16.51300       18.99700       4.107000
C   17.41000       17.97200       3.810000
C   17.48600       17.37200       2.527000
C   16.60600       17.86100       1.547000
C   15.70300       18.89200       1.854000
C   15.65400       19.47000       3.123000
H   14.96900       20.29200       3.327000
H   15.03200       19.27500       1.092000
H   16.61400       17.45600       0.523000
H   18.19800       16.55200       2.242000
H   18.01270       17.77884       4.709692
H   16.61088       19.36781       5.137980
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '14')] = qcdb.Molecule("""
0 1
C   10.37400       21.44300       10.31100
H   9.934000       21.48700       9.305000
C   9.280000       21.86400       11.28300
H   8.982000       22.89200       11.04700
H   8.438000       21.18400       11.17100
H   9.688000       21.81300       12.29800
C   11.56200       22.42400       10.37400
H   12.42800       21.94400       9.922000
H   11.32000       23.31700       9.801000
C   11.94400       22.85500       11.78700
H   12.79300       23.53700       11.72100
H   11.08500       23.36200       12.24600
H   12.21000       21.96700       12.37000
H   10.70966       20.41678       10.52123
--
0 1
C   10.60200       24.81800       6.466000
O   10.95600       23.84000       7.103000
N   10.17800       25.94300       7.070000
C   10.09100       26.25600       8.476000
C   9.372000       27.59000       8.640000
C   11.44600       26.35600       9.091000
C   9.333000       25.25000       9.282000
H   9.874000       26.68900       6.497000
H   9.908000       28.37100       8.093000
H   8.364000       27.46400       8.233000
H   9.317000       27.84600       9.706000
H   9.807000       24.28200       9.160000
H   9.371000       25.57400       10.32900
H   8.328000       25.26700       8.900000
H   11.28800       26.57600       10.14400
H   11.97000       27.14900       8.585000
H   11.93200       25.39300       8.957000
H   10.61998       24.85900       5.366911
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '15')] = qcdb.Molecule("""
0 1
C   12.74600       22.16800      -1.090000
H   12.81800       23.17800      -0.699000
H   11.91700       22.12800      -1.801000
C   12.43800       21.21300       0.067000
H   11.49600       21.52600       0.536000
H   12.33700       20.19000      -0.320000
H   13.25200       21.25600       0.799000
H   13.67406       21.92095      -1.626354
--
0 1
N   9.254000       23.85400       3.012000
C   11.78200       23.78100       2.948000
C   10.44700       23.17200       2.478000
N   11.86200       24.24500       4.359000
H   10.43300       22.11300       2.752000
H   10.40400       23.24700       1.396000
H   12.53400       23.01100       2.841000
H   12.00900       24.60600       2.286000
H   8.356916       23.29360       2.710019
H   9.400716       23.94588       4.098293
H   10.95781       24.70625       4.782907
H   12.80321       24.79031       4.522592
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '16')] = qcdb.Molecule("""
0 1
C   20.37500       27.13100       10.18700
H   20.77800       26.39300       10.88300
H   19.30700       27.23700       10.38300
C   20.54000       26.57700       8.773000
H   20.37800       27.33600       8.009000
C   21.95300       26.03000       8.616000
H   22.05800       25.62600       7.609000
H   22.66100       26.85500       8.771000
H   22.12100       25.24800       9.363000
C   19.49000       25.48000       8.621000
H   19.56600       25.05200       7.621000
H   19.65600       24.70700       9.381000
H   18.48800       25.92200       8.763000
H   20.85217       28.10683       10.36039
--
0 1
C   16.73700       21.75300       8.985000
C   18.06300       21.93100       8.570000
C   19.04400       21.01900       8.966000
C   18.68400       19.93600       9.775000
C   17.35900       19.76800       10.19300
C   16.38200       20.67900       9.796000
H   15.33000       20.56400       10.09500
H   17.07400       18.92500       10.82100
H   19.43700       19.21300       10.07200
H   20.08100       21.14800       8.627000
H   18.32800       22.76900       7.913000
C   15.63700       22.68100       8.524000
H   15.79700       23.65100       8.994000
H   14.68200       22.28100       8.905000
H   15.59011       22.79893       7.431346
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '17')] = qcdb.Molecule("""
0 1
C   5.679000       20.88000       0.749000
H   6.221000       21.83800       0.688000
H   4.631000       21.09600       0.905000
C   6.254000       19.97400       1.854000
H   6.621000       20.53900       2.695000
H   5.514000       19.25100       2.203000
H   7.080433       19.44257       1.359442
H   5.830258       20.31081      -0.1800558
--
0 1
C   6.715000       23.76400       2.930000
C   6.218000       24.82200       2.171000
N   4.990000       25.38300       2.391000
C   4.226000       24.84700       3.389000
C   4.663000       23.78000       4.166000
C   5.918000       23.22600       3.934000
C   8.039000       23.09500       2.603000
H   8.026000       22.11200       3.087000
H   8.096000       22.92700       1.519000
H   6.274000       22.37100       4.504000
H   4.018000       23.38100       4.943000
H   3.211000       25.24700       3.575000
H   6.798000       25.19700       1.318000
H   8.936083       23.65540       2.904981
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '18')] = qcdb.Molecule("""
0 1
C   11.54100       27.68600       13.69600
H   12.45900       27.15000       13.44600
C   10.79000       27.96500       12.40600
H   10.55700       27.01400       11.92400
H   9.879000       28.51400       12.64300
H   11.44300       28.56800       11.76200
H   10.90337       27.06487       14.34224
H   11.78789       28.62476       14.21347
--
0 1
C   10.60200       24.81800       6.466000
O   10.95600       23.84000       7.103000
N   10.17800       25.94300       7.070000
C   10.09100       26.25600       8.476000
C   9.372000       27.59000       8.640000
C   11.44600       26.35600       9.091000
C   9.333000       25.25000       9.282000
H   9.874000       26.68900       6.497000
H   9.908000       28.37100       8.093000
H   8.364000       27.46400       8.233000
H   9.317000       27.84600       9.706000
H   9.807000       24.28200       9.160000
H   9.371000       25.57400       10.32900
H   8.328000       25.26700       8.900000
H   11.28800       26.57600       10.14400
H   11.97000       27.14900       8.585000
H   11.93200       25.39300       8.957000
H   10.61998       24.85900       5.366911
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '19')] = qcdb.Molecule("""
0 1
C   18.95600       23.00600       13.42400
H   19.58200       22.11000       13.49500
C   17.99700       22.76900       12.25600
H   18.56900       22.61800       11.32900
H   17.41200       21.87100       12.48100
H   17.34100       23.63600       12.15300
C   19.86100       24.16500       13.07400
H   20.34500       23.96700       12.11700
H   19.24600       25.06600       13.01700
H   20.59800       24.25000       13.87300
H   18.44867       23.19359       14.38182
--
0 1
C   16.73700       21.75300       8.985000
C   18.06300       21.93100       8.570000
C   19.04400       21.01900       8.966000
C   18.68400       19.93600       9.775000
C   17.35900       19.76800       10.19300
C   16.38200       20.67900       9.796000
H   15.33000       20.56400       10.09500
H   17.07400       18.92500       10.82100
H   19.43700       19.21300       10.07200
H   20.08100       21.14800       8.627000
H   18.32800       22.76900       7.913000
H   15.93631       22.42849       8.649437
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '20')] = qcdb.Molecule("""
0 1
C   13.79100       17.02500      -2.243000
H   12.90600       17.67000      -2.286000
C   13.29600       15.57900      -2.178000
H   12.68500       15.45200      -1.281000
H   12.69000       15.37300      -3.075000
H   14.15900       14.91000      -2.152000
C   14.52700       17.45100      -0.990000
H   13.87200       17.32300      -0.127000
H   15.43100       16.85400      -0.884000
H   14.78900       18.51600      -1.107000
H   14.41225       17.12759      -3.144954
--
0 1
C   16.51300       18.99700       4.107000
C   17.41000       17.97200       3.810000
C   17.48600       17.37200       2.527000
C   16.60600       17.86100       1.547000
C   15.70300       18.89200       1.854000
C   15.65400       19.47000       3.123000
H   14.96900       20.29200       3.327000
H   15.03200       19.27500       1.092000
H   16.61400       17.45600       0.523000
H   18.19800       16.55200       2.242000
H   18.01270       17.77884       4.709692
H   16.61088       19.36781       5.137980
units angstrom
""")

GEOS['%s-%s-dimer' % (dbse, '21')] = qcdb.Molecule("""
-1 1
O   8.976000       28.18400       5.336000
H   9.797000       28.02600       4.860000
H   8.600000       28.95300       4.860000
C   7.707000       31.02100       3.964000
O   8.101000       31.66300       2.963000
O   7.068000       29.95800       3.828000
C   8.014000       31.53300       5.387000
H   7.465000       32.46700       5.552000
H   7.660000       30.81200       6.145000
H   9.081843       31.72975       5.563073
--
0 1
C   10.62700       24.87500       4.938000
C   10.60200       24.81800       6.466000
O   10.95600       23.84000       7.103000
N   10.17800       25.94300       7.070000
C   10.09100       26.25600       8.476000
H   9.874000       26.68900       6.497000
H   10.49500       25.96800       4.694000
H   9.572130       27.21869       8.594352
H   9.533425       25.51600       9.068883
H   11.09040       26.32976       8.929603
H   11.53119       24.41376       4.514093
H   9.792144       24.23767       4.611163
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
DATA['NUCLEAR REPULSION ENERGY']['HSG-1-dimer'                    ] =     409.61526850
DATA['NUCLEAR REPULSION ENERGY']['HSG-1-monoA-unCP'               ] =      13.33595232
DATA['NUCLEAR REPULSION ENERGY']['HSG-1-monoB-unCP'               ] =     332.12261009
DATA['NUCLEAR REPULSION ENERGY']['HSG-2-dimer'                    ] =     693.84322132
DATA['NUCLEAR REPULSION ENERGY']['HSG-2-monoA-unCP'               ] =      41.89071165
DATA['NUCLEAR REPULSION ENERGY']['HSG-2-monoB-unCP'               ] =     501.75349414
DATA['NUCLEAR REPULSION ENERGY']['HSG-3-dimer'                    ] =     578.08454963
DATA['NUCLEAR REPULSION ENERGY']['HSG-3-monoA-unCP'               ] =     194.80446994
DATA['NUCLEAR REPULSION ENERGY']['HSG-3-monoB-unCP'               ] =     202.93507303
DATA['NUCLEAR REPULSION ENERGY']['HSG-4-dimer'                    ] =     536.02111700
DATA['NUCLEAR REPULSION ENERGY']['HSG-4-monoA-unCP'               ] =     336.06029689
DATA['NUCLEAR REPULSION ENERGY']['HSG-4-monoB-unCP'               ] =      40.09418196
DATA['NUCLEAR REPULSION ENERGY']['HSG-5-dimer'                    ] =     641.07583890
DATA['NUCLEAR REPULSION ENERGY']['HSG-5-monoA-unCP'               ] =     440.48402439
DATA['NUCLEAR REPULSION ENERGY']['HSG-5-monoB-unCP'               ] =      39.79355972
DATA['NUCLEAR REPULSION ENERGY']['HSG-6-dimer'                    ] =     440.32913479
DATA['NUCLEAR REPULSION ENERGY']['HSG-6-monoA-unCP'               ] =     112.25425669
DATA['NUCLEAR REPULSION ENERGY']['HSG-6-monoB-unCP'               ] =     202.38032057
DATA['NUCLEAR REPULSION ENERGY']['HSG-7-dimer'                    ] =     825.37483209
DATA['NUCLEAR REPULSION ENERGY']['HSG-7-monoA-unCP'               ] =     302.68630925
DATA['NUCLEAR REPULSION ENERGY']['HSG-7-monoB-unCP'               ] =     256.12378323
DATA['NUCLEAR REPULSION ENERGY']['HSG-8-dimer'                    ] =     721.36437027
DATA['NUCLEAR REPULSION ENERGY']['HSG-8-monoA-unCP'               ] =     298.54657988
DATA['NUCLEAR REPULSION ENERGY']['HSG-8-monoB-unCP'               ] =     204.68604075
DATA['NUCLEAR REPULSION ENERGY']['HSG-9-dimer'                    ] =     699.77856295
DATA['NUCLEAR REPULSION ENERGY']['HSG-9-monoA-unCP'               ] =     298.58992071
DATA['NUCLEAR REPULSION ENERGY']['HSG-9-monoB-unCP'               ] =     179.49546339
DATA['NUCLEAR REPULSION ENERGY']['HSG-10-dimer'                   ] =     538.20524151
DATA['NUCLEAR REPULSION ENERGY']['HSG-10-monoA-unCP'              ] =     179.66798724
DATA['NUCLEAR REPULSION ENERGY']['HSG-10-monoB-unCP'              ] =     180.34079666
DATA['NUCLEAR REPULSION ENERGY']['HSG-11-dimer'                   ] =     697.51311416
DATA['NUCLEAR REPULSION ENERGY']['HSG-11-monoA-unCP'              ] =     296.89990217
DATA['NUCLEAR REPULSION ENERGY']['HSG-11-monoB-unCP'              ] =     180.34079666
DATA['NUCLEAR REPULSION ENERGY']['HSG-12-dimer'                   ] =     553.87245309
DATA['NUCLEAR REPULSION ENERGY']['HSG-12-monoA-unCP'              ] =      82.71734142
DATA['NUCLEAR REPULSION ENERGY']['HSG-12-monoB-unCP'              ] =     332.12261009
DATA['NUCLEAR REPULSION ENERGY']['HSG-13-dimer'                   ] =     492.23285254
DATA['NUCLEAR REPULSION ENERGY']['HSG-13-monoA-unCP'              ] =     134.28280330
DATA['NUCLEAR REPULSION ENERGY']['HSG-13-monoB-unCP'              ] =     202.38032057
DATA['NUCLEAR REPULSION ENERGY']['HSG-14-dimer'                   ] =     670.02074299
DATA['NUCLEAR REPULSION ENERGY']['HSG-14-monoA-unCP'              ] =     134.10189365
DATA['NUCLEAR REPULSION ENERGY']['HSG-14-monoB-unCP'              ] =     332.12261009
DATA['NUCLEAR REPULSION ENERGY']['HSG-15-dimer'                   ] =     242.88545739
DATA['NUCLEAR REPULSION ENERGY']['HSG-15-monoA-unCP'              ] =      42.22202660
DATA['NUCLEAR REPULSION ENERGY']['HSG-15-monoB-unCP'              ] =     131.69625678
DATA['NUCLEAR REPULSION ENERGY']['HSG-16-dimer'                   ] =     551.59382982
DATA['NUCLEAR REPULSION ENERGY']['HSG-16-monoA-unCP'              ] =     135.70381177
DATA['NUCLEAR REPULSION ENERGY']['HSG-16-monoB-unCP'              ] =     269.04078448
DATA['NUCLEAR REPULSION ENERGY']['HSG-17-dimer'                   ] =     421.73710621
DATA['NUCLEAR REPULSION ENERGY']['HSG-17-monoA-unCP'              ] =      42.20972067
DATA['NUCLEAR REPULSION ENERGY']['HSG-17-monoB-unCP'              ] =     270.70970086
DATA['NUCLEAR REPULSION ENERGY']['HSG-18-dimer'                   ] =     474.74808030
DATA['NUCLEAR REPULSION ENERGY']['HSG-18-monoA-unCP'              ] =      42.43370398
DATA['NUCLEAR REPULSION ENERGY']['HSG-18-monoB-unCP'              ] =     332.12261009
DATA['NUCLEAR REPULSION ENERGY']['HSG-19-dimer'                   ] =     410.08888873
DATA['NUCLEAR REPULSION ENERGY']['HSG-19-monoA-unCP'              ] =      83.35857717
DATA['NUCLEAR REPULSION ENERGY']['HSG-19-monoB-unCP'              ] =     202.93507303
DATA['NUCLEAR REPULSION ENERGY']['HSG-20-dimer'                   ] =     392.20505391
DATA['NUCLEAR REPULSION ENERGY']['HSG-20-monoA-unCP'              ] =      82.90559609
DATA['NUCLEAR REPULSION ENERGY']['HSG-20-monoB-unCP'              ] =     202.38032057
DATA['NUCLEAR REPULSION ENERGY']['HSG-21-dimer'                   ] =     495.71409832
DATA['NUCLEAR REPULSION ENERGY']['HSG-21-monoA-unCP'              ] =     169.11593456
DATA['NUCLEAR REPULSION ENERGY']['HSG-21-monoB-unCP'              ] =     179.49546339
DATA['NUCLEAR REPULSION ENERGY']['HSG-1-monoA-CP'                 ] =      13.33595232
DATA['NUCLEAR REPULSION ENERGY']['HSG-1-monoB-CP'                 ] =     332.12261009
DATA['NUCLEAR REPULSION ENERGY']['HSG-2-monoA-CP'                 ] =      41.89071165
DATA['NUCLEAR REPULSION ENERGY']['HSG-2-monoB-CP'                 ] =     501.75349414
DATA['NUCLEAR REPULSION ENERGY']['HSG-3-monoA-CP'                 ] =     194.80446994
DATA['NUCLEAR REPULSION ENERGY']['HSG-3-monoB-CP'                 ] =     202.93507303
DATA['NUCLEAR REPULSION ENERGY']['HSG-4-monoA-CP'                 ] =     336.06029689
DATA['NUCLEAR REPULSION ENERGY']['HSG-4-monoB-CP'                 ] =      40.09418196
DATA['NUCLEAR REPULSION ENERGY']['HSG-5-monoA-CP'                 ] =     440.48402439
DATA['NUCLEAR REPULSION ENERGY']['HSG-5-monoB-CP'                 ] =      39.79355972
DATA['NUCLEAR REPULSION ENERGY']['HSG-6-monoA-CP'                 ] =     112.25425669
DATA['NUCLEAR REPULSION ENERGY']['HSG-6-monoB-CP'                 ] =     202.38032057
DATA['NUCLEAR REPULSION ENERGY']['HSG-7-monoA-CP'                 ] =     302.68630925
DATA['NUCLEAR REPULSION ENERGY']['HSG-7-monoB-CP'                 ] =     256.12378323
DATA['NUCLEAR REPULSION ENERGY']['HSG-8-monoA-CP'                 ] =     298.54657988
DATA['NUCLEAR REPULSION ENERGY']['HSG-8-monoB-CP'                 ] =     204.68604075
DATA['NUCLEAR REPULSION ENERGY']['HSG-9-monoA-CP'                 ] =     298.58992071
DATA['NUCLEAR REPULSION ENERGY']['HSG-9-monoB-CP'                 ] =     179.49546339
DATA['NUCLEAR REPULSION ENERGY']['HSG-10-monoA-CP'                ] =     179.66798724
DATA['NUCLEAR REPULSION ENERGY']['HSG-10-monoB-CP'                ] =     180.34079666
DATA['NUCLEAR REPULSION ENERGY']['HSG-11-monoA-CP'                ] =     296.89990217
DATA['NUCLEAR REPULSION ENERGY']['HSG-11-monoB-CP'                ] =     180.34079666
DATA['NUCLEAR REPULSION ENERGY']['HSG-12-monoA-CP'                ] =      82.71734142
DATA['NUCLEAR REPULSION ENERGY']['HSG-12-monoB-CP'                ] =     332.12261009
DATA['NUCLEAR REPULSION ENERGY']['HSG-13-monoA-CP'                ] =     134.28280330
DATA['NUCLEAR REPULSION ENERGY']['HSG-13-monoB-CP'                ] =     202.38032057
DATA['NUCLEAR REPULSION ENERGY']['HSG-14-monoA-CP'                ] =     134.10189365
DATA['NUCLEAR REPULSION ENERGY']['HSG-14-monoB-CP'                ] =     332.12261009
DATA['NUCLEAR REPULSION ENERGY']['HSG-15-monoA-CP'                ] =      42.22202660
DATA['NUCLEAR REPULSION ENERGY']['HSG-15-monoB-CP'                ] =     131.69625678
DATA['NUCLEAR REPULSION ENERGY']['HSG-16-monoA-CP'                ] =     135.70381177
DATA['NUCLEAR REPULSION ENERGY']['HSG-16-monoB-CP'                ] =     269.04078448
DATA['NUCLEAR REPULSION ENERGY']['HSG-17-monoA-CP'                ] =      42.20972067
DATA['NUCLEAR REPULSION ENERGY']['HSG-17-monoB-CP'                ] =     270.70970086
DATA['NUCLEAR REPULSION ENERGY']['HSG-18-monoA-CP'                ] =      42.43370398
DATA['NUCLEAR REPULSION ENERGY']['HSG-18-monoB-CP'                ] =     332.12261009
DATA['NUCLEAR REPULSION ENERGY']['HSG-19-monoA-CP'                ] =      83.35857717
DATA['NUCLEAR REPULSION ENERGY']['HSG-19-monoB-CP'                ] =     202.93507303
DATA['NUCLEAR REPULSION ENERGY']['HSG-20-monoA-CP'                ] =      82.90559609
DATA['NUCLEAR REPULSION ENERGY']['HSG-20-monoB-CP'                ] =     202.38032057
DATA['NUCLEAR REPULSION ENERGY']['HSG-21-monoA-CP'                ] =     169.11593456
DATA['NUCLEAR REPULSION ENERGY']['HSG-21-monoB-CP'                ] =     179.49546339
