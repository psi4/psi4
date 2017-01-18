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
| Database of Hydrogen transfer reactions.
| Geometries from Bozkaya and Sherrill.
| Reference energies from Bozkaya and Sherrill.


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

# <<< HTR40 Database Module >>>
dbse = 'HTR40'
isOS = 'True'

# <<< Database Members >>>
HRXN = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
HRXN_SM = []
HRXN_LG = []
HTR15 = [1, 2, 3, 4, 11, 12, 13, 20, 21, 34, 35, 36, 37, 38,  39]

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV['%s-%s'            % (dbse, 1                     )] = ['%s-%s-reagent'      % (dbse, 'ch3'),
                                                             '%s-%s-reagent'      % (dbse, 'h2'),
                                                             '%s-%s-reagent'      % (dbse, 'ch4'),
                                                             '%s-%s-reagent'      % (dbse, 'h') ]
RXNM['%s-%s'            % (dbse, 1                     )] = dict(zip(ACTV['%s-%s' % (dbse, 1)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 2                     )] = ['%s-%s-reagent'      % (dbse, 'c2h'),
                                                              '%s-%s-reagent'      % (dbse, 'h2'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h2'),
                                                              '%s-%s-reagent'      % (dbse, 'h') ]
RXNM['%s-%s'            % (dbse, 2                     )] = dict(zip(ACTV['%s-%s' % (dbse, 2)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 3                     )] = ['%s-%s-reagent'      % (dbse, 'c2h3'),
                                                              '%s-%s-reagent'      % (dbse, 'h2'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h4'),
                                                              '%s-%s-reagent'      % (dbse, 'h') ]
RXNM['%s-%s'            % (dbse, 3                     )] = dict(zip(ACTV['%s-%s' % (dbse, 3)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 4                     )] = ['%s-%s-reagent'      % (dbse, 't-butyl'),
                                                              '%s-%s-reagent'      % (dbse, 'h2'),
                                                              '%s-%s-reagent'      % (dbse, 'isobutane'),
                                                              '%s-%s-reagent'      % (dbse, 'h') ]
RXNM['%s-%s'            % (dbse, 4                     )] = dict(zip(ACTV['%s-%s' % (dbse, 4)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 5                     )] = ['%s-%s-reagent'      % (dbse, 'cfch2'),
                                                              '%s-%s-reagent'      % (dbse, 'h2'),
                                                              '%s-%s-reagent'      % (dbse, 'chfch2'),
                                                              '%s-%s-reagent'      % (dbse, 'h') ]
RXNM['%s-%s'            % (dbse, 5                     )] = dict(zip(ACTV['%s-%s' % (dbse, 5)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 6                     )] = ['%s-%s-reagent'      % (dbse, 'ch2cho'),
                                                              '%s-%s-reagent'      % (dbse, 'h2'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3cho'),
                                                              '%s-%s-reagent'      % (dbse, 'h') ]
RXNM['%s-%s'            % (dbse, 6                     )] = dict(zip(ACTV['%s-%s' % (dbse, 6)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 7                     )] = ['%s-%s-reagent'      % (dbse, 'ch2cn'),
                                                              '%s-%s-reagent'      % (dbse, 'h2'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3cn'),
                                                              '%s-%s-reagent'      % (dbse, 'h') ]
RXNM['%s-%s'            % (dbse, 7                     )] = dict(zip(ACTV['%s-%s' % (dbse, 7)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 8                     )] = ['%s-%s-reagent'      % (dbse, 'ch2cch'),
                                                              '%s-%s-reagent'      % (dbse, 'h2'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3cch'),
                                                              '%s-%s-reagent'      % (dbse, 'h') ]
RXNM['%s-%s'            % (dbse, 8                     )] = dict(zip(ACTV['%s-%s' % (dbse, 8)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 9                     )] = ['%s-%s-reagent'      % (dbse, 'ch2ccn'),
                                                              '%s-%s-reagent'      % (dbse, 'h2'),
                                                              '%s-%s-reagent'      % (dbse, 'ch2chcn'),
                                                              '%s-%s-reagent'      % (dbse, 'h') ]
RXNM['%s-%s'            % (dbse, 9                     )] = dict(zip(ACTV['%s-%s' % (dbse, 9)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 10                     )] = ['%s-%s-reagent'      % (dbse, 'allyl'),
                                                              '%s-%s-reagent'      % (dbse, 'h2'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3chch2'),
                                                              '%s-%s-reagent'      % (dbse, 'h') ]
RXNM['%s-%s'            % (dbse, 10                     )] = dict(zip(ACTV['%s-%s' % (dbse, 10)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 11                     )] = ['%s-%s-reagent'      % (dbse, 'c2h'),
                                                              '%s-%s-reagent'      % (dbse, 'ch4'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h2'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3') ]
RXNM['%s-%s'            % (dbse, 11                     )] = dict(zip(ACTV['%s-%s' % (dbse, 11)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 12                     )] = ['%s-%s-reagent'      % (dbse, 'c2h3'),
                                                              '%s-%s-reagent'      % (dbse, 'ch4'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h4'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3') ]
RXNM['%s-%s'            % (dbse, 12                     )] = dict(zip(ACTV['%s-%s' % (dbse, 12)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 13                     )] = ['%s-%s-reagent'      % (dbse, 't-butyl'),
                                                              '%s-%s-reagent'      % (dbse, 'ch4'),
                                                              '%s-%s-reagent'      % (dbse, 'isobutane'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3') ]
RXNM['%s-%s'            % (dbse, 13                     )] = dict(zip(ACTV['%s-%s' % (dbse, 13)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 14                     )] = ['%s-%s-reagent'      % (dbse, 'cfch2'),
                                                              '%s-%s-reagent'      % (dbse, 'ch4'),
                                                              '%s-%s-reagent'      % (dbse, 'chfch2'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3') ]
RXNM['%s-%s'            % (dbse, 14                     )] = dict(zip(ACTV['%s-%s' % (dbse, 14)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 15                     )] = ['%s-%s-reagent'      % (dbse, 'ch2cho'),
                                                              '%s-%s-reagent'      % (dbse, 'ch4'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3cho'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3') ]
RXNM['%s-%s'            % (dbse, 15                     )] = dict(zip(ACTV['%s-%s' % (dbse, 15)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 16                     )] = ['%s-%s-reagent'      % (dbse, 'ch2cn'),
                                                              '%s-%s-reagent'      % (dbse, 'ch4'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3cn'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3') ]
RXNM['%s-%s'            % (dbse, 16                     )] = dict(zip(ACTV['%s-%s' % (dbse, 16)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 17                     )] = ['%s-%s-reagent'      % (dbse, 'ch2cch'),
                                                              '%s-%s-reagent'      % (dbse, 'ch4'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3cch'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3') ]
RXNM['%s-%s'            % (dbse, 17                     )] = dict(zip(ACTV['%s-%s' % (dbse, 17)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 18                     )] = ['%s-%s-reagent'      % (dbse, 'ch2ccn'),
                                                              '%s-%s-reagent'      % (dbse, 'ch4'),
                                                              '%s-%s-reagent'      % (dbse, 'ch2chcn'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3') ]
RXNM['%s-%s'            % (dbse, 18                     )] = dict(zip(ACTV['%s-%s' % (dbse, 18)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 19                     )] = ['%s-%s-reagent'      % (dbse, 'allyl'),
                                                              '%s-%s-reagent'      % (dbse, 'ch4'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3chch2'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3') ]
RXNM['%s-%s'            % (dbse, 19                     )] = dict(zip(ACTV['%s-%s' % (dbse, 19)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 20                     )] = ['%s-%s-reagent'      % (dbse, 'c2h'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h4'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h2'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h3') ]
RXNM['%s-%s'            % (dbse, 20                     )] = dict(zip(ACTV['%s-%s' % (dbse, 20)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 21                     )] = ['%s-%s-reagent'      % (dbse, 't-butyl'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h4'),
                                                              '%s-%s-reagent'      % (dbse, 'isobutane'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h3') ]
RXNM['%s-%s'            % (dbse, 21                     )] = dict(zip(ACTV['%s-%s' % (dbse, 21)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 22                     )] = ['%s-%s-reagent'      % (dbse, 'cfch2'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h4'),
                                                              '%s-%s-reagent'      % (dbse, 'chfch2'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h3') ]
RXNM['%s-%s'            % (dbse, 22                     )] = dict(zip(ACTV['%s-%s' % (dbse, 22)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 23                     )] = ['%s-%s-reagent'      % (dbse, 'ch2cho'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h4'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3cho'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h3') ]
RXNM['%s-%s'            % (dbse, 23                     )] = dict(zip(ACTV['%s-%s' % (dbse, 23)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 24                     )] = ['%s-%s-reagent'      % (dbse, 'ch2cn'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h4'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3cn'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h3') ]
RXNM['%s-%s'            % (dbse, 24                     )] = dict(zip(ACTV['%s-%s' % (dbse, 24)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 25                     )] = ['%s-%s-reagent'      % (dbse, 'ch2cch'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h4'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3cch'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h3') ]
RXNM['%s-%s'            % (dbse, 25                     )] = dict(zip(ACTV['%s-%s' % (dbse, 25)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 26                     )] = ['%s-%s-reagent'      % (dbse, 'ch2ccn'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h4'),
                                                              '%s-%s-reagent'      % (dbse, 'ch2chcn'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h3') ]
RXNM['%s-%s'            % (dbse, 26                     )] = dict(zip(ACTV['%s-%s' % (dbse, 26)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 27                     )] = ['%s-%s-reagent'      % (dbse, 'allyl'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h4'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3chch2'),
                                                              '%s-%s-reagent'      % (dbse, 'h') ]
RXNM['%s-%s'            % (dbse, 27                     )] = dict(zip(ACTV['%s-%s' % (dbse, 27)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 28                     )] = ['%s-%s-reagent'      % (dbse, 'cfch2'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h2'),
                                                              '%s-%s-reagent'      % (dbse, 'chfch2'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h') ]
RXNM['%s-%s'            % (dbse, 28                     )] = dict(zip(ACTV['%s-%s' % (dbse, 28)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 29                     )] = ['%s-%s-reagent'      % (dbse, 'ch2cho'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h2'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3cho'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h') ]
RXNM['%s-%s'            % (dbse, 29                     )] = dict(zip(ACTV['%s-%s' % (dbse, 29)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 30                     )] = ['%s-%s-reagent'      % (dbse, 'ch2cn'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h2'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3cn'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h') ]
RXNM['%s-%s'            % (dbse, 30                     )] = dict(zip(ACTV['%s-%s' % (dbse, 30)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 31                     )] = ['%s-%s-reagent'      % (dbse, 'ch2cch'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h2'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3cch'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h') ]
RXNM['%s-%s'            % (dbse, 31                     )] = dict(zip(ACTV['%s-%s' % (dbse, 31)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 32                     )] = ['%s-%s-reagent'      % (dbse, 'ch2ccn'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h2'),
                                                              '%s-%s-reagent'      % (dbse, 'ch2chcn'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h') ]
RXNM['%s-%s'            % (dbse, 32                     )] = dict(zip(ACTV['%s-%s' % (dbse, 32)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 33                     )] = ['%s-%s-reagent'      % (dbse, 'allyl'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h2'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3chch2'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h') ]
RXNM['%s-%s'            % (dbse, 33                     )] = dict(zip(ACTV['%s-%s' % (dbse, 33)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 34                     )] = ['%s-%s-reagent'      % (dbse, 'c2h'),
                                                              '%s-%s-reagent'      % (dbse, 'isobutane'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h2'),
                                                              '%s-%s-reagent'      % (dbse, 't-butyl') ]
RXNM['%s-%s'            % (dbse, 34                     )] = dict(zip(ACTV['%s-%s' % (dbse, 34)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 35                     )] = ['%s-%s-reagent'      % (dbse, 'c6h5'),
                                                              '%s-%s-reagent'      % (dbse, 'h2'),
                                                              '%s-%s-reagent'      % (dbse, 'c6h6'),
                                                              '%s-%s-reagent'      % (dbse, 'h') ]
RXNM['%s-%s'            % (dbse, 35                     )] = dict(zip(ACTV['%s-%s' % (dbse, 35)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 36                     )] = ['%s-%s-reagent'      % (dbse, 'c6h5'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h4'),
                                                              '%s-%s-reagent'      % (dbse, 'c6h6'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h3') ]
RXNM['%s-%s'            % (dbse, 36                     )] = dict(zip(ACTV['%s-%s' % (dbse, 36)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 37                     )] = ['%s-%s-reagent'      % (dbse, 'c6h5'),
                                                              '%s-%s-reagent'      % (dbse, 'isobutane'),
                                                              '%s-%s-reagent'      % (dbse, 'c6h6'),
                                                              '%s-%s-reagent'      % (dbse, 't-butyl') ]
RXNM['%s-%s'            % (dbse, 37                     )] = dict(zip(ACTV['%s-%s' % (dbse, 37)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 38                     )] = ['%s-%s-reagent'      % (dbse, 'c2h'),
                                                              '%s-%s-reagent'      % (dbse, 'c6h6'),
                                                              '%s-%s-reagent'      % (dbse, 'c2h2'),
                                                              '%s-%s-reagent'      % (dbse, 'c6h5') ]
RXNM['%s-%s'            % (dbse, 38                     )] = dict(zip(ACTV['%s-%s' % (dbse, 38)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 39                     )] = ['%s-%s-reagent'      % (dbse, 'c6h5'),
                                                              '%s-%s-reagent'      % (dbse, 'ch4'),
                                                              '%s-%s-reagent'      % (dbse, 'c6h6'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3') ]
RXNM['%s-%s'            % (dbse, 39                     )] = dict(zip(ACTV['%s-%s' % (dbse, 39)], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, 40                     )] = ['%s-%s-reagent'      % (dbse, 'c6h5'),
                                                              '%s-%s-reagent'      % (dbse, 'ch3cn'),
                                                              '%s-%s-reagent'      % (dbse, 'c6h6'),
                                                              '%s-%s-reagent'      % (dbse, 'ch2cn') ]
RXNM['%s-%s'            % (dbse, 40                     )] = dict(zip(ACTV['%s-%s' % (dbse, 40)], [-1,-1,1,1]))

# <<< Reference Values [kcal/mol] >>>
BIND = {}
BIND['%s-%s'            % (dbse, 1                    )] =    0.000
BIND['%s-%s'            % (dbse, 2                    )] =    0.000
BIND['%s-%s'            % (dbse, 3                    )] =    0.000
BIND['%s-%s'            % (dbse, 4                    )] =    0.000
BIND['%s-%s'            % (dbse, 5                    )] =    0.000
BIND['%s-%s'            % (dbse, 6                    )] =    0.000
BIND['%s-%s'            % (dbse, 7                    )] =    0.000
BIND['%s-%s'            % (dbse, 8                    )] =    0.000
BIND['%s-%s'            % (dbse, 9                    )] =    0.000
BIND['%s-%s'            % (dbse, 10                    )] =    0.000
BIND['%s-%s'            % (dbse, 11                    )] =    0.000
BIND['%s-%s'            % (dbse, 12                    )] =    0.000
BIND['%s-%s'            % (dbse, 13                    )] =    0.000
BIND['%s-%s'            % (dbse, 14                    )] =    0.000
BIND['%s-%s'            % (dbse, 15                    )] =    0.000
BIND['%s-%s'            % (dbse, 16                    )] =    0.000
BIND['%s-%s'            % (dbse, 17                    )] =    0.000
BIND['%s-%s'            % (dbse, 18                    )] =    0.000
BIND['%s-%s'            % (dbse, 19                    )] =    0.000
BIND['%s-%s'            % (dbse, 20                    )] =    0.000
BIND['%s-%s'            % (dbse, 21                    )] =    0.000
BIND['%s-%s'            % (dbse, 22                    )] =    0.000
BIND['%s-%s'            % (dbse, 23                    )] =    0.000
BIND['%s-%s'            % (dbse, 24                    )] =    0.000
BIND['%s-%s'            % (dbse, 25                    )] =    0.000
BIND['%s-%s'            % (dbse, 27                    )] =    0.000
BIND['%s-%s'            % (dbse, 28                    )] =    0.000
BIND['%s-%s'            % (dbse, 29                    )] =    0.000
BIND['%s-%s'            % (dbse, 30                    )] =    0.000
BIND['%s-%s'            % (dbse, 31                    )] =    0.000
BIND['%s-%s'            % (dbse, 32                    )] =    0.000
BIND['%s-%s'            % (dbse, 33                    )] =    0.000
BIND['%s-%s'            % (dbse, 34                    )] =    0.000
BIND['%s-%s'            % (dbse, 35                    )] =    0.000
BIND['%s-%s'            % (dbse, 36                    )] =    0.000
BIND['%s-%s'            % (dbse, 37                    )] =    0.000
BIND['%s-%s'            % (dbse, 38                    )] =    0.000
BIND['%s-%s'            % (dbse, 39                    )] =    0.000
BIND['%s-%s'            % (dbse, 40                    )] =    0.000


# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, 1                     )] = """Reaction 1 """
TAGL['%s-%s'            % (dbse, 2                     )] = """Reaction 2 """
TAGL['%s-%s'            % (dbse, 3                     )] = """Reaction 3 """
TAGL['%s-%s'            % (dbse, 4                     )] = """Reaction 4 """
TAGL['%s-%s'            % (dbse, 5                     )] = """Reaction 5 """
TAGL['%s-%s'            % (dbse, 6                     )] = """Reaction 6 """
TAGL['%s-%s'            % (dbse, 7                     )] = """Reaction 7 """
TAGL['%s-%s'            % (dbse, 8                     )] = """Reaction 8 """
TAGL['%s-%s'            % (dbse, 9                     )] = """Reaction 9 """
TAGL['%s-%s'            % (dbse, 10                    )] = """Reaction 10 """
TAGL['%s-%s'            % (dbse, 11                    )] = """Reaction 11 """
TAGL['%s-%s'            % (dbse, 12                    )] = """Reaction 12 """
TAGL['%s-%s'            % (dbse, 13                    )] = """Reaction 13 """
TAGL['%s-%s'            % (dbse, 14                    )] = """Reaction 14 """
TAGL['%s-%s'            % (dbse, 15                    )] = """Reaction 15 """
TAGL['%s-%s'            % (dbse, 16                    )] = """Reaction 16 """
TAGL['%s-%s'            % (dbse, 17                    )] = """Reaction 17 """
TAGL['%s-%s'            % (dbse, 18                    )] = """Reaction 18 """
TAGL['%s-%s'            % (dbse, 19                    )] = """Reaction 19 """
TAGL['%s-%s'            % (dbse, 20                    )] = """Reaction 20 """
TAGL['%s-%s'            % (dbse, 21                    )] = """Reaction 21 """
TAGL['%s-%s'            % (dbse, 22                    )] = """Reaction 22 """
TAGL['%s-%s'            % (dbse, 23                    )] = """Reaction 23 """
TAGL['%s-%s'            % (dbse, 24                    )] = """Reaction 24 """
TAGL['%s-%s'            % (dbse, 25                    )] = """Reaction 25 """
TAGL['%s-%s'            % (dbse, 26                    )] = """Reaction 26 """
TAGL['%s-%s'            % (dbse, 27                    )] = """Reaction 27 """
TAGL['%s-%s'            % (dbse, 28                    )] = """Reaction 28 """
TAGL['%s-%s'            % (dbse, 29                    )] = """Reaction 29 """
TAGL['%s-%s'            % (dbse, 30                    )] = """Reaction 30 """
TAGL['%s-%s'            % (dbse, 31                    )] = """Reaction 31 """
TAGL['%s-%s'            % (dbse, 32                    )] = """Reaction 32 """
TAGL['%s-%s'            % (dbse, 33                    )] = """Reaction 33 """
TAGL['%s-%s'            % (dbse, 34                    )] = """Reaction 34 """
TAGL['%s-%s'            % (dbse, 35                    )] = """Reaction 35 """
TAGL['%s-%s'            % (dbse, 36                    )] = """Reaction 36 """
TAGL['%s-%s'            % (dbse, 37                    )] = """Reaction 37 """
TAGL['%s-%s'            % (dbse, 38                    )] = """Reaction 38 """
TAGL['%s-%s'            % (dbse, 39                    )] = """Reaction 39 """
TAGL['%s-%s'            % (dbse, 40                    )] = """Reaction 40 """
TAGL['%s-%s-reagent'            % (dbse, 't-butyl'               )] = """ """
TAGL['%s-%s-reagent'            % (dbse, 'cfch2'                 )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'c2h2'                  )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'ch3cho'                )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'ch2cn'                 )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'c2h4'                  )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'ch2chcn'               )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'c2h'                   )] = """ """  
TAGL['%s-%s-reagent'            % (dbse, 'ch4'                   )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'c2h3'                  )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'ch3cn'                 )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'allyl'                 )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'ch3cch'                )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'ch3'                   )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'h2'                    )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'ch3chch2'              )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'chfch2'                )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'h'                     )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'isobutane'             )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'ch2cho'                )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'ch2cch'                )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'ch2ccn'                )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'c6h6'                  )] = """ """   
TAGL['%s-%s-reagent'            % (dbse, 'c6h5'                  )] = """ """   

# <<< Geometry Specification Strings >>>
GEOS = {}

GEOS['%s-%s-%s' % (dbse, 't-butyl', 'reagent')] = qcdb.Molecule("""
0 2
C          0.00000073      -0.17324401       0.00000000
C         -1.49313671       0.03293053       0.00000000
C          0.74656816       0.03293930       1.29309432
C          0.74656816       0.03293930      -1.29309432
H         -1.75412591       1.11750215       0.00000000
H          0.87705536       1.11751245       1.51911455
H          0.87705536       1.11751245      -1.51911455
H         -1.96447858      -0.41104829       0.89724346
H         -1.96447858      -0.41104829      -0.89724346
H          1.75927771      -0.41103277       1.25266784
H          0.20520649      -0.41104076       2.14991136
H          0.20520649      -0.41104076      -2.14991136
H          1.75927771      -0.41103277      -1.25266784
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'cfch2', 'reagent')] = qcdb.Molecule("""
0 2
F          1.18236000       0.12888000       0.00008000
C         -0.00743000      -0.42937000      -0.00017000
C         -1.20008000       0.11938000      -0.00004000
H         -2.08423000      -0.50156000       0.00065000
H         -1.31194000       1.20159000      -0.00015000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'c2h2', 'reagent')] = qcdb.Molecule("""
0 1
C         -0.00000000       0.00000000       0.60249005
C          0.00000000       0.00000000      -0.60249005
H          0.00000000      -0.00000000       1.66141025
H          0.00000000      -0.00000000      -1.66141025
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3cho', 'reagent')] = qcdb.Molecule("""
0 1
H          1.15548000      -1.23749000      -0.00007000
C          1.16885000      -0.14776000      -0.00001000
C         -0.23563000       0.39721000      -0.00001000
H         -0.30509000       1.50872000       0.00003000
H          1.70764000       0.22227000       0.87912000
H          1.70778000       0.22245000      -0.87897000
O         -1.23314000      -0.27658000       0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2cn', 'reagent')] = qcdb.Molecule("""
0 2
C          0.18728000       0.00000000      -0.00001000
N          1.35587000       0.00000000       0.00000000
C         -1.19103000       0.00000000       0.00000000
H         -1.73431000       0.93522000       0.00001000
H         -1.73432000      -0.93521000       0.00001000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'c2h4', 'reagent')] = qcdb.Molecule("""
0 1
C          0.00000000       0.00000000       0.66741206
C          0.00000000       0.00000000      -0.66741206
H         -0.00000000       0.92046521       1.22998610
H          0.00000000      -0.92046521       1.22998610
H         -0.00000000       0.92046521      -1.22998610
H          0.00000000      -0.92046521      -1.22998610
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2chcn', 'reagent')] = qcdb.Molecule("""
0 1
H         -2.62698000       0.00343000       0.00128000
C         -1.60538000      -0.35565000       0.00019000
C         -0.58353000       0.50215000      -0.00012000
C          0.78296000       0.08953000      -0.00091000
N          1.89484000      -0.22376000       0.00060000
H         -1.45042000      -1.42771000      -0.00098000
H         -0.75076000       1.57438000       0.00050000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'c2h', 'reagent')] = qcdb.Molecule("""
0 2
C          0.00000000       0.00000000       0.00000000
C          1.21283562       0.00000000       0.00000000
H         -1.05818189       0.00000000       0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch4', 'reagent')] = qcdb.Molecule("""
0 1
C          0.00000000       0.00000000       0.00000000
H          1.08613677       0.00000000       0.00000000
H         -0.36204538       1.02401965       0.00000000
H         -0.36204538      -0.51200982      -0.88682703
H         -0.36204538      -0.51200982       0.88682703
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'c2h3', 'reagent')] = qcdb.Molecule("""
0 2
C          0.02607811       0.69573733       0.00000000
C          0.02857044      -0.62089246       0.00000000
H         -0.70268837       1.48374496       0.00000000
H         -0.89678124      -1.18847800       0.00000000
H          0.94877872      -1.18643193       0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3cn', 'reagent')] = qcdb.Molecule("""
0 1
H          1.55303000      -0.19362000       1.00615000
C          1.17602000       0.00000000       0.00000000
C         -0.28096000       0.00000000      -0.00001000
N         -1.43280000       0.00000000       0.00000000
H          1.55307000       0.96815000      -0.33536000
H          1.55311000      -0.77453000      -0.67072000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'allyl', 'reagent')] = qcdb.Molecule("""
0 2
H          1.29620000       1.27815000       0.00032000
C          1.22748000       0.19570000       0.00008000
C          0.00002000      -0.44151000      -0.00012000
C         -1.22747000       0.19573000      -0.00006000
H         -2.15453000      -0.36307000       0.00092000
H          2.15460000      -0.36295000      -0.00001000
H         -0.00013000      -1.52978000      -0.00012000
H         -1.29630000       1.27818000      -0.00046000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3cch', 'reagent')] = qcdb.Molecule("""
0 1
H          1.62987000      -0.22419000       0.99609000
C          1.23812000       0.00000000      -0.00001000
C         -0.21922000      -0.00003000       0.00000000
C         -1.42012000       0.00006000       0.00007000
H         -2.48217000      -0.00018000      -0.00029000
H          1.62978000       0.97478000      -0.30390000
H          1.62986000      -0.75054000      -0.69223000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3', 'reagent')] = qcdb.Molecule("""
0 2
C          0.00000000       0.00000000       0.00000031
H         -0.00000000       0.00000000       1.07554864
H         -0.00000000       0.93145124      -0.53777618
H          0.00000000      -0.93145124      -0.53777618
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'h2', 'reagent')] = qcdb.Molecule("""
0 1
H          0.00000000       0.00000000      -0.37169941
H          0.00000000       0.00000000       0.37169941
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3chch2', 'reagent')] = qcdb.Molecule("""
0 1
H         -1.80765000      -0.15368000      -0.87813000
C         -1.23352000       0.16234000       0.00003000
C          0.13465000      -0.45362000      -0.00010000
C          1.28048000       0.22043000      -0.00006000
H          2.23888000      -0.28641000       0.00025000
H         -1.18210000       1.25365000      -0.00053000
H         -1.80711000      -0.15289000       0.87881000
H          0.16668000      -1.54201000       0.00023000
H          1.30167000       1.30642000       0.00021000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'chfch2', 'reagent')] = qcdb.Molecule("""
0 1
F          1.15776000      -0.22307000       0.00001000
C         -0.02076000       0.43302000      -0.00002000
C         -1.18005000      -0.19843000       0.00000000
H          0.11660000       1.50810000       0.00003000
H         -2.09896000       0.37141000       0.00003000
H         -1.23262000      -1.27944000      -0.00001000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'h', 'reagent')] = qcdb.Molecule("""
0 2
H          0.00000000       0.00000000       0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'isobutane', 'reagent')] = qcdb.Molecule("""
0 1
C          0.00000090      -0.36984652       0.00000000
H          0.00000133      -1.47997543       0.00000000
C         -1.46235985       0.10584355       0.00000000
C          0.73117989       0.10584959       1.26644108
C          0.73117989       0.10584959      -1.26644108
H         -1.50901508       1.21298165       0.00000000
H          0.75449342       1.21298769       1.30684890
H          0.75449342       1.21298769      -1.30684890
H         -2.00230292      -0.25604001       0.89527086
H         -2.00230292      -0.25604001      -0.89527086
H          1.77648250      -0.25602315       1.28640747
H          0.22582896      -0.25604150       2.18168063
H          0.22582896      -0.25604150      -2.18168063
H          1.77648250      -0.25602315      -1.28640747
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2cho', 'reagent')] = qcdb.Molecule("""
0 2
H          1.26423000      -1.25075000      -0.00022000
C          1.16786000      -0.17142000      -0.00004000
C         -0.13387000       0.40647000      -0.00006000
H         -0.18501000       1.51210000      -0.00018000
H          2.05914000       0.44514000       0.00051000
O         -1.16779000      -0.26460000       0.00006000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2cch', 'reagent')] = qcdb.Molecule("""
0 2
C          0.11561000      -0.00003000      -0.00001000
C          1.33791000       0.00003000      -0.00001000
H          2.40011000      -0.00005000       0.00008000
C         -1.25132000       0.00001000      -0.00001000
H         -1.80663000       0.93004000       0.00004000
H         -1.80669000      -0.93000000       0.00004000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2ccn', 'reagent')] = qcdb.Molecule("""
0 2
C         -1.80740000       0.05630000       0.00027000
C         -0.52077000      -0.14829000      -0.00057000
C          0.80722000      -0.01017000       0.00028000
N          1.98158000       0.04542000       0.00000000
H         -2.22720000       1.06314000      -0.00029000
H         -2.51811000      -0.76812000       0.00035000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'c6h6', 'reagent')] = qcdb.Molecule("""
0 1
C         -0.00000000       1.41066086       0.00000000
C          1.22166790       0.70533038      -0.00000000
C          1.22166790      -0.70533038       0.00000000
C          0.00000000      -1.41066086       0.00000000
C         -1.22166790      -0.70533038       0.00000000
C         -1.22166790       0.70533038       0.00000000
H         -0.00000000       2.50726822       0.00000000
H          2.17135800       1.25363362       0.00000000
H          2.17135800      -1.25363362       0.00000000
H          0.00000000      -2.50726822       0.00000000
H         -2.17135800      -1.25363362       0.00000000
H         -2.17135800       1.25363362       0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'c6h5', 'reagent')] = qcdb.Molecule("""
0 2
C         -0.02911138       1.44968932       0.00000000
C          1.19145514       0.72918056       0.00000000
C          1.18256095      -0.68275137       0.00000000
C         -0.03576037      -1.39642075       0.00000000
C         -1.27002298      -0.69963430       0.00000000
C         -1.20569095       0.69610546       0.00000000
H         -0.03893619       2.54557277       0.00000000
H          2.14438286       1.27284078       0.00000000
H          2.13198633      -1.23090189       0.00000000
H         -0.03011941      -2.49351046       0.00000000
H         -2.22399795      -1.23906809       0.00000000
units angstrom
""")

