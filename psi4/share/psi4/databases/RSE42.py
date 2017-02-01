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
| Database of radical stabilization energies.
| Geometries from [E. Soydas and U. Bozkaya, JCTC, 9, 1452-1460 (2013)].
| Reference radical stabilization energies from [E. Soydas and U. Bozkaya, JCTC, 9, 1452-1460 (2013)] at CCSD(T)/cc-pCVTZ level.


- **benchmark**

  - ``'RSE42'`` [E. Soydas and U. Bozkaya, JCTC, 9, 1452-1460 (2013)].
  - |dl| ``'RSE42'`` |dr| [E. Soydas and U. Bozkaya, JCTC, 9, 1452-1460 (2013)].

- **subset**

  - ``'small'`` <members_description>
  - ``'large'`` <members_description>
  - ``'RSE30'`` smaller systems in RSE42
  - ``'<subset>'`` <members_description>

"""
import re
import qcdb

# <<< RSE42 Database Module >>>
dbse = 'RSE42'
isOS = 'True'

# <<< Database Members >>>
HRXN = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42']
HRXN_SM = []
HRXN_LG = []
RSE30 = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30']

# <<< Chemical Systems Involved >>>
RXNM = {}     # reaction matrix of reagent contributions per reaction
ACTV = {}     # order of active reagents per reaction
ACTV['%s-%s'            % (dbse, '1'                     )] = ['%s-%s-reagent'      % (dbse, 'ch3no2'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2no2'),
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '1'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '1')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '2'                     )] = ['%s-%s-reagent'      % (dbse, 'ch3ocho'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2ocho'),
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '2'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '2')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '3'                     )] = ['%s-%s-reagent'      % (dbse, 'ch3sch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2sch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '3'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '3')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '4'                     )] = ['%s-%s-reagent'      % (dbse, 'cfhch2'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'cfch2'),
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '4'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '4')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '5'                     )] = ['%s-%s-reagent'      % (dbse, 'ch3ch2f'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2ch2f'),
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '5'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '5')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '6'                     )] = ['%s-%s-reagent'      % (dbse, 'ch3cho'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2cho'),
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '6'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '6')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '7'                     )] = ['%s-%s-reagent'      % (dbse, 'ch3cn'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2cn'),
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '7'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '7')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '8'                     )] = ['%s-%s-reagent'      % (dbse, 'ch3f'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2f'),
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '8'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '8')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '9'                     )] = ['%s-%s-reagent'      % (dbse, 'ch3nh2'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2nh2'),
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '9'                     )] = dict(zip(ACTV['%s-%s' % (dbse, '9')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '10'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3nh3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2nh3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '10'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '10')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '11'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3nhoh'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2nhoh'),
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '11'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '11')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '12'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3oh'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2oh'),
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '12'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '12')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '13'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3ph3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2ph3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '13'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '13')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '14'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3sh2'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2sh2'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '14'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '14')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '15'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3sh'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2sh'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '15'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '15')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '16'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3cch'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2cch'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '16'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '16')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '17'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2ch3'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '17'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '17')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '18'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3cl'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2cl'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '18'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '18')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '19'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3bh2'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2bh2'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '19'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '19')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '20'                    )] = ['%s-%s-reagent'      % (dbse, 'ch2o'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'cho'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '20'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '20')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '21'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3ph2'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2ph2'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '21'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '21')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '22'                    )] = ['%s-%s-reagent'      % (dbse, 'ch2clf'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'chclf'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '22'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '22')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '23'                    )] = ['%s-%s-reagent'      % (dbse, 'ch2fch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'chfch3'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '23'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '23')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '24'                    )] = ['%s-%s-reagent'      % (dbse, 'ch2ohoh'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'chohoh'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '24'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '24')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '25'                    )] = ['%s-%s-reagent'      % (dbse, 'ch2cl2'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'chcl2'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '25'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '25')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '26'                    )] = ['%s-%s-reagent'      % (dbse, 'ch2f2'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'chf2'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '26'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '26')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '27'                    )] = ['%s-%s-reagent'      % (dbse, 'ch2chcn'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2ccn'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '27'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '27')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '28'                    )] = ['%s-%s-reagent'      % (dbse, 'c2h2'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'hcc'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '28'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '28')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '29'                    )] = ['%s-%s-reagent'      % (dbse, 'c2h4'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'c2h3'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '29'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '29')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '30'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3chch2'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2chch2'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '30'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '30')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '31'                    )] = ['%s-%s-reagent'      % (dbse, 'cyclopropane'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'cyclopropyl'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '31'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '31')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '32'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3ch2cl'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2ch2cl'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '32'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '32')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '33'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3ch2oh'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2ch2oh'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '33'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '33')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '34'                    )] = ['%s-%s-reagent'      % (dbse, 'methylcyclopropane'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'cyclopropylmethyl'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '34'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '34')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '35'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3coch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2coch3'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '35'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '35')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '36'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3conh2'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2conh2'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '36'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '36')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '37'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3cooh'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2cooh'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '37'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '37')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '38'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3nhch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2nhch3'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '38'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '38')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '39'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3nhcho'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2nhcho'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '39'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '39')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '40'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3och3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch2och3'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '40'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '40')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '41'                    )] = ['%s-%s-reagent'      % (dbse, 'nh2ch2cn'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'nh2chcn'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '41'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '41')], [-1,-1,1,1]))

ACTV['%s-%s'            % (dbse, '42'                    )] = ['%s-%s-reagent'      % (dbse, 'ch3chch2'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3'),
                                                               '%s-%s-reagent'      % (dbse, 'ch3cch2'), 
                                                               '%s-%s-reagent'      % (dbse, 'ch4') ]
RXNM['%s-%s'            % (dbse, '42'                    )] = dict(zip(ACTV['%s-%s' % (dbse, '42')], [-1,-1,1,1]))



# <<< Reference Values [kcal/mol] >>>
BIND = {}
BIND['%s-%s'            % (dbse, '1'                     )] =   -3.18
BIND['%s-%s'            % (dbse, '2'                     )] =   -4.74
BIND['%s-%s'            % (dbse, '3'                     )] =  -10.86
BIND['%s-%s'            % (dbse, '4'                     )] =    6.72
BIND['%s-%s'            % (dbse, '5'                     )] =   -1.49
BIND['%s-%s'            % (dbse, '6'                     )] =   -9.65
BIND['%s-%s'            % (dbse, '7'                     )] =   -8.20
BIND['%s-%s'            % (dbse, '8'                     )] =   -4.17
BIND['%s-%s'            % (dbse, '9'                     )] =  -12.00
BIND['%s-%s'            % (dbse, '10'                    )] =    4.59
BIND['%s-%s'            % (dbse, '11'                    )] =   -8.67
BIND['%s-%s'            % (dbse, '12'                    )] =   -9.21
BIND['%s-%s'            % (dbse, '13'                    )] =    0.55
BIND['%s-%s'            % (dbse, '14'                    )] =    2.35
BIND['%s-%s'            % (dbse, '15'                    )] =   -9.56 
BIND['%s-%s'            % (dbse, '16'                    )] =  -12.66
BIND['%s-%s'            % (dbse, '17'                    )] =   -3.32
BIND['%s-%s'            % (dbse, '18'                    )] =   -5.58
BIND['%s-%s'            % (dbse, '19'                    )] =  -11.61
BIND['%s-%s'            % (dbse, '20'                    )] =  -17.41
BIND['%s-%s'            % (dbse, '21'                    )] =   -6.39
BIND['%s-%s'            % (dbse, '22'                    )] =   -6.48
BIND['%s-%s'            % (dbse, '23'                    )] =   -5.79
BIND['%s-%s'            % (dbse, '24'                    )] =   -6.59
BIND['%s-%s'            % (dbse, '25'                    )] =   -9.41
BIND['%s-%s'            % (dbse, '26'                    )] =   -3.97
BIND['%s-%s'            % (dbse, '27'                    )] =    3.78
BIND['%s-%s'            % (dbse, '28'                    )] =   27.08
BIND['%s-%s'            % (dbse, '29'                    )] =    5.87
BIND['%s-%s'            % (dbse, '30'                    )] =  -17.05
BIND['%s-%s'            % (dbse, '31'                    )] =    0.000
BIND['%s-%s'            % (dbse, '32'                    )] =    0.000
BIND['%s-%s'            % (dbse, '33'                    )] =    0.000
BIND['%s-%s'            % (dbse, '34'                    )] =    0.000
BIND['%s-%s'            % (dbse, '35'                    )] =    0.000
BIND['%s-%s'            % (dbse, '36'                    )] =    0.000
BIND['%s-%s'            % (dbse, '37'                    )] =    0.000
BIND['%s-%s'            % (dbse, '38'                    )] =    0.000
BIND['%s-%s'            % (dbse, '39'                    )] =    0.000
BIND['%s-%s'            % (dbse, '40'                    )] =    0.000
BIND['%s-%s'            % (dbse, '41'                    )] =    0.000
BIND['%s-%s'            % (dbse, '42'                    )] =    0.000

# <<< Comment Lines >>>
TAGL = {}
TAGL['%s-%s'            % (dbse, '1'                     )] = """Reaction 1 """
TAGL['%s-%s'            % (dbse, '2'                     )] = """Reaction 2 """
TAGL['%s-%s'            % (dbse, '3'                     )] = """Reaction 3 """
TAGL['%s-%s'            % (dbse, '4'                     )] = """Reaction 4 """
TAGL['%s-%s'            % (dbse, '5'                     )] = """Reaction 5 """
TAGL['%s-%s'            % (dbse, '6'                     )] = """Reaction 6 """
TAGL['%s-%s'            % (dbse, '7'                     )] = """Reaction 7 """
TAGL['%s-%s'            % (dbse, '8'                     )] = """Reaction 8 """
TAGL['%s-%s'            % (dbse, '9'                     )] = """Reaction 9 """
TAGL['%s-%s'            % (dbse, '10'                    )] = """Reaction 10 """
TAGL['%s-%s'            % (dbse, '11'                    )] = """Reaction 11 """
TAGL['%s-%s'            % (dbse, '12'                    )] = """Reaction 12 """
TAGL['%s-%s'            % (dbse, '13'                    )] = """Reaction 13 """
TAGL['%s-%s'            % (dbse, '14'                    )] = """Reaction 14 """
TAGL['%s-%s'            % (dbse, '15'                    )] = """Reaction 15 """
TAGL['%s-%s'            % (dbse, '16'                    )] = """Reaction 16 """
TAGL['%s-%s'            % (dbse, '17'                    )] = """Reaction 17 """
TAGL['%s-%s'            % (dbse, '18'                    )] = """Reaction 18 """
TAGL['%s-%s'            % (dbse, '19'                    )] = """Reaction 19 """
TAGL['%s-%s'            % (dbse, '20'                    )] = """Reaction 20 """
TAGL['%s-%s'            % (dbse, '21'                    )] = """Reaction 21 """
TAGL['%s-%s'            % (dbse, '22'                    )] = """Reaction 22 """
TAGL['%s-%s'            % (dbse, '23'                    )] = """Reaction 23 """
TAGL['%s-%s'            % (dbse, '24'                    )] = """Reaction 24 """
TAGL['%s-%s'            % (dbse, '25'                    )] = """Reaction 25 """
TAGL['%s-%s'            % (dbse, '26'                    )] = """Reaction 26 """
TAGL['%s-%s'            % (dbse, '27'                    )] = """Reaction 27 """
TAGL['%s-%s'            % (dbse, '28'                    )] = """Reaction 28 """
TAGL['%s-%s'            % (dbse, '29'                    )] = """Reaction 29 """
TAGL['%s-%s'            % (dbse, '30'                    )] = """Reaction 30 """
TAGL['%s-%s'            % (dbse, '31'                    )] = """Reaction 31 """
TAGL['%s-%s'            % (dbse, '32'                    )] = """Reaction 32 """
TAGL['%s-%s'            % (dbse, '33'                    )] = """Reaction 33 """
TAGL['%s-%s'            % (dbse, '34'                    )] = """Reaction 34 """
TAGL['%s-%s'            % (dbse, '35'                    )] = """Reaction 35 """
TAGL['%s-%s'            % (dbse, '36'                    )] = """Reaction 36 """
TAGL['%s-%s'            % (dbse, '37'                    )] = """Reaction 37 """
TAGL['%s-%s'            % (dbse, '38'                    )] = """Reaction 38 """
TAGL['%s-%s'            % (dbse, '39'                    )] = """Reaction 39 """
TAGL['%s-%s'            % (dbse, '40'                    )] = """Reaction 40 """
TAGL['%s-%s'            % (dbse, '41'                    )] = """Reaction 41 """
TAGL['%s-%s'            % (dbse, '42'                    )] = """Reaction 42 """
TAGL['%s-%s-reagent'    % (dbse, 'ch2clf'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2fch3'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2chch2'              )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3f'                  )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'nh2ch2cn'              )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2coch3'              )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3cl'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2ch2cl'              )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3ch2cl'              )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'cyclopropylmethyl'     )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'chohoh'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3nhch3'              )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3cch'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'cfch2'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3bh2'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2cl'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'cyclopropane'          )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3ocho'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2f'                  )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2ccn'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2nhoh'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'methylcyclopropane'    )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3ch2f'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2chcn'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3ph2'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'chfch3'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'c2h2'                  )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3cho'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3cch2'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2cho'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2nhch3'              )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2bh2'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3nh2'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3cn'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3nhcho'              )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2ph2'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2conh2'              )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2ocho'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3'                   )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3conh2'              )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'c2h3'                  )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2ch2oh'              )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3sh2'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'hcc'                   )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2ohoh'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2o'                  )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2cl2'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3sch3'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2cooh'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2ch2f'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2f2'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'chclf'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch4'                   )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2nh2'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'chf2'                  )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3nhoh'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3chch2'              )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3coch3'              )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'cfhch2'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3oh'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2sh2'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2ph3'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3sh'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3och3'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3ch2oh'              )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2cch'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'c2h4'                  )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2sh'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'cho'                   )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3cooh'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2cn'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2ch3'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2nhcho'              )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2oh'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2nh3'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3ph3'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2och3'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'nh2chcn'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2sch3'               )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'chcl2'                 )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3nh3'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3ch3'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch3no2'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'ch2no2'                )] = """ """
TAGL['%s-%s-reagent'    % (dbse, 'cyclopropyl'           )] = """ """

# <<< Geometry Specification Strings >>>
GEOS = {}

GEOS['%s-%s-%s' % (dbse, 'ch2clf', 'reagent')] = qcdb.Molecule("""
0 1
C          0.58067800       0.57509900       0.00000000
H          0.68300400       1.16697100       0.90686900
H          0.68300300       1.16698200      -0.90686000
F          1.50669700      -0.42314600       0.00000000
CL        -1.08296100      -0.11624900       0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2fch3', 'reagent')] = qcdb.Molecule("""
0 1
C         -0.10731600       0.55006900      -0.00001800
C          1.19344700      -0.22296700      -0.00006600
H         -0.20208100       1.17829900       0.89066500
H         -0.20250500       1.17945300      -0.88990400
H          1.26426900      -0.85757700      -0.88617900
H          2.04390700       0.46538200       0.00014000
H          1.26409900      -0.85777600       0.88602700
F         -1.18716300      -0.34115500      -0.00002700
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2chch2', 'reagent')] = qcdb.Molecule("""
0 2
H          1.29620400       1.27814700       0.00031500
C          1.22747600       0.19569600       0.00007500
C          0.00002000      -0.44151000      -0.00011900
C         -1.22747000       0.19572600      -0.00006200
H         -2.15453000      -0.36306700       0.00092200
H          2.15459900      -0.36294700      -0.00001300
H         -0.00012600      -1.52977700      -0.00012200
H         -1.29630000       1.27817500      -0.00046200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3f', 'reagent')] = qcdb.Molecule("""
0 1
C         -0.63474800       0.00000000      -0.00000500
H         -0.99345300       0.93604800      -0.43577500
H         -0.99347500      -0.09061900       1.02851700
H         -0.99345500      -0.84542400      -0.59273200
F          0.75431900       0.00000000       0.00000200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'nh2ch2cn', 'reagent')] = qcdb.Molecule("""
0 1
H         -2.42292000      -0.10756500      -0.17398600
N         -1.48513000      -0.48905400      -0.11388600
C         -0.52267600       0.60284100       0.04818300
C          0.85402700       0.09830600      -0.00779600
N          1.93691300      -0.29425300      -0.02269800
H         -1.45838500      -1.10667700       0.69080300
H         -0.64927500       1.30540400      -0.78074200
H         -0.62000100       1.18510500       0.97768400
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2coch3', 'reagent')] = qcdb.Molecule("""
0 2
C         -1.12950700      -0.87211300      -0.00005500
C         -0.08591400       0.11935900      -0.00029700
C          1.35884600      -0.35579400      -0.00012600
H         -2.15866300      -0.53573300       0.00114600
H         -0.91846000      -1.93531800      -0.00052700
O         -0.36674800       1.31947500      -0.00001600
H          2.02124300       0.50840600      -0.01441200
H          1.56766900      -0.95785600       0.88952900
H          1.56164700      -0.98401000      -0.87274100
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3cl', 'reagent')] = qcdb.Molecule("""
0 1
C          1.14172600       0.00001000      -0.00001100
CL        -0.66455600      -0.00000200       0.00000300
H          1.48242800      -0.15709100       1.02027600
H          1.48235400      -0.80512100      -0.64615700
H          1.48231600       0.96218500      -0.37410500
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2ch2cl', 'reagent')] = qcdb.Molecule("""
0 2
C          1.65215500      -0.36546200      -0.00002000
C          0.61456300       0.65491000       0.00017500
CL        -1.10301600      -0.15556000      -0.00007700
H          1.97726100      -0.81625800       0.92878800
H          1.97884400      -0.81446300      -0.92913600
H          0.59733200       1.26901200       0.89634600
H          0.59752100       1.26954500      -0.89561800
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3ch2cl', 'reagent')] = qcdb.Molecule("""
0 1
C         -1.61792700      -0.35828400       0.00000300
C         -0.49398800       0.65857000       0.00000500
CL         1.14228200      -0.15032300      -0.00000200
H         -2.57929200       0.16589400      -0.00142700
H         -1.57402700      -0.99560000      -0.88502900
H         -1.57568300      -0.99386300       0.88636700
H         -0.50916200       1.28871900      -0.88740600
H         -0.50913600       1.28862100       0.88748100
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'cyclopropylmethyl', 'reagent')] = qcdb.Molecule("""
0 2
C          1.57690500       0.00025700      -0.16395700
C          0.27580500      -0.00015600       0.48454100
C         -0.91190300       0.74706100      -0.13393400
C         -0.91183900      -0.74701500      -0.13432800
H          1.65509200       0.00064000      -1.24484900
H          2.49026900      -0.00112700       0.41491300
H          0.29545600      -0.00050400       1.56912100
H         -1.57797300       1.26783500       0.54415100
H         -0.72934000       1.26136600      -1.07000800
H         -1.57793900      -1.26822500       0.54341400
H         -0.72937600      -1.26087000      -1.07067900
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'chohoh', 'reagent')] = qcdb.Molecule("""
0 2
C          0.00688400       0.51186100      -0.14492200
O         -1.14672000      -0.22110200      -0.07812600
O          1.17471800      -0.16554500       0.08128400
H          0.02457000       1.50148200       0.30721800
H          1.07671900      -1.04466600      -0.30624400
H         -1.36657400      -0.43480400       0.84328900
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3nhch3', 'reagent')] = qcdb.Molecule("""
0 1
H          2.08967400       0.42785800      -0.06224700
C          1.21538900      -0.22241300       0.02028200
N          0.00000700       0.56360400      -0.14844600
C         -1.21542600      -0.22239500       0.02030000
H         -2.08960200       0.42801100      -0.06218200
H          1.28226700      -0.96647500      -0.77983000
H          1.27608700      -0.76437700       0.98168500
H          0.00021800       1.32940400       0.51640300
H         -1.28238800      -0.96637700      -0.77987600
H         -1.27608400      -0.76442700       0.98167300
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3cch', 'reagent')] = qcdb.Molecule("""
0 1
H          1.62986900      -0.22418700       0.99608700
C          1.23811600      -0.00000200      -0.00001400
C         -0.21921700      -0.00003100       0.00000300
C         -1.42012300       0.00005600       0.00006700
H         -2.48216700      -0.00018300      -0.00029300
H          1.62978000       0.97477500      -0.30390300
H          1.62985800      -0.75053900      -0.69223100
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'cfch2', 'reagent')] = qcdb.Molecule("""
0 2
F          1.18236000       0.12888300       0.00008300
C         -0.00743100      -0.42937300      -0.00017000
C         -1.20008000       0.11937700      -0.00003700
H         -2.08422800      -0.50156000       0.00065000
H         -1.31193900       1.20158500      -0.00015200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3bh2', 'reagent')] = qcdb.Molecule("""
0 1
C         -0.68214300      -0.00001300      -0.01636900
B          0.87205600      -0.00003500      -0.02090300
H         -1.14585100      -0.89605500      -0.43613100
H         -0.94870400       0.00006100       1.05694600
H         -1.14545300       0.89608700      -0.43641500
H          1.48609900       1.02439000       0.00926200
H          1.48648600      -1.02423100       0.00906800
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2cl', 'reagent')] = qcdb.Molecule("""
0 2
C          1.12668100       0.00000000      -0.00045200
CL        -0.58859200       0.00000000       0.00003600
H          1.62298700      -0.95630200       0.00105100
H          1.62298700       0.95630300       0.00105100
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'cyclopropane', 'reagent')] = qcdb.Molecule("""
0 1
C          0.26222500       0.82955100       0.00000100
C          0.58730900      -0.64189800       0.00000100
C         -0.84955700      -0.18773900      -0.00000300
H          0.43991900       1.39156800      -0.90914700
H          0.43990100       1.39158200       0.90914400
H         -1.42517500      -0.31472400      -0.90913600
H         -1.42515800      -0.31472000       0.90914200
H          0.98531900      -1.07660600       0.90916100
H          0.98533000      -1.07658600      -0.90916300
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3ocho', 'reagent')] = qcdb.Molecule("""
0 1
C          1.36739400       0.41209300       0.00000000
O          0.00000000       0.87018200       0.00000000
C         -0.93310600      -0.09275600       0.00000000
H         -1.92702500       0.37721300       0.00000000
H          1.97602500       1.31296200       0.00000000
H          1.56764700      -0.18831200       0.88859200
H          1.56764700      -0.18831200      -0.88859200
O         -0.72375300      -1.27387900       0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2f', 'reagent')] = qcdb.Molecule("""
0 2
C          0.65517400       0.00000000      -0.05820600
H          1.12122500      -0.95818100       0.12958600
H          1.12122500       0.95818100       0.12958600
F         -0.68594400       0.00000000       0.01000700
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2ccn', 'reagent')] = qcdb.Molecule("""
0 2
C         -1.80740300       0.05630000       0.00026900
C         -0.52077400      -0.14829300      -0.00056500
C          0.80722000      -0.01016700       0.00027900
N          1.98157900       0.04542000       0.00000400
H         -2.22720000       1.06313700      -0.00028800
H         -2.51811100      -0.76812000       0.00035400
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2nhoh', 'reagent')] = qcdb.Molecule("""
0 2
H         -2.10110400       0.31208400      -0.21944500
C         -1.19615900      -0.26035100      -0.06713800
N         -0.05044800       0.47182000       0.18732900
O          1.15271000      -0.24272900      -0.16818500
H          1.66921300      -0.19874400       0.64455800
H         -1.23356900      -1.26401500       0.33245000
H         -0.02612900       1.35187600      -0.32055700
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'methylcyclopropane', 'reagent')] = qcdb.Molecule("""
0 1
H          2.12855200       0.88463900       0.15688700
C          1.55441500       0.00000200      -0.13663300
C          0.18108800      -0.00019200       0.49815600
C         -0.95723500       0.75588200      -0.13990700
C         -0.95721600      -0.75573000      -0.14010200
H          1.48094300       0.00056200      -1.22847700
H          2.12860600      -0.88478000       0.15624200
H          0.18845400      -0.00024300       1.58445700
H         -1.66438600       1.26919000       0.50089000
H         -0.76177200       1.26046600      -1.07978800
H         -1.66416000      -1.26924700       0.50085800
H         -0.76255500      -1.26035500      -1.08014900
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3ch2f', 'reagent')] = qcdb.Molecule("""
0 1
H          2.04393800       0.46531400      -0.00025500
C          1.19346600      -0.22298200       0.00000300
C         -0.10736700       0.55014500       0.00003600
F         -1.18710800      -0.34119200      -0.00002700
H          1.26381800      -0.85762100      -0.88615200
H          1.26425400      -0.85768800       0.88610200
H         -0.20233100       1.17862500       0.89046700
H         -0.20230000       1.17911900      -0.89014800
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2chcn', 'reagent')] = qcdb.Molecule("""
0 1
H         -2.62698000       0.00343300       0.00127500
C         -1.60537600      -0.35564600       0.00019100
C         -0.58353400       0.50215200      -0.00011700
C          0.78295500       0.08952700      -0.00090900
N          1.89484100      -0.22375700       0.00060300
H         -1.45041800      -1.42770800      -0.00098200
H         -0.75075800       1.57437800       0.00049800
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3ph2', 'reagent')] = qcdb.Molecule("""
0 1
C         -1.19705600      -0.00003100       0.02529600
P          0.66982000       0.00000200      -0.12443800
H         -1.58957300      -0.88074500      -0.48670100
H         -1.55478100      -0.00028500       1.05489100
H         -1.58939300       0.88107700      -0.48619200
H          0.93418900       1.03642300       0.81636300
H          0.93459700      -1.03630900       0.81643200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'chfch3', 'reagent')] = qcdb.Molecule("""
0 2
C          0.11105500       0.51607700      -0.10093000
C         -1.19373800      -0.17214000       0.01228500
H          0.29639400       1.52732100       0.24560500
H         -1.35692900      -0.58834000       1.01911100
H         -2.00676900       0.52649400      -0.19767600
H         -1.25878800      -1.00450200      -0.69563600
F          1.20246600      -0.28051100       0.01783000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'c2h2', 'reagent')] = qcdb.Molecule("""
0 1
C         -0.00000000       0.00000000       0.60499861
C          0.00000000       0.00000000      -0.60499861
H          0.00000000      -0.00000000       1.66488377
H          0.00000000       0.00000000      -1.66488377
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3cho', 'reagent')] = qcdb.Molecule("""
0 1
H          1.15548300      -1.23749400      -0.00007400
C          1.16884900      -0.14775500      -0.00000800
C         -0.23563200       0.39720900      -0.00001200
H         -0.30508600       1.50871600       0.00003100
H          1.70764100       0.22227400       0.87912400
H          1.70778300       0.22244500      -0.87896500
O         -1.23314000      -0.27658300       0.00000100
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3cch2', 'reagent')] = qcdb.Molecule("""
0 2
H         -1.80737800      -0.28751200      -0.88203600
C         -1.27092200       0.07589700       0.00003200
C          0.12553700      -0.37933600      -0.00002800
C          1.34313900       0.09759700      -0.00018800
H          2.22165600      -0.54200900       0.00065800
H         -1.32823400       1.17576200      -0.00010500
H         -1.80727200      -0.28728400       0.88226200
H          1.53470900       1.17609600       0.00032000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2cho', 'reagent')] = qcdb.Molecule("""
0 2
H          1.26423000      -1.25074500      -0.00022300
C          1.16786100      -0.17141500      -0.00003900
C         -0.13387000       0.40646800      -0.00006100
H         -0.18501200       1.51209700      -0.00017600
H          2.05914100       0.44514300       0.00050900
O         -1.16778800      -0.26460200       0.00006200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2nhch3', 'reagent')] = qcdb.Molecule("""
0 2
C          1.25665700      -0.25296800       0.09101800
N          0.09359200       0.47138600      -0.13643100
C         -1.19132800      -0.18644800       0.03339200
H         -1.99251100       0.48192100      -0.28866900
H          1.25395100      -1.28237900      -0.24805200
H          2.18576000       0.30121100       0.04601900
H          0.11541700       1.41397000       0.22790100
H         -1.22447800      -1.07807600      -0.59874400
H         -1.38526200      -0.49985500       1.07010400
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2bh2', 'reagent')] = qcdb.Molecule("""
0 2
C          0.71330600       0.00000000       0.00000200
B         -0.81296800      -0.00000200      -0.00002100
H          1.30746700      -0.90739700      -0.08361900
H          1.30746400       0.90739800       0.08363400
H         -1.41496300       1.02725300      -0.07382900
H         -1.41496300      -1.02724800       0.07390900
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3nh2', 'reagent')] = qcdb.Molecule("""
0 1
H         -1.11810600      -0.87790900      -0.48746900
C         -0.70669300      -0.00000300       0.01773400
N          0.75219200      -0.00000100      -0.12374800
H         -1.11797400       0.87906800      -0.48553000
H         -1.08372700      -0.00111200       1.05295300
H          1.14730200       0.81201200       0.33996300
H          1.14731900      -0.81203100       0.33991600
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3cn', 'reagent')] = qcdb.Molecule("""
0 1
H          1.55303200      -0.19362100       1.00615000
C          1.17602200       0.00000000      -0.00000200
C         -0.28095500       0.00000200      -0.00001400
N         -1.43280200      -0.00000200       0.00000400
H          1.55307200       0.96815300      -0.33536300
H          1.55311000      -0.77452900      -0.67072400
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3nhcho', 'reagent')] = qcdb.Molecule("""
0 1
C          1.45060800      -0.43388200       0.00000600
N          0.47537300       0.64484300      -0.00000900
C         -0.86810900       0.42639900       0.00000000
H          2.08490300      -0.39459000      -0.89019800
H          2.08477400      -0.39465900       0.89029400
H          0.89750900      -1.37158600      -0.00008900
H          0.79723800       1.59973400       0.00000200
H         -1.43727300       1.37613100       0.00002100
O         -1.40621900      -0.66050400       0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2ph2', 'reagent')] = qcdb.Molecule("""
0 2
C          1.18935200       0.00005700       0.07043000
P         -0.58298700      -0.00011500      -0.12323900
H          1.74314200       0.92389300      -0.04461000
H          1.74250200      -0.92404900      -0.04565700
H         -0.93834000      -1.05920700       0.75926200
H         -0.93860700       1.06075100       0.75701500
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2conh2', 'reagent')] = qcdb.Molecule("""
0 2
H          2.18625900      -0.02366000       0.02261200
C          1.27446100      -0.60287200       0.00168800
C          0.01415400       0.12760900      -0.00043300
N         -1.13042200      -0.64022200      -0.03845600
H         -2.00633300      -0.15662700       0.08092200
H          1.32929100      -1.68477600      -0.02718100
H         -1.10934400      -1.62800200       0.15069600
O         -0.02732600       1.35327500       0.00432700
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2ocho', 'reagent')] = qcdb.Molecule("""
0 2
C          1.39696200      -0.44693100      -0.00027000
O          0.57750500       0.64559700       0.00004800
C         -0.77182900       0.44666800       0.00003400
H         -1.26113100       1.42735100      -0.00011500
H          0.95655300      -1.43028600       0.00054600
H          2.44597000      -0.20760900       0.00061800
O         -1.31402900      -0.61908100      -0.00000200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3', 'reagent')] = qcdb.Molecule("""
0 2
C          0.00000000       0.00000000       0.00021400
H          0.00000000       1.08040900      -0.00042800
H          0.93566200      -0.54020500      -0.00042800
H         -0.93566200      -0.54020500      -0.00042800
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3conh2', 'reagent')] = qcdb.Molecule("""
0 1
H          1.90246300       0.16251900      -0.79962900
C          1.36161600      -0.34591400      -0.00025200
C         -0.07735100       0.14833800      -0.00371100
N         -1.03441600      -0.82787400      -0.00411700
H         -2.00130700      -0.54487700       0.01827000
H          1.46157800      -1.42571100      -0.12782200
H          1.82569300      -0.05631300       0.94525000
H         -0.81281100      -1.80771000       0.01025300
O         -0.35503700       1.33158200       0.00078500
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'c2h3', 'reagent')] = qcdb.Molecule("""
0 2
C          0.04825900      -0.58537400       0.00000000
C          0.04825900       0.71905900       0.00000000
H         -0.88001800      -1.16408100       0.00000000
H          0.96858300      -1.16591900       0.00000000
H         -0.66767700       1.52789000       0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2ch2oh', 'reagent')] = qcdb.Molecule("""
0 2
H         -2.18970300       0.23835000      -0.24005400
C         -1.25191200      -0.25471300      -0.01994900
C          0.00851200       0.51937100       0.03731500
O          1.11275700      -0.37712300      -0.06380000
H          1.91487400       0.11923700       0.12180500
H         -1.25306000      -1.28679500       0.30475400
H          0.02620700       1.27214500      -0.76868100
H          0.06002400       1.08609600       0.98837600
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3sh2', 'reagent')] = qcdb.Molecule("""
1 1
H         -1.51535400      -0.00004800       1.07239600
C         -1.21297000      -0.00000100       0.02856200
S          0.62777800       0.00000100      -0.11140700
H         -1.55326900      -0.89268700      -0.49425800
H         -1.55326500       0.89273500      -0.49417700
H          0.92762600       0.99336200       0.76359300
H          0.92763100      -0.99337500       0.76357900
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'hcc', 'reagent')] = qcdb.Molecule("""
0 2
C          0.00000000       0.00000000       0.00000000
C          1.21283562       0.00000000       0.00000000
H         -1.05818189       0.00000000       0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2ohoh', 'reagent')] = qcdb.Molecule("""
0 1
C         -0.00005100       0.53181800       0.00000000
O          1.16882200      -0.24677600       0.09147700
O         -1.16881300      -0.24681500      -0.09148400
H          0.00848900       1.16260500       0.89514600
H         -0.00857300       1.16266000      -0.89507300
H         -1.21965900      -0.78378900       0.70749200
H          1.21998200      -0.78365400      -0.70751500
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2o', 'reagent')] = qcdb.Molecule("""
0 1
C          0.00000100      -0.52592000       0.00000000
O          0.00000100       0.67402700       0.00000000
H          0.93860800      -1.11835000       0.00000000
H         -0.93861700      -1.11835000       0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2cl2', 'reagent')] = qcdb.Molecule("""
0 1
C          0.00000000       0.76838400       0.00000000
CL         1.49583700      -0.21639700       0.00000000
CL        -1.49583600      -0.21641300       0.00000000
H         -0.00000700       1.37372600       0.89921100
H         -0.00000700       1.37372600      -0.89921100
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3sch3', 'reagent')] = qcdb.Molecule("""
0 1
C          1.39450900      -0.51481700       0.00000000
S          0.00000000       0.66160700       0.00000000
C         -1.39450800      -0.51481700       0.00000000
H          2.30921400       0.07884900       0.00000000
H          1.38175000      -1.14140200       0.89387100
H          1.38175000      -1.14140200      -0.89387100
H         -2.30921400       0.07885100       0.00000000
H         -1.38175200      -1.14140000      -0.89387100
H         -1.38175200      -1.14140000       0.89387100
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2cooh', 'reagent')] = qcdb.Molecule("""
0 2
C          1.37652700      -0.27885400       0.00007400
C         -0.01480300       0.11096500       0.00004500
O         -0.86030200      -0.95687700      -0.00004800
H         -1.75682400      -0.59202100       0.00082400
H          1.66483000      -1.32063100      -0.00085900
H          2.12298100       0.50163400       0.00071600
O         -0.41486400       1.25917100      -0.00012700
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2ch2f', 'reagent')] = qcdb.Molecule("""
0 2
C         -1.22285000      -0.25059600      -0.02923000
C          0.03294900       0.52095200       0.06492100
F          1.14514700      -0.32233900      -0.06018000
H         -1.26498700      -1.24902400       0.38614400
H         -2.13334000       0.21667000      -0.38082800
H          0.09765700       1.28601300      -0.71691200
H          0.13375700       1.02525400       1.03907300
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2f2', 'reagent')] = qcdb.Molecule("""
0 1
C          0.00000400       0.50314900       0.00000000
F         -1.10821200      -0.29035800      -0.00000100
H          0.00000100       1.10379200      -0.91279600
H          0.00000100       1.10377600       0.91280700
F          1.10820900      -0.29036000      -0.00000100
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'chclf', 'reagent')] = qcdb.Molecule("""
0 2
C         -0.55167800       0.55013700      -0.13353900
H         -0.72586900       1.49900600       0.36337900
F         -1.52740700      -0.34238000       0.02747200
CL         1.04603500      -0.10108300       0.01121200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch4', 'reagent')] = qcdb.Molecule("""
0 1
C          0.00000000       0.00000000       0.00000000
H          0.62958700       0.62958700       0.62958700
H         -0.62958700      -0.62958700       0.62958700
H         -0.62958700       0.62958700      -0.62958700
H          0.62958700      -0.62958700      -0.62958700
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2nh2', 'reagent')] = qcdb.Molecule("""
0 2
H          1.24137600       0.92884200      -0.13169100
C          0.72878100      -0.00007100       0.08564100
N         -0.65633800      -0.00005500      -0.09925400
H          1.24163700      -0.92867500      -0.13216700
H         -1.13115900      -0.83256900       0.22273500
H         -1.13017600       0.83321100       0.22205700
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'chf2', 'reagent')] = qcdb.Molecule("""
0 2
C          0.00000000       0.48713600      -0.15110400
F         -1.10027700      -0.24185600       0.02841400
H          0.00000000       1.43059400       0.39518300
F          1.10027700      -0.24185600       0.02841400
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3nhoh', 'reagent')] = qcdb.Molecule("""
0 1
H          2.02509400       0.44185800      -0.06388800
C          1.16605500      -0.23232900       0.00456500
N         -0.04441600       0.57005800      -0.15285600
O         -1.17243400      -0.29041100       0.12920900
H         -1.68759600      -0.22935400      -0.68128200
H          1.22725100      -0.95421100      -0.81193400
H          1.20629000      -0.77610500       0.95778200
H         -0.07697800       1.24466600       0.60825200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3chch2', 'reagent')] = qcdb.Molecule("""
0 1
H         -1.80765000      -0.15368300      -0.87813100
C         -1.23352400       0.16233900       0.00002500
C          0.13465000      -0.45361900      -0.00010100
C          1.28048100       0.22043300      -0.00006200
H          2.23887900      -0.28641300       0.00024900
H         -1.18210100       1.25365400      -0.00052600
H         -1.80710900      -0.15289100       0.87880500
H          0.16667600      -1.54200800       0.00022700
H          1.30166500       1.30642200       0.00020600
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3coch3', 'reagent')] = qcdb.Molecule("""
0 1
C         -1.29144800      -0.61344900       0.00162800
C         -0.00000100       0.18639500      -0.00001300
C          1.29144900      -0.61344500      -0.00163200
H         -2.14125300       0.06176600       0.09074300
H         -1.37741200      -1.18620100      -0.92773100
H         -1.30085900      -1.33611700       0.82327900
O         -0.00000300       1.39550900       0.00000200
H          2.14124800       0.06175100      -0.09095500
H          1.37750300      -1.18597500       0.92785600
H          1.30079800      -1.33630500      -0.82311200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'cfhch2', 'reagent')] = qcdb.Molecule("""
0 1
F          1.15775800      -0.22307000       0.00000600
C         -0.02075700       0.43302300      -0.00002200
C         -1.18005100      -0.19842800       0.00000300
H          0.11660400       1.50809800       0.00003400
H         -2.09895700       0.37140500       0.00003100
H         -1.23261900      -1.27943700      -0.00000900
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3oh', 'reagent')] = qcdb.Molecule("""
0 1
H          1.08276700       0.98760000      -0.00028200
C          0.66362300      -0.01954600       0.00000100
O         -0.75009300       0.12169100       0.00000000
H         -1.13398500      -0.75944800       0.00001200
H          1.03509500      -0.54244500      -0.89192500
H          1.03513400      -0.54195800       0.89219400
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2sh2', 'reagent')] = qcdb.Molecule("""
1 2
C         -1.20662600       0.02170400       0.02219700
H         -1.68207400       0.97499900       0.20542400
H         -1.73384700      -0.90995700      -0.12920000
S          0.55425100      -0.07913500      -0.07976000
H          0.86827600       1.23974300      -0.17152600
H          0.91939200      -0.16885000       1.23828900
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2ph3', 'reagent')] = qcdb.Molecule("""
1 2
C         -1.24118800      -0.00000500      -0.00420900
P          0.52397100       0.00000300      -0.00112200
H         -1.77980500      -0.93960200       0.02212500
H         -1.77982900       0.93957700       0.02212200
H          1.03398600       1.13558400      -0.64768800
H          1.07918000      -0.00072100       1.29444700
H          1.03403900      -1.13485700      -0.64892800
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3sh', 'reagent')] = qcdb.Molecule("""
0 1
H         -1.52774200      -1.00777000      -0.00169900
C         -1.16516100       0.01984900       0.00004300
S          0.66721700      -0.08708800       0.00001100
H          0.90833000       1.23982400      -0.00009200
H         -1.53240700       0.52241100      -0.89383000
H         -1.53269700       0.51984700       0.89518300
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3och3', 'reagent')] = qcdb.Molecule("""
0 1
H          2.02193500      -0.48923200      -0.00071000
C          1.17294200       0.19529900       0.00000000
O          0.00000100      -0.58914400       0.00001900
C         -1.17294400       0.19529800       0.00000000
H         -2.02193900      -0.48923000      -0.00106500
H          1.23374500       0.83653400       0.89241600
H          1.23300600       0.83748100      -0.89177200
H         -1.23391300       0.83628100       0.89258300
H         -1.23283200       0.83773500      -0.89160200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3ch2oh', 'reagent')] = qcdb.Molecule("""
0 1
H          2.07592400       0.45870400      -0.00139100
C          1.22080500      -0.22279500       0.00033800
C         -0.08478500       0.55050000      -0.00015000
O         -1.15070400      -0.39819000      -0.00039200
H         -1.98159000       0.08549700       0.00300000
H          1.28351300      -0.86183700      -0.88367400
H          1.28455500      -0.85892900       0.88610000
H         -0.13737200       1.19906600       0.88662100
H         -0.13551500       1.19679300      -0.88864500
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2cch', 'reagent')] = qcdb.Molecule("""
0 2
C          0.11561300      -0.00003000      -0.00001300
C          1.33791000       0.00002600      -0.00000900
H          2.40011400      -0.00004900       0.00007800
C         -1.25132200       0.00000500      -0.00000500
H         -1.80663200       0.93004300       0.00004100
H         -1.80668700      -0.92999900       0.00004100
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'c2h4', 'reagent')] = qcdb.Molecule("""
0 1
C          0.00000000       0.66349000       0.00000000
C          0.00000000      -0.66349000       0.00000000
H          0.00000000       1.23461800       0.92256200
H          0.00000000       1.23461800      -0.92256200
H          0.00000000      -1.23461800      -0.92256200
H          0.00000000      -1.23461800       0.92256200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2sh', 'reagent')] = qcdb.Molecule("""
0 2
H          1.68655500      -0.90997600       0.00055800
C          1.14729300       0.02417600      -0.03176100
S         -0.58501000      -0.08911600       0.00912900
H         -0.86056300       1.22722500      -0.09184000
H          1.65041300       0.96355500       0.13578900
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'cho', 'reagent')] = qcdb.Molecule("""
0 2
C          0.06228700       0.58421000       0.00000000
O          0.06228700      -0.59005800       0.00000000
H         -0.87201400       1.21520500       0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3cooh', 'reagent')] = qcdb.Molecule("""
0 1
H         -1.91951200       0.83381100      -0.00001900
C         -1.39368200      -0.11802900      -0.00000300
C          0.09176300       0.12824400      -0.00001100
O          0.78506300      -1.03886100       0.00000700
H          1.72359800      -0.79889600       0.00001600
H         -1.67259500      -0.70149700       0.88023200
H         -1.67259900      -0.70153300      -0.88021200
O          0.63401500       1.20221400       0.00000100
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2cn', 'reagent')] = qcdb.Molecule("""
0 2
C          0.18727900      -0.00000100      -0.00000500
N          1.35587400       0.00000100       0.00000300
C         -1.19102700       0.00000000       0.00000000
H         -1.73431100       0.93521500       0.00000700
H         -1.73431500      -0.93521200       0.00000700
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2ch3', 'reagent')] = qcdb.Molecule("""
0 2
C          0.79404000      -0.00001000      -0.01898900
C         -0.69354400       0.00011200      -0.00171900
H          1.35165100      -0.92648600       0.04178700
H          1.35211600       0.92615600       0.04184900
H         -1.10700000       0.89016100      -0.48573900
H         -1.10738100      -0.88205600      -0.49991300
H         -1.09236600      -0.00838900       1.02626600
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2nhcho', 'reagent')] = qcdb.Molecule("""
0 2
C         -1.44716800      -0.49182200      -0.00025100
N         -0.57509300       0.57056600      -0.00005000
C          0.79897400       0.44788700       0.00008900
H         -2.50479900      -0.29017100       0.00146300
H         -1.02384400      -1.48129200       0.00008900
H         -0.95098500       1.50778300      -0.00055400
H          1.30076800       1.43105100       0.00039600
O          1.38670900      -0.61221600      -0.00000900
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2oh', 'reagent')] = qcdb.Molecule("""
0 2
H         -1.23081700      -0.88927600       0.10207100
C         -0.68449600       0.02764800      -0.06653200
O          0.67116100      -0.12516300       0.02350300
H          1.09168100       0.73511600      -0.06829800
H         -1.12317500       0.98957200       0.17739600
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2nh3', 'reagent')] = qcdb.Molecule("""
1 2
H          1.31094400       0.95958500       0.06254700
C          0.82661100      -0.00000600      -0.02780600
N         -0.64330400       0.00002700       0.00195700
H          1.31073200      -0.95972200       0.06250000
H         -1.01795500      -0.00400200       0.96403900
H         -1.03003400       0.82965200      -0.46454100
H         -1.03022800      -0.82566500      -0.47141100
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3ph3', 'reagent')] = qcdb.Molecule("""
1 1
H          1.58306100       0.47144600      -0.91637400
C          1.22355800      -0.00000600       0.00000000
P         -0.58438900       0.00000100       0.00000200
H          1.58314400      -1.02930900       0.04987600
H          1.58312400       0.55787600       0.86644700
H         -1.10827300      -0.59602900       1.15572600
H         -1.10823400       1.29890900      -0.06168600
H         -1.10833400      -0.70287400      -1.09402100
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2och3', 'reagent')] = qcdb.Molecule("""
0 2
C          1.20368100       0.22731000       0.07010300
O          0.09172500      -0.54123500      -0.03988900
C         -1.13925700       0.16821300       0.01374700
H         -1.93367000      -0.56449600      -0.12107400
H          2.12393300      -0.33163800      -0.02661000
H          1.14335000       1.27002200      -0.23072900
H         -1.25993800       0.66692800       0.98160700
H         -1.19401900       0.91592200      -0.78718500
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'nh2chcn', 'reagent')] = qcdb.Molecule("""
0 2
H         -2.45610400      -0.01444300       0.23937200
N         -1.55523600      -0.34528300      -0.07015500
C         -0.49903400       0.53856300       0.01366400
C          0.81999200       0.11101400       0.00107900
N          1.93043500      -0.24950800       0.00385500
H         -1.37038800      -1.31387700       0.14493400
H         -0.72565300       1.59439500      -0.00866900
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2sch3', 'reagent')] = qcdb.Molecule("""
0 2
C         -1.38974300       0.43014600       0.02091500
S          0.11017200      -0.60914300      -0.01904600
C          1.35489400       0.58185000       0.04179000
H         -2.24570300      -0.24284500      -0.02480000
H         -1.42650300       1.00472000       0.94638700
H         -1.41259000       1.09892200      -0.84044700
H          1.16313900       1.59665400      -0.27952900
H          2.36799900       0.21686200       0.12689900
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'chcl2', 'reagent')] = qcdb.Molecule("""
0 2
C          0.00000000       0.68873000      -0.09737200
CL        -1.48435900      -0.17167400       0.00930300
CL         1.48435900      -0.17167400       0.00930300
H          0.00000000       1.70452500       0.26792100
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3nh3', 'reagent')] = qcdb.Molecule("""
1 1
H          1.14365300      -0.16818600       1.01921200
C          0.80331800      -0.00000100      -0.00000600
N         -0.71232800       0.00000000       0.00000700
H          1.14367100       0.96675700      -0.36395500
H          1.14364200      -0.79859600      -0.65523500
H         -1.08821900      -0.89243500       0.33598100
H         -1.08817300       0.73722000       0.60485500
H         -1.08818700       0.15524400      -0.94087200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3ch3', 'reagent')] = qcdb.Molecule("""
0 1
C          0.00000000      -0.76522200       0.00000000
C          0.00000000       0.76522200       0.00000000
H          1.01839600      -1.16365400       0.00000000
H         -0.50924700      -1.16358500       0.88196300
H         -0.50924700      -1.16358500      -0.88196300
H         -1.01839600       1.16365300       0.00000000
H          0.50924900       1.16358400       0.88196200
H          0.50924900       1.16358400      -0.88196200
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch3no2', 'reagent')] = qcdb.Molecule("""
0 1
C                 -0.03891200   -1.32592100    0.00000000
H                  0.49240200   -1.65949700    0.88961400
H                  0.49240200   -1.65949700   -0.88961400
H                 -1.07696300   -1.64238400    0.00000000
O                 -1.06579500    0.76938100    0.00000000
O                  1.10649900    0.69073400    0.00000000
N                  0.00000000    0.17657000    0.00000000
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'ch2no2', 'reagent')] = qcdb.Molecule("""
0 2
C                 -1.30974800   -0.00015200   -0.00012900
H                 -1.79634400    0.96237700    0.00044400
H                 -1.79649800   -0.96262600    0.00011800
O                  0.66955200   -1.09726400    0.00001500
O                  0.66919400    1.09753900   -0.00002000
N                  0.10590800   -0.00014800    0.00003600
units angstrom
""")

GEOS['%s-%s-%s' % (dbse, 'cyclopropyl', 'reagent')] = qcdb.Molecule("""
0 2
C                     0.00000   0.87112  -0.16751 
C                     0.76811  -0.36353   0.0305 
C                    -0.76811  -0.36353   0.0305 
H                     0.00000   1.81695   0.35474 
H                    -1.25447  -0.82499  -0.82552 
H                    -1.29862  -0.51568   0.9677 
H                     1.29862  -0.51568   0.9677 
H                     1.25447  -0.82499  -0.82552 
units angstrom
""")

