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

try:
    from collections import OrderedDict
except ImportError:
    from .oldpymodules import OrderedDict

def sapt_psivars():
    """Returns dictionary of PsiVariable definitions.

    """
    pv1 = OrderedDict()
    pv1['SAPT EXCHSCAL1'] = {'func': lambda x: 1.0 if x[0] < 1.0e-5 else x[0] / x[1], 'args': ['SAPT EXCH10 ENERGY', 'SAPT EXCH10(S^2) ENERGY']}  # special treatment in pandas
    pv1['SAPT EXCHSCAL3'] = {'func': lambda x: x[0] ** 3, 'args': ['SAPT EXCHSCAL1']}
    pv1['SAPT EXCHSCAL'] = {'func': lambda x: x[0] ** x[1], 'args': ['SAPT EXCHSCAL1', 'SAPT ALPHA']}
    pv1['SAPT HF(2) ALPHA=0.0 ENERGY'] = {'func': lambda x: x[0] - (x[1] + x[2] + x[3] + x[4]),
                                          'args': ['SAPT HF TOTAL ENERGY', 'SAPT ELST10,R ENERGY', 'SAPT EXCH10 ENERGY',
                                                   'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY']}

    pv1['SAPT HF(2) ENERGY'] = {'func': lambda x: x[1] + (1.0 - x[0]) * x[2],
                                'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ALPHA=0.0 ENERGY', 'SAPT EXCH-IND20,R ENERGY']}
    pv1['SAPT HF(3) ENERGY'] = {'func': lambda x: x[1] - (x[2] + x[0] * x[3]),
                                'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ENERGY', 'SAPT IND30,R ENERGY', 'SAPT EXCH-IND30,R ENERGY']}
    pv1['SAPT MP2(2) ENERGY'] = {'func': lambda x: x[1] - (x[2] + x[3] + x[4] + x[0] * (x[5] + x[6] + x[7] + x[8])),
                                 'args': ['SAPT EXCHSCAL', 'SAPT MP2 CORRELATION ENERGY', 'SAPT ELST12,R ENERGY',  # MP2 CORRELATION ENERGY renamed here from pandas since this is IE  # renamed again SA --> SAPT
                                          'SAPT IND22 ENERGY', 'SAPT DISP20 ENERGY', 'SAPT EXCH11(S^2) ENERGY',
                                          'SAPT EXCH12(S^2) ENERGY', 'SAPT EXCH-IND22 ENERGY', 'SAPT EXCH-DISP20 ENERGY']}
    pv1['SAPT MP2(3) ENERGY'] = {'func': lambda x: x[1] - (x[2] + x[0] * x[3]),
                                 'args': ['SAPT EXCHSCAL', 'SAPT MP2(2) ENERGY', 'SAPT IND-DISP30 ENERGY', 'SAPT EXCH-IND-DISP30 ENERGY']}
    pv1['SAPT MP4 DISP'] = {'func': lambda x: x[0] * x[1] + x[2] + x[3] + x[4] + x[5],
                            'args': ['SAPT EXCHSCAL', 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP20 ENERGY',
                                     'SAPT DISP21 ENERGY', 'SAPT DISP22(SDQ) ENERGY', 'SAPT EST.DISP22(T) ENERGY']}
    pv1['SAPT CCD DISP'] = {'func': lambda x: x[0] * x[1] + x[2] + x[3] + x[4],
                            'args': ['SAPT EXCHSCAL', 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP2(CCD) ENERGY',
                                     'SAPT DISP22(S)(CCD) ENERGY', 'SAPT EST.DISP22(T)(CCD) ENERGY']}
    pv1['SAPT0 ELST ENERGY'] = {'func': sum, 'args': ['SAPT ELST10,R ENERGY']}
    pv1['SAPT0 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT EXCH10 ENERGY']}
    pv1['SAPT0 IND ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3],
                                'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ENERGY', 'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY']}
    pv1['SAPT0 DISP ENERGY'] = {'func': lambda x: x[0] * x[1] + x[2],
                                'args': ['SAPT EXCHSCAL', 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP20 ENERGY']}
    pv1['SAPT0 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT0 ELST ENERGY', 'SAPT0 EXCH ENERGY', 'SAPT0 IND ENERGY', 'SAPT0 DISP ENERGY']}
    pv1['SSAPT0 ELST ENERGY'] = {'func': sum, 'args': ['SAPT0 ELST ENERGY']}
    pv1['SSAPT0 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT0 EXCH ENERGY']}
    pv1['SSAPT0 IND ENERGY'] = {'func': lambda x: x[1] + (x[0] - 1.0) * x[2],
                                 'args': ['SAPT EXCHSCAL3', 'SAPT0 IND ENERGY', 'SAPT EXCH-IND20,R ENERGY']}
    pv1['SSAPT0 DISP ENERGY'] = {'func': lambda x: x[0] * x[1] + x[2],
                                 'args': ['SAPT EXCHSCAL3', 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP20 ENERGY']}
    pv1['SSAPT0 TOTAL ENERGY'] = {'func': sum, 'args': ['SSAPT0 ELST ENERGY', 'SSAPT0 EXCH ENERGY', 'SSAPT0 IND ENERGY', 'SSAPT0 DISP ENERGY']}
    pv1['SCS-SAPT0 ELST ENERGY'] = {'func': sum, 'args': ['SAPT0 ELST ENERGY']}
    pv1['SCS-SAPT0 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT0 EXCH ENERGY']}
    pv1['SCS-SAPT0 IND ENERGY'] = {'func': sum, 'args': ['SAPT0 IND ENERGY']}
    pv1['SCS-SAPT0 DISP ENERGY'] = {'func': lambda x: (x[0] - x[3]) * (x[1] + x[2]) + x[3] * (x[4] + x[5]),
                                    'args': [0.66, 'SAPT SAME-SPIN EXCH-DISP20 ENERGY', 'SAPT SAME-SPIN DISP20 ENERGY',
                                             1.2, 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP20 ENERGY']}  # note no xs for SCS disp
    pv1['SCS-SAPT0 TOTAL ENERGY'] = {'func': sum, 'args': ['SCS-SAPT0 ELST ENERGY', 'SCS-SAPT0 EXCH ENERGY', 'SCS-SAPT0 IND ENERGY', 'SCS-SAPT0 DISP ENERGY']}
    pv1['SAPT2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT ELST10,R ENERGY', 'SAPT ELST12,R ENERGY']}
    pv1['SAPT2 EXCH ENERGY'] = {'func': lambda x: x[1] + x[0] * (x[2] + x[3]),
                                'args': ['SAPT EXCHSCAL', 'SAPT EXCH10 ENERGY', 'SAPT EXCH11(S^2) ENERGY', 'SAPT EXCH12(S^2) ENERGY']}
    pv1['SAPT2 IND ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5],
                                'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ENERGY', 'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY',
                                         'SAPT IND22 ENERGY', 'SAPT EXCH-IND22 ENERGY']}
    pv1['SAPT2 DISP ENERGY'] = {'func': lambda x: x[0] * x[1] + x[2],
                                'args': ['SAPT EXCHSCAL', 'SAPT EXCH-DISP20 ENERGY', 'SAPT DISP20 ENERGY']}
    pv1['SAPT2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2 ELST ENERGY', 'SAPT2 EXCH ENERGY', 'SAPT2 IND ENERGY', 'SAPT2 DISP ENERGY']}
    pv1['SAPT2+ ELST ENERGY'] = {'func': sum, 'args': ['SAPT ELST10,R ENERGY', 'SAPT ELST12,R ENERGY']}
    pv1['SAPT2+ EXCH ENERGY'] = {'func': lambda x: x[1] + x[0] * (x[2] + x[3]),
                                 'args': ['SAPT EXCHSCAL', 'SAPT EXCH10 ENERGY', 'SAPT EXCH11(S^2) ENERGY', 'SAPT EXCH12(S^2) ENERGY']}
    pv1['SAPT2+ IND ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5],
                                 'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ENERGY', 'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY',
                                          'SAPT IND22 ENERGY', 'SAPT EXCH-IND22 ENERGY']}
    pv1['SAPT2+ DISP ENERGY'] = {'func': sum, 'args': ['SAPT MP4 DISP']}
    pv1['SAPT2+ TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+ ELST ENERGY', 'SAPT2+ EXCH ENERGY', 'SAPT2+ IND ENERGY', 'SAPT2+ DISP ENERGY']}
    pv1['SAPT2+(CCD) ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+ ELST ENERGY']}
    pv1['SAPT2+(CCD) EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+ EXCH ENERGY']}
    pv1['SAPT2+(CCD) IND ENERGY'] = {'func': sum, 'args': ['SAPT2+ IND ENERGY']}
    pv1['SAPT2+(CCD) DISP ENERGY'] = {'func': sum, 'args': ['SAPT CCD DISP']}
    pv1['SAPT2+(CCD) TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(CCD) ELST ENERGY', 'SAPT2+(CCD) EXCH ENERGY', 'SAPT2+(CCD) IND ENERGY', 'SAPT2+(CCD) DISP ENERGY']}
    pv1['SAPT2+DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+ ELST ENERGY']}
    pv1['SAPT2+DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+ EXCH ENERGY']}
    pv1['SAPT2+DMP2 IND ENERGY'] = {'func': sum, 'args': ['SAPT2+ IND ENERGY', 'SAPT MP2(2) ENERGY']}
    pv1['SAPT2+DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+ DISP ENERGY']}
    pv1['SAPT2+DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+DMP2 ELST ENERGY', 'SAPT2+DMP2 EXCH ENERGY', 'SAPT2+DMP2 IND ENERGY', 'SAPT2+DMP2 DISP ENERGY']}
    pv1['SAPT2+(CCD)DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+ ELST ENERGY']}
    pv1['SAPT2+(CCD)DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+ EXCH ENERGY']}
    pv1['SAPT2+(CCD)DMP2 IND ENERGY'] = {'func': sum, 'args': ['SAPT2+DMP2 IND ENERGY']}
    pv1['SAPT2+(CCD)DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+(CCD) DISP ENERGY']}
    pv1['SAPT2+(CCD)DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(CCD)DMP2 ELST ENERGY', 'SAPT2+(CCD)DMP2 EXCH ENERGY', 'SAPT2+(CCD)DMP2 IND ENERGY', 'SAPT2+(CCD)DMP2 DISP ENERGY']}
    pv1['SAPT2+(3) ELST ENERGY'] = {'func': sum, 'args': ['SAPT ELST10,R ENERGY', 'SAPT ELST12,R ENERGY', 'SAPT ELST13,R ENERGY']}
    pv1['SAPT2+(3) EXCH ENERGY'] = {'func': lambda x: x[1] + x[0] * (x[2] + x[3]),
                                    'args': ['SAPT EXCHSCAL', 'SAPT EXCH10 ENERGY', 'SAPT EXCH11(S^2) ENERGY', 'SAPT EXCH12(S^2) ENERGY']}
    pv1['SAPT2+(3) IND ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5],
                                    'args': ['SAPT EXCHSCAL', 'SAPT HF(2) ENERGY', 'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY',
                                             'SAPT IND22 ENERGY', 'SAPT EXCH-IND22 ENERGY']}
    pv1['SAPT2+(3) DISP ENERGY'] = {'func': sum, 'args': ['SAPT MP4 DISP', 'SAPT DISP30 ENERGY']}
    pv1['SAPT2+(3) TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) ELST ENERGY', 'SAPT2+(3) EXCH ENERGY', 'SAPT2+(3) IND ENERGY', 'SAPT2+(3) DISP ENERGY']}
    pv1['SAPT2+(3)(CCD) ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) ELST ENERGY']}
    pv1['SAPT2+(3)(CCD) EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) EXCH ENERGY']}
    pv1['SAPT2+(3)(CCD) IND ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) IND ENERGY']}
    pv1['SAPT2+(3)(CCD) DISP ENERGY'] = {'func': sum, 'args': ['SAPT CCD DISP', 'SAPT DISP30 ENERGY']}
    pv1['SAPT2+(3)(CCD) TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(3)(CCD) ELST ENERGY', 'SAPT2+(3)(CCD) EXCH ENERGY', 'SAPT2+(3)(CCD) IND ENERGY', 'SAPT2+(3)(CCD) DISP ENERGY']}
    pv1['SAPT2+(3)DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) ELST ENERGY']}
    pv1['SAPT2+(3)DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) EXCH ENERGY']}
    pv1['SAPT2+(3)DMP2 IND ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) IND ENERGY', 'SAPT MP2(2) ENERGY']}
    pv1['SAPT2+(3)DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) DISP ENERGY']}
    pv1['SAPT2+(3)DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(3)DMP2 ELST ENERGY', 'SAPT2+(3)DMP2 EXCH ENERGY', 'SAPT2+(3)DMP2 IND ENERGY', 'SAPT2+(3)DMP2 DISP ENERGY']}
    pv1['SAPT2+(3)(CCD)DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) ELST ENERGY']}
    pv1['SAPT2+(3)(CCD)DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+(3) EXCH ENERGY']}
    pv1['SAPT2+(3)(CCD)DMP2 IND ENERGY'] = {'func': sum, 'args': ['SAPT2+(3)DMP2 IND ENERGY']}
    pv1['SAPT2+(3)(CCD)DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+(3)(CCD) DISP ENERGY']}
    pv1['SAPT2+(3)(CCD)DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+(3)(CCD)DMP2 ELST ENERGY', 'SAPT2+(3)(CCD)DMP2 EXCH ENERGY', 'SAPT2+(3)(CCD)DMP2 IND ENERGY', 'SAPT2+(3)(CCD)DMP2 DISP ENERGY']}
    pv1['SAPT2+3 ELST ENERGY'] = {'func': sum, 'args': ['SAPT ELST10,R ENERGY', 'SAPT ELST12,R ENERGY', 'SAPT ELST13,R ENERGY']}
    pv1['SAPT2+3 EXCH ENERGY'] = {'func': lambda x: x[1] + x[0] * (x[2] + x[3]),
                                  'args': ['SAPT EXCHSCAL', 'SAPT EXCH10 ENERGY', 'SAPT EXCH11(S^2) ENERGY', 'SAPT EXCH12(S^2) ENERGY']}
    pv1['SAPT2+3 IND ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5] + x[6] + x[0] * x[7],
                                  'args': ['SAPT EXCHSCAL', 'SAPT HF(3) ENERGY', 'SAPT IND20,R ENERGY', 'SAPT EXCH-IND20,R ENERGY',
                                           'SAPT IND22 ENERGY', 'SAPT EXCH-IND22 ENERGY', 'SAPT IND30,R ENERGY', 'SAPT EXCH-IND30,R ENERGY']}
    pv1['SAPT2+3 DISP ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5],
                                  'args': ['SAPT EXCHSCAL', 'SAPT MP4 DISP', 'SAPT DISP30 ENERGY', 'SAPT EXCH-DISP30 ENERGY',
                                           'SAPT IND-DISP30 ENERGY', 'SAPT EXCH-IND-DISP30 ENERGY']}
    pv1['SAPT2+3 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+3 ELST ENERGY', 'SAPT2+3 EXCH ENERGY', 'SAPT2+3 IND ENERGY', 'SAPT2+3 DISP ENERGY']}
    pv1['SAPT2+3(CCD) ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+3 ELST ENERGY']}
    pv1['SAPT2+3(CCD) EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+3 EXCH ENERGY']}
    pv1['SAPT2+3(CCD) IND ENERGY'] = {'func': sum, 'args': ['SAPT2+3 IND ENERGY']}
    pv1['SAPT2+3(CCD) DISP ENERGY'] = {'func': lambda x: x[1] + x[2] + x[0] * x[3] + x[4] + x[0] * x[5],
                                       'args': ['SAPT EXCHSCAL', 'SAPT CCD DISP', 'SAPT DISP30 ENERGY', 'SAPT EXCH-DISP30 ENERGY',
                                                'SAPT IND-DISP30 ENERGY', 'SAPT EXCH-IND-DISP30 ENERGY']}
    pv1['SAPT2+3(CCD) TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+3(CCD) ELST ENERGY', 'SAPT2+3(CCD) EXCH ENERGY', 'SAPT2+3(CCD) IND ENERGY', 'SAPT2+3(CCD) DISP ENERGY']}
    pv1['SAPT2+3DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+3 ELST ENERGY']}
    pv1['SAPT2+3DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+3 EXCH ENERGY']}
    pv1['SAPT2+3DMP2 IND ENERGY'] = {'func': sum, 'args': ['SAPT2+3 IND ENERGY', 'SAPT MP2(3) ENERGY']}
    pv1['SAPT2+3DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+3 DISP ENERGY']}
    pv1['SAPT2+3DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+3DMP2 ELST ENERGY', 'SAPT2+3DMP2 EXCH ENERGY', 'SAPT2+3DMP2 IND ENERGY', 'SAPT2+3DMP2 DISP ENERGY']}
    pv1['SAPT2+3(CCD)DMP2 ELST ENERGY'] = {'func': sum, 'args': ['SAPT2+3 ELST ENERGY']}
    pv1['SAPT2+3(CCD)DMP2 EXCH ENERGY'] = {'func': sum, 'args': ['SAPT2+3 EXCH ENERGY']}
    pv1['SAPT2+3(CCD)DMP2 IND ENERGY'] = {'func': sum, 'args': ['SAPT2+3DMP2 IND ENERGY']}
    pv1['SAPT2+3(CCD)DMP2 DISP ENERGY'] = {'func': sum, 'args': ['SAPT2+3(CCD) DISP ENERGY']}
    pv1['SAPT2+3(CCD)DMP2 TOTAL ENERGY'] = {'func': sum, 'args': ['SAPT2+3(CCD)DMP2 ELST ENERGY', 'SAPT2+3(CCD)DMP2 EXCH ENERGY', 'SAPT2+3(CCD)DMP2 IND ENERGY', 'SAPT2+3(CCD)DMP2 DISP ENERGY']}

    return pv1
