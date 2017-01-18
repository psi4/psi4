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

from __future__ import absolute_import
from __future__ import print_function
from decimal import Decimal, ROUND_FLOOR, ROUND_CEILING
from .exceptions import *


class PreservingDict(dict):
    """Class to store quantum chemical quantities extracted from output
    files. Extends the dictionary object to (1) store key as all-caps
    version of itself and (2) validate value for duplicate values for the
    same key by testing which has more decimal places and whether value
    the same within a plausing rounding error. Allows consistency checks
    when parsing output files without loss of precision.

    """

    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def __setitem__(self, key, value):
        try:
            key = key.upper()
        except AttributeError:
            raise AttributeError('Keys stored as upper-case strings: %s unsuitable' % (key))
        value = Decimal(value)
        if key in self.keys() and not 'CURRENT' in key:
            # Validate choosing more detailed value for variable
            existing_exp = self[key].as_tuple().exponent  # 0.1111 --> -4
            candidate_exp = value.as_tuple().exponent
            if existing_exp > candidate_exp:  # candidate has more digits
                places = Decimal(10) ** (existing_exp + 1)  # exp+1 permits slack in rounding
                best_value = value
            else:                             # existing has more digits
                places = Decimal(10) ** (candidate_exp + 1)
                best_value = self[key]
            # Validate values are the same
            places = max(places, Decimal('1E-11'))  # for computed psivars
            #print('FLOOR: ', self[key].quantize(places, rounding=ROUND_FLOOR) - value.quantize(places, rounding=ROUND_FLOOR))
            #print('CEIL:  ', self[key].quantize(places, rounding=ROUND_CEILING) - value.quantize(places, rounding=ROUND_CEILING))
            if (self[key].quantize(places, rounding=ROUND_CEILING).compare(value.quantize(places, rounding=ROUND_CEILING)) != 0) and \
               (self[key].quantize(places, rounding=ROUND_FLOOR).compare(value.quantize(places, rounding=ROUND_FLOOR)) != 0):
                raise ParsingValidationError(
                    """Output file yielded both %s and %s as values for quantity %s.""" %
                    (self[key].to_eng_string(), value.to_eng_string(), key))
            #print 'Resetting variable %s to %s' % (key, best_value.to_eng_string())
        else:
            best_value = value
            #print 'Setting   variable %s to %s' % (key, best_value.to_eng_string())
        super(PreservingDict, self).__setitem__(key, best_value)

    def update(self, *args, **kwargs):
        if args:
            if len(args) > 1:
                raise TypeError("update expected at most 1 arguments, "
                                "got %d" % len(args))
            other = dict(args[0])
            for key in other:
                self[key] = other[key]
        for key in kwargs:
            self[key] = kwargs[key]

    def setdefault(self, key, value=None):
        if key not in self:
            self[key] = value
        return self[key]


if __name__ == '__main__':
    c4info = PreservingDict()
    c4info['scf 4.5e0 total energy'] = '-1.e-4'
    c4info['1.3'] = '.4'
    c4info['curl'] = '-437.12345678'
    c4info['curl'] = '-437.12345677'
    c4info['curl'] = '-437.123456'
    c4info['curl'] = '-437.123457'
    c4info['curl'] = '-437.1234444'  # fails
    c4info['curl'] = '-437.123456789'
    #c4info['curl'] = '-437.1234567779'  # fails
    print(c4info)
