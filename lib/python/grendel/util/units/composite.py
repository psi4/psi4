from __future__ import division
from collections import defaultdict
from copy import copy
from fractions import Fraction
from itertools import groupby
from numbers import Number, Real
from grendel import type_checking_enabled
import math
from grendel.util.decorators import typechecked
from grendel.util.misc import is_near_integer, is_pretty_fraction
from grendel.util.strings import classname

__all__ = [
    "CompositeUnit"
]

class CompositeUnit(object):
    """ A class for composite units, e.g. Joules/Angstrom^2
    """

    ##############
    # Attributes #
    ##############

    composed = None
    """ dict of Unit -> exponent pairs """

    coefficient = 1.0


    ##################
    # Initialization #
    ##################

    def __init__(self,
            numer_units_or_composed,
            denom_units=None,
            coeff=1.0,
        ):
        self.composed = defaultdict(lambda: 0)
        if isinstance(numer_units_or_composed, dict):
            self.composed.update(numer_units_or_composed)
        else:
            for unit in numer_units_or_composed:
                self.composed[unit] += 1
            for unit in denom_units:
                self.composed[unit] -= 1
        self.coefficient = coeff

    ##############
    # Properties #
    ##############

    @property
    def genre(self):
        """ Analog of `Unit.genre`.  Always returns `CompositeUnit`
        """
        return CompositeUnit

    @property
    def name(self):
        """ The name of the composite unit, as a product of base units

        """
        numsorted = sorted(
            [(unit, exponent) for unit, exponent in self.composed.iteritems() if exponent > 0],
            key=lambda x: (x[0].genre.__name__, x[0].__name__)
        )
        numstrs = []
        for unit, exp in numsorted:
            if exp == 1:
                numstrs.append(unit.name)
            elif is_near_integer(exp):
                numstrs.append(unit.name + "**" + str(int(round(exp))))
            elif is_pretty_fraction(exp, largest_denominator=9):
                numstrs.append(unit.name + "**(" + str(Fraction(exp).limit_denominator(9)) + ")")
            else:
                numstrs.append('{}**{:.5f}'.format(unit.name, exp))
        #----------------------------------------#
        denomsorted = sorted(
            [(unit, exponent) for unit, exponent in self.composed.iteritems() if exponent < 0],
            key=lambda x: (x[0].genre.__name__, x[0].__name__)
        )
        denomstrs = []
        for unit, exp in denomsorted:
            if exp == -1:
                denomstrs.append(unit.name)
            elif is_near_integer(exp):
                denomstrs.append(unit.name + "**" + str(int(round(abs(exp)))))
            elif is_pretty_fraction(exp, largest_denominator=9):
                denomstrs.append(unit.name + "**(" + str(Fraction(abs(exp)).limit_denominator(9)) + ")")
            else:
                denomstrs.append('{}**{:.5f}'.format(unit.name, abs(exp)))
        #----------------------------------------#
        ret_val = " * ".join(numstrs)
        if len(numstrs) == 0:
            ret_val += str(self.coefficient)
        if len(denomstrs) > 1:
            ret_val += " / (" + " * ".join(denomstrs) + ")"
        elif len(denomstrs) == 1:
            ret_val += " / " + denomstrs[0]
        if self.coefficient != 1.0 and not len(numsorted) == 0:
            ret_val = "(" + str(self.coefficient) + " " + ret_val + ")"
        return ret_val


    ###################
    # Special Methods #
    ###################

    def __contains__(self, item):
        """
        Mimics the behavior of `__contains__` in `Unit`
        """
        if isinstance(item, ValueWithUnits):
            if item.units == self:
                return True
            else:
                return False
        else:
            raise TypeError()

    #----------------------#
    # Comparison Operators #
    #----------------------#

    def __eq__(self, other):
        try:
            return self.to(other) == 1.0
        except IncompatibleUnitsError:
            return False

    #----------------------#
    # Arithmetic Operators #
    #----------------------#

    def __mul__(self, other):
        new_composed = copy(self.composed)
        if isinstance(other, CompositeUnit):
            for unit, exp in other.composed.iteritems():
                new_composed[unit] += exp
            new_coeff = self.coefficient * other.coefficient
            return CompositeUnit(new_composed, coeff=new_coeff)
        elif isinstance(other, Unit):
            new_composed[other] += 1
            return CompositeUnit(new_composed, coeff=self.coefficient)
        elif isinstance(other, Number):
            return ValueWithUnits(self.coefficient * other, CompositeUnit(new_composed))
        else: # pragma: no cover
            return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, Number):
            return ValueWithUnits(self.coefficient * other, CompositeUnit(copy(self.composed)))
        if isinstance(other, Unit):
            new_composed = copy(self.composed)
            new_composed[other] += 1
            return CompositeUnit(new_composed, coeff=self.coefficient)
        else: # pragma: no cover
            return NotImplemented

    def __div__(self, other):
        new_composed = copy(self.composed)
        if isinstance(other, CompositeUnit):
            for unit, exp in other.composed.iteritems():
                new_composed[unit] -= exp
            new_coeff = self.coefficient / other.coefficient
            return CompositeUnit(new_composed, coeff=new_coeff)
        elif isinstance(other, Unit):
            new_composed[other] -= 1
            return CompositeUnit(new_composed, coeff=self.coefficient)
        elif isinstance(other, Number):
            return ValueWithUnits(self.coefficient / other, CompositeUnit(new_composed))
        else: # pragma: no cover
            return NotImplemented
    __truediv__ = __div__

    def __rdiv__(self, other):
        new_composed = defaultdict(lambda: 0)
        for unit, exp in self.composed.iteritems():
            new_composed[unit] = -exp
        if isinstance(other, Unit):
            new_composed[other] += 1
            return CompositeUnit(new_composed, coeff=1.0/self.coefficient)
        elif isinstance(other, Number):
            return ValueWithUnits(other / self.coefficient, CompositeUnit(new_composed))
        else: # pragma: no cover
            return NotImplemented
    __rtruediv__ = __rdiv__

    def __pow__(self, power):
        if isinstance(power, Real):
            new_composed = copy(self.composed)
            for unit, exp in self.composed.iteritems():
                new_composed[unit] = exp * power
            new_coeff = self.coefficient ** power
            return CompositeUnit(new_composed, coeff=new_coeff)
        else: # pragma: no cover
            return NotImplemented

    #------------------------#
    # Output Representations #
    #------------------------#

    def __repr__(self):
        return self.name
    __str__ = __repr__


    ###########
    # Methods #
    ###########

    def reduced(self, using_units=None):
        """ Reduce the composite unit into the fewest base units possible.

        Examples
        --------

        >>> from grendel.util.units import *
        >>> (Angstroms**2/Bohr).reduced()
        (1.88972612457 Angstrom)
        >>> (Degree / (Bohr**2)).reduced()
        Degree / Bohr**2
        >>> (Degree / (Bohr*Meter)).reduced()
        (5.2917721092e-11 Degree / Bohr**2)

        """
        using_units = using_units or {}
        genre_groups = defaultdict(lambda: [])
        for unit, exp in self.composed.iteritems():
            genre_groups[unit.genre].append((unit, exp))
        new_coeff = self.coefficient
        new_composed = defaultdict(lambda: 0)
        for genre, group in genre_groups.iteritems():
            # Decide which unit to reduce each genre to
            if genre in using_units:
                to_unit = using_units[genre]
            elif any(pair[0] is genre.default for pair in group):
                to_unit = genre.default
            elif any(not isinstance(pair[0], PrefixedUnit) for pair in group):
                # arbitrarily choose the first non-prefixed unit alphabetically
                to_unit = min(
                    [pair[0] for pair in group if not isinstance(pair[0], PrefixedUnit)],
                    key=lambda x: x.name
                )
            else:
                # if they are all prefixed units, choose the one that has the smallest (absolute)
                #   prefix exponent, then the non-prefixed unit corresponding to the genre default,
                #   and finally, if all else fails, the first alphabetically
                to_unit = min(
                    [pair[0] for pair in group],
                    key=lambda x: (abs(math.log10(x.prefix.multiplier)),
                                   not x.base_unit is genre.default,
                                   x.name)
                )
            for unit, exp in group:
                new_composed[to_unit] += exp
                new_coeff *= unit.to(to_unit) ** exp
        # clean out 0-exponent keys
        for unit in [k for k in new_composed.keys()]:
            if new_composed[unit] == 0:
                del new_composed[unit]
        return CompositeUnit(new_composed, coeff=new_coeff)

    @typechecked(other_in='isunit')
    def to(self, other_in):
        self_red = self.reduced(UnitGenre.GenreDefaultDict())
        #----------------------------------------#
        if not isinstance(other_in, CompositeUnit):
            if len(self_red.composed) == 1 and self_red.composed.values()[0] == 1.0:
                return self_red.coefficient/other_in.to(self_red.composed.keys()[0])
            else:
                raise IncompatibleUnitsError(self, other_in)
        #----------------------------------------#
        other = other_in.reduced(UnitGenre.GenreDefaultDict())
        #----------------------------------------#
        if len(self_red.composed) != len(other.composed):
            raise IncompatibleUnitsError(self, other_in)
        #----------------------------------------#
        conv_factor = self_red.coefficient / other.coefficient
        for (my_unit, my_exp), (o_unit, o_exp) in zip(
                sorted(self_red.composed.iteritems(), key=lambda x: (x[0].genre, x[0].__name__)),
                sorted(other.composed.iteritems(), key=lambda x: (x[0].genre, x[0].__name__))):
            if my_exp != o_exp:
                raise IncompatibleUnitsError(self, other_in)
            conv_factor *= my_unit.to(o_unit) ** my_exp
        #----------------------------------------#
        return conv_factor


#####################
# Dependent Imports #
#####################

from grendel.util.units.errors import IncompatibleUnitsError
from grendel.util.units.unit import Unit
from grendel.util.units.value_with_units import ValueWithUnits
from grendel.util.units.unit import isunit, PrefixedUnit
from grendel.util.units.unit import UnitGenre

