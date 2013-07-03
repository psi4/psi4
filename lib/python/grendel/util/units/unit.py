from __future__ import absolute_import
from collections import defaultdict
import math
from numbers import Number, Real
from grendel.util.aliasing import function_alias
from grendel.util.strings import classname
import sys
# Version 3 compatibility
if sys.version_info[0] == 3:
    basestring = str


__all__ = [
    'Unit',
    'isunit', 'is_unit',
    'convert', 'convert_units',
    'compatible_units', 'iscompatible',

    # Unit genres:
    'DistanceUnit',
    'EnergyUnit',
    'AngularUnit',
    'ElectricChargeUnit',
    'MassUnit',
    'TimeUnit'
]

#############
# Utilities #
#############

def isunit(unit):
    if isinstance(unit, Unit) or isinstance(unit, CompositeUnit):
        return True
    else:
        return False
is_unit = function_alias('is_unit', isunit)

def plural(unit): # pragma: no cover
    if not isunit(unit):
        raise TypeError
    return unit.__plural__

def convert_units(val, from_unit, to_unit):
    if not isunit(from_unit):
        raise UnknownUnitError(from_unit)
    if not isunit(to_unit):
        raise UnknownUnitError(to_unit)
    if from_unit == to_unit:
        return val
    return val * from_unit.to(to_unit)
convert = function_alias('convert', convert_units)

def compatible_units(unit1, unit2):
    try:
        unit1.to(unit2)
        return True
    except IncompatibleUnitsError:
        return False
iscompatible = function_alias('iscompatible', compatible_units)

########################
# Metaclasses and such #
########################

class Prefix(object):
    """ The prefix for a unit, e.g. centi-, kilo-, mega-, etc.
    """

    ##############
    # Attributes #
    ##############

    in_words = None
    abbrev = None
    multiplier = None

    ##################
    # Initialization #
    ##################

    def __init__(self, in_words, abbrev, multiplier):
        self.in_words = in_words
        self.abbrev = abbrev
        self.multiplier = multiplier


class Unit(type):
    """ Metaclass for a general unit of something.
    """

    ########################
    # Metaclass Attributes #
    ########################

    known_units = []
    prefixes = [
        Prefix("Yotta", "Y", 1.0e24),
        Prefix("Zetta", "Z", 1.0e21),
        Prefix("Exa", "E", 1.0e18),
        Prefix("Peta", "P", 1.0e15),
        Prefix("Tera", "T", 1.0e12),
        Prefix("Giga", "G", 1.0e9),
        Prefix("Mega", "M", 1.0e6),
        Prefix("Kilo", "k", 1.0e3),
        Prefix("Hecto", "h", 100.0),
        Prefix("Deca", "da", 10.0),
        Prefix("Deci", "d", 1.0e-1),
        Prefix("Centi", "c", 1.0e-2),
        Prefix("Milli", "m", 1.0e-3),
        Prefix("Micro", "u", 1.0e-6),
        Prefix("Nano", "n", 1.0e-9),
        Prefix("Pico", "p", 1.0e-12),
        Prefix("Femto", "f", 1.0e-15),
        Prefix("Atto", "a", 1.0e-18),
        Prefix("Zepto", "z", 1.0e-21),
        Prefix("Yocto", "y", 1.0e-24)
    ]


    ####################
    # Class Attributes #
    ####################

    __plural__ = None
    __aliases__ = None
    __abbrev__ = None
    __prefixed__ = True


    ############################
    # Metaclass Initialization #
    ############################

    def __init__(cls, name, bases, dct):
        Unit.known_units.append(name)
        if not all(issubclass(base, UnitGenre) for base in bases) or not len(bases) == 1:
            raise TypeError("Units must inherit from a single class with the UnitGenre superclass.")
        super(Unit, cls).__init__(name, bases, dct)
        globals()['__all__'].append(str(cls))
        # Automatically create a plural alias for the unit if one is not given
        if cls.__plural__ is None:
            cls.__plural__ = str(cls) + "s"
        if not cls.__plural__ == name:
            globals()[cls.__plural__] = cls
            globals()['__all__'].append(cls.__plural__)
            Unit.known_units.append(cls.__plural__)
            # Automatically create prefixed versions of the unit
        if cls.__prefixed__:
            for prefix in Unit.prefixes:
                d = {'prefix': prefix, 'base_unit': cls}
                name1 = prefix.in_words + name
                pre = PrefixedUnit.__new__(PrefixedUnit, name1, (cls,), d)
                globals()[name1] = pre
                globals()['__all__'].append(name1)
                Unit.known_units.append(pre)
                name2 = prefix.in_words + cls.__plural__
                globals()[name2] = pre
                globals()['__all__'].append(name1)
                # If the name is not CamelCase or UPPERCASE, append uncapitalized versions (e.g. Kilogram as well
                # as KiloGram, but not KiloaMU, only KiloAMU)
                if not any(letter.isupper() for letter in name[1:]):
                    name3 = prefix.in_words + name[0].lower() + name[1:]
                    globals()[name3] = pre
                    globals()['__all__'].append(name3)
                    name4 = prefix.in_words + cls.__plural__[0].lower() + cls.__plural__[1:]
                    globals()[name4] = pre
                    globals()['__all__'].append(name4)

    ####################
    # Class Properties #
    ####################

    @property
    def genre(cls):
        return cls.__mro__[1]

    @property
    def name(cls):
        return cls.__name__

    #########################
    # Special Class Methods #
    #########################

    def __contains__(self, item):
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

    def __eq__(cls, other):
        try:
            return other.to(cls) == 1.0
        except IncompatibleUnitsError:
            return False
        except AttributeError:
            # Other doesn't even have a 'to()' method...
            return NotImplemented

    def __ne__(self, other):
        eq_val = self.__eq__(other)
        if eq_val is NotImplemented:
            return NotImplemented
        else:
            return not eq_val

    #----------------------#
    # Arithmetic Operators #
    #----------------------#

    def __mul__(cls, other):
        if isinstance(other, Number):
            return ValueWithUnits(other, cls)
        elif isinstance(other, Unit):
            return CompositeUnit({cls: 1, other: 1})
        else:
            return NotImplemented

    def __rmul__(cls, other):
        if isinstance(other, Number):
            return ValueWithUnits(other, cls)
        else:
            return NotImplemented

    def __div__(cls, other):
        if isinstance(other, Unit):
            return CompositeUnit({cls: 1, other:-1})
        else:
            return NotImplemented
    __truediv__ = __div__

    def __rdiv__(cls, other):
        if isinstance(other, Number):
            return ValueWithUnits(other, CompositeUnit({cls: -1}))
        else: # pragma: no cover
            return NotImplemented
    __rtruediv__ = __rdiv__

    def __pow__(cls, power):
        if isinstance(power, Real):
            return CompositeUnit({cls: power})
        else: # pragma: no cover
            return NotImplemented

    #------------------------#
    # Output Representations #
    #------------------------#

    def __repr__(cls):
        return classname(super(Unit, cls).__repr__())
    __str__ = __repr__


    #################
    # Class Methods #
    #################

    def genre_check(cls, other):
        if not issubclass(other, cls.genre):
            raise IncompatibleUnitsError(cls, other)

    def prefix_factor(cls, other):
        other_fact = 1.0
        if isinstance(other, PrefixedUnit):
            other_fact = other.prefix.multiplier
        my_fact = 1.0
        if isinstance(cls, PrefixedUnit):
            my_fact = cls.prefix.multiplier
        return my_fact / other_fact

    def to(cls, other):
        cls.genre_check(other)
        if other is cls:
            return 1.0
        elif issubclass(other, cls) or issubclass(cls, other):
            return cls.prefix_factor(other)
        else:
            return (1.0 / cls.genre.reference_unit.to(cls)) * cls.genre.reference_unit.to(other)


class PrefixedUnit(Unit):
    """ Metaclass for a unit with a prefix, such as a Kilogram, Centimeter, etc.
    """

    ####################
    # Class Attributes #
    ####################

    base_unit = None
    prefix = None

    ############################
    # Metaclass Initialization #
    ############################

    def __init__(cls, name, bases, dct):
        cls.known_units.append(name)
        if not 'to' in dct:
            dct['to'] = PrefixedUnit.to
        if not all(isinstance(base, Unit) for base in bases) or not len(bases) == 1: # pragma: no cover
            raise TypeError("PrefixedUnits must inherit from a single class which is a Unit.")
        super(Unit, cls).__init__(name, bases, dct)

    @property
    def genre(cls):
        return cls.base_unit.genre


class UnitGenre(object):
    """ Superclass for classes of things that can be measured by units.
    For instance, DistanceUnit, AngularUnit, EnergyUnit, etc.
    """
    default = None
    reference_unit = None

    class GenreDefaultDict(defaultdict):
        def __missing__(self, key):
            return key.genre


#--------------------------------------------------------------------------------#


####################
# Helper functions #
####################

# Can be called as either def_unit_alias(alias, unit) or def_unit_alias(unit, alias) (as long as alias is a str and
# is_unit(unit) is True)
def def_unit_alias(arg1, arg2, plural=True, prefixed=True): # pragma: no cover
    alias = None
    unit = None
    if isinstance(arg1, basestring) and is_unit(arg2):
        alias = arg1
        unit = arg2
    elif isinstance(arg2, basestring) and is_unit(arg1):
        alias = arg2
        unit = arg1
    else:
        raise TypeError()
    globals()[alias] = unit
    globals()['__all__'].append(alias)
    my_plural = None
    if plural is True:
        # Automatically add plural with 's' unless the user specifies a specific plural or if the user specifies 'False'
        globals()[alias + 's'] = unit
        globals()['__all__'].append(alias + 's')
        my_plural = alias + 's'
    elif plural is not False and not str(plural) == alias:
        my_plural = str(plural)
        globals()[my_plural] = unit
        globals()['__all__'].append(my_plural)
        # Automatically create prefixed versions of the unit alias
    if prefixed:
        for prefix in Unit.prefixes:
            d = {'prefix': prefix, 'base_unit': unit}
            name = prefix.in_words + alias
            pre = PrefixedUnit.__new__(PrefixedUnit, name, (unit,), d)
            PrefixedUnit.__init__(pre, name, (unit,), d)
            globals()[name] = pre
            globals()['__all__'].append(name)
            Unit.known_units.append(pre)
            if not plural is False:
                name = prefix.in_words + my_plural
                globals()[name] = pre
                globals()['__all__'].append(name)
            if not any(letter.isupper() for letter in alias[1:]):
                # If the name is not CamelCase or UPPERCASE, append uncapitalized versions
                #   (e.g. Kilogram as well as KiloGram, but not KiloaMU, only KiloAMU)
                name = prefix.in_words + alias[0].lower() + alias[1:]
                globals()[name] = pre
                globals()['__all__'].append(name)
                if not plural is False:
                    name = prefix.in_words + my_plural[0].lower() + my_plural[1:]
                    globals()[name] = pre
                    globals()['__all__'].append(name)

def def_unit_aliases(unit, *args, **kwargs): # pragma: no cover
    for al in args:
        alias = str(al)
        plural = kwargs.pop(al + "_plural", True)
        prefixed = kwargs.pop(al + "_prefixed", True)
        def_unit_alias(unit, alias, plural, prefixed)

def def_unit(genre, unit, plural=True, prefixed=True):
    d = {} #{'to': Unit.to}
    if plural is False: # pragma: no cover
        # Define a plural that is the same as the unit to prevent plural from being defined
        d['__plural__'] = unit
    elif plural is not True:
        # When plural is True, use the default plural.  Otherwise, define it
        d['__plural__'] = str(plural)
    new_cls = globals()[unit] = Unit.__new__(Unit, unit, (genre,), d)
    new_cls.__prefixed__ = prefixed
    Unit.__init__(globals()[unit], unit, (genre,), d)
    globals()['__all__'].append(unit)

def def_units(genre, *args, **kwargs): # pragma: no cover
    for unit in args:
        prefixed = kwargs.pop(unit + "_prefixed", True)
        plural = kwargs.pop(unit + "_plural", True)
        def_unit(genre, unit, plural, prefixed)
        if (unit + "_alias") in kwargs:
            if (unit + "_alias_plural") in kwargs:
                def_unit_alias(kwargs[unit + "_alias"], eval(unit, globals()), kwargs[unit + "_alias_plural"])
            elif kwargs[unit + 'alias'] + "_alias" in kwargs:
                def_unit_alias(kwargs[unit + "_alias"], eval(unit, globals()), kwargs[kwargs[unit + "_alias"] + '_alias'])
            else:
                def_unit_alias(kwargs[unit + "_alias"], eval(unit, globals()))
        elif (unit + "_aliases") in kwargs:
            for alias in kwargs[unit + "_aliases"]:
                aplural = kwargs.pop(alias + "_plural", True)
                aprefixed = kwargs.pop(alias + "_prefixed", prefixed)
                def_unit_alias(alias, eval(unit, globals()), aplural, aprefixed)


#--------------------------------------------------------------------------------#


##################
# Distance Units #
##################

class DistanceUnit(UnitGenre):
    """ General superclass for all distance units
    """

class Angstrom(DistanceUnit):
    __metaclass__ = Unit
    @classmethod
    def to(cls, other):
        cls.genre_check(other)
        pf = cls.prefix_factor(other)
        if issubclass(other, Bohr):
            return pf / BohrRadius.value
        elif issubclass(other, Meter):
            return 1e-10 * pf
        elif issubclass(other, Angstrom):
            return pf
        else: # pragma: no cover
            raise NotImplementedError("Conversion from units " + classname(cls) + " to units " + classname(other) + " is not implemented.")
DistanceUnit.reference_unit = Angstrom

class Bohr(DistanceUnit):
    __metaclass__ = Unit
def_unit_alias('AtomicUnitOfDistance', Bohr, plural='AtomicUnitsOfDistance')

class Meter(DistanceUnit):
    __metaclass__ = Unit

DistanceUnit.default = Angstrom
#DistanceUnit.default = Bohr

#################
# Angular Units #
#################

class AngularUnit(UnitGenre):
    """ General superclass for all angular units
    """

class Degree(AngularUnit):
    __metaclass__ = Unit
    @classmethod
    def to(cls, other):
        cls.genre_check(other)
        pf = cls.prefix_factor(other)
        if issubclass(other, Radian):
            return pf * math.pi / 180.0
        elif issubclass(other, Degree):
            return pf
        else: # pragma: no cover
            raise NotImplementedError("Conversion from units " + classname(cls) + " to units " + classname(other) + " is not implemented.")
AngularUnit.reference_unit = Degree

class Radian(AngularUnit):
    __metaclass__ = Unit

# For now, using default units of Radians causes some unit tests to fail
#AngularUnit.default = Radian
AngularUnit.default = Degree

################
# Energy Units #
################

class EnergyUnit(UnitGenre):
    """ General superclass for all energy units
    """

class Joule(EnergyUnit):
    __metaclass__ = Unit
    @classmethod
    def to(cls, other):
        cls.genre_check(other)
        pf = cls.prefix_factor(other)
        if issubclass(other, Joule):
            return pf
        elif issubclass(other, Hartree):
            return pf / 4.35974434e-18
        elif issubclass(other, Wavenumbers):
            return pf / (PlanckConstant.value * SpeedOfLight.in_units(Centimeters/Second).value)
        elif issubclass(other, KiloCaloriePerMol):
            return pf * AvogadrosNumber / 1000.0 / 4.184
        elif issubclass(other, KiloJoulePerMol):
            return pf * AvogadrosNumber / 1000.0
        elif issubclass(other, Hertz):
            return pf / PlanckConstant.value
        else: # pragma: no cover
            raise NotImplementedError("Conversion from units " + classname(cls) + " to units " + classname(other) + " is not implemented.")
EnergyUnit.reference_unit = Joule


class Wavenumber(EnergyUnit):
    __metaclass__ = Unit
EnergyUnit.default = Wavenumber

# TODO Molar energy unit?
def_units(EnergyUnit,
    #'ElectronVolt',
    'Hertz',
    'Hartree',
    'KiloCaloriePerMol',
    'KiloJoulePerMol',
    #------------------#
    Hartree_alias = 'AtomicUnitOfEnergy',
    Hartree_alias_plural = 'AtomicUnitsOfEnergy',
    #------------------#
    KiloCaloriePerMol_prefixed = False,   # Don't create prefixes, since e.g. MicroKCalPerMol doesn't make sense
    KiloCaloriePerMol_aliases = [
        'KiloCaloriePerMole',
        'KCalPerMol',
        'KcalPerMol',
    ],
    KiloCaloriePerMole_plural = 'KiloCaloriesPerMol',
    KcalPerMol_plural = 'KcalsPerMol',
    KCalPerMol_plural = 'KCalsPerMol',
    #------------------#
    KiloJoulePerMol_plural = 'KiloJoulesPerMol',
    KiloJoulesPerMol_prefixed = False,
    KiloJoulesPerMol_aliases = [
        'KJPerMol',
    ],
    KJPerMol_plural = False,
    #------------------#
)


##############
# Time Units #
##############

class TimeUnit(UnitGenre):
    """ General superclass for all time units
    """

class Second(TimeUnit):
    __metaclass__ = Unit
    @classmethod
    def to(cls, other):
        cls.genre_check(other)
        pf = cls.prefix_factor(other)
        if issubclass(other, Second):
            return pf
        elif issubclass(other, AtomicUnitOfTime):
            return pf / 2.418884326502e-17
        elif issubclass(other, Minute):
            return pf / 60.0
        elif issubclass(other, Hour):
            return pf / 3600.0
        elif issubclass(other, Day):
            return pf / 86400.0
        elif issubclass(other, Week):
            return pf / 604800.0
        elif issubclass(other, Year):
            return pf / 31556925.445
        elif issubclass(other, Decade):
            return pf / (31556925.445 * 10)
        elif issubclass(other, Century):
            return pf / (31556925.445 * 100)
        elif issubclass(other, Millennium):
            return pf / (31556925.445 * 1000)
        else: # pragma: no cover
            raise NotImplementedError("Conversion from units " + classname(cls) + " to units " + classname(other) + " is not implemented.")
TimeUnit.default = Second
TimeUnit.reference_unit = Second

# Just to demonstrate how the process works...
def_units(TimeUnit,
    'AtomicUnitOfTime',
    'Minute',
    'Hour',
    'Day',
    'Week',
    'Year',
    'Decade',
    'Century',
    'Millennium',
    AtomicUnitOfTime_plural = "AtomicUnitsOfTime",
    Century_plural = "Centuries",
    Millennium_plural = 'Millennia')

#########################
# Electric Charge Units #
#########################

class ElectricChargeUnit(UnitGenre):
    """ General superclass for all units of electric charge
    """

class Coulomb(ElectricChargeUnit):
    __metaclass__ = Unit
    @classmethod
    def to(cls, other):
        cls.genre_check(other)
        pf = cls.prefix_factor(other)
        if issubclass(other, Coulomb):
            return pf
        elif issubclass(other, AtomicUnitOfElectricCharge):
            return pf / ElementaryCharge.in_units(Coulomb).value
        else: # pragma: no cover
            raise NotImplementedError("Conversion from units " + classname(cls) + " to units " + classname(other) + " is not implemented.")
ElectricChargeUnit.default = Coulomb
ElectricChargeUnit.reference_unit = Coulomb

def_units(ElectricChargeUnit,
    'AtomicUnitOfElectricCharge',
    AtomicUnitOfElectricCharge_plural = 'AtomicUnitsOfElectricCharge',
    AtomicUnitOfElectricCharge_alias = 'AtomicUnitOfCharge',
    AtomicUnitOfElectricCharge_alias_plural = 'AtomicUnitsOfCharge',
)

##############
# Mass Units #
##############

class MassUnit(UnitGenre):
    """ General superclass for all units of mass
    """

class Gram(MassUnit):
    __metaclass__ = Unit
    @classmethod
    def to(cls, other):
        cls.genre_check(other)
        pf = cls.prefix_factor(other)
        if issubclass(other, Gram):
            return pf
        if issubclass(other, AtomicMassUnit):
            # NIST
            return pf / 1.660538921e-24
            # IUPAC
            #return pf / 1.6605402e-24
        elif issubclass(other, AtomicUnitOfMass):
            return pf / ElectronMass.in_units(Gram).value
        else: # pragma: no cover
            raise NotImplementedError("Conversion from units " + classname(cls) + " to units " + classname(other) + " is not implemented.")
MassUnit.reference_unit = Gram

class AtomicMassUnit(MassUnit):
    __metaclass__ = Unit
def_unit_alias('AMU', AtomicMassUnit)
MassUnit.default = AtomicMassUnit

class AtomicUnitOfMass(MassUnit):
    __metaclass__ = Unit
    __plural__ = 'AtomicUnitsOfMass'


#####################
# Dependent Imports #
#####################

from grendel.util.units.composite import CompositeUnit
from grendel.util.units.errors import IncompatibleUnitsError, UnknownUnitError
from grendel.util.units.value_with_units import ValueWithUnits
from grendel.util.units.physical_constants import ElectronMass, ElementaryCharge, PlanckConstant, SpeedOfLight, AvogadrosNumber, BohrRadius

