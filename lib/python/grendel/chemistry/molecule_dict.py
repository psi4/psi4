from __future__ import print_function
from collections import defaultdict
from copy import copy
import sys
import math

import numpy as np

from grendel import sanity_checking_enabled, type_checking_enabled
from grendel.util.containers import CustomHashDict
from grendel.util.decorators import typechecked, IterableOf
from grendel.util.exceptions import ProgrammerNeedsMoreCoffeeError
from grendel.util.metaprogramming import ReadOnlyAttribute
from grendel.util.strings import indented
from grendel.util.units import DistanceUnit
from grendel.util.units.value_with_units import strip_units


class MoleculeDict(object):

    #########################
    # Private Inner Classes #
    #########################

    class _HashValueContainer(object):
        value = None
        def __init__(self, val):
            self.value = val
        def __str__(self):
            return "_HashValueContainer({})".format(self.value)
        __repr__=__str__

    ####################
    # Class Attributes #
    ####################

    hash_digits = 6
    equality_digits = 7

    ##############
    # Attributes #
    ##############

    default_factory = None
    default = None
    key_representation = None

    ######################
    # Private Attributes #
    ######################

    _heads = None
    _len = None
    _keys = None
    _hash_mult = None

    ##################
    # Initialization #
    ##################

    @typechecked(
        key_representation='InternalRepresentation',
        pairs=(IterableOf(IterableOf(object)), None),
        default=(callable, None),
        default_factory=(callable, None),
        hash_digits=(int, None),
        equality_digits=(int, None)
    )
    def __init__(self,
            key_representation,
            pairs=None,
            default=None,
            default_factory=None,
            hash_digits=None,
            equality_digits=None,
            **kwargs):
        self.key_representation = key_representation
        if hash_digits is not None:
            # Only override the class-level default if an argument is given...
            self.hash_digits = hash_digits
        if equality_digits is not None:
            # Only override the class-level default if an argument is given...
            self.equality_digits = equality_digits
        self._heads = defaultdict(lambda: dict())
        self._hash_mult = 10**self.hash_digits
        self._keys = []
        self.default_factory = default_factory
        self.default = default
        if pairs is not None:
            for pair in pairs:
                self[pair[0]] = pair[1]
        self.update(**kwargs)

    ###################
    # Special Methods #
    ###################

    #---------------------#
    # Container Emulation #
    #---------------------#

    def __len__(self):
        return len(self._keys)

    def __contains__(self, item):
        if type_checking_enabled and not isinstance(item, Molecule):
            raise TypeError("only Molecule instances can be keys in a MoleculeDict")
        if sanity_checking_enabled:
            if item.ninternals != len(self.key_representation):
                raise ValueError(
                    "dimension mismatch.  Molecule with {} internals cannot be used as a key"
                    " with representation that has {} coordinates.".format(
                        item.ninternals,
                        len(self.key_representation)
                    ))
        #----------------------------------------#
        item = copy(item)
        #if item.cartesian_units is not DistanceUnit.default:
        item.convert_units(DistanceUnit.default)
        possibilities = self._get_possibilities(item, create=False)
        if possibilities is False:
            return False
        for key, valc in possibilities:
            if self._is_valid_key_for(key, item):
                return True
        return False

    def __getitem__(self, item):
        if type_checking_enabled and not isinstance(item, Molecule):
            raise TypeError("only Molecule instances can be keys in a MoleculeDict")
        if sanity_checking_enabled:
            if item.ninternals != len(self.key_representation):
                raise ValueError(
                    "dimension mismatch.  Molecule with {} internals cannot be used as a key"
                    " with representation that has {} coordinates.".format(
                        item.ninternals,
                        len(self.key_representation)
                    ))
        #----------------------------------------#
        # First, see if we already have the item we're looking for
        item = copy(item)
        #if item.cartesian_units is not DistanceUnit.default:
        item.convert_units(DistanceUnit.default)
        possibilities = self._get_possibilities(item)
        for key, valc in possibilities:
            if self._is_valid_key_for(key, item):
                return valc.value
        #----------------------------------------#
        # Not already here, see if we can generate an item
        if isinstance(item, MoleculeStub):
            key = item
        else:
            key = MoleculeStub.stub_for(item)
        val = self._get_default(item)
        self._keys.append(key)
        valc = MoleculeDict._HashValueContainer(val)
        possibilities.append((key, valc))
        return val

    def __setitem__(self, item, value):
        if type_checking_enabled and not isinstance(item, Molecule):
            raise TypeError("only Molecule instances can be keys in a MoleculeDict")
        if sanity_checking_enabled:
            if item.ninternals != len(self.key_representation):
                raise KeyError(
                    "dimension mismatch.  Molecule with {} internals cannot be used as a key"
                    " with representation that has {} coordinates.".format(
                        item.ninternals,
                        len(self.key_representation)
                    ))
        #----------------------------------------#
        # First, see if the key already exists
        item = copy(item)
        item.convert_units(DistanceUnit.default)
        possibilities = self._get_possibilities(item)
        for key, valc in possibilities:
            if self._is_valid_key_for(key, item):
                return valc.value
        #----------------------------------------#
        # Not already here, so add it
        if isinstance(item, MoleculeStub):
            key = item
        else:
            key = MoleculeStub.stub_for(item)
        valc = MoleculeDict._HashValueContainer(value)
        possibilities.append((key,valc))
        self._keys.append(key)

    def __delitem__(self, key):
        raise NotImplementedError

    def __iter__(self):
        for key in self._keys:
            yield key


    ###########
    # Methods #
    ###########

    def keys(self):
        for key in self:
            yield key

    def keys(self):
        return list(self.keys())

    def items(self):
        for key in self:
            yield key, self[key]

    def items(self):
        return list(self.items())

    def values(self):
        for k, v in self.items():
            yield v

    def values(self):
        return list(self.values())

    def update(self, **kwargs):
        for key, val in kwargs.items():
            self[key] = val

    ###################
    # Private Methods #
    ###################

    def _is_valid_key_for(self, key, item):
        if key is item:
            return True
        keysyms = tuple(a.full_symbol for a in key)
        itsyms = tuple(a.full_symbol for a in item)
        if keysyms != itsyms:
            return False
        kvals = strip_units(self.key_representation.values_for_molecule(key),
            DistanceUnit.default,
            assume_units=key.cartesian_units)
        itvals = strip_units(self.key_representation.values_for_molecule(item),
            DistanceUnit.default,
            assume_units=item.cartesian_units)
        cutoff = 10**(-self.equality_digits)
        return np.all(abs(kvals - itvals) < abs(cutoff*kvals))

    def _get_possibilities(self, item, create=True):
        itsyms = tuple(a.full_symbol for a in item)
        items = self._heads[itsyms]
        itvals = strip_units(self.key_representation.values_for_molecule(item),
            DistanceUnit.default,
            assume_units=item.cartesian_units)
        itvals = [i * self._hash_mult for i in itvals]
        return self._r_get_possibilities(items, itvals, create)

    # TODO this could be optimized by passing just a position in itvals rather than duplicating the list each recursion
    def _r_get_possibilities(self, currdict, itvals, create):
        if len(itvals) == 1:
            cl, flr = int(math.ceil(itvals[0])), int(math.floor(itvals[0]))
            if cl in currdict:
                if flr not in currdict:
                    # This could be done more efficiently than just connecting neighboring nodes...but shouldn't be a problem
                    currdict[flr] = currdict[cl]
                elif currdict[flr] is not currdict[cl]:
                    # We have a neighboring value collision.  I could handle this, or just do this for now:
                    raise NotImplementedError("MoleculeDict hash collision.  Make Molecule.hash_digits a larger number (but not too big)")
                return currdict[cl]
            elif flr in currdict:
                if cl not in currdict:
                    currdict[cl] = currdict[flr]
                elif currdict[flr] is not currdict[cl]:
                    # We have a neighboring value collision.  I could handle this, or just do this for now:
                    raise NotImplementedError("MoleculeDict hash collision.  Make Molecule.hash_digits a larger number (but not too big)")
                return currdict[flr]
            else:
                if create:
                    rv = []
                    currdict[cl] = rv
                    currdict[flr] = rv
                    return rv
                else:
                    return False
        else:
            cl, flr = int(math.ceil(itvals[0])), int(math.floor(itvals[0]))
            if cl in currdict:
                if flr not in currdict:
                    currdict[flr] = currdict[cl]
                elif currdict[flr] is not currdict[cl]:
                    # We have a neighboring value collision.  I could handle this, or just do this for now:
                    raise NotImplementedError("MoleculeDict hash collision.  Make Molecule.hash_digits a larger number (but not too big)")
                recurs_dict = currdict[cl]
            elif flr in currdict:
                if cl not in currdict:
                    currdict[cl] = currdict[flr]
                elif currdict[flr] is not currdict[cl]:
                    # We have a neighboring value collision.  I could handle this, or just do this for now:
                    raise NotImplementedError("MoleculeDict hash collision.  Make Molecule.hash_digits a larger number (but not too big)")
                recurs_dict = currdict[flr]
            else:
                if create:
                    recurs_dict = dict()
                    currdict[cl] = currdict[flr] = recurs_dict
                else:
                    return False
            return self._r_get_possibilities(recurs_dict, itvals[1:], create)

    def _get_default(self, item):
        valset = False
        val = None
        if self.default_factory is not None:
            try:
                val = self.default_factory(item)
                valset = True
            except:
                if self.default is not None:
                    val = self.default()
                    valset = True
        elif self.default is not None:
            val = self.default()
            valset=True
        if not valset:
            raise KeyError(item)
        else:
            return val

    def _key_chain_str(self, item):
        itsyms = tuple(a.full_symbol for a in item)
        items = self._heads[itsyms]
        itvals = strip_units(self.key_representation.values_for_molecule(item),
            DistanceUnit.default,
            assume_units=item.cartesian_units)
        itvals = [i * self._hash_mult for i in itvals]
        return self._r_key_chain_str(items, itvals, "")

    def _r_key_chain_str(self, currdict, itvals, currstr):
        if len(itvals) > 0:
            cl, flr = int(math.ceil(itvals[0])), int(math.floor(itvals[0]))
            if cl in currdict:
                recur_dict = currdict[cl]
                currstr += str(cl) + " > "
            elif flr in currdict:
                recur_dict = currdict[flr]
                currstr += str(flr) + " > "
            else:
                currstr += "missing ({} and {})".format(cl, flr)
                return currstr
            return self._r_key_chain_str(recur_dict, itvals[1:], currstr)
        else:
            return currstr + "<end> => {}".format(currdict)

#####################
# Dependent Imports #
#####################

from grendel.chemistry.molecule import MoleculeStub, Molecule
if type_checking_enabled:
    from grendel.representations.internal_representation import InternalRepresentation

