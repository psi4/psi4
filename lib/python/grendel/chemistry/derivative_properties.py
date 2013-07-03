from grendel.chemistry.molecular_properties import TensorProperty, MolecularProperty, Energy
from grendel.util.aliasing import function_alias
from grendel.representations.cartesian_representation import CartesianRepresentation
from grendel.util.exceptions import raises_error, ArgTypeError
from grendel.util.overloading import PartiallyConstructed
from grendel.util.parsing import MatrixRegexSequence
from grendel.util.units.errors import IncompatibleUnitsError

class RepresentationDependentProperty(MolecularProperty):
    """ A property whose value is dependent on a given representation of the molecule.
    """

    ##############
    # Attributes #
    ##############

    representation_getter = None

    ######################
    # Private Attributes #
    ######################

    _representation = None
    _needs_representation_instance = None

    ##############
    # Properties #
    ##############

    @property
    def representation(self):
        return self._representation

    @representation.setter
    def representation(self, new_val):
        if isinstance(new_val, type) and issubclass(new_val, Representation):
            if self.molecule is not None:
                if issubclass(new_val, CartesianRepresentation):
                    self._representation = self.molecule.cartesian_representation.frozen_copy()
                elif issubclass(new_val, InternalRepresentation):
                    raise NotImplementedError("need to be able to freeze the values in an"
                                              " InternalRepresentation in order to construct"
                                              " a property dependent on it.  This functionality"
                                              " is not yet implemented, but shouldn't be hard.")
                    #self._representation = self.molecule.internal_representation
                    ## Now make sure the molecule had an internal representation:
                    #if self._representation is None:
                    #    self._representation = new_val
                    #    self._needs_representation_instance = True
            else:
                self._representation = new_val
                self._needs_representation_instance = True
        elif not isinstance(new_val, Representation):
            raise ArgTypeError("new_val", new_val, Representation)
        else:
            self._representation = new_val
        if not self._needs_representation_instance:
            # Freeze it once it's assigned.  We don't want it changing on us...
            self._representation.freeze()

    @property
    def molecule(self):
        return self._molecule

    @molecule.setter
    def molecule(self, new_mol):
        self._molecule = new_mol
        if self._needs_representation_instance and self.representation is not None:
            # trigger the setter, which retrieves the Representation instance from the Molecule instance
            self.representation = self.representation

    # DEPRECATED
    def get_value(self, *args, **kwargs):
        pass


class DerivativeProperty(TensorProperty, RepresentationDependentProperty):

    ####################
    # Class Attributes #
    ####################

    known_subclasses = {}

    ##############
    # Attributes #
    ##############

    order = None

    ##################
    # Initialization #
    ##################

    def __new__(cls, *args, **kwargs):
        if cls is DerivativeProperty:
            return MolecularPropertyDerivative(*args, **kwargs)
        else:
            return super(DerivativeProperty, cls).__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        if type(self) is DerivativeProperty:
            raise NotImplementedError
        super(DerivativeProperty, self).__init__(*args, **kwargs)

    #################
    # Class Methods #
    #################

    @classmethod
    def with_respect_to(cls, rep):
        ret_val = cls()
        ret_val.representation = rep
        return ret_val

    ###########
    # Methods #
    ###########

    # DEPRECATED!
    def from_sequence(self, seq, units):
        """ Expects a MatrixRegexSequence which returns a Nx3 matrix
        """
        if self.order == 1:
            if isinstance(seq, MatrixRegexSequence):
                if units is not None:
                    val = (
                        RepresentationDependentTensor(seq.get_matrix().flatten('C'),
                                         rep=self.representation)
                        * units).in_units(self.units)
                    # Strip the units for the contents of _value
                    val.units = None
                    self._value = val
                    # Calling the property descriptor puts the units back on
                    return self.value_with_units
                elif self.units is None:
                    return seq.get_matrix()
                else:
                    raise IncompatibleUnitsError(self.units, units)
        else:
            raise NotImplementedError

    # deprecated
    def from_groups(self, groups, units):
        """ Placeholder implementation to satisfy abstract class requirements.  Everything
        is already taken care of in `from_sequence`
        """
        return NotImplemented


def MolecularPropertyDerivative(prop, representation=CartesianRepresentation, order=1):
    """ Factory for creating `MolecularProperty` subclasses that represent the nth derivative
     of a `MolecularProperty` with respect to some set of coordinates.

     Several possibilities for `representation`:  If it's a `Representation` subclass, then we want that representation
     of the molecule in question (used when requesting the derivative).  If it's a representation instance, then we
     want the Derivative with respect to that specific representation (raises an error if it is not compatible).
     If it's an arbitrary callable, then the callable is called with the contents of the file and should return
     a representation, which the (typically parsed-from-a-file) derivative is considered to be with respect to.
     The default is `CartesianRepresentation`.

     Derivative properties can be automatically or manually transformed to other representations.
    """
    if not issubclass(prop, MolecularProperty):
        raise TypeError
    if raises_error(int, order):
        raise TypeError
    name = '_{1}Derivative{0}'.format(order, prop.__name__)
    if name in DerivativeProperty.known_subclasses:
        ret_val = PartiallyConstructed(DerivativeProperty.known_subclasses[name])
    else:
        ret = type(name, (DerivativeProperty,), {})
        DerivativeProperty.known_subclasses[name] = ret
        ret_val = PartiallyConstructed(ret)
    ret_val.with_attributes(order=order, representation=representation)
    return ret_val

PropertyDerivative = function_alias('PropertyDerivative', MolecularPropertyDerivative)

Gradient = EnergyGradient = Forces = PropertyDerivative(Energy, order=1)
Hessian = EnergyHessian = PropertyDerivative(Energy, order=2)

#####################
# Dependent Imports #
#####################

from grendel.differentiation.derivative_tensor import RepresentationDependentTensor
from grendel.representations.internal_representation import InternalRepresentation
from grendel.representations.representation import Representation

