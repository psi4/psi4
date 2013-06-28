from grendel.util.decorators import CachedProperty
import grendel as g

class Derivative(object):
    """ Encapsulates the concept of an (energy) derivative with respect to a coordinate.
    This may later be expanded (in a backwards compatible manner) to include other types of derivatives,
    e.g. dipole moment derivatives

    Attributes
    ----------
    of : subclass of Differentiable
        The entity that the derivative is being taken of (usually the energy of the molecule, but can be anything
        that is a subclass of the abstract base class Differentiable
    value : float
        The value of the derivative, in the units given by `units`
    units : `CompositeUnit` or some class with a `Unit` metaclass
        The units of the derivative
    coords : list of Coordinate
        The coordinate(s) that the derivative is taken with respect to (multiple if the derivative is higher than first order)


    """



    ##############
    # Attributes #
    ##############

    of = None
    units = None
    coords = None


    ##################
    # Initialization #
    ##################

    def __init__(self, coords):
        """
        """
        # Sanity check:  All coordinates must be from the same molecule
        if g.sanity_checking_enabled:
            for i in xrange(1, len(coords)):
                if not coords[i].molecule is coords[0].molecule:
                    raise ValueError("All coordinates a derivative is taken with respect to must refer to the"
                                     " same Molecule object (not even a copy thereof)")
        if g.type_checking_enabled:
            for coord in coords:
                if not isinstance(coord, Coordinate):
                    raise TypeError("Argument 'coords' given to Derivative constructor must be a list of objects that "
                                    "are subclasses of Coordinate")
        self.coords = coords



    ##############
    # Properties #
    ##############

    @property
    def order(self):
        """The order of the derivative.  '1' indicates first derivative, '2' indicates second derivative, etc.
        """
        return len(self.coords)

    @CachedProperty
    def value(self):
        """ The value of the derivative (as a float), in the units given by `units`
        """
        return NotImplemented

    @property
    def molecule(self):
        """ The `Molecule` object to which the coordinate
        """
        return self.coords[0].molecule

#####################
# Dependent Imports #
#####################

from grendel.coordinates.coordinate import Coordinate

