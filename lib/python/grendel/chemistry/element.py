from collections import Iterable
from functools import total_ordering
import math
from grendel.util.units.value_with_units import has_units
from grendel.util.units.composite import CompositeUnit

__all__ = [
    "Element",
    "Isotope"
]

# TODO __str__ for Element and Isotope


# TODO Nuclear spin attribute
@total_ordering
class Isotope(object):
    """
    Struct for containing the information about an isotope.
    For abundance, ``None`` means trace abundance
    """

    ##############
    # Attributes #
    ##############

    abundance = None
    element = None
    mass = None
    spin = None

    ##################
    # Initialization #
    ##################

    def __init__(self,
            mass,
            mass_uncertainty,
            abundance,
            abundance_uncertainty,
            mass_number,
            special_symbol=None
    ):
        self.mass = mass
        self.mass_uncertainty = mass_uncertainty
        self.abundance = abundance
        self.abundance_uncertainty = abundance_uncertainty
        self.element = None
        self.mass_number = mass_number
        self.special_symbol = special_symbol

    ##############
    # Properties #
    ##############

    @property
    def symbol(self):
        if self.special_symbol is not None:
            return self.special_symbol
        else:
            if self.element is not None:
                return str(self.mass_number) + self.element.symbol
            else:
                return str(self.mass_number) + "(unknown element)"

    ###################
    # Special Methods #
    ###################

    def __repr__(self):
        return "Isotope(" + repr(self.mass) + ", " + repr(self.abundance) + ")"

    def __str__(self):
        return self.symbol

    def __eq__(self, other):
        return (self.mass, self.spin) == (other.mass, other.spin)

    def __lt__(self, other):
        return (self.mass, self.spin) < (other.mass, other.spin)

    def __hash__(self):
        return hash(self.mass) + hash(self.spin)

    ###########
    # Methods #
    ###########

    #-------------------------------------#
    # Inquiry methods (which return bool) #
    #-------------------------------------#

    def is_principal(self):
        return self.element.principal_isotope is self

class Element(object):
    """Encapsulates an element.  Contains attributes that are constant for a given element.

    """

    #--------------------------------------------------------------------------------#

    #region | Attributes                                                                {{{1 |

    symbol = None
    """
    TODO Document this attribute
    """

    atomic_number = None
    """
    TODO Document this attribute
    """

    isotopes = None
    """
    TODO Document this attribute
    """

    atomic_weight = None
    """
    TODO Document this attribute
    """

    atomic_weight_uncertainty = None
    """
    TODO Document this attribute
    """

    is_synthetic = None
    """
    TODO Document this attribute
    """

    ionization_energies = None
    """
    TODO Document this attribute
    """

    vdw_radius = None
    """
    TODO Document this attribute
    """

    electron_affinity = None
    """
    TODO Document this attribute
    """

    electronegativity = None
    """
    TODO Document this attribute
    """

    #endregion }}}1

    #--------------------------------------------------------------------------------#

    #region | Initialization                                                            {{{1 |

    def __init__(self,
            symbol,
            atomic_number,
            isotopes,
            atomic_weight,
            atomic_weight_uncertainty=None,
            is_synthetic=False,
            **kwargs
    ):
        """
        """
        self.symbol = symbol
        self.atomic_number = atomic_number
        self.isotopes = []
        self.atomic_weight = atomic_weight
        self.atomic_weight_uncertainty = atomic_weight_uncertainty
        self.is_synthetic = is_synthetic
        for iso in isotopes:
            self.add_isotope(iso)
        for kw in kwargs:
            if kw in self.__class__.__dict__:
                val = kwargs[kw]
                if has_units(val) and not isinstance(val.units, CompositeUnit):
                    val = val.in_units(val.units.genre.default)
                elif isinstance(val, Iterable):
                    new_vals = []
                    val = list(val)
                    for v in val:
                        if has_units(v) and not isinstance(v.units, CompositeUnit):
                            new_vals.append(v.in_units(v.units.genre.default))
                        else:
                            new_vals.append(v)
                    val = new_vals
                setattr(self, kw, val)
            else:
                raise ValueError("Unknown Element constructor keyword '{}'".format(kw))

    #endregion }}}1

    #--------------------------------------------------------------------------------#
    ##############
    # Properties #
    ##############

    @property
    def principal_isotope(self):
        """The most abundant isotope.  This contains the default mass for the element."""
        def abundance(iso): return iso.abundance if not iso.abundance is None else 0.0
        return max(self.isotopes, key=abundance)

    ###################
    # Special Methods #
    ###################

    def __repr__(self):
        return "Element(" + repr(self.symbol) + ", " + repr(self.atomic_number) \
            + ", [" + ", ".join(map(repr, self.isotopes)) + "])"

    def __str__(self):
        return self.symbol

    def __eq__(self, other):
        """
        True if the symbols are the same.  Does not check isotope list!  You should not be creating copies of elements
        anyway.  Modify the instance created in ElementData if you need to add an isotope.
        """
        if not isinstance(other, Element):
            raise TypeError("Element object compared to non-Element type.")
        return self.symbol == other.symbol

    def __hash__(self):
        return hash(self.symbol)

    def add_isotope(self, iso):
        iso.element = self
        self.isotopes.append(iso)

