from functools import total_ordering
import math

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

    def __init__(self, mass, abundance, element):
        self.mass = mass
        self.abundance = abundance
        self.element = element

    ##############
    # Properties #
    ##############

    @property
    def symbol(self):
        return str(int(round(self.mass))) + self.element.symbol

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

    ##############
    # Attributes #
    ##############

    symbol = None
    atomic_number = None

    def __init__(self, symbol, atomic_number, isotopes):
        """
            The following example is pretty self explanitory hopefully:

            Element("H", 1,
                (
                    ( 1.00782503207, 0.999885 ),
                    ( 2.0141017778, 0.000115 ),
                    ( 3.0160492777, None )
                )
            )
        """
        self.symbol = symbol
        self.atomic_number = atomic_number
        self.isotopes = []
        for tup in isotopes:
            if isinstance(tup, tuple):
                self.isotopes.append(Isotope(tup[0], tup[1], self))
            elif isinstance(tup, Isotope):
                self.isotopes.append(tup)
            else:
                raise TypeError

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


