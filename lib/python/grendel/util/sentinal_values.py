""" A module containing blank (or almost blank) sentinal value classes to be used throughout the program.
"""

from grendel.util.strings import classname


__all__ = [
]

class Keyword(object):
    """ Simple container for keyword-like objects ("sentinal values").
    Can also be used as a container for sub-keywords via instance objects.
    Sub-keywords are case-insensative.  Keywords cannot be assigned values using
    the `value` attribute (note that this means you can't create a sub-keyword called 'value').

    .. WARNING::
      Since there is exactly one instance of every unique keyword (with
      two keywords that differ by capitalization *not* being unique), storing a value can have some
      unexpected consequences.  For instance::

          >>> Hello = Keyword('Hello')
          >>> foo = Hello.World
          >>> foo.value = "good bye"
          >>> # ...some time later...
          >>> Hello.WORLD.value = "good night"
          >>> foo.value
          'good night'

    .. note::
      Since `value`, `parent`, `_parent`, `parents`, `_name`, and any of the double-underscore methods inherited from
      `object` are already attributes of `Keyword` objects, you can't make sub-keywords with these names.  Everything
      else is fair game.


    :Examples:

    >>> Hello = Keyword('Hello')
    >>> Hello.World
    Hello.World
    >>> # Any future changes in case simply refer back to the original
    >>> Hello.WoRlD
    Hello.World
    >>> # But, of course, not for the first level
    >>> hEllO.World
    Traceback (most recent call last):
       ...
    NameError: name 'hEllO' is not defined
    >>> # Multiple levels can be created at once
    >>> Hello.there.world
    Hello.there.world
    >>> Hello.there == Hello.ThErE
    True
    >>> Hello.there == Hello.ThErE.world
    False
    >>> Hello.there.WORLD == Hello.ThErE.world
    True
    >>> Hello.There.World.value = "Hello World"
    >>> Hello.THERE.WORLD.value
    'Hello World'


    """

    ####################
    # Class Attributes #
    ####################

    _known_roots = []


    ##############
    # Attributes #
    ##############

    value = None
    """ Used for assignment of values to keywords. """


    ######################
    # Private Attributes #
    ######################

    _name = None
    _parent = None


    ##################
    # Initialization #
    ##################

    # I would want to implement this if for some reason I wanted root keywords to be unique
    #def __new__(cls, name, parent = None):
    #    self = object.__new__(cls)
    #    if
    #    self.__init__(name, parent)
    #    return self


    def __init__(self, name, parent = None):
        self._name = name
        self._parent = parent
        if parent is None:
            Keyword._known_roots.append(self)


    ##############
    # Properties #
    ##############

    @property
    def parent(self):
        """ The keyword immediately preceding `self` in the heirarchical keyword notation.
        For instance:
            >>> Hello = Keyword('Hello')
            >>> Hello.World.parent
            Hello
            >>> Hello.There.World.parent
            Hello.There

        """
        return self._parent

    @parent.setter
    def parent(self, new_parent):
        if new_parent is None:
            if self not in Keyword._known_roots:
                Keyword._known_roots.append(self)
            self._parent = new_parent
        else:
            if self in Keyword._known_roots:
                Keyword._known_roots.remove(self)
            self._parent = new_parent

    @property
    def parents(self):
        """ list of all parent keywords, ordered from immediate parent to most distant
        """
        spot = self
        ret_val = []
        while spot._parent is not None:
            ret_val.append(spot._parent)
            spot = spot._parent
        return ret_val


    ###################
    # Special Methods #
    ###################

    def __str__(self):
        if self._parent is not None:
            return ".".join(p._name for p in reversed(self.parents)) + "." + self._name
        else:
            return self._name
    __repr__ = __str__

    def __eq__(self, other):
        if isinstance(other, basestring):
            return self._name == other
        if not isinstance(other, self.__class__):
            raise TypeError("Invalid equality comparison of {0} to {1}".format(
                self.__class__.__name__,
                other.__class__.__name__))
        elif self.value is not None or other.value is not None:
            return (self._name.upper(), self.value) == (other._name.upper(), other.value)
        else:
            return self._name.upper() == other._name.upper()

    def __hash__(self):
        return self._name.upper().__hash__()

    def __getattr__(self, item):
        if item.upper() not in self.__dict__:
            self.__dict__[item.upper()] = Keyword(item, self)
            return self.__dict__[item.upper()]
        else:
            return self.__dict__[item.upper()]

    def __setattr__(self, key, value):
        if key in self.__class__.__dict__:
            self.__dict__[key] = value
        else:
            raise SyntaxError("Sentinal value dummy objects can't be assigned values.  "
                              "Check your code; you're doing something that doesn't make sense.")


# TODO Remove this class
class SentinalValue(object):
    """
    DON'T USE THIS; IT BREAKS PARALLELISM
    SentinalValue objects mimic builtins like None and NotImplemented.  There is exactly
    one instance of a SentinalValue with a given name, so "is" and "is not" can be
    safely used, just as with None and NotImplemented.  SentinalValues can also be pickled
    and unpickled without changing this behavior.
    """
    known_values = dict()

    def __new__(cls, name):
        if name in cls.known_values:
            return cls.known_values[name]
        else:
            rv = object.__new__(cls)
            rv.name = name
            cls.known_values[rv] = rv
            return rv

    def __repr__(self):
        return self.name

    def __hash__(self):
        return hash(self.name)

    #def __eq__(self, other):
    #    if isinstance(other, SentinalValue):
    #        return self.name == other.name
    #    else:
    #        return False

    class _SVUnloader(object):
        def __init__(self, name):
            self.name = name
        def __call__(self):
            return SentinalValue(self.name)

    def __reduce__(self):
        return SentinalValue._SVUnloader(self.name), tuple()

###############################################################
# Some keywords that make sense to define for use many places #
###############################################################

# DEPRECATED: Don't use these.  They don't work as expected 
#   in parallel contexts
All = SentinalValue("All")
ArgumentNotGiven = SentinalValue("ArgumentNotGiven")
