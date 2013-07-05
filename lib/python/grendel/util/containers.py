""" Special containers that don't logically belong somewhere else
"""
import __builtin__
from collections import defaultdict
from grendel.util.decorators import typechecked, IterableOf
from grendel.util.sentinal_values import ArgumentNotGiven


class CustomHashDict(object):

    ##############
    # Attributes #
    ##############

    default_factory = None
    factory_takes_key_argument = None
    hash_function = None
    equality_function = None

    ######################
    # Private Attributes #
    ######################

    _dict = None
    _len = None

    ##################
    # Initialization #
    ##################

    #TODO document this
    @typechecked(
        hash_function=(callable, None),
        pairs=(IterableOf(IterableOf(object)), None),
        default_factory=(callable, None),
        factory_takes_key_argument=(bool, None),
        equality_function=callable)
    def __init__(self,
            hash_function=None,
            pairs=None,
            default_factory=None,
            factory_takes_key_argument=False,
            equality_function=lambda x, y: x == y,
            **kwargs):
        """
        """
        self.hash_function = hash_function or __builtin__.hash
        self._dict = defaultdict(lambda: [])
        self.default_factory = default_factory
        self.factory_takes_key_argument = factory_takes_key_argument
        self.equality_function = equality_function
        self._len = 0
        if pairs is not None:
            for pair in pairs:
                self[pair[0]] = pair[1]
        self.update(kwargs)


    ###################
    # Special Methods #
    ###################

    def __getitem__(self, item):
        items = self._get_possibilites(item)
        for k, v in items:
            if self.__equality_function__(item, k):
                return v
        ret_val = self.__missing__(item)
        self[item] = ret_val
        return ret_val

    def __setitem__(self, key, value):
        pairs = self._get_possibilites(key)
        for i, (k, v) in enumerate(pairs):
            if self.__equality_function__(key, k):
                pairs[i] = (key, value)
                return
        # We need to add it...
        pairs_prim = self._get_primary_possibilities(key)
        pairs_prim.append((key, value))
        self._len += 1

    def __delitem__(self, key):
        pairs = self._get_possibilites(key)
        for i, (k, v) in enumerate(pairs):
            if self.__equality_function__(key, k):
                del pairs[i]
                self._len -= 1
                return
        raise KeyError(key)

    def __len__(self):
        return self._len

    def __contains__(self, item):
        items = self._get_possibilites(item)
        for k, v in items:
            if self.__equality_function__(item, k):
                return True
        return False

    def __iter__(self):
        for pairs in self._dict.values():
            for key, value in pairs:
                yield key

    def __missing__(self, key):
        """
        Allows CustomHashDict to act like a `collections.defaultdict`.  The default implementation
        works just like the `defaultdict` implementation.  Can be overridden in subclasses for
        a bit more control.
        """
        if self.default_factory is None:
            raise KeyError(key)
        elif self.factory_takes_key_argument:
            return self.default_factory(key)
        else:
            return self.default_factory()

    def __hash_function__(self, key):
        """
        Default implementation that can be overridden in subclasses.  By default, it calls the
        hash_function passed in to the constructor.  If no hash function was passed in, it calls
        Python's standard `hash`.  The hash function can also return a tuple, the first item
        of which is the hash to which the value will be assigned if the key is not equal to any
        of the values found with hashes corresponding to the hashes in the returned tuple.
        """
        return self.hash_function(key)

    def __equality_function__(self, in_key, found_key):
        """
        Default implementation that can be overridden in subclasses.  By default, it calls the
        equality_function passed in to the constructor.  If no equality function was passed in,
        it calls uses `key_a == key_b`
        """
        self.equality_function(in_key, found_key)

    ###########
    # Methods #
    ###########

    def items(self):
        for pairs in self._dict.values():
            for pair in pairs:
                yield pair

    def items(self):
        return [i for i in self.items()]

    def keys(self):
        for pair in self.items():
            yield pair[0]
    iter = keys

    def keys(self):
        return [k for k in self.keys()]

    def values(self):
        for pair in self.items():
            yield pair[1]

    def values(self):
        return [v for v in self.values()]

    def has_key(self, key):
        return key in self

    def clear(self):
        self._dict.clear()
        self._len = 0

    def pop(self, key, default=ArgumentNotGiven):
        try:
            ret_val = self[key]
            del self[key]
            return ret_val
        except KeyError as k:
            if default is not ArgumentNotGiven:
                return default
            else:
                raise k

    def update(self, other=None, **kwargs):
        # basically defer to super:
        to_add = dict(other, **kwargs) if other else dict(**kwargs)
        # then assign
        for key, val in to_add.items():
            self[key] = val
        return

    ###################
    # Private Methods #
    ###################

    def _get_possibilites(self, key):
        hsh = self.__hash_function__(key)
        if isinstance(hsh, tuple):
            return sum((self._dict[h] for h in hsh), [])
        else:
            return self._dict[hsh]

    def _get_primary_possibilities(self, key):
        hsh = self.__hash_function__(key)
        if isinstance(hsh, tuple):
            return self._dict[hsh[0]]
        else:
            return self._dict[hsh]




#--------------------------------------------------------------------------------#
#                   Least-recently-used Dictionary                               #
#--------------------------------------------------------------------------------#


## {{{ http://code.activestate.com/recipes/252524/ (r3)
class Node(object):
    __slots__ = ['prev', 'next', 'me']
    def __init__(self, prev, me):
        self.prev = prev
        self.me = me
        self.next = None

class LRUDict:
    """
    Implementation of a length-limited O(1) LRUDict queue.
    Built for and used by PyPE:
    http://pype.sourceforge.net
    Copyright 2003 Josiah Carlson.
    """
    def __init__(self, count, pairs=[]):
        self.count = max(count, 1)
        self.d = {}
        self.first = None
        self.last = None
        for key, value in pairs:
            self[key] = value

    def __contains__(self, obj):
        return obj in self.d

    def __getitem__(self, obj):
        a = self.d[obj].me
        self[a[0]] = a[1]
        return a[1]

    def __setitem__(self, obj, val):
        if obj in self.d:
            del self[obj]
        nobj = Node(self.last, (obj, val))
        if self.first is None:
            self.first = nobj
        if self.last:
            self.last.next = nobj
        self.last = nobj
        self.d[obj] = nobj
        if len(self.d) > self.count:
            if self.first == self.last:
                self.first = None
                self.last = None
                return
            a = self.first
            a.next.prev = None
            self.first = a.next
            a.next = None
            del self.d[a.me[0]]
            del a

    def __delitem__(self, obj):
        nobj = self.d[obj]
        if nobj.prev:
            nobj.prev.next = nobj.next
        else:
            self.first = nobj.next
        if nobj.next:
            nobj.next.prev = nobj.prev
        else:
            self.last = nobj.prev
        del self.d[obj]

    def __iter__(self):
        cur = self.first
        while cur is not None:
            cur2 = cur.next
            yield cur.me[1]
            cur = cur2

    def items(self):
        cur = self.first
        while cur is not None:
            cur2 = cur.next
            yield cur.me
            cur = cur2

    def keys(self):
        return iter(self.d)

    def values(self):
        for i,j in self.items():
            yield j

    def keys(self):
        return self.d.keys()

# end of http://code.activestate.com/recipes/252524/ }}}


