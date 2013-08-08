""" Special containers that don't logically belong somewhere else
"""
from __future__ import print_function
import __builtin__
from collections import defaultdict
from grendel.util.decorators import typechecked, IterableOf
from grendel.util.misc import have_python3
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

    def iteritems(self):
        for pairs in self._dict.itervalues():
            for pair in pairs:
                yield pair

    def items(self):
        return [i for i in self.iteritems()]

    def iterkeys(self):
        for pair in self.iteritems():
            yield pair[0]
    iter = iterkeys

    def keys(self):
        return [k for k in self.iterkeys()]

    def itervalues(self):
        for pair in self.iteritems():
            yield pair[1]

    def values(self):
        return [v for v in self.itervalues()]

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
        for key, val in to_add.iteritems():
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

    def iteritems(self):
        cur = self.first
        while cur is not None:
            cur2 = cur.next
            yield cur.me
            cur = cur2

    def iterkeys(self):
        return iter(self.d)

    def itervalues(self):
        for i,j in self.iteritems():
            yield j

    def keys(self):
        return self.d.keys()
# end of http://code.activestate.com/recipes/252524/ }}}


#--------------------------------------------------------------------------------#
#                              AliasedDict                                       #
#--------------------------------------------------------------------------------#

class AliasedDict(dict):
    """
    A dictionary with multiple keys pointing to a given value
    """

    _sequence_key_types = (tuple, set, frozenset)

    #region | Initialization |

    def __init__(self, initial_dict=None):
        self._equivalent_keys = {}
        self._key_sets = set()
        self._first_keys = dict()
        super(AliasedDict, self).__init__()
        if initial_dict is not None:
            for keys, value in initial_dict.items():
                if isinstance(keys, AliasedDict._sequence_key_types):
                    keyset = frozenset(keys)
                    if isinstance(keys, tuple):
                        self._first_keys[keyset] = keys[0]
                    self._key_sets.add(keyset)
                    for key in keys:
                        if key in self._equivalent_keys:
                            raise KeyError("Key '{0}' is in multiple key lists".format(key))
                        if isinstance(key, AliasedDict._sequence_key_types):
                            raise KeyError("For simplicity, a key in an AliasedDict cannot be one of {0}".format(
                                AliasedDict._sequence_key_types
                            ))
                        self._equivalent_keys[key] = keyset
                        super(AliasedDict, self).__setitem__(key, value)
                else:
                    key = keys
                    keyset = frozenset((key,))
                    self._first_keys[keyset] = key
                    self._equivalent_keys[key] = keyset
                    self._key_sets.add(keyset)
                    super(AliasedDict, self).__setitem__(key, value)

    #endregion

    #========================================#

    #region | Special Methods |

    # The methods in this section and the next correspond to the
    #   order in which the mapping type operations are introduced in
    #   the python language documentation,
    #   http://docs.python.org/2/library/stdtypes.html#mapping-types-dict

    class _AliasedDictUnloader(object):
        def __init__(self, adict):
            self.pickled_state = dict()
            for keyset in adict._key_sets:
                if keyset in adict._first_keys:
                    key = [adict._first_keys[keyset]]
                    key += [k for k in keyset if k is not key[0]]
                    key = tuple(key)
                else:
                    key = keyset
                self.pickled_state[key] = adict[keyset]
        def __call__(self):
            return AliasedDict(self.pickled_state)

    def __reduce__(self):
        return AliasedDict._AliasedDictUnloader(self)

    def __len__(self):
        return len(self._key_sets)

    def __getitem__(self, item):
        if isinstance(item, AliasedDict._sequence_key_types):
            item = frozenset(item)
            if item not in self._key_sets:
                raise KeyError("Key '{0}' not found".format(item))
                # Just retrieve the first value
            return super(AliasedDict, self).__getitem__(list(item)[0])
        else:
            return super(AliasedDict, self).__getitem__(item)

    def __setitem__(self, key, value):
        if isinstance(key, AliasedDict._sequence_key_types):
            keyset = frozenset(key)
            if isinstance(key, tuple):
                self._first_keys[keyset] = key[0]
            if keyset not in self._key_sets:
                # It's a new set of aliases
                self._key_sets.add(keyset)
                for k in key:
                    if k in self._equivalent_keys:
                        raise KeyError("Key '{0}' is already in key list {1}, and"
                                         " {1} is not the same as given key list {2}".format(
                            k, self._equivalent_keys[k], key
                        ))
                    if isinstance(k, AliasedDict._sequence_key_types):
                        raise KeyError("For simplicity, a key in an AliasedDict cannot be one of {0}".format(
                            AliasedDict._sequence_key_types
                        ))
                    self._equivalent_keys[k] = keyset
                    super(AliasedDict, self).__setitem__(k, value)
            else:
                # we already have a value for the keyset, just assign each alias
                for k in keyset:
                    super(AliasedDict, self).__setitem__(k, value)
        elif key in self._equivalent_keys:
            # It's an existing key.  Set all of the aliases
            for k in self._equivalent_keys[key]:
                super(AliasedDict, self).__setitem__(k, value)
        else:
            # It's a single key that we don't already have.  Just set the value
            keyset = frozenset((key,))
            self._equivalent_keys[key] = keyset
            self._key_sets.add(keyset)
            self._first_keys[keyset] = key
            super(AliasedDict, self).__setitem__(key, value)

    def __delitem__(self, key):
        if isinstance(key, AliasedDict._sequence_key_types):
            keyset = frozenset(key)
        elif key in self._equivalent_keys:
            keyset = self._equivalent_keys[key]
        else:
            keyset = frozenset((key,))
        if keyset not in self._key_sets:
            raise KeyError("Key '{0}' not found.".format(key))
        for k in keyset:
            super(AliasedDict, self).__delitem__(k)
            del self._equivalent_keys[k]
        if keyset in self._first_keys:
            del self._first_keys[keyset]
        self._key_sets.remove(keyset)

    def __contains__(self, item):
        if isinstance(item, AliasedDict._sequence_key_types):
            return frozenset(item) in self._key_sets
        else:
            return item in self._equivalent_keys

    def __repr__(self):
        rv = type(self).__name__ + "("
        rv += repr(dict(
            (tuple(ks), self[ks]) for ks in self._key_sets
        ))
        return rv + ")"

    #endregion

    #========================================#

    #region | Methods |

    def iterkeys(self):
        for keyset in self._key_sets:
            yield keyset
    __iter__ = iterkeys

    def clear(self):
        self._key_sets.clear()
        self._equivalent_keys.clear()
        super(AliasedDict, self).clear()

    def copy(self):
        init_dict = dict([
            (keyset, self[keyset]) for keyset in self._key_sets
        ])
        return type(self)(init_dict)

    def get(self, k, d=None):
        if k in self:
            return self[k]
        else:
            return d

    def has_key(self, k):
        return k in self

    # TODO Return dictview objects for items(), keys(), and values() if have_python3

    def items(self):
        return list(self.iteritems())

    def iteritems(self):
        for keyset in self.iterkeys():
            yield keyset, self[keyset]

    def iterkeys(self):
        return iter(self._key_sets)

    def itervalues(self):
        for keyset, value in self.iteritems():
            yield value

    def keys(self):
        return list(self.iterkeys())

    def pop(self, k, d=ArgumentNotGiven):
        if k not in self:
            if d is ArgumentNotGiven:
                # This will raise an error.  Let it:
                return self[k]
            else:
                return d
        else:
            rv = self[k]
            del self[k]
            return rv

    def popitem(self):
        if len(self) == 0:
            raise KeyError("popitem(): dictionary is empty")
        key, val = self.items()[0]
        del self[key]
        return key, val

    def setdefault(self, k, d=None):
        if k in self:
            return self[k]
        else:
            self[k] = d
            return d

    def update(self, E=None, **F):
        if E is not None:
            for key, val in E.items():
                self[key] = val
        for key, val in F.items():
            self[key] = val

    def values(self):
        return list(self.itervalues())

    # TODO viewitems(), viewkeys(), viewvalues()

    def viewkeys(self):
        raise NotImplementedError()

    def viewitems(self):
        raise NotImplementedError()

    def viewvalues(self):
        raise NotImplementedError()

    def firstkey(self, key):
        if isinstance(key, AliasedDict._sequence_key_types):
            keyset = frozenset(key)
        else:
            if key not in self._equivalent_keys:
                raise KeyError("Key not found: {0}".format(key))
            keyset = self._equivalent_keys[key]
        if keyset in self._first_keys:
            return self._first_keys[keyset]
        elif keyset in self._key_sets:
            raise KeyError("Set of keys {0} was not ordered when it was passed in".format(keyset))
        else:
            raise KeyError("Keyset not found: {0}".format(keyset))

    def get_simple_dict(self):
        rv = dict()
        for keyset, value in self.items():
            if keyset in self._first_keys:
                rv[self._first_keys[keyset]] = value
            else:
                # Just choose one
                rv[iter(keyset).next()] = value
        return rv



    #endregion

    #========================================#

    # End of AliasedDict class
    pass



#TODO AccessibleSet as a dict of item->item acting as a subclass of set
