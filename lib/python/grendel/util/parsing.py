""" Tools for parsing files.
Be a bit careful with this...there are still some edge cases that fail for RegexSequence.  I've tried to create
a test set of the test cases that fail (marked with skips in the test suite), but there may be some that I've missed.
"""
from abc import ABCMeta
from collections import Iterable, Sequence, Sized, namedtuple
from copy import copy, deepcopy
from itertools import chain
import re
import types
from warnings import warn

from grendel.gmath.matrix import Matrix
from grendel.util.decorators import with_flexible_arguments
from grendel.util.metaprogramming import ImportableClass
from grendel.util.overloading import pop_kwarg

__all__ = [
    "RegexSequence",
    "Regex",
    "Duplicated",
    "Repeated",
    "RepeatedRegex"
]

re_float = r'-?\d+\.\d+'
re_int = r'-?\d+'
re_scientific = r'-?\d+\.\d+[dDeE][-\+]?\d+'
re_float_or_int = r'(?:{float}|{int})'.format(float=re_float, int=re_int)
re_float_or_int_or_scientific = r'(?:{float}|{int}|{scientific})'.format(float=re_float, int=re_int, scientific=re_scientific)
re_float_or_scientific = r'(?:{float}|{scientific})'.format(float=re_float, scientific=re_scientific)




class Regex(object):
    """ Abstract base class for regular expressions.
    For now, anything that responds to 'match', 'search', and has the attributes 'flags' and 'pattern' is
    considered a Regex object.  You can't instantiate a `Regex` object; the `__new__()` just returns `re.compile()`
    with the same arguments.

    Examples
    --------
    >>> import re
    >>> isinstance(re.compile(r'test123'), Regex)
    True
    >>> isinstance(r'test123', Regex)
    False

    """
    __metaclass__ = ABCMeta

    def __new__(cls, *args, **kwargs):
        return re.compile(*args, **kwargs)

    @classmethod
    def __subclasshook__(cls, C):
        if cls is Regex:
            if any('match' in B.__dict__ and callable(B.__dict__['match']) for B in C.__mro__) and \
                    any('search' in B.__dict__ and callable(B.__dict__['search']) for B in C.__mro__) and \
                    any(hasattr(B, 'flags') for B in C.__mro__) and any(hasattr(B, 'pattern') for B in C.__mro__):
                return True
        return NotImplemented

CompiledRegexType = type(re.compile(r'.'))

#TODO complete coverage
#TODO Multiline mode
#TODO retrieval of named capture groups
#TODO More thorough documentation
#TODO arbitrarily nestable repeated sequences (NEEDED FOR MATRIX PARSING)
class RegexSequence(object):
    """ A sequence of regular expressions for parsing a piece of data from a file.
    The sequence can have optional activation ("start") and deactivation ("stop") sequences, and matches are collected
    from every match of the main regex (or main list of regexes) given.  The Ellipsis object may be given
    in the main sequence (or the start/stop sequences) to indicate that any number of lines may come in between the
    preceeding and succeding regexes.  A RegexSequence.RepeatedRegex object can be used to indicate that a regex can be
    repeated without invalidating the sequence.

    The RegexSequence object is used by passing a file in line by line to the `parse_line()` method.

    If the `last_only` keyword argument is given to the initialization, only the last match (or sequence of matches,
     if a list of regexes is given for the main sequence) is kept

    """

    class RepeatedRegex(object):
        """ A regex that may be repeated on multiple consecutive lines without invalidating the sequence.
        The RepeatedRegex may be repeated an arbitrary (repeats = None, the default) or specified number of times.
        Currently, the behavior of repeated regexes with exact=False is not stable when multiple (different) RepeatedRegex
        objects are given consecutively and the number of intervening lines to be matched is fewer than the number of
        consecutive RepeatedRegex objects minus 1.  This is a fairly extreme edge case that is demonstrated in the
        test suite in test_expected_failures() in the file "test_parsing.py".
        """

        ##############
        # Attributes #
        ##############

        repeats = None
        regexes = None
        exactly = None
        min_repeats = None
        greedy = None


        ######################
        # Private Attributes #
        ######################

        _counter = None
        _spot = None
        _was_okay_to_finish = None


        ##################
        # Initialization #
        ##################

        def __init__(self, regexes, repeats=None, min_repeats=1, flags=0, exactly=False, greedy=True):
            if isinstance(regexes, basestring):
                self.regexes = [re.compile(regexes, flags)]
            elif isinstance(regexes, Regex):
                # Note that flags are ignored if regexes are passed in like this
                self.regexes = [regexes]
            elif isinstance(regexes, RepeatedRegex):
                raise NotImplementedError("Nested RepeatedRegexes are not yet implemented because they would be horrific to implement.  It's on my TODO list, though...")
            elif isinstance(regexes, Iterable):
                self.regexes = []
                for i in regexes:
                    if isinstance(i, basestring):
                        self.regexes.append(re.compile(i, flags))
                    elif isinstance(i, RepeatedRegex):
                        raise NotImplementedError("Nested RepeatedRegexes are not yet implemented because they would be horrific to implement.  It's on my TODO list, though...")
                    elif isinstance(i, Regex):
                        # Note that flags are ignored for regexes passed in like this
                        self.regexes.append(i)
            else:
                raise TypeError("Invalid type for regex given to RepeatedRegex constructor: " + str(regexes))
            self.repeats = repeats
            if repeats is not None:
                if not isinstance(repeats, int) or repeats < 1:
                    raise TypeError("Number of repeats in repeated regex must be a positive integer (if any number is specified)")
            self.exactly = exactly
            self.min_repeats = min_repeats
            self.greedy = greedy
            self.reset()


        ###################
        # Special Methods #
        ###################

        def __deepcopy__(self, memo):
            # no need to copy the individual regexes, they're immutable
            # but RegexSequence.Duplicate objects need to be copied;
            # do call deepcopy on everything but compiled regexes
            new_regexes = [deepcopy(r, memo) if not isinstance(r, CompiledRegexType) else r for r in self.regexes]
            # everything else is immutable
            return self.__class__(new_regexes, self.repeats, self.min_repeats, 0, self.exactly, self.greedy)

        def __copy__(self):
            # just need to copy the self.regexes container
            # note: this probably is not what you meant to do, since RegexSequence.Duplicate objects will not be copied
            return self.__class__(copy(self.regexes), self.repeats, self.min_repeats, 0, self.exactly, self.greedy)

        ###########
        # Methods #
        ###########

        def match(self, string):
            return self._search_or_match(string, re.match)

        def search(self, string):
            return self._search_or_match(string, re.search)

        def reset(self):
            self._counter = self.repeats
            self._spot = 0
            self._was_okay_to_finish = False

        def okay_to_finish(self):
            if not self._spot % len(self.regexes) == 0:
                self._was_okay_to_finish = False
                return False
            elif self._spot / len(self.regexes) < self.min_repeats:
                self._was_okay_to_finish = False
                return False
            elif self.repeats is None:
                self._was_okay_to_finish = True
                return True
            elif self.exactly and self._counter != 0:
                self._was_okay_to_finish = False
                return False
            else:
                if self._counter >= 0:
                    self._was_okay_to_finish = True
                    return True
                else:
                    self._was_okay_to_finish = False
                    return False

        def in_reset_state(self):
            """ True if reset() would change nothing
            """
            return self._counter == self.repeats and self._spot == 0 and self._was_okay_to_finish is False

        ###################
        # Private Methods #
        ###################

        def _search_or_match(self, string, func):
            if func is re.search:
                m = self.regexes[self._spot % len(self.regexes)].search(string)
            elif func is re.match:
                m = self.regexes[self._spot % len(self.regexes)].match(string)
            else:
                m = func(self.regexes[self._spot % len(self.regexes)], string)
            if m:
                self._spot += 1
                if self._counter is not None:
                    self._counter -= 1
                    if self._counter < 0:
                        self.reset()
                        return None
                    else:
                        return m
                else:
                    return m
            else:
                self.reset()
                return None


    class Duplicated(object):
        """ A mechanism for allowing a certain part of a pattern to be duplicated an arbitrary number of times in a line.

        """

        class Interlude(object):
            """ A super-lightweight class to hold an interlude """
            regex = None
            frequency = None
            max_repeats = None
            offset = None
            _repeats = 0

            def __init__(self, regex, frequency, max_repeats=None, offset=None):
                self.regex, self.frequency, self.max_repeats, self.offset = regex, frequency, max_repeats, offset
                self.offset = self.offset or 0
                self.reset()

            def __copy__(self):
                # basically we just need to create another instance which will allow _repeats to reset.
                return self.__class__(self.regex, self.frequency, self.max_repeats, self.offset)

            def reset(self):
                self._repeats = 0

        ##############
        # Attributes #
        ##############

        duplicated = None
        prelude = None
        postlude = None
        interludes = None
        min_dups = None
        max_dups = None
        flags = None

        ##################
        # Initialization #
        ##################

        def __init__(self, duplicated, prelude=None, postlude=None, interludes=None, min_dups=None, max_dups=None, flags=None):
            """
            """
            if isinstance(duplicated, basestring):
                self.duplicated = duplicated
            elif isinstance(duplicated, Regex):
                # "Decompile" the pattern, since we need to be able to concatinate it
                self.duplicated = duplicated.pattern
                if duplicated.flags != 0:
                    warn("Flags given to argument 'duplicated' of RegexSequence.Duplicated constructor will be ignored")
            else:
                raise ArgTypeError('duplicated', 'string', 'compiled regular expression', caller='RegexSequence.Duplicated constructor')
            #--------------------------------------------------------------------------------#
            if prelude is None:
                self.prelude = ''
            elif isinstance(prelude, basestring):
                self.prelude = prelude
            elif isinstance(prelude, Regex):
                # "Decompile" the pattern, since we need to be able to concatinate it
                self.prelude = prelude.pattern
                if prelude.flags != 0:
                    warn("Flags given to argument 'prelude' of RegexSequence.Duplicated constructor will be ignored")
            else:
                raise ArgTypeError('prelude', 'string', 'compiled regular expression', caller='RegexSequence.Duplicated constructor')
            #--------------------------------------------------------------------------------#
            if postlude is None:
                self.postlude = ''
            elif isinstance(postlude, basestring):
                self.postlude = postlude
            elif isinstance(postlude, Regex):
                # "Decompile" the pattern, since we need to be able to concatinate it
                self.postlude = postlude.pattern
                if postlude.flags != 0:
                    warn("Flags given to argument 'postlude' of RegexSequence.Duplicated constructor will be ignored")
            else:
                raise ArgTypeError('postlude', 'string', 'compiled regular expression', caller='RegexSequence.Duplicated constructor')
            #--------------------------------------------------------------------------------#
            # Interlude parsing helper functions
            def _invalid_interlude():
                raise TypeError("Interludes defined in RepeatedReges.Duplicated constructor must be of the form (str_or_regex, repeat_frequency[, max_repeats_or_none[, offset]])")
            def _check_interlude_seq(seq):
                # must have between 2 and 4 elements
                if not (2 <= len(interludes) <= 4): _invalid_interlude()
                # first entry must be a string
                if isinstance(interludes[0], (basestring, Regex)): _invalid_interlude()
                # second entry must be an int
                if not isinstance(interludes[1], int): _invalid_interlude()
                # if third entry is present, must be an int
                if len(interludes) > 2 and not isinstance(interludes[2], int): _invalid_interlude()
                # if fourth entry is present, must be an int
                if len(interludes) > 3 and not isinstance(interludes[3], int): _invalid_interlude()
            #--------------------------------------------------------------------------------#
            # Now check and assemble the 'interludes', things that can go between the main repeated parts
            if interludes is not None:
                if isinstance(interludes, Sequence) and not isinstance(interludes, basestring):
                    # First check to see if we have 0, 1 or multiple interludes
                    if len(interludes) == 0:
                        # If we have no interludes, this is the same as None
                        self.interludes = None
                    elif isinstance(interludes[0], (basestring, Regex)):
                        # this means we should only have 1 interlude...
                        # make sure it's valid
                        _check_interlude_seq(interludes)
                        # It's valid, let's add it
                        self.interludes = []
                        # We need to "fix" the pattern if it's a Regex object to be a string so we can concatinate it later.
                        fixed=[]
                        if isinstance(interludes[0], Regex):
                            fixed[0] = interludes[0].pattern
                            if interludes[0].flags != 0:
                                warn("Flags given to iterlude regex in RegexSequence.Duplicated constructor will be ignored")
                        else:
                            fixed[0] = interludes[0]
                        # everything else should be okay, just append
                        fixed.extend(interludes[1:])
                        # Now that we've fixed the pattern part, add it to the interludes...
                        self.interludes.append(RegexSequence.Duplicated.Interlude(*tuple(fixed)))
                    else:
                        # We have more than one interlude
                        self.interludes = []
                        for inter in interludes:
                            if isinstance(inter, Duplicated.Interlude):
                                # Mostly for calls from copy, but also could come from a really informed user...
                                self.interludes.append(inter)
                            elif isinstance(inter, Sequence):
                                # It's a tuple definition of an interlude; we need to unpack it into a Interlude object
                                # make sure it's valid
                                _check_interlude_seq(inter)
                                # It's valid, let's add it
                                # We need to "fix" the pattern if it's a Regex object to be a string so we can concatinate it later.
                                fixed=[]
                                if isinstance(inter[0], Regex):
                                    fixed[0] = inter[0].pattern
                                    # Check to make sure no flags are set, and warn the user that they will be ignored if this is the case...
                                    if inter[0].flags != 0:
                                        warn("Flags given to iterlude regex in RegexSequence.Duplicated constructor will be ignored")
                                else:
                                    fixed[0] = inter[0]
                                # everything else should be okay
                                fixed.extend(inter[1:])
                                # Now that we've fixed the pattern part, add it to the interludes...
                                self.interludes.append(RegexSequence.Duplicated.Interlude(*tuple(fixed)))
                            else:
                                # it's not an Interlude or an Sequence, so it's invalid
                                _invalid_interlude()
                elif isinstance(interludes, Duplicated.Interlude):
                    # wow...this user must really know what they are doing...
                    self.interludes = [interludes]
                else:
                    _invalid_interlude()
            #--------------------------------------------------------------------------------#
            self.min_dups = min_dups or 1
            self.max_dups = max_dups
            self.flags = flags or 0

        ##############
        # Properties #
        ##############

        @property
        def pattern(self):
            """ Returns the main regex, in order to conform with the `Regex` class protocol
            """
            return self.duplicated

        ###################
        # Special Methods #
        ###################

        def __copy__(self):
            # everything here is immutable except for interlude objects
            # technically, this is probably more of a deepcopy, but it doesn't make sense to do a shallow copy, so don't allow it
            # especially since the Interlude class isn't really supposed to be exposed to the user...
            new_interludes = [copy(inter) for inter in self.interludes] if self.interludes is not None else None
            return self.__class__(self.duplicated, self.prelude, self.postlude, new_interludes, self.min_dups, self.max_dups, self.flags)

        def __deepcopy__(self, memo):
            # same as __copy__, see comment in __copy__
            return copy(self)

        ###########
        # Methods #
        ###########

        def search(self, string):
            return self._search_or_match(string, re.search)

        def match(self, string):
            return self._search_or_match(string, re.match)

        ###################
        # Private Methods #
        ###################

        def _search_or_match(self, string, match_func):
            def _build_regex(dups, starred=False):
                if self.interludes is None:
                    main_part = self.duplicated * n_dups
                else:
                    main_part = self.duplicated
                    for inter in self.interludes: inter.reset()
                    for i in xrange(1, n_dups):
                        for inter in self.interludes:
                            if i >= inter.offset:
                                if (i + inter.offset) % inter.frequency == 0:
                                    inter._repeats += 1
                                    if inter.max_repeats is None or inter._repeats <= inter.max_repeats:
                                        main_part += inter.regex
                        main_part += self.duplicated
                return self.prelude + main_part + ('.*' if starred else '') + self.postlude
            def _do_match(rexp, str):
                return match_func(curr_regex, string, self.flags)
            n_dups = self.min_dups
            curr_regex = _build_regex(n_dups, True)
            m = _do_match(curr_regex, str)
            # First, search upward for the largest possible match...
            while m is not None and (self.max_dups is None or n_dups < self.max_dups):
                n_dups += 1
                curr_regex = _build_regex(n_dups, True)
                m = _do_match(curr_regex, str)
            # We've found the most that there can be, now try without the '.*' part...
            curr_regex = _build_regex(n_dups)
            m = _do_match(curr_regex, str)
            # Now scan downwards
            while m is None and n_dups > self.min_dups:
                n_dups -= 1
                curr_regex = _build_regex(n_dups)
                m = _do_match(curr_regex, str)
            return m


    ##############
    # Attributes #
    ##############

    main = None
    start = None
    stop = None
    last_only = None
    require_stop = None
    match_function = None

    ######################
    # Private Attributes #
    ######################

    _start_spot = 0
    _stop_spot = 0
    _main_spot = 0
    _active = False
    _curr_match_set = None
    _matches = None
    _last_match_okay = None

    ##################
    # Initialization #
    ##################

    # TODO multiline option that collects lines between start and stop and performs matching only once
    # TODO allow callables to be passed in
    # TODO distribute flags across repeated regexes (by adding flags argument to _add_repeated)
    # TODO !!! Greedy RepeatedSequences (currently they aren't greedy, and greediness seems to be the expected behavior)
    # TODO first_only kwarg
    # NOTE!!! make sure if you change this/add options to add them to __copy__ and __deepcopy__
    def __init__(self, main, main_flags=0,
                 start=None, start_flags=0,
                 stop=None, stop_flags=0,
                 last_only=False,
                 require_stop=False,
                 match_function=re.search):
        self.last_only = last_only
        self.require_stop = require_stop
        if self.require_stop: raise NotImplementedError("required presenece of the stop sequence is not yet implemented.")
        self.match_function = match_function
        # makes sure last part of  fixes the behavior of
        # repeated regexes with a given length to work
        def _add_repeated(ary, size, rexp, aryname):
            if len(ary) == 0:  # check if we're at the beginning
                ary.append(rexp)
            elif len(ary) == size - 1:  #check if we're at the end
                if aryname == 'start' and rexp.greedy:
                    raise NotImplementedError("Start sequences ending in greedy RepeatedRegexes are not implemented.")
                ary.append(rexp)
            else: # we're in the middle.  Check for an identical regex, and raise an exception because that's probably not what the user wanted to do
                if isinstance(ary[-1], RegexSequence.RepeatedRegex) and (rexp.regexes, rexp.exactly) == (ary[-1].regexes, ary[-1].exactly):
                    raise TypeError("Multiple identical RepeatedRegex objects given consecutively."
                                    "  This is probably not what you meant to do; if it is, please combine them.")
                else:
                    ary.append(rexp)
        #--------------------------------------------------------------------------------#
        # Get the main regex or regex sequence
        if isinstance(main, basestring):
            self.main = re.compile(main, main_flags)
        elif isinstance(main, RepeatedRegex):
            # always handle RepeatedRegex objects as a part of a list, even if there is only one
            self.main = [main]
        elif hasattr(main, 'match'):
            self.main = main
        elif isinstance(main, Sequence):
            self.main = []
            if isinstance(main_flags, Sequence):
                flags = main_flags
                if len(main) != len(main_flags):
                    raise TypeError("Dimension mismatch in RegexSequence constructor: main_flags must either be an integer or a list of the same length as main.")
            else:
                flags = [main_flags] * len(main)
            for regex, flgs in zip(main, flags):
                if isinstance(regex, basestring):
                    self.main.append(re.compile(regex, flgs))
                elif isinstance(regex, RegexSequence.RepeatedRegex):
                    _add_repeated(self.main, len(main), regex, 'main')
                elif hasattr(regex, 'match'):
                    self.main.append(regex)
                elif regex is Ellipsis:
                    # This only makes sense in Python 3, where the ellipsis literal '...' can be used outside of the context of a slice
                    _add_repeated(self.main, len(main), RegexSequence.RepeatedRegex(r''), 'main')
                else:
                    raise TypeError("Invalid item given in regex list 'main' to RegexSequence constructor: " + str(regex))
        else:
            raise TypeError("Invalid type for regex given for argument 'main' of RegexSequence constructor: " + str(main))
        #--------------------------------------------------------------------------------#
        # Get the start regex or regex sequence
        if not start is None:
            if isinstance(start, basestring):
                self.start = re.compile(start, start_flags)
            elif isinstance(start, RepeatedRegex):
                # always handle RepeatedRegex objects as a part of a list, even if there is only one
                self.start = [start]
            elif hasattr(start, 'match'):
                self.start = start
            elif isinstance(start, Sequence):
                self.start = []
                if isinstance(start_flags, Sequence):
                    flags = start_flags
                    if len(start) != len(start_flags):
                        raise TypeError("Dimension mismatch in RegexSequence constructor: start_flags must either be an integer or a list of the same length as start.")
                else:
                    flags = [start_flags] * len(start)
                for regex, flgs in zip(start, flags):
                    if isinstance(regex, basestring):
                        self.start.append(re.compile(regex, flgs))
                    elif isinstance(regex, RegexSequence.RepeatedRegex):
                        _add_repeated(self.start, len(start), regex, 'start')
                    elif hasattr(regex, 'match'):
                        self.start.append(regex)
                    elif regex is Ellipsis:
                        # This only makes sense in Python 3, where the ellipsis literal '...' can be used outside of the context of a slice
                        _add_repeated(self.start, len(start), RegexSequence.RepeatedRegex(r''), 'start')
                    else:
                        raise TypeError("Invalid item given in regex list 'start' to RegexSequence constructor: " + str(regex))
            else:
                raise TypeError("Invalid type for regex given for argument 'start' of RegexSequence constructor: " + str(start))
        else:
            self.start = None
        #--------------------------------------------------------------------------------#
        # Get the stop regex or regex sequence
        if not stop is None:
            if isinstance(stop, basestring):
                self.stop = re.compile(stop, stop_flags)
            elif isinstance(stop, RepeatedRegex):
                # always handle RepeatedRegex objects as a part of a list, even if there is only one
                self.stop = [stop]
            elif hasattr(stop, 'match'):
                self.stop = stop
            elif isinstance(stop, Sequence):
                self.stop = []
                if isinstance(stop_flags, Sequence):
                    flags = stop_flags
                    if len(stop) != len(stop_flags):
                        raise TypeError("Dimension mismatch in RegexSequence constructor: stop_flags must either be an integer or a list of the same length as stop.")
                else:
                    flags = [stop_flags] * len(stop)
                for regex, flgs in zip(stop, flags):
                    if isinstance(regex, basestring):
                        self.stop.append(re.compile(regex, flgs))
                    elif isinstance(regex, RegexSequence.RepeatedRegex):
                        _add_repeated(self.stop, len(stop), regex, 'stop')
                    elif hasattr(regex, 'match'):
                        self.stop.append(regex)
                    elif regex is Ellipsis:
                        # This only makes sense in Python 3, where the ellipsis literal '...' can be used outside of the context of a slice
                        _add_repeated(self.stop, len(stop), RegexSequence.RepeatedRegex(r''), 'stop')
                    else:
                        raise TypeError("Invalid item given in regex list 'stop' to RegexSequence constructor: " + str(regex))
            else:
                raise TypeError("Invalid type for regex given for argument 'stop' of RegexSequence constructor: " + str(stop))
        else:
            self.stop = None
        #--------------------------------------------------------------------------------#
        self._check_for_edge_cases()
        self.reset()

    ##############
    # Properties #
    ##############

    @property
    def matches(self):
        """ A (possibly nested) list of the matches collected thus far.
        """
        if self._last_match_okay:
            # This only applies to cases where we have a sequence at the end
            # thus, we should expect a nexted list in matches...
            if not self.last_only:
                return self._matches + [self._curr_match_set]
            else:
                return self._curr_match_set
        else:
            return self._matches


    ###################
    # Special Methods #
    ###################

    def __copy__(self):
        # Just shallow copy all of the arrays.  It rarely makes sense to do this; you should almost certainly be
        # using deepcopy.
        warn("Calling copy with a RegexSequence object doesn't make much sense, and this is a really advanced feature.  Are you sure you didn't mean to use 'deepcopy'?")
        # note that flags don't need to be copied since they are already compiled in to the compiled regular expressions
        return self.__class__(
            main = copy(self.main) if isinstance(self.main, Sequence) else self.main, main_flags = 0,
            start = copy(self.start) if isinstance(self.start, Sequence) else self.start, start_flags = 0,
            stop = copy(self.stop) if isinstance(self.stop, Sequence) else self.stop, stop_flags = 0,
            last_only=self.last_only, require_stop=self.require_stop, match_function=self.match_function
        )

    def __deepcopy__(self, memo):
        if isinstance(self.main, Sequence):
            new_main = [deepcopy(item, memo) if not isinstance(item, CompiledRegexType) else item for item in self.main]
        else:
            # Single regular expression, which is immutable so don't copy it
            new_main = self.main
        # do the same thing for start:
        if isinstance(self.start, Sequence):
            new_start = [deepcopy(item, memo) if not isinstance(item, CompiledRegexType) else item for item in self.start]
        else:
            new_start = self.start
        # and for stop:
        if isinstance(self.stop, Sequence):
            new_stop = [deepcopy(item, memo) if not isinstance(item, CompiledRegexType) else item for item in self.stop]
        else:
            new_stop = self.stop
        # now construct the new instance and return
        # note that flags don't need to be copied since they are already compiled in to the compiled regular expressions
        return self.__class__(
            main=new_main, main_flags=0,
            start=new_start, start_flags=0,
            stop=new_stop, stop_flags=0,
            last_only=self.last_only,
            require_stop=self.require_stop,
            match_function=self.match_function
        )

    ###########
    # Methods #
    ###########

    def reset(self):
        """ Reset the regex sequence.
        Do this, for instance, when you've finished parsing one file and need to start parsing another.
        This gets called by the constructor.  Most internal uses of this class should already do this
        when it is supposed to be done.
        """
        self._start_spot = 0
        self._stop_spot = 0
        self._main_spot = 0
        self._matches = []
        self._curr_match_set = []
        self._active = False
        self._last_match_okay = False
        if self.start is None:
            self._active = True
        for items in [self.start, self.main, self.stop]:
            if isinstance(items, Iterable):
                for item in items:
                    if hasattr(item, 'reset'):
                        item.reset()

    def parse_line(self, line):
        """ Parse a line from a file.
        """
        # Handles greedy matches at the end of the main sequence.
        # Irrelevant if
        self._last_match_okay = False
        if self._active:
            #-------------------------#
            # Check for stop sequence #
            #-------------------------#
            if self.stop is not None:
                if not isinstance(self.stop, Sequence):
                    if self._match(self.stop, line):
                        self._active = False
                        return
                else:  # We have a list, and our life becomes much more difficult
                    spot = self.stop[self._stop_spot]
                    if isinstance(spot, RegexSequence.RepeatedRegex):
                        # See if we're at the end...
                        if self._stop_spot + 1 == len(self.stop):
                            # See if we're okay to finish
                            if spot.okay_to_finish():
                                # Then finish...
                                self._stop_spot = 0
                                self._active = False
                                spot.reset()
                                # and we want to make sure we're not starting another start sequence right away, so now
                                # that we've reset everything, recurse
                                self.parse_line(line)
                                return
                            else:
                                # Try to keep matching
                                if self._match(spot, line):
                                    # Then we're okay to keep going...
                                    return
                                else:
                                    # we've failed.
                                    # if we were at spot 1, perhaps we just started too soon
                                    if self._stop_spot == 1:
                                        # go back and try spot 0 again
                                        self._stop_spot = 0
                                        spot.reset()
                                        self.parse_line(line)
                                    else:
                                        # reset self._stop_spot...
                                        self._stop_spot = 0
                                        spot.reset()
                        else:
                            # Not at the end
                            next = self.stop[self._stop_spot + 1]
                            if self._match(next, line):
                                if spot.okay_to_finish():
                                    # Move on to the next spot in the sequence
                                    self._stop_spot += 1
                                    # Make sure that the RepeatedRegex in spot is ready for action next time
                                    spot.reset()
                                    # Now increment self._stop_spot *again* unless the next regex is also a RepeatedRegex
                                    # (since matching next "gobbles up" the next part of the sequence)
                                    if not isinstance(next, RegexSequence.RepeatedRegex):
                                        self._stop_spot += 1
                                    else:
                                        # otherwise, we want to try again, since the next set of repeats are allowed to be 0 length...
                                        # first reset to avoid double counting
                                        next.reset()
                                        # and recurse...
                                        self.parse_line(line)
                                    # and check if we're done...
                                    if self._stop_spot == len(self.stop):
                                        # The whole stop sequence has been matched, so stop
                                        self._active = False
                                        self._stop_spot = 0
                                        return
                                else:
                                    # spot wasn't an okay repeated regex, and we moved were ready to move on
                                    # so the stop sequence fails
                                    spot.reset()
                                    # reset next if it's a RepeatedRegex, since we called match on it and it succeeded,
                                    # thus decrementing the _counter in the repeated regex
                                    if isinstance(next, RegexSequence.RepeatedRegex): next.reset()
                                    self._stop_spot = 0
                            elif self._match(spot, line):
                                # Don't increment, since the repeated sequence is okay...
                                # if spot has run out of repeats, self._match(spot, line) will return None even if the
                                #   regex itself matches
                                pass
                            else:
                                # spot doesn't match, and the line doesn't match the next part of the stop sequence,
                                # so the stop sequence fails
                                self._stop_spot = 0
                                spot.reset()
                    else:  # Not a RepeatedRegex
                        if self._match(spot, line):
                            # Move on to the next spot in the sequence
                            self._stop_spot += 1
                            if self._stop_spot == len(self.stop):
                                # The whole stop sequence has been matched, so stop
                                self._active = False
                                self._stop_spot = 0
                                return
                        else:
                            # Match failure, reset stop spot...
                            self._stop_spot = 0
            #--------------------------------------------------------------------------------#
            #----------------------------------#
            # Check main sequence, get matches #
            #----------------------------------#
            if not isinstance(self.main, Sequence):
                m = self._match(self.main, line)
                if m:
                    if not self.last_only:
                        self._matches.append(m)
                    else:
                        self._matches = [m]
                    return
            else:  # We have a list, and our life becomes much more difficult
                spot = self.main[self._main_spot]
                if isinstance(spot, RegexSequence.RepeatedRegex):
                    if self._main_spot + 1 == len(self.main):
                        if spot.okay_to_finish():
                            # Note that greedy doesn't apply to exact sequences
                            if not spot.greedy or spot.exactly:
                                # We've been cleared to finish.  This means we can reset spot
                                spot.reset()
                                if not self.last_only:
                                    self._matches.append(self._curr_match_set)
                                else:
                                    self._matches = self._curr_match_set
                                self._curr_match_set = []
                                self._main_spot = 0
                                # and we want to make sure we're not starting another sequence right away, so now
                                # that we've reset everything, recurse
                                self.parse_line(line)
                            else:
                                # Try to keep matching
                                self._last_match_okay = True
                                was_okay = spot._was_okay_to_finish
                                was_in_reset_state = spot.in_reset_state()
                                match = self._match(spot, line)
                                if match:
                                    # Then we're okay to keep going...
                                    self._curr_match_set.append(match)
                                elif was_okay:
                                    # It was okay to finish last time, but we were greedy.  That's okay, finish now...
                                    #TODO wrap-up after last line...
                                    # We've been cleared to finish.  This means we can reset spot
                                    spot.reset()
                                    # And add matches...
                                    if not self.last_only:
                                        self._matches.append(self._curr_match_set)
                                    else:
                                        self._matches = self._curr_match_set
                                    self._curr_match_set = []
                                    self._main_spot = 0
                                    # and we want to make sure we're not starting another sequence right away, so now
                                    # that we've reset everything, recurse
                                    self.parse_line(line)
                                else:
                                    # All is lost.  We came to an unmatched line and we weren't okay to finish previously.
                                    # If this was the first regex (and also the last), maybe we started too soon
                                    # recurse only once then...
                                    if self._main_spot == 0 and not was_in_reset_state:
                                        # reset things...
                                        spot.reset()
                                        self._main_spot = 0
                                        self._curr_match_set = []
                                        # and try again with spot reset.  Perhaps the first Regex in spot will work
                                        # if spot is a multi-regex repeated sequence
                                        self.parse_line(line)
                                    else:
                                        spot.reset()
                                        self._main_spot = 0
                                        self._curr_match_set = []
                                        # if we were at spot 1, perhaps we just started too soon
                                        if self._main_spot == 1:
                                            # go back and try spot 0 again
                                            self.parse_line(line)
                        else:
                            # Try to keep matching
                            match = self._match(spot, line)
                            if match:
                                # Then we're okay to keep going...
                                self._curr_match_set.append(match)
                                # check to see if *now* we are okay to finish...
                                if spot.okay_to_finish():
                                    # if we don't get run again, we still want to pick up the last match...
                                    self._last_match_okay = True
                            else:
                                # we've failed.
                                # if we were at spot 1, perhaps we just started too soon
                                if self._main_spot == 1:
                                    # go back and try spot 0 again
                                    self._main_spot = 0
                                    spot.reset()
                                    self.parse_line(line)
                                else:
                                    # reset self._stop_spot...
                                    self._main_spot = 0
                                    spot.reset()
                    else:
                        # We're not at the end; try again the next one
                        next = self.main[self._main_spot + 1]
                        nextmatch = self._match(next, line)
                        if nextmatch:
                            # now see if we are okay to finish...
                            if spot.okay_to_finish():
                                # Move on to the next spot in the sequence
                                self._main_spot += 1
                                # Make sure that the RepeatedRegex in spot is ready for action next time
                                spot.reset()
                                # Now increment self._main_spot *again* unless the next regex is also a RepeatedRegex
                                if not isinstance(next, RegexSequence.RepeatedRegex):
                                    self._main_spot += 1
                                    # And add the match from next to the current match set...
                                    self._curr_match_set.append(nextmatch)
                                else:
                                    # otherwise, we want to try again, since the next set of repeats are allowed to be 0 length...
                                    # first reset to avoid double counting
                                    next.reset()
                                    # don't add the match; it will get added next time if the following sequence doesn't match
                                    # and recurse...
                                    self.parse_line(line)
                                # and check if we're done...
                                if self._main_spot == len(self.main):
                                    # The whole main sequence has been matched, so add the match set to the matches and start over
                                    if not self.last_only:
                                        self._matches.append(self._curr_match_set)
                                    else:
                                        self._matches = [self._curr_match_set]
                                    self._curr_match_set = []
                                    self._main_spot = 0
                            else:
                                # spot wasn't an okay repeated regex, and we moved were ready to move on
                                # so the main sequence fails
                                spot.reset()
                                # reset next if it's a RepeatedRegex, since we called match on it and it succeeded,
                                # thus decrementing the _counter in the repeated regex
                                if isinstance(next, RegexSequence.RepeatedRegex): next.reset()
                                self._main_spot = 0
                                # empty the main match set; these weren't the droids we were looking for
                                self._curr_match_set = []
                        else:
                            spotmatch = self._match(spot, line)
                            if spotmatch:
                                # Don't increment, since the repeated sequence is okay...
                                # if spot has run out of repeats, self._match(spot, line) will return None even if the
                                #   regex itself matches
                                self._curr_match_set.append(spotmatch)
                            else:
                                # spot doesn't match, and the line doesn't match the next part of the main sequence,
                                # so the main sequence match fails
                                self._main_spot = 0
                                spot.reset()
                                # empty the main match set; this wasn't the droids we were looking for
                                self._curr_match_set = []
                else:  # Not a RepeatedRegex
                    match = self._match(spot, line)
                    if match:
                        # add the match to the current match set
                        self._curr_match_set.append(match)
                        # and move on to the next spot in the sequence
                        self._main_spot += 1
                        if self._main_spot == len(self.main):
                            # The whole main sequence has been matched, so add the match set to the matches and start over
                            if not self.last_only:
                                self._matches.append(self._curr_match_set)
                            else:
                                self._matches = self._curr_match_set
                            self._curr_match_set = []
                            self._main_spot = 0
                    else:
                        # If we were at spot 1, perhaps we just started too soon.  See if everything will be okay if we just restart.
                        if self._main_spot == 1:
                            self._main_spot = 0
                            self._curr_match_set = []
                            self.parse_line(line)
                        # Match failure, reset main spot...
                        self._main_spot = 0
                        # and flush the current match set...
                        self._curr_match_set = []
        #--------------------------------------------------------------------------------#
        else:
            #---------------------------------------------#
            # Check start sequence, activate if necessary #
            #---------------------------------------------#
            if self.start is not None:
                if not isinstance(self.start, Sequence):
                    if self._match(self.start, line):
                        self._active = True
                        return
                else:  # We have a list, and our life becomes much more difficult
                    spot = self.start[self._start_spot]
                    if isinstance(spot, RegexSequence.RepeatedRegex):
                        # See if we're at the end...
                        if self._start_spot + 1 == len(self.start):
                            # See if we're okay to finish
                            if spot.okay_to_finish():
                                # Then finish...
                                self._start_spot = 0
                                self._active = True
                                spot.reset()
                                # and we want to make sure we're not starting the main sequence right away, so now
                                # that we've reset everything, recurse
                                self.parse_line(line)
                                return
                            else:
                                # Try to keep matching
                                if self._match(spot, line):
                                    # Then we're okay to keep going...
                                    return
                                else:
                                    # we've failed.
                                    # if we were at spot 1, perhaps we just started too soon
                                    if self._start_spot == 1:
                                        # go back and try spot 0 again
                                        self._start_spot = 0
                                        spot.reset()
                                        self.parse_line(line)
                                    else:
                                        # reset self._stop_spot...
                                        self._start_spot = 0
                                        spot.reset()
                        else:
                            # We're not at the end, we can go ahead and check the next spot
                            next = self.start[self._start_spot + 1]
                            if self._match(next, line):
                                if spot.okay_to_finish():
                                    # Move on to the next spot in the sequence
                                    self._start_spot += 1
                                    # Make sure that the RepeatedRegex in spot is ready for action next time
                                    spot.reset()
                                    # Now increment self._start_spot *again* unless the next regex is also a RepeatedRegex
                                    # (since matching next "gobbles up" the next part of the sequence)
                                    if not isinstance(next, RegexSequence.RepeatedRegex):
                                        self._start_spot += 1
                                        # and check if we're done...
                                    else:
                                        # otherwise, we want to try again, since the next set of repeats are allowed to be 0 length...
                                        # first reset to avoid double counting
                                        next.reset()
                                        # and recurse...
                                        self.parse_line(line)
                                    if self._start_spot == len(self.start):
                                        # The whole start sequence has been matched, so start
                                        self._active = True
                                        self._start_spot = 0
                                        return
                                else:
                                    # spot wasn't an okay repeated regex, and we're ready to move on
                                    # so the start sequence fails
                                    spot.reset()
                                    # reset next if it's a RepeatedRegex, since we called match on it and it succeeded,
                                    # thus decrementing the _counter in the repeated regex
                                    if isinstance(next, RegexSequence.RepeatedRegex): next.reset()
                                    self._start_spot = 0
                            elif self._match(spot, line):
                                # Don't increment, since the repeated sequence is okay...
                                # if spot has run out of repeats, self._match(spot, line) will return None even if the
                                #   regex itself matches
                                pass
                            else:
                                # spot doesn't match, and the line doesn't match the next part of the start sequence,
                                # so the start sequence fails
                                self._start_spot = 0
                                spot.reset()
                    else:  # Not a RepeatedRegex
                        if self._match(spot, line):
                            # Move on to the next spot in the sequence
                            self._start_spot += 1
                            if self._start_spot == len(self.start):
                                # The whole start sequence has been matched, so start
                                self._active = True
                                self._start_spot = 0
                                return
                        else:
                            # Match failure, reset start spot...
                            self._start_spot = 0
        #--------------------------------------------------------------------------------#
        return

    def groups(self, flat=True):
        """ Returns a tuple containing all the match groups from all the matches.
        If the `flat` argument is True, (default) the groups are flattened across all matches.  If flat is False,
        the groups are only flattened within a given series (this makes no difference for cases where the `main`
        regex is a single expression rather than a list of regexes; in this case a flat list of groups is always
        returned).

        The name "groups" is meant to correspond roughly with the naming convention in Python's re module.
        """
        if not isinstance(self.main, Iterable):
            return tuple(chain.from_iterable(m.groups() for m in self.matches))
        else:
            if flat:
                return tuple(chain.from_iterable(chain.from_iterable(m.groups() for m in seq) for seq in self.matches))
            else:
                return tuple(tuple(chain.from_iterable(m.groups() for m in seq)) for seq in self.matches)

    ###################
    # Private Methods #
    ###################

    def _match(self, regex_like, line):
        if self.match_function == re.match:
            return regex_like.match(line)
        elif self.match_function == re.search:
            return regex_like.search(line)
        else:
            #assume its a callable that takes a regex_like thing (compiled regex, Duplicated, or Repeated) and a line...
            # for super-advanced use only
            # TODO raise an error if self.match_function isn't a callable
            return self.match_function(regex_like, line)

    def _check_for_edge_cases(self):
        # Check the main sequence...
        # If only one item is given, and it is an exact RepeatedRegex, things won't work.  You can always do this
        # another way anyhow...
        if isinstance(self.main, Sized) and len(self.main) == 1 \
                and isinstance(self.main[0], RegexSequence.RepeatedRegex) and self.main[0].exactly:
            raise ValueError('a single, exactly repeated regex in the main sequence is ambiguous.\n'
                'Consider the set of lines:\n    a1\n    a2\n    a3\nwith the repeated regex: Repeated(r\'a(\\d+)\', repeats=2, exactly=True)\n'
                'Technically, there are three repeats, so the RepeatedRegex should not match.  But there are actually\n'
                'two different cases that match: the first two lines and the last two lines.  To avoid any confusion,\n'
                'the use of this form has been prohibited.  Please try something else (this exact RepeatedRegex can always be \n'
                'written out explicitly if this is really what you want to do)'
            )

# Convenient Aliases
Duplicated = RegexSequence.Duplicated
Repeated = RepeatedRegex = RegexSequence.RepeatedRegex

class MatrixRegexSequence(RegexSequence):
    """ Special subclass of `RegexSequence` that make parsing matrices easier.
    """

    ##############
    # Attributes #
    ##############

    transpose = None

    ##################
    # Initialization #
    ##################

    @with_flexible_arguments(
        optional = [
            ('transpose',),
        ]
    )
    def __init__(self, main=Repeated(Duplicated('(-?\d+\.\d+)\s*')), *args, **kwargs):
        self.transpose = pop_kwarg(kwargs, 'transpose') or False
        super(MatrixRegexSequence, self).__init__(main, *args, **kwargs)

    ###################
    # Special Methods #
    ###################

    def __deepcopy__(self, memo):
        ret_val = super(MatrixRegexSequence, self).__deepcopy__(memo)
        ret_val.transpose = self.transpose
        return ret_val

    def __copy__(self):
        ret_val = super(MatrixRegexSequence, self).__copy__()
        ret_val.transpose = self.transpose
        return ret_val

    ###########
    # Methods #
    ###########

    def get_matrix(self):
        # Assumes all of the rows are given together.  If not, things get to be more difficult...
        rows = []
        row_num = 0
        first_time = True
        for set in self.matches:
            for match in set:
                if first_time:
                    rows.append([])
                for item in match.groups():
                    rows[row_num].append(float(item))
                row_num += 1
            first_time = False
            row_num = 0
        if self.transpose:
            return Matrix(rows).T
        else:
            return Matrix(rows)

MatrixRegex = MatrixRegexSequence

#####################
# Dependent Imports #
#####################

from .exceptions import ArgTypeError

