from __future__ import print_function
import sys
import os
import unittest

# Add the directory containing 'grendel' to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir))

# Add the directory containing the 'grendel_tests' package to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from grendel.util.overloading import overloaded, OverloadedFunctionCallError, PartiallyConstructed


class overloadedTest(unittest.TestCase):

    @overloaded
    def example(self, *args, **kwargs):
        raise OverloadedFunctionCallError

    @example.overload_with(
        number=(int, float, complex),
        index=int)
    def example(self, number, index):
        #print("Number and index: %d %d" % (int(number), index))
        return number + index

    @example.overload_with(
        word=basestring)
    def example(self, word):
        #print("%s is the word backwards!" % word[::-1])
        return word[::-1]

    def test_overload1(self):
        assert self.example(5.5, 10) == 15.5

    def test_overload2(self):
        ret_val = self.example("hello")
        assert ret_val == "olleh"


class partiallyConstructedTest(unittest.TestCase):

    class Custom1(object):
        def __init__(self, arg1, arg2, arg3, arg4='hello', arg5='world', *args, **kwargs):
            self.attr1, self.attr2, self.attr3 = 1, 2, 3
            self.string1 = ', '.join(str(i) for i in ((arg1, arg2, arg3, arg4, arg5) + args + tuple(sorted(kwargs.keys()))))
        def attr_string(self):
            return ', '.join(str(a) for a in [self.attr1, self.attr2, self.attr3])

    def test_custom1_1(self):
        p = PartiallyConstructed(partiallyConstructedTest.Custom1, 4, 5, some_kwarg='something')
        full = p(6)
        self.assertEqual(full.string1, '4, 5, 6, hello, world, some_kwarg')

    def test_custom1_2(self):
        p = PartiallyConstructed(partiallyConstructedTest.Custom1, 'a', g='something')
        p.with_attributes(attr1=4, attr2=5, attr3=6)
        full = p('b', 'c', 'd', 'e', 'f')
        self.assertEqual(full.string1, 'a, b, c, d, e, f, g')
        self.assertEqual(full.attr_string(), '4, 5, 6')

    def test_custom1_3(self):
        # Should be okay as long as keyword versions are used in both cases...
        p = PartiallyConstructed(
            partiallyConstructedTest.Custom1,
            arg1='a1',
            arg3='a3',
            another_kwarg='kwarg',
            arg5='a5',
        )
        full = p(arg2='a2', arg4='a4')
        self.assertEqual(full.string1, 'a1, a2, a3, a4, a5, another_kwarg')

    def test_custom_double(self):
        p = PartiallyConstructed(partiallyConstructedTest.Custom1, 'a', g='something')
        p.with_attributes(attr1='a1')
        p.with_args('b')
        q = PartiallyConstructed(p, arg5='e', h='blah')
        q.with_attributes(attr3='a3')
        q.with_kwargs(arg4='d', f="kwarg")
        q.with_args('c')
        full = q(i='something else')
        full.attr2='a2'
        self.assertEqual(full.string1, 'a, b, c, d, e, f, g, h, i')
        self.assertEqual(full.attr_string(), 'a1, a2, a3')

