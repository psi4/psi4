import unittest

import sys
import os

# Add the directory containing the 'grendel_tests' package to sys.path
#sys.path.insert(0, os.path.join(os.path.abspath(os.path.dirname(__file__)), os.pardir, os.pardir, os.pardir, "lib", "python"))
#sys.path.insert(0, "/Users/dhollman/Projects/PyGrendel")

from grendel.util.metaclasses import Immutable, mutator


class Foo(object):
    __metaclass__ = Immutable

    def __init__(self, a, b, c):
        self.a, self.b, self.c = a, b, c
        self.change_a(a)
        self.abc = a, b, c

    def __str__(self):
        return "Foo" + str((self.a, self.b, self.c))
    __repr__ = __str__

    @mutator
    def change_a(self, new_val):
        self.a = new_val

    @property
    def abc(self):
        return self.a, self.b, self.c

    @abc.setter
    def abc(self, val):
        self.a, self.b, self.c = val

class ImmutableMetaclassTest(unittest.TestCase):

    def test_construct_1(self):
        t = Foo(1, 'a', (123, "abc"))
        self.assertEqual(str(t), "Foo(1, 'a', (123, 'abc'))")
        self.assertEqual(t.b, 'a')
        self.assertEqual(t.a, 1)
        self.assertEqual(t.c, (123, "abc"))

    def test_construct_2(self):
        t = Foo(1, Foo(1, 2, 3), (123, "abc", Foo('a', 2, 4)))
        self.assertEqual(str(t), "Foo(1, Foo(1, 2, 3), (123, 'abc', Foo('a', 2, 4)))")

    def test_construct_3(self):
        with self.assertRaises(AttributeError):
            Foo(["hello", "world"], 'a', (123, "abc"))

    def test_no_set_1(self):
        t = Foo(1, 'a', (123, "abc"))
        with self.assertRaisesRegexp(TypeError, "Can't change attribute 'a' of object of immutable type 'Foo'"):
            t.a = 6

    def test_no_set_2(self):
        t = Foo(1, 'a', (123, "abc"))
        with self.assertRaisesRegexp(TypeError, "Can't change attribute 'a' of object of immutable type 'Foo'"):
            setattr(t, 'a', 6)

    def test_no_set_3(self):
        t = Foo(1, 'a', (123, "abc"))
        with self.assertRaisesRegexp(TypeError, "Can't change attribute 'a' of object of immutable type 'Foo'"):
            t.change_a(6)

    def test_property(self):
        t = Foo(1, 'a', (123, "abc"))
        self.assertEqual(t.abc, (1, 'a', (123, "abc")))

    def test_no_set_property(self):
        t = Foo(1, 'a', (123, "abc"))
        with self.assertRaisesRegexp(TypeError, "Can't change attribute 'abc' of object of immutable type 'Foo'"):
            t.abc = 1, 2, 3

    def test_subclass_1(self):
        with self.assertRaises(TypeError):
            class Bar(Foo):
                def __setattr__(self, key, value):
                    pass

    def test_subclass_2(self):
        with self.assertRaises(TypeError):
            class Bar(Foo):
                def __delattr__(self, item):
                    pass

    def test_subclass_3(self):
        with self.assertRaises(TypeError):
            class FooBar(object):
                pass
            class Bar(FooBar):
                __metaclass__ = Immutable
                pass

    def test_subclass_4(self):
        class Bar(Foo):
            __metaclass__ = Immutable
            def __init__(self, a, b, c, d):
                self.d = d
                self.change_a(a)
                super(Bar, self).__init__(a, b, c)
        bar = Bar(1, 2, 3, 4)
        self.assertEqual(str(bar), "Foo(1, 2, 3)")
