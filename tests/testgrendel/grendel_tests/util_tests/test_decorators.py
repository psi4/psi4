import unittest

import sys
import os

# Add the directory containing the 'grendel_tests' package to sys.path
from grendel.util.metaclasses import Immutable

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from grendel.util.decorators import with_flexible_arguments, CachedMethod

class FlexibleArgsTest(unittest.TestCase):

    def test_flexible_args_trivial(self):
        @with_flexible_arguments()
        def flexible():
            return 'hello'
        self.assertEqual(flexible(), 'hello')
        with self.assertRaisesRegexp(TypeError, r"flexible.* takes no arguments \(1 given\)"):
            flexible('goodbye')

    def test_flexible_args_basic(self):
        @with_flexible_arguments(
            required=[
                ('test', 'test_arg')
            ]
        )
        def flexible(test):
            return test
        self.assertEqual(flexible('hello'), 'hello')
        self.assertEqual(flexible(test='hello'), 'hello')
        self.assertEqual(flexible(test_arg='hello'), 'hello')
        # Exceptions that can be raised...
        with self.assertRaisesRegexp(TypeError, r"flexible\(\) got multiple values for keyword argument 'test'"):
            flexible('hello', test='hello')
        with self.assertRaisesRegexp(TypeError, r"flexible\(\) got multiple values for keyword argument 'test'"):
            flexible('hello', test_arg='hello')
        with self.assertRaisesRegexp(TypeError, r"flexible\(\) got multiple values for keyword argument 'test'"):
            flexible(test_arg='hello', test='hello')
        with self.assertRaisesRegexp(TypeError, r"flexible.* missing required keyword argument 'test'"):
            flexible()
        with self.assertRaisesRegexp(TypeError, r"flexible.* got unexpected keyword argument 'not_an_arg'"):
            flexible('hello', not_an_arg='goodbye')
        with self.assertRaisesRegexp(TypeError, r"flexible.* takes exactly 1 argument "):
            flexible('hello', 'goodbye')

    def test_flexible_args_optional(self):
        @with_flexible_arguments(
            required=[
                ('test', 'test_arg')
            ],
            optional=[
                ('optional', 'foo', 'foobar')
            ]
        )
        def flexible(test, optional=''):
            return test + optional
        self.assertEqual(flexible(test_arg='hello', foobar=' world'), 'hello world')
        #make sure the rest of the tests still work...
        self.assertEqual(flexible('hello'), 'hello')
        self.assertEqual(flexible(test='hello'), 'hello')
        self.assertEqual(flexible(test_arg='hello'), 'hello')
        # A few exceptions that can be raised...
        with self.assertRaisesRegexp(TypeError, r"flexible\(\) got multiple values for keyword argument 'optional'"):
            flexible('hello', 'world', optional='hello')
        with self.assertRaisesRegexp(TypeError, r"flexible\(\) got multiple values for keyword argument 'optional'"):
            flexible('hello', optional='world', foo='test')

    def test_flexible_args_vargs(self):
        @with_flexible_arguments(
            required=[
                ('test', 'test_arg')
            ],
            optional=[
                ('optional', 'foo', 'foobar')
            ]
        )
        def flexible(test, optional='', *args):
            return test + optional + ''.join(args)
        self.assertEqual(flexible('hello', ' world', 'a', 'b', 'c'), 'hello worldabc')
        # Make sure the rest of the tests still work...
        self.assertEqual(flexible('hello'), 'hello')
        self.assertEqual(flexible(test='hello'), 'hello')
        self.assertEqual(flexible(test_arg='hello'), 'hello')
        self.assertEqual(flexible(test_arg='hello', foobar=' world'), 'hello world')

    def test_flexible_args_vargs_kwargs(self):
        @with_flexible_arguments(
            required=[
                ('test', 'test_arg')
            ],
            optional=[
                ('optional', 'foo', 'foobar')
            ]
        )
        def flexible(test, optional='', *args, **kwargs):
            return test + optional + ''.join(args) + ','.join(k+'='+kwargs[k] for k in kwargs)
        self.assertEqual(flexible('hello', ' world', 'a', 'b', 'c', jkl='def'), 'hello worldabcjkl=def')
        # Make sure the rest of the tests still work...
        self.assertEqual(flexible('hello'), 'hello')
        self.assertEqual(flexible(test='hello'), 'hello')
        self.assertEqual(flexible(test_arg='hello'), 'hello')
        self.assertEqual(flexible(test_arg='hello', foobar=' world'), 'hello world')
        self.assertEqual(flexible('hello', ' world', 'a', 'b', 'c'), 'hello worldabc')
        # A few exceptions that can be raised...
        with self.assertRaisesRegexp(TypeError, r"flexible\(\) got multiple values for keyword argument 'optional'"):
            flexible('hello', 'world', optional='hello')
        with self.assertRaisesRegexp(TypeError, r"flexible\(\) got multiple values for keyword argument 'optional'"):
            flexible('hello', optional='world', foo='test')

    def test_flexible_kwargs(self):
        class Foo(object):
            @with_flexible_arguments(
                required=[
                    ('test', 'test_arg')
                ],
                optional=[
                    ('optional', 'foo', 'foobar')
                ]
            )
            def flexible(self, **kwargs):
                return kwargs['test'] + (kwargs['optional'] if 'optional' in kwargs else '')
        foo = Foo()
        self.assertEqual(foo.flexible(test='hello'), 'hello')
        self.assertEqual(foo.flexible(test_arg='hello'), 'hello')
        self.assertEqual(foo.flexible(test_arg='hello', foobar=' world'), 'hello world')
        with self.assertRaisesRegexp(TypeError, r'flexible\(\) takes exactly 1 argument \(2 given\)'):
            foo.flexible(2, test=5)

    def test_flecible_example(self):
        @with_flexible_arguments(
            required=[
                ('name', 'greet', 'name_to_greet'),
                ('from_where', 'source')
             ],
            optional=[
                ('greeting', 'hello_phrase', 'salutation'),
                ('message',)
            ]
        )
        def say_hello(name, from_where, greeting='Hello', message=''):
            return greeting + ', ' + name + ', from ' + from_where + '!' + message
        self.assertEqual(say_hello('moon', 'Earth'), "Hello, moon, from Earth!")
        self.assertEqual(say_hello('world', source='the moon'), "Hello, world, from the moon!")
        self.assertEqual(
            say_hello(source='France', name_to_greet='visitor', salutation='Bonjour'),
            "Bonjour, visitor, from France!"
        )
        self.assertEqual(
            say_hello('earthlings', 'outer space', 'Greetings', message='  We come in peace!'),
            "Greetings, earthlings, from outer space!  We come in peace!"
        )


class FooCached(object):

    def __init__(self, d):
        self.d = d

    @CachedMethod
    def bar(self, a, b, c):
        return ((a - b) // c) * self.d

    @CachedMethod
    def baz(self, a, b=2, c=4):
        return (a + 2*b + 4*c) * self.d

class FooImmutable(object):
    __metaclass__ = Immutable

    def __init__(self, d):
        self.d = d

    @CachedMethod
    def bar(self, a, b):
        return (2 * a + 4 * b) * self.d


class CachedMethodTest(unittest.TestCase):

    def test_order(self):
        f = FooCached(3)
        self.assertEqual(f.bar(14, 2, 4), 9)
        self.assertEqual(f.bar(2, 14, 4), -9)
        self.assertEqual(f.bar(4, 14, 2), -15)
        self.assertEqual(f.bar(4, 2, 14), 0)
        self.assertEqual(f.bar(2, 4, 14), -3)
        self.assertEqual(len(f.bar.cache), 5)

    def test_kwargs(self):
        f = FooCached(4)
        self.assertEqual(f.bar(2, 12, 5), -8)
        self.assertEqual(f.bar(2, 12, c=5), -8)
        self.assertEqual(f.bar(b=12, a=2, c=5), -8)
        self.assertEqual(f.bar(a=2, c=5, b=12), -8)
        self.assertEqual(len(f.bar.cache), 1)

    def test_defaults(self):
        f = FooCached(3)
        self.assertEqual(f.baz(1), 63)
        self.assertEqual(f.baz(3, c=1), 33)
        self.assertEqual(f.baz(3, 2, 1), 33)
        self.assertEqual(len(f.baz.cache), 2)

    def test_immutable(self):
        f = FooImmutable(2)
        self.assertEqual(f.bar(3, 6), 60)
        self.assertEqual(f.bar(a=3, b=6), 60)
        self.assertEqual(len(f.bar.cache), 1)

    def test_errors_1(self):
        f = FooCached(2)
        with self.assertRaisesRegexp(TypeError, r'bar\(\) takes exactly 4 arguments \(3 given\)'):
            f.bar(1,2)

    def test_errors_2(self):
        f = FooCached(2)
        with self.assertRaisesRegexp(TypeError, r'bar\(\) takes exactly 4 arguments \(5 given\)'):
            f.bar(1,2,3,4)

    def test_errors_3(self):
        f = FooCached(2)
        with self.assertRaisesRegexp(TypeError, r'baz\(\) takes at least 2 arguments \(1 given\)'):
            f.baz()

    def test_errors_4(self):
        f = FooCached(2)
        with self.assertRaisesRegexp(TypeError, r'baz\(\) takes at most 4 arguments \(5 given\)'):
            f.baz(1,2,3,4)

    def test_errors_5(self):
        f = FooCached(2)
        with self.assertRaisesRegexp(TypeError, r"unexpected keyword argument 'f'"):
            f.bar(1, 2, f=4)

    def test_errors_6(self):
        f = FooCached(2)
        with self.assertRaisesRegexp(TypeError, r"multiple values .* 'b'"):
            f.bar(1, 2, b=4)

    def test_errors_7(self):
        f = FooCached(2)
        with self.assertRaisesRegexp(TypeError, r"multiple values .* 'c'"):
            f.baz(1, 2, 3, c=4)
