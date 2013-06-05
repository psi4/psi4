import unittest

import sys
import os

# Add the directory containing 'grendel' to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir))

# Add the directory containing the 'grendel_tests' package to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from grendel.util.decorators import with_flexible_arguments

class DecoratorTest(unittest.TestCase):

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

    def test_example(self):
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



