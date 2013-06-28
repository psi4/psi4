import unittest
from copy import deepcopy
import sys
import os

# Add the directory containing 'grendel' to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir))

# Add the directory containing the 'grendel_tests' package to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from grendel import *
from grendel_tests import expected_failure

class ParsingTest(unittest.TestCase):

    # A few test strings that get reused...
    test_strings = [
        """
        something before the start sequence
        ---------------
        irrelevent part
        ---------------
        ...
        ...
        something else before the start sequence
        -------------------
        vvvvvvvvvvvvvvvvvvv
        relevant
        data:
        3.14159
        is in this
        part of the
        file:
        2.71828
        and some other stuff
        ^^^^^^^^^^^^^^^^^^^
        -------------------
        something after the end sequence
        numbers that shouldn't be matched:
        3.1
        2.1
        1.2
        1.3
        4.5.6
        """,
        #--------------------------------------------------------------------------------#
        # for the test_multi_repeated series
        """

        1 2 3
        a b c
        4 5 6
        d e f
        7 8 9
        g h i

        """,
        #--------------------------------------------------------------------------------#
    ]

    def test_basic(self):
        test_str =  """
         -------
         oooooo
         hello
         oooooo
         world
         oooooo
         -------
        """
        rep = RegexSequence(
            start = r'^\s*--+',
            main = r'\s+([^aeiou\s]+).+?([^aeiou\s]+)',  # match consonent sets
            stop = r'^\s*--+'
        )
        for line in test_str.splitlines():
            rep.parse_line(line)
        self.assertEqual(len(rep.matches), 2)
        self.assertEqual(len(rep.groups()), 4)
        self.assertSequenceEqual(rep.groups(), ('h', 'll', 'w', 'rld'))
        #--------------------------------------------------------------------------------#
        # last match only
        rep = RegexSequence(
            start = r'^\s*--+',
            main = r'\s+([^aeiou\s]+).+?([^aeiou\s]+)',  # match consonent sets
            stop = r'^\s*--+',
            last_only = True
        )
        for line in test_str.splitlines():
            rep.parse_line(line)
        self.assertEqual(len(rep.matches), 1)
        self.assertEqual(len(rep.groups()), 2)
        self.assertSequenceEqual(rep.groups(), ('w', 'rld'))

    def test_sequences_1(self):
        test_str = self.test_strings[0]
        rseq = RegexSequence(
            start = [
                r'^\s*-+\s*$',
                r'^\s*v+\s*$',
            ],
            main = [
                r'^\s*(\w+):\s*$',
                r'^\s*(\d+)\.(\d+)$',
            ],
            stop = [
                r'^\s*\^+\s*$',
                r'^\s*-+\s*$',
            ],
        )
        for line in test_str.splitlines():
            rseq.parse_line(line)
        self.assertEqual(len(rseq.matches), 2)
        self.assertSequenceEqual(rseq.groups(), ('data', '3', '14159', 'file', '2', '71828'))
        self.assertSequenceEqual(rseq.groups(False), (('data', '3', '14159'), ('file', '2', '71828')))

    def test_ellipsis(self):
        # Test the ellipsis/plain repeated regex
        test_str = self.test_strings[0]
        rseq = RegexSequence(
            start = [
                r'^\s*-+\s*$',
                r'^\s*v+\s*$',
            ],
            main = [
                r'^\s*(\d+)\.(\d+)$',
                Ellipsis,
                r'^\s*(\d+)\.(\d+)$',
            ],
            stop = [
                r'^\s*\^+\s*$',
                r'^\s*-+\s*$',
            ],
        )
        for line in test_str.splitlines():
            rseq.parse_line(line)
        self.assertEqual(len(rseq.matches), 1)
        self.assertSequenceEqual(rseq.groups(), ('3', '14159', '2', '71828'))
        self.assertSequenceEqual(rseq.groups(False), (('3', '14159', '2', '71828'),))

    def test_should_not_match(self):
        # slightly different version of the test_ellipsis that should not match...
        test_str = self.test_strings[0]
        rseq = RegexSequence(
            start = [
                r'^\s*-+\s*$',
                r'^\s*v+\s*$',
            ],
            main = [
                r'^\s*(\d+)\.(\d+)$',
                RegexSequence.RepeatedRegex(r'is'),
                r'^\s*(\d+)\.(\d+)$',
                r'.*other stuff'
            ],
            stop = [
                r'^\s*\^+\s*$',
                r'^\s*-+\s*$',
            ],
        )
        for line in test_str.splitlines():
            rseq.parse_line(line)
        self.assertEqual(len(rseq.matches), 0)

    def test_sequences_2(self):
        test_str = """
        -------------------
        vvvvvvvvvvvvvvvvvvv
        numbers that should be matched:
        3.14159
        ^^^^^^^^^^^^^^^^^^^
        -------------------
        numbers that shouldn't be matched:
        2.71828
        """
        #--------------------------------------------------------------------------------#
        # Length 1 sequences and beginning/end repeats:
        rseq = RegexSequence(
            start = [
                RegexSequence.RepeatedRegex(r'^\s*-+\s*$', 2),
                r'^\s*v+\s*$',
            ],
            main = [
                r'^\s*(\d+)\.(\d+)$',
            ],
            stop = [
                RegexSequence.RepeatedRegex(r'^\s*\^+\s*$', 1),
                RegexSequence.RepeatedRegex(r'^\s*-+\s*$', 2)
            ],
        )
        for line in test_str.splitlines():
            rseq.parse_line(line)
        self.assertEqual(len(rseq.matches), 1)
        self.assertSequenceEqual(rseq.groups(), ('3', '14159'))

    def test_multi_repeated_1(self):
        test_str = self.test_strings[1]
        rseq = RegexSequence(
            main = [
                RegexSequence.RepeatedRegex([
                    r'\d \d (\d)',
                    r'\w \w \w'
                ]),
                r'^\s*$'
            ]
        )
        for line in test_str.splitlines():
            rseq.parse_line(line)
        self.assertSequenceEqual(rseq.groups(), ('3','6','9'))

    def test_multi_repeated_2(self):
        test_str = self.test_strings[1]
        rseq = RegexSequence(
            main = [
                RegexSequence.RepeatedRegex([
                    r'\d \d (\d)',
                    r'\w \w \w'
                ]),
            ]
        )
        for line in test_str.splitlines():
            rseq.parse_line(line)
        self.assertSequenceEqual(rseq.groups(), ('3','6','9'))

    def test_multi_repeated_3(self):
        test_str = self.test_strings[1]
        rseq = RegexSequence(
            main = [
                Repeated([
                    Duplicated(r'(\d)\s*', min_dups=0),
                    Duplicated(r'\w\s*')
                ]),
                r'^\s*$',
                ]
        )
        for line in test_str.splitlines():
            rseq.parse_line(line)
        self.assertSequenceEqual(rseq.groups(), ('1', '2', '3', '4', '5', '6', '7', '8', '9'))

    def test_repeated_special(self):
        test_str = """
        -------------------
        vvvvvvvvvvvvvvvvvvv
        numbers that should be matched:
        3.14159
        ^^^^^^^^^^^^^^^^^^^
        -------------------
        numbers that shouldn't be matched:
        2.71828
        """
        #--------------------------------------------------------------------------------#
        # Make sure starts:
        rseq = RegexSequence(
            start = [
                r'^\s*-+\s*$',
                RegexSequence.RepeatedRegex(r'', min_repeats=0),
                RegexSequence.RepeatedRegex(r'.*', min_repeats=0),
                r'^\s*v+\s*$',
            ],
            main = [
                r'^\s*(\d+)\.(\d+)$',
            ],
            stop = [
                r'^\s*\^+\s*$',
                r'^\s*-+\s*$',
            ],
        )
        for line in test_str.splitlines():
            rseq.parse_line(line)
        self.assertEqual(len(rseq.matches), 1)
        self.assertSequenceEqual(rseq.groups(), ('3', '14159'))
        #--------------------------------------------------------------------------------#
        # Make sure this stops:
        rseq = RegexSequence(
            main = [
                r'^\s*(\d+)\.(\d+)$',
                ],
            stop = [
                r'^\s*\^+\s*$',
                RegexSequence.RepeatedRegex(r'', 10, min_repeats=0),
                RegexSequence.RepeatedRegex(r'.*', 3, min_repeats=0),
                r'^\s*-+\s*$',
            ],
        )
        for line in test_str.splitlines():
            rseq.parse_line(line)
        self.assertEqual(len(rseq.matches), 1)
        self.assertSequenceEqual(rseq.groups(), ('3', '14159'))

    def test_deepcopy(self):
        test_str = self.test_strings[0]
        # initially just the same test as the test_sequences test, above
        rseq = RegexSequence(
            start = [
                r'^\s*-+\s*$',
                r'^\s*v+\s*$',
                ],
            main = [
                r'^\s*(\w+):\s*$',
                r'^\s*(\d+)\.(\d+)$',
                ],
            stop = [
                r'^\s*\^+\s*$',
                r'^\s*-+\s*$',
                ],
        )
        # Do some parsing on the first instance
        for line in test_str.splitlines():
            rseq.parse_line(line)
        # Now do a deepcopy
        rseq2 = deepcopy(rseq)
        # The original should remain unchanged
        self.assertEqual(len(rseq.matches), 2)
        self.assertSequenceEqual(rseq.groups(), ('data', '3', '14159', 'file', '2', '71828'))
        self.assertSequenceEqual(rseq.groups(False), (('data', '3', '14159'), ('file', '2', '71828')))
        # The copy should act like nothing has been done
        self.assertEqual(len(rseq2.matches), 0)
        # and the copy should work as expected
        for line in test_str.splitlines():
            rseq2.parse_line(line)
        # and pass the same tests as the original:
        self.assertEqual(len(rseq2.matches), 2)
        self.assertSequenceEqual(rseq2.groups(), ('data', '3', '14159', 'file', '2', '71828'))
        self.assertSequenceEqual(rseq2.groups(False), (('data', '3', '14159'), ('file', '2', '71828')))

    def test_deepcopy_multi(self):
        # same test as test_multi_repeated_3
        test_str = self.test_strings[1]
        rseq = RegexSequence(
            main = [
                Repeated([
                    Duplicated(r'(\d)\s*', min_dups=0),
                    Duplicated(r'\w\s*')
                ]),
                r'^\s*$',
                ]
        )
        for line in test_str.splitlines():
            rseq.parse_line(line)
        rseq2 = deepcopy(rseq)
        # The original should remain unchanged
        self.assertEqual(len(rseq.matches), 1)
        self.assertSequenceEqual(rseq.groups(), ('1', '2', '3', '4', '5', '6', '7', '8', '9'))
        # The copy should act like nothing has been done
        self.assertEqual(len(rseq2.matches), 0)
        # and the copy should work as expected
        for line in test_str.splitlines():
            rseq2.parse_line(line)
            # and pass the same tests as the original:
        self.assertEqual(len(rseq2.matches), 1)
        self.assertSequenceEqual(rseq2.groups(), ('1', '2', '3', '4', '5', '6', '7', '8', '9'))

    def test_greedy_1(self):
        test_str = """
        aaa1
        aaa2
        aaa3
        bbbb
        """
        rseq = RegexSequence(
            main = [
                Repeated(r'(a)aa', greedy=True)
            ]
        )
        for line in test_str.splitlines():
            rseq.parse_line(line)
        self.assertEqual(len(rseq.matches), 1)
        self.assertEqual(len(rseq.matches[0]), 3)
        self.assertSequenceEqual(rseq.groups(False), (('a','a','a'),))

    def test_greedy_2(self):
        test_str = """
        aaa1
        aaa2
        aaa3
        bbbb
        aaa4
        aaa5
        bbbbbb
        aaa6"""
        rseq = RegexSequence(
            main = [
                Repeated(r'aaa(\d)')
            ]
        )
        for line in test_str.splitlines():
            rseq.parse_line(line)
        self.assertEqual(len(rseq.matches), 3)
        self.assertEqual(len(rseq.matches[0]), 3)
        self.assertEqual(len(rseq.matches[1]), 2)
        self.assertSequenceEqual(rseq.groups(False), (('1','2','3'), ('4','5'), ('6',)))
        # Now try the same thing with greed off
        rseq = RegexSequence(
            main = [
                Repeated(r'aaa(\d)', greedy=False)
            ]
        )
        for line in test_str.splitlines():
            rseq.parse_line(line)
        self.assertEqual(len(rseq.matches), 6)
        self.assertEqual(len(rseq.matches[0]), 1)
        self.assertSequenceEqual(rseq.groups(False), (('1',),('2',),('3',),('4',),('5',),('6',)))

    @expected_failure(allow_exceptions=[ValueError])
    def test_exact_1(self):
        # Should raise a value error.
        test_str = """
        aaa1
        aaa2
        aaa3
        bbbb
        aaa4
        aaa5
        bbbbbb
        aaa6
        """
        rseq = RegexSequence(
            main = [
                Repeated(r'aaa(\d)', repeats=2, exactly=True)
            ]
        )
        for line in test_str.splitlines():
            rseq.parse_line(line)
        self.assertEqual(len(rseq.matches), 1)
        self.assertEqual(len(rseq.matches[0]), 2)
        self.assertSequenceEqual(rseq.groups(False), (('4','5'),))

    def test_repeats(self):
        # Not entirely sure what the most intuitive thing to do with this case is...
        # all three possibilities could be matched with a "look-behind" like feature
        # that keeps a temporary array of lines matched since the start and goes back to the
        # second line of that array when the match is finished...that's pretty complicated though
        test_str = """
        aaa1
        aaa2
        aaa3
        bbbb
        aaa4
        aaa5
        bbbbbb
        aaa6
        """
        rseq = RegexSequence(
            main = [
                r'aaa(\d)',
                r'aaa(\d)'
            ]
        )
        for line in test_str.splitlines():
            rseq.parse_line(line)
        self.assertEqual(len(rseq.matches), 2)
        self.assertEqual(len(rseq.matches[0]), 2)
        self.assertSequenceEqual(rseq.groups(False), (('1','2'), ('4','5')))

