from collections import defaultdict
from itertools import groupby, chain, combinations_with_replacement
from copy import copy
#from decimal import Decimal
#import decimal
from math import factorial
from fractions import Fraction
from operator import add
import re
import string
import textwrap
import unittest
import sys
import os

import math

# Add the directory containing 'grendel' to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, os.pardir))

# Add the directory containing the 'grendel_tests' package to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from grendel_tests import long_test
from grendel import *
from grendel.util.strings import indented
from grendel.differentiation.finite_difference import Differentiable, FiniteDifferenceVariable, FiniteDifferenceFunction, FiniteDifferenceDerivative

#--------------------------------------------------------------------------------#
#         Polynomial class to test finite difference derivatives with            #
#--------------------------------------------------------------------------------#

Differentiable.register(Fraction)
FiniteDifferenceVariable.register(str)

class Polynomial(FiniteDifferenceFunction):

    #################
    # Inner classes #
    #################

    class Term(object):

        ##############
        # Attributes #
        ##############

        powers = None
        coeff = None

        ##################
        # Initialization #
        ##################

        def __init__(self, str_or_coeff, pows=None):
            if isinstance(str_or_coeff, str):
                m = re.match(r'\s*([-+]?)\s*(-?\d+\.?\d*(?:/\d+)?)\s*(.*)\s*$', str_or_coeff)
                self.powers = {}
                if m:
                    self.coeff = Fraction((m.group(1) or '') + m.group(2))
                    for part in re.finditer(r'([a-zA-Z])(?:\^?(\d+))?', m.group(3)):
                        if part.group(2) is not None:
                            self.powers[part.group(1)] = int(part.group(2))
                        else:
                            self.powers[part.group(1)] = 1
                else:
                    raise ValueError('{} does not match the syntax of a polynomial term'.format(
                        str_or_coeff
                    ))
            else:
                self.coeff = str_or_coeff
                self.powers = pows

        ##############
        # Properties #
        ##############

        @property
        def vars(self):
            return sorted(self.powers.keys())

        ###################
        # Special Methods #
        ###################

        def __str__(self):
            if self.coeff < 0:
                return '(' + str(self.coeff) + ')' + ''.join((var if self.powers[var] > 0 else '') + ('' if self.powers[var] <= 1 else '^'+str(self.powers[var])) for var in self.vars)
            else:
                return str(self.coeff) + ''.join((var if self.powers[var] > 0 else '') + ('' if self.powers[var] <= 1 else '^'+str(self.powers[var])) for var in self.vars)

        def __call__(self, *args):
            tot = self.coeff
            args = [Fraction(num) for num in args]
            s = sorted(self.powers.keys())
            for var, exp in self.powers.items():
                tot *= args[s.index(var)] ** exp
            return tot

        ###########
        # Methods #
        ###########

        def derivative(self, *vars):
            ret_val = Polynomial.Term(self.coeff, copy(self.powers))
            for var in vars:
                if var not in ret_val.powers or ret_val.powers[var] == 0:
                    new_pows = copy(ret_val.powers)
                    for key in new_pows:
                        new_pows[key] = 0
                    return Polynomial.Term(0.0, new_pows)
                else:
                    new_pows = copy(ret_val.powers)
                    mult = new_pows[var]
                    new_pows[var] -= 1
                    ret_val = Polynomial.Term(mult * ret_val.coeff, new_pows)
            return ret_val

    ####################
    # Class Attributes #
    ####################

    delta = Fraction('0.0001')


    ##############
    # Attributes #
    ##############

    terms = None
    initial_vals = None
    variables = None
    use_fractions = None

    ##################
    # Initialization #
    ##################

    def __init__(self, in_str_or_terms, initial_vals, use_fractions=True):
        self.terms = []
        self.use_fractions = use_fractions
        if isinstance(in_str_or_terms, str):
            for m in re.finditer(r'[-+]?\s*[^\-\+]+', in_str_or_terms):
                self.terms.append(Polynomial.Term(m.group(0)))
        else:
            self.terms = in_str_or_terms[:]
        self.variables = sorted(list(set(chain(*[t.vars for t in self.terms]))))
        if isinstance(initial_vals, dict):
            self.initial_vals = initial_vals
            self.variables = sorted(list(set(chain(*[self.variables, self.initial_vals.keys()]))))
        else:
            self.initial_vals = dict(zip(self.variables, initial_vals))

    ###################
    # Special Methods #
    ###################

    def __str__(self):
        ret_val = ' +#'.join(str(t) for t in self.terms if t.coeff != 0.0)
        tw = textwrap.TextWrapper(subsequent_indent='    ', width=76)
        ret_val = tw.fill(ret_val)
        return ret_val.replace('#', ' ')


    def __call__(self, *args):
        total = Fraction(0) if self.use_fractions else 0.0
        for term in self.terms:
            if not term.coeff == 0:
                idxs = [self.variables.index(v) for v in self.variables if v in term.vars]
                total += term(*[args[idx] for idx in idxs])
        return total

    def __repr__(self):
        return 'Polynomial(' + str(self) + ', ' + str(self.initial_vals) + ')'

    ###########
    # Methods #
    ###########

    def value_for_displacements(self, pairs):
        total = Fraction(0)
        ndeltas = defaultdict(lambda:0)
        ndeltas.update(dict(pairs))
        for term in self.terms:
            total += term(*[self.initial_vals[var] + (ndeltas[var]*self.delta) for var in term.vars])
        return total

    #def deltas_for_variables(self, vars):
    #    return (self.delta,) * len(set(vars))

    # For testing purposes
    def derivative_polynomial(self, *dvars):
        ret_val = Polynomial([], self.initial_vals, use_fractions=self.use_fractions)
        for term in self.terms:
            ret_val.terms.append(term.derivative(*dvars))
        return ret_val

    # For testing purposes
    def actual_derivative(self, *dvars):
        return self.derivative_polynomial(*dvars)(*[self.initial_vals[k] for k in self.variables])

#--------------------------------------------------------------------------------#

class FiniteDifferenceTest(unittest.TestCase):

    def assertOkayFdiffDerivative(self, poly, vars, places=4, robustness=2, forward=False):
        """ Assert almost equal with significant figures rather than digits after the decimal
        """
        df = FiniteDifferenceDerivative(poly, *tuple(vars), robustness=robustness, forward=forward).value
        exp = poly.actual_derivative(*tuple(vars))
        if abs(exp) > 0:
            div = 10.0 ** float(math.ceil(math.log(abs(exp), 10.0)))
            ract = round(df/div, places)
            rexp = round(exp/div, places)
        else:
            ract = round(df, places)
            rexp = round(exp, places)
        tmp = self.longMessage
        self.longMessage = True
        try:
            #("%"+str(places+1)+"f") % rexp, ("%"+str(places+1)+"f") % ract,
            self.assertAlmostEqual(ract, rexp, places,
                msg="\nfor derivative {} with robustness {} of function:\n{}".format(
                    vars, robustness, indented(str(poly)))
            )
        finally:
            self.longMessage = tmp


    def test_one_var(self):
        p = Polynomial('3x^2 + 2x + 1', [5.5])
        f = FiniteDifferenceDerivative(p, 'x')
        self.assertAlmostEqual(f.value, p.actual_derivative('x'))

    def test_three_var(self):
        p = Polynomial('5x^5y^4z^4 + 1x^4 + 1y^4 + 1z^4', [3.2, 2.1, -1.6])
        for i in range(1,5):
            for v in combinations_with_replacement('xyz', i):
                self.assertOkayFdiffDerivative(p, v, places=2)

    @long_test
    def test_seven_var(self):
        p = Polynomial('3.0a^5b^5c^5d^5e^5f^5g^5 + 2a^6 + 3b^6 + 4c^6 + 5d^2 + 6e^4 +7f^4',
            [3.2, 2.1, -1.6, -9.5, -1.2, 3.1, 3.0])
        for i in range(1,5):
            for v in combinations_with_replacement('abcdefg', i):
                self.assertOkayFdiffDerivative(p, v, 2)

    @long_test
    def test_one_var_huge_forward(self):
        terms = []
        for exp in range(1, 45):
            terms.append('{}x^{}'.format(45-exp,exp))
        p = Polynomial(' + '.join(terms), [1.25])
        self.assertOkayFdiffDerivative(p, 'x', robustness=2, places=2, forward=True)
        self.assertOkayFdiffDerivative(p, 'x', robustness=4, places=3, forward=True)
        for rob in range(6, 18):
            self.assertOkayFdiffDerivative(p, 'x', robustness=rob, places=6, forward=True)
        for rob in range(19, 58):
            self.assertOkayFdiffDerivative(p, 'x', robustness=rob, places=10, forward=True)
        self.assertOkayFdiffDerivative(p, 'xx', robustness=1, places=1, forward=True)
        self.assertOkayFdiffDerivative(p, 'xx', robustness=2, places=2, forward=True)
        self.assertOkayFdiffDerivative(p, 'xx', robustness=3, places=3, forward=True)
        self.assertOkayFdiffDerivative(p, 'xx', robustness=4, places=4, forward=True)
        for rob in range(6, 12):
            self.assertOkayFdiffDerivative(p, 'xx', robustness=rob, places=6, forward=True)
        for rob in range(12, 58):
            self.assertOkayFdiffDerivative(p, 'xx', robustness=rob, places=10, forward=True)
        for rob in range(6, 12):
            self.assertOkayFdiffDerivative(p, 'xxx', robustness=rob, places=6, forward=True)
        for rob in range(12, 58):
            self.assertOkayFdiffDerivative(p, 'xxx', robustness=rob, places=10, forward=True)

    @long_test
    def test_one_var_huge(self):
        FiniteDifferenceDerivative.precompute_single_variable(3, 60)
        terms = []
        for exp in range(1, 45):
            terms.append('{}x^{}'.format(45-exp,exp))
        p = Polynomial(' + '.join(terms), [1.25])
        #self.assertOkayFdiffDerivative(p, 'x', robustness=2, places=2)
        #self.assertOkayFdiffDerivative(p, 'x', robustness=4, places=3)
        #for rob in range(6, 12, 2):
        #    self.assertOkayFdiffDerivative(p, 'x', robustness=rob, places=6)
        #for rob in range(12, 58, 2):
        #    self.assertOkayFdiffDerivative(p, 'x', robustness=rob, places=10)
        self.assertOkayFdiffDerivative(p, 'xx', robustness=2, places=2)
        self.assertOkayFdiffDerivative(p, 'xx', robustness=4, places=3)
        for rob in range(6, 12, 1):
            self.assertOkayFdiffDerivative(p, 'xx', robustness=rob, places=6)
        for rob in range(12, 58, 1):
            self.assertOkayFdiffDerivative(p, 'xx', robustness=rob, places=10)
        for rob in range(6, 12, 2):
            self.assertOkayFdiffDerivative(p, 'xxx', robustness=rob, places=6)
        for rob in range(12, 58, 2):
            self.assertOkayFdiffDerivative(p, 'xxx', robustness=rob, places=10)

    @long_test
    def test_one_var_huger(self):
        terms = []
        for exp in range(1, 90):
            terms.append('{}x^{}'.format(120-exp,exp))
        p = Polynomial(' + '.join(terms), [2.25])
        self.assertOkayFdiffDerivative(p, 'x', robustness=2, places=2)
        self.assertOkayFdiffDerivative(p, 'x', robustness=4, places=3)
        for rob in range(6, 12, 2):
            self.assertOkayFdiffDerivative(p, 'x', robustness=rob, places=6)
        for rob in range(12, 48, 2):
            self.assertOkayFdiffDerivative(p, 'x', robustness=rob, places=10)
        self.assertOkayFdiffDerivative(p, 'xx', robustness=2, places=2)
        self.assertOkayFdiffDerivative(p, 'xx', robustness=4, places=3)
        for rob in range(6, 12, 1):
            self.assertOkayFdiffDerivative(p, 'xx', robustness=rob, places=6)
        for rob in range(12, 48, 1):
            self.assertOkayFdiffDerivative(p, 'xx', robustness=rob, places=10)



    #def test_twenty_six_var(self):
    #    terms = []
    #    for coeff in range(1,10):
    #        for letter, exponent in zip(string.ascii_lowercase, range(coeff+10, 37)):
    #            terms.append("{}{}^{}".format(coeff, letter, exponent))
    #    p = Polynomial(' + '.join(terms), [0.56 * i for i in range(1,27)])
    #    for i in range(1,15):
    #        self.assertOkayFdiffDerivative(p, , 2)

