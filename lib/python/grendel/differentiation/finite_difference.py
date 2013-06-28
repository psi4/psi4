"""
"""
from abc import ABCMeta, abstractproperty, abstractmethod
from fractions import Fraction
from itertools import groupby
from operator import add, attrgetter, mul

from grendel import type_checking_enabled, sanity_checking_enabled
import string
import math
from grendel.gmath.tensor import Tensor
from grendel.util.decorators import with_flexible_arguments, typechecked, IterableOf
from grendel.util.metaprogramming import ReadOnlyAttribute
from grendel.util.overloading import listify_args
from grendel.util.strings import indented


#class FiniteDifferenceDerivativeCollection(object):
#    """
#    Note:  Not a subclass of `Tensor`
#    """
#
#    ######################
#    # Private Attributes #
#    ######################
#
#    _value_tens = None
#
#    ##################
#    # Initialization #
#    ##################
#
#    @typechecked(
#        function="FiniteDifferenceFunction",
#        variable_list=IterableOf('FiniteDifferenceVariable'),
#        max_order=int
#    )
#    def __init__(self, function, variable_list, max_order):
#        pass


#noinspection PyTypeChecker
class FiniteDifferenceDerivative(object):
    """ An arbitrary finite difference derivative with respect to k `FiniteDifferenceVariable` instances (not necessarily distinct)
    of a `FiniteDifferenceFunction` of an arbitrary number of `FiniteDifferenceVariable` instances whose output is a `Differentiable` instance.
    Ideally, the derivative can be calculated to an arbitrary order of robustness in the displacement, though in
    practice not all orders may be implemented yet.
    """

    ####################
    # Class Attributes #
    ####################

    formulas = []
    generated_single_variable_formulas = {}

    ##############
    # Attributes #
    ##############

    function = ReadOnlyAttribute('function',
        doc="""The `FiniteDifferenceFunction` instance to differentiate."""
    )

    variables = ReadOnlyAttribute('variables',
        doc="""The list of `FiniteDifferenceVariable` instances to be displaced for computation of the finite difference derivative."""
    )

    target_robustness = ReadOnlyAttribute('target_robustness',
        doc="""The minimum order of the error in the displacement.  Defaults to 2."""
    )

    formula = ReadOnlyAttribute('formula',
        doc="""The FiniteDifferenceFormula object we need to use, based on the input parameters."""
    )

    orders = ReadOnlyAttribute('orders')

    ######################
    # Private Attributes #
    ######################

    _value = None
    _value_function = None
    _delta_function = None
    _delta = None
    _forward = None


    ##################
    # Initialization #
    ##################

    @with_flexible_arguments(
        optional = [
            ('target_robustness', 'robustness', 'accuracy', 'order', 'correct_to_order')
        ]
    )
    @typechecked(function='FiniteDifferenceFunction')
    def __init__(self, function, *variables, **kwargs):
        """
        """
        if len(FiniteDifferenceDerivative.formulas) == 0:
            # Load the formulas generated "by hand", which (for now, anyway) require fewer
            #   displacements than the automatically generated formulas if we also need to
            #   compute the lower order derivatives as well, as is the case with the computation
            #   of quartic forcefields. (But not, for instance, the B tensor.  So the
            #   FiniteDifferenceDerivative constructor could be optimized to take a parameter
            #   which specifies whether we should choose the formula with the fewest overall
            #   displacements or the fewest "new" displacements not needed for smaller derivatives)
            load_formulas()
        #--------------------------------------------------------------------------------#
        # miscellanea
        self._target_robustness = kwargs.pop('target_robustness', 2)
        self._value_function = kwargs.pop('value_function', None)
        self._delta_function = kwargs.pop('delta_function', None)
        self._delta = kwargs.pop('delta', None)
        self._forward = kwargs.pop('forward', False)
        self._function = function
        #--------------------------------------------------------------------------------#
        # type checking
        if type_checking_enabled:
            if not all(isinstance(v, FiniteDifferenceVariable) for v in variables):
                raise TypeError
            if not isinstance(self.target_robustness, int):
                raise TypeError
        #--------------------------------------------------------------------------------#
        # Get the variables and the orders....
        vars = listify_args(*variables)
        # Determine which formula we need
        vars = sorted(vars, key=id)
        # This is nasty, but it works...The zip(*list_of_lists) effectively "unzips"
        self._orders, self._variables = zip(
            *sorted(
                [(len(list(g)), k) for k, g in groupby(vars)],
                reverse=True)
        )
        #--------------------------------------------------------------------------------#
        # Determine which formula to use
        # This gets reused, so define a quicky function...
        def get_possibilities(formula_list):
            return [f for f in formula_list
                if f.orders == list(self.orders)
                        and f.robustness >= self.target_robustness
                        and (f.is_forward() if self._forward else f.is_central())
            ]
        #----------------------------------------#
        # First, try and get a "hand-generated" formula
        possibilities = get_possibilities(FiniteDifferenceDerivative.formulas)
        if len(possibilities) == 0:
            # We know how to generate single variable formulas to arbitrary order, so let's do it
            n_derivative_vars = len(self.orders)
            derivative_order = sum(self.orders)
            if n_derivative_vars == 1:
                # This long name is unweildy...
                gen_dict = FiniteDifferenceDerivative.generated_single_variable_formulas
                # See if we've already generated it...
                formula = gen_dict.get(
                    (
                        derivative_order,
                        self.target_robustness
                            + (1 if not self._forward and self.target_robustness % 2 == 1 else 0),
                        self._forward
                    ),
                    None)
                if formula is None:
                    # okay, we can generate it.
                    generate_single_variable_formulas(
                        derivative_order,
                        self.target_robustness
                            + (1 if not self._forward and self.target_robustness % 2 == 1 else 0),
                        self._forward)
                    formula = gen_dict[(
                        derivative_order,
                        self.target_robustness
                            + (1 if not self._forward and self.target_robustness % 2 == 1 else 0),
                        self._forward)]
                possibilities.append(formula)
                if sanity_checking_enabled:
                    possibilities = get_possibilities(possibilities)
            else:
                # we don't know how to generate these...yet...but I'm working on it!
                raise RuntimeError("Can't find formula for orders {0} and"
                                   " robustness {1}".format(
                    self.orders, self.target_robustness))
        # Use the minimum robustness for now.  Later we can make it use
        #   the best possible without additional calculations.
        self._formula = sorted(possibilities, key=attrgetter('robustness'))[0]

    ##############
    # Properties #
    ##############

    @property
    def value(self):
        if self._value is None:
            self.compute()
        return self._value

    @property
    def needed_increments(self):
        return self.formula.coefficients.keys()

    #################
    # Class Methods #
    #################

    @classmethod
    def precompute_single_variable(cls, max_derivative, max_order, forward=False):
        """ Save a little bit of time by prepopulating the single variable displacement
        formulas dictionary up to `max_derivative` and `max_order`.  If `forward` is True,
        the forward formulas are precomputed instead of the central ones.
        """
        generate_single_variable_formulas(max_derivative, max_order, forward)



    ###########
    # Methods #
    ###########

    def compute(self):
        #TODO handle units (efficiently!!!)
        if self._value is not None:
            return self._value
        total = None
        for increments, coeff in self.formula.coefficients.items():
            if self._value_function:
                tmp = self._value_function(zip(self.variables, increments))
            else:
                tmp = self.function.value_for_displacements(zip(self.variables, increments))
            if tmp is None:
                raise ValueError("the value_for_displacements method of FiniteDifferenceFunction"
                                 " '{}' returned `None` for increments {}".format(
                    self.function, increments
                ))
            if hasattr(tmp, 'value'):
                val = tmp.value * coeff
            else:
                val = tmp * coeff
            if total is not None:
                total += val
            else:
                total = val
        if self._delta is not None:
            deltas = (self._delta,) * len(set(self.variables))
        elif self._delta_function:
            deltas = self._delta_function(self.variables)
        else:
            deltas = self.function.deltas_for_variables(self.variables)
        if isinstance(total, Fraction):
            # Try and keep it that way
            denom = reduce(mul, [d**exp for d, exp in zip(deltas, self.orders)])
            total /= denom
        else:
            total /= reduce(mul, Tensor(deltas)**Tensor(self.orders))
        self._value = total




class FiniteDifferenceFunction(object):
    """
    """
    __metaclass__ = ABCMeta

    @property
    def variables(self):
        """The list of `FiniteDifferenceVariable` instances on which the function depends."""
        raise NotImplementedError

    @abstractmethod
    def value_for_displacements(self, pairs):
        """ Get a value for the displacement corresponding to the (variable, number of deltas) pairs given in the argument `pairs`
        """
        return NotImplemented

    def deltas_for_variables(self, vars):
        """ The displacement amounts for each of the variables.
        """
        if hasattr(self, 'delta'):
            return (self.delta,) * len(set(vars))
        else:
            raise NotImplementedError


class FiniteDifferenceVariable(object):
    """
    """
    __metaclass__ = ABCMeta


class Differentiable(object):
    """ Abstract base class for things that you are allowed to take derivatives of.
    """
    __metaclass__ = ABCMeta

    @property
    def value(self):
        """ The value of the differentiable property.  If this function is not
        overridden, assume that `self` can be added and subtracted, as well
        as multiplied by a float.
        """
        return self

    @property
    def shape(self):
        return tuple()

# Simplest possible subclass of Differentiable that acts like a float
class FloatShell(float, Differentiable): pass


class FiniteDifferenceFormula(object):
    """ A formula for a given finite difference derivative.
    """

    #############
    # Constants #
    #############

    CENTRAL = 0
    FORWARD = 1

    ##############
    # Attributes #
    ##############

    orders = ReadOnlyAttribute('orders')
    robustness = ReadOnlyAttribute('robustness')
    coefficients = ReadOnlyAttribute('coefficients',
        """ displacement => coefficient dictionary, where
        the elements of the displacement tuple correspond
        to the orders tuple
        """
    )
    direction = None

    ##################
    # Initialization #
    ##################

    def __init__(self, orders, robustness, pairs, forward=False):
        self._orders = orders
        self._robustness = robustness
        # pairs is a list of (coefficient, displacement tuple) tuples
        self._coefficients = dict((p[1], p[0]) for p in pairs)
        if forward:
            self.direction = FiniteDifferenceFormula.FORWARD
        else:
            self.direction = FiniteDifferenceFormula.CENTRAL

    ###################
    # Special Methods #
    ###################

    #------------------------#
    # Output Representations #
    #------------------------#

    #def __eq__(self, other):
    #    return self._robustness

    def __str__(self):
        max_width=80
        indent_size = 4
        function_name = 'F'
        disp_name = 'h'
        #----------------------------------------#
        def var(idx):
            vars = 'xyz' + str(reversed(string.ascii_lowercase[:-3]))
            vars = vars.replace(disp_name, '')
            vars = vars.replace(function_name, '')
            try:
                return vars[idx]
            except IndexError:
                return 'x_' + idx
        #----------------------------------------#
        flines = []
        curr_line = '['
        for dispnum, (disp, coeff) in enumerate(sorted(self.coefficients.items())):
            if coeff == 0:
                continue
            if coeff.denominator == 1:
                # Only print non-unit coefficients
                if abs(coeff) != 1:
                    disp_f = str(abs(coeff))
                else:
                    disp_f = ''
            else:
                disp_f = '(' + str(abs(coeff)) + ')'
            disp_f += function_name + '('
            for idx, d in enumerate(disp):
                disp_f += var(idx)
                if d == -1:
                    disp_f += ' - ' + disp_name + '_' + var(idx)
                elif d == 1:
                    disp_f += ' + ' + disp_name + '_' + var(idx)
                elif d < 0:
                    disp_f += ' - ' + str(abs(d)) + disp_name + '_' + var(idx)
                elif d > 0:
                    disp_f += ' + ' + str(d) + disp_name + '_' + var(idx)
                if d is not disp[-1]:
                    disp_f += ', '
            disp_f += ')'
            if dispnum != 0:
                if coeff < 0:
                    disp_f = ' - ' + disp_f
                else:
                    disp_f = ' + ' + disp_f
            elif coeff < 0:  # and dispnum == 0
                disp_f = '-' + disp_f
            if len(curr_line + disp_f) > (max_width if len(flines) == 0 else max_width - indent_size):
                if curr_line != '':
                    flines.append(curr_line)
                    curr_line = disp_f
                else:
                    # one term per line is the best we can do, and we still overflow...yikes...
                    flines.append(disp_f)
            else:
                curr_line += disp_f
        denom = '] / '
        if len(self.orders) > 1:
            denom += '('
        for idx, exp in enumerate(self.orders):
            denom += disp_name + '_' + var(idx)
            if exp > 1:
                denom += '^' + str(exp)
        if len(self.orders) > 1:
            denom += ')'
        if len(curr_line + denom) > (max_width if len(flines) == 0 else max_width - indent_size):
            if curr_line != '':
                flines.append(curr_line)
                curr_line = denom
            else:
                curr_line = denom
        else:
            curr_line += denom
        flines.append(curr_line)

        return '{cent} finite difference formula for df/{dvars},' \
             ' correct to order {order} in displacement ({terms} terms):\n{formula}'.format(
            cent='Central' if self.direction == FiniteDifferenceFormula.CENTRAL else 'Forward',
            order = self.robustness,
            dvars=''.join('d' + var(idx) + ('^'+str(exp) if exp > 1 else '')
                              for idx, exp in enumerate(self.orders)),
            formula=indented(('\n' + ' ' * indent_size).join(flines), indent_size),
            terms=len([c for c in self.coefficients.values() if c != 0])
        )


    def __repr__(self):
        return "FiniteDifferenceFormula({dord}, {rob}, [\n{terms}\n])".format(
            dord=repr(self.orders),
            rob=self.robustness,
            terms=indented(
                ',\n'.join(repr((coeff, deltas)) for deltas, coeff in
                    sorted(self.coefficients.items(), key=lambda x: ' '.join(str(i) for i in x[0]))
                        if coeff != 0)
                )
        )

    ##############
    # Properties #
    ##############

    @property
    def order(self):
        return reduce(add, self.orders)

    @property
    def needed_displacements(self):
        return self.coefficients.items()

    ###########
    # Methods #
    ###########

    #-----------------#
    # Inquiry methods #
    #-----------------#

    def is_forward(self):
        return self.direction == FiniteDifferenceFormula.FORWARD

    def is_central(self):
        return self.direction == FiniteDifferenceFormula.CENTRAL



def generate_single_variable_formulas(max_derivative, max_order, forward=False):
    # generate the 'number of deltas' list
    fornberg_N = max_derivative + max_order - 1
    if forward:
        disps = range(fornberg_N + 1)
    else:
        disps = sum(([i, -i] for i in range(1, int(math.ceil(float(fornberg_N)/2.0)) + 1)), [0])
    # Generate the formulas!
    deltas = fornberg_coefficient_generator(max_derivative, fornberg_N, 0, *disps)
    # Long name is too unweildy...
    gen_formulas = FiniteDifferenceDerivative.generated_single_variable_formulas
    # Parse out the coefficients from the three-dimensional array returned by Fornberg's algorithm
    #  into FiniteDifferenceFormula instances
    for derivative in xrange(1, max_derivative+1):
        for n in xrange(derivative, max_derivative + max_order):
            robustness = n - derivative + 1
            if not forward and robustness % 2 == 1:
                continue
            terms = []
            coeffs = deltas[derivative][n]
            for nu, coeff in enumerate(coeffs):
                terms.append((coeff, (disps[nu],)))
            # Just overwrite what was there (even though it is exactly the same,
            #    it takes longer and is more complicated  to check for it's existence
            #    than to just overwrite what's already there, since we have it anyway)
            gen_formulas[(derivative, robustness, forward)] = \
                FiniteDifferenceFormula([derivative], robustness, terms, forward=forward)


def fornberg_coefficient_generator(M, N, x0, *alpha):
    """
    From Fornberg, Bengt. Mathematics of Computation 51 (1988), p. 699-706
    M: "the order of the highest derivative we wish to approximate"
    N: given "a set of N + 1 grid points"
    x0: the point at which we wish to take the derivative
    alpha: alpha_i is the ith grid point
    """
    delta = [[[0 for nu in xrange(n+1)] for n in xrange(N+1)] for m in xrange(M+1)]
    delta[0][0][0] = Fraction(1)
    c1 = Fraction(1)
    for n in xrange(1, N+1):
        c2 = Fraction(1)
        alpha = tuple(Fraction(a) for a in alpha)
        for nu in xrange(n):
            c3 = alpha[n] - alpha[nu]
            c2 = c2 * c3
            for m in range(min(n, M)+1):
                delta[m][n][nu] = ((alpha[n] - x0) * delta[m][n-1][nu]
                                       - (m * delta[m-1][n-1][nu] if m >= 0 else 0))/c3
        for m in range(min(n, M)+1):
            delta[m][n][n] = c1/c2 * (m * (delta[m-1][n-1][n-1] if m >= 0 else 0)
                                          - (alpha[n-1] - x0) * delta[m][n-1][n-1])
        c1 = c2
    return delta


#####################
# Dependent Imports #
#####################

from grendel.differentiation.fdiff_formulas import load_formulas

