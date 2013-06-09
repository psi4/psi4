
from fractions import Fraction
from grendel.differentiation.finite_difference import FiniteDifferenceFormula, FiniteDifferenceDerivative


def load_formulas():
    #--------------------------------------------------------------------------------#
    #                          begin script generated code                           #
    #--------------------------------------------------------------------------------#
    # Finite difference formula for i-type derivatives
    FiniteDifferenceDerivative.formulas.append(
        FiniteDifferenceFormula([1], 2, [
            (Fraction(-1, 2), (-1,)),
            (Fraction(1, 2), (1,)),
        ])
    )
    #--------------------------------------------------------------------------------#
    # Finite difference formula for i-type derivatives with robustness = 3
    FiniteDifferenceDerivative.formulas.append(
        FiniteDifferenceFormula([1], 3, [
            (Fraction(-2, 3), (-1,)),
            (Fraction(2, 3), (1,)),
            (Fraction(1, 12), (-2,)),
            (Fraction(-1, 12), (2,)),
        ])
    )
    #--------------------------------------------------------------------------------#
    # Finite difference formula for ii-type derivatives
    FiniteDifferenceDerivative.formulas.append(
        FiniteDifferenceFormula([2], 2, [
            (Fraction(-2, 1), (0,)),
            (Fraction(1, 1), (-1,)),
            (Fraction(1, 1), (1,)),
        ])
    )
    #--------------------------------------------------------------------------------#
    # Finite difference formula for ii-type derivatives with robustness = 3
    FiniteDifferenceDerivative.formulas.append(
        FiniteDifferenceFormula([2], 3, [
            (Fraction(-5, 2), (0,)),
            (Fraction(4, 3), (-1,)),
            (Fraction(4, 3), (1,)),
            (Fraction(-1, 12), (-2,)),
            (Fraction(-1, 12), (2,)),
        ])
    )
    #--------------------------------------------------------------------------------#
    # Finite difference formula for iii-type derivatives
    FiniteDifferenceDerivative.formulas.append(
        FiniteDifferenceFormula([3], 2, [
            (Fraction(1, 1), (-1, 0, 0)),
            (Fraction(-1, 1), (1, 0, 0)),
            (Fraction(-1, 2), (-2, 0, 0)),
            (Fraction(1, 2), (2, 0, 0)),
        ])
    )
    #--------------------------------------------------------------------------------#
    # Finite difference formula for iiii-type derivatives
    FiniteDifferenceDerivative.formulas.append(
        FiniteDifferenceFormula([4], 2, [
            (Fraction(6, 1), (0, 0)),
            (Fraction(-4, 1), (-1, 0)),
            (Fraction(-4, 1), (1, 0)),
            (Fraction(1, 1), (-2, 0)),
            (Fraction(1, 1), (2, 0)),
        ])
    )
    #--------------------------------------------------------------------------------#
    # Finite difference formula for iiij-type derivatives
    FiniteDifferenceDerivative.formulas.append(
        FiniteDifferenceFormula([3, 1], 2, [
            (Fraction(-3, 1), (0, 0)),
            (Fraction(2, 1), (-1, 0)),
            (Fraction(3, 2), (0, -1)),
            (Fraction(3, 2), (0, 1)),
            (Fraction(2, 1), (1, 0)),
            (Fraction(-1, 2), (-2, 0)),
            (Fraction(-1, 2), (2, 0)),
            (Fraction(-3, 2), (-1, -1)),
            (Fraction(-1, 2), (-1, 1)),
            (Fraction(-1, 2), (1, -1)),
            (Fraction(-3, 2), (1, 1)),
            (Fraction(1, 2), (2, 1)),
            (Fraction(1, 2), (-2, -1)),
        ])
    )
    #--------------------------------------------------------------------------------#
    # Finite difference formula for iij-type derivatives
    FiniteDifferenceDerivative.formulas.append(
        FiniteDifferenceFormula([2, 1], 2, [
            (Fraction(1, 1), (0, -1, 0)),
            (Fraction(-1, 1), (0, 1, 0)),
            (Fraction(-1, 2), (-1, -1, 0)),
            (Fraction(1, 2), (-1, 1, 0)),
            (Fraction(-1, 2), (1, -1, 0)),
            (Fraction(1, 2), (1, 1, 0)),
        ])
    )
    #--------------------------------------------------------------------------------#
    # Finite difference formula for iijj-type derivatives
    FiniteDifferenceDerivative.formulas.append(
        FiniteDifferenceFormula([2, 2], 2, [
            (Fraction(4, 1), (0, 0)),
            (Fraction(-2, 1), (-1, 0)),
            (Fraction(-2, 1), (0, -1)),
            (Fraction(-2, 1), (0, 1)),
            (Fraction(-2, 1), (1, 0)),
            (Fraction(1, 1), (-1, -1)),
            (Fraction(1, 1), (-1, 1)),
            (Fraction(1, 1), (1, -1)),
            (Fraction(1, 1), (1, 1)),
        ])
    )
    #--------------------------------------------------------------------------------#
    # Finite difference formula for iijk-type derivatives
    FiniteDifferenceDerivative.formulas.append(
        FiniteDifferenceFormula([2, 1, 1], 2, [
            (Fraction(-1, 2), (0, -1, -1)),
            (Fraction(1, 2), (0, -1, 1)),
            (Fraction(1, 2), (0, 1, -1)),
            (Fraction(-1, 2), (0, 1, 1)),
            (Fraction(1, 4), (-1, -1, -1)),
            (Fraction(1, 4), (1, 1, 1)),
            (Fraction(-1, 4), (-1, -1, 1)),
            (Fraction(-1, 4), (-1, 1, -1)),
            (Fraction(1, 4), (1, -1, -1)),
            (Fraction(1, 4), (-1, 1, 1)),
            (Fraction(-1, 4), (1, -1, 1)),
            (Fraction(-1, 4), (1, 1, -1)),
        ])
    )
    #--------------------------------------------------------------------------------#
    # Finite difference formula for ij-type derivatives
    FiniteDifferenceDerivative.formulas.append(
        FiniteDifferenceFormula([1, 1], 2, [
            (Fraction(-1, 2), (-1, 0)),
            (Fraction(-1, 2), (0, -1)),
            (Fraction(-1, 2), (0, 1)),
            (Fraction(-1, 2), (1, 0)),
            (Fraction(1, 2), (-1, -1)),
            (Fraction(1, 2), (1, 1)),
            (Fraction(1, 1), (0, 0)),
        ])
    )
    #--------------------------------------------------------------------------------#
    # Finite difference formula for ij-type derivatives with robustness = 3
    FiniteDifferenceDerivative.formulas.append(
        FiniteDifferenceFormula([1, 1], 3, [
            (Fraction(1, 1), (0, 0)),
            (Fraction(-7, 12), (-1, 0)),
            (Fraction(-7, 12), (0, -1)),
            (Fraction(-7, 12), (0, 1)),
            (Fraction(-7, 12), (1, 0)),
            (Fraction(1, 12), (-2, 0)),
            (Fraction(1, 12), (0, -2)),
            (Fraction(1, 12), (0, 2)),
            (Fraction(1, 12), (2, 0)),
            (Fraction(3, 4), (-1, -1)),
            (Fraction(-1, 12), (-1, 1)),
            (Fraction(-1, 12), (1, -1)),
            (Fraction(3, 4), (1, 1)),
            (Fraction(-1, 12), (1, 2)),
            (Fraction(-1, 12), (2, 1)),
            (Fraction(-1, 12), (-2, -1)),
            (Fraction(-1, 12), (-1, -2)),
        ])
    )
    #--------------------------------------------------------------------------------#
    # Finite difference formula for ijk-type derivatives
    FiniteDifferenceDerivative.formulas.append(
        FiniteDifferenceFormula([1, 1, 1], 2, [
            (Fraction(-1, 2), (-1, 0, 0)),
            (Fraction(-1, 2), (0, -1, 0)),
            (Fraction(-1, 2), (0, 0, -1)),
            (Fraction(1, 2), (0, 0, 1)),
            (Fraction(1, 2), (0, 1, 0)),
            (Fraction(1, 2), (1, 0, 0)),
            (Fraction(1, 2), (-1, -1, 0)),
            (Fraction(1, 2), (-1, 0, -1)),
            (Fraction(1, 2), (0, -1, -1)),
            (Fraction(-1, 2), (0, 1, 1)),
            (Fraction(-1, 2), (1, 0, 1)),
            (Fraction(-1, 2), (1, 1, 0)),
            (Fraction(-1, 2), (-1, -1, -1)),
            (Fraction(1, 2), (1, 1, 1)),
        ])
    )
    #--------------------------------------------------------------------------------#
    # Finite difference formula for ijkl-type derivatives
    FiniteDifferenceDerivative.formulas.append(
        FiniteDifferenceFormula([1, 1, 1, 1], 2, [
            (Fraction(-1, 3), (0, 0, 0, 0)),
            (Fraction(1, 3), (-1, -1, 0, 0)),
            (Fraction(1, 3), (-1, 0, -1, 0)),
            (Fraction(1, 3), (-1, 0, 0, -1)),
            (Fraction(1, 3), (0, -1, -1, 0)),
            (Fraction(1, 3), (0, -1, 0, -1)),
            (Fraction(1, 3), (0, 0, -1, -1)),
            (Fraction(-1, 6), (-1, 0, 0, 1)),
            (Fraction(-1, 6), (-1, 0, 1, 0)),
            (Fraction(-1, 6), (-1, 1, 0, 0)),
            (Fraction(-1, 6), (0, -1, 0, 1)),
            (Fraction(-1, 6), (0, -1, 1, 0)),
            (Fraction(-1, 6), (0, 0, -1, 1)),
            (Fraction(-1, 6), (0, 0, 1, -1)),
            (Fraction(-1, 6), (0, 1, -1, 0)),
            (Fraction(-1, 6), (0, 1, 0, -1)),
            (Fraction(-1, 6), (1, -1, 0, 0)),
            (Fraction(-1, 6), (1, 0, -1, 0)),
            (Fraction(-1, 6), (1, 0, 0, -1)),
            (Fraction(1, 3), (0, 0, 1, 1)),
            (Fraction(1, 3), (0, 1, 0, 1)),
            (Fraction(1, 3), (0, 1, 1, 0)),
            (Fraction(1, 3), (1, 0, 0, 1)),
            (Fraction(1, 3), (1, 0, 1, 0)),
            (Fraction(1, 3), (1, 1, 0, 0)),
            (Fraction(-11, 24), (-1, -1, -1, 0)),
            (Fraction(-11, 24), (-1, -1, 0, -1)),
            (Fraction(-11, 24), (-1, 0, -1, -1)),
            (Fraction(-11, 24), (0, -1, -1, -1)),
            (Fraction(-11, 24), (0, 1, 1, 1)),
            (Fraction(-11, 24), (1, 0, 1, 1)),
            (Fraction(-11, 24), (1, 1, 0, 1)),
            (Fraction(-11, 24), (1, 1, 1, 0)),
            (Fraction(1, 24), (-1, -1, 0, 1)),
            (Fraction(1, 24), (-1, -1, 1, 0)),
            (Fraction(1, 24), (-1, 0, -1, 1)),
            (Fraction(1, 24), (-1, 0, 1, -1)),
            (Fraction(1, 24), (-1, 1, -1, 0)),
            (Fraction(1, 24), (-1, 1, 0, -1)),
            (Fraction(1, 24), (0, -1, -1, 1)),
            (Fraction(1, 24), (0, -1, 1, -1)),
            (Fraction(1, 24), (0, 1, -1, -1)),
            (Fraction(1, 24), (1, -1, -1, 0)),
            (Fraction(1, 24), (1, -1, 0, -1)),
            (Fraction(1, 24), (1, 0, -1, -1)),
            (Fraction(1, 24), (-1, 0, 1, 1)),
            (Fraction(1, 24), (-1, 1, 0, 1)),
            (Fraction(1, 24), (-1, 1, 1, 0)),
            (Fraction(1, 24), (0, -1, 1, 1)),
            (Fraction(1, 24), (0, 1, -1, 1)),
            (Fraction(1, 24), (0, 1, 1, -1)),
            (Fraction(1, 24), (1, -1, 0, 1)),
            (Fraction(1, 24), (1, -1, 1, 0)),
            (Fraction(1, 24), (1, 0, -1, 1)),
            (Fraction(1, 24), (1, 0, 1, -1)),
            (Fraction(1, 24), (1, 1, -1, 0)),
            (Fraction(1, 24), (1, 1, 0, -1)),
            (Fraction(1, 2), (1, 1, 1, 1)),
            (Fraction(1, 2), (-1, -1, -1, -1)),
        ])
    )
    #--------------------------------------------------------------------------------#

