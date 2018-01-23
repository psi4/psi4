import re

import numpy as np

from .. import periodictable
from ..exceptions import *



def reconcile_nucleus(A=None, Z=None, E=None, mass=None, label=None, nonphysical=False, mtol=1.e-3, verbose=1):
    """

    Considers all the atom identity information available from arguments,
    supplemented by the periodic table. At the least, must provide element
    information somehow.

    Parameters
    ----------
    A : int, optional
        Mass number, number of protons and neutrons.
    Z : int, optional
        Atomic number, number of protons.
    E : str, optional
        Element symbol from periodic table.
    mass : float, optional
        Atomic mass [u].
    label : str, optional
        TODO
    nonphysical : bool, optional
        When `True`, turns off sanity checking that prevents periodic table
        violations (e.g, light uranium: 1U@1.007).
    mtol : float, optional
        How different `mass` can be from a known nuclide mass and still
        merit the mass number assignment. Note that for elements dominated
        by a single isotope, the default may not be tight enough to
        prevent standard atomic weight (abundance-average of isotopes)
        from being labeled as the dominant isotope for A.
    verbose : int, optional
        Quantity of printing.

    Returns
    -------
    A, Z, E, mass : int or '(none)', int, str, float
        mass number, unless clues don't point to a known nuclide, in which case '(none)'.
        atomic number
        element symbol
        mass [u]

    """
    def reconcile(exact, tests, feature):
        """Returns a member from `exact` that passes all `tests`, else raises error for `feature`."""

        for candidate in exact:
            assessment = [fn(candidate) for fn in tests]
            text.append("""Assess {} candidate {}: {} --> {}""".format(feature, candidate, assessment, all(assessment)))
            if all(assessment):
                return candidate
        else:
            err = """Inconsistent or unspecified {}: A: {}, Z: {}, E: {}, mass: {}, label: {}""".format(
                                  feature, A, Z, E, mass, label)
            if verbose > -1:
                print('\n\n' + '\n'.join(text))
            raise ValidationError(err)

    def offer_atomic_number(z):
        """Given an element, what can be suggested and asserted about Z, A, mass?"""

        z_symbol = periodictable.z2el[z]
        z_mass = periodictable.z2mass[z]
        re_eliso = re.compile(z_symbol + '[0-9]{1,3}')  # lone symbol (val equals most common isotope) will not match
        z_a2mass = {int(k[len(z_symbol):]): v for k, v in periodictable.eliso2mass.items() if re_eliso.match(k)}
        z_a2mass_min = min(z_a2mass.keys())
        z_a2mass_max = max(z_a2mass.keys())
        z_mass2a = {v: k for k, v in z_a2mass.items()}
        z_mass2a_min = min(z_mass2a.keys())
        z_mass2a_max = max(z_mass2a.keys())
        z_a = z_mass2a[z_mass]

        Z_exact.append(z)
        Z_range.append(lambda x: x == z)

        A_exact.append(z_a)
        if nonphysical:
            A_range.append(lambda x: x == '(none)' or x >= 1)
            text.append("""For A, input Z: {}, requires 1 < A or (none), nonphysical""".format(z))
        else:
            A_range.append(lambda x: x == '(none)' or (x >= z_a2mass_min and x <= z_a2mass_max))
            text.append("""For A, input Z: {} requires {} < A < {} or (none), the known mass number range for element""".format(z, z_a2mass_min, z_a2mass_max))

        m_exact.append(z_mass)
        if nonphysical:
            m_range.append(lambda x: x > 0.5)
            text.append("""For mass, input Z: {}, requires 0.5 < mass, nonphysical""".format(z))
        else:
            mmtol = 0.5
            m_range.append(lambda x: x >= z_mass2a_min - mmtol and x <= z_mass2a_max + mmtol)
            text.append("""For mass, input Z: {} requires {} < mass < {} +/-{}, the known mass range for element""".format(z, z_mass2a_min, z_mass2a_max, mmtol))

    def offer_mass_number(z, a):
        """Given a mass number and element, what can be suggested and asserted about A, mass?"""

        a_eliso = periodictable.z2el[z] + str(a)
        try:
            a_mass = periodictable.eliso2mass[a_eliso]
        except KeyError as err:
            raise ValidationError('Invalid nuclide: {}'.format(a_eliso))

        A_exact.append(a)
        A_range.append(lambda x: x == a)
        text.append("""For A, inputs Z: {}, A: {} require A == {}.""".format(z, a, a))

        m_exact.append(a_mass)
        m_range.append(lambda x: abs(x - a_mass) < mtol)
        text.append("""For mass, inputs Z: {}, A: {} require abs(mass - {}) < {}""".format(z, a, a_mass, mtol))

    def offer_mass(z, m):
        """Given a mass and element, what can be suggested and asserted about A, mass?"""

        m_a = int(round(m, 0))
        m_eliso = periodictable.z2el[z] + str(m_a)

        try:
            if abs(periodictable.eliso2mass[m_eliso] - m) > mtol:
                # only offer A if known nuclide. C@12.4 != 12C
                m_a = '(none)'
        except KeyError:
            m_a = '(none)'

        A_exact.append(m_a)
        A_range.append(lambda x: x == m_a)
        text.append("""For A, inputs Z: {}, m: {} require A == {}""".format(z, m, m_a))

        m_exact.append(m)
        m_range.append(lambda x: x == m)
        text.append("""For mass, inputs Z: {}, m: {} require m == {}""".format(z, m, m))

    # initialize

    text = ['', """--> Inp: A={}, Z={}, E={}, mass={}""".format(A, Z, E, mass)]

    Z_exact = []  # *_exact are candidates for the final value
    Z_range = []  # *_range are tests that the final value must pass to be valid
    A_exact = []
    A_range = []
    m_exact = []
    m_range = []

    # collect evidence for Z/A/m, then reconcile Z

    if Z is not None:
        _Z = int(Z)
        if _Z not in periodictable.z2el:
            raise ValidationError('Invalid element: {}'.format(Z))
        offer_atomic_number(_Z)

    if E is not None:
        try:
            _Z = periodictable.el2z[E.upper()]
        except KeyError as err:
            raise ValidationError('Invalid element: {}'.format(E))
        offer_atomic_number(_Z)

    Z_final = reconcile(Z_exact, Z_range, 'atomic number')
    E_final = periodictable.z2el[Z_final].capitalize()

    # collect more evidence for A/m, then reconcile them

    if A is not None:
        offer_mass_number(Z_final, int(A))

    if mass is not None:
        offer_mass(Z_final, mass * 1.)

    if label is not None:
        pass
        # TODO label not yet used
    # TODO ghosting not yet used

    mass_final = reconcile(m_exact, m_range, 'mass')
    A_final = reconcile(A_exact, A_range, 'mass number')

    text.append("""<-- Out: A={}, Z={}, E={}, mass={}""".format(A_final, Z_final, E_final, mass_final))

    if verbose >= 2:
        print('\n'.join(text))

    return (A_final, Z_final, E_final, mass_final, E_final)


if __name__ == "__main__":

    import qcdb
    from qcdb import periodictable
    from qcdb.exceptions import *

    co_dominant = (59, 27, 'Co', 58.933195048, 'Co')
    assert(co_dominant == reconcile_nucleus(E='co'))
    assert(co_dominant == reconcile_nucleus(Z=27))
    assert(co_dominant == reconcile_nucleus(A=59, Z=27))
    assert(co_dominant == reconcile_nucleus(E='cO', mass=58.933195048))
    assert(co_dominant == reconcile_nucleus(A=59, Z=27, E='CO'))
    assert(co_dominant == reconcile_nucleus(A=59, E='cO', mass=58.933195048))

    co_dominant_shortmass = (59, 27, 'Co', 58.933, 'Co')
    assert(co_dominant_shortmass == reconcile_nucleus(E='cO', mass=58.933))
    try:
        assert(co_dominant_shortmass == reconcile_nucleus(E='cO', mass=58.933, mtol=1.e-4))
    except AssertionError:
        pass

    co60 = (60, 27, 'Co', 59.933817059, 'Co')
    assert(co60 == reconcile_nucleus(E='Co', A=60))
    assert(co60 == reconcile_nucleus(Z=27, A=60))
    assert(co60 == reconcile_nucleus(E='Co', A=60))
    assert(co60 == reconcile_nucleus(Z=27, mass=59.933817059))
    assert(co60 == reconcile_nucleus(A=60, Z=27, mass=59.933817059))

    co_unspecified = ('(none)', 27, 'Co', 60.6, 'Co')
    assert(co_unspecified == reconcile_nucleus(mass=60.6, Z=27))
    assert(co_unspecified == reconcile_nucleus(mass=60.6, E='Co'))
    try:
        assert(co_unspecified == reconcile_nucleus(mass=60.6, Z=27, A=61))
    except ValidationError:
        pass

    try:
        reconcile_nucleus(A=80, Z=27)
    except ValidationError:
        pass

    try:
        reconcile_nucleus(Z=27, mass=200)
    except ValidationError:
        pass

    reconcile_nucleus(Z=27, mass=200, nonphysical=True)

    try:
        reconcile_nucleus(Z=27, mass=-200, nonphysical=True)
    except ValidationError:
        pass

    try:
        reconcile_nucleus(Z=-27, mass=200, nonphysical=True)
    except ValidationError:
        pass

    print('All reconcile_atom tests pass')
