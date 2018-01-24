import re

import numpy as np

from .. import periodictable
from ..exceptions import *
from . import regex


def reconcile_nucleus(A=None,
                      Z=None,
                      E=None,
                      mass=None,
                      real=None,
                      label=None,
                      nonphysical=False,
                      mtol=1.e-3,
                      verbose=1):
    """Considers all nuclear identity information available from
    arguments, supplemented by the periodic table. At the least, must
    provide element identity somehow. Returns consistent set of nucleus
    descriptors, defaulting to most-abundant isotope.

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
    real : bool, optional
        Whether real or ghost/absent.
    label : str, optional
        Atom label according to :py:attribute:`qcdb.molparse.NUCLEUS`.
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
    A, Z, E, mass, real, userlabel : int, int, str, float, bool, str
        mass number, unless clues don't point to a known nuclide, in which case `-1`.
        atomic number.
        element symbol, capitalized.
        mass value [u].
        real/ghost.
        user portion of `label` if present, else ''.

    """

    # <<< define functions

    def reconcile(exact, tests, feature):
        """Returns a member from `exact` that passes all `tests`, else raises error for `feature`."""

        for candidate in exact:
            assessment = [fn(candidate) for fn in tests]
            text.append(
                """Assess {} candidate {}: {} --> {}""".format(feature, candidate, assessment, all(assessment)))
            if all(assessment):
                return candidate
        else:
            err = """Inconsistent or unspecified {}: A: {}, Z: {}, E: {}, mass: {}, real: {}, label: {}""".format(
                feature, A, Z, E, mass, real, label)
            if verbose > -1:
                print('\n\n' + '\n'.join(text))
            raise ValidationError(err)

    def offer_element_symbol(e):
        """Given an element, what can be suggested and asserted about Z, A, mass?"""

        try:
            _Z = periodictable.el2z[e.upper()]
        except KeyError as err:
            raise ValidationError('Invalid element: {}'.format(e))
        offer_atomic_number(_Z)

    def offer_atomic_number(z):
        """Given an atomic number, what can be suggested and asserted about Z, A, mass?"""

        z = int(z)
        if z not in periodictable.z2el:
            raise ValidationError('Invalid element: {}'.format(z))

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
        Z_range.append(lambda x, z=z: x == z)

        A_exact.append(z_a)
        if nonphysical:
            A_range.append(lambda x: x == -1 or x >= 1)
            text.append("""For A, input Z: {}, requires 1 < A or -1, nonphysical""".format(z))
        else:
            A_range.append(lambda x, amin=z_a2mass_min, amax=z_a2mass_max: x == -1 or (x >= amin and x <= amax))
            text.append("""For A, input Z: {} requires {} < A < {} or -1, the known mass number range for element""".
                        format(z, z_a2mass_min, z_a2mass_max))

        m_exact.append(z_mass)
        if nonphysical:
            m_range.append(lambda x: x > 0.5)
            text.append("""For mass, input Z: {}, requires 0.5 < mass, nonphysical""".format(z))
        else:
            m_range.append(lambda x, mmin=z_mass2a_min, mmax=z_mass2a_max: x >= mmin - mmtol and x <= mmax + mmtol)
            text.append("""For mass, input Z: {} requires {} < mass < {} +/-{}, the known mass range for element""".
                        format(z, z_mass2a_min, z_mass2a_max, mmtol))

    def offer_mass_number(z, a):
        """Given a mass number and element, what can be suggested and asserted about A, mass?"""

        a = int(a)
        a_eliso = periodictable.z2el[z] + str(a)
        try:
            a_mass = periodictable.eliso2mass[a_eliso]
        except KeyError as err:
            raise ValidationError('Invalid nuclide: {}'.format(a_eliso))

        A_exact.append(a)
        A_range.append(lambda x, a=a: x == a)
        text.append("""For A, inputs Z: {}, A: {} require A == {}.""".format(z, a, a))

        m_exact.append(a_mass)
        m_range.append(lambda x, a_mass=a_mass: abs(x - a_mass) < mtol)
        text.append("""For mass, inputs Z: {}, A: {} require abs(mass - {}) < {}""".format(z, a, a_mass, mtol))

    def offer_mass_value(z, m):
        """Given a mass and element, what can be suggested and asserted about A, mass?"""

        m = float(m)
        m_a = int(round(m, 0))
        m_eliso = periodictable.z2el[z] + str(m_a)

        try:
            if abs(periodictable.eliso2mass[m_eliso] - m) > mtol:
                # only offer A if known nuclide. C@12.4 != 12C
                m_a = -1
        except KeyError:
            m_a = -1

        A_exact.append(m_a)
        A_range.append(lambda x, m_a=m_a: x == m_a)
        text.append("""For A, inputs Z: {}, m: {} require A == {}""".format(z, m, m_a))

        m_exact.append(m)
        m_range.append(lambda x, m=m: x == m)
        text.append("""For mass, inputs Z: {}, m: {} require m == {}""".format(z, m, m))

    def offer_reality(rgh):
        r_exact.append(rgh)
        r_range.append(lambda x, rgh=rgh: x == rgh)
        text.append("""For real/ghost, input rgh: {} requires rgh == {}""".format(rgh, rgh))

    def offer_user_label(lbl):
        lbl = lbl.lower()
        l_exact.append(lbl)
        l_range.append(lambda x, lbl=lbl: x == lbl)
        text.append("""For user label, input lbl: {} requires lbl == {}""".format(lbl, lbl))

    # <<< initialize

    text = ['', """--> Inp: A={}, Z={}, E={}, mass={}, real={}, label={}""".format(A, Z, E, mass, real, label)]

    Z_exact = []  # *_exact are candidates for the final value
    Z_range = []  # *_range are tests that the final value must pass to be valid
    A_exact = []
    A_range = []
    m_exact = []
    m_range = []

    r_exact = [True]  # default real/ghost is real
    r_range = []
    l_exact = ['']  # default user label is empty string
    l_range = []
    mmtol = 0.5  # tolerance for mass outside known masses for element

    # <<< collect evidence for Z/A/m, then reconcile Z

    if Z is not None:
        offer_atomic_number(Z)

    if E is not None:
        offer_element_symbol(E)

    if label is not None:
        lbl_A, lbl_Z, lbl_E, lbl_mass, lbl_real, lbl_user = parse_nucleus_label(label)

        if lbl_Z is not None:
            offer_atomic_number(lbl_Z)
        if lbl_E is not None:
            offer_element_symbol(lbl_E)

    Z_final = reconcile(Z_exact, Z_range, 'atomic number')
    E_final = periodictable.z2el[Z_final].capitalize()

    # <<< collect more evidence for A/m, then reconcile them

    if A is not None:
        offer_mass_number(Z_final, A)

    if mass is not None:
        offer_mass_value(Z_final, mass)

    if real is not None:
        offer_reality(real)

    if label is not None:
        offer_reality(lbl_real)
        if lbl_A is not None:
            offer_mass_number(Z_final, lbl_A)
        if lbl_mass is not None:
            offer_mass_value(Z_final, lbl_mass)
        if lbl_user is not None:
            offer_user_label(lbl_user)

    mass_final = reconcile(m_exact, m_range, 'mass')
    A_final = reconcile(A_exact, A_range, 'mass number')
    real_final = reconcile(r_exact, r_range, 'real/ghost')
    user_final = reconcile(l_exact, l_range, 'user label')

    text.append("""<-- Out: A={}, Z={}, E={}, mass={}, real={}, user={}""".format(A_final, Z_final, E_final,
                                                                                  mass_final, real_final, user_final))

    if verbose >= 2:
        print('\n'.join(text))

    return (A_final, Z_final, E_final, mass_final, real_final, user_final)


def parse_nucleus_label(label):
    """Separate the contents of molecule nucleus string into fields.

    Parameters
    ----------
    label : str
        Conveys at least element and ghostedness and possibly isotope, mass, and
        user info in accordance with :py:attribute:`qcdb.molparse.NUCLEUS`.

    Returns
    -------
    A, Z, E, mass, real, user : int or None, int or None, str or None, float or None, bool, str or None
        Field breakdown of `label`.

    Raises
    ------
    ValidationError
        If `label` does not match NUCLEUS.

    """
    nucleus = re.compile(r'\A' + regex.NUCLEUS + r'\Z', re.IGNORECASE | re.VERBOSE)
    matchobj = nucleus.match(label)

    if matchobj:
        real = not (matchobj.group('gh1') or matchobj.group('gh2'))

        if matchobj.group('A'):
            A = int(matchobj.group('A'))
        else:
            A = None

        if matchobj.group('Z'):
            Z = int(matchobj.group('Z'))
        else:
            Z = None

        E = matchobj.group('E')

        if matchobj.group('user1'):
            user = matchobj.group('user1')
        elif matchobj.group('user2'):
            user = matchobj.group('user2')
        else:
            user = None

        if matchobj.group('mass'):
            mass = float(matchobj.group('mass'))
        else:
            mass = None
    else:
        raise ValidationError("""Nucleus label is not parseable: {}""".format(label))

    return A, Z, E, mass, real, user


if __name__ == "__main__":

    # NOTE: to run tests as python module, comment out relative imports at top
    import qcdb
    from qcdb import periodictable
    from qcdb.exceptions import *
    from qcdb.molparse import regex

    co_dominant = (59, 27, 'Co', 58.933195048, True, '')
    co_dominant_mine = (59, 27, 'Co', 58.933195048, True, '_mine')
    assert (co_dominant == reconcile_nucleus(E='co'))
    assert (co_dominant == reconcile_nucleus(Z=27))
    assert (co_dominant == reconcile_nucleus(A=59, Z=27))
    assert (co_dominant == reconcile_nucleus(E='cO', mass=58.933195048))
    assert (co_dominant == reconcile_nucleus(A=59, Z=27, E='CO'))
    assert (co_dominant == reconcile_nucleus(A=59, E='cO', mass=58.933195048))
    assert (co_dominant == reconcile_nucleus(label='co'))
    assert (co_dominant == reconcile_nucleus(label='59co'))
    assert (co_dominant == reconcile_nucleus(label='co@58.933195048'))
    assert (co_dominant == reconcile_nucleus(A=59, Z=27, E='cO', mass=58.933195048, label='co@58.933195048'))
    assert (co_dominant == reconcile_nucleus(A=59, Z=27, E='cO', mass=58.933195048, label='27@58.933195048'))
    assert (co_dominant == reconcile_nucleus(label='27'))
    assert (co_dominant_mine == reconcile_nucleus(label='co_miNe'))
    assert (co_dominant_mine == reconcile_nucleus(label='co_mIne@58.933195048'))

    co_dominant_shortmass = (59, 27, 'Co', 58.933, True, '')
    assert (co_dominant_shortmass == reconcile_nucleus(E='cO', mass=58.933))
    assert (co_dominant_shortmass == reconcile_nucleus(label='cO@58.933'))
    try:
        assert (co_dominant_shortmass == reconcile_nucleus(E='cO', mass=58.933, mtol=1.e-4))
    except AssertionError:
        pass
    try:
        assert (co_dominant_shortmass == reconcile_nucleus(label='27@58.933', mtol=1.e-4))
    except AssertionError:
        pass

    co60 = (60, 27, 'Co', 59.933817059, True, '')
    assert (co60 == reconcile_nucleus(E='Co', A=60))
    assert (co60 == reconcile_nucleus(Z=27, A=60, real=True))
    assert (co60 == reconcile_nucleus(E='Co', A=60))
    assert (co60 == reconcile_nucleus(Z=27, mass=59.933817059))
    assert (co60 == reconcile_nucleus(A=60, Z=27, mass=59.933817059))
    assert (co60 == reconcile_nucleus(label='60Co'))
    assert (co60 == reconcile_nucleus(label='27', mass=59.933817059))
    assert (co60 == reconcile_nucleus(label='Co', mass=59.933817059))
    assert (co60 == reconcile_nucleus(A=60, label='Co'))

    co60ghost = (60, 27, 'Co', 59.933817059, False, '')
    assert (co60ghost == reconcile_nucleus(E='Co', A=60, real=False))
    assert (co60ghost == reconcile_nucleus(A=60, Z=27, mass=59.933817059, real=0))
    assert (co60ghost == reconcile_nucleus(label='@60Co'))
    assert (co60ghost == reconcile_nucleus(label='Gh(27)', mass=59.933817059))
    assert (co60ghost == reconcile_nucleus(label='@Co', mass=59.933817059))
    assert (co60ghost == reconcile_nucleus(A=60, label='Gh(Co)'))

    co_unspecified = (-1, 27, 'Co', 60.6, True, '')
    assert (co_unspecified == reconcile_nucleus(mass=60.6, Z=27))
    assert (co_unspecified == reconcile_nucleus(mass=60.6, E='Co'))
    assert (co_unspecified == reconcile_nucleus(mass=60.6, label='27'))
    assert (co_unspecified == reconcile_nucleus(label='Co@60.6'))
    try:
        assert (co_unspecified == reconcile_nucleus(mass=60.6, Z=27, A=61))
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

    try:
        reconcile_nucleus(Z=1, label='he')
    except ValidationError:
        pass

    try:
        reconcile_nucleus(A=4, label='3he')
    except ValidationError:
        pass

    try:
        reconcile_nucleus(label='@U', real=True)
    except ValidationError:
        pass

    try:
        reconcile_nucleus(label='U', real=False)
    except ValidationError:
        pass

    print('All reconcile_nucleus tests pass')

    label_tests = {
        '@ca_miNe': {
            'E': 'ca',
            'Z': None,
            'user': '_miNe',
            'A': None,
            'real': False,
            'mass': None
        },
        'Gh(Ca_mine)': {
            'E': 'Ca',
            'Z': None,
            'user': '_mine',
            'A': None,
            'real': False,
            'mass': None
        },
        '@Ca_mine@1.07': {
            'E': 'Ca',
            'Z': None,
            'user': '_mine',
            'A': None,
            'real': False,
            'mass': 1.07
        },
        'Gh(cA_MINE@1.07)': {
            'E': 'cA',
            'Z': None,
            'user': '_MINE',
            'A': None,
            'real': False,
            'mass': 1.07
        },
        '@40Ca_mine@1.07': {
            'E': 'Ca',
            'Z': None,
            'user': '_mine',
            'A': 40,
            'real': False,
            'mass': 1.07
        },
        'Gh(40Ca_mine@1.07)': {
            'E': 'Ca',
            'Z': None,
            'user': '_mine',
            'A': 40,
            'real': False,
            'mass': 1.07
        },
        '444lu333@4.0': {
            'E': 'lu',
            'Z': None,
            'user': '333',
            'A': 444,
            'real': True,
            'mass': 4.0
        },
        '@444lu333@4.4': {
            'E': 'lu',
            'Z': None,
            'user': '333',
            'A': 444,
            'real': False,
            'mass': 4.4
        },
        '8i': {
            'E': 'i',
            'Z': None,
            'user': None,
            'A': 8,
            'real': True,
            'mass': None
        },
        '53_mI4': {
            'Z': 53,
            'E': None,
            'user': '_mI4',
            'A': None,
            'real': True,
            'mass': None
        },
        '@5_MINEs3@4.4': {
            'Z': 5,
            'E': None,
            'user': '_MINEs3',
            'A': None,
            'real': False,
            'mass': 4.4
        },
        'Gh(555_mines3@0.1)': {
            'Z': 555,
            'E': None,
            'user': '_mines3',
            'A': None,
            'real': False,
            'mass': 0.1
        },
    }

    for test, ans in label_tests.items():
        lbl_A, lbl_Z, lbl_E, lbl_mass, lbl_real, lbl_user = parse_nucleus_label(test)

        qcdb.compare_integers(ans['real'], lbl_real, test + " real")
        qcdb.compare_integers(ans['A'], lbl_A, test + " A")
        qcdb.compare_integers(ans['Z'], lbl_Z, test + " Z")
        qcdb.compare_strings(ans['E'], lbl_E, test + " symbol")
        qcdb.compare_strings(ans['user'], lbl_user, test + " user")
        qcdb.compare_values(ans['mass'], lbl_mass, 6, test + " mass", passnone=True)
    print('All parse_nucleus_label tests pass')
