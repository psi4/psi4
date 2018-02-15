from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import sys
import math
import collections

import numpy as np

from .physconst import *
from .psiutil import *
from .util import *
from .libmintsmolecule import compute_atom_map

try:
    from itertools import izip_longest as zip_longest  # py2
except ImportError:
    from itertools import zip_longest  # py3

VibrationAspect = collections.namedtuple('VibrationAspect', 'lbl unit data comment')

LINEAR_A_TOL = 1.0E-2  # tolerance (roughly max dev) for TR space


def compare_vibinfos(expected, computed, tol, label, verbose=1, forgive=None, required=None):
    """Compare two dictionaries of VibrationAspect objects. All items in
    `expected` must be present in `computed` and agree to `tol` (outright
    (<1) or decimal places (>=1)). At minimum, both must contain items in
    `required` ('omega' recc. for vibs).  Failed items in `forgive` will
    not count against passing.

    """
    np.set_printoptions(formatter={'float': '{: 0.4f}'.format})

    def _success(label):
        """Function to print a '*label*...PASSED' line to screen.
        Used by :py:func:`util.compare_values` family when functions pass.
        """
        msg = '\t{0:.<66}PASSED'.format(label)
        print(msg)
        sys.stdout.flush()

    def print_stuff(asp, same, ref, val, space=''):
        if verbose >= 1:
            print(asp, ':', same)
        if (verbose >= 2) or (not same and verbose >= 1):
            print('\texp:', space, ref)
            print('\tobs:', space, val)
        if verbose >= 1:
            if not same:
                try:
                    print('\tdif:', space, val - ref)
                except TypeError:
                    print('\tdif: Different, inspect arrays')

    if tol >= 1.0:
        tol = 10.**-tol

    if forgive is None:
        forgive = []

    summsame = []
    if required is None:
        checkkeys = []
    else:
        checkkeys = required
    checkkeys.extend(expected.keys())

    for asp in checkkeys:
        if asp in ['q', 'w', 'x']:
            ccnc = _phase_cols_to_max_element(computed[asp].data)
            eenc = _phase_cols_to_max_element(expected[asp].data)
            ccnc = _check_degen_modes(ccnc, computed['omega'].data)
            eenc = _check_degen_modes(eenc, expected['omega'].data)
            same = np.allclose(eenc, ccnc, atol=tol)
            print_stuff(asp=asp, same=same, ref=eenc, val=ccnc, space='\n')

        elif asp in ['gamma', 'TRV']:
            same = all([computed[asp].data[idx] == val for idx, val in enumerate(expected[asp].data)])
            print_stuff(asp=asp, same=same, ref=expected[asp].data, val=computed[asp].data)

        elif isinstance(expected[asp].data, float):
            same = abs(expected[asp].data - computed[asp].data) < tol
            print_stuff(asp=asp, same=same, ref=expected[asp].data, val=computed[asp].data)
        else:
            same = np.allclose(expected[asp].data, computed[asp].data, atol=tol) and \
                   (expected[asp].data.shape == computed[asp].data.shape)
            print_stuff(asp=asp, same=same, ref=expected[asp].data, val=computed[asp].data)

        if asp not in forgive:
            summsame.append(same)

    passed = all(summsame)
    if passed:
        _success(label)
    return passed


def hessian_symmetrize(hess, mol):
    """Apply Abelian symmetry of `mol` to Hessian `hess`.

    Parameters
    ----------
    hess : ndarray of float
        (3 * nat, 3 * nat) Hessian array perhaps with jitter unbecoming a symmetric molecule.
    mol : qcdb.Molecule
        Molecule at which Hessian computed.

    Returns
    -------
    ndarray
        (3 * nat, 3 * nat) symmetrized Hessian array.

    """
    ct = mol.point_group().char_table()

    # Obtain atom mapping of atom * symm op to atom
    atom_map = compute_atom_map(mol)

    syms = []
    smap = []
    for g in range(ct.order()):
        syms.append(np.asarray(ct.symm_operation(g).d))
        smap.append([atom_map[at][g] for at in range(mol.natom())])

    np.set_printoptions(formatter={'float': '{: 16.12f}'.format})
    b_hess = blockwise_expand(hess, (3, 3), False)

    bDG = []
    nat = b_hess.shape[0]
    for iat in range(nat):
        for jat in range(nat):
            for sym in range(len(syms)):
                bDG.append(np.zeros_like(b_hess))
                bDG[sym][iat, jat] = syms[sym].dot(b_hess[iat, jat].dot(syms[sym]))
                # Note that tested syms all diagonal, so above may be off by some transposes

    for sym in range(len(syms)):
        bDG[sym] = bDG[sym][:, smap[sym]]
        bDG[sym] = bDG[sym][smap[sym], :]
    tot = np.sum(bDG, axis=0)
    tot = np.divide(tot, len(syms))

    print('symmetrization diff:', np.linalg.norm(tot - b_hess))
    m_tot = blockwise_contract(tot)
    return m_tot


def print_molden_vibs(vibinfo, atom_symbol, geom, standalone=True):
    """Format vibrational analysis for Molden.

    Parameters
    ----------
    vibinfo : dict of VibrationAspect
        Holds results of vibrational analysis.
    atom_symbol : array-like of str
        (nat,) element symbols for geometry of vibrational analysis.
    geom : array-like of float
        (nat, 3) geometry of vibrational analysis [a0].
    standalone : bool, optional
        Whether returned string prefixed "[Molden Format]" for standalone rather than append.

    Returns
    -------
    str
        `vibinfo` formatted for Molden, including FREQ, FR-COORD, & FR-NORM-COORD fields.

    Notes
    -----
    Molden format spec from http://www.cmbi.ru.nl/molden/molden_format.html
    Specifies "atomic coordinates x,y,z and atomic displacements dx,dy,dz are all in Bohr (Atomic Unit of length)"

    Despite it being quite wrong, imaginary modes are represented by a negative frequency.

    """
    nat = int(len(vibinfo['q'].data[:, 0]) / 3)
    active = [idx for idx, trv in enumerate(vibinfo['TRV'].data) if trv == 'V']

    text = ''
    if standalone:
        text += """[Molden Format]\n"""

    text += """\n[FREQ]\n"""
    for vib in active:
        if vibinfo['omega'].data[vib].imag > vibinfo['omega'].data[vib].real:
            freq = -1.0 * vibinfo['omega'].data[vib].imag
        else:
            freq = vibinfo['omega'].data[vib].real
        text += """   {:20.10f}\n""".format(freq)

    text += """\n[FR-COORD]\n"""
    for at in range(nat):
        text += ("""{:3}""" + """{:20.10f}""" * 3 + '\n').format(atom_symbol[at], *geom[at])

    text += """\n[FR-NORM-COORD]\n"""
    for idx, vib in enumerate(active):
        text += """vibration {}\n""".format(idx + 1)
        for at in range(nat):
            text += ('   ' + """{:20.10f}""" * 3 + '\n').format(*(vibinfo['x'].data[:, vib].reshape(nat, 3)[at].real))


#     text += """\n[INT]\n"""
#     for vib in active:
#         text += """1.0\n"""

    return text


def _check_degen_modes(arr, freq, verbose=1):
    """Use `freq` to identify degenerate columns of eigenvectors `arr` and
    sort into std order for comparison. Returns eigenvectors back sorted.

    """
    arr2 = np.zeros_like(arr)
    dfreq, didx, dinv, dcts = np.unique(np.around(freq, 2), return_index=True, return_inverse=True, return_counts=True)

    # judging degen normco to only 2 decimals is probably sign need to resolve evec
    idx_max_elem_each_normco = np.argmax(np.around(arr, 2), axis=0)
    max_elem_each_normco = np.amax(np.around(arr, 2), axis=0)
    idx_vib_reordering = np.empty_like(idx_max_elem_each_normco)

    for idegen, istart in enumerate(didx):
        degree = dcts[idegen]

        # sort degen evec
        #   primarily (last arg) by value of extreme element
        #             (sep evec that in this coord sys have diff elements)
        #   & secondarily (2nd-to-last arg) by index of extreme element
        #             (order evec with same elements in diff (xyz) arrangements)
        idx_sort_wi_degen = np.lexsort((idx_max_elem_each_normco[istart:istart + degree],
                                        max_elem_each_normco[istart:istart + degree]))
        idx_vib_reordering[istart:istart + degree] = np.arange(istart, istart + degree)[idx_sort_wi_degen]

    arr2 = arr[:, idx_vib_reordering]

    reorderings = ['{}-->{}'.format(i, v) for i, v in enumerate(idx_vib_reordering) if (i != v)]
    if reorderings and verbose >= 2:
        print('Degenerate modes reordered:', ', '.join(reorderings))

    return arr2


def _phase_cols_to_max_element(arr, tol=1.e-2, verbose=1):
    """Returns copy of 2D `arr` scaled such that, within cols, max(fabs)
    element is positive. If max(fabs) is pos/neg pair, scales so first
    element (within `tol`) is positive.

    """
    arr2 = np.copy(arr)

    rephasing = []
    for v in range(arr.shape[1]):
        vextreme = 0.0
        iextreme = None

        # find most extreme value
        for varr in arr[:, v]:
            vextreme = max(np.absolute(varr), vextreme)

        # find the first index whose fabs equals that value, w/i tolerance
        for iarr, varr in enumerate(arr[:, v]):
            if (vextreme - np.absolute(varr)) < tol:
                iextreme = iarr
                break

        sign = np.sign(arr[iextreme, v])
        if sign == -1.:
            rephasing.append(str(v))
        arr2[:, v] *= sign

    if rephasing and verbose >= 2:
        print('Negative modes rephased:', ', '.join(rephasing))

    return arr2


def harmonic_analysis(hess, geom, mass, basisset, irrep_labels, project_trans=True, project_rot=True):
    """Like so much other Psi4 goodness, originally by @andysim

    Parameters
    ----------
    hess : ndarray of float
        (3*nat, 3*nat) non-mass-weighted Hessian in atomic units, [Eh/a0/a0].
    geom : ndarray of float
        (nat, 3) geometry [a0] at which Hessian computed.
    mass : ndarray of float
        (nat,) atomic masses [u].
    basisset : psi4.core.BasisSet
        Basis set object (can be dummy, e.g., STO-3G) for SALCs.
    irrep_labels : list of str
        Irreducible representation labels.
    project_trans : bool, optional
        Idealized translations projected out of final vibrational analysis.
    project_rot : bool, optional
        Idealized rotations projected out of final vibrational analysis.

    Returns
    -------
    dict, text
        Returns dictionary of VibrationAspect objects (fields: lbl unit data comment).
        Also returns text suitable for printing.

    .. _`table:vibaspectinfo`:

    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | key           | description (lbl & comment)                | units     | data (real/imaginary modes)                          |
    +===============+============================================+===========+======================================================+
    | omega         | frequency                                  | cm^-1     | nd.array(ndof) complex (real/imag)                   |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | q             | normal mode, normalized mass-weighted      | a0 u^1/2  | ndarray(ndof, ndof) float                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | w             | normal mode, un-mass-weighted              | a0        | ndarray(ndof, ndof) float                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | x             | normal mode, normalized un-mass-weighted   | a0        | ndarray(ndof, ndof) float                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | degeneracy    | degree of degeneracy                       |           | ndarray(ndof) int                                    |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | TRV           | translation/rotation/vibration             |           | ndarray(ndof) str 'TR' or 'V' or '-' for partial     |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | gamma         | irreducible representation                 |           | ndarray(ndof) str irrep or None if unclassifiable    |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | mu            | reduced mass                               | u         | ndarray(ndof) float (+/+)                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | k             | force constant                             | mDyne/A   | ndarray(ndof) float (+/-)                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | DQ0           | RMS deviation v=0                          | a0 u^1/2  | ndarray(ndof) float (+/0)                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | Qtp0          | Turning point v=0                          | a0 u^1/2  | ndarray(ndof) float (+/0)                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | Xtp0          | Turning point v=0                          | a0        | ndarray(ndof) float (+/0)                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | theta_vib     | char temp                                  | K         | ndarray(ndof) float (+/0)                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+

    Examples
    --------
    # displacement of first atom in highest energy mode
    >>> vibinfo['x'].data[:, -1].reshape(nat, 3)[0]

    # remove translations & rotations
    >>> vibonly = filter_nonvib(vibinfo)

    """
    from psi4 import core

    if (mass.shape[0] == geom.shape[0] == (hess.shape[0] // 3) == (hess.shape[1] // 3)) and (geom.shape[1] == 3):
        pass
    else:
        raise ValidationError("""Dimension mismatch among mass ({}), geometry ({}), and Hessian ({})""".format(
            mass.shape, geom.shape, hess.shape))

    def mat_symm_info(a, atol=1e-14, lbl='array', stol=None):
        symm = np.allclose(a, a.T, atol=atol)
        herm = np.allclose(a, a.conj().T, atol=atol)
        ivrt = a.shape[0] - np.linalg.matrix_rank(a, tol=stol)
        return """  {:32} Symmetric? {}   Hermitian? {}   Lin Dep Dim? {:2}""".format(lbl + ':', symm, herm, ivrt)

    def vec_in_space(vec, space, tol=1.0e-4):
        merged = np.vstack((space, vec))
        u, s, v = np.linalg.svd(merged)
        return (s[-1] < tol)

    vibinfo = {}
    text = []

    nat = len(mass)
    text.append("""\n\n  ==> Harmonic Vibrational Analysis <==\n""")

    if nat == 1:
        nrt_expected = 3
    elif np.linalg.matrix_rank(geom) == 1:
        nrt_expected = 5
    else:
        nrt_expected = 6

    nmwhess = hess.copy()
    text.append(mat_symm_info(nmwhess, lbl='non-mass-weighted Hessian') + ' (0)')

    # get SALC object, possibly w/o trans & rot
    mints = core.MintsHelper(basisset)
    cdsalcs = mints.cdsalcs(0xFF, project_trans, project_rot)

    Uh = collections.OrderedDict()
    for h, lbl in enumerate(irrep_labels):
        tmp = np.asarray(cdsalcs.matrix_irrep(h))
        if tmp.size > 0:
            Uh[lbl] = tmp

    # form projector of translations and rotations
    space = ('T' if project_trans else '') + ('R' if project_rot else '')
    TRspace = _get_TR_space(mass, geom, space=space, tol=LINEAR_A_TOL)
    nrt = TRspace.shape[0]
    text.append(
        '  projection of translations ({}) and rotations ({}) removed {} degrees of freedom ({})'.
        format(project_trans, project_rot, nrt, nrt_expected))

    P = np.identity(3 * nat)
    for irt in TRspace:
        P -= np.outer(irt, irt)
    text.append(mat_symm_info(P, lbl='total projector') + ' ({})'.format(nrt))

    # mass-weight & solve
    sqrtmmm = np.repeat(np.sqrt(mass), 3)
    sqrtmmminv = np.divide(1.0, sqrtmmm)
    mwhess = np.einsum('i,ij,j->ij', sqrtmmminv, nmwhess, sqrtmmminv)
    text.append(mat_symm_info(mwhess, lbl='mass-weighted Hessian') + ' (0)')

    pre_force_constant_au = np.linalg.eigvalsh(mwhess)

    idx = np.argsort(pre_force_constant_au)
    pre_force_constant_au = pre_force_constant_au[idx]
    uconv_cm_1 = np.sqrt(psi_na * psi_hartree2J * 1.0e19) / (2 * np.pi * psi_c * psi_bohr2angstroms)
    pre_frequency_cm_1 = np.lib.scimath.sqrt(pre_force_constant_au) * uconv_cm_1

    pre_lowfreq = np.where(np.real(pre_frequency_cm_1) < 100.0)[0]
    pre_lowfreq = np.append(pre_lowfreq, np.arange(nrt_expected))  # catch at least nrt modes
    for lf in set(pre_lowfreq):
        vlf = pre_frequency_cm_1[lf]
        if vlf.imag > vlf.real:
            text.append('  pre-proj  low-frequency mode: {:9.4f}i [cm^-1]'.format(vlf.real, vlf.imag))
        else:
            text.append('  pre-proj  low-frequency mode: {:9.4f}  [cm^-1]'.format(vlf.real, ''))
    text.append('  pre-proj  all modes:' + str(_format_omega(pre_frequency_cm_1, 4)))

    # project & solve
    mwhess_proj = np.dot(P.T, mwhess).dot(P)
    text.append(mat_symm_info(mwhess_proj, lbl='projected mass-weighted Hessian') + ' ({})'.format(nrt))

    #print('projhess = ', np.array_repr(mwhess_proj))
    force_constant_au, qL = np.linalg.eigh(mwhess_proj)

    # expected order for vibrations is steepest downhill to steepest uphill
    idx = np.argsort(force_constant_au)
    force_constant_au = force_constant_au[idx]
    qL = qL[:, idx]
    qL = _phase_cols_to_max_element(qL)
    vibinfo['q'] = VibrationAspect('normal mode', '[a0 u^1/2]', qL, 'normalized mass-weighted')

    # frequency, LAB II.17
    frequency_cm_1 = np.lib.scimath.sqrt(force_constant_au) * uconv_cm_1
    vibinfo['omega'] = VibrationAspect('frequency', '[cm^-1]', frequency_cm_1, '')

    # degeneracies
    ufreq, uinv, ucts = np.unique(np.around(frequency_cm_1, 2), return_inverse=True, return_counts=True)
    vibinfo['degeneracy'] = VibrationAspect('degeneracy', '', ucts[uinv], '')

    # look among the symmetry subspaces h for one to which the normco
    #   of vib does *not* add an extra dof to the vector space
    active = []
    irrep_classification = []
    for idx, vib in enumerate(frequency_cm_1):

        if vec_in_space(qL[:, idx], TRspace, 1.0e-4):
            active.append('TR')
            irrep_classification.append(None)

        else:
            active.append('V')

            for h in Uh.keys():
                if vec_in_space(qL[:, idx], Uh[h], 1.0e-4):
                    irrep_classification.append(h)
                    break
            else:
                irrep_classification.append(None)

                # catch partial Hessians
                if np.linalg.norm(vib) < 1.e-3:
                    active[-1] = '-'

    vibinfo['TRV'] = VibrationAspect('translation/rotation/vibration', '', active, '')
    vibinfo['gamma'] = VibrationAspect('irreducible representation', '', irrep_classification, '')

    lowfreq = np.where(np.real(frequency_cm_1) < 100.0)[0]
    lowfreq = np.append(lowfreq, np.arange(nrt_expected))  # catch at least nrt modes
    for lf in set(lowfreq):
        vlf = frequency_cm_1[lf]
        if vlf.imag > vlf.real:
            text.append('  post-proj low-frequency mode: {:9.4f}i [cm^-1] ({})'.format(vlf.imag, active[lf]))
        else:
            text.append('  post-proj low-frequency mode: {:9.4f}  [cm^-1] ({})'.format(vlf.real, active[lf]))
    text.append('  post-proj  all modes:' + str(_format_omega(frequency_cm_1, 4)) + '\n')
    if project_trans and not project_rot:
        text.append('  Note that "Vibration"s include {} un-projected rotation-like modes.'.format(nrt_expected - 3))
    elif not project_trans and not project_rot:
        text.append('  Note that "Vibration"s include {} un-projected rotation-like and translation-like modes.'.format(nrt_expected))

    # general conversion factors, LAB II.11
    uconv_K = (psi_h * psi_na * 1.0e21) / (8 * np.pi * np.pi * psi_c)
    uconv_S = np.sqrt((psi_c * (2 * np.pi * psi_bohr2angstroms)**2) / (psi_h * psi_na * 1.0e21))

    # normco & reduced mass, LAB II.14 & II.15
    wL = np.einsum('i,ij->ij', sqrtmmminv, qL)
    vibinfo['w'] = VibrationAspect('normal mode', '[a0]', wL, 'un-mass-weighted')

    reduced_mass_u = np.divide(1.0, np.linalg.norm(wL, axis=0)**2)
    vibinfo['mu'] = VibrationAspect('reduced mass', '[u]', reduced_mass_u, '')

    xL = np.sqrt(reduced_mass_u) * wL
    vibinfo['x'] = VibrationAspect('normal mode', '[a0]', xL, 'normalized un-mass-weighted')

    # force constants, LAB II.16 (real compensates for earlier sqrt)
    uconv_mdyne_a = (0.1 * (2 * np.pi * psi_c)**2) / psi_na
    force_constant_mdyne_a = reduced_mass_u * (frequency_cm_1 * frequency_cm_1).real * uconv_mdyne_a
    vibinfo['k'] = VibrationAspect('force constant', '[mDyne/A]', force_constant_mdyne_a, '')

    force_constant_cm_1_bb = reduced_mass_u * (frequency_cm_1 * frequency_cm_1).real * uconv_S * uconv_S
    VibrationAspect('force constant', '[cm^-1/a0^2]', force_constant_cm_1_bb, "Hooke's Law")

    # turning points, LAB II.20 (real & zero since turning point silly for imag modes)
    nu = 0
    turning_point_rnc = np.sqrt(2.0 * nu + 1.0)

    with np.errstate(divide='ignore'):
        turning_point_bohr_u = turning_point_rnc / (np.sqrt(frequency_cm_1.real) * uconv_S)
    turning_point_bohr_u[turning_point_bohr_u == np.inf] = 0.
    vibinfo['Qtp0'] = VibrationAspect('Turning point v=0', '[a0 u^1/2]', turning_point_bohr_u, '')

    with np.errstate(divide='ignore'):
        turning_point_bohr = turning_point_rnc / (np.sqrt(frequency_cm_1.real * reduced_mass_u) * uconv_S)
    turning_point_bohr[turning_point_bohr == np.inf] = 0.
    vibinfo['Xtp0'] = VibrationAspect('Turning point v=0', '[a0]', turning_point_bohr, '')

    rms_deviation_bohr_u = turning_point_bohr_u / np.sqrt(2.0)
    vibinfo['DQ0'] = VibrationAspect('RMS deviation v=0', '[a0 u^1/2]', rms_deviation_bohr_u, '')

    # characteristic vibrational temperature, RAK thermo & https://en.wikipedia.org/wiki/Vibrational_temperature
    #   (imag freq zeroed)
    uconv_K = 100 * psi_h * psi_c / psi_kb
    vib_temperature_K = frequency_cm_1.real * uconv_K
    vibinfo['theta_vib'] = VibrationAspect('char temp', '[K]', vib_temperature_K, '')

    return vibinfo, '\n'.join(text)


def _format_omega(omega, decimals):
    """Return complex frequencies in `omega` into strings showing only real or imag ("i"-labeled)
    to `decimals` precision.

    """
    freqs = []
    for fr in omega:
        if fr.imag > fr.real:
            freqs.append("""{:.{prec}f}i""".format(fr.imag, prec=decimals))
        else:
            freqs.append("""{:.{prec}f}""".format(fr.real, prec=decimals))
    return np.array(freqs)


def print_vibs(vibinfo, atom_lbl=None, normco='x', shortlong=True, **kwargs):
    """Pretty printer for vibrational analysis.

    Parameters
    ----------
    vibinfo : dict of VibrationAspect
        Results of a Hessian solution.
    atom_lbl : list of str, optional
        Atomic symbols for printing. If None, integers used.
    normco : {'q', 'w', 'x'}, optional
        Which normal coordinate definition to print (reduced mass, etc. unaffected by this parameter).
        Must be `q` [a0 u^1/2], the mass-weighted normalized eigenvectors of the Hessian,
                `w` [a0], the un-mass-weighted (Cartesian) of q, or
                `x` [a0], the normalized w.
    shortlong : bool, optional
        Whether normal coordinates should be in (nat, 3) `True` or (nat * 3, 1) `False` format.
    groupby : int, optional
        How many normal coordinates per row. Defaults to 3/6 for shortlong=T/F. Value of `-1` uses one row.
    prec : int, optional
        Number of decimal places for frequencies, reduced masses, etc.
    ncprec : int, optional
        Number of decimal places for normal coordinates. Defaults to 2 for shortlong=short and 4 for shortlong=long.

    Returns
    -------
    str
        String suitable for printing.

    """

    def grouper(iterable, n, fillvalue=None):
        "Collect data into fixed-length chunks or blocks"
        # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
        args = [iter(iterable)] * n
        return zip_longest(*args, fillvalue=fillvalue)

    if normco not in ['q', 'w', 'x']:
        raise ValidationError("""Requested normal coordinates not among allowed q/w/x: """ + normco)

    nat = int(len(vibinfo['q'].data[:, 0]) / 3)
    if atom_lbl is None:
        atom_lbl = [''] * nat

    active = [idx for idx, trv in enumerate(vibinfo['TRV'].data) if trv == 'V']

    presp = 2
    colsp = 2
    prec = kwargs.get('prec', 4)
    if shortlong:
        groupby = kwargs.get('groupby', 3)
        ncprec = kwargs.get('ncprec', 2)
        width = (ncprec + 4) * 3
        prewidth = 24
    else:
        groupby = kwargs.get('groupby', 6)
        ncprec = kwargs.get('ncprec', 4)
        width = ncprec + 8
        prewidth = 24
    if groupby == -1:
        groupby = len(active)

    omega_str = _format_omega(vibinfo['omega'].data, decimals=prec)

    text = ''
    for row in grouper(active, groupby):

        text += """\n{:{presp}}{:{prewidth}}""".format('', 'Vibration', prewidth=prewidth, presp=presp)
        for vib in row:
            if vib is None:
                # ran out of vibrations in this row
                break
            text += """{:^{width}d}{:{colsp}}""".format(vib + 1, '', width=width, colsp=colsp)
        text += '\n'

        text += """{:{presp}}{:{prewidth}}""".format('', 'Freq [cm^-1]', prewidth=prewidth, presp=presp)
        for vib in row:
            if vib is None:
                break
            text += """{:^{width}}  """.format(omega_str[vib], width=width)
        text += '\n'

        text += """{:{presp}}{:{prewidth}}""".format('', 'Irrep', prewidth=prewidth, presp=presp)
        for vib in row:
            if vib is None:
                break
            val = vibinfo['gamma'].data[vib]
            if val is None:
                val = ''
            text += """{:^{width}}{:{colsp}}""".format(val, '', width=width, colsp=colsp)
        text += '\n'

        text += """{:{presp}}{:{prewidth}}""".format(
            '', 'Reduced mass ' + vibinfo['mu'].unit, prewidth=prewidth, presp=presp)
        for vib in row:
            if vib is None:
                break
            text += """{:^{width}.{prec}f}{:{colsp}}""".format(
                vibinfo['mu'].data[vib], '', width=width, prec=prec, colsp=colsp)
        text += '\n'

        text += """{:{presp}}{:{prewidth}}""".format(
            '', 'Force const ' + vibinfo['k'].unit, prewidth=prewidth, presp=presp)
        for vib in row:
            if vib is None:
                break
            text += """{:^{width}.{prec}f}{:{colsp}}""".format(
                vibinfo['k'].data[vib], '', width=width, prec=prec, colsp=colsp)
        text += '\n'

        text += """{:{presp}}{:{prewidth}}""".format(
            '', 'Turning point v=0 ' + vibinfo['Xtp0'].unit, prewidth=prewidth, presp=presp)
        for vib in row:
            if vib is None:
                break
            text += """{:^{width}.{prec}f}{:{colsp}}""".format(
                vibinfo['Xtp0'].data[vib], '', width=width, prec=prec, colsp=colsp)
        text += '\n'

        text += """{:{presp}}{:{prewidth}}""".format(
            '', 'RMS dev v=0 ' + vibinfo['DQ0'].unit, prewidth=prewidth, presp=presp)
        for vib in row:
            if vib is None:
                break
            text += """{:^{width}.{prec}f}{:{colsp}}""".format(
                vibinfo['DQ0'].data[vib], '', width=width, prec=prec, colsp=colsp)
        text += '\n'

        if 'theta_vib' in vibinfo:
            text += """{:{presp}}{:{prewidth}}""".format(
                '', 'Char temp ' + vibinfo['theta_vib'].unit, prewidth=prewidth, presp=presp)
            for vib in row:
                if vib is None:
                    break
                text += """{:^{width}.{prec}f}{:{colsp}}""".format(
                    vibinfo['theta_vib'].data[vib], '', width=width, prec=prec, colsp=colsp)
            text += '\n'

        #text += 'IR activ [KM/mol]\n'
        #text += 'Raman activ [A^4/u]\n'
        text += ' ' * presp + '-' * (prewidth + groupby * (width + colsp) - colsp) + '\n'

        if shortlong:
            for at in range(nat):
                text += """{:{presp}}{:5d}   {:{width}}""".format(
                    '', at + 1, atom_lbl[at], width=prewidth - 8, presp=presp)
                for vib in row:
                    if vib is None:
                        break
                    text += ("""{:^{width}.{prec}f}""" * 3).format(
                        *(vibinfo[normco].data[:, vib].reshape(nat, 3)[at]), width=int(width / 3), prec=ncprec)
                    text += """{:{colsp}}""".format('', colsp=colsp)
                text += '\n'
        else:
            for at in range(nat):
                for xyz in range(3):
                    text += """{:{presp}}{:5d}    {}    {:{width}}""".format(
                        '', at + 1, 'XYZ' [xyz], atom_lbl[at], width=prewidth - 14, presp=presp)
                    for vib in row:
                        if vib is None:
                            break
                        text += """{:^{width}.{prec}f}""".format(
                            (vibinfo[normco].data[3 * at + xyz, vib]), width=width, prec=ncprec)
                        text += """{:{colsp}}""".format('', colsp=colsp)
                    text += '\n'

    return text


def thermo(vibinfo, T, P, multiplicity, molecular_mass, E0, sigma, rot_const, rotor_type=None):
    """Perform thermochemical analysis from vibrational output.

    Parameters
    ----------
    E0 : float
        Electronic energy [Eh] at well bottom at 0 [K], :psivar:`CURRENT ENERGY`.
    molecular_mass : float
        Mass in [u] of molecule under analysis.
    multiplicity: int
        Spin multiplicity of molecule under analysis.
    rot_const : ndarray of floats
        (3,) rotational constants in [cm^-1] of molecule under analysis.
    sigma : int
        The rotational or external symmetry number determined from the point group.
    rotor_type : str
        The rotor type for rotational stat mech purposes: RT_ATOM, RT_LINEAR, other.
    T : float
        Temperature in [K]. Psi default 298.15. Note that 273.15 is IUPAC STP.
    P : float
        Pressure in [Pa]. Psi default 101325. Note that 100000 is IUPAC STP.

    Returns
    -------
    dict of VibrationAspect, str
        First is every thermochemistry component in atomic units.
        Second is formatted presentation of analysis.

    """
    sm = collections.defaultdict(float)

    # electronic
    q_elec = multiplicity
    sm[('S', 'elec')] = math.log(q_elec)

    # translational
    beta = 1 / (psi_kb * T)
    q_trans = (2.0 * np.pi * molecular_mass * psi_amu2kg / (beta * psi_h * psi_h))**1.5 * psi_na / (beta * P)
    sm[('S', 'trans')] = 5 / 2 + math.log(q_trans / psi_na)
    sm[('Cv', 'trans')] = 3 / 2
    sm[('Cp', 'trans')] = 5 / 2
    sm[('E', 'trans')] = 3 / 2 * T
    sm[('H', 'trans')] = 5 / 2 * T

    # rotational
    if rotor_type == "RT_ATOM":
        pass
    elif rotor_type == "RT_LINEAR":
        q_rot = 1. / (beta * sigma * 100 * psi_c * psi_h * rot_const[1])
        sm[('S', 'rot')] = 1.0 + math.log(q_rot)
        sm[('Cv', 'rot')] = 1
        sm[('Cp', 'rot')] = 1
        sm[('E', 'rot')] = T
    else:
        phi_A, phi_B, phi_C = rot_const * 100 * psi_c * psi_h / psi_kb
        q_rot = math.sqrt(math.pi) * T**1.5 / (sigma * math.sqrt(phi_A * phi_B * phi_C))
        sm[('S', 'rot')] = 3 / 2 + math.log(q_rot)
        sm[('Cv', 'rot')] = 3 / 2
        sm[('Cp', 'rot')] = 3 / 2
        sm[('E', 'rot')] = 3 / 2 * T
    sm[('H', 'rot')] = sm[('E', 'rot')]

    # vibrational
    vibonly = filter_nonvib(vibinfo)
    ZPE_cm_1 = 1 / 2 * np.sum(vibonly['omega'].data.real)
    omega_str = _format_omega(vibonly['omega'].data, decimals=4)

    imagfreqidx = np.where(vibonly['omega'].data.imag > vibonly['omega'].data.real)[0]
    if len(imagfreqidx):
        print("Warning: thermodynamics relations excluded imaginary frequencies: {}".format(omega_str[imagfreqidx]))

    filtered_theta_vib = np.delete(vibonly['theta_vib'].data, imagfreqidx, None)
    filtered_omega_str = np.delete(omega_str, imagfreqidx, None)
    rT = filtered_theta_vib / T  # reduced temperature

    lowfreqidx = np.where(filtered_theta_vib < 900.)[0]
    if len(lowfreqidx):
        print("Warning: used thermodynamics relations inappropriate for low-frequency modes: {}".format(
            filtered_omega_str[lowfreqidx]))

    sm[('S', 'vib')] = np.sum(rT / np.expm1(rT) - np.log(1 - np.exp(-rT)))
    sm[('Cv', 'vib')] = np.sum(np.exp(rT) * (rT / np.expm1(rT))**2)
    sm[('Cp', 'vib')] = sm[('Cv', 'vib')]
    sm[('ZPE', 'vib')] = np.sum(rT) * T / 2
    sm[('E', 'vib')] = sm[('ZPE', 'vib')] + np.sum(rT * T / np.expm1(rT))
    sm[('H', 'vib')] = sm[('E', 'vib')]

    assert (abs(ZPE_cm_1 - sm[('ZPE', 'vib')] * psi_R * psi_hartree2wavenumbers * 0.001 / psi_hartree2kJmol) < 0.1)

    #real_vibs = np.ma.masked_where(vibinfo['omega'].data.imag > vibinfo['omega'].data.real, vibinfo['omega'].data)

    # compute Gibbs
    for term in ['elec', 'trans', 'rot', 'vib']:
        sm[('G', term)] = sm[('H', term)] - T * sm[('S', term)]

    # convert to atomic units
    for term in ['elec', 'trans', 'rot', 'vib']:
        # terms above are unitless (S, Cv, Cp) or in units of temperature (ZPE, E, H, G) as expressions are divided by R.
        # R [Eh/K], computed as below, slightly diff in 7th sigfig from 3.1668114e-6 (k_B in [Eh/K])
        #    value listed https://en.wikipedia.org/wiki/Boltzmann_constant
        uconv_R_EhK = psi_R / psi_hartree2kJmol
        for piece in ['S', 'Cv', 'Cp']:
            sm[(piece, term)] *= uconv_R_EhK  # [mEh/K] <-- []
        for piece in ['ZPE', 'E', 'H', 'G']:
            sm[(piece, term)] *= uconv_R_EhK * 0.001  # [Eh] <-- [K]

    # sum corrections and totals
    for piece in ['S', 'Cv', 'Cp']:
        for term in ['elec', 'trans', 'rot', 'vib']:
            sm[(piece, 'tot')] += sm[(piece, term)]
    for piece in ['ZPE', 'E', 'H', 'G']:
        for term in ['elec', 'trans', 'rot', 'vib']:
            sm[(piece, 'corr')] += sm[(piece, term)]
        sm[(piece, 'tot')] = E0 + sm[(piece, 'corr')]

    # package results for export
    therminfo = {}
    for entry in sm:
        if entry[0] in ['S', 'Cv', 'Cp']:
            unit = '[mEh/K]'
        elif entry[0] in ['ZPE', 'E', 'H', 'G']:
            unit = '[Eh]'
        therminfo['_'.join(entry)] = VibrationAspect('XXX', unit, sm[entry], '')

    # display
    terms = collections.OrderedDict()
    terms['elec'] = '  Electronic'
    terms['trans'] = '  Translational'
    terms['rot'] = '  Rotational'
    terms['vib'] = '  Vibrational'
    terms['tot'] = 'Total'

    format_S_Cv_Cp = """\n  {:19} {:11.3f} [cal/(mol K)]  {:11.3f} [J/(mol K)]  {:15.8f} [mEh/K]"""
    format_ZPE_E_H_G = """\n  {:19} {:11.3f} [kcal/mol]  {:11.3f} [kJ/mol]  {:15.8f} [Eh]"""
    uconv = np.asarray([psi_hartree2kcalmol, psi_hartree2kJmol, 1.])

    # TODO rot_const, rotor_type
    text = ''
    text += """\n  ==> Thermochemistry Components <=="""

    text += """\n\n  Entropy, S"""
    for term in terms:
        text += format_S_Cv_Cp.format(terms[term] + ' S', *sm[('S', term)] * uconv)
        if term == 'elec':
            text += """ (multiplicity = {})""".format(multiplicity)
        elif term == 'trans':
            text += """ (mol. weight = {:.4f} [u], P = {:.2f} [Pa])""".format(molecular_mass, P)
        elif term == 'rot':
            text += """ (symmetry no. = {})""".format(sigma)

    text += """\n\n  Constant volume heat capacity, Cv"""
    for term in terms:
        text += format_S_Cv_Cp.format(terms[term] + ' Cv', *sm[('Cv', term)] * uconv)

    text += """\n\n  Constant pressure heat capacity, Cp"""
    for term in terms:
        text += format_S_Cv_Cp.format(terms[term] + ' Cp', *sm[('Cp', term)] * uconv)

    del terms['tot']
    terms['corr'] = 'Correction'

    text += """\n\n  ==> Thermochemistry Energy Analysis <=="""

    text += """\n\n  Raw electronic energy, E0"""
    text += """\n  Total E0, Electronic energy at well bottom at 0 [K]               {:15.8f} [Eh]""".format(E0)

    text += """\n\n  Zero-point energy, ZPE_vib = Sum_i nu_i / 2"""
    for term in terms:
        text += format_ZPE_E_H_G.format(terms[term] + ' ZPE', *sm[('ZPE', term)] * uconv)
        if term in ['vib', 'corr']:
            text += """ {:15.3f} [cm^-1]""".format(sm[('ZPE', term)] * psi_hartree2wavenumbers)
    text += """\n  Total ZPE, Electronic energy at 0 [K]                             {:15.8f} [Eh]""".format(
        sm[('ZPE', 'tot')])

    text += """\n\n  Thermal Energy, E (includes ZPE)"""
    for term in terms:
        text += format_ZPE_E_H_G.format(terms[term] + ' E', *sm[('E', term)] * uconv)
    text += """\n  Total E, Electronic energy at {:7.2f} [K]                         {:15.8f} [Eh]""".format(
        T, sm[('E', 'tot')])

    text += """\n\n  Enthalpy, H_trans = E_trans + k_B * T"""
    for term in terms:
        text += format_ZPE_E_H_G.format(terms[term] + ' H', *sm[('H', term)] * uconv)
    text += """\n  Total H, Enthalpy at {:7.2f} [K]                                  {:15.8f} [Eh]""".format(
        T, sm[('H', 'tot')])

    text += """\n\n  Gibbs free energy, G = H - T * S"""
    for term in terms:
        text += format_ZPE_E_H_G.format(terms[term] + ' G', *sm[('G', term)] * uconv)
    text += """\n  Total G, Free enthalpy at {:7.2f} [K]                             {:15.8f} [Eh]\n""".format(
        T, sm[('G', 'tot')])

    return therminfo, text


def filter_nonvib(vibinfo, remove=None):
    """From a dictionary of VibrationAspect, remove normal coordinates.

    Parameters
    ----------
    vibinfo : dict of VibrationAspect
        Results of Hessian analysis.
    remove : list of int, optional
	    0-indexed indices of normal modes to remove from `vibinfo`. If
	    None, non-vibrations (R, T, or TR as labeled in `vibinfo['TRV']`)
	    will be removed.

    Returns
    -------
    dict of VibrationAspect
        Copy of input `vibinfo` with the specified modes removed from all
        dictionary entries.

    Examples
    --------
    # after a harmonic analysis, remove first translations and rotations and then all non-A1 vibs
    >>> allnormco = harmonic_analysis(...)
    >>> allvibs = filter_nonvib(allnormco)
    >>> a1vibs = filter_nonvib(allvibs, remove=[i for i, d in enumerate(allvibs['gamma'].data) if d != 'A1'])

    """
    work = {}
    if remove is None:
        remove = [idx for idx, dat in enumerate(vibinfo['TRV'].data) if dat != 'V']
    for asp, oasp in vibinfo.items():
        if asp in ['q', 'w', 'x']:
            axis = 1
        else:
            axis = 0
        work[asp] = VibrationAspect(oasp.lbl, oasp.unit, np.delete(oasp.data, remove, axis=axis), oasp.comment)

    return work


def filter_omega_to_real(omega):
    """Returns ndarray (float) of `omega` (complex) where imaginary entries are converted to negative reals."""
    freqs = []
    for fr in omega:
        if fr.imag > fr.real:
            freqs.append(-1 * fr.imag)
        else:
            freqs.append(fr.real)
    return np.asarray(freqs)


def _get_TR_space(m, geom, space='TR', tol=None, verbose=1):
    """Form the idealized translation and rotation dof from geometry `geom` and masses `m`.
    Remove any linear dependencies and return an array of shape (3, 3) for atoms, (5, 3 * nat) for linear `geom`,
    or (6, 3 * nat) otherwise. To handle noisy linear geometries, pass `tol` on the order of max deviation.

    m1 = np.asarray([1.])
    m2 = np.asarray([1., 1.])
    m3 = np.asarray([1., 1., 1.])
    m4 = np.asarray([1., 1., 1., 1.])
    g4 = np.asarray([[ 1.,  1., 0.],
          [-1.,  1., 0.],
          [-1., -1., 0.],
          [ 1., -1., 0.]])
    g2 = np.asarray([[ 1.,  1., 0.],
          [-1., -1., 0.]])
    g3 = np.asarray([[3., 3., 3.],
          [4., 4., 4.,],
          [5., 5., 5.]])
    g3noisy = np.asarray([[3., 3.001, 3.],
          [4., 4.001, 4.,],
          [5., 5., 5.01]])
    g33 = np.asarray([[0., 0., 0.],
                     [1., 0., 0.],
                     [-1., 0., 0.]])
    g1 = np.asarray([[0., 0., 0.]])
    g11 = np.asarray([[1., 2., 3.]])
    noise = np.random.normal(0, 1, 9).reshape(3, 3)
    noise = np.divide(noise, np.max(noise))

    assert(_get_TR_space(m4, g4).shape == (6, 12))
    assert(_get_TR_space(m2, g2).shape == (5, 6))
    assert(_get_TR_space(m3, g3).shape == (5, 9))
    assert(_get_TR_space(m3, g33).shape == (5, 9))
    assert(_get_TR_space(m1, g1).shape == (3, 3))
    assert(_get_TR_space(m1, g11).shape == (3, 3))
    assert(_get_TR_space(m3, g3noisy, tol=1.e-2).shape == (5, 9))
    for ns in range(2, 6):
        tol = 10. ** -ns
        gnoisy = g3 + tol * noise
        assert(_get_TR_space(m3, gnoisy, tol=10*tol).shape == (5, 9))

    """
    sqrtmmm = np.repeat(np.sqrt(m), 3)
    xxx = np.repeat(geom[:, 0], 3)
    yyy = np.repeat(geom[:, 1], 3)
    zzz = np.repeat(geom[:, 2], 3)

    z = np.zeros_like(m)
    i = np.ones_like(m)
    ux = np.ravel([i, z, z], order='F')
    uy = np.ravel([z, i, z], order='F')
    uz = np.ravel([z, z, i], order='F')

    # form translation and rotation unit vectors
    T1 = sqrtmmm * ux
    T2 = sqrtmmm * uy
    T3 = sqrtmmm * uz
    R4 = sqrtmmm * (yyy * uz - zzz * uy)
    R5 = sqrtmmm * (zzz * ux - xxx * uz)
    R6 = sqrtmmm * (xxx * uy - yyy * ux)

    TRspace = []
    if 'T' in space:
        TRspace.append([T1, T2, T3])
    if 'R' in space:
        TRspace.append([R4, R5, R6])
    if not TRspace:
        # not sure about this, but it runs
        ZZ = np.zeros_like(T1)
        TRspace.append([ZZ])

    TRspace = np.vstack(TRspace)

    def orth(A, tol=tol):
        u, s, vh = np.linalg.svd(A, full_matrices=False)
        if verbose >= 2:
            print(s)
        M, N = A.shape
        eps = np.finfo(float).eps
        if tol is None:
            tol = max(M, N) * np.amax(s) * eps
        num = np.sum(s > tol, dtype=int)
        Q = u[:, :num]
        return Q

    TRindep = orth(TRspace.T)
    TRindep = TRindep.T

    if verbose >= 2:
        print(TRindep.shape, '<--', TRspace.shape)
        print(np.linalg.norm(TRindep, axis=1))
        print('-' * 80)

    return TRindep
