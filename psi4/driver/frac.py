#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2023 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

__all__ = [
    "frac_nuke",
    "frac_traverse",
    "ip_fitting",
]

from typing import Callable, Dict, Union

from psi4 import core
from psi4.driver import p4util
from psi4.driver import driver
from psi4.driver.p4util.exceptions import *


def frac_traverse(name: Union[str, Callable], **kwargs) -> Dict[float, float]:
    """Scan electron occupancy from +1 electron to -1.

    Parameters
    ----------
    name
        DFT functional string name or function defining functional
        whose omega is to be optimized.
    molecule : :ref:`molecule <op_py_molecule>`, optional
        Target molecule (neutral) for which omega is to be tuned, if not last defined.
    cation_mult : Optional[int]
        Multiplicity of cation, if not neutral multiplicity + 1.
    anion_mult : Optional[int]
        Multiplicity of anion, if not neutral multiplicity + 1.
    frac_start : Optional[int]
        Iteration at which to start frac procedure when not reading previous
        guess. Defaults to 25.
    HOMO_occs : Optional[List]
        Occupations to step through for cation, by default `[1 - 0.1 * x for x in range(11)]`.
    LUMO_occs : Optional[List]
        Occupations to step through for anion, by default `[1 - 0.1 * x for x in range(11)]`.
    HOMO : Optional[int]
        Index of HOMO.
    LUMO : Optional[int]
        Index of LUMO.
    frac_diis : Optional[bool]
        Do use DIIS for non-1.0-occupied points?
    neutral_guess : Optional[bool]
        Do use neutral orbitals as guess for the anion?
    hf_guess: Optional[bool]
        Do use UHF guess before UKS?
    continuous_guess : Optional[bool]
        Do carry along guess rather than reguessing at each occupation?
    filename : Optional[str]
        Result filename, if not name of molecule.

    Returns
    -------
    ~typing.Dict[float, float]
        Dictionary associating SCF energies with occupations.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'GUESS'],
        ['SCF', 'DF_INTS_IO'],
        ['SCF', 'REFERENCE'],
        ["SCF", "FRAC_START"],
        ["SCF", "FRAC_RENORMALIZE"],
        #["SCF", "FRAC_LOAD"],
        ["SCF", "FRAC_OCC"],
        ["SCF", "FRAC_VAL"],
        ["SCF", "FRAC_DIIS"])
    kwargs = p4util.kwargs_lower(kwargs)

    # Make sure the molecule the user provided is the active one, and neutral
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()

    if molecule.molecular_charge() != 0:
        raise ValidationError("""frac_traverse requires neutral molecule to start.""")
    if molecule.schoenflies_symbol() != 'c1':
        core.print_out("""  Requested procedure `frac_traverse` does not make use of molecular symmetry: """
                       """further calculations in C1 point group.\n""")
    molecule = molecule.clone()
    molecule.reset_point_group('c1')
    molecule.update_geometry()

    charge0 = molecule.molecular_charge()
    mult0 = molecule.multiplicity()

    chargep = charge0 + 1
    chargem = charge0 - 1

    multp = kwargs.get('cation_mult', mult0 + 1)
    multm = kwargs.get('anion_mult', mult0 + 1)

    # By default, we start the frac procedure on the 25th iteration
    # when not reading a previous guess
    frac_start = kwargs.get('frac_start', 25)

    # By default, we occupy by tenths of electrons
    HOMO_occs = kwargs.get('HOMO_occs', [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])
    LUMO_occs = kwargs.get('LUMO_occs', [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])

    # By default, HOMO and LUMO are both in alpha
    Z = 0
    for A in range(molecule.natom()):
        Z += molecule.Z(A)
    Z -= charge0
    HOMO = kwargs.get('HOMO', (Z / 2 + 1 if (Z % 2) else Z / 2))
    LUMO = kwargs.get('LUMO', HOMO + 1)

    # By default, DIIS in FRAC (1.0 occupation is always DIIS'd)
    frac_diis = kwargs.get('frac_diis', True)

    # By default, use the neutral orbitals as a guess for the anion
    neutral_guess = kwargs.get('neutral_guess', True)

    # By default, burn-in with UHF first, if UKS
    hf_guess = False
    if core.get_local_option('SCF', 'REFERENCE') == 'UKS':
        hf_guess = kwargs.get('hf_guess', True)

    # By default, re-guess at each N
    continuous_guess = kwargs.get('continuous_guess', False)

    # By default, drop the files to the molecule's name
    root = kwargs.get('filename', molecule.name())
    traverse_filename = root + '.traverse.dat'
    # => Traverse <= #
    occs = []
    energies = []
    potentials = []
    convs = []

    # => Run the neutral for its orbitals, if requested <= #

    core.set_local_option("SCF", "DF_INTS_IO", "SAVE")

    old_guess = core.get_local_option("SCF", "GUESS")
    if (neutral_guess):
        if (hf_guess):
            core.set_local_option("SCF", "REFERENCE", "UHF")
        driver.energy('scf', dft_functional=name, molecule=molecule, **kwargs)
        core.set_local_option("SCF", "GUESS", "READ")
        core.set_local_option("SCF", "DF_INTS_IO", "LOAD")

    # => Run the anion first <= #

    molecule.set_molecular_charge(chargem)
    molecule.set_multiplicity(multm)

    # => Burn the anion in with hf, if requested <= #
    if hf_guess:
        core.set_local_option("SCF", "REFERENCE","UHF")
        driver.energy('scf', dft_functional=name, molecule=molecule, **kwargs)
        core.set_local_option("SCF", "REFERENCE", "UKS")
        core.set_local_option("SCF", "GUESS", "READ")
        core.set_local_option("SCF", "DF_INTS_IO", "SAVE")

    core.set_local_option("SCF", "FRAC_START", frac_start)
    core.set_local_option("SCF", "FRAC_RENORMALIZE", True)
    # NYI core.set_local_option("SCF", "FRAC_LOAD", False)

    for occ in LUMO_occs:

        core.set_local_option("SCF", "FRAC_OCC", [LUMO])
        core.set_local_option("SCF", "FRAC_VAL", [occ])

        E, wfn = driver.energy('scf', dft_functional=name, return_wfn=True, molecule=molecule, **kwargs)
        C = 1
        if E == 0.0:
            E = core.variable('SCF ITERATION ENERGY')
            C = 0

        if LUMO > 0:
            eps = wfn.epsilon_a()
            potentials.append(eps.get(int(LUMO) - 1))
        else:
            eps = wfn.epsilon_b()
            potentials.append(eps.get(-int(LUMO) - 1))

        occs.append(occ)
        energies.append(E)
        convs.append(C)

        core.set_local_option("SCF", "FRAC_START", 2)
        #core.set_local_option("SCF", "FRAC_LOAD", True)
        core.set_local_option("SCF", "GUESS", "READ")
        core.set_local_option("SCF", "FRAC_DIIS", frac_diis)
        core.set_local_option("SCF", "DF_INTS_IO", "LOAD")


    # => Run the neutral next <= #

    molecule.set_molecular_charge(charge0)
    molecule.set_multiplicity(mult0)

    # Burn the neutral in with hf, if requested <= #

    if not continuous_guess:
        core.set_local_option("SCF", "GUESS", old_guess)
        if hf_guess:
            core.set_local_option("SCF", "FRAC_START", 0)
            core.set_local_option("SCF", "REFERENCE", "UHF")
            driver.energy('scf', dft_functional=name, molecule=molecule, **kwargs)
            core.set_local_option("SCF", "REFERENCE", "UKS")
            core.set_local_option("SCF", "GUESS", "READ")
        # NYI core.set_local_option("SCF", "FRAC_LOAD", False)

    core.set_local_option("SCF", "FRAC_START", frac_start)
    core.set_local_option("SCF", "FRAC_RENORMALIZE", True)

    for occ in HOMO_occs:

        core.set_local_option("SCF", "FRAC_OCC", [HOMO])
        core.set_local_option("SCF", "FRAC_VAL", [occ])

        E, wfn = driver.energy('scf', dft_functional=name, return_wfn=True, molecule=molecule, **kwargs)
        C = 1
        if E == 0.0:
            E = core.variable('SCF ITERATION ENERGY')
            C = 0

        if LUMO > 0:
            eps = wfn.epsilon_a()
            potentials.append(eps.get(int(HOMO) - 1))
        else:
            eps = wfn.epsilon_b()
            potentials.append(eps.get(-int(HOMO) - 1))

        occs.append(occ - 1.0)
        energies.append(E)
        convs.append(C)

        core.set_local_option("SCF", "FRAC_START", 2)
        # NYI core.set_local_option("SCF", "FRAC_LOAD", True)
        core.set_local_option("SCF", "GUESS", "READ")
        core.set_local_option("SCF", "FRAC_DIIS", frac_diis)
        core.set_local_option("SCF", "DF_INTS_IO", "LOAD")

    # => Print the results out <= #
    E = {}
    core.print_out("""\n    ==> Fractional Occupation Traverse Results <==\n\n""")
    core.print_out("""    %-11s %-24s %-24s %11s\n""" % ('N', 'Energy', 'HOMO Energy', 'Converged'))
    for k in range(len(occs)):
        core.print_out("""    %11.3E %24.16E %24.16E %11d\n""" % (occs[k], energies[k], potentials[k], convs[k]))
        E[occs[k]] = energies[k]

    core.print_out("""
    You trying to be a hero Watkins?
    Just trying to kill some bugs sir!
            -Starship Troopers""")

    # Drop the files out
    with open(traverse_filename, 'w') as fh:
        fh.write("""    %-11s %-24s %-24s %11s\n""" % ('N', 'Energy', 'HOMO Energy', 'Converged'))
        for k in range(len(occs)):
            fh.write("""    %11.3E %24.16E %24.16E %11d\n""" % (occs[k], energies[k], potentials[k], convs[k]))

    optstash.restore()
    return E


def frac_nuke(name: Union[str, Callable], **kwargs) -> Dict[float, float]:
    """Pull all the electrons out, one at a time"""
    optstash = p4util.OptionsState(
        ['SCF', 'GUESS'],
        ['SCF', 'DF_INTS_IO'],
        ["SCF", "FRAC_START"],
        ["SCF", "FRAC_RENORMALIZE"],
        # NYI ["SCF", "FRAC_LOAD"],
        ["SCF", "FRAC_OCC"],
        ["SCF", "FRAC_VAL"],
        ["SCF", "FRAC_DIIS"])

    kwargs = p4util.kwargs_lower(kwargs)

    # Make sure the molecule the user provided is the active one, and neutral
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()

    if molecule.molecular_charge() != 0:
        raise ValidationError("""frac_nuke requires neutral molecule to start.""")
    if molecule.schoenflies_symbol() != 'c1':
        core.print_out("""  Requested procedure `frac_nuke` does not make use of molecular symmetry: """
                       """further calculations in C1 point group.\n""")
    molecule = molecule.clone()
    molecule.reset_point_group('c1')
    molecule.update_geometry()

    charge0 = molecule.molecular_charge()
    mult0 = molecule.multiplicity()

    # By default, we start the frac procedure on the 25th iteration
    # when not reading a previous guess
    frac_start = kwargs.get('frac_start', 25)

    # By default, we occupy by tenths of electrons
    foccs = kwargs.get('foccs', [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])

    # By default, HOMO and LUMO are both in alpha
    N = 0
    for A in range(molecule.natom()):
        N += molecule.Z(A)
    N -= charge0
    N = int(N)
    Nb = int((N - mult0 + 1) / 2)
    Na = int(N - Nb)

    charge = charge0
    mult = mult0

    # By default, nuke all the electrons
    Nmin = 0
    if 'nmax' in kwargs:
        Nmin = N - int(kwargs['nmax'])

    # By default, DIIS in FRAC (1.0 occupation is always DIIS'd)
    frac_diis = kwargs.get('frac_diis', True)

    # By default, drop the files to the molecule's name
    root = kwargs.get('filename', molecule.name())
    traverse_filename = root + '.traverse.dat'
    stats_filename = root + '.stats.dat'

    # => Traverse <= #
    core.set_local_option("SCF", "DF_INTS_IO", "SAVE")

    Ns = []
    energies = []
    potentials = []
    convs = []
    stats = []

    # Run one SCF to burn things in
    E, wfn = driver.energy('scf', dft_functional=name, return_wfn=True, molecule=molecule, **kwargs)

    # Determine HOMO
    eps_a = wfn.epsilon_a()
    eps_b = wfn.epsilon_b()
    eps_a.print_out()
    if Na == Nb:
        HOMO = -Nb
    elif Nb == 0:
        HOMO = Na
    else:
        E_a = eps_a.get(int(Na - 1))
        E_b = eps_b.get(int(Nb - 1))
        if E_a >= E_b:
            HOMO = Na
        else:
            HOMO = -Nb

    stats.append("""    %6d %6d %6d %6d %6d %6d\n""" % (N, Na, Nb, charge, mult, HOMO))

    if HOMO > 0:
        Na -= 1
    else:
        Nb -= 1
    charge += 1
    mult = Na - Nb + 1

    core.set_local_option("SCF", "DF_INTS_IO", "LOAD")
    core.set_local_option("SCF", "FRAC_START", frac_start)
    core.set_local_option("SCF", "FRAC_RENORMALIZE", True)

    # Nuke 'em Rico!
    for Nintegral in range(N, Nmin, -1):

        # Nuke the current HOMO
        for occ in foccs:

            core.set_local_option("SCF", "FRAC_OCC", [HOMO])
            core.set_local_option("SCF", "FRAC_VAL", [occ])

            E, wfn = driver.energy('scf', dft_functional=name, return_wfn=True, molecule=molecule, **kwargs)
            C = 1
            if E == 0.0:
                E = core.variable('SCF ITERATION ENERGY')
                C = 0

            if HOMO > 0:
                eps = wfn.epsilon_a()
                potentials.append(eps.np[HOMO - 1])
            else:
                eps = wfn.epsilon_b()
                potentials.append(eps.np[-HOMO - 1])

            Ns.append(Nintegral + occ - 1.0)
            energies.append(E)
            convs.append(C)

            core.set_local_option("SCF", "FRAC_START", 2)
            # NYI core.set_local_option("SCF", "FRAC_LOAD", True)
            core.set_local_option("SCF", "FRAC_DIIS", frac_diis)
            core.set_local_option("SCF", "GUESS", "READ")

        # Set the next charge/mult
        molecule.set_molecular_charge(charge)
        molecule.set_multiplicity(mult)

        # Determine HOMO
        print('DGAS: What ref should this point to?')
        eps_a = wfn.epsilon_a()
        eps_b = wfn.epsilon_b()
        if Na == Nb:
            HOMO = -Nb
        elif Nb == 0:
            HOMO = Na
        else:
            E_a = eps_a.np[int(Na - 1)]
            E_b = eps_b.np[int(Nb - 1)]
            if E_a >= E_b:
                HOMO = Na
            else:
                HOMO = -Nb

        stats.append("""    %6d %6d %6d %6d %6d %6d\n""" % (Nintegral-1, Na, Nb, charge, mult, HOMO))

        if HOMO > 0:
            Na -= 1
        else:
            Nb -= 1
        charge += 1
        mult = Na - Nb + 1

    core.set_local_option("SCF", "DF_INTS_IO", "NONE")

    # => Print the results out <= #
    E = {}
    core.print_out("""\n    ==> Fractional Occupation Nuke Results <==\n\n""")
    core.print_out("""    %-11s %-24s %-24s %11s\n""" % ('N', 'Energy', 'HOMO Energy', 'Converged'))
    for k in range(len(Ns)):
        core.print_out("""    %11.3E %24.16E %24.16E %11d\n""" % (Ns[k], energies[k], potentials[k], convs[k]))
        E[Ns[k]] = energies[k]

    core.print_out('\n')
    core.print_out("""    %6s %6s %6s %6s %6s %6s\n""" % ('N', 'Na', 'Nb', 'Charge', 'Mult', 'HOMO'))
    for line in stats:
        core.print_out(line)

    core.print_out('\n    "You shoot a nuke down a bug hole, you got a lot of dead bugs"\n')
    core.print_out('            -Starship Troopers\n')

    # Drop the files out
    with open(traverse_filename, 'w') as fh:
        fh.write("""    %-11s %-24s %-24s %11s\n""" % ('N', 'Energy', 'HOMO Energy', 'Converged'))
        for k in range(len(Ns)):
            fh.write("""    %11.3E %24.16E %24.16E %11d\n""" % (Ns[k], energies[k], potentials[k], convs[k]))

    with open(stats_filename, 'w') as fh:
        fh.write("""    %6s %6s %6s %6s %6s %6s\n""" % ('N', 'Na', 'Nb', 'Charge', 'Mult', 'HOMO'))
        for line in stats:
            fh.write(line)

    optstash.restore()
    return E


def ip_fitting(name: Union[str, Callable], omega_l: float = 0.05, omega_r: float = 2.5, omega_convergence: float = 1.0e-3, maxiter: int = 20, **kwargs) -> float:
    """Optimize DFT omega parameter for molecular system.

    Parameters
    ----------
    name
        DFT functional string name or function defining functional
        whose omega is to be optimized.
    omega_l
        Minimum omega to be considered during fitting.
    omega_r
        Maximum omega to be considered during fitting.
    molecule : :ref:`molecule <op_py_molecule>`, optional
        Target molecule (neutral) for which omega is to be tuned, if not last defined.
    omega_convergence
        Threshold below which to consider omega converged. (formerly omega_tolerance)
    maxiter
        Maximum number of iterations towards omega convergence.

    Returns
    -------
    float
        Optimal omega parameter.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'REFERENCE'],
        ['SCF', 'GUESS'],
        ['SCF', 'DF_INTS_IO'],
        ['SCF', 'DFT_OMEGA'],
        ['DOCC'],
        ['SOCC'])

    kwargs = p4util.kwargs_lower(kwargs)

    # By default, do not read previous 180 orbitals file
    read = False
    read180 = ''
    if 'read' in kwargs:
        read = True
        read180 = kwargs['read']

    if core.get_option('SCF', 'REFERENCE') != 'UKS':
        core.print_out("""  Requested procedure `ip_fitting` runs further calculations with UKS reference.\n""")
        core.set_local_option('SCF', 'REFERENCE', 'UKS')

    # Make sure the molecule the user provided is the active one, and neutral
    molecule = kwargs.pop('molecule', core.get_active_molecule())
    molecule.update_geometry()

    if molecule.molecular_charge() != 0:
        raise ValidationError("""IP Fitting requires neutral molecule to start.""")
    if molecule.schoenflies_symbol() != 'c1':
        core.print_out("""  Requested procedure `ip_fitting` does not make use of molecular symmetry: """
                       """further calculations in C1 point group.\n""")
    molecule = molecule.clone()
    molecule.reset_point_group('c1')
    molecule.update_geometry()

    charge0 = molecule.molecular_charge()
    mult0 = molecule.multiplicity()

    # How many electrons are there?
    N = 0
    for A in range(molecule.natom()):
        N += molecule.Z(A)
    N -= charge0
    N = int(N)
    Nb = int((N - mult0 + 1) / 2)
    Na = int(N - Nb)

    # Work in the ot namespace for this procedure
    core.IO.set_default_namespace("ot")

    # Burn in to determine orbital eigenvalues
    if read:
        core.set_local_option("SCF", "GUESS", "READ")
        copy_file_to_scratch(read180, 'psi', 'ot', 180)
    core.set_local_option("SCF", "DF_INTS_IO", "SAVE")
    E, wfn = driver.energy('scf', dft_functional=name, return_wfn=True, molecule=molecule,
                           banner='IP Fitting SCF: Burn-in', **kwargs)
    core.set_local_option("SCF", "DF_INTS_IO", "LOAD")

    if not wfn.functional().is_x_lrc():
        raise ValidationError("""Not sensible to optimize omega for non-long-range-correction functional.""")

    # Determine HOMO, to determine mult1
    eps_a = wfn.epsilon_a()
    eps_b = wfn.epsilon_b()
    if Na == Nb:
        HOMO = -Nb
    elif Nb == 0:
        HOMO = Na
    else:
        E_a = eps_a.np[int(Na - 1)]
        E_b = eps_b.np[int(Nb - 1)]
        if E_a >= E_b:
            HOMO = Na
        else:
            HOMO = -Nb

    Na1 = Na
    Nb1 = Nb
    if HOMO > 0:
        Na1 -= 1
    else:
        Nb1 -= 1

    charge1 = charge0 + 1
    mult1 = Na1 - Nb1 + 1

    omegas = []
    E0s = []
    E1s = []
    kIPs = []
    IPs = []
    types = []

    # Right endpoint
    core.set_local_option('SCF', 'DFT_OMEGA', omega_r)

    # Neutral
    if read:
        core.set_local_option("SCF", "GUESS", "READ")
        p4util.copy_file_to_scratch(read180, 'psi', 'ot', 180)

    molecule.set_molecular_charge(charge0)
    molecule.set_multiplicity(mult0)
    E0r, wfn = driver.energy('scf', dft_functional=name, return_wfn=True, molecule=molecule,
                             banner='IP Fitting SCF: Neutral, Right Endpoint', **kwargs)
    eps_a = wfn.epsilon_a()
    eps_b = wfn.epsilon_b()
    if Nb == 0:
        E_HOMO = eps_a.np[int(Na - 1)]
    else:
        E_a = eps_a.np[int(Na - 1)]
        E_b = eps_b.np[int(Nb - 1)]
        E_HOMO = max(E_a, E_b)
    E_HOMOr = E_HOMO
    core.IO.change_file_namespace(180, "ot", "neutral")

    # Cation
    if read:
        core.set_local_option("SCF", "GUESS", "READ")
        p4util.copy_file_to_scratch(read180, 'psi', 'ot', 180)

    molecule.set_molecular_charge(charge1)
    molecule.set_multiplicity(mult1)
    E1r = driver.energy('scf', dft_functional=name, molecule=molecule,
                               banner='IP Fitting SCF: Cation, Right Endpoint', **kwargs)
    core.IO.change_file_namespace(180, "ot", "cation")

    IPr = E1r - E0r
    kIPr = -E_HOMOr
    delta_r = IPr - kIPr

    if IPr > kIPr:
        raise ValidationError("""\n***IP Fitting Error: Right Omega limit should have kIP > IP: {} !> {}""".format(kIPr, IPr))

    omegas.append(omega_r)
    types.append('Right Limit')
    E0s.append(E0r)
    E1s.append(E1r)
    IPs.append(IPr)
    kIPs.append(kIPr)

    # Use previous orbitals from here out
    core.set_local_option("SCF", "GUESS", "READ")

    # Left endpoint
    core.set_local_option('SCF', 'DFT_OMEGA', omega_l)

    # Neutral
    core.IO.change_file_namespace(180, "neutral", "ot")
    molecule.set_molecular_charge(charge0)
    molecule.set_multiplicity(mult0)
    core.set_global_option("DOCC", [Nb])
    core.set_global_option("SOCC", [Na - Nb])
    E0l, wfn = driver.energy('scf', dft_functional=name, return_wfn=True, molecule=molecule,
                             banner='IP Fitting SCF: Neutral, Left Endpoint', **kwargs)
    eps_a = wfn.epsilon_a()
    eps_b = wfn.epsilon_b()
    if Nb == 0:
        E_HOMO = eps_a.np[int(Na - 1)]
    else:
        E_a = eps_a.np[int(Na - 1)]
        E_b = eps_b.np[int(Nb - 1)]
        E_HOMO = max(E_a, E_b)
    E_HOMOl = E_HOMO
    core.IO.change_file_namespace(180, "ot", "neutral")

    # Cation
    core.IO.change_file_namespace(180, "cation", "ot")
    molecule.set_molecular_charge(charge1)
    molecule.set_multiplicity(mult1)
    core.set_global_option("DOCC", [Nb1])
    core.set_global_option("SOCC", [Na1 - Nb1])
    E1l = driver.energy('scf', dft_functional=name, molecule=molecule,
                        banner='IP Fitting SCF: Cation, Left Endpoint', **kwargs)
    core.IO.change_file_namespace(180, "ot", "cation")

    IPl = E1l - E0l
    kIPl = -E_HOMOl
    delta_l = IPl - kIPl

    if IPl < kIPl:
        raise ValidationError("""\n***IP Fitting Error: Left Omega limit should have kIP < IP: {} !< {}""".format(kIPl, IPl))

    omegas.append(omega_l)
    types.append('Left Limit')
    E0s.append(E0l)
    E1s.append(E1l)
    IPs.append(IPl)
    kIPs.append(kIPl)

    converged = False
    repeat_l = 0
    repeat_r = 0
    for step in range(maxiter):

        # Regula Falsi (modified)
        if repeat_l > 1:
            delta_l /= 2.0
        if repeat_r > 1:
            delta_r /= 2.0
        omega = - (omega_r - omega_l) / (delta_r - delta_l) * delta_l + omega_l
        core.set_local_option('SCF', 'DFT_OMEGA', omega)

        # Neutral
        core.IO.change_file_namespace(180, "neutral", "ot")
        molecule.set_molecular_charge(charge0)
        molecule.set_multiplicity(mult0)
        core.set_global_option("DOCC", [Nb])
        core.set_global_option("SOCC", [Na - Nb])
        E0, wfn = driver.energy('scf', dft_functional=name, return_wfn=True, molecule=molecule,
                                banner='IP Fitting SCF: Neutral, Omega = {:11.3E}'.format(omega), **kwargs)
        eps_a = wfn.epsilon_a()
        eps_b = wfn.epsilon_b()
        if Nb == 0:
            E_HOMO = eps_a.np[int(Na - 1)]
        else:
            E_a = eps_a.np[int(Na - 1)]
            E_b = eps_b.np[int(Nb - 1)]
            E_HOMO = max(E_a, E_b)
        core.IO.change_file_namespace(180, "ot", "neutral")

        # Cation
        core.IO.change_file_namespace(180, "cation", "ot")
        molecule.set_molecular_charge(charge1)
        molecule.set_multiplicity(mult1)
        core.set_global_option("DOCC", [Nb1])
        core.set_global_option("SOCC", [Na1 - Nb1])
        E1 = driver.energy('scf', dft_functional=name, molecule=molecule,
                           banner='IP Fitting SCF: Cation, Omega = {:11.3E}'.format(omega), **kwargs)
        core.IO.change_file_namespace(180, "ot", "cation")

        IP = E1 - E0
        kIP = -E_HOMO
        delta = IP - kIP

        if kIP > IP:
            omega_r = omega
            E0r = E0
            E1r = E1
            IPr = IP
            kIPr = kIP
            delta_r = delta
            repeat_r = 0
            repeat_l += 1
        else:
            omega_l = omega
            E0l = E0
            E1l = E1
            IPl = IP
            kIPl = kIP
            delta_l = delta
            repeat_l = 0
            repeat_r += 1

        omegas.append(omega)
        types.append('Regula-Falsi')
        E0s.append(E0)
        E1s.append(E1)
        IPs.append(IP)
        kIPs.append(kIP)

        # Termination
        if abs(omega_l - omega_r) < omega_convergence:
            converged = True
            break

    core.IO.set_default_namespace("")
    core.print_out("""\n    ==> IP Fitting Results <==\n\n""")

    core.print_out("""     => Occupation Determination <= \n\n""")
    core.print_out("""              %6s %6s %6s %6s %6s %6s\n""" % ('N', 'Na', 'Nb', 'Charge', 'Mult', 'HOMO'))
    core.print_out("""     Neutral: %6d %6d %6d %6d %6d %6d\n""" % (N, Na, Nb, charge0, mult0, HOMO))
    core.print_out("""     Cation:  %6d %6d %6d %6d %6d\n\n""" % (N - 1, Na1, Nb1, charge1, mult1))

    core.print_out("""     => Regula Falsi Iterations <=\n\n""")
    core.print_out("""    %3s %11s %14s %14s %14s %s\n""" % ('N','Omega','IP','kIP','Delta','Type'))
    for k in range(len(omegas)):
        core.print_out("""    %3d %11.3E %14.6E %14.6E %14.6E %s\n""" %
                       (k + 1, omegas[k], IPs[k], kIPs[k], IPs[k] - kIPs[k], types[k]))

    optstash.restore()
    if converged:
        core.print_out("""\n    IP Fitting Converged\n""")
        core.print_out("""    Final omega = %14.6E\n""" % ((omega_l + omega_r) / 2))
        core.print_out("""\n    "M,I. does the dying. Fleet just does the flying."\n""")
        core.print_out("""            -Starship Troopers\n""")

    else:
        raise ConvergenceError("""IP Fitting """, step + 1)

    return ((omega_l + omega_r) / 2)
