#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

from __future__ import absolute_import

import core
import os
import math
import p4util
from molutil import *
from driver import *
from p4util.exceptions import *

# Scan from +1 electron to -1 electron
def frac_traverse(molecule, **kwargs):
    kwargs = p4util.kwargs_lower(kwargs)

    # The molecule is required, and should be the neutral species
    molecule.update_geometry()
    charge0 = molecule.molecular_charge()
    mult0 = molecule.multiplicity()

    chargep = charge0 + 1
    chargem = charge0 - 1

    # By default, the multiplicity of the cation/anion are mult0 + 1
    # These are overridden with the cation_mult and anion_mult kwargs
    multp = kwargs.get('cation_mult', mult0 + 1)
    multm = kwargs.get('anion_mult', mult0 + 1)

    # By default, we start the frac procedure on the 25th iteration
    # when not reading a previous guess
    frac_start = kwargs.get('frac_start', 25)

    # By default, we occupy by tenths of electrons
    HOMO_occs = kwargs.get('HOMO_occs', [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])
    LUMO_occs = kwargs.get('LUMO_occs', [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])

    # By default, HOMO and LUMO are both in alpha
    Z = 0;
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
    if core.get_global_option('REFERENCE') == 'UKS':
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

    old_df_ints_io = core.get_global_option("DF_INTS_IO")
    core.set_global_option("DF_INTS_IO", "SAVE")

    old_guess = core.get_global_option("GUESS")
    if (neutral_guess):
        if (hf_guess):
            core.set_global_option("REFERENCE","UHF")
        energy('scf')
        core.set_global_option("GUESS", "READ")
        core.set_global_option("DF_INTS_IO", "LOAD")

    # => Run the anion first <= #

    molecule.set_molecular_charge(chargem)
    molecule.set_multiplicity(multm)

    # => Burn the anion in with hf, if requested <= #
    if hf_guess:
        core.set_global_option("REFERENCE","UHF")
        energy('scf', molecule=molecule, **kwargs)
        core.set_global_option("REFERENCE","UKS")
        core.set_global_option("GUESS", "READ")
        core.set_global_option("DF_INTS_IO", "SAVE")

    core.set_global_option("FRAC_START", frac_start)
    core.set_global_option("FRAC_RENORMALIZE", True)
    core.set_global_option("FRAC_LOAD", False)

    for occ in LUMO_occs:

        core.set_global_option("FRAC_OCC", [LUMO])
        core.set_global_option("FRAC_VAL", [occ])

        E, wfn = energy('scf', return_wfn=True, molecule=molecule, **kwargs)
        C = 1
        if E == 0.0:
            E = core.get_variable('SCF ITERATION ENERGY')
            C = 0

        if LUMO > 0:
            eps = wfn.epsilon_a()
            potentials.append(eps[int(LUMO) - 1])
        else:
            eps = wfn.epsilon_b()
            potentials.append(eps[-int(LUMO) - 1])

        occs.append(occ)
        energies.append(E)
        convs.append(C)

        core.set_global_option("FRAC_START", 2)
        core.set_global_option("FRAC_LOAD", True)
        core.set_global_option("GUESS", "READ")
        core.set_global_option("FRAC_DIIS", frac_diis)
        core.set_global_option("DF_INTS_IO", "LOAD")


    # => Run the neutral next <= #

    molecule.set_molecular_charge(charge0)
    molecule.set_multiplicity(mult0)

    # Burn the neutral in with hf, if requested <= #

    if not continuous_guess:
        core.set_global_option("GUESS", old_guess)
        if hf_guess:
            core.set_global_option("FRAC_START", 0)
            core.set_global_option("REFERENCE", "UHF")
            energy('scf', molecule=molecule, **kwargs)
            core.set_global_option("REFERENCE", "UKS")
            core.set_global_option("GUESS", "READ")
        core.set_global_option("FRAC_LOAD", False)

    core.set_global_option("FRAC_START", frac_start)
    core.set_global_option("FRAC_RENORMALIZE", True)

    for occ in HOMO_occs:

        core.set_global_option("FRAC_OCC", [HOMO])
        core.set_global_option("FRAC_VAL", [occ])

        E, wfn = energy('scf', return_wfn=True, molecule=molecule, **kwargs)
        C = 1
        if E == 0.0:
            E = core.get_variable('SCF ITERATION ENERGY')
            C = 0

        if LUMO > 0:
            eps = wfn.epsilon_a()
            potentials.append(eps[int(HOMO) - 1])
        else:
            eps = wfn.epsilon_b()
            potentials.append(eps[-int(HOMO) - 1])

        occs.append(occ - 1.0)
        energies.append(E)
        convs.append(C)

        core.set_global_option("FRAC_START", 2)
        core.set_global_option("FRAC_LOAD", True)
        core.set_global_option("GUESS", "READ")
        core.set_global_option("FRAC_DIIS", frac_diis)
        core.set_global_option("DF_INTS_IO", "LOAD")

    core.set_global_option("DF_INTS_IO", old_df_ints_io)

    # => Print the results out <= #
    E = {}
    core.print_out("""\n    ==> Fractional Occupation Traverse Results <==\n\n""")
    core.print_out("""\t%-11s %-24s %-24s %11s\n""" % ('N', 'Energy', 'HOMO Energy', 'Converged'))
    for k in range(len(occs)):
        core.print_out("""\t%11.3E %24.16E %24.16E %11d\n""" % (occs[k], energies[k], potentials[k], convs[k]))
        E[occs[k]] = energies[k]

    core.print_out('\n\t"You trying to be a hero Watkins?"\n')
    core.print_out('\t"Just trying to kill some bugs sir!"\n')
    core.print_out('\t\t\t-Starship Troopers\n')

    # Drop the files out
    fh = open(traverse_filename, 'w')
    fh.write("""\t%-11s %-24s %-24s %11s\n""" % ('N', 'Energy', 'HOMO Energy', 'Converged'))
    for k in range(len(occs)):
        fh.write("""\t%11.3E %24.16E %24.16E %11d\n""" % (occs[k], energies[k], potentials[k], convs[k]))
    fh.close()

    # Properly, should clone molecule but since not returned and easy to unblemish,
    molecule.set_molecular_charge(charge0)
    molecule.set_multiplicity(mult0)

    return E

# Pull all the electrons out, one at a time
def frac_nuke(molecule, **kwargs):
    kwargs = p4util.kwargs_lower(kwargs)

    # The molecule is required, and should be the neutral species
    molecule.update_geometry()
    charge0 = molecule.molecular_charge()
    mult0 = molecule.multiplicity()

    # By default, we start the frac procedure on the 25th iteration
    # when not reading a previous guess
    frac_start = kwargs.get('frac_start', 25)

    # By default, we occupy by tenths of electrons
    foccs = kwargs.get('foccs', [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])

    # By default, HOMO and LUMO are both in alpha
    N = 0;
    for A in range(molecule.natom()):
        N += molecule.Z(A)
    N -= charge0
    N = int(N)
    Nb = int((N - mult0 + 1) / 2)
    Na = int(N - Nb)

    charge = charge0
    mult = mult0

    # By default, nuke all the electrons
    Nmin = 0;
    if ('nmax' in kwargs):
        Nmin = N - int(kwargs['nmax'])

    # By default, DIIS in FRAC (1.0 occupation is always DIIS'd)
    frac_diis = kwargs.get('frac_diis', True)

    # By default, drop the files to the molecule's name
    root = kwargs.get('filename', molecule.name())
    traverse_filename = root + '.traverse.dat'
    stats_filename = root + '.stats.dat'

    # => Traverse <= #
    core.set_global_option("DF_INTS_IO", "SAVE")

    Ns = []
    energies = []
    potentials = []
    convs = []
    stats = []

    # Run one SCF to burn things in
    E, wfn= energy('scf', return_wfn=True, molecule=molecule, **kwargs)

    # Determine HOMO
    eps_a = wfn.epsilon_a()
    eps_b = wfn.epsilon_b()
    if Na == Nb:
        HOMO = -Nb
    elif Nb == 0:
        HOMO = Na
    else:
        E_a = eps_a[int(Na - 1)]
        E_b = eps_b[int(Nb - 1)]
        if E_a >= E_b:
            HOMO = Na
        else:
            HOMO = -Nb

    stats.append("""\t%6d %6d %6d %6d %6d %6d\n""" % (N, Na, Nb, charge, mult, HOMO))

    if HOMO > 0:
        Na = Na - 1
    else:
        Nb = Nb - 1
    charge = charge + 1
    mult = Na - Nb + 1

    core.set_global_option("DF_INTS_IO", "LOAD")
    core.set_global_option("FRAC_START", frac_start)
    core.set_global_option("FRAC_RENORMALIZE", True)

    # Nuke 'em Rico!
    for Nintegral in range(N, Nmin, -1):

        # Nuke the current HOMO
        for occ in foccs:

            core.set_global_option("FRAC_OCC", [HOMO])
            core.set_global_option("FRAC_VAL", [occ])

            E, wfn = energy('scf', return_wfn=True, molecule=molecule, **kwargs)
            C = 1
            if E == 0.0:
                E = core.get_variable('SCF ITERATION ENERGY')
                C = 0

            if HOMO > 0:
                eps = wfn.epsilon_a()
                potentials.append(eps[HOMO - 1])
            else:
                eps = wfn.epsilon_b()
                potentials.append(eps[-HOMO - 1])

            Ns.append(Nintegral + occ - 1.0)
            energies.append(E)
            convs.append(C)

            core.set_global_option("FRAC_START", 2)
            core.set_global_option("FRAC_LOAD", True)
            core.set_global_option("FRAC_DIIS", frac_diis)
            core.set_global_option("GUESS", "READ")

        # Set the next charge/mult
        molecule.set_molecular_charge(charge)
        molecule.set_multiplicity(mult)

        # Determine HOMO
        print('DGAS: What ref should this point to?')
        #ref = core.legacy_wavefunction()
        eps_a = wfn.epsilon_a()
        eps_b = wfn.epsilon_b()
        if Na == Nb:
            HOMO = -Nb
        elif Nb == 0:
            HOMO = Na
        else:
            E_a = eps_a[int(Na - 1)]
            E_b = eps_b[int(Nb - 1)]
            if E_a >= E_b:
                HOMO = Na
            else:
                HOMO = -Nb

        stats.append("""\t%6d %6d %6d %6d %6d %6d\n""" % (Nintegral-1, Na, Nb, charge, mult, HOMO))

        if HOMO > 0:
            Na = Na - 1
        else:
            Nb = Nb - 1
        charge = charge + 1
        mult = Na - Nb + 1

    core.set_global_option("DF_INTS_IO", "NONE")

    # => Print the results out <= #
    E = {}
    core.print_out("""\n    ==> Fractional Occupation Nuke Results <==\n\n""")
    core.print_out("""\t%-11s %-24s %-24s %11s\n""" % ('N', 'Energy', 'HOMO Energy', 'Converged'))
    for k in range(len(Ns)):
        core.print_out("""\t%11.3E %24.16E %24.16E %11d\n""" % (Ns[k], energies[k], potentials[k], convs[k]))
        E[Ns[k]] = energies[k]

    core.print_out('\n')
    core.print_out("""\t%6s %6s %6s %6s %6s %6s\n""" % ('N', 'Na', 'Nb', 'Charge', 'Mult', 'HOMO'))
    for line in stats:
        core.print_out(line)

    core.print_out('\n\t"You shoot a nuke down a bug hole, you got a lot of dead bugs"\n')
    core.print_out('\t\t\t-Starship Troopers\n')

    # Drop the files out
    fh = open(traverse_filename, 'w')
    fh.write("""\t%-11s %-24s %-24s %11s\n""" % ('N', 'Energy', 'HOMO Energy', 'Converged'))
    for k in range(len(Ns)):
        fh.write("""\t%11.3E %24.16E %24.16E %11d\n""" % (Ns[k], energies[k], potentials[k], convs[k]))
    fh.close()

    fh = open(stats_filename, 'w')
    fh.write("""\t%6s %6s %6s %6s %6s %6s\n""" % ('N', 'Na', 'Nb', 'Charge', 'Mult', 'HOMO'))
    for line in stats:
        fh.write(line)
    fh.close()

    # Properly, should clone molecule but since not returned and easy to unblemish,
    molecule.set_molecular_charge(charge0)
    molecule.set_multiplicity(mult0)

    return E

def ip_fitting(molecule, omega_l, omega_r, **kwargs):
    kwargs = p4util.kwargs_lower(kwargs)

    # By default, zero the omega to 3 digits
    omega_tol = kwargs.get('omega_tolerance', 1.0E-3)

    # By default, do up to twenty iterations
    maxiter = kwargs.get('maxiter', 20)

    # By default, do not read previous 180 orbitals file
    read = False
    read180 = ''
    if 'read' in kwargs:
        read = True
        read180 = kwargs['read']

    # The molecule is required, and should be the neutral species
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
        core.set_global_option("GUESS", "READ")
        copy_file_to_scratch(read180, 'psi', 'ot', 180)
    old_guess = core.get_global_option("GUESS")
    core.set_global_option("DF_INTS_IO", "SAVE")
    core.print_out("""\n\t==> IP Fitting SCF: Burn-in <==\n""")
    E, wfn = energy('scf', return_wfn=True, molecule=molecule, **kwargs)
    core.set_global_option("DF_INTS_IO", "LOAD")

    # Determine HOMO, to determine mult1
    eps_a = wfn.epsilon_a()
    eps_b = wfn.epsilon_b()
    if Na == Nb:
        HOMO = -Nb
    elif Nb == 0:
        HOMO = Na
    else:
        E_a = eps_a[int(Na - 1)]
        E_b = eps_b[int(Nb - 1)]
        if E_a >= E_b:
            HOMO = Na
        else:
            HOMO = -Nb

    Na1 = Na;
    Nb1 = Nb;
    if HOMO > 0:
        Na1 = Na1 - 1;
    else:
        Nb1 = Nb1 - 1;

    charge1 = charge0 + 1;
    mult1 = Na1 - Nb1 + 1

    omegas = []
    E0s = []
    E1s = []
    kIPs = []
    IPs = []
    types = []

    # Right endpoint
    core.set_global_option('DFT_OMEGA', omega_r)

    # Neutral
    if read:
        core.set_global_option("GUESS", "READ")
        p4util.copy_file_to_scratch(read180, 'psi', 'ot', 180)

    molecule.set_molecular_charge(charge0)
    molecule.set_multiplicity(mult0)
    core.print_out("""\n\t==> IP Fitting SCF: Neutral, Right Endpoint <==\n""")
    E0r, wfn = energy('scf', return_wfn=True, molecule=molecule, **kwargs)
    eps_a = wfn.epsilon_a()
    eps_b = wfn.epsilon_b()
    E_HOMO = 0.0;
    if Nb == 0:
        E_HOMO = eps_a[int(Na - 1)]
    else:
        E_a = eps_a[int(Na - 1)]
        E_b = eps_b[int(Nb - 1)]
        if E_a >= E_b:
            E_HOMO = E_a
        else:
            E_HOMO = E_b
    E_HOMOr = E_HOMO
    core.IO.change_file_namespace(180, "ot", "neutral")

    # Cation
    if read:
        core.set_global_option("GUESS", "READ")
        p4util.copy_file_to_scratch(read180, 'psi', 'ot', 180)

    molecule.set_molecular_charge(charge1)
    molecule.set_multiplicity(mult1)
    core.print_out("""\n\t==> IP Fitting SCF: Cation, Right Endpoint <==\n""")
    E1r = energy('scf', molecule=molecule, **kwargs)
    core.IO.change_file_namespace(180, "ot", "cation")

    IPr = E1r - E0r;
    kIPr = -E_HOMOr;
    delta_r = IPr - kIPr;

    if IPr > kIPr:
        message = ("""\n***IP Fitting Error: Right Omega limit should have kIP > IP""")
        raise ValidationError(message)

    omegas.append(omega_r)
    types.append('Right Limit')
    E0s.append(E0r)
    E1s.append(E1r)
    IPs.append(IPr)
    kIPs.append(kIPr)

    # Use previous orbitals from here out
    core.set_global_option("GUESS", "READ")

    # Left endpoint
    core.set_global_option('DFT_OMEGA', omega_l)

    # Neutral
    core.IO.change_file_namespace(180, "neutral", "ot")
    molecule.set_molecular_charge(charge0)
    molecule.set_multiplicity(mult0)
    core.print_out("""\n\t==> IP Fitting SCF: Neutral, Left Endpoint <==\n""")
    E0l, wfn = energy('scf', return_wfn=True, molecule=molecule, **kwargs)
    eps_a = wfn.epsilon_a()
    eps_b = wfn.epsilon_b()
    E_HOMO = 0.0
    if Nb == 0:
        E_HOMO = eps_a[int(Na - 1)]
    else:
        E_a = eps_a[int(Na - 1)]
        E_b = eps_b[int(Nb - 1)]
        if E_a >= E_b:
            E_HOMO = E_a
        else:
            E_HOMO = E_b
    E_HOMOl = E_HOMO
    core.IO.change_file_namespace(180, "ot", "neutral")

    # Cation
    core.IO.change_file_namespace(180, "cation", "ot")
    molecule.set_molecular_charge(charge1)
    molecule.set_multiplicity(mult1)
    core.print_out("""\n\t==> IP Fitting SCF: Cation, Left Endpoint <==\n""")
    E1l = energy('scf', molecule=molecule, **kwargs)
    core.IO.change_file_namespace(180, "ot", "cation")

    IPl = E1l - E0l
    kIPl = -E_HOMOl
    delta_l = IPl - kIPl

    if IPl < kIPl:
        message = ("""\n***IP Fitting Error: Left Omega limit should have kIP < IP""")
        raise ValidationError(message)

    omegas.append(omega_l)
    types.append('Left Limit')
    E0s.append(E0l)
    E1s.append(E1l)
    IPs.append(IPl)
    kIPs.append(kIPl)

    converged = False
    repeat_l = 0
    repeat_r = 0
    step = 0
    while True:

        step = step + 1

        # Regula Falsi (modified)
        if repeat_l > 1:
            delta_l = delta_l / 2.0
        if repeat_r > 1:
            delta_r = delta_r / 2.0
        omega = - (omega_r - omega_l) / (delta_r - delta_l) * delta_l + omega_l
        core.set_global_option('DFT_OMEGA', omega)

        # Neutral
        core.IO.change_file_namespace(180, "neutral", "ot")
        molecule.set_molecular_charge(charge0)
        molecule.set_multiplicity(mult0)
        core.print_out("""\n\t==> IP Fitting SCF: Neutral, Omega = %11.3E <==\n""" % omega)
        E0, wfn = energy('scf', return_wfn=True, molecule=molecule, **kwargs)
        eps_a = wfn.epsilon_a()
        eps_b = wfn.epsilon_b()
        E_HOMO = 0.0
        if Nb == 0:
            E_HOMO = eps_a[int(Na - 1)]
        else:
            E_a = eps_a[int(Na - 1)]
            E_b = eps_b[int(Nb - 1)]
            if E_a >= E_b:
                E_HOMO = E_a
            else:
                E_HOMO = E_b
        core.IO.change_file_namespace(180, "ot", "neutral")

        # Cation
        core.IO.change_file_namespace(180, "cation", "ot")
        molecule.set_molecular_charge(charge1)
        molecule.set_multiplicity(mult1)
        core.print_out("""\n\t==> IP Fitting SCF: Cation, Omega = %11.3E <==\n""" % omega)
        E1 = energy('scf', molecule=molecule, **kwargs)
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
            repeat_l = repeat_l + 1
        else:
            omega_l = omega
            E0l = E0
            E1l = E1
            IPl = IP
            kIPl = kIP
            delta_l = delta
            repeat_l = 0;
            repeat_r = repeat_r + 1

        omegas.append(omega)
        types.append('Regula-Falsi')
        E0s.append(E0)
        E1s.append(E1)
        IPs.append(IP)
        kIPs.append(kIP)

        # Termination
        if (abs(omega_l - omega_r) < omega_tol or step > maxiter):
            converged = True
            break

    # Properly, should clone molecule but since not returned and easy to unblemish,
    molecule.set_molecular_charge(charge0)
    molecule.set_multiplicity(mult0)
    core.IO.set_default_namespace("")

    core.print_out("""\n\t==> IP Fitting Results <==\n\n""")

    core.print_out("""\t => Occupation Determination <= \n\n""")
    core.print_out("""\t          %6s %6s %6s %6s %6s %6s\n""" % ('N', 'Na', 'Nb', 'Charge', 'Mult', 'HOMO'))
    core.print_out("""\t Neutral: %6d %6d %6d %6d %6d %6d\n""" % (N, Na, Nb, charge0, mult0, HOMO))
    core.print_out("""\t Cation:  %6d %6d %6d %6d %6d\n\n""" % (N - 1, Na1, Nb1, charge1, mult1))

    core.print_out("""\t => Regula Falsi Iterations <=\n\n""")
    core.print_out("""\t%3s %11s %14s %14s %14s %s\n""" % ('N','Omega','IP','kIP','Delta','Type'))
    for k in range(len(omegas)):
        core.print_out("""\t%3d %11.3E %14.6E %14.6E %14.6E %s\n""" % 
                       (k + 1, omegas[k], IPs[k], kIPs[k], IPs[k] - kIPs[k], types[k]))
    if converged:
        core.print_out("""\n\tIP Fitting Converged\n""")
        core.print_out("""\tFinal omega = %14.6E\n""" % ((omega_l + omega_r) / 2))
        core.print_out("""\n\t"M,I. does the dying. Fleet just does the flying."\n""")
        core.print_out("""\t\t\t-Starship Troopers\n""")

    else:
        core.print_out("""\n\tIP Fitting did not converge!\n""")

    core.set_global_option("DF_INTS_IO", "NONE")
    core.set_global_option("GUESS", old_guess)
