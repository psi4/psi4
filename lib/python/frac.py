#
#@BEGIN LICENSE
#
# PSI4: an ab initio quantum chemistry software package
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
#@END LICENSE
#

import psi4
import os
import math
from molutil import *
from driver import *
from procutil import *
from util import *

# Scan from +1 electron to -1 electron
def frac_traverse(mol, **kwargs):
    kwargs = kwargs_lower(kwargs)

    # The molecule is required, and should be the neutral species
    mol.update_geometry()
    activate(mol)
    charge0 = mol.molecular_charge()
    mult0   = mol.multiplicity()

    chargep = charge0 + 1
    chargem = charge0 - 1

    # By default, the multiplicity of the cation/anion are mult0 + 1
    # These are overridden with the cation_mult and anion_mult kwargs
    multp = mult0 + 1
    multm = mult0 + 1
    if kwargs.has_key('cation_mult'):
        multp = kwargs['cation_mult']
    if kwargs.has_key('anion_mult'):
        multm = kwargs['anion_mult']

    # By default, we start the frac procedure on the 25th iteration
    # when not reading a previous guess
    frac_start = 25
    if kwargs.has_key('frac_start'):
        frac_start = kwargs['frac_start']

    # By default, we occupy by tenths of electrons
    LUMO_occs = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
    HOMO_occs  = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
    if kwargs.has_key('HOMO_occs'):
        HOMO_occs = kwargs['HOMO_occs']
    if kwargs.has_key('LUMO_occs'):
        LUMO_occs = kwargs['LUMO_occs']

    # By default, HOMO and LUMO are both in alpha
    Z = 0;
    for A in range(mol.natom()):
        Z += mol.Z(A)
    Z -= charge0
    if (Z%2):
        HOMO = Z/2+1
    else:
        HOMO = Z/2
    LUMO = HOMO+1
    if kwargs.has_key('HOMO'):
        HOMO = kwargs['HOMO']
    if kwargs.has_key('LUMO'):
        LUMO = kwargs['LUMO']

    # By default, DIIS in FRAC (1.0 occupation is always DIIS'd)
    frac_diis = True
    if kwargs.has_key('frac_diis'):
        frac_diis = kwargs['frac_diis']

    # By default, use the neutral orbitals as a guess for the anion
    neutral_guess = True
    if kwargs.has_key('neutral_guess'):
        neutral_guess = kwargs['neutral_guess']

    # By default, burn-in with UHF first, if UKS
    hf_guess = False
    if psi4.get_global_option('REFERENCE') == 'UKS':
        hf_guess = True
        if kwargs.has_key('hf_guess'):
            hf_guess = kwargs['hf_guess']

    # By default, re-guess at each N
    continuous_guess = False
    if kwargs.has_key('continuous_guess'):
         continuous_guess = kwargs['continuous_guess']

    # By default, drop the files to the molecule's name
    root = mol.name()
    if kwargs.has_key('filename'):
        root = kwargs['filename']
    traverse_filename = root + '.traverse.dat'
    # => Traverse <= #
    occs = []
    energies = []
    potentials = []
    convs = []

    # => Run the neutral for its orbitals, if requested <= #

    old_df_ints_io = psi4.get_global_option("DF_INTS_IO")
    psi4.set_global_option("DF_INTS_IO", "SAVE")

    old_guess = psi4.get_global_option("GUESS")
    if (neutral_guess):
        if (hf_guess):
            psi4.set_global_option("REFERENCE","UHF")
        energy('scf')
        psi4.set_global_option("GUESS", "READ")
        psi4.set_global_option("DF_INTS_IO", "LOAD")

    # => Run the anion first <= #

    mol.set_molecular_charge(chargem)
    mol.set_multiplicity(multm)

    # => Burn the anion in with hf, if requested <= #
    if (hf_guess):
        psi4.set_global_option("REFERENCE","UHF")
        energy('scf')
        psi4.set_global_option("REFERENCE","UKS")
        psi4.set_global_option("GUESS", "READ")
        psi4.set_global_option("DF_INTS_IO", "SAVE")

    psi4.set_global_option("FRAC_START", frac_start)
    psi4.set_global_option("FRAC_RENORMALIZE", True)
    psi4.set_global_option("FRAC_LOAD", False)

    for occ in LUMO_occs:

        psi4.set_global_option("FRAC_OCC", [LUMO])
        psi4.set_global_option("FRAC_VAL", [occ])

        E = energy('scf')
        C = 1
        if (E == 0.0):
            E = psi4.get_variable('SCF ITERATION ENERGY')
            C = 0

        if (LUMO > 0):
            ref = psi4.wavefunction()
            eps = ref.epsilon_a()
            potentials.append(eps[int(LUMO)-1])
        else:
            ref = psi4.wavefunction()
            eps = ref.epsilon_b()
            potentials.append(eps[-int(LUMO)-1])

        occs.append(occ)
        energies.append(E)
        convs.append(C)

        psi4.set_global_option("FRAC_START", 2)
        psi4.set_global_option("FRAC_LOAD", True)
        psi4.set_global_option("GUESS", "READ")
        psi4.set_global_option("FRAC_DIIS", frac_diis)
        psi4.set_global_option("DF_INTS_IO", "LOAD")


    # => Run the neutral next <= #

    mol.set_molecular_charge(charge0)
    mol.set_multiplicity(mult0)

    # Burn the neutral in with hf, if requested <= #

    if (not continuous_guess):
        psi4.set_global_option("GUESS", old_guess)
        if (hf_guess):
            psi4.set_global_option("FRAC_START", 0)
            psi4.set_global_option("REFERENCE","UHF")
            energy('scf')
            psi4.set_global_option("REFERENCE","UKS")
            psi4.set_global_option("GUESS", "READ")
        psi4.set_global_option("FRAC_LOAD", False)

    psi4.set_global_option("FRAC_START", frac_start)
    psi4.set_global_option("FRAC_RENORMALIZE", True)

    for occ in HOMO_occs:

        psi4.set_global_option("FRAC_OCC", [HOMO])
        psi4.set_global_option("FRAC_VAL", [occ])

        E = energy('scf')
        C = 1
        if (E == 0.0):
            E = psi4.get_variable('SCF ITERATION ENERGY')
            C = 0

        if (LUMO > 0):
            ref = psi4.wavefunction()
            eps = ref.epsilon_a()
            potentials.append(eps[int(HOMO)-1])
        else:
            ref = psi4.wavefunction()
            eps = ref.epsilon_b()
            potentials.append(eps[-int(HOMO)-1])

        occs.append(occ - 1.0)
        energies.append(E)
        convs.append(C)

        psi4.set_global_option("FRAC_START", 2)
        psi4.set_global_option("FRAC_LOAD", True)
        psi4.set_global_option("GUESS", "READ")
        psi4.set_global_option("FRAC_DIIS", frac_diis)
        psi4.set_global_option("DF_INTS_IO", "LOAD")

    psi4.set_global_option("DF_INTS_IO", old_df_ints_io)

    # => Print the results out <= #
    E = {}
    psi4.print_out('\n    ==> Fractional Occupation Traverse Results <==\n\n')
    psi4.print_out('\t%-11s %-24s %-24s %11s\n' %('N', 'Energy', 'HOMO Energy', 'Converged'))
    for k in range(len(occs)):
        psi4.print_out('\t%11.3E %24.16E %24.16E %11d\n' % (occs[k], energies[k], potentials[k], convs[k]))
        E[occs[k]] = energies[k]

    psi4.print_out('\n\t"You trying to be a hero Watkins?"\n')
    psi4.print_out('\t"Just trying to kill some bugs sir!"\n')
    psi4.print_out('\t\t\t-Starship Troopers\n')

    # Drop the files out
    fh = open(traverse_filename, 'w')
    fh.write('\t%-11s %-24s %-24s %11s\n' %('N', 'Energy', 'HOMO Energy', 'Converged'))
    for k in range(len(occs)):
        fh.write('\t%11.3E %24.16E %24.16E %11d\n' % (occs[k], energies[k], potentials[k], convs[k]))
    fh.close()

    return E

# Pull all the electrons out, one at a time
def frac_nuke(mol, **kwargs):
    kwargs = kwargs_lower(kwargs)

    # The molecule is required, and should be the neutral species
    mol.update_geometry()
    activate(mol)
    charge0 = mol.molecular_charge()
    mult0   = mol.multiplicity()

    # By default, we start the frac procedure on the 25th iteration
    # when not reading a previous guess
    frac_start = 25
    if kwargs.has_key('frac_start'):
        frac_start = kwargs['frac_start']

    # By default, we occupy by tenths of electrons
    foccs = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
    if kwargs.has_key('foccs'):
        foccs = kwargs['foccs']

    # By default, HOMO and LUMO are both in alpha
    N = 0;
    for A in range(mol.natom()):
        N += mol.Z(A)
    N -= charge0
    N  = int(N)
    Nb = int((N - mult0 + 1)/2)
    Na = int(N - Nb)

    charge = charge0
    mult = mult0

    # By default, nuke all the electrons
    Nmin = 0;
    if (kwargs.has_key('nmax')):
        Nmin = N - int(kwargs['nmax'])

    # By default, DIIS in FRAC (1.0 occupation is always DIIS'd)
    frac_diis = True
    if kwargs.has_key('frac_diis'):
        frac_diis = kwargs['frac_diis']

    # By default, drop the files to the molecule's name
    root = mol.name()
    if kwargs.has_key('filename'):
        root = kwargs['filename']
    traverse_filename = root + '.traverse.dat'
    stats_filename = root + '.stats.dat'

    # => Traverse <= #
    psi4.set_global_option("DF_INTS_IO", "SAVE")

    Ns = []
    energies = []
    potentials = []
    convs = []
    stats = []

    # Run one SCF to burn things in
    energy('scf')

    # Determine HOMO
    ref = psi4.wavefunction()
    eps_a = ref.epsilon_a()
    eps_b = ref.epsilon_b()
    if (Na == Nb):
        HOMO = -Nb
    elif (Nb == 0):
        HOMO = Na
    else:
        E_a = eps_a[int(Na - 1)]
        E_b = eps_b[int(Nb - 1)]
        if (E_a >= E_b):
            HOMO = Na
        else:
            HOMO = -Nb

    stats.append('\t%6d %6d %6d %6d %6d %6d\n' %(N, Na, Nb, charge, mult, HOMO))

    if (HOMO > 0):
        Na = Na - 1
    else:
        Nb = Nb - 1
    charge = charge + 1
    mult = Na - Nb + 1

    psi4.set_global_option("DF_INTS_IO", "LOAD")
    psi4.set_global_option("FRAC_START", frac_start)
    psi4.set_global_option("FRAC_RENORMALIZE", True)

    # Nuke 'em Rico!
    for Nintegral in range(N,Nmin,-1):

        # Nuke the current HOMO
        for occ in foccs:

            psi4.set_global_option("FRAC_OCC", [HOMO])
            psi4.set_global_option("FRAC_VAL", [occ])

            E = energy('scf')
            C = 1
            if (E == 0.0):
                E = psi4.get_variable('SCF ITERATION ENERGY')
                C = 0

            if (HOMO > 0):
                ref = psi4.wavefunction()
                eps = ref.epsilon_a()
                potentials.append(eps[HOMO-1])
            else:
                ref = psi4.wavefunction()
                eps = ref.epsilon_b()
                potentials.append(eps[-HOMO-1])

            Ns.append(Nintegral + occ - 1.0)
            energies.append(E)
            convs.append(C)

            psi4.set_global_option("FRAC_START", 2)
            psi4.set_global_option("FRAC_LOAD", True)
            psi4.set_global_option("FRAC_DIIS", frac_diis)
            psi4.set_global_option("GUESS", "READ")

        # Set the next charge/mult
        mol.set_molecular_charge(charge)
        mol.set_multiplicity(mult)

        # Determine HOMO
        ref = psi4.wavefunction()
        eps_a = ref.epsilon_a()
        eps_b = ref.epsilon_b()
        if (Na == Nb):
            HOMO = -Nb
        elif (Nb == 0):
            HOMO = Na
        else:
            E_a = eps_a[int(Na - 1)]
            E_b = eps_b[int(Nb - 1)]
            if (E_a >= E_b):
                HOMO = Na
            else:
                HOMO = -Nb

        stats.append('\t%6d %6d %6d %6d %6d %6d\n' %(Nintegral-1, Na, Nb, charge, mult, HOMO))

        if (HOMO > 0):
            Na = Na - 1
        else:
            Nb = Nb - 1
        charge = charge + 1
        mult = Na - Nb + 1

    psi4.set_global_option("DF_INTS_IO", "NONE")

    # => Print the results out <= #
    E = {}
    psi4.print_out('\n    ==> Fractional Occupation Nuke Results <==\n\n')
    psi4.print_out('\t%-11s %-24s %-24s %11s\n' %('N', 'Energy', 'HOMO Energy', 'Converged'))
    for k in range(len(Ns)):
        psi4.print_out('\t%11.3E %24.16E %24.16E %11d\n' % (Ns[k], energies[k], potentials[k], convs[k]))
        E[Ns[k]] = energies[k]

    psi4.print_out('\n')
    psi4.print_out('\t%6s %6s %6s %6s %6s %6s\n' %('N', 'Na', 'Nb', 'Charge', 'Mult', 'HOMO'))
    for line in stats:
        psi4.print_out(line)

    psi4.print_out('\n\t"You shoot a nuke down a bug hole, you got a lot of dead bugs"\n')
    psi4.print_out('\t\t\t-Starship Troopers\n')

    # Drop the files out
    fh = open(traverse_filename, 'w')
    fh.write('\t%-11s %-24s %-24s %11s\n' %('N', 'Energy', 'HOMO Energy', 'Converged'))
    for k in range(len(Ns)):
        fh.write('\t%11.3E %24.16E %24.16E %11d\n' % (Ns[k], energies[k], potentials[k], convs[k]))
    fh.close()

    fh = open(stats_filename, 'w')
    fh.write('\t%6s %6s %6s %6s %6s %6s\n' %('N', 'Na', 'Nb', 'Charge', 'Mult', 'HOMO'))
    for line in stats:
        fh.write(line)
    fh.close()

    return E

def ip_fitting(mol, omega_l, omega_r, **kwargs):
    kwargs = kwargs_lower(kwargs)

    # By default, zero the omega to 3 digits
    omega_tol = 1.0E-3;
    if (kwargs.has_key('omega_tolerance')):
        omega_tol = kwargs['omega_tolerance']

    # By default, do up to twenty iterations
    maxiter = 20;
    if (kwargs.has_key('maxiter')):
        maxiter = kwargs['maxiter']

    # By default, do not read previous 180 orbitals file
    read = False;
    read180 = ''
    if (kwargs.has_key('read')):
        read = True;
        read180 = kwargs['read']

    # The molecule is required, and should be the neutral species
    mol.update_geometry()
    activate(mol)
    charge0 = mol.molecular_charge()
    mult0   = mol.multiplicity()

    # How many electrons are there?
    N = 0;
    for A in range(mol.natom()):
        N += mol.Z(A)
    N -= charge0
    N  = int(N)
    Nb = int((N - mult0 + 1)/2)
    Na = int(N - Nb)

    # Work in the ot namespace for this procedure
    psi4.IO.set_default_namespace("ot")

    # Burn in to determine orbital eigenvalues
    if (read):
        psi4.set_global_option("GUESS", "READ")
        copy_file_to_scratch(read180, 'psi', 'ot', 180)
    old_guess = psi4.get_global_option("GUESS")
    psi4.set_global_option("DF_INTS_IO", "SAVE")
    psi4.print_out('\n\t==> IP Fitting SCF: Burn-in <==\n')
    energy('scf')
    psi4.set_global_option("DF_INTS_IO", "LOAD")

    # Determine HOMO, to determine mult1
    ref = psi4.wavefunction()
    eps_a = ref.epsilon_a()
    eps_b = ref.epsilon_b()
    if (Na == Nb):
        HOMO = -Nb
    elif (Nb == 0):
        HOMO = Na
    else:
        E_a = eps_a[int(Na - 1)]
        E_b = eps_b[int(Nb - 1)]
        if (E_a >= E_b):
            HOMO = Na
        else:
            HOMO = -Nb

    Na1 = Na;
    Nb1 = Nb;
    if (HOMO > 0):
        Na1 = Na1-1;
    else:
        Nb1 = Nb1-1;

    charge1 = charge0 + 1;
    mult1 = Na1 - Nb1 + 1

    omegas = [];
    E0s = [];
    E1s = [];
    kIPs = [];
    IPs = [];
    types = [];

    # Right endpoint
    psi4.set_global_option('DFT_OMEGA',omega_r)

    # Neutral
    if (read):
        psi4.set_global_option("GUESS", "READ")
        copy_file_to_scratch(read180, 'psi', 'ot', 180)

    mol.set_molecular_charge(charge0)
    mol.set_multiplicity(mult0)
    psi4.print_out('\n\t==> IP Fitting SCF: Neutral, Right Endpoint <==\n')
    E0r = energy('scf')
    ref = psi4.wavefunction()
    eps_a = ref.epsilon_a()
    eps_b = ref.epsilon_b()
    E_HOMO = 0.0;
    if (Nb == 0):
        E_HOMO = eps_a[int(Na-1)]
    else:
        E_a = eps_a[int(Na - 1)]
        E_b = eps_b[int(Nb - 1)]
        if (E_a >= E_b):
            E_HOMO = E_a;
        else:
            E_HOMO = E_b;
    E_HOMOr = E_HOMO;
    psi4.IO.change_file_namespace(180,"ot","neutral")

    # Cation
    if (read):
        psi4.set_global_option("GUESS", "READ")
        copy_file_to_scratch(read180, 'psi', 'ot', 180)

    mol.set_molecular_charge(charge1)
    mol.set_multiplicity(mult1)
    psi4.print_out('\n\t==> IP Fitting SCF: Cation, Right Endpoint <==\n')
    E1r = energy('scf')
    psi4.IO.change_file_namespace(180,"ot","cation")

    IPr = E1r - E0r;
    kIPr = -E_HOMOr;
    delta_r = IPr - kIPr;

    if (IPr > kIPr):
        psi4.print_out('\n***IP Fitting Error: Right Omega limit should have kIP > IP')
        sys.exit(1)

    omegas.append(omega_r)
    types.append('Right Limit')
    E0s.append(E0r)
    E1s.append(E1r)
    IPs.append(IPr)
    kIPs.append(kIPr)

    # Use previous orbitals from here out
    psi4.set_global_option("GUESS","READ")

    # Left endpoint
    psi4.set_global_option('DFT_OMEGA',omega_l)

    # Neutral
    psi4.IO.change_file_namespace(180,"neutral","ot")
    mol.set_molecular_charge(charge0)
    mol.set_multiplicity(mult0)
    psi4.print_out('\n\t==> IP Fitting SCF: Neutral, Left Endpoint <==\n')
    E0l = energy('scf')
    ref = psi4.wavefunction()
    eps_a = ref.epsilon_a()
    eps_b = ref.epsilon_b()
    E_HOMO = 0.0;
    if (Nb == 0):
        E_HOMO = eps_a[int(Na-1)]
    else:
        E_a = eps_a[int(Na - 1)]
        E_b = eps_b[int(Nb - 1)]
        if (E_a >= E_b):
            E_HOMO = E_a;
        else:
            E_HOMO = E_b;
    E_HOMOl = E_HOMO;
    psi4.IO.change_file_namespace(180,"ot","neutral")

    # Cation
    psi4.IO.change_file_namespace(180,"cation","ot")
    mol.set_molecular_charge(charge1)
    mol.set_multiplicity(mult1)
    psi4.print_out('\n\t==> IP Fitting SCF: Cation, Left Endpoint <==\n')
    E1l = energy('scf')
    psi4.IO.change_file_namespace(180,"ot","cation")

    IPl = E1l - E0l;
    kIPl = -E_HOMOl;
    delta_l = IPl - kIPl;

    if (IPl < kIPl):
        psi4.print_out('\n***IP Fitting Error: Left Omega limit should have kIP < IP')
        sys.exit(1)

    omegas.append(omega_l)
    types.append('Left Limit')
    E0s.append(E0l)
    E1s.append(E1l)
    IPs.append(IPl)
    kIPs.append(kIPl)

    converged = False
    repeat_l = 0;
    repeat_r = 0;
    step = 0;
    while True:

        step = step + 1;

        # Regula Falsi (modified)
        if (repeat_l > 1):
            delta_l = delta_l / 2.0;
        if (repeat_r > 1):
            delta_r = delta_r / 2.0;
        omega = - (omega_r - omega_l) / (delta_r - delta_l) * delta_l + omega_l;
        psi4.set_global_option('DFT_OMEGA',omega)

        # Neutral
        psi4.IO.change_file_namespace(180,"neutral","ot")
        mol.set_molecular_charge(charge0)
        mol.set_multiplicity(mult0)
        psi4.print_out('\n\t==> IP Fitting SCF: Neutral, Omega = %11.3E <==\n' % omega)
        E0 = energy('scf')
        ref = psi4.wavefunction()
        eps_a = ref.epsilon_a()
        eps_b = ref.epsilon_b()
        E_HOMO = 0.0;
        if (Nb == 0):
            E_HOMO = eps_a[int(Na-1)]
        else:
            E_a = eps_a[int(Na - 1)]
            E_b = eps_b[int(Nb - 1)]
            if (E_a >= E_b):
                E_HOMO = E_a;
            else:
                E_HOMO = E_b;
        psi4.IO.change_file_namespace(180,"ot","neutral")

        # Cation
        psi4.IO.change_file_namespace(180,"cation","ot")
        mol.set_molecular_charge(charge1)
        mol.set_multiplicity(mult1)
        psi4.print_out('\n\t==> IP Fitting SCF: Cation, Omega = %11.3E <==\n' % omega)
        E1 = energy('scf')
        psi4.IO.change_file_namespace(180,"ot","cation")

        IP = E1 - E0;
        kIP = -E_HOMO;
        delta = IP - kIP;

        if (kIP > IP):
            omega_r = omega
            E0r = E0
            E1r = E1
            IPr = IP
            kIPr = kIP
            delta_r = delta
            repeat_r = 0;
            repeat_l = repeat_l + 1;
        else:
            omega_l = omega
            E0l = E0
            E1l = E1
            IPl = IP
            kIPl = kIP
            delta_l = delta
            repeat_l = 0;
            repeat_r = repeat_r + 1;

        omegas.append(omega)
        types.append('Regula-Falsi')
        E0s.append(E0)
        E1s.append(E1)
        IPs.append(IP)
        kIPs.append(kIP)

        # Termination
        if (abs(omega_l - omega_r) < omega_tol or step > maxiter):
            converged = True;
            break

    psi4.IO.set_default_namespace("")

    psi4.print_out('\n\t==> IP Fitting Results <==\n\n')

    psi4.print_out('\t => Occupation Determination <= \n\n')
    psi4.print_out('\t          %6s %6s %6s %6s %6s %6s\n' %('N', 'Na', 'Nb', 'Charge', 'Mult', 'HOMO'))
    psi4.print_out('\t Neutral: %6d %6d %6d %6d %6d %6d\n' %(N, Na, Nb, charge0, mult0, HOMO))
    psi4.print_out('\t Cation:  %6d %6d %6d %6d %6d\n\n' %(N-1, Na1, Nb1, charge1, mult1))

    psi4.print_out('\t => Regula Falsi Iterations <=\n\n')
    psi4.print_out('\t%3s %11s %14s %14s %14s %s\n' % ('N','Omega','IP','kIP','Delta','Type'))
    for k in range(len(omegas)):
        psi4.print_out('\t%3d %11.3E %14.6E %14.6E %14.6E %s\n' % (k+1,omegas[k],IPs[k],kIPs[k],IPs[k] - kIPs[k], types[k]))
    if (converged):
        psi4.print_out('\n\tIP Fitting Converged\n')
        psi4.print_out('\tFinal omega = %14.6E\n' % ((omega_l + omega_r) / 2))
        psi4.print_out('\n\t"M,I. does the dying. Fleet just does the flying."\n')
        psi4.print_out('\t\t\t-Starship Troopers\n')

    else:
        psi4.print_out('\n\tIP Fitting did not converge!\n')

    psi4.set_global_option("DF_INTS_IO", "NONE")
    psi4.set_global_option("GUESS", old_guess)

