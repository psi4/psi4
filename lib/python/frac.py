import PsiMod
import os
import input
import math
from molecule import * 
from driver import * 
from procutil import *

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
    if PsiMod.get_global_option('REFERENCE') == 'UKS':
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
    convs = []

    # => Run the neutral for its orbitals, if requested <= #

    old_df_ints_io = PsiMod.get_global_option("DF_INTS_IO")
    PsiMod.set_global_option("DF_INTS_IO", "SAVE")

    old_guess = PsiMod.get_global_option("GUESS")
    if (neutral_guess):
        if (hf_guess):
            PsiMod.set_global_option("REFERENCE","UHF")
        energy('scf')
        PsiMod.set_global_option("GUESS", "READ")
        PsiMod.set_global_option("DF_INTS_IO", "LOAD")

    # => Run the anion first <= #

    mol.set_molecular_charge(chargem)
    mol.set_multiplicity(multm)
    
    # => Burn the anion in with hf, if requested <= #
    if (hf_guess):
        PsiMod.set_global_option("REFERENCE","UHF")
        energy('scf')
        PsiMod.set_global_option("REFERENCE","UKS")
        PsiMod.set_global_option("GUESS", "READ")
        PsiMod.set_global_option("DF_INTS_IO", "SAVE")

    PsiMod.set_global_option("FRAC_START", frac_start)
    PsiMod.set_global_option("FRAC_RENORMALIZE", True)
    PsiMod.set_global_option("FRAC_LOAD", False)

    for occ in LUMO_occs:
       
        PsiMod.set_global_option("FRAC_OCC", [LUMO])
        PsiMod.set_global_option("FRAC_VAL", [occ]) 

        E = energy('scf')
        C = 1
        if (E == 0.0):
            E = PsiMod.get_variable('SCF ITERATION ENERGY')
            C = 0

        occs.append(occ)
        energies.append(E)
        convs.append(C)

        PsiMod.set_global_option("FRAC_START", 2)
        PsiMod.set_global_option("FRAC_LOAD", True)
        PsiMod.set_global_option("GUESS", "READ")
        PsiMod.set_global_option("FRAC_DIIS", frac_diis)
        PsiMod.set_global_option("DF_INTS_IO", "LOAD")


    # => Run the neutral next <= #

    mol.set_molecular_charge(charge0)
    mol.set_multiplicity(mult0)

    # Burn the neutral in with hf, if requested <= #

    if (not continuous_guess):
        PsiMod.set_global_option("GUESS", old_guess)
        if (hf_guess):
            PsiMod.set_global_option("FRAC_START", 0)
            PsiMod.set_global_option("REFERENCE","UHF")
            energy('scf')
            PsiMod.set_global_option("REFERENCE","UKS")
            PsiMod.set_global_option("GUESS", "READ")
        PsiMod.set_global_option("FRAC_LOAD", False)

    PsiMod.set_global_option("FRAC_START", frac_start)
    PsiMod.set_global_option("FRAC_RENORMALIZE", True)

    for occ in HOMO_occs:
       
        PsiMod.set_global_option("FRAC_OCC", [HOMO])
        PsiMod.set_global_option("FRAC_VAL", [occ]) 

        E = energy('scf')
        C = 1
        if (E == 0.0):
            E = PsiMod.get_variable('SCF ITERATION ENERGY')
            C = 0

        occs.append(occ - 1.0)
        energies.append(E)
        convs.append(C)

        PsiMod.set_global_option("FRAC_START", 2)
        PsiMod.set_global_option("FRAC_LOAD", True)
        PsiMod.set_global_option("GUESS", "READ")
        PsiMod.set_global_option("FRAC_DIIS", frac_diis)
        PsiMod.set_global_option("DF_INTS_IO", "LOAD")

    PsiMod.set_global_option("DF_INTS_IO", old_df_ints_io)

    # => Print the results out <= #
    E = {}
    PsiMod.print_out('\n    ==> Fractional Occupation Traverse Results <==\n\n')
    PsiMod.print_out('\t%-11s %-24s %11s\n' %('Delta', 'Energy', 'Converged'))
    for k in range(len(occs)):
        PsiMod.print_out('\t%11.3E %24.16E %11d\n' % (occs[k], energies[k], convs[k]))
        E[occs[k]] = energies[k]

    PsiMod.print_out('\n\t"You trying to be a hero Watkins?"\n')
    PsiMod.print_out('\t"Just trying to kill some bugs sir!"\n')
    PsiMod.print_out('\t\t\t-Starship Troopers\n')

    # Drop the files out
    fh = open(traverse_filename, 'w')
    fh.write('\t%-11s %-24s %11s\n' %('N', 'Energy', 'Converged'))
    for k in range(len(occs)):
        fh.write('\t%11.3E %24.16E %11d\n' % (occs[k], energies[k], convs[k]))
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
    PsiMod.set_global_option("DF_INTS_IO", "SAVE")

    Ns = []
    energies = []
    convs = []
    stats = []

    # Run one SCF to burn things in
    energy('scf')

    # Determine HOMO
    ref = PsiMod.reference_wavefunction()
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

    PsiMod.set_global_option("DF_INTS_IO", "LOAD")
    PsiMod.set_global_option("FRAC_START", frac_start)
    PsiMod.set_global_option("FRAC_RENORMALIZE", True)

    # Nuke 'em Rico! 
    for Nintegral in range(N,0,-1):

        # Nuke the current HOMO
        for occ in foccs:

            PsiMod.set_global_option("FRAC_OCC", [HOMO])
            PsiMod.set_global_option("FRAC_VAL", [occ]) 

            E = energy('scf')
            C = 1
            if (E == 0.0):
                E = PsiMod.get_variable('SCF ITERATION ENERGY')
                C = 0

            Ns.append(Nintegral + occ - 1.0)
            energies.append(E)
            convs.append(C)

            PsiMod.set_global_option("FRAC_START", 2)
            PsiMod.set_global_option("FRAC_LOAD", True)
            PsiMod.set_global_option("FRAC_DIIS", frac_diis)
            PsiMod.set_global_option("GUESS", "READ")

        # Set the next charge/mult
        mol.set_molecular_charge(charge)    
        mol.set_multiplicity(mult)

        # Determine HOMO
        ref = PsiMod.reference_wavefunction()
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
    
    PsiMod.set_global_option("DF_INTS_IO", "NONE")

    # => Print the results out <= #
    E = {}
    PsiMod.print_out('\n    ==> Fractional Occupation Nuke Results <==\n\n')
    PsiMod.print_out('\t%-11s %-24s %11s\n' %('N', 'Energy', 'Converged'))
    for k in range(len(Ns)):
        PsiMod.print_out('\t%11.3E %24.16E %11d\n' % (Ns[k], energies[k], convs[k]))
        E[Ns[k]] = energies[k]

    PsiMod.print_out('\n')
    PsiMod.print_out('\t%6s %6s %6s %6s %6s %6s\n' %('N', 'Na', 'Nb', 'Charge', 'Mult', 'HOMO'))
    for line in stats: 
        PsiMod.print_out(line)

    PsiMod.print_out('\n\t"You shoot a nuke down a bug hole, you got a lot of dead bugs"\n')
    PsiMod.print_out('\t\t\t-Starship Troopers\n')

    # Drop the files out
    fh = open(traverse_filename, 'w')
    fh.write('\t%-11s %-24s %11s\n' %('N', 'Energy', 'Converged'))
    for k in range(len(Ns)):
        fh.write('\t%11.3E %24.16E %11d\n' % (Ns[k], energies[k], convs[k]))
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
    PsiMod.IO.set_default_namespace("ot")

    # Burn in to determine orbital eigenvalues
    old_guess = PsiMod.get_global_option("GUESS")
    PsiMod.set_global_option("DF_INTS_IO", "SAVE")
    PsiMod.print_out('\n\t==> IP Fitting SCF: Burn-in <==\n')
    energy('scf')
    PsiMod.set_global_option("DF_INTS_IO", "LOAD")

    # Determine HOMO, to determine mult1
    ref = PsiMod.reference_wavefunction()
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
    PsiMod.set_global_option('DFT_OMEGA',omega_r)

    # Neutral
    mol.set_molecular_charge(charge0)
    mol.set_multiplicity(mult0)
    PsiMod.print_out('\n\t==> IP Fitting SCF: Neutral, Right Endpoint <==\n')
    E0r = energy('scf')
    ref = PsiMod.reference_wavefunction()
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
    PsiMod.IO.change_file_namespace(180,"ot","neutral")
    
    # Cation
    mol.set_molecular_charge(charge1)
    mol.set_multiplicity(mult1)
    PsiMod.print_out('\n\t==> IP Fitting SCF: Cation, Right Endpoint <==\n')
    E1r = energy('scf')
    PsiMod.IO.change_file_namespace(180,"ot","cation")

    IPr = E1r - E0r;
    kIPr = -E_HOMOr;
    delta_r = IPr - kIPr;

    if (IPr > kIPr):
        PsiMod.print_out('\n***IP Fitting Error: Right Omega limit should have kIP > IP')
        sys,exit(1)

    omegas.append(omega_r)
    types.append('Right Limit')
    E0s.append(E0r)
    E1s.append(E1r)
    IPs.append(IPr)
    kIPs.append(kIPr)

    # Use previous orbitals from here out
    PsiMod.set_global_option("GUESS","READ")
    
    # Left endpoint
    PsiMod.set_global_option('DFT_OMEGA',omega_l)

    # Neutral
    PsiMod.IO.change_file_namespace(180,"neutral","ot")
    mol.set_molecular_charge(charge0)
    mol.set_multiplicity(mult0)
    PsiMod.print_out('\n\t==> IP Fitting SCF: Neutral, Left Endpoint <==\n')
    E0l = energy('scf')
    ref = PsiMod.reference_wavefunction()
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
    PsiMod.IO.change_file_namespace(180,"ot","neutral")
    
    # Cation
    PsiMod.IO.change_file_namespace(180,"cation","ot")
    mol.set_molecular_charge(charge1)
    mol.set_multiplicity(mult1)
    PsiMod.print_out('\n\t==> IP Fitting SCF: Cation, Left Endpoint <==\n')
    E1l = energy('scf')
    PsiMod.IO.change_file_namespace(180,"ot","cation")

    IPl = E1l - E0l;
    kIPl = -E_HOMOl;
    delta_l = IPl - kIPl;

    if (IPl < kIPl):
        PsiMod.print_out('\n***IP Fitting Error: Left Omega limit should have kIP < IP')
        sys,exit(1)

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
        PsiMod.set_global_option('DFT_OMEGA',omega)

        # Neutral
        PsiMod.IO.change_file_namespace(180,"neutral","ot")
        mol.set_molecular_charge(charge0)
        mol.set_multiplicity(mult0)
        PsiMod.print_out('\n\t==> IP Fitting SCF: Neutral, Omega = %11.3E <==\n' % omega)
        E0 = energy('scf')
        ref = PsiMod.reference_wavefunction()
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
        PsiMod.IO.change_file_namespace(180,"ot","neutral")
        
        # Cation
        PsiMod.IO.change_file_namespace(180,"cation","ot")
        mol.set_molecular_charge(charge1)
        mol.set_multiplicity(mult1)
        PsiMod.print_out('\n\t==> IP Fitting SCF: Cation, Omega = %11.3E <==\n' % omega)
        E1 = energy('scf')
        PsiMod.IO.change_file_namespace(180,"ot","cation")

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

    PsiMod.IO.set_default_namespace("")

    PsiMod.print_out('\n\t==> IP Fitting Results <==\n\n')

    PsiMod.print_out('\t => Occupation Determination <= \n\n')
    PsiMod.print_out('\t          %6s %6s %6s %6s %6s %6s\n' %('N', 'Na', 'Nb', 'Charge', 'Mult', 'HOMO'))
    PsiMod.print_out('\t Neutral: %6d %6d %6d %6d %6d %6d\n' %(N, Na, Nb, charge0, mult0, HOMO))
    PsiMod.print_out('\t Cation:  %6d %6d %6d %6d %6d\n\n' %(N-1, Na1, Nb1, charge1, mult1))

    PsiMod.print_out('\t => Regula Falsi Iterations <=\n\n')
    PsiMod.print_out('\t%3s %11s %14s %14s %14s %s\n' % ('N','Omega','IP','kIP','Delta','Type'))
    for k in range(len(omegas)):
        PsiMod.print_out('\t%3d %11.3E %14.6E %14.6E %14.6E %s\n' % (k+1,omegas[k],IPs[k],kIPs[k],IPs[k] - kIPs[k], types[k]))
    if (converged):
        PsiMod.print_out('\n\tIP Fitting Converged\n')
        PsiMod.print_out('\tFinal omega = %14.6E\n' % ((omega_l + omega_r) / 2))
    else:
        PsiMod.print_out('\n\tIP Fitting did not converge!\n')
    
    PsiMod.set_global_option("DF_INTS_IO", "NONE")
    PsiMod.set_global_option("GUESS", old_guess)
