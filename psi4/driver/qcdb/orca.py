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
from __future__ import print_function
from collections import defaultdict
from .pdict import PreservingDict
from .molecule import Molecule
from .physconst import *


def harvest(p4Mol, orca_out, **largs):
    """Harvest variables, gradient, and the molecule from the output and other
    files
    """

    # Split into lines as it is much easier to find what is needed
    out_lines = orca_out.split('\n')

    mol = harvest_molecule_from_outfile(out_lines)

    file_name = "NONE"
    grad = harvest_engrad(file_name)

    # Harvest energies and properties from the output file
    psivar = PreservingDict()
    harvest_scf_from_outfile(out_lines, psivar)
    harvest_dipole(out_lines, psivar)
    harvest_mp2(out_lines, psivar)
    harvest_coupled_cluster(out_lines, psivar)

    return psivar, grad, mol


def muster_memory(mem):
    """Transform input *mem* in MB into psi4-type options for orca.

    """
    # prepare memory keywords to be set as c-side keywords
    options = defaultdict(lambda: defaultdict(dict))
    options['ORCA']['ORCA_MAXCORE']['value'] = int(mem)
    text = "%MaxCore {}\n".format(options['ORCA']['ORCA_MAXCORE']['value'])

    for item in options['ORCA']:
        options['ORCA'][item]['clobber'] = True

    return text, options


def muster_modelchem(name, dertype):
    """Transform calculation method *name* and derivative level *dertype*
    into options for orca. While deliberately requested pieces,
    generally orca__orca_deriv_level and orca__orca_calc_level,
    are set to complain if contradicted ('clobber' set to True), other
    'recommended' settings, can be countermanded by keywords in input file
    ('clobber' set to False). Occasionally, we want these pieces to actually
    overcome keywords in input file ('superclobber' set to True).
    """
    text = ''
    lowername = name.lower()
    options = defaultdict(lambda: defaultdict(dict))

    if dertype == 0:
        options['ORCA']['ORCA_RUNTYP']['value'] = 'ENERGY'
    elif dertype == 1:
        options['ORCA']['ORCA_RUNTYP']['value'] = 'ENGRAD'
    else:
        raise ValidationError("Requested Orca dertype {} is not available."
                              .format(dertype))

    if lowername == 'orca':
        pass
    elif lowername == 'orca-b3lyp':
        options['ORCA']['ORCA_FUNCTIONAL']['value'] = 'B3LYP_G'
    elif lowername == 'orca-mp2':
        options['ORCA']['ORCA_CALC_LEVEL']['value'] = 'MP2'
    elif lowername == 'orca-ccsd':
        options['ORCA']['ORCA_CALC_LEVEL']['value'] = 'CCSD'
    elif lowername == 'orca-ccsd(t)':
        options['ORCA']['ORCA_CALC_LEVEL']['value'] = 'CCSD(T)'
    else:
        raise ValidationError("Requested Orca computational methods {} is not "
                              "available." .format(lowername))

    # Set clobbering
    if 'ORCA_RUNTYP' in options['ORCA']:
        options['ORCA']['ORCA_RUNTYP']['clobber'] = True
        options['ORCA']['ORCA_RUNTYP']['superclobber'] = True
    if 'ORCA_FUNCTIONAL' in options['ORCA']:
        options['ORCA']['ORCA_FUNCTIONAL']['clobber'] = True
        options['ORCA']['ORCA_FUNCTIONAL']['superclobber'] = True

    return text, options


def orca_list():
    """Return an array of Orca methods with energies. Appended
    to procedures['energy'].

    """
    val = []
    val.append('orca')
    val.append('orca-b3lyp')
    return val


def orca_gradient_list():
    """Return an array of Orca methods with analytical gradients.
    Appended to procedures['gradient'].

    """
    val = []
    val.append('oc-b3lyp')
    return val


def harvest_molecule_from_outfile(lines):
    """Return a molecule of the last geometry"""
    """Sample molecule block"""
    #----------------------------
    #CARTESIAN COORDINATES (A.U.)
    #----------------------------
    #  NO LB      ZA    FRAG    MASS        X           Y           Z
    #   0 O     8.0000    0    15.999         -0.043407801307192         -0.055556028344352          0.000000000000000
    #   1 H     1.0000    0     1.008          1.780497256508764         -0.017018089151928          0.000000000000000
    #   2 H     1.0000    0     1.008         -0.462170608038134          1.719154625261312          0.000000000000000
    #

    geom_start = find_start(lines, 'CARTESIAN COORDINATES (A.U.)')
    if geom_start == -1:
        return Molecule()

    # Geometry starts 3 lines after header and ends with a blank line
    geom_start += 3
    end = ''
    mol_str = ''
    for i, line in enumerate(lines[geom_start:], start=geom_start):
        if line == end:
            break
        num, atom, z, frag, mass, x, y, z = line.split()
        mol_str += '{} {} {} {}\n'.format(atom, x, y, z)

    return Molecule(mol_str)


def harvest_scf_from_outfile(lines, psivar):
    """Harvest SCF results from the SCF section of the output file"""
    """Sample SCF results block"""
    #----------------
    #TOTAL SCF ENERGY
    #----------------
    #
    #Total Energy       :          -76.02602169 Eh           -2068.77322 eV
    #
    #Components:
    #Nuclear Repulsion  :            9.12509697 Eh             248.30651 eV
    #Electronic Energy  :          -85.15111867 Eh           -2317.07974 eV
    #
    #One Electron Energy:         -123.01434123 Eh           -3347.39040 eV
    #Two Electron Energy:           37.86322256 Eh            1030.31067 eV
    #
    #Virial components:
    #Potential Energy   :         -151.99262033 Eh           -4135.92947 eV
    #Kinetic Energy     :           75.96659864 Eh            2067.15624 eV
    #Virial Ratio       :            2.00078223
    #
    #

    scf_start = find_start(lines, 'TOTAL SCF ENERGY')
    if scf_start == -1:
        return ''

    # Energies in SCF block
    psivar['SCF TOTAL ENERGY'] = float(lines[scf_start + 3].split()[3])
    psivar['NUCLEAR REPULSION ENERGY'] = float(lines[scf_start + 6].split()[3])


def harvest_dipole(lines, psivar):
    """Harvest the dipole, and return as a tuple (x, y, z)
    Multiple different dipole moments can be output if post-HF calculations are
    run and their dipoles are requested resulting in highly similar blocks.
    It by default collects the last which appears to always be the one requested

    TODO: collect all the different types of dipole moments
    """
    """Sample dipole moment results block"""
    #-------------
    #DIPOLE MOMENT
    #-------------
    #                                X             Y             Z
    #Electronic contribution:     -0.11359      -0.14669      -0.00000
    #Nuclear contribution   :      0.61892       0.79867       0.00000
    #                        -----------------------------------------
    #Total Dipole Moment    :      0.50533       0.65198      -0.00000
    #                        -----------------------------------------
    #Magnitude (a.u.)       :      0.82489
    #Magnitude (Debye)      :      2.09670
    #
    #

    dipole_start = find_start(lines, 'DIPOLE MOMENT')

    if dipole_start != -1:
        # Dipole x, y, z are the last items 6 lines down in the dipole block
        dipole_str_list = lines[dipole_start + 6].split()[-3:]
        # Convert the dipole to debye
        dipole = [float(i)*psi_dipmom_au2debye for i in dipole_str_list]
        psivar['CURRENT DIPOLE X'] = dipole[0]
        psivar['CURRENT DIPOLE Y'] = dipole[1]
        psivar['CURRENT DIPOLE Z'] = dipole[2]

        # Dipole magnitude is 8 line down in the dipole block
        magnitude = float(lines[dipole_start + 8][-1])


def harvest_mp2(lines, psivar):
    """Harvest the MP2 results"""
    """Sample MP2 energy line (works for both MP2 and RI-MP2)"""
    #---------------------------------------
    #MP2 TOTAL ENERGY:      -76.226803665 Eh
    #---------------------------------------

    """Sample MP2 correlation energy line (yes there is a space)"""
    #-----------------------------------------------
    # MP2 CORRELATION ENERGY   :     -0.125436532 Eh
    #-----------------------------------------------

    """Sample RI-MP2 Correlation energy line (yes there is a space)"""
    #-----------------------------------------------
    # RI-MP2 CORRELATION ENERGY:     -0.125496692 Eh
    #-----------------------------------------------

    for line in reversed(lines):
        if line[:16] == 'MP2 TOTAL ENERGY':
            psivar['MP2 TOTAL ENERGY'] = line.split()[-2]
            break
    for line in reversed(lines):
        if line[:23] == ' MP2 CORRELATION ENERGY' or\
                line[:26] == ' RI-MP2 CORRELATION ENERGY':
            psivar['MP2 CORRELATION ENERGY'] = line.split()[-2]
            break


def harvest_coupled_cluster(lines, psivar):
    """Harvest coupled cluster results
    WARNING: Canonical and DLPNO print the coupled cluster results differently
    """
    """Sample (canonical) CCSD results block"""
    #----------------------
    #COUPLED CLUSTER ENERGY
    #----------------------
    #
    #E(0)                                       ...    -76.063720080
    #E(CORR)                                    ...     -0.288938791
    #E(TOT)                                     ...    -76.352658871
    #Singles Norm <S|S>**1/2                    ...      0.021106262
    #T1 diagnostic                              ...      0.007462191
    #

    """Sample DLPNO coupled cluster block (CCSD)"""
    #----------------------
    #COUPLED CLUSTER ENERGY
    #----------------------
    #
    #E(0)                                       ...    -76.026019996
    #E(CORR)(strong-pairs)                      ...     -0.211953159
    #E(CORR)(weak-pairs)                        ...     -0.000007244
    #E(CORR)(corrected)                         ...     -0.211960403
    #E(TOT)                                     ...    -76.237980399
    #Singles Norm <S|S>**1/2                    ...      0.014443573
    #T1 diagnostic                              ...      0.005106574
    #

    """Sample CCSD(T) block (same for DLPNO and canonical)"""
    #
    #Triples Correction (T)                     ...     -0.001544381
    #Final correlation energy                   ...     -0.134770265
    #E(CCSD)                                    ...    -75.709548429
    #E(CCSD(T))                                 ...    -75.711092810
    #

    cc_start = find_start(lines, 'COUPLED CLUSTER ENERGY')
    if cc_start == -1:
        return

    #psivar["CC REFERENCE"] = float(lines[cc_start + 3].split()[-1])

    # CCSD energy block is less than 20 lines
    for i, line in enumerate(lines[cc_start:cc_start + 20], start=cc_start):
        if line[:6] == "E(TOT)":
            psivar["CCSD TOTAL ENERGY"] = line.split()[-1]
            psivar["CCSD CORRELATION ENERGY"] = lines[i-1].split()[-1]
            #psivar["SINGLES NORM"] = lines[i+1].split()[-1]
            #psivar["T1 DIAGNOSTIC"] = lines[i+2].split()[-1]
            break

    # CCSD(T) energy block
    for i, line in enumerate(lines[cc_start:], start=cc_start):
        if line[:22] == "Triples Correction (T)":
            #psivar["TRIPLES CORRELATION ENERGY"] = line.split()[-1]
            psivar["CCSD(T) CORRELATION ENERGY"] = lines[i+1].split()[-1]
            psivar["CCSD TOTAL ENERGY"] = lines[i+2].split()[-1]
            psivar["CCSD(T) TOTAL ENERGY"] = lines[i+3].split()[-1]
            break


def harvest_engrad(engrad):
    """Parse the engrad file for the gradient"""
    try:
        lines = open(engrad).readlines()
    except IOError:
        return []
    num_atoms = int(lines[3].strip())
    energy = lines[7].strip()
    grad = []
    for i in range(12, 13 + num_atoms*3, 3):
        grad.append(list(map(float, lines[i:i + 3])))
    return grad


def find_start(lines, start_str, reverse=True):
    """Find the start of a block, iterate backwards by default,
    Usually the last one is wanted
    If not found, return -1
    """
    start = -1
    # Iterate backwards until the last value is found
    if reverse:
        for i, line in reversed(list(enumerate(lines))):
            if start_str == line:
                return i
    else:
        for i, line in enumerate(lines):
            if start_str == line:
                return i
    return start
