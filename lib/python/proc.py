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

from __future__ import print_function
"""Module with functions that encode the sequence of PSI module
calls for each of the *name* values of the energy(), optimize(),
response(), and frequency() function.

"""
import shutil
import os
import subprocess
import re
import psi4
import p4const
import p4util
from p4regex import *
#from extend_Molecule import *
from molutil import *
from functional import *
# never import driver, wrappers, or aliases into this file

# ATTN NEW ADDITIONS!
# consult http://sirius.chem.vt.edu/psi4manual/master/proc_py.html

def run_lmp2(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an LMP2 theory calculation.

    """

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
        scf_helper(name, **kwargs)
    psi4.lmp2()


def run_dcft(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density cumulant functional theory calculation.

    """

    flag = False
    if psi4.get_option('SCF', 'REFERENCE') == 'RHF':
        optstash = p4util.OptionsState(
            ['SCF', 'REFERENCE'],
            ['DCFT', 'REFERENCE'])
        psi4.set_local_option('SCF', 'REFERENCE', 'UHF')
        psi4.set_local_option('DCFT', 'REFERENCE', 'UHF')
        flag = True

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
        scf_helper(name, **kwargs)

    psi4.dcft()

    if flag:
        optstash.restore()


def run_dcft_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    DCFT gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['GLOBALS', 'DERTYPE'])

    psi4.set_global_option('DERTYPE', 'FIRST')
    run_dcft(name, **kwargs)
    psi4.deriv()

    optstash.restore()


def run_omp2(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an orbital-optimized MP2 computation

    """
    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
        scf_helper(name, **kwargs)

    return psi4.occ()


def run_omp2_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    OMP2 gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['REFERENCE'],
        ['GLOBALS', 'DERTYPE'])

    psi4.set_global_option('DERTYPE', 'FIRST')
    run_omp2(name, **kwargs)
    psi4.deriv()

    optstash.restore()


def run_mp2(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a MP2 calculation.

    """
    optstash = p4util.OptionsState(
        ['OCC', 'ORB_OPT'])

    # If the scf type is DF/CD, then the AO integrals were never written to disk
    if psi4.get_option('SCF', 'SCF_TYPE') == 'DF' or psi4.get_option('SCF', 'SCF_TYPE') == 'CD':
            mints = psi4.MintsHelper()
            mints.integrals()

    psi4.set_local_option('OCC', 'ORB_OPT', 'FALSE')
    run_omp2(name, **kwargs)

    optstash.restore()


def run_oldmp2(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a MP2 calculation.

    """
    optstash = p4util.OptionsState(
        ['TRANSQT2', 'WFN'],
        ['CCSORT', 'WFN'],
        ['MP2', 'WFN'])

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
        scf_helper(name, **kwargs)

        # If the scf type is DF/CD, then the AO integrals were never written to disk
        if psi4.get_option('SCF', 'SCF_TYPE') == 'DF' or psi4.get_option('SCF', 'SCF_TYPE') == 'CD':
            mints = psi4.MintsHelper()
            mints.integrals()

    psi4.set_local_option('TRANSQT2', 'WFN', 'MP2')
    psi4.set_local_option('CCSORT', 'WFN', 'MP2')
    psi4.set_local_option('MP2', 'WFN', 'MP2')

    psi4.transqt2()
    psi4.ccsort()
    psi4.mp2()

    optstash.restore()


def run_mp2_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a MP2 gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['REFERENCE'],
        ['GLOBALS', 'DERTYPE'])

    psi4.set_global_option('DERTYPE', 'FIRST')
    psi4.set_local_option('OCC', 'ORB_OPT', 'FALSE')
    run_omp2(name, **kwargs)
    psi4.deriv()

    optstash.restore()


def run_scs_omp2(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a spin-component scaled OMP2 computation

    """
    lowername = name.lower()

    optstash = p4util.OptionsState(
        ['OCC', 'SCS_TYPE'],
        ['OCC', 'DO_SCS'])

    # what type of scs?
    if (lowername == 'scs-omp2'):
        psi4.set_local_option('OCC', 'SCS_TYPE', 'SCS')
    elif (lowername == 'scsn-omp2'):
        psi4.set_local_option('OCC', 'SCS_TYPE', 'SCSN')
    #elif (lowername == 'scs-mi-omp2'):
    #    psi4.set_local_option('OCC', 'SCS_TYPE', 'SCSMI')
    elif (lowername == 'scs-omp2-vdw'):
        psi4.set_local_option('OCC', 'SCS_TYPE', 'SCSVDW')

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
        scf_helper(name, **kwargs)

    psi4.set_local_option('OCC', 'DO_SCS', 'TRUE')
    psi4.occ()

    optstash.restore()


def run_sos_omp2(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a spin-opposite scaled OMP2 computation

    """
    lowername = name.lower()

    optstash = p4util.OptionsState(
        ['OCC', 'SOS_TYPE'],
        ['OCC', 'DO_SOS'])

    # what type of sos?
    if (lowername == 'sos-omp2'):
        psi4.set_local_option('OCC', 'SOS_TYPE', 'SOS')
    elif (lowername == 'sos-pi-omp2'):
        psi4.set_local_option('OCC', 'SOS_TYPE', 'SOSPI')

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
        scf_helper(name, **kwargs)

    psi4.set_local_option('OCC', 'DO_SOS', 'TRUE')
    psi4.occ()

    optstash.restore()


def run_omp3(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an orbital-optimized MP3 computation

    """
    optstash = p4util.OptionsState(
        ['OCC', 'WFN_TYPE'])

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
        scf_helper(name, **kwargs)

    psi4.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
    psi4.occ()

    optstash.restore()


def run_omp3_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    OMP3 gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['OCC', 'WFN_TYPE'],
        ['GLOBALS', 'DERTYPE'])

    psi4.set_global_option('DERTYPE', 'FIRST')
    psi4.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
    run_omp3(name, **kwargs)
    psi4.deriv()

    optstash.restore()


def run_mp3(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a MP3 calculation.

    """
    optstash = p4util.OptionsState(
        ['OCC', 'ORB_OPT'])

    psi4.set_local_option('OCC', 'ORB_OPT', 'FALSE')
    run_omp3(name, **kwargs)

    optstash.restore()


def run_mp3_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a MP3 gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['GLOBALS', 'DERTYPE'],
        ['OCC', 'WFN_TYPE'],
        ['OCC', 'ORB_OPT'])

    psi4.set_global_option('DERTYPE', 'FIRST')
    psi4.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
    psi4.set_local_option('OCC', 'ORB_OPT', 'FALSE')
    run_omp3(name, **kwargs)
    psi4.deriv()

    optstash.restore()


def run_scs_omp3(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a spin-component scaled OMP3 computation

    """
    lowername = name.lower()

    optstash = p4util.OptionsState(
        ['OCC', 'SCS_TYPE'],
        ['OCC', 'DO_SCS'],
        ['OCC', 'WFN_TYPE'])

    # what type of scs?
    if (lowername == 'scs-omp3'):
        psi4.set_local_option('OCC', 'SCS_TYPE', 'SCS')
    elif (lowername == 'scsn-omp3'):
        psi4.set_local_option('OCC', 'SCS_TYPE', 'SCSN')
    #elif (lowername == 'scs-mi-omp3'):
    #    psi4.set_local_option('OCC', 'SCS_TYPE', 'SCSMI')
    elif (lowername == 'scs-omp3-vdw'):
        psi4.set_local_option('OCC', 'SCS_TYPE', 'SCSVDW')

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
        scf_helper(name, **kwargs)

    psi4.set_local_option('OCC', 'DO_SCS', 'TRUE')
    psi4.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
    psi4.occ()

    optstash.restore()


def run_sos_omp3(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a spin-opposite scaled OMP3 computation

    """
    lowername = name.lower()

    optstash = p4util.OptionsState(
        ['OCC', 'SOS_TYPE'],
        ['OCC', 'DO_SOS'],
        ['OCC', 'WFN_TYPE'])

    # what type of sos?
    if (lowername == 'sos-omp3'):
        psi4.set_local_option('OCC', 'SOS_TYPE', 'SOS')
    elif (lowername == 'sos-pi-omp3'):
        psi4.set_local_option('OCC', 'SOS_TYPE', 'SOSPI')

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
        scf_helper(name, **kwargs)

    psi4.set_local_option('OCC', 'DO_SOS', 'TRUE')
    psi4.set_local_option('OCC', 'WFN_TYPE', 'OMP3')
    psi4.occ()

    optstash.restore()


def run_ocepa(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an orbital-optimized CEPA computation

    """
    optstash = p4util.OptionsState(
        ['OCC', 'WFN_TYPE'])

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
        scf_helper(name, **kwargs)

    psi4.set_local_option('OCC', 'WFN_TYPE', 'OCEPA')
    psi4.occ()

    optstash.restore()


def run_ocepa_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    OCEPA gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['GLOBALS', 'DERTYPE'])

    psi4.set_global_option('DERTYPE', 'FIRST')
    run_ocepa(name, **kwargs)
    psi4.deriv()

    optstash.restore()


def run_cepa0(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a CEPA (LCCD) computation

    """
    optstash = p4util.OptionsState(
        ['OCC', 'WFN_TYPE'],
        ['OCC', 'ORB_OPT'])

    psi4.set_local_option('OCC', 'WFN_TYPE', 'OCEPA')
    psi4.set_local_option('OCC', 'ORB_OPT', 'FALSE')
    run_ocepa(name, **kwargs)

    optstash.restore()


def run_cepa0_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a CEPA(0) gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['GLOBALS', 'DERTYPE'],
        ['OCC', 'WFN_TYPE'],
        ['OCC', 'ORB_OPT'])

    psi4.set_global_option('DERTYPE', 'FIRST')
    psi4.set_local_option('OCC', 'WFN_TYPE', 'OCEPA')
    psi4.set_local_option('OCC', 'ORB_OPT', 'FALSE')
    run_ocepa(name, **kwargs)
    psi4.deriv()

    optstash.restore()


def run_omp2_5(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an orbital-optimized MP2.5 computation

    """
    optstash = p4util.OptionsState(
        ['OCC', 'WFN_TYPE'])

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
        scf_helper(name, **kwargs)

    psi4.set_local_option('OCC', 'WFN_TYPE', 'OMP2.5')
    psi4.occ()

    optstash.restore()


def run_omp2_5_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    OMP2.5 gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['GLOBALS', 'DERTYPE'],
        ['OCC', 'WFN_TYPE'])

    psi4.set_global_option('DERTYPE', 'FIRST')
    psi4.set_local_option('OCC', 'WFN_TYPE', 'OMP2.5')
    run_omp2_5(name, **kwargs)
    psi4.deriv()

    optstash.restore()


def run_mp2_5(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a MP2.5 calculation.

    """
    optstash = p4util.OptionsState(
        ['OCC', 'ORB_OPT'])

    psi4.set_local_option('OCC', 'ORB_OPT', 'FALSE')
    run_omp2_5(name, **kwargs)

    optstash.restore()


def run_mp2_5_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a MP3 gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['GLOBALS', 'DERTYPE'],
        ['OCC', 'WFN_TYPE'],
        ['OCC', 'ORB_OPT'])

    psi4.set_global_option('DERTYPE', 'FIRST')
    psi4.set_local_option('OCC', 'WFN_TYPE', 'OMP2.5')
    psi4.set_local_option('OCC', 'ORB_OPT', 'FALSE')
    run_omp2_5(name, **kwargs)
    psi4.deriv()

    optstash.restore()


def run_scf(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a self-consistent-field theory (HF & DFT) calculation.

    """
    lowername = name.lower()

    optstash = p4util.OptionsState(
        ['SCF', 'DFT_FUNCTIONAL'],
        ['SCF', 'SCF_TYPE'],
        ['SCF', 'REFERENCE'])

    # Alter default algorithm
    if not psi4.has_option_changed('SCF', 'SCF_TYPE'):
        psi4.set_local_option('SCF', 'SCF_TYPE', 'DF')

    if lowername == 'hf':
        if psi4.get_option('SCF', 'REFERENCE') == 'RKS':
            psi4.set_local_option('SCF', 'REFERENCE', 'RHF')
        elif psi4.get_option('SCF', 'REFERENCE') == 'UKS':
            psi4.set_local_option('SCF', 'REFERENCE', 'UHF')
        else:
            pass
    elif lowername == 'rhf':
        psi4.set_local_option('SCF', 'REFERENCE', 'RHF')
    elif lowername == 'uhf':
        psi4.set_local_option('SCF', 'REFERENCE', 'UHF')
    elif lowername == 'rohf':
        psi4.set_local_option('SCF', 'REFERENCE', 'ROHF')
    elif lowername == 'rscf':
        if (len(psi4.get_option('SCF', 'DFT_FUNCTIONAL')) > 0) or psi4.get_option('SCF', 'DFT_CUSTOM_FUNCTIONAL') is not None:
            psi4.set_local_option('SCF', 'REFERENCE', 'RKS')
        else:
            psi4.set_local_option('SCF', 'REFERENCE', 'RHF')
    elif lowername == 'uscf':
        if (len(psi4.get_option('SCF', 'DFT_FUNCTIONAL')) > 0) or psi4.get_option('SCF', 'DFT_CUSTOM_FUNCTIONAL') is not None:
            psi4.set_local_option('SCF', 'REFERENCE', 'UKS')
        else:
            psi4.set_local_option('SCF', 'REFERENCE', 'UHF')
    elif lowername == 'roscf':
        if (len(psi4.get_option('SCF', 'DFT_FUNCTIONAL')) > 0) or psi4.get_option('SCF', 'DFT_CUSTOM_FUNCTIONAL') is not None:
            raise ValidationError('ROHF reference for DFT is not available.')
        else:
            psi4.set_local_option('SCF', 'REFERENCE', 'ROHF')

    scf_helper(name, **kwargs)

    optstash.restore()


def run_scf_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a SCF gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['DF_BASIS_SCF'],
        ['SCF', 'SCF_TYPE'])

    # Alter default algorithm
    if not psi4.has_option_changed('SCF', 'SCF_TYPE'):
        psi4.set_local_option('SCF', 'SCF_TYPE', 'DF')

    run_scf(name, **kwargs)

    if (psi4.get_option('SCF', 'SCF_TYPE') == 'DF'):

        # if the df_basis_scf basis is not set, pick a sensible one.
        if psi4.get_global_option('DF_BASIS_SCF') == '':
            jkbasis = p4util.corresponding_jkfit(psi4.get_global_option('BASIS'))
            if jkbasis:
                psi4.set_global_option('DF_BASIS_SCF', jkbasis)
                psi4.print_out('\n  No DF_BASIS_SCF auxiliary basis selected, defaulting to %s\n\n' % (jkbasis))
            else:
                raise ValidationError('Keyword DF_BASIS_SCF is required.')

    psi4.scfgrad()
    optstash.restore()


def run_libfock(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a calculation through libfock, namely RCPHF,
    RCIS, RTDHF, RTDA, and RTDDFT.

    """
    if (name.lower() == 'cphf'):
        psi4.set_global_option('MODULE', 'RCPHF')
    if (name.lower() == 'cis'):
        psi4.set_global_option('MODULE', 'RCIS')
    if (name.lower() == 'tdhf'):
        psi4.set_global_option('MODULE', 'RTDHF')
    if (name.lower() == 'cpks'):
        psi4.set_global_option('MODULE', 'RCPKS')
    if (name.lower() == 'tda'):
        psi4.set_global_option('MODULE', 'RTDA')
    if (name.lower() == 'tddft'):
        psi4.set_global_option('MODULE', 'RTDDFT')

    psi4.libfock()


def run_mcscf(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a multiconfigurational self-consistent-field calculation.

    """
    return psi4.mcscf()


def scf_helper(name, **kwargs):
    """Function serving as helper to SCF, choosing whether to cast
    up or just run SCF with a standard guess. This preserves
    previous SCF options set by other procedures (e.g., SAPT
    output file types for SCF).

    """
    optstash = p4util.OptionsState(
        ['PUREAM'],
        ['BASIS'],
        ['DF_BASIS_SCF'],
        ['SCF', 'SCF_TYPE'],
        ['SCF', 'GUESS'],
        ['SCF', 'DF_INTS_IO'],
        ['SCF', 'SCF_TYPE'] # Hack: scope gets changed internally with the Andy trick
    )

    # if the df_basis_scf basis is not set, pick a sensible one.
    if psi4.get_option('SCF', 'SCF_TYPE') == 'DF':
        if psi4.get_global_option('DF_BASIS_SCF') == '':
            jkbasis = p4util.corresponding_jkfit(psi4.get_global_option('BASIS'))
            if jkbasis:
                psi4.set_global_option('DF_BASIS_SCF', jkbasis)
                psi4.print_out('\n  No DF_BASIS_SCF auxiliary basis selected, defaulting to %s\n\n' % (jkbasis))
            else:
                raise ValidationError('Keyword DF_BASIS_SCF is required.')

    optstash2 = p4util.OptionsState(
        ['BASIS'],
        ['DF_BASIS_SCF'],
        ['SCF', 'SCF_TYPE'],
        ['SCF', 'DF_INTS_IO'])

    # sort out cast_up settings. no need to stash these since only read, never reset
    cast = False
    if psi4.has_option_changed('SCF', 'BASIS_GUESS'):
        cast = psi4.get_option('SCF', 'BASIS_GUESS')
        if yes.match(str(cast)):
            cast = True
        elif no.match(str(cast)):
            cast = False

        if psi4.get_option('SCF', 'SCF_TYPE') == 'DF':
            castdf = True
        else:
            castdf = False

        if psi4.has_option_changed('SCF', 'DF_BASIS_GUESS'):
            castdf = psi4.get_option('SCF', 'DF_BASIS_GUESS')
            if yes.match(str(castdf)):
                castdf = True
            elif no.match(str(castdf)):
                castdf = False

    # sort out broken_symmetry settings.
    if 'brokensymmetry' in kwargs:
        molecule = psi4.get_active_molecule()
        multp = molecule.multiplicity()
        if multp != 1:
            raise ValidationError('Broken symmetry is only for singlets.')
        if psi4.get_option('SCF','REFERENCE') != 'UHF' and psi4.get_option('SCF','REFERENCE') != 'UKS':
            raise ValidationError('You must specify "set reference uhf" to use broken symmetry.')
        do_broken = True
    else:
        do_broken = False

    precallback = None
    if 'precallback' in kwargs:
        precallback = kwargs.pop('precallback')

    postcallback = None
    if 'postcallback' in kwargs:
        postcallback = kwargs.pop('postcallback')

    # Hack to ensure cartesian or pure are used throughout
    # Note that can't query PUREAM option directly, as it only
    #   reflects user changes to value, so load basis and
    #   read effective PUREAM setting off of it
    psi4.set_global_option('BASIS', psi4.get_global_option('BASIS'))
    psi4.set_global_option('PUREAM', psi4.MintsHelper().basisset().has_puream())

    # broken set-up
    if do_broken:
        molecule.set_multiplicity(3)
        psi4.print_out('\n')
        p4util.banner('  Computing high-spin triplet guess  ')
        psi4.print_out('\n')

    # cast set-up
    if (cast):

        if yes.match(str(cast)):
            guessbasis = '3-21G'
        else:
            guessbasis = cast

        if (castdf):
            if yes.match(str(castdf)):
                guessbasisdf = p4util.corresponding_jkfit(guessbasis)
            else:
                guessbasisdf = castdf

        # Switch to the guess namespace
        namespace = psi4.IO.get_default_namespace()
        psi4.IO.set_default_namespace((namespace + '.guess'))

        # Setup initial SCF
        psi4.set_global_option('BASIS', guessbasis)
        if (castdf):
            psi4.set_local_option('SCF', 'SCF_TYPE', 'DF')
            psi4.set_local_option('SCF', 'DF_INTS_IO', 'none')
            psi4.set_global_option('DF_BASIS_SCF', guessbasisdf)

        # Print some info about the guess
        psi4.print_out('\n')
        p4util.banner('Guess SCF, %s Basis' % (guessbasis))
        psi4.print_out('\n')

    # the FIRST scf call
    if cast or do_broken:
        # Perform the guess scf
        psi4.scf()

    # broken clean-up
    if do_broken:
        molecule.set_multiplicity(1)
        psi4.set_local_option('SCF', 'GUESS', 'READ')
        psi4.print_out('\n')
        p4util.banner('  Computing broken symmetry solution from high-spin triplet guess  ')
        psi4.print_out('\n')

    # cast clean-up
    if (cast):

        # Move files to proper namespace
        psi4.IO.change_file_namespace(180, (namespace + '.guess'), namespace)
        psi4.IO.set_default_namespace(namespace)

        # Set to read and project, and reset bases to final ones
        optstash2.restore()
        psi4.set_local_option('SCF', 'GUESS', 'READ')

        # Print the banner for the standard operation
        psi4.print_out('\n')
        p4util.banner(name.upper())
        psi4.print_out('\n')


    # the SECOND scf call
    e_scf = psi4.scf(precallback, postcallback)

    optstash.restore()
    return e_scf


def run_mp2_select(name, **kwargs):
    """Function selecting the algorithm for a MP2 energy call
    and directing toward the OCC (conv MP2) or the DFMP2 modules.

    """
    if (psi4.get_option("DFMP2", "MP2_TYPE") == "CONV") or (psi4.get_option("OCC", "MP2_TYPE") == "CONV"):
        return run_mp2(name, **kwargs)
    else:
        return run_dfmp2(name, **kwargs)


def run_mp2_select_gradient(name, **kwargs):
    """Function selecting the algorithm for a MP2 gradient call
    and directing toward the OCC (conv MP2) or the DFMP2 modules.

    """
    if (psi4.get_option("DFMP2", "MP2_TYPE") == "CONV") or (psi4.get_option("OCC", "MP2_TYPE") == "CONV"):
        return run_mp2_gradient(name, **kwargs)
    else:
        return run_dfmp2_gradient(name, **kwargs)


def run_dfmp2_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a DFMP2 gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['DF_BASIS_SCF'],
        ['DF_BASIS_MP2'],
        ['SCF_TYPE'])

    # Alter default algorithm
    if not psi4.has_option_changed('SCF', 'SCF_TYPE'):
        #psi4.set_local_option('SCF', 'SCF_TYPE', 'DF')  # insufficient b/c SCF option read in DFMP2
        psi4.set_global_option('SCF_TYPE', 'DF')

    if not psi4.get_option('SCF', 'SCF_TYPE') == 'DF':
        raise ValidationError('DF-MP2 gradients need DF-SCF reference, for now.')

    if 'restart_file' in kwargs:
        restartfile = kwargs.pop('restart_file')
        # Rename the checkpoint file to be consistent with psi4's file system
        psioh = psi4.IOManager.shared_object()
        psio = psi4.IO.shared_object()
        filepath = psioh.get_file_path(32)
        namespace = psio.get_default_namespace()
        pid = str(os.getpid())
        prefix = 'psi'
        targetfile = filepath + prefix + '.' + pid + '.' + namespace + '.32'
        if(psi4.me() == 0):
            shutil.copy(restartfile, targetfile)
    else:
        # if the df_basis_scf basis is not set, pick a sensible one.
        if psi4.get_global_option('DF_BASIS_SCF') == '':
            jkbasis = p4util.corresponding_jkfit(psi4.get_global_option('BASIS'))
            if jkbasis:
                psi4.set_global_option('DF_BASIS_SCF', jkbasis)
                psi4.print_out('\n  No DF_BASIS_SCF auxiliary basis selected, defaulting to %s\n\n' % (jkbasis))
            else:
                raise ValidationError('Keyword DF_BASIS_SCF is required.')

        scf_helper(name, **kwargs)

    psi4.print_out('\n')
    p4util.banner('DFMP2')
    psi4.print_out('\n')

    # if the df_basis_mp2 basis is not set, pick a sensible one.
    if psi4.get_global_option('DF_BASIS_MP2') == '':
        ribasis = p4util.corresponding_rifit(psi4.get_global_option('BASIS'))
        if ribasis:
            psi4.set_global_option('DF_BASIS_MP2', ribasis)
            psi4.print_out('  No DF_BASIS_MP2 auxiliary basis selected, defaulting to %s\n' % (ribasis))
        else:
            raise ValidationError('Keyword DF_BASIS_MP2 is required.')

    psi4.dfmp2grad()
    e_dfmp2 = psi4.get_variable('MP2 TOTAL ENERGY')
    e_scs_dfmp2 = psi4.get_variable('SCS-MP2 TOTAL ENERGY')

    optstash.restore()

    if (name.upper() == 'SCS-MP2'):
        return e_scs_dfmp2
    elif (name.upper() == 'DF-MP2') or (name.upper() == 'DFMP2') or (name.upper() == 'MP2'):
        return e_dfmp2


def run_ccenergy(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a CCSD, CC2, and CC3 calculation.

    """
    lowername = name.lower()

    optstash = p4util.OptionsState(
        ['TRANSQT2', 'WFN'],
        ['CCSORT', 'WFN'],
        ['CCENERGY', 'WFN'])

    if (lowername == 'ccsd'):
        psi4.set_local_option('TRANSQT2', 'WFN', 'CCSD')
        psi4.set_local_option('CCSORT', 'WFN', 'CCSD')
        psi4.set_local_option('CCENERGY', 'WFN', 'CCSD')
    elif (lowername == 'ccsd(t)'):
        psi4.set_local_option('TRANSQT2', 'WFN', 'CCSD_T')
        psi4.set_local_option('CCSORT', 'WFN', 'CCSD_T')
        psi4.set_local_option('CCENERGY', 'WFN', 'CCSD_T')
    elif (lowername == 'ccsd(at)' or lowername == 'a-ccsd(t)'):
        psi4.set_local_option('TRANSQT2', 'WFN', 'CCSD_AT')
        psi4.set_local_option('CCSORT', 'WFN', 'CCSD_AT')
        psi4.set_local_option('CCENERGY', 'WFN', 'CCSD_AT')
        psi4.set_local_option('CCHBAR', 'WFN', 'CCSD_AT')
        psi4.set_local_option('CCLAMBDA', 'WFN', 'CCSD_AT')
    elif (lowername == 'cc2'):
        psi4.set_local_option('TRANSQT2', 'WFN', 'CC2')
        psi4.set_local_option('CCSORT', 'WFN', 'CC2')
        psi4.set_local_option('CCENERGY', 'WFN', 'CC2')
    elif (lowername == 'cc3'):
        psi4.set_local_option('TRANSQT2', 'WFN', 'CC3')
        psi4.set_local_option('CCSORT', 'WFN', 'CC3')
        psi4.set_local_option('CCENERGY', 'WFN', 'CC3')
    elif (lowername == 'eom-cc2'):
        psi4.set_local_option('TRANSQT2', 'WFN', 'EOM_CC2')
        psi4.set_local_option('CCSORT', 'WFN', 'EOM_CC2')
        psi4.set_local_option('CCENERGY', 'WFN', 'EOM_CC2')
    elif (lowername == 'eom-ccsd'):
        psi4.set_local_option('TRANSQT2', 'WFN', 'EOM_CCSD')
        psi4.set_local_option('CCSORT', 'WFN', 'EOM_CCSD')
        psi4.set_local_option('CCENERGY', 'WFN', 'EOM_CCSD')
    # Call a plain energy('ccenergy') and have full control over options, incl. wfn
    elif(lowername == 'ccenergy'):
        pass

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
        scf_helper(name, **kwargs)

        # If the scf type is DF/CD, then the AO integrals were never written to disk
        if psi4.get_option('SCF', 'SCF_TYPE') == 'DF' or psi4.get_option('SCF', 'SCF_TYPE') == 'CD':
            mints = psi4.MintsHelper()
            mints.integrals()

    psi4.transqt2()
    psi4.ccsort()
    psi4.ccenergy()

    if (lowername == 'ccsd(at)' or lowername == 'a-ccsd(t)'):
        psi4.cchbar()
        psi4.cclambda()

    optstash.restore()


def run_cc_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a CCSD and CCSD(T) gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['GLOBALS', 'DERTYPE'],
        ['CCLAMBDA', 'WFN'],
        ['CCDENSITY', 'WFN'])

    psi4.set_global_option('DERTYPE', 'FIRST')

    run_ccenergy(name, **kwargs)
    if (name.lower() == 'ccsd'):
        psi4.set_local_option('CCLAMBDA', 'WFN', 'CCSD')
        psi4.set_local_option('CCDENSITY', 'WFN', 'CCSD')
    elif (name.lower() == 'ccsd(t)'):
        psi4.set_local_option('CCLAMBDA', 'WFN', 'CCSD_T')
        psi4.set_local_option('CCDENSITY', 'WFN', 'CCSD_T')

        user_ref = psi4.get_option('CCENERGY', 'REFERENCE')
        if (user_ref != 'UHF'):
            raise ValidationError('Reference %s for CCSD(T) gradients is not available.' % user_ref)

    psi4.cchbar()
    psi4.cclambda()
    psi4.ccdensity()
    psi4.deriv()

    optstash.restore()


def run_bccd(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a Brueckner CCD calculation.

    """
    optstash = p4util.OptionsState(
        ['TRANSQT2', 'DELETE_TEI'],
        ['TRANSQT2', 'WFN'],
        ['CCSORT', 'WFN'],
        ['CCENERGY', 'WFN'])

    if (name.lower() == 'bccd'):
        psi4.set_local_option('TRANSQT2', 'WFN', 'BCCD')
        psi4.set_local_option('CCSORT', 'WFN', 'BCCD')
        psi4.set_local_option('CCENERGY', 'WFN', 'BCCD')

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
        scf_helper(name, **kwargs)

        # If the scf type is DF/CD, then the AO integrals were never written to disk
        if psi4.get_option('SCF', 'SCF_TYPE') == 'DF' or psi4.get_option('SCF', 'SCF_TYPE') == 'CD':
            mints = psi4.MintsHelper()
            mints.integrals()

    psi4.set_local_option('TRANSQT2', 'DELETE_TEI', 'false')

    while True:
        psi4.transqt2()
        psi4.ccsort()
        psi4.ccenergy()
        psi4.print_out('Brueckner convergence check: %d\n' % psi4.get_variable('BRUECKNER CONVERGED'))
        if (psi4.get_variable('BRUECKNER CONVERGED') == True):
            break

    optstash.restore()


def run_bccd_t(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a Brueckner CCD(T) calculation.

    """
    optstash = p4util.OptionsState(
        ['TRANSQT2', 'WFN'],
        ['CCSORT', 'WFN'],
        ['CCENERGY', 'WFN'],
        ['CCTRIPLES', 'WFN'])

    psi4.set_local_option('TRANSQT2', 'WFN', 'BCCD_T')
    psi4.set_local_option('CCSORT', 'WFN', 'BCCD_T')
    psi4.set_local_option('CCENERGY', 'WFN', 'BCCD_T')
    psi4.set_local_option('CCTRIPLES', 'WFN', 'BCCD_T')
    run_bccd(name, **kwargs)
    psi4.cctriples()

    optstash.restore()


def run_scf_property(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    SCF calculations. This is a simple alias to :py:func:`~proc.run_scf`
    since SCF properties all handled through oeprop.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'SCF_TYPE'])

    # Alter default algorithm
    if not psi4.has_option_changed('SCF', 'SCF_TYPE'):
        psi4.set_local_option('SCF', 'SCF_TYPE', 'DF')

    run_scf(name, **kwargs)

    optstash.restore()


def run_cc_property(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    all CC property calculations.

    """
    oneel_properties = ['dipole', 'quadrupole']
    twoel_properties = []
    response_properties = ['polarizability', 'rotation', 'roa']
    excited_properties = ['oscillator_strength', 'rotational_strength']

    one = []
    two = []
    response = []
    excited = []
    invalid = []

    if 'properties' in kwargs:
        properties = kwargs.pop('properties')
        properties = p4util.drop_duplicates(properties)

        for prop in properties:
            if prop in oneel_properties:
                one.append(prop)
            elif prop in twoel_properties:
                two.append(prop)
            elif prop in response_properties:
                response.append(prop)
            elif prop in excited_properties:
                excited.append(prop)
            else:
                invalid.append(prop)
    else:
        raise ValidationError("The \"properties\" keyword is required with the property() function.")

    n_one = len(one)
    n_two = len(two)
    n_response = len(response)
    n_excited = len(excited)
    n_invalid = len(invalid)

    if (n_invalid > 0):
        print("The following properties are not currently supported: %s" % invalid)

    if (n_excited > 0 and (name.lower() != 'eom-ccsd' and name.lower() != 'eom-cc2')):
        raise ValidationError("Excited state CC properties require EOM-CC2 or EOM-CCSD.")

    if ((name.lower() == 'eom-ccsd' or name.lower() == 'eom-cc2') and n_response > 0):
        raise ValidationError("Cannot (yet) compute response properties for excited states.")

    if (n_one > 0 or n_two > 0) and (n_response > 0):
        print("Computing both density- and response-based properties.")

    if (name.lower() == 'ccsd'):
        psi4.set_global_option('WFN', 'CCSD')
        run_ccenergy('ccsd', **kwargs)
        psi4.set_global_option('WFN', 'CCSD')
    elif (name.lower() == 'cc2'):
        psi4.set_global_option('WFN', 'CC2')
        run_ccenergy('cc2', **kwargs)
        psi4.set_global_option('WFN', 'CC2')
    elif (name.lower() == 'eom-ccsd'):
        psi4.set_global_option('WFN', 'EOM_CCSD')
        run_ccenergy('eom-ccsd', **kwargs)
        psi4.set_global_option('WFN', 'EOM_CCSD')
    elif (name.lower() == 'eom-cc2'):
        psi4.set_global_option('WFN', 'EOM_CC2')
        run_ccenergy('eom-cc2', **kwargs)
        psi4.set_global_option('WFN', 'EOM_CC2')

    # Need cchbar for everything
    psi4.cchbar()

    # Need ccdensity at this point only for density-based props
    if (n_one > 0 or n_two > 0):
        if (name.lower() == 'eom-ccsd'):
            psi4.set_global_option('WFN', 'EOM_CCSD')
            psi4.set_global_option('DERTYPE', 'NONE')
            psi4.set_global_option('ONEPDM', 'TRUE')
            psi4.cceom()
        elif (name.lower() == 'eom-cc2'):
            psi4.set_global_option('WFN', 'EOM_CC2')
            psi4.set_global_option('DERTYPE', 'NONE')
            psi4.set_global_option('ONEPDM', 'TRUE')
            psi4.cceom()
        psi4.set_global_option('DERTYPE', 'NONE')
        psi4.set_global_option('ONEPDM', 'TRUE')
        psi4.cclambda()
        psi4.ccdensity()

    # Need ccresponse only for response-type props
    if (n_response > 0):
        psi4.set_global_option('DERTYPE', 'RESPONSE')
        psi4.cclambda()
        for prop in response:
            psi4.set_global_option('PROPERTY', prop)
            psi4.ccresponse()

    # Excited-state transition properties
    if (n_excited > 0):
        if (name.lower() == 'eom-ccsd'):
            psi4.set_global_option('WFN', 'EOM_CCSD')
        elif (name.lower() == 'eom-cc2'):
            psi4.set_global_option('WFN', 'EOM_CC2')
        else:
            raise ValidationError("Unknown excited-state CC wave function.")
        psi4.set_global_option('DERTYPE', 'NONE')
        psi4.set_global_option('ONEPDM', 'TRUE')
        psi4.cceom()
        psi4.cclambda()
        psi4.ccdensity()

    psi4.set_global_option('WFN', 'SCF')
    psi4.revoke_global_option_changed('WFN')
    psi4.set_global_option('DERTYPE', 'NONE')
    psi4.revoke_global_option_changed('DERTYPE')


def run_dfmp2_property(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a DFMP2 property calculation.

    """
    optstash = p4util.OptionsState(
        ['DF_BASIS_SCF'],
        ['DF_BASIS_MP2'],
        ['SCF_TYPE'])

    psi4.set_global_option('ONEPDM', 'TRUE')
    psi4.set_global_option('OPDM_RELAX', 'TRUE')

    # Alter default algorithm
    if not psi4.has_option_changed('SCF', 'SCF_TYPE'):
        #psi4.set_local_option('SCF', 'SCF_TYPE', 'DF')  # insufficient b/c SCF option read in DFMP2
        psi4.set_global_option('SCF_TYPE', 'DF')

    if not psi4.get_option('SCF', 'SCF_TYPE') == 'DF':
        raise ValidationError('DF-MP2 properties need DF-SCF reference, for now.')

    if 'restart_file' in kwargs:
        restartfile = kwargs.pop('restart_file')
        # Rename the checkpoint file to be consistent with psi4's file system
        psioh = psi4.IOManager.shared_object()
        psio = psi4.IO.shared_object()
        filepath = psioh.get_file_path(32)
        namespace = psio.get_default_namespace()
        pid = str(os.getpid())
        prefix = 'psi'
        targetfile = filepath + prefix + '.' + pid + '.' + namespace + '.32'
        if(psi4.me() == 0):
            shutil.copy(restartfile, targetfile)
    else:
        # if the df_basis_scf basis is not set, pick a sensible one.
        if psi4.get_global_option('DF_BASIS_SCF') == '':
            jkbasis = p4util.corresponding_jkfit(psi4.get_global_option('BASIS'))
            if jkbasis:
                psi4.set_global_option('DF_BASIS_SCF', jkbasis)
                psi4.print_out('\nNo DF_BASIS_SCF auxiliary basis selected, defaulting to %s\n\n' % (jkbasis))
            else:
                raise ValidationError('Keyword DF_BASIS_SCF is required.')

        scf_helper(name, **kwargs)

    psi4.print_out('\n')
    p4util.banner('DFMP2')
    psi4.print_out('\n')

    # if the df_basis_mp2 basis is not set, pick a sensible one.
    if psi4.get_global_option('DF_BASIS_MP2') == '':
        ribasis = p4util.corresponding_rifit(psi4.get_global_option('BASIS'))
        if ribasis:
            psi4.set_global_option('DF_BASIS_MP2', ribasis)
            psi4.print_out('No DF_BASIS_MP2 auxiliary basis selected, defaulting to %s\n' % (ribasis))
        else:
            raise ValidationError('Keyword DF_BASIS_MP2 is required.')

    psi4.dfmp2grad()
    e_dfmp2 = psi4.get_variable('MP2 TOTAL ENERGY')
    e_scs_dfmp2 = psi4.get_variable('SCS-MP2 TOTAL ENERGY')

    optstash.restore()

    if (name.upper() == 'SCS-MP2'):
        return e_scs_dfmp2
    elif (name.upper() == 'DF-MP2') or (name.upper() == 'DFMP2') or (name.upper() == 'MP2'):
        return e_dfmp2


def run_eom_cc(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an EOM-CC calculation, namely EOM-CC2, EOM-CCSD, and EOM-CC3.

    """
    optstash = p4util.OptionsState(
        ['TRANSQT2', 'WFN'],
        ['CCSORT', 'WFN'],
        ['CCENERGY', 'WFN'],
        ['CCHBAR', 'WFN'],
        ['CCEOM', 'WFN'])

    if (name.lower() == 'eom-ccsd'):
        psi4.set_local_option('TRANSQT2', 'WFN', 'EOM_CCSD')
        psi4.set_local_option('CCSORT', 'WFN', 'EOM_CCSD')
        psi4.set_local_option('CCENERGY', 'WFN', 'EOM_CCSD')
        psi4.set_local_option('CCHBAR', 'WFN', 'EOM_CCSD')
        psi4.set_local_option('CCEOM', 'WFN', 'EOM_CCSD')
        run_ccenergy('ccsd', **kwargs)
    elif (name.lower() == 'eom-cc2'):

        user_ref = psi4.get_option('CCENERGY', 'REFERENCE')
        if (user_ref != 'RHF') and (user_ref != 'UHF'):
            raise ValidationError('Reference %s for EOM-CC2 is not available.' % user_ref)

        psi4.set_local_option('TRANSQT2', 'WFN', 'EOM_CC2')
        psi4.set_local_option('CCSORT', 'WFN', 'EOM_CC2')
        psi4.set_local_option('CCENERGY', 'WFN', 'EOM_CC2')
        psi4.set_local_option('CCHBAR', 'WFN', 'EOM_CC2')
        psi4.set_local_option('CCEOM', 'WFN', 'EOM_CC2')
        run_ccenergy('cc2', **kwargs)
    elif (name.lower() == 'eom-cc3'):
        psi4.set_local_option('TRANSQT2', 'WFN', 'EOM_CC3')
        psi4.set_local_option('CCSORT', 'WFN', 'EOM_CC3')
        psi4.set_local_option('CCENERGY', 'WFN', 'EOM_CC3')
        psi4.set_local_option('CCHBAR', 'WFN', 'EOM_CC3')
        psi4.set_local_option('CCEOM', 'WFN', 'EOM_CC3')
        run_ccenergy('cc3', **kwargs)

    psi4.cchbar()
    psi4.cceom()

    optstash.restore()


def run_eom_cc_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an EOM-CCSD gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['CCDENSITY', 'XI'],
        ['CCDENSITY', 'ZETA'],
        ['CCLAMBDA', 'ZETA'],
        ['DERTYPE'],
        ['CCDENSITY', 'WFN'],
        ['CCLAMBDA', 'WFN'])

    psi4.set_global_option('DERTYPE', 'FIRST')

    if (name.lower() == 'eom-ccsd'):
        psi4.set_local_option('CCLAMBDA', 'WFN', 'EOM_CCSD')
        psi4.set_local_option('CCDENSITY', 'WFN', 'EOM_CCSD')
        run_eom_cc(name, **kwargs)

    psi4.set_local_option('CCLAMBDA', 'ZETA', 'FALSE')
    psi4.set_local_option('CCDENSITY', 'ZETA', 'FALSE')
    psi4.set_local_option('CCDENSITY', 'XI', 'TRUE')
    psi4.cclambda()
    psi4.ccdensity()
    psi4.set_local_option('CCLAMBDA', 'ZETA', 'TRUE')
    psi4.set_local_option('CCDENSITY', 'ZETA', 'TRUE')
    psi4.set_local_option('CCDENSITY', 'XI', 'FALSE')
    psi4.cclambda()
    psi4.ccdensity()
    psi4.deriv()

    optstash.restore()


def run_adc(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    an algebraic diagrammatic construction calculation.

    .. caution:: Get rid of active molecule lines- should be handled in energy.

    """
    if (psi4.get_option('ADC', 'REFERENCE') != 'RHF'):
        raise ValidationError('ADC requires reference RHF')

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
        scf_helper(name, **kwargs)

    return psi4.adc()


def run_dft(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-functional-theory calculation.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'DFT_FUNCTIONAL'],
        ['SCF', 'REFERENCE'],
        ['SCF', 'SCF_TYPE'],
        ['DF_BASIS_MP2'],
        ['DFMP2', 'MP2_OS_SCALE'],
        ['DFMP2', 'MP2_SS_SCALE'])

    # Alter default algorithm
    if not psi4.has_option_changed('SCF', 'SCF_TYPE'):
        psi4.set_local_option('SCF', 'SCF_TYPE', 'DF')

    psi4.set_local_option('SCF', 'DFT_FUNCTIONAL', name)

    user_ref = psi4.get_option('SCF', 'REFERENCE')
    if (user_ref == 'RHF'):
        psi4.set_local_option('SCF', 'REFERENCE', 'RKS')
    elif (user_ref == 'UHF'):
        psi4.set_local_option('SCF', 'REFERENCE', 'UKS')
    elif (user_ref == 'ROHF'):
        raise ValidationError('ROHF reference for DFT is not available.')
    elif (user_ref == 'CUHF'):
        raise ValidationError('CUHF reference for DFT is not available.')

    run_scf(name, **kwargs)
    returnvalue = psi4.get_variable('CURRENT ENERGY')

    for ssuper in superfunctional_list():
        if ssuper.name().lower() == name.lower():
            dfun = ssuper

    if dfun.is_c_hybrid():

        # if the df_basis_mp2 basis is not set, pick a sensible one.
        if psi4.get_global_option('DF_BASIS_MP2') == '':
            ribasis = p4util.corresponding_rifit(psi4.get_global_option('BASIS'))
            if ribasis:
                psi4.set_global_option('DF_BASIS_MP2', ribasis)
                psi4.print_out('  No DF_BASIS_MP2 auxiliary basis selected, defaulting to %s\n' % (ribasis))
            else:
                raise ValidationError('Keyword DF_BASIS_MP2 is required.')

        if dfun.is_c_scs_hybrid():
            psi4.set_local_option('DFMP2', 'MP2_OS_SCALE', dfun.c_os_alpha())
            psi4.set_local_option('DFMP2', 'MP2_SS_SCALE', dfun.c_ss_alpha())
            psi4.dfmp2()
            vdh = dfun.c_alpha() * psi4.get_variable('SCS-MP2 CORRELATION ENERGY')

        else:
            psi4.dfmp2()
            vdh = dfun.c_alpha() * psi4.get_variable('MP2 CORRELATION ENERGY')

        psi4.set_variable('DOUBLE-HYBRID CORRECTION ENERGY', vdh)
        returnvalue += vdh
        psi4.set_variable('DFT TOTAL ENERGY', returnvalue)
        psi4.set_variable('CURRENT ENERGY', returnvalue)

    optstash.restore()


def run_dft_gradient(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-functional-theory gradient calculation.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'DFT_FUNCTIONAL'],
        ['SCF', 'REFERENCE'],
        ['SCF', 'SCF_TYPE'])

    # Alter default algorithm
    if not psi4.has_option_changed('SCF', 'SCF_TYPE'):
        psi4.set_local_option('SCF', 'SCF_TYPE', 'DF')

    psi4.set_local_option('SCF', 'DFT_FUNCTIONAL', name)

    user_ref = psi4.get_option('SCF', 'REFERENCE')
    if (user_ref == 'RHF'):
        psi4.set_local_option('SCF', 'REFERENCE', 'RKS')
    elif (user_ref == 'UHF'):
        psi4.set_local_option('SCF', 'REFERENCE', 'UKS')
    elif (user_ref == 'ROHF'):
        raise ValidationError('ROHF reference for DFT is not available.')
    elif (user_ref == 'CUHF'):
        raise ValidationError('CUHF reference for DFT is not available.')

    run_scf_gradient(name, **kwargs)

    optstash.restore()


def run_detci(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a configuration interaction calculation, namely FCI,
    CIn, MPn, and ZAPTn.

    """
    optstash = p4util.OptionsState(
        ['TRANSQT2', 'WFN'],
        ['DETCI', 'WFN'],
        ['DETCI', 'MAX_NUM_VECS'],
        ['DETCI', 'MPN_ORDER_SAVE'],
        ['DETCI', 'MPN'],
        ['DETCI', 'FCI'],
        ['DETCI', 'EX_LEVEL'])

    user_ref = psi4.get_option('DETCI', 'REFERENCE')
    if (user_ref != 'RHF') and (user_ref != 'ROHF'):
        raise ValidationError('Reference %s for DETCI is not available.' % user_ref)

    if (name.lower() == 'zapt'):
        psi4.set_local_option('TRANSQT2', 'WFN', 'ZAPTN')
        psi4.set_local_option('DETCI', 'WFN', 'ZAPTN')
        level = kwargs['level']
        maxnvect = int((level + 1) / 2) + (level + 1) % 2
        psi4.set_local_option('DETCI', 'MAX_NUM_VECS', maxnvect)
        if ((level + 1) % 2):
            psi4.set_local_option('DETCI', 'MPN_ORDER_SAVE', 2)
        else:
            psi4.set_local_option('DETCI', 'MPN_ORDER_SAVE', 1)
    elif (name.lower() == 'detci-mp') or (name.lower() == 'mp'):
        psi4.set_local_option('TRANSQT2', 'WFN', 'DETCI')
        psi4.set_local_option('DETCI', 'WFN', 'DETCI')
        psi4.set_local_option('DETCI', 'MPN', 'TRUE')

        level = kwargs['level']
        maxnvect = int((level + 1) / 2) + (level + 1) % 2
        psi4.set_local_option('DETCI', 'MAX_NUM_VECS', maxnvect)
        if ((level + 1) % 2):
            psi4.set_local_option('DETCI', 'MPN_ORDER_SAVE', 2)
        else:
            psi4.set_local_option('DETCI', 'MPN_ORDER_SAVE', 1)
    elif (name.lower() == 'fci'):
            psi4.set_local_option('TRANSQT2', 'WFN', 'DETCI')
            psi4.set_local_option('DETCI', 'WFN', 'DETCI')
            psi4.set_local_option('DETCI', 'FCI', 'TRUE')
    elif (name.lower() == 'cisd'):
            psi4.set_local_option('TRANSQT2', 'WFN', 'DETCI')
            psi4.set_local_option('DETCI', 'WFN', 'DETCI')
            psi4.set_local_option('DETCI', 'EX_LEVEL', 2)
    elif (name.lower() == 'cisdt'):
            psi4.set_local_option('TRANSQT2', 'WFN', 'DETCI')
            psi4.set_local_option('DETCI', 'WFN', 'DETCI')
            psi4.set_local_option('DETCI', 'EX_LEVEL', 3)
    elif (name.lower() == 'cisdtq'):
            psi4.set_local_option('TRANSQT2', 'WFN', 'DETCI')
            psi4.set_local_option('DETCI', 'WFN', 'DETCI')
            psi4.set_local_option('DETCI', 'EX_LEVEL', 4)
    elif (name.lower() == 'ci'):
        psi4.set_local_option('TRANSQT2', 'WFN', 'DETCI')
        psi4.set_local_option('DETCI', 'WFN', 'DETCI')
        level = kwargs['level']
        psi4.set_local_option('DETCI', 'EX_LEVEL', level)
    # Call a plain energy('detci') and have full control over options
    elif(name.lower() == 'detci'):
        pass

    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
        scf_helper(name, **kwargs)

        # If the scf type is DF/CD, then the AO integrals were never written to disk
        if psi4.get_option('SCF', 'SCF_TYPE') == 'DF' or psi4.get_option('SCF', 'SCF_TYPE') == 'CD':
            psi4.MintsHelper().integrals()

    psi4.transqt2()
    psi4.detci()

    optstash.restore()


def run_dfmp2(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a density-fitted MP2 calculation.

    """
    optstash = p4util.OptionsState(
        ['DF_BASIS_MP2'],
        ['SCF', 'SCF_TYPE'])

    # Alter default algorithm
    if not psi4.has_option_changed('SCF', 'SCF_TYPE'):
        psi4.set_local_option('SCF', 'SCF_TYPE', 'DF')

    if 'restart_file' in kwargs:
        restartfile = kwargs.pop('restart_file')
        # Rename the checkpoint file to be consistent with psi4's file system
        psioh = psi4.IOManager.shared_object()
        psio = psi4.IO.shared_object()
        filepath = psioh.get_file_path(32)
        namespace = psio.get_default_namespace()
        pid = str(os.getpid())
        prefix = 'psi'
        targetfile = filepath + prefix + '.' + pid + '.' + namespace + '.32'
        if(psi4.me() == 0):
            shutil.copy(restartfile, targetfile)
    else:
        scf_helper(name, **kwargs)

    psi4.print_out('\n')
    p4util.banner('DFMP2')
    psi4.print_out('\n')

    # if the df_basis_mp2 basis is not set, pick a sensible one.
    if psi4.get_global_option('DF_BASIS_MP2') == '':
        ribasis = p4util.corresponding_rifit(psi4.get_global_option('BASIS'))
        if ribasis:
            psi4.set_global_option('DF_BASIS_MP2', ribasis)
            psi4.print_out('  No DF_BASIS_MP2 auxiliary basis selected, defaulting to %s\n' % (ribasis))
        else:
            raise ValidationError('Keyword DF_BASIS_MP2 is required.')

    e_dfmp2 = psi4.dfmp2()
    e_scs_dfmp2 = psi4.get_variable('SCS-MP2 TOTAL ENERGY')

    optstash.restore()

    if (name.upper() == 'SCS-MP2'):
        return e_scs_dfmp2
    elif (name.upper() == 'DF-MP2') or (name.upper() == 'DFMP2') or (name.upper() == 'MP2'):
        return e_dfmp2


def run_psimrcc(name, **kwargs):
    """Function encoding sequence of PSI module calls for a PSIMRCC computation
     using a reference from the MCSCF module

    """
    run_mcscf(name, **kwargs)
    psi4.psimrcc()
    return psi4.get_variable("CURRENT ENERGY")


def run_psimrcc_scf(name, **kwargs):
    """Function encoding sequence of PSI module calls for a PSIMRCC computation
     using a reference from the SCF module

    """
    # Bypass routine scf if user did something special to get it to converge
    if not (('bypass_scf' in kwargs) and yes.match(str(kwargs['bypass_scf']))):
        scf_helper(name, **kwargs)

    psi4.psimrcc()
    return psi4.get_variable("CURRENT ENERGY")


def run_mp2c(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a coupled MP2 calculation.

    """
    optstash = p4util.OptionsState(
        ['DF_BASIS_MP2'])

    molecule = psi4.get_active_molecule()
    molecule.update_geometry()
    monomerA = molecule.extract_subsets(1, 2)
    monomerA.set_name('monomerA')
    monomerB = molecule.extract_subsets(2, 1)
    monomerB.set_name('monomerB')

    # if the df_basis_mp2 basis is not set, pick a sensible one.
    if psi4.get_global_option('DF_BASIS_MP2') == '':
        ribasis = p4util.corresponding_rifit(psi4.get_global_option('BASIS'))
        if ribasis:
            psi4.set_global_option('DF_BASIS_MP2', ribasis)
            psi4.print_out('  No DF_BASIS_MP2 auxiliary basis selected, defaulting to %s\n' % (ribasis))
        else:
            raise ValidationError('Keyword DF_BASIS_MP2 is required.')

    ri = psi4.get_option('SCF', 'SCF_TYPE')
    df_ints_io = psi4.get_option('SCF', 'DF_INTS_IO')
    # inquire if above at all applies to dfmp2

    psi4.IO.set_default_namespace('dimer')
    psi4.set_local_option('SCF', 'SAPT', '2-dimer')
    psi4.print_out('\n')
    p4util.banner('Dimer HF')
    psi4.print_out('\n')
    psi4.set_global_option('DF_INTS_IO', 'SAVE')
    e_dimer = scf_helper('RHF', **kwargs)
    psi4.print_out('\n')
    p4util.banner('Dimer DFMP2')
    psi4.print_out('\n')
    e_dimer_mp2 = psi4.dfmp2()
    psi4.set_global_option('DF_INTS_IO', 'LOAD')

    activate(monomerA)
    if (ri == 'DF'):
        psi4.IO.change_file_namespace(97, 'dimer', 'monomerA')
    psi4.IO.set_default_namespace('monomerA')
    psi4.set_local_option('SCF', 'SAPT', '2-monomer_A')
    psi4.print_out('\n')
    p4util.banner('Monomer A HF')
    psi4.print_out('\n')
    e_monomerA = scf_helper('RHF', **kwargs)
    psi4.print_out('\n')
    p4util.banner('Monomer A DFMP2')
    psi4.print_out('\n')
    e_monomerA_mp2 = psi4.dfmp2()

    activate(monomerB)
    if (ri == 'DF'):
        psi4.IO.change_file_namespace(97, 'monomerA', 'monomerB')
    psi4.IO.set_default_namespace('monomerB')
    psi4.set_local_option('SCF', 'SAPT', '2-monomer_B')
    psi4.print_out('\n')
    p4util.banner('Monomer B HF')
    psi4.print_out('\n')
    e_monomerB = scf_helper('RHF', **kwargs)
    psi4.print_out('\n')
    p4util.banner('Monomer B DFMP2')
    psi4.print_out('\n')
    e_monomerB_mp2 = psi4.dfmp2()
    psi4.set_global_option('DF_INTS_IO', df_ints_io)

    psi4.IO.change_file_namespace(p4const.PSIF_SAPT_MONOMERA, 'monomerA', 'dimer')
    psi4.IO.change_file_namespace(p4const.PSIF_SAPT_MONOMERB, 'monomerB', 'dimer')

    activate(molecule)
    psi4.IO.set_default_namespace('dimer')
    psi4.set_local_option('SAPT', 'E_CONVERGENCE', 10e-10)
    psi4.set_local_option('SAPT', 'D_CONVERGENCE', 10e-10)
    psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'MP2C')
    psi4.print_out('\n')
    p4util.banner('MP2C')
    psi4.print_out('\n')

    psi4.set_variable('MP2C DIMER MP2 ENERGY', e_dimer_mp2)
    psi4.set_variable('MP2C MONOMER A MP2 ENERGY', e_monomerA_mp2)
    psi4.set_variable('MP2C MONOMER B MP2 ENERGY', e_monomerB_mp2)

    e_sapt = psi4.sapt()

    optstash.restore()
    return e_sapt


def run_sapt(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a SAPT calculation of any level.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'SCF_TYPE'])

    # Alter default algorithm
    if not psi4.has_option_changed('SCF', 'SCF_TYPE'):
        psi4.set_local_option('SCF', 'SCF_TYPE', 'DF')

    molecule = psi4.get_active_molecule()
    user_pg = molecule.schoenflies_symbol()
    molecule.reset_point_group('c1')
    molecule.fix_orientation(True)
    molecule.fix_com(True) # This should always have been set, very dangerous bug here
    molecule.update_geometry()
    if user_pg != 'c1':
        psi4.print_out('  SAPT does not make use of molecular symmetry, further calculations in C1 point group.\n')

    if (psi4.get_option('SCF','REFERENCE') != 'RHF' ):
        raise ValidationError('SAPT requires requires \"reference rhf\".')

    nfrag = molecule.nfragments()
    if nfrag != 2:
        raise ValidationError('SAPT requires active molecule to have 2 fragments, not %s.' % (nfrag))

    sapt_basis = 'dimer'
    if 'sapt_basis' in kwargs:
        sapt_basis = kwargs.pop('sapt_basis')
    sapt_basis = sapt_basis.lower()

    if (sapt_basis == 'dimer'):
        #molecule.update_geometry()
        monomerA = molecule.extract_subsets(1, 2)
        monomerA.set_name('monomerA')
        monomerB = molecule.extract_subsets(2, 1)
        monomerB.set_name('monomerB')
    elif (sapt_basis == 'monomer'):
        #molecule.update_geometry()
        monomerA = molecule.extract_subsets(1)
        monomerA.set_name('monomerA')
        monomerB = molecule.extract_subsets(2)
        monomerB.set_name('monomerB')

    ri = psi4.get_option('SCF', 'SCF_TYPE')
    df_ints_io = psi4.get_option('SCF', 'DF_INTS_IO')
    # inquire if above at all applies to dfmp2

    psi4.IO.set_default_namespace('dimer')
    psi4.set_local_option('SCF', 'SAPT', '2-dimer')
    psi4.print_out('\n')
    p4util.banner('Dimer HF')
    psi4.print_out('\n')
    if (sapt_basis == 'dimer'):
        psi4.set_global_option('DF_INTS_IO', 'SAVE')
    e_dimer = scf_helper('RHF', **kwargs)
    if (sapt_basis == 'dimer'):
        psi4.set_global_option('DF_INTS_IO', 'LOAD')

    activate(monomerA)
    if (ri == 'DF' and sapt_basis == 'dimer'):
        psi4.IO.change_file_namespace(97, 'dimer', 'monomerA')
    psi4.IO.set_default_namespace('monomerA')
    psi4.set_local_option('SCF', 'SAPT', '2-monomer_A')
    psi4.print_out('\n')
    p4util.banner('Monomer A HF')
    psi4.print_out('\n')
    e_monomerA = scf_helper('RHF', **kwargs)

    activate(monomerB)
    if (ri == 'DF' and sapt_basis == 'dimer'):
        psi4.IO.change_file_namespace(97, 'monomerA', 'monomerB')
    psi4.IO.set_default_namespace('monomerB')
    psi4.set_local_option('SCF', 'SAPT', '2-monomer_B')
    psi4.print_out('\n')
    p4util.banner('Monomer B HF')
    psi4.print_out('\n')
    e_monomerB = scf_helper('RHF', **kwargs)
    psi4.set_global_option('DF_INTS_IO', df_ints_io)

    psi4.IO.change_file_namespace(p4const.PSIF_SAPT_MONOMERA, 'monomerA', 'dimer')
    psi4.IO.change_file_namespace(p4const.PSIF_SAPT_MONOMERB, 'monomerB', 'dimer')

    activate(molecule)
    psi4.IO.set_default_namespace('dimer')
    psi4.set_local_option('SAPT', 'E_CONVERGENCE', 10e-10)
    psi4.set_local_option('SAPT', 'D_CONVERGENCE', 10e-10)
    if (name.lower() == 'sapt0'):
        psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT0')
    elif (name.lower() == 'sapt2'):
        psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2')
    elif (name.lower() == 'sapt2+'):
        psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+')
        psi4.set_local_option('SAPT', 'DO_CCD_DISP', False)
    elif (name.lower() == 'sapt2+(3)'):
        psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        psi4.set_local_option('SAPT', 'DO_THIRD_ORDER', False)
        psi4.set_local_option('SAPT', 'DO_CCD_DISP', False)
    elif (name.lower() == 'sapt2+3'):
        psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        psi4.set_local_option('SAPT', 'DO_THIRD_ORDER', True)
        psi4.set_local_option('SAPT', 'DO_CCD_DISP', False)
    elif (name.lower() == 'sapt2+(ccd)'):
        psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+')
        psi4.set_local_option('SAPT', 'DO_CCD_DISP', True)
    elif (name.lower() == 'sapt2+(3)(ccd)'):
        psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        psi4.set_local_option('SAPT', 'DO_THIRD_ORDER', False)
        psi4.set_local_option('SAPT', 'DO_CCD_DISP', True)
    elif (name.lower() == 'sapt2+3(ccd)'):
        psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        psi4.set_local_option('SAPT', 'DO_THIRD_ORDER', True)
        psi4.set_local_option('SAPT', 'DO_CCD_DISP', True)

    # if the df_basis_sapt basis is not set, pick a sensible one.
    if psi4.get_global_option('DF_BASIS_SAPT') == '':
        ribasis = p4util.corresponding_rifit(psi4.get_global_option('BASIS'))
        if ribasis:
            psi4.set_global_option('DF_BASIS_SAPT', ribasis)
            psi4.print_out('  No DF_BASIS_SAPT auxiliary basis selected, defaulting to %s\n' % (ribasis))
        else:
            raise ValidationError('Keyword DF_BASIS_SAPT is required.')

    psi4.print_out('\n')
    p4util.banner(name.upper())
    psi4.print_out('\n')
    e_sapt = psi4.sapt()

    molecule.reset_point_group(user_pg)
    molecule.update_geometry()

    optstash.restore()
    return e_sapt

def run_sapt_ct(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a charge-transfer SAPT calcuation of any level.

    """
    optstash = p4util.OptionsState(
        ['SCF', 'SCF_TYPE'])

    # Alter default algorithm
    if not psi4.has_option_changed('SCF', 'SCF_TYPE'):
        psi4.set_local_option('SCF', 'SCF_TYPE', 'DF')

    molecule = psi4.get_active_molecule()
    user_pg = molecule.schoenflies_symbol()
    molecule.reset_point_group('c1')
    molecule.fix_orientation(True)
    molecule.update_geometry()
    if user_pg != 'c1':
        psi4.print_out('  SAPT does not make use of molecular symmetry, further calculations in C1 point group.\n')

    if (psi4.get_option('SCF','REFERENCE') != 'RHF' ):
        raise ValidationError('SAPT requires requires \"reference rhf\".')

    nfrag = molecule.nfragments()
    if nfrag != 2:
        raise ValidationError('SAPT requires active molecule to have 2 fragments, not %s.' % (nfrag))

    monomerA = molecule.extract_subsets(1, 2)
    monomerA.set_name('monomerA')
    monomerB = molecule.extract_subsets(2, 1)
    monomerB.set_name('monomerB')
    molecule.update_geometry()
    monomerAm = molecule.extract_subsets(1)
    monomerAm.set_name('monomerAm')
    monomerBm = molecule.extract_subsets(2)
    monomerBm.set_name('monomerBm')

    ri = psi4.get_option('SCF', 'SCF_TYPE')
    df_ints_io = psi4.get_option('SCF', 'DF_INTS_IO')
    # inquire if above at all applies to dfmp2

    psi4.IO.set_default_namespace('dimer')
    psi4.set_local_option('SCF', 'SAPT', '2-dimer')
    psi4.print_out('\n')
    p4util.banner('Dimer HF')
    psi4.print_out('\n')
    psi4.set_global_option('DF_INTS_IO', 'SAVE')
    e_dimer = scf_helper('RHF', **kwargs)
    psi4.set_global_option('DF_INTS_IO', 'LOAD')

    activate(monomerA)
    if (ri == 'DF'):
        psi4.IO.change_file_namespace(97, 'dimer', 'monomerA')
    psi4.IO.set_default_namespace('monomerA')
    psi4.set_local_option('SCF', 'SAPT', '2-monomer_A')
    psi4.print_out('\n')
    p4util.banner('Monomer A HF (Dimer Basis)')
    psi4.print_out('\n')
    e_monomerA = scf_helper('RHF', **kwargs)

    activate(monomerB)
    if (ri == 'DF'):
        psi4.IO.change_file_namespace(97, 'monomerA', 'monomerB')
    psi4.IO.set_default_namespace('monomerB')
    psi4.set_local_option('SCF', 'SAPT', '2-monomer_B')
    psi4.print_out('\n')
    p4util.banner('Monomer B HF (Dimer Basis)')
    psi4.print_out('\n')
    e_monomerB = scf_helper('RHF', **kwargs)
    psi4.set_global_option('DF_INTS_IO', df_ints_io)

    activate(monomerAm)
    psi4.IO.set_default_namespace('monomerAm')
    psi4.set_local_option('SCF', 'SAPT', '2-monomer_A')
    psi4.print_out('\n')
    p4util.banner('Monomer A HF (Monomer Basis)')
    psi4.print_out('\n')
    e_monomerA = scf_helper('RHF', **kwargs)

    activate(monomerBm)
    psi4.IO.set_default_namespace('monomerBm')
    psi4.set_local_option('SCF', 'SAPT', '2-monomer_B')
    psi4.print_out('\n')
    p4util.banner('Monomer B HF (Monomer Basis)')
    psi4.print_out('\n')
    e_monomerB = scf_helper('RHF', **kwargs)

    activate(molecule)
    psi4.IO.set_default_namespace('dimer')
    psi4.set_local_option('SAPT', 'E_CONVERGENCE', 10e-10)
    psi4.set_local_option('SAPT', 'D_CONVERGENCE', 10e-10)
    if (name.lower() == 'sapt0-ct'):
        psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT0')
    elif (name.lower() == 'sapt2-ct'):
        psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2')
    elif (name.lower() == 'sapt2+-ct'):
        psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+')
    elif (name.lower() == 'sapt2+(3)-ct'):
        psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        psi4.set_local_option('SAPT', 'DO_THIRD_ORDER', False)
    elif (name.lower() == 'sapt2+3-ct'):
        psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        psi4.set_local_option('SAPT', 'DO_THIRD_ORDER', True)
    elif (name.lower() == 'sapt2+(ccd)-ct'):
        psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+')
        psi4.set_local_option('SAPT', 'DO_CCD_DISP', True)
    elif (name.lower() == 'sapt2+(3)(ccd)-ct'):
        psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        psi4.set_local_option('SAPT', 'DO_THIRD_ORDER', False)
        psi4.set_local_option('SAPT', 'DO_CCD_DISP', True)
    elif (name.lower() == 'sapt2+3(ccd)-ct'):
        psi4.set_local_option('SAPT', 'SAPT_LEVEL', 'SAPT2+3')
        psi4.set_local_option('SAPT', 'DO_THIRD_ORDER', True)
        psi4.set_local_option('SAPT', 'DO_CCD_DISP', True)
    psi4.print_out('\n')
    p4util.banner('SAPT Charge Transfer')
    psi4.print_out('\n')

    # if the df_basis_sapt basis is not set, pick a sensible one.
    if psi4.get_global_option('DF_BASIS_SAPT') == '':
        ribasis = p4util.corresponding_rifit(psi4.get_global_option('BASIS'))
        if ribasis:
            psi4.set_global_option('DF_BASIS_SAPT', ribasis)
            psi4.print_out('  No DF_BASIS_SAPT auxiliary basis selected, defaulting to %s\n' % (ribasis))
        else:
            raise ValidationError('Keyword DF_BASIS_SAPT is required.')

    psi4.print_out('\n')
    p4util.banner('Dimer Basis SAPT')
    psi4.print_out('\n')
    psi4.IO.change_file_namespace(p4const.PSIF_SAPT_MONOMERA, 'monomerA', 'dimer')
    psi4.IO.change_file_namespace(p4const.PSIF_SAPT_MONOMERB, 'monomerB', 'dimer')
    e_sapt = psi4.sapt()
    CTd = psi4.get_variable('SAPT CT ENERGY')

    psi4.print_out('\n')
    p4util.banner('Monomer Basis SAPT')
    psi4.print_out('\n')
    psi4.IO.change_file_namespace(p4const.PSIF_SAPT_MONOMERA, 'monomerAm', 'dimer')
    psi4.IO.change_file_namespace(p4const.PSIF_SAPT_MONOMERB, 'monomerBm', 'dimer')
    e_sapt = psi4.sapt()
    CTm = psi4.get_variable('SAPT CT ENERGY')
    CT = CTd - CTm

    psi4.print_out('\n\n')
    psi4.print_out('    SAPT Charge Transfer Analysis\n')
    psi4.print_out('  -----------------------------------------------------------------------------\n')
    line1 = '    SAPT Induction (Dimer Basis)      %10.4lf mH    %10.4lf kcal mol^-1\n' % (CTd * 1000.0, CTd * p4const.psi_hartree2kcalmol)
    line2 = '    SAPT Induction (Monomer Basis)    %10.4lf mH    %10.4lf kcal mol^-1\n' % (CTm * 1000.0, CTm * p4const.psi_hartree2kcalmol)
    line3 = '    SAPT Charge Transfer              %10.4lf mH    %10.4lf kcal mol^-1\n\n' % (CT * 1000.0, CT * p4const.psi_hartree2kcalmol)
    psi4.print_out(line1)
    psi4.print_out(line2)
    psi4.print_out(line3)
    psi4.set_variable('SAPT CT ENERGY', CT)

    molecule.reset_point_group(user_pg)
    molecule.update_geometry()

    optstash.restore()
    return e_sapt


def run_mrcc(name, **kwargs):
    """Function that prepares environment and input files
    for a calculation calling Kallay's MRCC code.

    """
    # TODO: Check to see if we really need to run the SCF code.
    scf_helper(name, **kwargs)
    vscf = psi4.get_variable('SCF TOTAL ENERGY')

    # The parse_arbitrary_order method provides us the following information
    # We require that level be provided. level is a dictionary
    # of settings to be passed to psi4.mrcc
    if not('level' in kwargs):
        raise ValidationError('level parameter was not provided.')

    level = kwargs['level']

    # Fullname is the string we need to search for in iface
    fullname = level['fullname']

    # User can provide 'keep' to the method.
    # When provided, do not delete the MRCC scratch directory.
    keep = False
    if 'keep' in kwargs:
        keep = kwargs['keep']

    # Save current directory location
    current_directory = os.getcwd()

    # Need to move to the scratch directory, perferrably into a separate directory in that location
    psi_io = psi4.IOManager.shared_object()
    os.chdir(psi_io.get_default_path())

    # Make new directory specifically for mrcc
    mrcc_tmpdir = 'mrcc_' + str(os.getpid())
    if 'path' in kwargs:
        mrcc_tmpdir = kwargs['path']

    # Check to see if directory already exists, if not, create.
    if os.path.exists(mrcc_tmpdir) == False:
        os.mkdir(mrcc_tmpdir)

    # Move into the new directory
    os.chdir(mrcc_tmpdir)

    # Generate integrals and input file (dumps files to the current directory)
    psi4.mrcc_generate_input(level)

    # Load the fort.56 file
    # and dump a copy into the outfile
    psi4.print_out('\n===== Begin fort.56 input for MRCC ======\n')
    psi4.print_out(open('fort.56', 'r').read())
    psi4.print_out('===== End   fort.56 input for MRCC ======\n')

    # Close output file
    psi4.close_outfile()

    # Modify the environment:
    #    PGI Fortan prints warning to screen if STOP is used
    os.environ['NO_STOP_MESSAGE'] = '1'

    # Obtain user's OMP_NUM_THREADS so that we don't blow it away.
    omp_num_threads_found = 'OMP_NUM_THREADS' in os.environ
    if omp_num_threads_found == True:
        omp_num_threads_user = os.environ['OMP_NUM_THREADS']

    # If the user provided MRCC_OMP_NUM_THREADS set the environ to it
    if psi4.has_option_changed('MRCC', 'MRCC_OMP_NUM_THREADS') == True:
        os.environ['OMP_NUM_THREADS'] = str(psi4.get_option('MRCC', 'MRCC_OMP_NUM_THREADS'))

    # Call dmrcc, directing all screen output to the output file
    try:
        if psi4.outfile_name() == 'stdout':
            retcode = subprocess.call('dmrcc', shell=True)
        else:
            retcode = subprocess.call('dmrcc >> ' + current_directory + '/' + psi4.outfile_name(), shell=True)

        if retcode < 0:
            print('MRCC was terminated by signal %d' % -retcode, file=sys.stderr)
            exit(1)
        elif retcode > 0:
            print('MRCC errored %d' % retcode, file=sys.stderr)
            exit(1)

    except OSError as e:
        print('Execution failed: %s' % e, file=sys.stderr)
        exit(1)

    # Restore the OMP_NUM_THREADS that the user set.
    if omp_num_threads_found == True:
        if psi4.has_option_changed('MRCC', 'MRCC_OMP_NUM_THREADS') == True:
            os.environ['OMP_NUM_THREADS'] = omp_num_threads_user

    # Scan iface file and grab the file energy.
    e = 0.0
    for line in file('iface'):
        fields = line.split()
        m = fields[1]
        try:
            e = float(fields[5])
            if m == "MP(2)":
                m = "MP2"
            psi4.set_variable(m + ' TOTAL ENERGY', e)
            psi4.set_variable(m + ' CORRELATION ENERGY', e - vscf)
        except ValueError:
            continue

    # The last 'e' in iface is the one the user requested.
    psi4.set_variable('CURRENT ENERGY', e)
    psi4.set_variable('CURRENT CORRELATION ENERGY', e - vscf)

    # Load the iface file
    iface = open('iface', 'r')
    iface_contents = iface.read()

    # Delete mrcc tempdir
    os.chdir('..')
    try:
        # Delete unless we're told not to
        if (keep == False and not('path' in kwargs)):
            shutil.rmtree(mrcc_tmpdir)
    except OSError as e:
        print('Unable to remove MRCC temporary directory %s' % e, file=sys.stderr)
        exit(1)

    # Revert to previous current directory location
    os.chdir(current_directory)

    # Reopen output file
    psi4.reopen_outfile()

    # If we're told to keep the files or the user provided a path, do nothing.
    if (keep != False or ('path' in kwargs)):
        psi4.print_out('\nMRCC scratch files have been kept.\n')
        psi4.print_out('They can be found in ' + mrcc_tmpdir)

    # Dump iface contents to output
    psi4.print_out('\n')
    p4util.banner('Full results from MRCC')
    psi4.print_out('\n')
    psi4.print_out(iface_contents)

    return e


def run_fnodfcc(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a DF-CCSD(T) computation.

    >>> energy('df-ccsd(t)')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # stash user options
    optstash = p4util.OptionsState(
        ['FNOCC','COMPUTE_TRIPLES'],
        ['FNOCC','DFCC'],
        ['FNOCC','NAT_ORBS'],
        ['FNOCC','RUN_CEPA'],
        ['SCF','DF_BASIS_SCF'],
        ['SCF','DF_INTS_IO'],
        ['SCF','SCF_TYPE'])

    psi4.set_local_option('FNOCC', 'DFCC', True)
    psi4.set_local_option('FNOCC', 'RUN_CEPA', False)

    # throw an exception for open-shells
    if (psi4.get_option('SCF','REFERENCE') != 'RHF' ):
        raise ValidationError("Error: %s requires \"reference rhf\"." % lowername)

    # override symmetry:
    molecule = psi4.get_active_molecule()
    user_pg = molecule.schoenflies_symbol()
    molecule.reset_point_group('c1')
    molecule.fix_orientation(1)
    molecule.update_geometry()
    if user_pg != 'c1':
        psi4.print_out('  FNOCC does not make use of molecular symmetry, further calculations in C1 point group.\n')

    # hack to ensure puream (or not) throughout
    psi4.set_global_option('PUREAM', psi4.MintsHelper().basisset().has_puream())

    # triples?
    if (lowername == 'df-ccsd'):
        psi4.set_local_option('FNOCC', 'COMPUTE_TRIPLES', False)
    if (lowername == 'df-ccsd(t)'):
        psi4.set_local_option('FNOCC', 'COMPUTE_TRIPLES', True)
    if (lowername == 'fno-df-ccsd'):
        psi4.set_local_option('FNOCC', 'COMPUTE_TRIPLES', False)
        psi4.set_local_option('FNOCC', 'NAT_ORBS', True)
    if (lowername == 'fno-df-ccsd(t)'):
        psi4.set_local_option('FNOCC', 'COMPUTE_TRIPLES', True)
        psi4.set_local_option('FNOCC', 'NAT_ORBS', True)

    # set scf-type to df unless the user wants something else
    if psi4.has_option_changed('SCF','SCF_TYPE') == False:
       psi4.set_global_option('SCF_TYPE', 'DF')

    scf_type = psi4.get_option('SCF','SCF_TYPE')
    if ( scf_type != 'CD' and scf_type != 'DF' ):
        raise ValidationError("Invalid scf_type for DFCC.")

    # save DF or CD ints generated by SCF for use in CC
    psi4.set_local_option('SCF','DF_INTS_IO', 'SAVE')

    # the default auxiliary basis
    if psi4.get_option('FNOCC','DF_BASIS_CC') == '':
       basis   = psi4.get_global_option('BASIS')
       dfbasis = p4util.corresponding_rifit(basis)
       psi4.set_local_option('FNOCC','DF_BASIS_CC',dfbasis)

    # make sure this module knows what df basis was used in the scf
    if ( psi4.get_option('SCF','SCF_TYPE') == "DF" ):
        df_basis_scf = psi4.get_option('SCF','DF_BASIS_SCF')
        if df_basis_scf == '':
           basis        = psi4.get_global_option('BASIS')
           df_basis_scf = p4util.corresponding_jkfit(basis)

        psi4.set_global_option('DF_BASIS_SCF',df_basis_scf)

    scf_helper(name,**kwargs)

    psi4.fnocc()

    molecule.reset_point_group(user_pg)
    molecule.update_geometry()

    # restore options
    optstash.restore()

    return psi4.get_variable("CURRENT ENERGY")


def run_fnocc(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a QCISD(T), CCSD(T), MP2.5, MP3, and MP4 computation.

    >>> energy('fno-ccsd(t)')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)
    if 'level' in kwargs:
        level = kwargs['level']
    else:
        level = 0

    # stash user options:
    optstash = p4util.OptionsState(
        ['TRANSQT2','WFN'],
        ['FNOCC','RUN_MP2'],
        ['FNOCC','RUN_MP3'],
        ['FNOCC','RUN_MP4'],
        ['FNOCC','RUN_CCSD'],
        ['FNOCC','COMPUTE_TRIPLES'],
        ['FNOCC','COMPUTE_MP4_TRIPLES'],
        ['FNOCC','DFCC'],
        ['FNOCC','RUN_CEPA'],
        ['FNOCC','USE_DF_INTS'],
        ['FNOCC','NAT_ORBS'])

    psi4.set_local_option('FNOCC','DFCC', False)
    psi4.set_local_option('FNOCC','RUN_CEPA', False)
    psi4.set_local_option('FNOCC','USE_DF_INTS', False)

    # which method?
    if (lowername == '_ccsd'):
        psi4.set_local_option('FNOCC','COMPUTE_TRIPLES', False)
        psi4.set_local_option('FNOCC','RUN_CCSD', True)
    elif (lowername == '_ccsd(t)'):
        psi4.set_local_option('FNOCC','COMPUTE_TRIPLES', True)
        psi4.set_local_option('FNOCC','RUN_CCSD', True)
    elif (lowername == 'fno-ccsd'):
        psi4.set_local_option('FNOCC','COMPUTE_TRIPLES', False)
        psi4.set_local_option('FNOCC','RUN_CCSD', True)
        psi4.set_local_option('FNOCC','NAT_ORBS', True)
    elif (lowername == 'fno-ccsd(t)'):
        psi4.set_local_option('FNOCC','COMPUTE_TRIPLES', True)
        psi4.set_local_option('FNOCC','RUN_CCSD', True)
        psi4.set_local_option('FNOCC','NAT_ORBS', True)
    elif (lowername == 'qcisd'):
        psi4.set_local_option('FNOCC','COMPUTE_TRIPLES', False)
        psi4.set_local_option('FNOCC','RUN_CCSD', False)
    elif (lowername == 'qcisd(t)'):
        psi4.set_local_option('FNOCC','COMPUTE_TRIPLES', True)
        psi4.set_local_option('FNOCC','RUN_CCSD', False)
    elif (lowername == 'fno-qcisd'):
        psi4.set_local_option('FNOCC','COMPUTE_TRIPLES', False)
        psi4.set_local_option('FNOCC','RUN_CCSD', False)
        psi4.set_local_option('FNOCC','NAT_ORBS', True)
    elif (lowername == 'fno-qcisd(t)'):
        psi4.set_local_option('FNOCC','COMPUTE_TRIPLES', True)
        psi4.set_local_option('FNOCC','NAT_ORBS', True)
        psi4.set_local_option('FNOCC','RUN_CCSD', False)
    elif (lowername == '_mp2'):
        psi4.set_local_option('FNOCC','RUN_MP2', True)
    elif (lowername == 'fno-mp3'):
        psi4.set_local_option('FNOCC','RUN_MP3', True)
        psi4.set_local_option('FNOCC','NAT_ORBS', True)
    elif (lowername == 'fno-mp4'):
        psi4.set_local_option('FNOCC','RUN_MP4', True)
        psi4.set_local_option('FNOCC','COMPUTE_MP4_TRIPLES', True)
        psi4.set_local_option('FNOCC','COMPUTE_TRIPLES', True)
        psi4.set_local_option('FNOCC','NAT_ORBS', True)
    elif (lowername == 'mp4(sdq)'):
        psi4.set_local_option('FNOCC','RUN_MP4', True)
        psi4.set_local_option('FNOCC','COMPUTE_MP4_TRIPLES', False)
        psi4.set_local_option('FNOCC','COMPUTE_TRIPLES', False)
    elif (lowername == 'fno-mp4(sdq)'):
        psi4.set_local_option('FNOCC','RUN_MP4', True)
        psi4.set_local_option('FNOCC','COMPUTE_MP4_TRIPLES', False)
        psi4.set_local_option('FNOCC','COMPUTE_TRIPLES', False)
        psi4.set_local_option('FNOCC','NAT_ORBS', True)
    elif (lowername == 'fnocc-mp') and (level == 3):
        psi4.set_local_option('FNOCC','RUN_MP3', True)
    elif (lowername == 'fnocc-mp') and (level == 4):
        psi4.set_local_option('FNOCC','RUN_MP4', True)
        psi4.set_local_option('FNOCC','COMPUTE_MP4_TRIPLES', True)
        psi4.set_local_option('FNOCC','COMPUTE_TRIPLES', True)

    # throw an exception for open-shells
    if (psi4.get_option('SCF','REFERENCE') != 'RHF' ):
        raise ValidationError("Error: %s requires \"reference rhf\"." % lowername)

    # scf
    scf_helper(name,**kwargs)

    # if the scf type is df/cd, then the ao integrals were never written to disk.
    if psi4.get_option('SCF', 'SCF_TYPE') == 'DF' or psi4.get_option('SCF', 'SCF_TYPE') == 'CD':
        # do we generate 4-index eri's with 3-index ones, or do we want conventional eri's?
        if psi4.get_option('FNOCC', 'USE_DF_INTS') == False:
            mints = psi4.MintsHelper()
            mints.integrals()

    # if this is not cim or FNO-CC, run transqt2.  otherwise, libtrans will be used
    if psi4.get_option('FNOCC','NAT_ORBS') == False and psi4.get_option('FNOCC','RUN_MP2') == False:
        if psi4.get_option('FNOCC', 'USE_DF_INTS') == False:
            psi4.set_local_option('TRANSQT2', 'WFN', 'CCSD')
            psi4.transqt2()

    # run ccsd
    psi4.fnocc()

    # set current correlation energy and total energy.  only need to treat mpn here.
    if (lowername == 'fnocc-mp') and (level == 3):
        emp3     = psi4.get_variable("MP3 TOTAL ENERGY")
        cemp3    = psi4.get_variable("MP3 CORRELATION ENERGY")
        psi4.set_variable("CURRENT ENERGY",emp3)
        psi4.set_variable("CURRENT CORRELATION ENERGY",cemp3)
    elif ( lowername == 'fno-mp3' ):
        emp3     = psi4.get_variable("MP3 TOTAL ENERGY")
        cemp3    = psi4.get_variable("MP3 CORRELATION ENERGY")
        psi4.set_variable("CURRENT ENERGY",emp3)
        psi4.set_variable("CURRENT CORRELATION ENERGY",cemp3)
    elif ( lowername == 'mp4(sdq)'):
        emp4sdq  = psi4.get_variable("MP4(SDQ) TOTAL ENERGY")
        cemp4sdq = psi4.get_variable("MP4(SDQ) CORRELATION ENERGY")
        psi4.set_variable("CURRENT ENERGY",emp4sdq)
        psi4.set_variable("CURRENT CORRELATION ENERGY",cemp4sdq)
    elif ( lowername == 'fno-mp4(sdq)'):
        emp4sdq  = psi4.get_variable("MP4(SDQ) TOTAL ENERGY")
        cemp4sdq = psi4.get_variable("MP4(SDQ) CORRELATION ENERGY")
        psi4.set_variable("CURRENT ENERGY",emp4sdq)
        psi4.set_variable("CURRENT CORRELATION ENERGY",cemp4sdq)
    elif ( lowername == 'fno-mp4'):
        emp4     = psi4.get_variable("MP4 TOTAL ENERGY")
        cemp4    = psi4.get_variable("MP4 CORRELATION ENERGY")
        psi4.set_variable("CURRENT ENERGY",emp4)
        psi4.set_variable("CURRENT CORRELATION ENERGY",cemp4)
    elif (lowername == 'fnocc-mp') and (level == 4):
        emp4     = psi4.get_variable("MP4 TOTAL ENERGY")
        cemp4    = psi4.get_variable("MP4 CORRELATION ENERGY")
        psi4.set_variable("CURRENT ENERGY",emp4)
        psi4.set_variable("CURRENT CORRELATION ENERGY",cemp4)

    # restore options
    optstash.restore()

    return psi4.get_variable("CURRENT ENERGY")


def run_cepa(name, **kwargs):
    """Function encoding sequence of PSI module calls for
    a cepa-like calculation.

    >>> energy('cepa(1)')

    """
    lowername = name.lower()
    uppername = name.upper()
    kwargs = p4util.kwargs_lower(kwargs)

    # save user options
    optstash = p4util.OptionsState(
        ['TRANSQT2', 'WFN'],
        ['FNOCC', 'NAT_ORBS'],
        ['FNOCC', 'RUN_CEPA'],
        ['FNOCC', 'USE_DF_INTS'],
        ['FNOCC', 'CEPA_NO_SINGLES'])

    psi4.set_local_option('FNOCC','RUN_CEPA', True)
    psi4.set_local_option('FNOCC','USE_DF_INTS', False)

    # what type of cepa?
    cepa_level = uppername
    if (lowername == 'cepa(2)'):
        raise ValidationError("Error: %s not implemented\n" % lowername)
    if (lowername == 'dci'):
        cepa_level = 'CISD'
    if (lowername == 'sdci'):
        cepa_level = 'CISD'

    if (lowername == 'fno-cepa(0)'):
        cepa_level = 'CEPA(0)'
        psi4.set_local_option('FNOCC', 'NAT_ORBS', True)
    if (lowername == 'fno-cepa(1)'):
        psi4.set_local_option('FNOCC', 'NAT_ORBS', True)
        cepa_level = 'CEPA(1)'
    if (lowername == 'fno-cepa(3)'):
        psi4.set_local_option('FNOCC', 'NAT_ORBS', True)
        cepa_level = 'CEPA(3)'
    if (lowername == 'fno-acpf'):
        psi4.set_local_option('FNOCC', 'NAT_ORBS', True)
        cepa_level = 'ACPF'
    if (lowername == 'fno-aqcc'):
        psi4.set_local_option('FNOCC', 'NAT_ORBS', True)
        cepa_level = 'AQCC'
    if (lowername == 'fno-sdci'):
        psi4.set_local_option('FNOCC', 'NAT_ORBS', True)
        cepa_level = 'CISD'
    if (lowername == 'fno-dci'):
        psi4.set_local_option('FNOCC', 'NAT_ORBS', True)
        cepa_level = 'CISD'

    psi4.set_local_option('FNOCC', 'CEPA_LEVEL', cepa_level)

    # throw an exception for open-shells
    if (psi4.get_option('SCF','REFERENCE') != 'RHF' ):
        raise ValidationError("Error: %s requires \"reference rhf\"." % lowername)

    psi4.set_local_option('TRANSQT2', 'WFN', 'CCSD')
    scf_helper(name, **kwargs)

    # If the scf type is DF/CD, then the AO integrals were never written to disk
    if psi4.get_option('SCF', 'SCF_TYPE') == 'DF' or psi4.get_option('SCF', 'SCF_TYPE') == 'CD':
        if psi4.get_option('FNOCC', 'USE_DF_INTS') == False:
            mints = psi4.MintsHelper()
            mints.integrals()

    if psi4.get_option('FNOCC','NAT_ORBS') == False:
        if psi4.get_option('FNOCC', 'USE_DF_INTS') == False:
            psi4.set_local_option('TRANSQT2', 'WFN', 'CCSD')
            psi4.transqt2()

    # run cepa
    psi4.fnocc()

    # one-electron properties
    if psi4.get_option('FNOCC', 'DIPMOM'):
        if cepa_level == "CEPA(1)" or cepa_level == "CEPA(3)":
            psi4.print_out("\n")
            psi4.print_out("    Error: one-electron properties not implemented for %s\n" % lowername)
            psi4.print_out("\n")
        elif psi4.get_option('FNOCC','NAT_ORBS'):
            psi4.print_out("\n")
            psi4.print_out("    Error: one-electron properties not implemented for %s\n" % lowername)
            psi4.print_out("\n")
        else:
            p4util.oeprop('DIPOLE','QUADRUPOLE','MULLIKEN_CHARGES','NO_OCCUPATIONS',title = cepa_level)

    # restore options
    optstash.restore()

    return psi4.get_variable("CURRENT ENERGY")

