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

"""
Module with database functionality.

"""

__all__ = [
    "database",
    "db",
    "DB_RGT",
    "DB_RXN",
]

import collections
import math
import os
import re
import sys

from psi4 import core
from .constants import constants
from . import p4util
from .driver import *
# never import aliases into this file


#########################
##  Start of Database  ##
#########################

DB_RGT = {}
DB_RXN = {}


def database(name, db_name, **kwargs):
    r"""Function to access the molecule objects and reference energies of
    popular chemical databases.

    :aliases: db()

    :returns: (*float*) Mean absolute deviation of the database in kcal/mol

    :PSI variables:

    .. hlist::
       :columns: 1

       * :psivar:`db_name DATABASE MEAN SIGNED DEVIATION`
       * :psivar:`db_name DATABASE MEAN ABSOLUTE DEVIATION`
       * :psivar:`db_name DATABASE ROOT-MEAN-SQUARE DEVIATION`
       * Python dictionaries of results accessible as ``DB_RGT`` and ``DB_RXN``.

    .. note:: It is very easy to make a database from a collection of xyz files
        using the script :source:`psi4/share/psi4/scripts/ixyz2database.py`.
        See :ref:`sec:createDatabase` for details.

    .. caution:: Some features are not yet implemented. Buy a developer some coffee.

       - In sow/reap mode, use only global options (e.g., the local option set by ``set scf scf_type df`` will not be respected).

    .. note:: To access a database that is not embedded in a |PSIfour|
        distribution, add the path to the directory containing the database
        to the environment variable :envvar:`PYTHONPATH`.

    :type name: str
    :param name: ``'scf'`` || ``'sapt0'`` || ``'ccsd(t)'`` || etc.

        First argument, usually unlabeled. Indicates the computational method
        to be applied to the database. May be any valid argument to
        :py:func:`psi4.driver.energy`.

    :type db_name: str
    :param db_name: ``'BASIC'`` || ``'S22'`` || ``'HTBH'`` || etc.

        Second argument, usually unlabeled. Indicates the requested database
        name, matching (case insensitive) the name of a python file in
        ``psi4/share/databases`` or :envvar:`PYTHONPATH`.  Consult that
        directory for available databases and literature citations.

    :type func: :ref:`function <op_py_function>`
    :param func: |dl| ``energy`` |dr| || ``optimize`` || ``cbs``

        Indicates the type of calculation to be performed on each database
        member. The default performs a single-point ``energy('name')``, while
        ``optimize`` perfoms a geometry optimization on each reagent, and
        ``cbs`` performs a compound single-point energy. If a nested series
        of python functions is intended (see :ref:`sec:intercalls`), use
        keyword ``db_func`` instead of ``func``.

    :type mode: str
    :param mode: |dl| ``'continuous'`` |dr| || ``'sow'`` || ``'reap'``

        Indicates whether the calculations required to complete the
        database are to be run in one file (``'continuous'``) or are to be
        farmed out in an embarrassingly parallel fashion
        (``'sow'``/``'reap'``).  For the latter, run an initial job with
        ``'sow'`` and follow instructions in its output file.

    :type cp: :ref:`boolean <op_py_boolean>`
    :param cp: ``'on'`` || |dl| ``'off'`` |dr|

        Indicates whether counterpoise correction is employed in computing
        interaction energies. Use this option and NOT the ``bsse_type="cp"``
        function for BSSE correction in database().  Option available
        (See :ref:`sec:availableDatabases`) only for databases of bimolecular complexes.

    :type rlxd: :ref:`boolean <op_py_boolean>`
    :param rlxd: ``'on'`` || |dl| ``'off'`` |dr|

        Indicates whether correction for deformation energy is
        employed in computing interaction energies.  Option available
        (See :ref:`sec:availableDatabases`) only for databases of bimolecular complexes
        with non-frozen monomers, e.g., HBC6.

    :type symm: :ref:`boolean <op_py_boolean>`
    :param symm: |dl| ``'on'`` |dr| || ``'off'``

        Indicates whether the native symmetry of the database reagents is
        employed (``'on'``) or whether it is forced to :math:`C_1` symmetry
        (``'off'``). Some computational methods (e.g., SAPT) require no
        symmetry, and this will be set by database().

    :type zpe: :ref:`boolean <op_py_boolean>`
    :param zpe: ``'on'`` || |dl| ``'off'`` |dr|

        Indicates whether zero-point-energy corrections are appended to
        single-point energy values. Option valid only for certain
        thermochemical databases. Disabled until Hessians ready.

    :type benchmark: str
    :param benchmark: |dl| ``'default'`` |dr| || ``'S22A'`` || etc.

        Indicates whether a non-default set of reference energies, if
        available (See :ref:`sec:availableDatabases`), are employed for the
        calculation of error statistics.

    :type tabulate: List[str]
    :param tabulate: |dl| ``[]`` |dr| || ``['scf total energy', 'natom']`` || etc.

        Indicates whether to form tables of variables other than the
        primary requested energy.  Available for any PSI variable.

    :type subset: Union[str, List[str]]
    :param subset:

        Indicates a subset of the full database to run. This is a very
        flexible option and can be used in three distinct ways, outlined
        below. Note that two take a string and the last takes an array.
        See :ref:`sec:availableDatabases` for available values.

        * ``'small'`` || ``'large'`` || ``'equilibrium'``
            Calls predefined subsets of the requested database, either
            ``'small'``, a few of the smallest database members,
            ``'large'``, the largest of the database members, or
            ``'equilibrium'``, the equilibrium geometries for a database
            composed of dissociation curves.
        * ``'BzBz_S'`` || ``'FaOOFaON'`` || ``'ArNe'`` ||  ``'HB'`` || etc.
            For databases composed of dissociation curves, or otherwise
            divided into subsets, individual curves and subsets can be
            called by name. Consult the database python files for available
            molecular systems (case insensitive).
        * ``[1,2,5]`` || ``['1','2','5']`` || ``['BzMe-3.5', 'MeMe-5.0']`` || etc.
            Specify a list of database members to run. Consult the
            database python files for available molecular systems.  This
            is the only portion of database input that is case sensitive;
            choices for this keyword must match the database python file.

    :examples:

    >>> # [1] Two-stage SCF calculation on short, equilibrium, and long helium dimer
    >>> db('scf','RGC10',cast_up='sto-3g',subset=['HeHe-0.85','HeHe-1.0','HeHe-1.5'], tabulate=['scf total energy','natom'])

    >>> # [2] Counterpoise-corrected interaction energies for three complexes in S22
    >>> #     Error statistics computed wrt an old benchmark, S22A
    >>> database('mp2','S22',cp=1,subset=[16,17,8],benchmark='S22A')

    >>> # [3] SAPT0 on the neon dimer dissociation curve
    >>> db('sapt0',subset='NeNe',cp=0,symm=0,db_name='RGC10')

    >>> # [4] Optimize system 1 in database S22, producing tables of scf and mp2 energy
    >>> db('mp2','S22',db_func=optimize,subset=[1], tabulate=['mp2 total energy','current energy'])

    >>> # [5] CCSD on the smallest systems of HTBH, a hydrogen-transfer database
    >>> database('ccsd','HTBH',subset='small', tabulate=['ccsd total energy', 'mp2 total energy'])

    """
    lowername = name  #TODO
    kwargs = p4util.kwargs_lower(kwargs)

    # Wrap any positional arguments into kwargs (for intercalls among wrappers)
    if not('name' in kwargs) and name:
        kwargs['name'] = name #.lower()
    if not('db_name' in kwargs) and db_name:
        kwargs['db_name'] = db_name

    # Establish function to call
    func = kwargs.pop('db_func', kwargs.pop('func', energy))
    kwargs['db_func'] = func
    # Bounce to CP if bsse kwarg (someday)
    if kwargs.get('bsse_type', None) is not None:
        raise ValidationError("""Database: Cannot specify bsse_type for database. Use the cp keyword withing database instead.""")

    allowoptexceeded = kwargs.get('allowoptexceeded', False)
    optstash = p4util.OptionsState(
        ['WRITER_FILE_LABEL'],
        ['SCF', 'REFERENCE'])

    # Wrapper wholly defines molecule. discard any passed-in
    kwargs.pop('molecule', None)

    # Paths to search for database files: here + PSIPATH + library + PYTHONPATH
    db_paths = []
    db_paths.append(os.getcwd())
    db_paths.extend(os.environ.get('PSIPATH', '').split(os.path.pathsep))
    db_paths.append(os.path.join(core.get_datadir(), 'databases'))
    db_paths.append(os.path.dirname(__file__))
    db_paths = list(map(os.path.abspath, db_paths))
    sys.path[1:1] = db_paths
    # TODO this should be modernized a la interface_cfour

    # Define path and load module for requested database
    database = p4util.import_ignorecase(db_name)
    if database is None:
        core.print_out('\nPython module for database %s failed to load\n\n' % (db_name))
        core.print_out('\nSearch path that was tried:\n')
        core.print_out(", ".join(map(str, sys.path)))
        raise ValidationError("Python module loading problem for database " + str(db_name))
    else:
        dbse = database.dbse
        HRXN = database.HRXN
        ACTV = database.ACTV
        RXNM = database.RXNM
        BIND = database.BIND
        TAGL = database.TAGL
        GEOS = database.GEOS
        try:
            DATA = database.DATA
        except AttributeError:
            DATA = {}

    user_writer_file_label = core.get_global_option('WRITER_FILE_LABEL')
    user_reference = core.get_global_option('REFERENCE')

    # Configuration based upon e_name & db_name options
    #   Force non-supramolecular if needed
    if not hasattr(lowername, '__call__') and re.match(r'^.*sapt', lowername):
        try:
            database.ACTV_SA
        except AttributeError:
            raise ValidationError('Database %s not suitable for non-supramolecular calculation.' % (db_name))
        else:
            ACTV = database.ACTV_SA
    #   Force open-shell if needed
    openshell_override = 0
    if user_reference in ['RHF', 'RKS']:
        try:
            database.isOS
        except AttributeError:
            pass
        else:
            if p4util.yes.match(str(database.isOS)):
                openshell_override = 1
                core.print_out('\nSome reagents in database %s require an open-shell reference; will be reset to UHF/UKS as needed.\n' % (db_name))

    # Configuration based upon database keyword options
    #   Option symmetry- whether symmetry treated normally or turned off (currently req'd for dfmp2 & dft)
    db_symm = kwargs.get('symm', True)

    symmetry_override = 0
    if db_symm is False:
        symmetry_override = 1
    elif db_symm is True:
        pass
    else:
        raise ValidationError("""Symmetry mode '%s' not valid.""" % (db_symm))

    #   Option mode of operation- whether db run in one job or files farmed out
    db_mode = kwargs.pop('db_mode', kwargs.pop('mode', 'continuous')).lower()
    kwargs['db_mode'] = db_mode

    if db_mode == 'continuous':
        pass
    elif db_mode == 'sow':
        pass
    elif db_mode == 'reap':
        db_linkage = kwargs.get('linkage', None)
        if db_linkage is None:
            raise ValidationError("""Database execution mode 'reap' requires a linkage option.""")
    else:
        raise ValidationError("""Database execution mode '%s' not valid.""" % (db_mode))

    #   Option counterpoise- whether for interaction energy databases run in bsse-corrected or not
    db_cp = kwargs.get('cp', False)

    if db_cp is True:
        try:
            database.ACTV_CP
        except AttributeError:
            raise ValidationError("""Counterpoise correction mode 'yes' invalid for database %s.""" % (db_name))
        else:
            ACTV = database.ACTV_CP
    elif db_cp is False:
        pass
    else:
        raise ValidationError("""Counterpoise correction mode '%s' not valid.""" % (db_cp))

    #   Option relaxed- whether for non-frozen-monomer interaction energy databases include deformation correction or not?
    db_rlxd = kwargs.get('rlxd', False)

    if db_rlxd is True:
        if db_cp is True:
            try:
                database.ACTV_CPRLX
                database.RXNM_CPRLX
            except AttributeError:
                raise ValidationError('Deformation and counterpoise correction mode \'yes\' invalid for database %s.' % (db_name))
            else:
                ACTV = database.ACTV_CPRLX
                RXNM = database.RXNM_CPRLX
        elif db_cp is False:
            try:
                database.ACTV_RLX
            except AttributeError:
                raise ValidationError('Deformation correction mode \'yes\' invalid for database %s.' % (db_name))
            else:
                ACTV = database.ACTV_RLX
    elif db_rlxd is False:
    #elif no.match(str(db_rlxd)):
        pass
    else:
        raise ValidationError('Deformation correction mode \'%s\' not valid.' % (db_rlxd))

    #   Option zero-point-correction- whether for thermochem databases jobs are corrected by zpe
    db_zpe = kwargs.get('zpe', False)

    if db_zpe is True:
        raise ValidationError('Zero-point-correction mode \'yes\' not yet implemented.')
    elif db_zpe is False:
        pass
    else:
        raise ValidationError('Zero-point-correction \'mode\' %s not valid.' % (db_zpe))

    #   Option benchmark- whether error statistics computed wrt alternate reference energies
    db_benchmark = 'default'
    if 'benchmark' in kwargs:
        db_benchmark = kwargs['benchmark']

        if db_benchmark.lower() == 'default':
            pass
        else:
            BIND = p4util.getattr_ignorecase(database, 'BIND_' + db_benchmark)
            if BIND is None:
                raise ValidationError('Special benchmark \'%s\' not available for database %s.' % (db_benchmark, db_name))

    #   Option tabulate- whether tables of variables other than primary energy method are formed
    # TODO db(func=cbs,tabulate=[non-current-energy])  # broken
    db_tabulate = []
    if 'tabulate' in kwargs:
        db_tabulate = kwargs['tabulate']

    #   Option subset- whether all of the database or just a portion is run
    db_subset = HRXN
    if 'subset' in kwargs:
        db_subset = kwargs['subset']

    if isinstance(db_subset, (str, bytes)):
        if db_subset.lower() == 'small':
            try:
                database.HRXN_SM
            except AttributeError:
                raise ValidationError("""Special subset 'small' not available for database %s.""" % (db_name))
            else:
                HRXN = database.HRXN_SM
        elif db_subset.lower() == 'large':
            try:
                database.HRXN_LG
            except AttributeError:
                raise ValidationError("""Special subset 'large' not available for database %s.""" % (db_name))
            else:
                HRXN = database.HRXN_LG
        elif db_subset.lower() == 'equilibrium':
            try:
                database.HRXN_EQ
            except AttributeError:
                raise ValidationError("""Special subset 'equilibrium' not available for database %s.""" % (db_name))
            else:
                HRXN = database.HRXN_EQ
        else:
            HRXN = p4util.getattr_ignorecase(database, db_subset)
            if HRXN is None:
                HRXN = p4util.getattr_ignorecase(database, 'HRXN_' + db_subset)
                if HRXN is None:
                    raise ValidationError("""Special subset '%s' not available for database %s.""" % (db_subset, db_name))
    else:
        temp = []
        for rxn in db_subset:
            if rxn in HRXN:
                temp.append(rxn)
            else:
                raise ValidationError("""Subset element '%s' not a member of database %s.""" % (str(rxn), db_name))
        HRXN = temp

    temp = []
    for rxn in HRXN:
        temp.append(ACTV['%s-%s' % (dbse, rxn)])
    HSYS = p4util.drop_duplicates(sum(temp, []))

    # Sow all the necessary reagent computations
    core.print_out("\n\n")
    p4util.banner(("Database %s Computation" % (db_name)))
    core.print_out("\n")

    #   write index of calcs to output file
    instructions = """\n    The database single-job procedure has been selected through mode='continuous'.\n"""
    instructions += """    Calculations for the reagents will proceed in the order below and will be followed\n"""
    instructions += """    by summary results for the database.\n\n"""
    for rgt in HSYS:
        instructions += """                    %-s\n""" % (rgt)
    core.print_out(instructions)

    #   Loop through chemical systems
    ERGT = {}
    ERXN = {}
    VRGT = {}
    VRXN = {}
    for rgt in HSYS:
        VRGT[rgt] = {}

        core.print_out('\n')
        p4util.banner(' Database {} Computation: Reagent {} \n   {}'.format(db_name, rgt, TAGL[rgt]))
        core.print_out('\n')

        molecule = core.Molecule.from_dict(GEOS[rgt].to_dict())
        molecule.set_name(rgt)
        molecule.update_geometry()

        if symmetry_override:
            molecule.reset_point_group('c1')
            molecule.fix_orientation(True)
            molecule.fix_com(True)
            molecule.update_geometry()

        if (openshell_override) and (molecule.multiplicity() != 1):
            if user_reference == 'RHF':
                core.set_global_option('REFERENCE', 'UHF')
            elif user_reference == 'RKS':
                core.set_global_option('REFERENCE', 'UKS')

        core.set_global_option('WRITER_FILE_LABEL', user_writer_file_label + ('' if user_writer_file_label == '' else '-') + rgt)

        if allowoptexceeded:
            try:
                ERGT[rgt] = func(molecule=molecule, **kwargs)
            except ConvergenceError:
                core.print_out(f"Optimization exceeded cycles for {rgt}")
                ERGT[rgt] = 0.0
        else:
            ERGT[rgt] = func(molecule=molecule, **kwargs)
        core.print_variables()
        core.print_out("   Database Contributions Map:\n   {}\n".format('-' * 75))
        for rxn in HRXN:
            db_rxn = dbse + '-' + str(rxn)
            if rgt in ACTV[db_rxn]:
                core.print_out('   reagent {} contributes by {:.4f} to reaction {}\n'.format(rgt, RXNM[db_rxn][rgt], db_rxn))
        core.print_out('\n')
        for envv in db_tabulate:
            VRGT[rgt][envv.upper()] = core.variable(envv)
        core.set_global_option("REFERENCE", user_reference)
        core.clean()
        #core.opt_clean()
        core.clean_variables()

    # Reap all the necessary reaction computations
    core.print_out("\n")
    p4util.banner(("Database %s Results" % (db_name)))
    core.print_out("\n")

    maxactv = []
    for rxn in HRXN:
        maxactv.append(len(ACTV[dbse + '-' + str(rxn)]))
    maxrgt = max(maxactv)
    table_delimit = '-' * (62 + 20 * maxrgt)
    tables = ''

    #   find any reactions that are incomplete
    FAIL = collections.defaultdict(int)
    for rxn in HRXN:
        db_rxn = dbse + '-' + str(rxn)
        for i in range(len(ACTV[db_rxn])):
            if abs(ERGT[ACTV[db_rxn][i]]) < 1.0e-12:
                if not allowoptexceeded:
                    FAIL[rxn] = 1

    #   tabulate requested process::environment variables
    tables += """   For each VARIABLE requested by tabulate, a 'Reaction Value' will be formed from\n"""
    tables += """   'Reagent' values according to weightings 'Wt', as for the REQUESTED ENERGY below.\n"""
    tables += """   Depending on the nature of the variable, this may or may not make any physical sense.\n"""
    for rxn in HRXN:
        db_rxn = dbse + '-' + str(rxn)
        VRXN[db_rxn] = {}

    for envv in db_tabulate:
        envv = envv.upper()
        tables += """\n   ==> %s <==\n\n""" % (envv.title())
        tables += _tblhead(maxrgt, table_delimit, 2)

        for rxn in HRXN:
            db_rxn = dbse + '-' + str(rxn)

            if FAIL[rxn]:
                tables += """\n%23s   %8s %8s %8s %8s""" % (db_rxn, '', '****', '', '')
                for i in range(len(ACTV[db_rxn])):
                    tables += """ %16.8f %2.0f""" % (VRGT[ACTV[db_rxn][i]][envv], RXNM[db_rxn][ACTV[db_rxn][i]])

            else:
                VRXN[db_rxn][envv] = 0.0
                for i in range(len(ACTV[db_rxn])):
                    VRXN[db_rxn][envv] += VRGT[ACTV[db_rxn][i]][envv] * RXNM[db_rxn][ACTV[db_rxn][i]]

                tables += """\n%23s        %16.8f                  """ % (db_rxn, VRXN[db_rxn][envv])
                for i in range(len(ACTV[db_rxn])):
                    tables += """ %16.8f %2.0f""" % (VRGT[ACTV[db_rxn][i]][envv], RXNM[db_rxn][ACTV[db_rxn][i]])
        tables += """\n   %s\n""" % (table_delimit)

    #   tabulate primary requested energy variable with statistics
    count_rxn = 0
    minDerror = 100000.0
    maxDerror = 0.0
    MSDerror = 0.0
    MADerror = 0.0
    RMSDerror = 0.0

    tables += """\n   ==> %s <==\n\n""" % ('Requested Energy')
    tables += _tblhead(maxrgt, table_delimit, 1)
    for rxn in HRXN:
        db_rxn = dbse + '-' + str(rxn)

        if FAIL[rxn]:
            tables += """\n%23s   %8.4f %8s %10s %10s""" % (db_rxn, BIND[db_rxn], '****', '****', '****')
            for i in range(len(ACTV[db_rxn])):
                tables += """ %16.8f %2.0f""" % (ERGT[ACTV[db_rxn][i]], RXNM[db_rxn][ACTV[db_rxn][i]])

        else:
            ERXN[db_rxn] = 0.0
            for i in range(len(ACTV[db_rxn])):
                ERXN[db_rxn] += ERGT[ACTV[db_rxn][i]] * RXNM[db_rxn][ACTV[db_rxn][i]]
            error = constants.hartree2kcalmol * ERXN[db_rxn] - BIND[db_rxn]

            tables += """\n%23s   %8.4f %8.4f %10.4f %10.4f""" % (db_rxn, BIND[db_rxn], constants.hartree2kcalmol * ERXN[db_rxn],
                error, error * constants.cal2J)
            for i in range(len(ACTV[db_rxn])):
                tables += """ %16.8f %2.0f""" % (ERGT[ACTV[db_rxn][i]], RXNM[db_rxn][ACTV[db_rxn][i]])

            if abs(error) < abs(minDerror):
                minDerror = error
            if abs(error) > abs(maxDerror):
                maxDerror = error
            MSDerror += error
            MADerror += abs(error)
            RMSDerror += error * error
            count_rxn += 1
    tables += """\n   %s\n""" % (table_delimit)

    if count_rxn:

        MSDerror /= float(count_rxn)
        MADerror /= float(count_rxn)
        RMSDerror = math.sqrt(RMSDerror / float(count_rxn))

        tables += """%23s %19s %10.4f %10.4f\n""" % ('Minimal Dev', '', minDerror, minDerror * constants.cal2J)
        tables += """%23s %19s %10.4f %10.4f\n""" % ('Maximal Dev', '', maxDerror, maxDerror * constants.cal2J)
        tables += """%23s %19s %10.4f %10.4f\n""" % ('Mean Signed Dev', '', MSDerror, MSDerror * constants.cal2J)
        tables += """%23s %19s %10.4f %10.4f\n""" % ('Mean Absolute Dev', '', MADerror, MADerror * constants.cal2J)
        tables += """%23s %19s %10.4f %10.4f\n""" % ('RMS Dev', '', RMSDerror, RMSDerror * constants.cal2J)
        tables += """   %s\n""" % (table_delimit)

        core.set_variable('%s DATABASE MEAN SIGNED DEVIATION' % (db_name), MSDerror)
        core.set_variable('%s DATABASE MEAN ABSOLUTE DEVIATION' % (db_name), MADerror)
        core.set_variable('%s DATABASE ROOT-MEAN-SQUARE DEVIATION' % (db_name), RMSDerror)

        core.print_out(tables)
        finalenergy = MADerror

    else:
        finalenergy = 0.0

    optstash.restore()

    DB_RGT.clear()
    DB_RGT.update(VRGT)
    DB_RXN.clear()
    DB_RXN.update(VRXN)
    return finalenergy


def _tblhead(tbl_maxrgt, tbl_delimit, ttype):
    r"""Function that prints the header for the changable-width results tables in db().
    *tbl_maxrgt* is the number of reagent columns the table must plan for. *tbl_delimit*
    is a string of dashes of the correct length to set off the table. *ttype* is 1 for
    tables comparing the computed values to the reference or 2 for simple tabulation
    and sum of the computed values.

    """
    tbl_str = ''
    tbl_str += """   %s""" % (tbl_delimit)
    if ttype == 1:
        tbl_str += """\n%23s %19s %21s""" % ('Reaction', 'Reaction Energy', 'Reaction Error')
    elif ttype == 2:
        tbl_str += """\n%23s     %19s %17s""" % ('Reaction', 'Reaction Value', '')
    for i in range(tbl_maxrgt):
        tbl_str += """%20s""" % ('Reagent ' + str(i + 1))
    if ttype == 1:
        tbl_str += """\n%23s   %8s %8s %10s %10s""" % ('', 'Ref', 'Calc', '[kcal/mol]', '[kJ/mol]')
    elif ttype == 2:
        tbl_str += """\n%65s""" % ('')
    for i in range(tbl_maxrgt):
        if ttype == 1:
            tbl_str += """%20s""" % ('[Eh] Wt')
        elif ttype == 2:
            tbl_str += """%20s""" % ('Value Wt')
    tbl_str += """\n   %s""" % (tbl_delimit)
    return tbl_str

##  Aliases  ##
db = database

#######################
##  End of Database  ##
#######################
