#!/usr/bin/env python

#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
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

import glob
import os
import sys

sys.path.append(os.path.dirname(__file__) + '/../../../driver')
try:
    import qcdb
except ImportError:
    print("""Cannot load qcdb python module. Run this script in situ or append the psi4/share/psi4/scripts directory to $PYTHONPATH.""")
    exit(1)

"""
Utility: This script converts a set of geometry files in XYZ format into
    a python database file for psi4 and qcdb scripts.
Instructions: Detailed instructions may be found at
    http://sirius.chem.vt.edu/psi4manual/latest/quickadddatabase.html .
    In short, move all XYZ files intended for a database into a directory
    and run this script from that directory. Answer a few questions about
    the intended database.  Edit the resulting database.py file if necessary,
    then copy it into psi4/share/psi4/databases/ . Its contents can be accessed as
    normal through the db() wrapper with no further configuration or recompiling.
Created: Monday, December 21, 2009, LAB
Last Modified: Tuesday, September 10, 2013, LAB
"""

# instructions
print("""
 Welcome to ixyz2database.
    Just fill in the variables when prompted. 
    Hit ENTER to accept default.
    Strings should not be in quotes.
    Elements in arrays should be space-delimited.
""")

# query database name
print("""
 Name your database.
    Recommend 3-8 characters, all capitalized letters.
    e.g., MINE or BZ6
""")
user_obedient = False
while not user_obedient:
    dbse = input('    dbse = ').strip()
    if dbse.isalnum():
        user_obedient = True

# query file extension
print("""
 XYZ file extension.
    All files with this extension in the current directory will be processed
    Additionally, all files with extension p4m in the current dir will be processed as psi4 mol format
""")
fext = input('    fext = [xyz] ').strip()
if fext == "":
    fext = 'xyz'

# query xyz file comment line
print("""
 What should line two of the XYZ file be used for (needn't be specially formatted in all files)
    [cgmp]      Treat first item in line as system charge, second as multiplicity, rest as comment
    [comment]   Treat content as text for the comment line
    [trash]     Ignore content
""")
user_obedient = False
while not user_obedient:
    line2 = input('    line2 = [cgmp] ').strip().lower()
    if line2 == "":
        line2 = 'cgmp'
    if line2 == 'comment' or line2 == 'cgmp' or line2 == 'trash':
        user_obedient = True

# query closed shell
print("""
 Are open-shell or non-singlets are present among your systems (or subsystems in the case of dimers)?
""")
isOS = qcdb.query_yes_no('    isOS =', False)

# query database type
print("""
 What is the nature of the systems in your incipient database?
    [1]         I have a bunch of plain molecules (no need to act on any subsystems)
                that I want to be able to act upon in parallel.
    [2]         I have a bunch of molecules that I want to form into a database
                whose reference quantity corresponds to various combinations thereof.
    [3]         I have a bunch of dimers (only dimer, no monomer, files should be present)
                that I want to form into a database whose reference quantity is interaction energy.
    Your final database may of course resemble any combination of these choices.
    This is but a humble script to get you started.
""")
user_obedient = False
while not user_obedient:
    route = input('    route = ').strip().lower()
    if route.isdigit():
        route = int(route)
        if route == 1 or route == 2 or route == 3:
            user_obedient = True

# query number of reactions
if route == 2:
    print("""
 How many reactions (things that have a reference quantity, as opposed
    to reagents that have a geometry) are in the database?
""")
    user_obedient = False
    while not user_obedient:
        Nrxn = input('    Nrxn = ').strip().lower()
        if Nrxn.isdigit():
            Nrxn = int(Nrxn)
            user_obedient = True
else:
    Nrxn = 1    # TODO really need?

# initialize containers
spy = ""
gpy = ""

HRXN = range(1, Nrxn + 1)
BINDRXN = {}
TAGLRXN = {}
for rxn in HRXN:
    BINDRXN[rxn] = None  # "nan" ? TODO
    TAGLRXN[rxn] = 'Reaction %s' % (rxn)

# reagent geometry section
gpy += "\n# <<< Geometry Specification Strings >>>\n"
gpy += "GEOS = {}\n\n"
count = 0
HRGT = []
TAGLRGT = {}
BINDRGT = {}

print("\n%-25s %6s %6s %6s %6s %6s\t\t%s\n" % ("system", "CHGsyst", "MLPsyst", "Natom", "Nmol1", "Nmol2", "Fragmentation Pattern"))

for xyzfile in (glob.glob('*.' + fext) + glob.glob('*.p4m')):

    # ascertain system name and open file
    system = os.path.splitext(xyzfile)[0]
    HRGT.append(system)

    f = open(xyzfile, 'r')
    text = f.readlines()
    f.close()

    # use Molecule object to read geometry in xyz file
    mol = qcdb.Molecule.from_string(''.join(text), fix_com=True, fix_orientation=True)
    Nsyst = mol.natom()

    # alter second line
    if line2 == 'cgmp':
        pass
    elif line2 == 'comment':
        mol.set_molecular_charge(0)
        mol.fragment_charges[0] = 0
        mol.set_multiplicity(1)
        mol.fragment_multiplicities[0] = 1
        mol.tagline = text[1].strip()
    elif line2 == 'trash':
        mol.set_molecular_charge(0)
        mol.fragment_charges[0] = 0
        mol.set_multiplicity(1)
        mol.fragment_multiplicities[0] = 1
        mol.tagline = ""

    CHGsyst = mol.molecular_charge()
    MLPsyst = mol.multiplicity()
    TAGLRGT[system] = mol.tagline
    BINDRGT[system] = None  # "nan" ?  # TODO

    if route == 3 and mol.nfragments() == 1:

        frag_pattern, mol = mol.BFS(return_molecule=True)
        Nmol1 = mol.fragments[0][1] - mol.fragments[0][0] + 1
        Nmol2 = mol.fragments[1][1] - mol.fragments[1][0] + 1

        print("%-25s %6d %6d %6d %6d %6d\t\t%s" % (system, CHGsyst, MLPsyst, Nsyst, Nmol1, Nmol2, frag_pattern))
        gpy += "GEOS['%%s-%%s-%%s' %% (dbse, '%s', 'dimer')] = qcdb.Molecule(\"\"\"\n" % (str(system))

        if mol.nfragments() != 2:
            print("ERROR: 2 fragments not detected for system %s." % (system))
            print("       If you really have trimers or above, contact LAB to modify this script.\n")
            sys.exit()

    else:

        print("%-25s %6d %6d %6d %6d %6d" % (system, CHGsyst, MLPsyst, Nsyst, Nsyst, 0))
        gpy += "GEOS['%%s-%%s-%%s' %% (dbse, '%s', 'reagent')] = qcdb.Molecule(\"\"\"\n" % (str(system))

    gpy += mol.create_psi4_string_from_molecule()
    gpy += """\"\"\")\n\n"""

    count += 1

Nrgt = len(HRGT)
if Nrgt != count:
    print("ERROR: discrepancy in counting systems $Nrgt vs $count!\n")
    sys.exit()


# python database file
docstring = """\"\"\"
| Database of <description of members and reference energy type>.
| Geometries from <Reference>.
| Reference interaction energies from <Reference>.

"""
if route == 3:
    docstring += """
- **cp**  ``'off'`` <erase this comment and after unless on is a valid option> || ``'on'``

- **rlxd** ``'off'`` <erase this comment and after unless on is valid option> || ``'on'``

"""
docstring += """
- **benchmark**

  - ``'<benchmark_name>'`` <Reference>.
  - |dl| ``'<default_benchmark_name>'`` |dr| <Reference>.

- **subset**

  - ``'small'`` <members_description>
  - ``'large'`` <members_description>
  - ``'<subset>'`` <members_description>

\"\"\"
"""

spy += docstring
spy += 'import re\n'
spy += 'import qcdb\n'

spy += "\n# <<< %s Database Module >>>\n" % (dbse)
spy += "dbse = %s\n" % ("'" + dbse + "'")
if isOS:
    spy += "isOS = '%s'\n" % (isOS)

spy += "\n# <<< Database Members >>>\n"
spy += "HRXN = ["
if route == 1:
    for rgt in HRGT:
        spy += "'%s', " % (rgt)
elif route == 2:
    for rxn in HRXN:
        spy += "'%s', " % (rxn)
elif route == 3:
    for rgt in HRGT:
        spy += "'%s', " % (rgt)
spy += "]\n"
spy += "HRXN_SM = []\n"
spy += "HRXN_LG = []\n"

spy += "\n# <<< Chemical Systems Involved >>>\n"
spy += "RXNM = {}     # reaction matrix of reagent contributions per reaction\n"
spy += "ACTV = {}     # order of active reagents per reaction\n"

if route == 1:
    for rgt in HRGT:
        spy += """ACTV['%%s-%%s'            %% (dbse, %-23s )] = """ % ("'" + rgt + "'")
        spy += """['%%s-%%s-reagent'      %% (dbse, %s)]\n""" % ("'" + rgt + "'")
        spy += """RXNM['%%s-%%s'            %% (dbse, %-23s )] = """ % ("'" + rgt + "'")
        spy += """dict(zip(ACTV['%%s-%%s' %% (dbse, %s)], [+1]))\n\n""" % ("'" + rgt + "'")

elif route == 2:
    for rxn in HRXN:
        spy += """ACTV['%%s-%%s'            %% (dbse, %-23s )] = """ % ("'" + str(rxn) + "'")
        spy += """['%%s-%%s-reagent'      %% (dbse, %s),\n""" % ("''")
        spy += """%62s '%%s-%%s-reagent'      %% (dbse, %s),\n""" % ("", "''")
        spy += """%62s '%%s-%%s-reagent'      %% (dbse, %s) ]\n""" % ("", "''")
        spy += """RXNM['%%s-%%s'            %% (dbse, %-23s )] = """ % ("'" + str(rxn) + "'")
        spy += """dict(zip(ACTV['%%s-%%s' %% (dbse, %s)], []))\n\n""" % ("'" + str(rxn) + "'")

elif route == 3:
    pass

    spy += "ACTV_CP = {}  # order of active reagents per counterpoise-corrected reaction\n"
    spy += "ACTV_SA = {}  # order of active reagents for non-supermolecular calculations\n"

    spy += "for rxn in HRXN:\n\n"

    spy += "    RXNM[   '%s-%s' % (dbse, rxn)] = {'%s-%s-dimer'      % (dbse, rxn) : +1,\n"
    spy += "                                      '%s-%s-monoA-CP'   % (dbse, rxn) : -1,\n"
    spy += "                                      '%s-%s-monoB-CP'   % (dbse, rxn) : -1,\n"
    spy += "                                      '%s-%s-monoA-unCP' % (dbse, rxn) : -1,\n"
    spy += "                                      '%s-%s-monoB-unCP' % (dbse, rxn) : -1 }\n\n"

    spy += "    ACTV_SA['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn) ]\n\n"

    spy += "    ACTV_CP['%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),\n"
    spy += "                                      '%s-%s-monoA-CP'   % (dbse, rxn),\n"
    spy += "                                      '%s-%s-monoB-CP'   % (dbse, rxn) ]\n\n"

    spy += "    ACTV[   '%s-%s' % (dbse, rxn)] = ['%s-%s-dimer'      % (dbse, rxn),\n"
    spy += "                                      '%s-%s-monoA-unCP' % (dbse, rxn),\n"
    spy += "                                      '%s-%s-monoB-unCP' % (dbse, rxn) ]\n\n"


spy += "# <<< Reference Values [kcal/mol] >>>\n"
spy += "BIND = {}\n"
#print SPY_OUT "nan = float('NaN')\n";
if route == 1:
    for rgt in HRGT:
        spy += """BIND['%%s-%%s'            %% (dbse, %-23s )] = """ % ("'" + rgt + "'")
        spy += """%8.3f\n""" % (0.0)  # TODO BINDRGT[rgt]))
elif route == 2:
    for rxn in HRXN:
        pass
        spy += """BIND['%%s-%%s'            %% (dbse, %-23s )] = """ % ("'" + str(rxn) + "'")
        spy += """%8.3f\n""" % (0.0)  # TODO BINDRGT[rxn]))
elif route == 3:
    for rgt in HRGT:
        pass
        spy += """BIND['%%s-%%s'            %% (dbse, %-23s )] = """ % ("'" + rgt + "'")
        spy += """%8.3f\n""" % (0.0)  # TODO BINDRGT[rgt]))

# write comment line section
spy += "\n# <<< Comment Lines >>>\n"
spy += "TAGL = {}\n"
if route == 1:
    for rgt in HRGT:
        spy += """TAGL['%%s-%%s'            %% (dbse, %-23s )] = """ % ("'" + rgt + "'")
        spy += """\"\"\"%s \"\"\"\n""" % (TAGLRGT[rgt])
        spy += """TAGL['%%s-%%s-reagent'    %% (dbse, %-23s )] = """ % ("'" + rgt + "'")
        spy += """\"\"\"%s \"\"\"\n""" % (TAGLRGT[rgt])

elif route == 2:
    for rxn in HRXN:
        spy += """TAGL['%%s-%%s'            %% (dbse, %-23s )] = """ % ("'" + str(rxn) + "'")
        spy += """\"\"\"%s \"\"\"\n""" % (TAGLRXN[rxn])
    for rgt in HRGT:
        spy += """TAGL['%%s-%%s-reagent'    %% (dbse, %-23s )] = """ % ("'" + rgt + "'")
        spy += """\"\"\"%s \"\"\"\n""" % (TAGLRGT[rgt])

elif route == 3:

    for rgt in HRGT:
        spy += """TAGL['%%s-%%s'            %% (dbse, %-23s )] = """ % ("'" + rgt + "'")
        spy += """\"\"\"%s \"\"\"\n""" % (TAGLRGT[rgt])
        spy += """TAGL['%%s-%%s-dimer'      %% (dbse, %-23s )] = """ % ("'" + rgt + "'")
        spy += """\"\"\"%s \"\"\"\n""" % ('Dimer from ' + TAGLRGT[rgt])
        spy += """TAGL['%%s-%%s-monoA-CP'   %% (dbse, %-23s )] = """ % ("'" + rgt + "'")
        spy += """\"\"\"%s \"\"\"\n""" % ('Monomer A from ' + TAGLRGT[rgt])
        spy += """TAGL['%%s-%%s-monoB-CP'   %% (dbse, %-23s )] = """ % ("'" + rgt + "'")
        spy += """\"\"\"%s \"\"\"\n""" % ('Monomer B from ' + TAGLRGT[rgt])
        spy += """TAGL['%%s-%%s-monoA-unCP' %% (dbse, %-23s )] = """ % ("'" + rgt + "'")
        spy += """\"\"\"%s \"\"\"\n""" % ('Monomer A from ' + TAGLRGT[rgt])
        spy += """TAGL['%%s-%%s-monoB-unCP' %% (dbse, %-23s )] = """ % ("'" + rgt + "'")
        spy += """\"\"\"%s \"\"\"\n""" % ('Monomer B from ' + TAGLRGT[rgt])

# write subset geometry section
if route == 3:
    gpy += "# <<< Derived Geometry Strings >>>\n"
    gpy += "for rxn in HRXN:\n"
    gpy += "    GEOS['%s-%s-monoA-unCP' % (dbse, rxn)] = "
    gpy += "GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(1)\n"
    gpy += "    GEOS['%s-%s-monoB-unCP' % (dbse, rxn)] = "
    gpy += "GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(2)\n"
    gpy += "    GEOS['%s-%s-monoA-CP'   % (dbse, rxn)] = "
    gpy += "GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(1, 2)\n"
    gpy += "    GEOS['%s-%s-monoB-CP'   % (dbse, rxn)] = "
    gpy += "GEOS['%s-%s-dimer' % (dbse, rxn)].extract_fragments(2, 1)\n"

# arrange intermediate strings into final database file
fpy = open('%s.py' % (dbse), 'w')
fpy.write(spy)
fpy.write(gpy)
fpy.close()

# display customized advice for finishing off the database
final = """
   **  Congratulations, your database file %s.py has been constructed!

   **  To have a minimally functioning database, do the following:
""" % (dbse)

if line2 == 'comment' and isOS:
    final += """
       *  If not all neutral singlets, fill in correct charge and
          multiplicity for all reagents.
    """

if line2 == 'comment' and not isOS:
    final += """
       *  If not all neutral, fill in correct charge for all reagents.
    """

if route == 3 and line2 == 'cgmp':
    final += """
       *  The charge and multiplicity read in from line2 of the xyz files
          has been assigned to fragmentA, leaving fragmentB as a neutral
          singlet. If this is incorrect for any reagents, reapportion the
          charge and multiplicity correctly between fragments A & B.
    """

if route == 3 and line2 == 'comment':
    final += """
       *  If dimer and both subsystems are not neutral singlets, fill in
          correct charge and multiplicity for each subsystem.
    """

if route == 2:
    final += """
       *  Define the reagents that contribute to reach reaction by filling
          in the empty single quotes in ACTV. Add or delete lines as
          necessary for each reaction if more or fewer than three reagents
          contribute. See NHTBH.py as an example.

       *  Define the mathematical contribution of reagents to reactions
          by filling in a number (most often +1 or -1) for each reagent to
          the RXNM of each reaction. See NHTBH.py as an example.
    """

final += """
       *  Make sure the Psi4 driver can find your new database.
          M %s.py into INSTALLED_DIRECTORY/share/psi4/databases .
          Alternatively, add the directory containing %s.py into PYTHONPATH .
""" % (dbse, dbse)

final += """
   **  To enhance the functionality/documentation of your database, do the following:

       *  Rearrange the order of reactions in HRXN, as this will define
          the order for the database.

       *  Fill in the skeleton docstring at top of file, adding sources
          for geometries and any reference data. This info will show up
          in the online documentation.

       *  Fill in the comment lines of TAGL in plain text. These show up
          as banners in job output files.

       *  Fill in reference values (in kcal/mol) into BIND.

       *  If multiple sets of reference values are available, define each
          in an array BIND_ALTREF so that they can be called in a psi4
          input file as benchmark='ALTREF'. Add the new reference to the
          docstring. See S22.py as an example.

       *  Fill in the least computationally expensive 2-3 reactions into
          HRXN_SM and the most expensive into HRXN_LG so that they can be
          called in a psi4 input file as subset='small' or subset='large'.

       *  Define subsets of reactions such as in an array
          SUBSETARRAY=['reaction', 'reaction'] so that they can be called
          in a psi4 input file as subset='SUBSETARRAY'. Add the new subset
          option to to the docstring. See NBC10.py for a simple example or
          CFLOW.py for a complex example.
"""
print(final)
