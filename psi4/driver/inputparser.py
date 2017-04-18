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

"""Module with functions to parse the input file and convert
Psithon into standard Python. Particularly, forms psi4
module calls that access the C++ side of Psi4.

"""

## Force Python 3 print syntax, if this is python 2.X
#if sys.hexversion < 0x03000000:
from __future__ import print_function
from __future__ import absolute_import

import re
import os
import sys
import uuid
from psi4.driver import pubchem
from psi4.driver.p4util.exceptions import *
from psi4.driver.p4util.util import set_memory
from psi4 import core

# globally available regex strings
pubchemre = re.compile(r'^(\s*pubchem\s*:\s*(.*)\n)$', re.MULTILINE | re.IGNORECASE)

# inputfile contents to be preserved from the processor
literals = {}

# experimental - whether to run py statements as they're parsed from psithon
runalso = False


def bad_option_syntax(line):
    """Function to report bad syntax to screen and output file."""
    message = ('Unsupported syntax:\n\n%s\n\n' % (line))
    raise TestComparisonError(message)


def process_word_quotes(matchobj):
    """Function to determine if argument needs wrapping in quotes as string."""
    dollar = matchobj.group(2)
    val = matchobj.group(3)
    if dollar:
        # This is a python variable, make sure that it starts with a letter
        if re.match(r'^[A-Za-z][\w]*', val):
            return val
        else:
            message = ("Invalid Python variable: %s" % (val))
            raise TestComparisonError(message)
    elif re.match(r'^-?\d+\.?\d*(?:[Ee]-?\d+)?$', val):
        # This must be a number, don't wrap it in quotes
        return val
    elif re.match(r'^\'.*\'$', val) or re.match(r'^\".*\"$', val):
        # This is already wrapped in quotes, do nothing
        return val
    else:
        # This must be a string
        return "\"%s\"" % (val)


def quotify(string, isbasis=False):
    """Function to wrap anything that looks like a string in quotes
    and to remove leading dollar signs from python variables. When *basis*
    is True, allows commas, since basis sets may have commas and are assured to
    not involve arrays.

    """
    # This wraps anything that looks like a string in quotes, and removes leading
    # dollar signs from python variables
    if isbasis:
        wordre = re.compile(r'(([$]?)([-+()*.,\w\"\'/\\]+))')
    else:
        wordre = re.compile(r'(([$]?)([-+()*.\w\"\'/\\]+))')
    string = wordre.sub(process_word_quotes, string)
    return string

def dequotify(string):
    if string[0] == '"' and string[-1] == '"':
        return string[1:-1]
    else:
        return string


def process_option(spaces, module, key, value, line):
    """Function to process a line with set or in a set block
    into global/local domain and keyword/value.

    """
    module = module.upper()
    key = key.upper()
    isbasis = True if 'BASIS' in key else False
    value = quotify(value.strip(), isbasis=isbasis)

    if module == "GLOBALS" or module == "GLOBAL" or module == "" or module.isspace():
        # If it's really a global, we need slightly different syntax
        if runalso:
            core.set_global_option(key, dequotify(value))
        return "%score.set_global_option(\"%s\", %s)\n" % (spaces, key, value)
    else:
        # It's a local option, so we need the module name in there too
        if runalso:
            core.set_local_option(module, key, dequotify(value))
        return "%score.set_local_option(\"%s\", \"%s\", %s)\n" % (spaces, module, key, value)


def process_set_command(matchobj):
    """Function to process match of all individual ``set (module_list)
    key {[value_list] or $value or value}``.

    """
    result = ""
    module_string = ""
    if matchobj.group(2):
        module_string = matchobj.group(2)
    for module in module_string.split(","):
        result += process_option(matchobj.group(1), module, matchobj.group(3), matchobj.group(4), matchobj.group(0))
    return result


def process_set_commands(matchobj):
    """Function to process match of ``set name? { ... }``."""
    spaces = matchobj.group(1)
    commands = matchobj.group(3)
    command_lines = re.split('\n', commands)
    # Remove trailing newline from each line
    map(lambda x: x.strip(), command_lines)
    result = ""
    module_string = ""
    command = ""
    if matchobj.group(2):
        module_string = matchobj.group(2)
    for module in module_string.split(","):
        for line in command_lines:
            # Chomp the trailing newline and accumulate
            command += line
            if not check_parentheses_and_brackets(command, 0):
                # If the brackets don't match up, we need to move on to the next line
                # and keep going, until they do match. Only then do we process the command
                continue
            # Ignore blank/empty lines
            if not line or line.isspace():
                continue
            matchobj = re.match(r'^\s*(\w+)[\s=]+(.*?)$', command)
            # Is the syntax correct? If so, process the line
            if matchobj:
                result += process_option(spaces, module, matchobj.group(1), matchobj.group(2), command)
                # Reset the string
                command = ""
            else:
                bad_option_syntax(command)
    return result


def process_from_file_command(matchobj):
    """Function that process a match of ``from_file`` in molecule block."""
    string = matchobj.group(2)
    mol=core.mol_from_file(string,1)
    tempmol=[line for line in mol.split('\n') if line.strip() != '']
    mol2=set(tempmol)
    mol=""
    for i in mol2:
           mol+=i
           mol+="\n"
    return mol


def process_pubchem_command(matchobj):
    """Function to process match of ``pubchem`` in molecule block."""
    string = matchobj.group(2)
    if re.match(r'^\s*[0-9]+\s*$', string):
        # This is just a number - must be a CID
        pcobj = pubchem.PubChemObj(int(string), '', '')
        try:
            return pcobj.getMoleculeString()
        except Exception as e:
            return e.message
    else:
        # Search pubchem for the provided string
        try:
            results = pubchem.getPubChemResults(string)
        except Exception as e:
            return e.message

        # N.B. Anything starting with PubchemError will be handled correctly by the molecule parser
        # in libmints, which will just print the rest of the string and exit gracefully.
        if not results:
            # Nothing!
            return "PubchemError\n\tNo results were found when searching PubChem for %s.\n" % (string)
        elif len(results) == 1:
            # There's only 1 result - use it
            return results[0].getMoleculeString()
        else:
            # There are multiple results. Print and exit
            msg = "\tPubchemError\n"
            msg += "\tMultiple pubchem results were found. Replace\n\n\t\tpubchem:%s\n\n" % (string)
            msg += "\twith the Chemical ID number or exact name from one of the following and re-run.\n\n"
            msg += "\t Chemical ID     IUPAC Name\n\n"
            for result in results:
                msg += "%s" % (result)
                if result.name().lower() == string.lower():
                    #We've found an exact match!
                    return result.getMoleculeString()
            return msg


def process_molecule_command(matchobj):
    """Function to process match of ``molecule name? { ... }``."""
    spaces = matchobj.group(1)
    name = matchobj.group(2)
    geometry = matchobj.group(3)
    geometry = pubchemre.sub(process_pubchem_command, geometry)
    from_filere = re.compile(r'^(\s*from_file\s*:\s*(.*)\n)$', re.MULTILINE | re.IGNORECASE)
    geometry = from_filere.sub(process_from_file_command,geometry)
    molecule = spaces

    if name != "":
        if sys.version_info >= (3, 0):
            if not name.isidentifier():
                raise ValidationError('Molecule name not valid Python identifier: ' + name)
        else:
            if not re.match(r'^[^\d\W]\w*\Z', name):
                raise ValidationError('Molecule name not valid Python identifier: ' + name)

    molecule += 'core.efp_init()\n'  # clear EFP object before Molecule read in
    molecule += spaces

    if name != "":
        molecule += '%s = ' % (name)

    molecule += 'geometry("""%s"""' % (geometry)
    if name != "":
        molecule += ',"%s"' % (name)

    molecule += ")\n"
    molecule += '%score.IO.set_default_namespace("%s")' % (spaces, name)

    return molecule


def process_literal_blocks(matchobj):
    """Function to process match of ``literals_psi4_yo-...``."""
    return literals[matchobj.group(1)]


def process_cfour_command(matchobj):
    """Function to process match of ``cfour name? { ... }``."""
    spaces = matchobj.group(1)
    name = matchobj.group(2)
    cfourblock = matchobj.group(3)

    literalkey = str(uuid.uuid4())[:8]
    literals[literalkey] = cfourblock
    return "%score.set_global_option(\"%s\", \"\"\"%s\n\"\"\")\n" % \
        (spaces, 'LITERAL_CFOUR', 'literals_psi4_yo-' + literalkey)


def process_extract_command(matchobj):
    """Function to process match of ``extract_subsets``."""
    spaces = matchobj.group(1)
    name = matchobj.group(2)
    result = matchobj.group(0)
    result += '%s%s.set_name("%s")' % (spaces, name, name)
    result += "\n%score.set_active_molecule(%s)" % (spaces, name)
    result += '\n%score.IO.set_default_namespace("%s")' % (spaces, name)

    return result


def process_print_command(matchobj):
    """Function to process match of ``print`` and transform
    it to ``core.print_out()``.

    """
    spaces = matchobj.group(1)
    string = matchobj.group(2)

    return "%score.print_out(str(%s))\n" % (spaces, str(string))


def process_memory_command(matchobj):
    """Function to process match of ``memory ...``."""
    spaces = str(matchobj.group(1))
    sig = str(matchobj.group(2))
    units = str(matchobj.group(3))

    mem_in_bytes = set_memory(sig + units, execute=False)

    return "%score.set_memory_bytes(%d)\n" % (spaces, mem_in_bytes)


def basname(name):
    """Imitates BasisSet.make_filename() without the gbs extension"""
    return name.lower().replace('+', 'p').replace('*', 's').replace('(', '_').replace(')', '_').replace(',', '_')


def process_basis_block(matchobj):
    """Function to process match of ``basis name? { ... }``."""
    spaces = matchobj.group(1)
    basistype = matchobj.group(2).upper()
    name = matchobj.group(3)
    name = ('anonymous' + str(uuid.uuid4())[:8]) if name == '' else name
    cleanbas = basname(name).replace('-', '')  # further remove hyphens so can be function name
    command_lines = re.split('\n', matchobj.group(4))

    symbol_re = re.compile(r'^\s*assign\s+(?P<symbol>[A-Z]{1,3})\s+(?P<basis>[-+*\(\)\w]+)\s*$', re.IGNORECASE)
    label_re = re.compile(r'^\s*assign\s+(?P<label>(?P<symbol>[A-Z]{1,3})(?:(_\w+)|(\d+))?)\s+(?P<basis>[-+*\(\)\w]+)\s*$', re.IGNORECASE)
    all_re = re.compile(r'^\s*assign\s+(?P<basis>[-+*\(\)\w]+)\s*$', re.IGNORECASE)
    basislabel = re.compile(r'\s*\[\s*([-*\(\)\w]+)\s*\]\s*')

    result = """%sdef basisspec_psi4_yo__%s(mol, role):\n""" % (spaces, cleanbas)
    result += """%s    basstrings = {}\n""" % (spaces)

    # Start by looking for assign lines, and remove them
    leftover_lines = []
    for line in command_lines:
        if symbol_re.match(line):
            m = symbol_re.match(line)
            result += """%s    mol.set_basis_by_symbol("%s", "%s", role=role)\n""" % \
                (spaces, m.group('symbol'), m.group('basis'))

        elif label_re.match(line):
            m = label_re.match(line)
            result += """%s    mol.set_basis_by_label("%s", "%s", role=role)\n""" % \
                (spaces, m.group('label'), m.group('basis'))

        elif all_re.match(line):
            m = all_re.match(line)
            result += """%s    mol.set_basis_all_atoms("%s", role=role)\n""" % \
                (spaces, m.group('basis'))

        else:
            # Ignore blank lines and accumulate remainder
            if line and not line.isspace():
                leftover_lines.append(line.strip())

    # Now look for regular basis set definitions
    basblock = list(filter(None, basislabel.split('\n'.join(leftover_lines))))
    if len(basblock) == 1:
        if len(result.split('\n')) == 3:
            # case with no [basname] markers where whole block is contents of gbs file
            result += """%s    mol.set_basis_all_atoms("%s", role=role)\n""" % \
                (spaces, name)
            result += """%s    basstrings['%s'] = \"\"\"\n%s\n\"\"\"\n""" % \
                (spaces, basname(name), basblock[0])
        else:
            message = ("Conflicting basis set specification: assign lines present but shells have no [basname] label.""")
            raise TestComparisonError(message)
    else:
        # case with specs separated by [basname] markers
        for idx in range(0, len(basblock), 2):
            result += """%s    basstrings['%s'] = \"\"\"\n%s\n\"\"\"\n""" % \
                (spaces, basname(basblock[idx]), basblock[idx + 1])

    result += """%s    return basstrings\n""" % (spaces)
    result += """{}qcdb.libmintsbasisset.basishorde['{}'] = {}\n""" \
              .format(spaces, name.upper(), 'basisspec_psi4_yo__' + cleanbas)
    result += """%score.set_global_option(\"%s\", \"%s\")""" % (spaces, basistype, name)
    return result


def process_pcm_command(matchobj):
    """Function to process match of ``pcm name? { ... }``."""
    spacing = str(matchobj.group(1)) # Ignore..
    name = str(matchobj.group(2)) # Ignore..
    block = str(matchobj.group(3)) # Get input to PCMSolver
    # Setup unique scratch directory and move in, as done for other add-ons
    psioh = core.IOManager.shared_object()
    psio = core.IO.shared_object()
    os.chdir(psioh.get_default_path())
    pcm_tmpdir = 'psi.' + str(os.getpid()) + '.' + psio.get_default_namespace() + \
        'pcmsolver.' + str(uuid.uuid4())[:8]
    if os.path.exists(pcm_tmpdir) is False:
        os.mkdir(pcm_tmpdir)
    os.chdir(pcm_tmpdir)
    with open('pcmsolver.inp', 'w') as handle:
        handle.write(block)
    import pcmsolver
    pcmsolver.parse_pcm_input('pcmsolver.inp')
    return "" # The file has been written to disk; nothing needed in Psi4 input


def process_external_command(matchobj):
    """Function to process match of ``external name? { ... }``."""
    spaces = str(matchobj.group(1))
    name = str(matchobj.group(2))
    if not name or name.isspace():
        name = "extern"
    block = str(matchobj.group(3))
    lines = re.split('\n', block)

    extern = "%sqmmm = QMMM()\n" % (spaces)

    NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

    # Comments are all removed by this point
    # 0. Remove blank lines
    re_blank = re.compile(r'^\s*$')
    lines2 = []
    for line in lines:
        mobj = re_blank.match(line)
        if mobj:
            pass
        else:
            lines2.append(line)
    lines = lines2

    # 1. Look for units [ang|bohr|au|a.u.] defaults to ang
    re_units = re.compile(r'^\s*units?[\s=]+((ang)|(angstrom)|(bohr)|(au)|(a\.u\.))$\s*', re.IGNORECASE)
    units = 'ang'
    lines2 = []
    for line in lines:
        mobj = re_units.match(line)
        if mobj:
            unit = mobj.group(1)
            if unit in ['bohr', 'au', 'a.u.']:
                units = 'bohr'
            else:
                units = 'ang'
        else:
            lines2.append(line)
    lines = lines2

    # 2. Look for basis basisname, defaults to cc-pvdz
    # 3. Look for df_basis_scf basisname, defaults to cc-pvdz-jkfit
    re_basis = re.compile(r'\s*basis[\s=]+(\S+)\s*$', re.IGNORECASE)
    re_df_basis = re.compile(r'\s*df_basis_scf[\s=]+(\S+)\s*$', re.IGNORECASE)
    basis = 'cc-pvdz'
    df_basis_scf = 'cc-pvdz-jkfit'
    lines2 = []
    for line in lines:
        mobj = re_basis.match(line)
        if mobj:
            basis = mobj.group(1)
        else:
            mobj = re_df_basis.match(line)
            if mobj:
                df_basis_scf = mobj.group(1)
            else:
                lines2.append(line)
    lines = lines2

    # 4. Look for charge lines Z x y z, convert according to unit convention
    charge_re = re.compile(r'^\s*' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$')
    lines2 = []
    for line in lines:
        mobj = charge_re.match(line)
        if mobj:
            if units == 'ang':
                extern += '%sqmmm.addChargeAngstrom(%s,%s,%s,%s)\n' % (spaces, mobj.group(1), mobj.group(2), mobj.group(3), mobj.group(4))
            if units == 'bohr':
                extern += '%sqmmm.addChargeBohr(%s,%s,%s,%s)\n' % (spaces, mobj.group(1), mobj.group(2), mobj.group(3), mobj.group(4))
        else:
            lines2.append(line)
    lines = lines2

    # 5. Look for diffuse regions, which are XYZ molecules seperated by the usual -- lines
    spacer_re = re.compile(r'^\s*--\s*$')
    frags = []
    frags.append([])
    for line in lines:
        mobj = spacer_re.match(line)
        if mobj:
            if len(frags[len(frags) - 1]):
                frags.append([])
        else:
            frags[len(frags) - 1].append(line)

    extern += '%sextern_mol_temp = core.get_active_molecule()\n' % (spaces)

    mol_re = re.compile(r'\s*\S+\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$')
    lines = []
    for frag in frags:

        if not len(frag):
            continue

        extern += '%sexternal_diffuse = geometry("""\n' % (spaces)
        extern += '%s0 1\n' % (spaces)

        for line in frag:
            if not mol_re.match(line):
                lines.append(line)
            else:
                extern += '%s%s\n' % (spaces, line)

        extern += '%sunits %s\n' % (spaces, units)
        extern += '%ssymmetry c1\n' % (spaces)
        extern += '%sno_reorient\n' % (spaces)
        extern += '%sno_com\n' % (spaces)
        extern += '%s""")\n' % (spaces)
        extern += "%sdiffuse = Diffuse(external_diffuse,'%s','%s')\n" % (spaces, basis, df_basis_scf)
        extern += '%sdiffuse.fitScf()\n' % (spaces)
        extern += '%sqmmm.addDiffuse(diffuse)\n' % (spaces)
        extern += '\n'

    extern += '%score.set_active_molecule(extern_mol_temp)\n' % (spaces)

    # 6. If there is anything left, the user messed up
    if len(lines):
        print('Input parsing for external {}: Extra line(s) present:')
        for line in lines:
            raise TestComparisonError(line)

    # Return is actually an ExternalPotential, not a QMMM
    extern += '%sqmmm.populateExtern()\n' % (spaces)
    extern += '%s%s = qmmm.extern\n' % (spaces, name)

    extern += '%score.set_global_option_python("EXTERN", extern)\n' % (spaces)

    return extern


def check_parentheses_and_brackets(input_string, exit_on_error):
    """Function to check that all parenthesis and brackets
    in *input_string* are paired. On that condition, *exit_on_error* =1,
    otherwise 0.

    """
    # This returns 1 if the string's all matched up, 0 otherwise
    import collections

    # create left to right parenthesis mappings
    lrmap = {"(": ")", "[": "]", "{": "}"}

    # derive sets of left and right parentheses
    lparens = set(lrmap.keys())
    rparens = set(lrmap.values())

    parenstack = collections.deque()
    all_matched = 1
    for ch in input_string:
        if ch in lparens:
            parenstack.append(ch)
        elif ch in rparens:
            opench = ""
            try:
                opench = parenstack.pop()
            except IndexError:
                # Run out of opening parens
                all_matched = 0
                if exit_on_error:
                    message = ("Input error: extra %s" % (ch))
                    raise TestComparisonError(message)
            if lrmap[opench] != ch:
                # wrong type of parenthesis popped from stack
                all_matched = 0
                if exit_on_error:
                    message = ("Input error: %s closed with a %s" % (opench, ch))
                    raise TestComparisonError(message)
    if len(parenstack) != 0:
        all_matched = 0
        if exit_on_error:
            message = ("Input error: Unmatched %s" % (parenstack.pop()))
            raise TestComparisonError(message)

    return all_matched


def parse_multiline_array(input_list):
    """Function to squash multiline arrays into a single line
    until all parentheses and brackets are fully paired.

    """
    line = input_list.pop(0)
    # Keep adding lines to the current one, until all parens match up
    while not check_parentheses_and_brackets(line, 0):
        thisline = input_list.pop(0).strip()
        line += thisline
    return "%s\n" % (line)


def process_multiline_arrays(inputfile):
    """Function to find array inputs that are spread across multiple
    lines and squash them into a single line.

    """
    # This function takes multiline array inputs, and puts them on a single line
    # Start by converting the input to a list, splitting at newlines
    input_list = inputfile.split("\n")
    set_re = re.compile(r'^(\s*?)set\s+(?:([-,\w]+)\s+)?(\w+)[\s=]+\[.*', re.IGNORECASE)
    newinput = ""
    while len(input_list):
        line = input_list[0]
        if set_re.match(line):
            # We've found the start of a set matrix [ .... line - hand it off for more checks
            newinput += parse_multiline_array(input_list)
        else:
            # Nothing to do - just add the line to the string
            newinput += "%s\n" % (input_list.pop(0))
    return newinput


def process_input(raw_input, print_level=1):
    """Function to preprocess *raw input*, the text of the input file, then
    parse it, validate it for format, and convert it into legitimate Python.
    *raw_input* is printed to the output file unless *print_level* =0. Does
    a series of regular expression filters, where the matching portion of the
    input is replaced by the output of the corresponding function (in this
    module) call. Returns a string concatenating module import lines, a copy
    of the user's .psi4rc files, a setting of the scratch directory, a dummy
    molecule, and the processed *raw_input*.

    """
    # Check if the infile is actually an outfile (yeah we did)
    psi4_id = re.compile(r'Psi4: An Open-Source Ab Initio Electronic Structure Package')
    if re.search(psi4_id, raw_input):
        input_lines = raw_input.split("\n")
        input_re = re.compile(r'^\s*?\=\=> Input File <\=\=')
        input_start = -1
        for line_count in range(len(input_lines)):
            line = input_lines[line_count]
            if re.match(input_re, line):
                input_start = line_count + 3
                break

        stop_re = re.compile(r'^-{74}')
        input_stop = -1
        for line_count in range(input_start, len(input_lines)):
            line = input_lines[line_count]
            if re.match(stop_re, line):
                input_stop = line_count
                break

        if input_start == -1 or input_stop == -1:
            message = ('Cannot extract infile from outfile.')
            raise TestComparisonError(message)

        raw_input = '\n'.join(input_lines[input_start:input_stop])
        raw_input += '\n'

    # Echo the infile on the outfile
    if print_level > 0:
        core.print_out("\n  ==> Input File <==\n\n")
        core.print_out("--------------------------------------------------------------------------\n")
        core.print_out(raw_input)
        core.print_out("--------------------------------------------------------------------------\n")
        core.flush_outfile()

    #NOTE: If adding mulitline data to the preprocessor, use ONLY the following syntax:
    #   function [objname] { ... }
    #   which has the regex capture group:
    #
    #   r'^(\s*?)FUNCTION\s*(\w*?)\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE
    #
    #   your function is in capture group #1
    #   your objname is in capture group #2
    #   your data is in capture group #3

    # Sections that are truly to be taken literally (spaces included)
    #   Must be stored then subbed in the end to escape the normal processing

    # Process "cfour name? { ... }"
    cfour = re.compile(r'^(\s*?)cfour[=\s]*(\w*?)\s*\{(.*?)\}',
                          re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(cfour, process_cfour_command, raw_input)

    # Return from handling literal blocks to normal processing

    # Nuke all comments
    comment = re.compile(r'(^|[^\\])#.*')
    temp = re.sub(comment, '', temp)
    # Now, nuke any escapes from comment lines
    comment = re.compile(r'\\#')
    temp = re.sub(comment, '#', temp)

    # Check the brackets and parentheses match up, as long as this is not a pickle input file
    #if not re.search(r'pickle_kw', temp):
    #    check_parentheses_and_brackets(temp, 1)

    # First, remove everything from lines containing only spaces
    blankline = re.compile(r'^\s*$')
    temp = re.sub(blankline, '', temp, re.MULTILINE)

    # Look for things like
    # set matrix [
    #              [ 1, 2 ],
    #              [ 3, 4 ]
    #            ]
    # and put them on a single line
    temp = process_multiline_arrays(temp)

    # Process all "set name? { ... }"
    set_commands = re.compile(r'^(\s*?)set\s*([-,\w]*?)[\s=]*\{(.*?)\}',
                              re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(set_commands, process_set_commands, temp)

    # Process all individual "set (module_list) key  {[value_list] or $value or value}"
    # N.B. We have to be careful here, because \s matches \n, leading to potential problems
    # with undesired multiline matches.  Better the double-negative [^\S\n] instead, which
    # will match any space, tab, etc., except a newline
    set_command = re.compile(r'^(\s*?)set\s+(?:([-,\w]+)[^\S\n]+)?(\w+)(?:[^\S\n]|=)+((\[.*\])|(\$?[-+,*()\.\w]+))\s*$',
                             re.MULTILINE | re.IGNORECASE)
    temp = re.sub(set_command, process_set_command, temp)

    # Process "molecule name? { ... }"
    molecule = re.compile(r'^(\s*?)molecule[=\s]*(\S*?)\s*\{(.*?)\}',
                          re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(molecule, process_molecule_command, temp)

    # Process "external name? { ... }"
    external = re.compile(r'^(\s*?)external[=\s]*(\w*?)\s*\{(.*?)\}',
                          re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(external, process_external_command, temp)

    # Process "pcm name? { ... }"
    pcm = re.compile(r'^(\s*?)pcm[=\s]*(\w*?)\s*\{(.*?)^\}',
                          re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(pcm, process_pcm_command, temp)

    # Then remove repeated newlines
    multiplenewlines = re.compile(r'\n+')
    temp = re.sub(multiplenewlines, '\n', temp)

    # Process " extract"
    extract = re.compile(r'(\s*?)(\w+)\s*=\s*\w+\.extract_subsets.*',
                         re.IGNORECASE)
    temp = re.sub(extract, process_extract_command, temp)

    # Process "print" and transform it to "core.print_out()"
    #print_string = re.compile(r'(\s*?)print\s+(.*)', re.IGNORECASE)
    #temp = re.sub(print_string, process_print_command, temp)

    # Process "memory ... "
    memory_string = re.compile(r'(\s*?)memory\s+(\d*\.?\d+)\s*([KMGTPBE]i?B)',
                               re.IGNORECASE)
    temp = re.sub(memory_string, process_memory_command, temp)

    # Process "basis name? { ... }"
    basis_block = re.compile(r'^(\s*?)(basis|df_basis_scf|df_basis_mp2|df_basis_cc|df_basis_sapt)[=\s]*(\w*?)\s*\{(.*?)\}',
                             re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(basis_block, process_basis_block, temp)

    # Process literal blocks by substituting back in
    lit_block = re.compile(r'literals_psi4_yo-(\d*\d)')
    temp = re.sub(lit_block, process_literal_blocks, temp)

    future_imports = []
    def future_replace(m):
        future_imports.append(m.group(0))
        return ''

    future_string = re.compile('^from __future__ import .*$', flags=re.MULTILINE)
    temp = re.sub(future_string, future_replace, temp)

    # imports
    imports = '\n'.join(future_imports) + '\n'
    imports += 'import psi4\n'
    imports += 'from psi4 import *\n'
    imports += 'from psi4.core import *\n'
    imports += 'from psi4.driver.diatomic import anharmonicity\n'
    imports += 'from psi4.driver.gaussian_n import *\n'
    imports += 'from psi4.driver.aliases import *\n'
    imports += 'from psi4.driver.driver_cbs import xtpl_highest_1, scf_xtpl_helgaker_2, scf_xtpl_helgaker_3, corl_xtpl_helgaker_2\n'
    imports += 'from psi4.driver.wrapper_database import database, db, DB_RGT, DB_RXN\n'
    imports += 'from psi4.driver.wrapper_autofrag import auto_fragments\n'
    imports += 'from psi4.driver.constants.physconst import *\n'
    imports += 'psi4_io = core.IOManager.shared_object()\n'

    # psirc (a baby PSIthon script that might live in ~/.psi4rc)
    psirc_file = os.path.expanduser('~') + os.path.sep + '.psi4rc'
    if os.path.isfile(psirc_file):
        fh = open(psirc_file)
        psirc = fh.read()
        fh.close()
        psirc = psirc.replace('psi4.IOManager', 'psi4.core.IOManager')
    else:
        psirc = ''

    # Override scratch directory if user specified via env_var
    scratch = ''
    scratch_env = core.get_environment('PSI_SCRATCH')
    if len(scratch_env):
        scratch += 'psi4_io.set_default_path("%s")\n' % (scratch_env)

    blank_mol = 'geometry("""\n'
    blank_mol += '0 1\nH\nH 1 0.74\n'
    blank_mol += '""","blank_molecule_psi4_yo")\n'

    temp = imports + psirc + scratch + blank_mol + temp

    # Move up the psi4.core namespace
    for func in dir(core):
        temp = temp.replace("psi4." + func, "psi4.core." + func)

    # Move pseudonamespace for physconst into proper namespace
    from psi4.driver.p4util import constants
    for pc in dir(constants.physconst):
        if not pc.startswith('__'):
            temp = temp.replace('psi_' + pc, 'psi4.constants.' + pc)

    return temp


if __name__ == "__main__":
    result = process_input("""
molecule h2 {
H
H 1 R

R = .9
}

set basis 6-31G**

""")

    print("Result\n==========================")
    print(result)
