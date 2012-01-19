import PsiMod
import re;
import os;

yes = re.compile(r'^(yes|true|on|1)', re.IGNORECASE)
no = re.compile(r'^(no|false|off|0)', re.IGNORECASE)
der0th = re.compile(r'^(0|none|energy)', re.IGNORECASE)
der1st = re.compile(r'^(1|first|gradient)', re.IGNORECASE)
der2nd = re.compile(r'^(2|second|hessian)', re.IGNORECASE)

def bad_option_syntax(line):
    print 'Unsupported syntax:\n\n%s\n\n' % (line)
    sys.exit(1)

def process_word_quotes(matchobj):
    dollar = matchobj.group(2)
    val    = matchobj.group(3)
    if(dollar):
        # This is a python variable, make sure that it starts with a letter
        if(re.match(r'^[A-Za-z][\w]*', val)):
            return val
        else:
            print "Invalid Python variable: %s" % val
            sys.exit(1)
    elif(re.match(r'^-?\d+\.?\d*(?:[Ee]-?\d+)?$', val)):
        # This must be a number, don't wrap it in quotes
        return val
    elif(re.match(r'^\'.*\'$', val) or re.match(r'^\".*\"$', val)):
        # This is already wrapped in quotes, do nothing
        return val
    else:
        # This must be a string
        return "\"%s\"" % val

def quotify(string):
    # This wraps anything that looks like a string in quotes, and removes leading
    # dollar signs from python variables
    wordre = re.compile(r'(([$]?)([-+()*.\w\"\']+))')
    string = wordre.sub(process_word_quotes, string)
    return string

def process_option(spaces, module, key, value, line):
    value  = quotify(value.strip())
    temp   = ""

    global_options = False
    module = module.upper()
    if(module == "GLOBALS" or module == "GLOBAL" or module == "" or module.isspace()):
        global_options = True

    if(global_options):
        # If it's really a global, we need slightly different syntax
        return spaces + "PsiMod.set_global_option(\"%s\", %s)\n" % (key, value)
    else:
        # It's a local option, so we need the module name in there too
        return spaces + "PsiMod.set_local_option(\"%s\", \"%s\", %s)\n" % (module, key, value)

def process_set_command(matchobj):
    result = ""
    module_string = ""
    if(matchobj.group(2)):
        module_string = matchobj.group(2)
    for module in module_string.split(","):
        result = result + process_option(matchobj.group(1), module, matchobj.group(3), matchobj.group(4), matchobj.group(0))
    return result

def process_set_commands(matchobj):
    spaces = matchobj.group(1)
    commands = matchobj.group(3)
    command_lines = re.split('\n', commands)
    # Remove trailing newline from each line
    map(lambda x: x.strip(), command_lines)
    result = ""
    module_string = ""
    command = ""
    if(matchobj.group(2)):
        module_string = matchobj.group(2)
    for module in module_string.split(","):
        for line in command_lines:
            # Chomp the trailing newline and accumulate
            command = command + line
            if not check_parentheses_and_brackets(command, 0):
                # If the brackets don't match up, we need to move on to the next line
                # and keep going, until they do match. Only then do we process the command
                continue
            # Ignore blank/empty lines
            if (not line or line.isspace()):
                continue
            matchobj = re.match(r'^\s*(\w+)[\s=]+(.*?)$', command)
            # Is the syntax correct? If so, process the line
            if matchobj:
                result = result + process_option(spaces, module, matchobj.group(1), matchobj.group(2), command)
                # Reset the string
                command = ""
            else:
                bad_option_syntax(command)
    return result


def process_molecule_command(matchobj):
    spaces = matchobj.group(1)
    name = matchobj.group(2)
    geometry = matchobj.group(3)
    molecule = spaces
    if name != "":
        molecule += '%s = ' % (name)

    molecule += 'geometry("""%s"""' % (geometry)
    if name != "":
        molecule += ',"%s"' % (name)

    molecule += ")\n"
    molecule += 'PsiMod.IO.set_default_namespace("%s")' % (name)

    return molecule

def process_extract_command(matchobj):
    spaces = matchobj.group(1)
    name = matchobj.group(2)
    extract = matchobj.group(0)
    extract += spaces + '%s.set_name("%s")' %(name,name)
    extract += "\n%sPsiMod.set_active_molecule(%s)" % (spaces,name)
    extract += '\n%sPsiMod.IO.set_default_namespace("%s")' % (spaces,name)

    return extract

def process_print_command(matchobj):
    spaces = matchobj.group(1)
    string = matchobj.group(2)

    printer = str(spaces)
    printer += "PsiMod.print_out(str(%s))\n" % str(string)

    return printer

def process_memory_command(matchobj):

    spacing = str(matchobj.group(1))
    sig = str(matchobj.group(2))
    units = str(matchobj.group(3))

    val = float(sig)
    memory_amount = val

    if (units.upper() == 'KB'):
        memory_amount = val*1000
    elif (units.upper() == 'MB'):
        memory_amount = val*1000000
    elif (units.upper() == 'GB'):
        memory_amount = val*1000000000

    command = "%sPsiMod.set_memory(%d)\n" % (spacing,int(memory_amount))
    return command

def process_basis_file(matchobj):
    spacing   = str(matchobj.group(1))
    basisfile = str(matchobj.group(2)).strip()
    command   = "%sPsiMod.add_user_basis_file(\"%s\")" % (spacing, basisfile)

    return command

def process_basis_block(matchobj):
    command_lines = re.split('\n', matchobj.group(2))
    spacing   = str(matchobj.group(1))
    result   = "%stemppsioman = PsiMod.IOManager.shared_object()" % spacing
    result  += "%spsi4tempscratchdir = temppsioman.get_file_path(100)" % spacing
    basislabel = re.compile(r'\s*\[([-*\(\)\w]+)\]\s*')

    # Start by looking for assign lines, and remove them
    label_re  = re.compile(r'^\s*assign\s*([A-Za-z]+\d+)\s+([-*\(\)\w]+)\s*(\w+)?\s*$')
    symbol_re = re.compile(r'^\s*assign\s*([A-Za-z]+)\s+([-*\(\)\w]+)\s*(\w+)?\s*$')
    all_re    = re.compile(r'^\s*assign\s*([-*\(\)\w]+)\s*(\w+)?\s*$')
    leftover_lines = []
    for line in command_lines:
        basistype = "BASIS"
        if(label_re.match(line)):
            m = label_re.match(line)
            if m.group(3): basistype = m.group(3).upper()
            result += "%sPsiMod.get_active_molecule().set_basis_by_label(\"%s\",\"%s\",\"%s\")" % (spacing, m.group(1), m.group(2), basistype)
        elif(symbol_re.match(line)):
            m = symbol_re.match(line)
            if m.group(3): basistype = m.group(3).upper()
            result += "%sPsiMod.get_active_molecule().set_basis_by_symbol(\"%s\",\"%s\",\"%s\")" % (spacing, m.group(1), m.group(2), basistype)
        elif(all_re.match(line)):
            m = all_re.match(line)
            if m.group(2): basistype = m.group(2).upper()
            result += "%sPsiMod.get_active_molecule().set_basis_all_atoms(\"%s\",\"%s\")" % (spacing, m.group(1), basistype)
        else:
            leftover_lines.append(line)

    # Now look for regular basis set definitions
    basisstring = ""
    for line in leftover_lines:
        # Ignore blank/empty lines
        if (not line or line.isspace()):
            continue
        m = re.match(basislabel, line)
        if(m):
            if(basisstring != ""):
                result += "%spsi4tempbasisfile = psi4tempscratchdir + \"%s\"" % (spacing, basisname)
                result += "%sPsiMod.add_user_basis_file(psi4tempbasisfile)" % (spacing)
                result += "%stemppsioman.write_scratch_file(psi4tempbasisfile, \"\"\"\n%s\"\"\")" % (spacing, basisstring)
                basisstring = ""
            basisname = PsiMod.BasisSet.make_filename(m.group(1))
        basisstring += line + "\n"
    if(basisstring != ""):
        result += "%spsi4tempbasisfile = psi4tempscratchdir + \"%s\"" % (spacing, basisname)
        result += "%sPsiMod.add_user_basis_file(psi4tempbasisfile)" % (spacing)
        result += "%stemppsioman.write_scratch_file(psi4tempbasisfile, \"\"\"\n%s\"\"\")" % (spacing, basisstring)
    return result


def process_external_command(matchobj):

    spacing = str(matchobj.group(1))
    name = str(matchobj.group(2))
    if (not name or name.isspace()):
        name = "extern"
    block = str(matchobj.group(3))
    lines = re.split('\n', block)

    extern =  "%sqmmm = QMMM()\n" % (spacing)

    NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

    # Comments are all removed by this point
    # 0. Remove blank lines
    re_blank = re.compile(r'^\s*$')
    lines2 = []
    for line in lines:
        mobj = re_blank.match(line)
        if (mobj):
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
        if (mobj):
            unit = mobj.group(1)
            if (unit == 'bohr' or unit == 'au' or unit == 'a.u.'):
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
        if (mobj):
            basis = mobj.group(1) 
        else:
            mobj = re_df_basis.match(line)
            if (mobj):
                df_basis_scf = mobj.group(1)
            else:
                lines2.append(line)
    lines = lines2

    # 4. Look for charge lines Z x y z, convert according to unit convention
    charge_re = re.compile(r'^\s*' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$')
    lines2 = []
    for line in lines: 
        mobj = charge_re.match(line)
        if (mobj):
            if (units == 'ang'):
                extern += '%sqmmm.addChargeAngstrom(%s,%s,%s,%s)\n' %(spacing,mobj.group(1),mobj.group(2),mobj.group(3),mobj.group(4))
            if (units == 'bohr'):
                extern += '%sqmmm.addChargeBohr(%s,%s,%s,%s)\n' %(spacing,mobj.group(1),mobj.group(2),mobj.group(3),mobj.group(4))
        else:
            lines2.append(line)
    lines = lines2

    # 5. Look for diffuse regions, which are XYZ molecules seperated by the usual -- lines
    spacer_re = re.compile(r'^\s*--\s*$')
    frags = []
    frags.append([])
    for line in lines:
        mobj = spacer_re.match(line)
        if (mobj):
            if (len(frags[len(frags)-1])):
                frags.append([])
        else:
            frags[len(frags)-1].append(line)

    extern += '%sextern_mol_temp = PsiMod.get_active_molecule()\n' %(spacing)

    mol_re = re.compile(r'\s*\S+\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$')
    lines = []
    for frag in frags:

        if not len(frag):
            continue

        extern += '%sexternal_diffuse = geometry("""\n' % (spacing)
        extern += '%s0 1\n' % (spacing) 

        for line in frag:
            if not mol_re.match(line):
                lines.append(line)
            else:
                extern += '%s%s\n' %(spacing,line)    

        extern += '%sunits %s\n' % (spacing, units) 
        extern += '%ssymmetry c1\n' % (spacing) 
        extern += '%sno_reorient\n' % (spacing) 
        extern += '%sno_com\n' % (spacing) 
        extern += '%s""")\n' % (spacing)
        extern += "%sdiffuse = Diffuse(external_diffuse,'%s','%s')\n" % (spacing, basis, df_basis_scf)
        extern += '%sdiffuse.fitScf()\n' % (spacing)
        extern += '%sqmmm.addDiffuse(diffuse)\n' % (spacing)
        extern += '\n'

    extern += '%sPsiMod.set_active_molecule(extern_mol_temp)\n' %(spacing)

    # 6. If there is anything left, the user messed up
    if (len(lines)):
        print 'Input parsing for external {}: Extra line(s) present:'
        for line in lines:
            print line
            sys.exit(1)

    # Return is actually an ExternalPotential, not a QMMM
    extern += '%sqmmm.populateExtern()\n' % (spacing)
    extern += '%s%s = qmmm.extern\n' % (spacing, name) 

    return extern

def check_parentheses_and_brackets(input_string, exit_on_error):
    # This returns 1 if the string's all matched up, 0 otherwise
    import collections

    # create left to right parenthesis mappings
    lrmap = {"(":")", "[":"]", "{":"}"}

    # derive sets of left and right parentheses
    lparens = set(lrmap.iterkeys())
    rparens = set(lrmap.itervalues())

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
                    print "Input error: extra %s" % ch
                    sys.exit(1)
            if lrmap[opench] != ch:
                # wrong type of parenthesis popped from stack
                all_matched = 0
                if exit_on_error:
                    print "Input error: %s closed with a %s" % (opench, ch)
                    sys.exit(1)
    if(len(parenstack) != 0):
        all_matched = 0
        if exit_on_error:
            print "Input error: Unmatched %s" % parenstack.pop()
            sys.exit(1)

    return all_matched


def parse_multiline_array(input_list):
    line = input_list.pop(0)
    # Keep adding lines to the current one, until all parens match up
    while not check_parentheses_and_brackets(line, 0):
        thisline = input_list.pop(0).strip()
        line += thisline
    return "%s\n" % line


def process_multiline_arrays(inputfile):
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
            newinput += "%s\n" % input_list.pop(0)
    return newinput


def process_input(raw_input):

    #NOTE: If adding mulitline data to the preprocessor, use ONLY the following syntax:
    #   function [objname] { ... }
    #   which has the regex capture group:
    #
    #   r'^(\s*?)FUNCTION\s*(\w*?)\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE
    #
    #   your function is in capture group #1
    #   your objname is in capture group #2
    #   your data is in capture group #3

    # Nuke all comments
    comment = re.compile(r'#.*')
    temp = re.sub(comment, '', raw_input)

    # Check the brackets and parentheses match up, as long as this is not a pickle input file
    if not re.search(r'pickle_kw', temp):
        check_parentheses_and_brackets(temp, 1)

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
    set_commands = re.compile(r'^(\s*?)set\s+([-,\w]*?)[\s=]*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(set_commands, process_set_commands, temp)

    # Process all individual "set (module_list) key  {[value_list] or $value or value}"
    set_command = re.compile(r'^(\s*?)set\s+(?:([-,\w]+)\s+)?(\w+)[\s=]+((\[.*\])|(\$?[-+,()*\.\w]+))\s*$', re.MULTILINE | re.IGNORECASE)
    temp = re.sub(set_command, process_set_command, temp)

    # Process "molecule name? { ... }"
    molecule = re.compile(r'^(\s*?)molecule[=\s]*(\w*?)\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(molecule, process_molecule_command, temp)

    # Process "external name? { ... }"
    external = re.compile(r'^(\s*?)external[=\s]*(\w*?)\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(external, process_external_command, temp)
    
    # Then remove repeated newlines
    multiplenewlines = re.compile(r'\n+')
    temp = re.sub(multiplenewlines, '\n', temp)

    # Process " extract"
    extract = re.compile(r'(\s*?)(\w+)\s*=\s*\w+\.extract_subsets.*', re.IGNORECASE)
    temp = re.sub(extract, process_extract_command, temp)

    # Process "print" and transform it to "PsiMod.print_out()"
    print_string = re.compile(r'(\s*?)print\s+(.*)',re.IGNORECASE)
    temp = re.sub(print_string,process_print_command,temp)

    # Process "memory ... "
    memory_string = re.compile(r'(\s*?)memory\s+([+-]?\d*\.?\d+)\s+([KMG]i?B)', re.IGNORECASE)
    temp = re.sub(memory_string,process_memory_command,temp)

    # Process "basis file ... "
    basis_file = re.compile(r'(\s*?)basis\s+file\s*(\b.*\b)\s*$', re.MULTILINE | re.IGNORECASE)
    temp = re.sub(basis_file,process_basis_file,temp)

    # Process "basis name { ... }"
    basis_block = re.compile(r'(\s*?)basis[=\s]*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(basis_block,process_basis_block,temp)

    # imports
    imports  = 'from PsiMod import *\n'
    imports += 'from physconst import *\n'
    imports += 'from molecule import *\n'
    imports += 'from driver import *\n'
    imports += 'from text import *\n'
    imports += 'from inpsight import *\n'
    imports += 'from wrappers import *\n'
    imports += 'from aliases import *\n'
    imports += 'from psiexceptions import *\n'
    imports += 'from util import *\n'
    imports += 'from qmmm import *\n'
    imports += 'from pubchem import *\n'
    imports += 'import pickle\n'
    imports += 'psi4_io = PsiMod.IOManager.shared_object()\n'

    # psirc (a baby PSIthon script that might live in ~/.psi4rc)
    psirc = ''
    homedir = os.path.expanduser('~')
    psirc_file = homedir + '/.psi4rc'
    if os.path.isfile(psirc_file):
        fh = open(psirc_file)
        psirc = fh.read()
        fh.close()

    blank_mol = 'geometry("""\n'
    blank_mol += '0 1\nH\nH 1 0.74\n'
    blank_mol += '""","blank_molecule_psi4_yo")\n'

    temp = imports + psirc + blank_mol + temp

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

    print "Result\n=========================="
    print result
