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
    basisstring = ""
    for line in command_lines:
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

def process_basis_assign_block(matchobj):
    spacing       = str(matchobj.group(1))
    basistype     = str(matchobj.group(2)).upper()
    block         = str(matchobj.group(3)).upper()
    command_lines = re.split('\n', block)
    label_re      = re.compile(r'^\s*([A-Za-z]+\d+)\s+([-*\(\)\w]+)\s*$')
    symbol_re     = re.compile(r'^\s*([A-Za-z]+)\s+([-*\(\)\w]+)\s*$')
    number_re     = re.compile(r'^\s*(\d+)\s+([-*\(\)\w]+)\s*$')
    all_re        = re.compile(r'^\s*([-*\(\)\w]+)\s*$')
    command = ""
    if(not basistype or basistype.isspace()):
        basistype = "BASIS"
    for line in command_lines:
        if (not line or line.isspace()):
            continue
    #void set_basis_all_atoms(const std::string& name, const std::string& type="BASIS");
    #void set_basis_by_symbol(const std::string& symbol, const std::string& name, const std::string& type="BASIS");
    #void set_basis_by_number(int number, const std::string& name, const std::string& type="BASIS");
    #void set_basis_by_label(const std::string& label, const std::string& name, const std::string& type="BASIS");
        if(label_re.match(line)):
            m = label_re.match(line)
            command += "%sPsiMod.get_active_molecule().set_basis_by_label(\"%s\",\"%s\",\"%s\")" % (spacing, m.group(1), m.group(2), basistype)
        elif(symbol_re.match(line)):
            m = symbol_re.match(line)
            command += "%sPsiMod.get_active_molecule().set_basis_by_symbol(\"%s\",\"%s\",\"%s\")" % (spacing, m.group(1), m.group(2), basistype)
        elif(number_re.match(line)):
            m = number_re.match(line)
            command += "%sPsiMod.get_active_molecule().set_basis_by_number(%s,\"%s\",\"%s\")" % (spacing, m.group(1), m.group(2), basistype)
        elif(all_re.match(line)):
            m = all_re.match(line)
            command += "%sPsiMod.get_active_molecule().set_basis_all_atoms(\"%s\",\"%s\")" % (spacing, m.group(1), basistype)
        else:
            bad_option_syntax(line)
    return command

def process_extern_block(matchobj):

    spacing = str(matchobj.group(1))
    name = str(matchobj.group(2))
    if (not name or name.isspace()):
        name = "EXTERN"
    block = str(matchobj.group(3))
    lines = re.split('\n', block)

    extern = "%s%s = PsiMod.ExternalPotential()" % (spacing, name)
    extern += '%sactivate_potential(%s)' %(spacing,name)
    extern += '%s%s.setName("%s")' % (spacing,name,name)

    addType = "None"

    NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

    charge_re = re.compile(r'^\s*(charges?)\s*$', re.IGNORECASE)
    dip_re = re.compile(r'^\s*(dipoles?)\s*$', re.IGNORECASE)
    quad_re = re.compile(r'^\s*(quadrupoles?)\s*$', re.IGNORECASE)
    spacer_re = re.compile(r'^\s*\*{4}\s*$')

    c_re = re.compile(r'^\s*' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$')
    d_re = re.compile(r'^\s*' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$')
    q_re = re.compile(r'^\s*' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$')

    for line in lines:
        if (not line or line.isspace()):
            continue

        if spacer_re.match(line):
            continue
        elif charge_re.match(line):
            addType = "Charge"
            continue
        else:
            if addType == "None":
                print "First Line of extern section must be a type of potential (charge, dipole, etc)\n"
                bad_option_syntax(line)

            if addType == "Charge":
                mobj = c_re.match(line)
                if not mobj:
                    bad_option_syntax(line)
                extern += "%s%s.addCharge(%s,%s,%s,%s)" % (spacing, name, mobj.group(1),\
                     mobj.group(2), mobj.group(3), mobj.group(4))

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

    # Process all "set name? { ... }"
    set_commands = re.compile(r'^(\s*?)set\s+([-,\w]*?)[\s=]*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(set_commands, process_set_commands, temp)

    # Process all individual "set (module_list) key  {[value_list] or $value or value}"
    set_command = re.compile(r'^(\s*?)set\s+(?:([-,\w]+)\s+)?(\w+)[\s=]+((\[.*\])|(\$?[-+,()*\.\w]+))\s*$', re.MULTILINE | re.IGNORECASE)
    temp = re.sub(set_command, process_set_command, temp)

    # Process "molecule name? { ... }"
    molecule = re.compile(r'^(\s*?)molecule[=\s]*(\w*?)\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(molecule, process_molecule_command, temp)

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

    # Process "basis assign (name) { ... }"
    basis_assign_block = re.compile(r'(\s*?)basis\s+assign\s*(\w*)\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(basis_assign_block,process_basis_assign_block,temp)

    # Process "basis name { ... }"
    basis_block = re.compile(r'(\s*?)basis[=\s]*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(basis_block,process_basis_block,temp)

    # Process "extern name { ... }"
    extern_block = re.compile(r'(\s*?)extern\s*(\w*?)[=\s]*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(extern_block,process_extern_block,temp)

    # imports
    imports  = 'from PsiMod import *\n'
    imports += 'from molecule import *\n'
    imports += 'from driver import *\n'
    imports += 'from text import *\n'
    imports += 'from inpsight import *\n'
    imports += 'from wrappers import *\n'
    imports += 'from psiexceptions import *\n'
    imports += 'from util import *\n'
    imports += 'import pickle\n'

    # psirc (a baby PSithon script that might live in ~/.psirc
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
