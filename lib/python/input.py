import re;
import os;

yes = re.compile(r'^(yes|true|on|1)', re.IGNORECASE)
no = re.compile(r'^(no|false|off|0)', re.IGNORECASE)

def bad_option_syntax(line):
    print 'Unsupported syntax:\n\n%s\n\n' % (line)
    exit()

def process_option(spaces, module, key, value, line):
    value  = value.strip()
    temp   = ""

    global_options = False
    module = module.upper()
    if(module == "GLOBALS" or module == "GLOBAL" or module == "" or module.isspace()):
        global_options = True

    # Start by assuming that we're setting a local option
    command_string = "PsiMod.set_local_option(\"%s\", \"%s\", " % (module, key)
    if(global_options):
        # If it's really a global, we need slightly different syntax
        command_string = "PsiMod.set_global_option(\"%s\", " % (key)

    if re.match(r'^[-]?\d*\.?\d+$', value) or re.match(r'^\[.*\]$', value):
        # This is a number
        temp ='%s %s)' % (command_string, value)
    elif re.match(r'^\$', value):
        # Assume the user has set the variable in value
        temp = '%s %s)' % (command_string, value[1:])
    else:
        # This must contain only alphanumeric (plus -,* for basis sets) or it's bad
        if re.match(r'^[-*\w]+$', value):
            temp = '%s "%s")' % (command_string, value)
        else:
            bad_option_syntax(line)
    return spaces + temp + "\n"

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
    map(lambda x: x.strip(), command_lines)
    result = ""
    module_string = ""
    if(matchobj.group(2)):
        module_string = matchobj.group(2)
    for module in module_string.split(","):
        for line in command_lines:
            # Ignore blank/empty lines
            if (not line or line.isspace()):
                continue
            matchobj = re.match(r'^\s*(\w+)[\s=]+(.*?)$', line)
            # Is the syntax correct? If so, process the line
            if matchobj:
                result = result + process_option(spaces, module, matchobj.group(1), matchobj.group(2), line)
            else:
                bad_option_syntax(line)
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

    spacing = str(matchobj.group(1))
    basisfile = str(matchobj.group(2)).strip()

    command = "%sPsiMod.add_user_basis_file(\"%s\")" % (spacing, basisfile)

    return command

def process_basis_block(matchobj):

    spacing = str(matchobj.group(1))
    filename = str(matchobj.group(2))
    block = str(matchobj.group(3))

    # TODO: Rob's updated PSIOManager is going to be updated to make this easy.
    command = "%s# do something with the block of basis set data for %s" % (spacing, filename)

    return command

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

    # First, remove everything from lines containing only spaces
    blankline = re.compile(r'^\s*$')
    temp = re.sub(blankline, '', temp, re.MULTILINE)

    # Process all "set name? { ... }"
    set_commands = re.compile(r'^(\s*?)set\s+([-,\w]*?)[\s=]*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(set_commands, process_set_commands, temp)

    # Process all individual "set (module_list) key  {[value_list] or $value or value}"
    set_command = re.compile(r'^(\s*?)set\s+(?:([-,\w]+)\s+)?(\w+)[\s=]+((\[.*\])|(\$?[-*\.\w]+))\s*$', re.MULTILINE | re.IGNORECASE)
    temp = re.sub(set_command, process_set_command, temp)

    # Process "molecule name? { ... }"
    molecule = re.compile(r'^(\s*?)molecule\s*(\w*?)\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
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
    memory_string = re.compile(r'(\s*?)memory\s+([+-]?\d*\.?\d+)\s+([KMG]B)', re.IGNORECASE)
    temp = re.sub(memory_string,process_memory_command,temp)

    # Process "basis file ... "
    basis_file = re.compile(r'(\s*?)basis\s+file\s*(\b.*\b)\s*$', re.MULTILINE | re.IGNORECASE)
    temp = re.sub(basis_file,process_basis_file,temp)

    # Process "basis name { ... }"
    basis_block = re.compile(r'(\s*?)basis\s+([-\(\)\+\*\w]*)\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(basis_block,process_basis_block,temp)

    # imports
    imports  = 'from PsiMod import *\n'
    imports += 'from opt import *\n'
    imports += 'from molecule import *\n'
    imports += 'from driver import *\n'
    imports += 'from text import *\n'
    imports += 'from wrappers import *\n'
    imports += 'from psiexceptions import *\n'
    imports += 'from util import *\n'

    # psirc (a baby PSithon script that might live in ~/.psirc
    psirc = ''
    homedir = os.path.expanduser('~')
    psirc_file = homedir + '/.psi4rc'
    if os.path.isfile(psirc_file):
        fh = open(psirc_file)
        psirc = fh.read()
        fh.close()

    temp = imports + psirc + temp

    return temp

if __name__ == "__main__":
    result = process_input("""
molecule h2 {
H
H 1 R

R = .9
}

set basis 6-31G**

#this is a comment
set globals {
    RI_BASIS_SCF STO-3G
}

set scf,ccsd  = {
    print 1
    DOCC [3, 0, 1, 1]
    DIIS on
}

    set globals freeze_core true
    set global freeze_core = true

    set global  docc   [2, 0, 1, 1]
    set  ss  = 3.0
    set  ss   3.0
    set scf socc  [23]
    set scf,ccsd docc  [34,43]
    set dostuff   1
    set globals,ccsd do_more_stuff $foo

    set mp2 {
        print  5
        print = 5
    }

basis file ~/basis/sto3g.gbs
basis file ~/basis sets/cc-pvdz.gbs

basis sto-3g {
 blah blah
}

""")

    print "Result\n=========================="
    print result
