import re;
import os;

yes = re.compile(r'^(yes|true|on|1)', re.IGNORECASE)
no = re.compile(r'^(no|false|off|0)', re.IGNORECASE)

def bad_option_syntax(line):
    print 'Unsupported syntax:\n\n%s\n\n' % (line)
    exit()

def process_option(matchobj, option_string, line):
    spaces = matchobj.group(1)
    key   = matchobj.group(2).upper()
    value = matchobj.group(3).strip()
    temp = ""
    if re.match(r'^[-]?\d*\.?\d+$', value) or re.match(r'^\[.*\]$', value):
        # This is a number
        temp ='PsiMod.%s("%s", %s)' % (option_string, key, value)
    elif re.match(r'^\$', value):
        # Assume the user has set the variable in value
        temp = 'PsiMod.%s("%s", %s)' % (option_string, key, value[1:])
    else:
        # This must contain only alphanumeric (plus -,* for basis sets) or it's bad
        if re.match(r'^[-*\w]+$', value):
            temp = 'PsiMod.%s("%s", "%s")' % (option_string, key, value)
        else:
            bad_option_syntax(line)
    return spaces + temp + "\n"

def process_global_command(matchobj):
    return process_option(matchobj, "set_global_option", matchobj.group(0))

def process_set_command(matchobj):
    return process_option(matchobj, "set_option", matchobj.group(0)) 

def process_options(matchobj, command_string):
    spaces = matchobj.group(1)
    module = matchobj.group(2)
    commands = matchobj.group(3)
    command_lines = re.split('\n', commands)
    map(lambda x: x.strip(), command_lines)
    set_command = "set_global_option"
    result = []
    if module != "":
        set_command = "set_option"
        result.append('%sPsiMod.set_default_options_for_module("%s")\n' % (spaces, module.upper()))
    for line in command_lines:
        # Ignore blank/empty lines
        if (not line or line.isspace()):
            continue
        matchobj = re.match(r'^\s*()(\w+)[\s=]+(.*?)$', line)
        # Is the syntax correct? If so, process the line 
        if matchobj:
            result.append(process_option(matchobj, command_string,  line))
        else:
            bad_option_syntax(line)
    return (spaces).join(result)

def process_global_commands(matchobj):
    return process_options(matchobj, "set_global_option")

def process_set_commands(matchobj):
    return process_options(matchobj, "set_option")

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
 
    # Then remove repeated newlines
    multiplenewlines = re.compile(r'\n+')
    temp = re.sub(multiplenewlines, '\n', temp)

    # Process all "set name? { ... }"
    #set_commands = re.compile(r'^(\s*?)set\s*(\w*?)\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    set_commands = re.compile(r'^(\s*?)set\s*([\w\.]*?)\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(set_commands, process_set_commands, temp)

    # Process all individual "set key value"
    set_command = re.compile(r'(\s*?)set\s+(\w+)[\s=]+(.*?)($|#.*)', re.MULTILINE | re.IGNORECASE)
    temp = re.sub(set_command, process_set_command, temp)

    # Process all "global(s) { ... }"
    global_commands = re.compile(r'^(\s*?)()globals?\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(global_commands, process_global_commands, temp)

    # Process all individual "global key value"
    global_command = re.compile(r'^(\s*?)globals?\s+(\w+)[\s=]+(.*?)(\s*)$', re.MULTILINE | re.IGNORECASE)
    temp = re.sub(global_command, process_global_command, temp)

    # Process "molecule name? { ... }"
    molecule = re.compile(r'^(\s*?)molecule\s*(\w*?)\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(molecule, process_molecule_command, temp)

    # Process " extract"
    extract = re.compile(r'(\s*?)(\w+)\s*=\s*\w+\.extract_subsets.*', re.IGNORECASE)
    temp = re.sub(extract, process_extract_command, temp)

    # Process "print" and transform it to "PsiMod.print_out()"
    print_string = re.compile(r'(\s*?)print\s+(.*)',re.IGNORECASE)
    temp = re.sub(print_string,process_print_command,temp)

    # Process "memory ... "
    memory_string = re.compile(r'(\s*?)memory\s+([+-]?\d*\.?\d+)\s+([KMG]B)', re.IGNORECASE)   
    temp = re.sub(memory_string,process_memory_command,temp)

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

#this is a comment
globals {
    RI_BASIS_SCF STO-3G
}

global ncycles = 5

set scf {
    SCF_TYPE DF # test this too
    DOCC [3, 0, 1, 1]
    DIIS on
}

    globals freeze_core true

set docc =  [2, 0, 1, 1]
set docc   [2, 0, 1, 1]
    set docc   [2, 0, 1, 1]

    set mp2 {
        print  5
        print = 5
    }

""")

    print "Result\n=========================="
    print result
