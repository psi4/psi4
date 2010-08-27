import re;

yes = re.compile(r'^(yes|true|on)', re.IGNORECASE)
no = re.compile(r'^(no|false|off)', re.IGNORECASE)

def process_global_command(matchobj):
    spaces = matchobj.group(1)
    key   = matchobj.group(2).upper()
    value = matchobj.group(3).strip()

    if re.match(r'^[-]?\d*\.?\d+$', value) or re.match(r'^\[.*\]$', value):
        return spaces + 'PsiMod.set_global_option("%s", %s)' % (key, value)
    elif re.match(r'^\$', value):
        # Assume the user as set the variable in value
        return spaces + 'PsiMod.set_global_option("%s", %s)' % (key, value[1:])
    elif re.match(yes, value):
        return spaces + 'PsiMod.set_global_option("%s", True)' % (key)
    elif re.match(no, value):
        return spaces + 'PsiMod.set_global_option("%s", False)' % (key)
    else:
        return spaces + 'PsiMod.set_global_option("%s", "%s")' % (key, value)

def process_globals_command(matchobj):
    spaces = matchobj.group(1)
    commands = matchobj.group(3)
    command_lines = re.split('\n', commands)

    result = []
    for line in command_lines:
        temp = re.sub(r'#.*', "", line)
        result.append(re.sub(r'^\s*()(\w+)\s+(.*)($|#.*)', process_global_command, temp))

    set_commands = spaces
    set_commands += (spaces).join(result)

    return set_commands

def process_set_command(matchobj):
    spaces = matchobj.group(1)
    key   = matchobj.group(2).upper()
    value = matchobj.group(3).strip()

    if re.match(r'^[-]?\d*\.?\d+$', value) or re.match(r'^\[.*\]$', value):
        return spaces + 'PsiMod.set_option("%s", %s)' % (key, value)
    elif re.match(r'^\$', value):
        # Assume the user as set the variable in value
        return spaces + 'PsiMod.set_option("%s", %s)' % (key, value[1:])
    elif re.match(yes, value):
        return spaces + 'PsiMod.set_option("%s", True)' % (key)
    elif re.match(no, value):
        return spaces + 'PsiMod.set_option("%s", False)' % (key)
    else:
        # Just place the value in quotes
        return spaces + 'PsiMod.set_option("%s", "%s")' % (key, value)

def process_set_commands(matchobj):
    spaces = matchobj.group(1)
    module = matchobj.group(2)
    commands = matchobj.group(3)
    command_lines = re.split('\n', commands)

    result = []
    for line in command_lines:
        temp = re.sub(r'#.*', "", line)
        result.append(re.sub(r'^\s*()(\w+)\s+(.*)($|#.*)', process_set_command, temp))

    if module != "":
        x = 'PsiMod.set_default_options_for_module("%s")' % (module.upper())
        result.insert(0, x)
        result.insert(1, 'PsiMod.set_option("NO_INPUT", True)')

    set_commands = spaces
    set_commands += (spaces).join(result)

    return set_commands

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
    extract += "\nPsiMod.set_active_molecule(%s)" % name
    extract += '\nPsiMod.IO.set_default_namespace("%s")' % name

    return extract

def process_print_command(matchobj):
    spaces = matchobj.group(1)
    string = matchobj.group(2)

    printer = "\npsi_string_print = str(%s)\n" % string
    printer += "PsiMod.print_out(psi_string_print)\n"

    return printer

def process_input(raw_input):
    # Process all "set name? { ... }"
    set_commands = re.compile(r'^(\s*?)set\s*(\w*?)\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(set_commands, process_set_commands, raw_input)

    # Process all individual "set key value"
    set_command = re.compile(r'(\s*?)set\s+(\w+)[\s=]+(.*?)($|#.*)', re.MULTILINE | re.IGNORECASE)
    temp = re.sub(set_command, process_set_command, temp)

    # Process all "global(s) { ... }"
    globals_command = re.compile(r'^(\s*?)(global|globals) \s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(globals_command, process_globals_command, temp)

    # Process all individual "global key value"
    global_command = re.compile(r'(\s*?)global\s+(\w+)[\s=]+(.*?)($|#.*)', re.MULTILINE | re.IGNORECASE)
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

    temp = "from PsiMod import *\nfrom psi import *\nfrom psiopt import *\n" + temp

    return temp

if __name__ == "__main__":
    result = process_input("""
molecule h2 {
H
H 1 R

R = .9
}

Rs = [ 0.9, 1.0, 1.1 ]

set scf {
    SCF_TYPE DF
    PRINT 2
    BASIS 3-21G
    RI_BASIS_SCF STO-3G
    DOCC [3, 0, 1, 1]
    DIIS on
}

globals {
    BASIS STO-3G
    RI_BASIS_SCF STO-3G
}

set docc [2, 0, 1, 1]

global print 1

for val in Rs:
    h2.R = val

    set basis cc-pv(T+d)Z

    for base in basissets:
        set basis $base

    set scf {
        SCF_TYPE DF
        PRINT 2
        BASIS 3-21G
        RI_BASIS_SCF STO-3G
        DOCC [3, 0, 1, 1]
    }

    escf = scf()

""")

    print "Result\n=========================="
    print result
