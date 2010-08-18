import re;

def process_set_command(matchobj):
    spaces = matchobj.group(1)
    key   = matchobj.group(2).upper()
    value = matchobj.group(3).strip()

    if re.match(r'^[-]?\d*\.?\d+$', value) or re.match(r'^\[.*\]$', value):
        return spaces + 'PsiMod.set_option("%s", %s)' % (key, value)
    else:
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
    set_commands += (spaces+"\n").join(result)

    return set_commands

def process_molecule_command(matchobj):
    spaces = matchobj.group(1)
    name = matchobj.group(2)
    geometry = matchobj.group(3)
    molecule = spaces
    if name != "":
        molecule += '%s = ' % (name)

    molecule += 'geometry("""%s""")' % (geometry)

    return molecule

def process_input(raw_input):
    # Process all "set name? { ... }"
    set_commands = re.compile(r'^(\s*?)set\s*(\w*?)\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(set_commands, process_set_commands, raw_input)

    # Process all individual "set key value"
    set_command = re.compile(r'(\s*?)set\s+(\w+)[\s=]+(.*?)($|#.*)', re.MULTILINE | re.IGNORECASE)
    temp = re.sub(set_command, process_set_command, temp)

    # Process "molecule name? { ... }"
    molecule = re.compile(r'^(\s*?)molecule\s*(\w*?)\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(molecule, process_molecule_command, temp)

    temp = "from PsiMod import *\nfrom psi import *\n" + temp

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
}

set docc [2, 0, 1, 1]

for val in Rs:
    h2.R = val

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
