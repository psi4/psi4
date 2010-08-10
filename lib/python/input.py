#import PsiMod;
#import psi;
import re;

def process_set_command(matchobj):
    key   = matchobj.group(1).upper()
    value = matchobj.group(2).strip()

    if re.match(r'^[-]?\d*\.?\d+$', value) or re.match(r'^\[.*\]$', value):
        return 'PsiMod.set_option("%s", %s)' % (key, value)
    else:
        return 'PsiMod.set_option("%s", "%s")' % (key, value)

def process_set_commands(matchobj):
    commands = matchobj.group(2)
    command_lines = re.split('\n', commands)

    result = []
    for line in command_lines:
        temp = re.sub(r'#.*', "", line)
        result.append(re.sub(r'^\s*(\w+)\s+(.*)($|#.*)', process_set_command, temp))

#        result.append(re.sub(r'^\s*(\w+)\s+([\w\.-]+)', process_set_command, line))

    set_commands = ""
    if matchobj.group(1) != "":
        set_commands += 'PsiMod.set_default_options_for_module("%s")\nPsiMod.set_option("NO_INPUT", True)\n' % (matchobj.group(1).upper())

    set_commands += "\n".join(result)

    return set_commands

def process_molecule_command(matchobj):
    molecule = ""
    if matchobj.group(1) != "":
        molecule += '%s = ' % (matchobj.group(1))

    molecule += 'geometry("""%s""")' % (matchobj.group(2))

    return molecule

def process_input(raw_input):
    # Process all "set name? { ... }"
    set_commands = re.compile(r'set\s*(\w*?)\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(set_commands, process_set_commands, raw_input)

    # Process all individual "set key value"
    temp = re.sub(r'set\s*(\w+)\s+([\w\.-]+)', process_set_command, temp)

    # Process "molecule name? { ... }"
    molecule = re.compile(r'molecule\s*(\w*?)\s*\{(.*?)\}', re.MULTILINE | re.DOTALL | re.IGNORECASE)
    temp = re.sub(molecule, process_molecule_command, temp)

    temp = "from PsiMod import *\nfrom psi import *\n" + temp

    return temp

if __name__ == "__main__":
    result = process_input("""
# Run: psi4 -s

molecule {
O 0.0 0.0 0.0
H 0.0 0.0 1.0
H 0.0 1.0 0.0
}

set scf {
#   GUESS   SAD
    SCF_TYPE DF
    PRINT 2
    BASIS 3-21G
    RI_BASIS_SCF STO-3G
    DOCC [3, 0, 1, 1]
}

scf();

set dfmp2 {
    BASIS  3-21G
    RI_BASIS_MP2 STO-3G
}

dfmp2();
""")

    print "Result\n=========================="
    print result
