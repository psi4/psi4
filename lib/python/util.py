import PsiMod
import sys
from psiexceptions import *

def set_memory(bytes):
    PsiMod.set_memory(bytes)

def get_memory():
    return PsiMod.get_memory()

def set_num_threads(nthread):
    PsiMod.set_n_threads(nthread)

def get_num_threads():
    return PsiMol.get_n_threads()

# Test functions
def compare_values(expected, computed, digits, label):
    if (abs(expected-computed) > 10**(-digits)):
        print "\t%s: computed value (%f) does not match (%f) to %d digits." % (label, expected, computed, digits)
        sys.exit(1)

    print "\t%s: matched." % (label)

def compare_strings(expected, computed, label):
    if(expected != computed):
        print "\t%s: computed value (%f) does not match (%f)." % (label, expected, computed)
        sys.exit(1)

    print "\t%s: matched." % (label)
