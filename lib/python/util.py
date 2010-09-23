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
