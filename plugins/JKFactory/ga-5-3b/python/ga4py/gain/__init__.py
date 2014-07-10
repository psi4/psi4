import mpi4py.MPI
from core import *
from misc import *
import random

if __name__ != '__main__':
    import inspect as _inspect
    import sys
    import numpy
    # imports from 'numpy' module every missing attribute into 'gain' module
    self_module = sys.modules[__name__]
    # import all classes not already overridden
    for name in dir(numpy):
        if not hasattr(self_module, name):
            attr = getattr(numpy, name)
            if _inspect.isclass(attr):
                setattr(self_module, name, attr)
    # import some other numpy functions directly
    from numpy import alen
    from numpy import newaxis

class PrintZero(object):
    def __init__(self):
        self.me = ga.nodeid()
        self.stdout = sys.stdout
    def write(self, something):
        if not self.me:
            self.stdout.write(something)
    def flush(self):
        if not self.me:
            self.stdout.flush()
import sys
_stdout = sys.stdout
sys.stdout = PrintZero()
