import sys
import os

# Modify this to the required path
if "PSICORE" in os.environ.keys():
    sys.path.insert(1, os.environ["PSICORE"])
    if not os.path.isfile(os.path.join(os.environ["PSICORE"], 'psi4core.so')):
        raise ImportError("No psi4core.so found in folder %s" % os.environ["PSICORE"])
    print("Found PSICORE=%s, attempting import." % os.environ["PSICORE"])
    import psi4core
else:
    try:
        import psi4core
    except ImportError:
        psi_path = os.path.abspath(__file__ + '/../../objdir/stage/usr/local/lib/')
        print("psi4core.so not found in local folder, attempting to guess relative location %s" % psi_path)
        sys.path.insert(1, psi_path)
        try:
            import psi4core
        except ImportError:
            raise ImportError("Could not find psi4core.so at %s" % psi_path)
        


# Init psi4core
psi4core.initialize()
psi4core.set_memory(int(512e6)) # Set to 512 MB

# Cleanup psi4core at exit
import atexit
atexit.register(psi4core.set_legacy_molecule, None)
atexit.register(psi4core.finalize)
from psi4core import *
