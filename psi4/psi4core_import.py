import sys
import os

try:
    from . import psi4core
except ImportError:
    print("psi4core.so not found in local folder. Psi4 is not installed, looking for the 'objdir' build directory...")
    base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    base_path += os.path.sep + 'objdir' + os.path.sep + 'stage'
    matches = []
    for root, dirnames, filenames in os.walk(base_path):
        if 'include' in root: continue
        if 'share' in root: continue
        if 'psi4core.so' in filenames:
            matches.append(root)

    if len(matches) == 0:
        raise ImportError("Could not find psi4core.so in basepath: %s" % base_path)

    sys.path.insert(1, matches[0])
    import psi4core


# Init psi4core
psi4core.initialize()
psi4core.set_memory(int(512e6)) # Set to 512 MB

# Set psidatadir
if "PSIDATADIR" not in os.environ.keys():
    datadir = os.path.dirname(os.path.abspath(__file__))
    datadir += os.path.sep + "share" + os.path.sep + "psi4"
    print datadir
    os.environ["PSIDATADIR"] = datadir
else:
    datadir = os.environ["PSIDATADIR"]

if not os.path.isdir(datadir):
     raise KeyError("Unable to read the Psi4 Python folder - check the PSIDATADIR environmental variable"
                    "      Current value of PSIDATADIR is %s" % datadir)

psi4core.set_environment("PSIDATADIR", datadir)


# Cleanup psi4core at exit
import atexit
atexit.register(psi4core.set_legacy_molecule, None)
atexit.register(psi4core.finalize)
atexit.register(psi4core.clean)

# Move up the namesapce
#from psi4core import *
