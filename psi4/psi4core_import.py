import sys
import os

try:
    from . import psi4core
except ImportError:
    print("Psi4 is not installed, looking for the 'objdir' build directory for ps4icore.so ...")
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

    print("Found psi4core.so at %s" % matches[0])
    sys.path.insert(1, matches[0])
    import psi4core

# Init psi4core
psi4core.initialize()

# Set psidatadir
if "PSIDATADIR" not in os.environ.keys():
    datadir = os.path.dirname(os.path.abspath(__file__))
    datadir += os.path.sep + "share" + os.path.sep + "psi4"
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
atexit.register(psi4core.clean)
atexit.register(psi4core.finalize)

# Numpy place holder for files and cleanup
numpy_files = []
def clean_numpy_files():
    for nfile in numpy_files:
        os.unlink(nfile)

atexit.register(clean_numpy_files)