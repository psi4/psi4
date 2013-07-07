
import os
import subprocess

from grendel.util.descriptors import RaiseOnAccessDescriptor
from grendel.util.overloading import get_kwarg
from grendel.util import SystemInfo

# OS check... currently only 'posix' is supported
if not os.name == 'posix':
    raise OSError("Q-Chem interface is currently only available on posix-like systems. "
                  " This system returns '" + os.name + "' for the Python expression 'os.name'")

SystemInfo.qchem_executable = get_kwarg(os.environ, 'QCHEM_EXECUTABLE', 'QCHEM_EXE', 'QCHEM')
if not SystemInfo.qchem_executable:
    try:
        # try to find Q-Chem using `which`, which should be always available on posix systems
        # the last part strips the newline from the output
        SystemInfo.qchem_executable = subprocess.check_output(["which", 'qchem'])[:-1]
        if not os.access(SystemInfo.qchem_executable, os.X_OK):
            SystemInfo.qchem_executable = RaiseOnAccessDescriptor(SystemError,
                'Q-Chem executable at {0} is not an executable file.  Check the contents of the '
                'QCHEM_EXECUTABLE environment variable (or, if it is not defined, check the '
                'QCHEM_EXE or QCHEM environment variables.)'.format(SystemInfo.qchem_executable))
    except subprocess.CalledProcessError:
        SystemInfo.qchem_executable = RaiseOnAccessDescriptor(SystemError,
            'Q-Chem executable not found.  Try a different interface or set the QCHEM_EXECUTABLE '
            'environment variable to the location of the executable.  (Equivalently, you can set'
            ' either of the QCHEM or QCHEM_EXE environment variables.')

# Check that the Q-Chem executable is in fact executable:

from grendel.interface.qchem.input_generator import *
from grendel.interface.qchem.runner import *
from grendel.interface.qchem.output_parser import *

