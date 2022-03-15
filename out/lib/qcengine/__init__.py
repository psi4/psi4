"""
Base file for the dqm_compute module.
"""

from . import config, exceptions
from .compute import compute, compute_procedure
from .config import get_config
from .extras import get_information
from .mdi_server import MDIServer
from .procedures import get_procedure, list_all_procedures, list_available_procedures
from .programs import get_program, list_all_programs, list_available_programs, register_program, unregister_program
from .stock_mols import get_molecule

# Handle versioneer
__version__ = get_information("version")
__git_revision__ = get_information("git_revision")
del get_information
