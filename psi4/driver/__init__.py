#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2021 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#
"""
Python driver for Psi4

isort:skip_file
"""

from qcelemental import constants

from . import dependency_check
from .inputparser import process_input
from . import psifiles as psif

from . import pluginutil
from .pluginutil import *

from . import ipi_broker
from .ipi_broker import *

from . import molutil
from .molutil import *

from . import wrapper_autofrag
from .wrapper_autofrag import *

from . import p4util
from .p4util import *

from . import mdi_engine
from .mdi_engine import *

from . import driver_findif
from .driver_findif import *

from . import procrouting
from .procrouting import *
from .procrouting.proc import scf_helper  # plugins use psi4.driver.scf_helper

from . import driver_util
from .driver_util import *

from . import driver_cbs
from .driver_cbs import *

from . import qmmm
from .qmmm import *

from . import driver_nbody_helper
from .driver_nbody_helper import *

from . import driver_nbody
from .driver_nbody import *

from . import driver
from .driver import *

from . import wrapper_database
from .wrapper_database import *

from . import prop_util
from .prop_util import *

from . import frac
from .frac import *

from . import aliases
from .aliases import *

from . import gaussian_n
from .gaussian_n import *

from . import diatomic
from .diatomic import *

from . import schema_wrapper
from . import schema_wrapper as json_wrapper # Deprecate in 1.4
from .schema_wrapper import *


__all__ = [
    "ConvergenceError",
    "DB_RGT",
    "DB_RXN",
    "ManagedMethodError",
    "MissingMethodError",
    "OptimizationConvergenceError",
    "QMMM",
    "SCFConvergenceError",
    "TDSCFConvergenceError",
    "TestComparisonError",
    "UpgradeHelper",
    "ValidationError",
    "activate",
    "allen_focal_point",
    "anharmonicity",
    "auto_fragments",
    "banner",
    "basis_helper",
    "cbs",
    "compare",
    "compare_arrays",
    "compare_cubes",
    "compare_fchkfiles",
    "compare_fcidumps",
    "compare_integers",
    "compare_matrices",
    "compare_molrecs",
    "compare_recursive",
    "compare_strings",
    "compare_values",
    "compare_vectors",
    "compare_wavefunctions",
    "constants",
    "corl_xtpl_helgaker_2",
    "cubeprop",
    "database",
    "db",
    "diatomic",  # module  # used in tests but should encourage anharmonic instead
    "driver_cbs",  # module  # used in tests but should encourage xtpl instead
    "energies_from_fcidump",
    "energy",
    "fake_file11",
    "fchk",
    "fcidump",
    "fcidump_from_file",
    "frac_nuke",
    "frac_traverse",
    "freq",
    "frequencies",
    "frequency",
    "gdma",
    "geometry",
    "get_memory",
    "gradient",
    "hessian",
    "ip_fitting",
    "ipi_broker",
    "json_wrapper",  # module  # kill off in 1.5
    "mdi_init",
    "mdi_run",
    "MDIEngine",
    "message_box",
    "molden",
    "oeprop",
    "opt",
    "optimize",
    "optimize_geometric",
    "p4util",  # module  # used in tests but psi4.driver.p4util preferred
    "pcm_helper",
    "procedures",  # old plugins using this as `psi4.procedures`. current plugins use `psi4.driver.procedures` so not needed in all
    "prop",
    "properties",
    "qcdb",  # module
    "scf_xtpl_helgaker_2",
    "scf_xtpl_helgaker_3",
    "scf_xtpl_karton_2",
    "scf_xtpl_truhlar_2",
    "schema_wrapper",  # module
    "set_memory",
    "set_module_options",  # kill off in 1.5
    "set_options",
    "sherrill_gold_standard",
    "tdscf",
    "vibanal_wfn",
    "xtpl_highest_1",
]
