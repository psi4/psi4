#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2022 The Psi4 Developers.
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

import pickle

from . import dependency_check

# printing and logging formatting niceties
import pprint
from functools import partial
import numpy as np
pp = pprint.PrettyPrinter(width=120, compact=True, indent=1)
nppp = partial(np.array_str, max_line_width=120, precision=8, suppress_small=True)
nppp10 = partial(np.array_str, max_line_width=120, precision=10, suppress_small=True)
del partial

from qcelemental import constants
from psi4.driver import psifiles as psif

from psi4.driver.ipi_broker import ipi_broker
from psi4.driver.molutil import *
from psi4.driver.inputparser import process_input
from psi4.driver.p4util.util import *
from psi4.driver.p4util.testing import *
from psi4.driver.p4util.fcidump import *
from psi4.driver.p4util.fchk import *
from psi4.driver.p4util.text import *
from psi4.driver.qmmm import QMMM, QMMMbohr
from psi4.driver.pluginutil import *

from psi4.driver import gaussian_n
from psi4.driver import aliases
from psi4.driver import diatomic
from psi4.driver import wrapper_database
from psi4.driver import wrapper_autofrag
from psi4.driver import schema_wrapper
from psi4.driver import schema_wrapper as json_wrapper # Deprecate in 1.4
from psi4.driver import frac

from psi4.driver.driver import *

# Single functions
from psi4.driver.driver_cbs import cbs  # remove in v1.8 when UpgradeHelper expires
from psi4.driver.p4util.python_helpers import set_options, set_module_options, pcm_helper, basis_helper
