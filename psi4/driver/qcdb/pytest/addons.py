#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2024 The Psi4 Developers.
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

import pytest


def _plugin_import(plug):
    import sys
    if sys.version_info >= (3, 4):
        from importlib import util
        plug_spec = util.find_spec(plug)
    else:
        import pkgutil
        plug_spec = pkgutil.find_loader(plug)
    if plug_spec is None:
        return False
    else:
        return True


def is_psi4_new_enough(version_feature_introduced):
    if not _plugin_import('psi4'):
        return False
    import psi4
    from qcelemental.util import parse_version
    return parse_version(psi4.__version__) >= parse_version(version_feature_introduced)


#def is_numpy_new_enough(version_feature_introduced):
#    if not _plugin_import('numpy'):
#        return False
#    import numpy
#    from pkg_resources import parse_version
#    return parse_version(numpy.version.version) >= parse_version(version_feature_introduced)
#
#
#using_scipy = pytest.mark.skipif(_plugin_import('scipy') is False,
#                                reason='Not detecting module scipy. Install package if necessary and add to envvar PYTHONPATH')

using_psi4 = pytest.mark.skipif(_plugin_import('psi4') is False,
                                 reason='Not detecting module psi4. Install package and add to envvar PYTHONPATH')

#using_psi4_libxc = pytest.mark.skipif(is_psi4_new_enough("1.2a1.dev100") is False,
#                                reason="Psi4 does not include DFT rewrite to use Libxc. Update to development head")
#
#using_psi4_efpmints = pytest.mark.skipif(is_psi4_new_enough("1.2a1.dev507") is False,
#                                reason="Psi4 does not include EFP integrals in mints. Update to development head")
#
#using_psi4_python_integral_deriv = pytest.mark.skipif(is_psi4_new_enough("1000") is False,
#                                reason="Psi4 does not include derivatives of integrals exported to python. Update to development head")

using_psi4_molrec = pytest.mark.skipif(is_psi4_new_enough("1.2a1.dev999") is False,
                                reason="Psi4 does not use the new Molecule parsing. Update to development head")

#using_numpy_113 = pytest.mark.skipif(is_numpy_new_enough("1.13.0") is False,
#                                reason='NumPy does not include 1.13 features. Update package and add to envvar PYTHONPATH')
#
#using_matplotlib = pytest.mark.skipif(_plugin_import('matplotlib') is False,
#                                reason='Note detecting module matplotlib. Install package if necessary and add to envvar PYTHONPATH')

using_pylibefp = pytest.mark.skipif(_plugin_import('pylibefp') is False,
                                reason='Not detecting module pylibefp. Install package if necessary and add to envvar PYTHONPATH')

