.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2017 The Psi4 Developers.
.. #
.. # The copyrights for code used from other parties are included in
.. # the corresponding files.
.. #
.. # This program is free software; you can redistribute it and/or modify
.. # it under the terms of the GNU General Public License as published by
.. # the Free Software Foundation; either version 2 of the License, or
.. # (at your option) any later version.
.. #
.. # This program is distributed in the hope that it will be useful,
.. # but WITHOUT ANY WARRANTY; without even the implied warranty of
.. # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. # GNU General Public License for more details.
.. #
.. # You should have received a copy of the GNU General Public License along
.. # with this program; if not, write to the Free Software Foundation, Inc.,
.. # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
.. #
.. # @END LICENSE
.. #

.. include:: autodoc_abbr_options_c.rst

.. _`sec:addAddOns`:

Adding Add-Ons
==============

* Project Name

  * Select a name. May be mixed case with numerals and underscores
    (*e.g.*, CheMPS2, libefp, PCMSolver, v2rdm_casscf). Shouldn't start with a
    numeral. Needn't start with "lib", even if a library.

  * GitHub repository name should be :samp:`{addon_name}` or
    :samp:`{addon_name}.lower()`. For example: CheMPS2, libefp, pcmsolver,
    v2rdm_casscf.

  * CMake project name should be :samp:`{addon_name}`. For example:
    ``project(libefp)``, ``project(CheMPS2)``, ``project(PCMSolver)``,
    ``project(v2rdm_casscf)``. Namespacing in the directory structure used
    to detect the addon should have this name (*e.g.*,
    ``share/cmake/CheMPS2``).

  * Restricted by the CMake project name, add-ons return CMake variables
    and compile definitions of :samp:`FOUND_{addon_name}` and
    :samp:`USING_{addon_name}`. For example: ``FOUND_libefp``,
    ``USING_CheMPS2``, ``PCMSolver_LIBRARIES``, ``USING_v2rdm_casscf``.

  * The CMake target(s) formed use the full add-on name as the namespace,
    :samp:`{addon_name}::{lib_name_without_lib}.lower()`. For example:
    ``libefp::efp``, ``CheMPS2::chemps2``, ``PCMSolver::pcm``,
    ``v2rdm_casscf::v2rdm_casscf``.

  * Following the CMake project name (though not restricted to it --
    |PSIfour| managment could change the pattern), the user flag to enable
    an add-on is :samp:`ENABLE_{addon_name}`. Note that binary-only
    add-ons don't go through this enabling process.

  * Internally, the ExternalProject_Add and dummy libraries as well as any
    tests/ and external/ subdirectories should all be lowercase,
    :samp:`{addon_name}.lower()`.

  * The `conda package <https://anaconda.org/psi4/repo>`_ and internal to
    |PSIfour| (that is, the ExternalProject_Add, dummy libraries, and any
    tests/ and external/ subdirectories) should all be lowercase,
    :samp:`{addon_name}.lower()`.

  * Alternatively, you can do everything mentioned here lowercase and just
    have a different capitalization for an advertising name. After all,
    that's what |PSIfour| does.



* Debugging

  >>> lldb -- python stage/CMAKE_INSTALL_PREFIX/bin/psi4 ../tests/tu1-h2o-energy/input.dat -l stage/CMAKE_INSTALL_PREFIX/share/psi4/
  >>> (lldb) run


