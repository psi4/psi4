.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2024 The Psi4 Developers.
.. #
.. # The copyrights for code used from other parties are included in
.. # the corresponding files.
.. #
.. # This file is part of Psi4.
.. #
.. # Psi4 is free software; you can redistribute it and/or modify
.. # it under the terms of the GNU Lesser General Public License as published by
.. # the Free Software Foundation, version 3.
.. #
.. # Psi4 is distributed in the hope that it will be useful,
.. # but WITHOUT ANY WARRANTY; without even the implied warranty of
.. # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. # GNU Lesser General Public License for more details.
.. #
.. # You should have received a copy of the GNU Lesser General Public License along
.. # with Psi4; if not, write to the Free Software Foundation, Inc.,
.. # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
.. #
.. # @END LICENSE
.. #

.. include:: autodoc_abbr_options_c.rst

.. _`sec:prog_tour`:

==========================================================
General layout of the core: where new C++ code should live
==========================================================

New integral operators
----------------------

Most of the heavy lifting in |PSIfour| is handled by libmints, which can be
found in the :source:`psi4/src/psi4/libmints` directory.  New types of one- and
two-electron integral operators should be added here.  The Wavefunction class
:source:`psi4/src/psi4/libmints/wavefunction.h` is also found here and is a key
part of the infrastructure.  Every energy calculation is performed by a class
that derives from Wavefunction and is accessible as a return value in the
Python layer.  The Wavefunction class contains all pertinent calculation
results, such as one-particle densities, molecular orbitals and gradients.

Completely new methods
----------------------

A new method that is not a modification of existing code belongs in its own
folder in :source:`psi4/src/psi4`; see other folders in that exist in that
location for examples of setting up CMake, and make sure that the new folder is
added to :source:`psi4/src/psi4/CMakeLists.txt`.  There are also a number of
variables that can be exported to be available to the user, as detailed in
:ref:`sec:psiVariables`.  To set these variables, the following member of
Wavefunction should be called::

    set_variable("Variable Name", variable_value);

The new variable should also be documented in
:source:`doc/sphinxman/source/glossary_psivariables.rst`.  There are a
number of different helpers to export various quantities from the wavefunction
to external formats such as FCHK and MOLDEN.  Because the Wavefunction makes
its members available to the Python layer, any other similar export functions
should be written in python.

Integral consuming technologies
-------------------------------

The general philosophy in |PSIfour| is to try write two-electron integral
driven tasks in methods like Hartree-Fock, CIS and CPHF in terms of generalized
Fock-like matrices.  From here, a single class can be used to construct these
generalized Fock matrices, which is what libFock
(:source:`psi4/src/psi4/libfock`) accomplishes.  A number of integral
technologies -- such as integral-direct, disk-based and density fitting -- are
supported in libFock, making them generally available to all elements of the
code that use the generalized Fock matrix strategy.

.. _`sec:prog_tour-exposing`:

Exposing C++ code to Python
---------------------------

The recent push to move sections of the code that are not a bottleneck into the
Python layer requires that the C++ code is callable from Python and that its
results are accessible.  The result accessibility is addressed by populating
the appropriate variables in the Wavefunction object.  To make the code
callable from Python, we rely on the excellent `PyBind11 <https://pybind11.readthedocs.io/en/stable/>`_ library
to create the bindings.  Existing code to export various |PSIfour| classes can
be found in :source:`psi4/src` in the files whose name begins with `export_`.
The code to export functions that run entire calculations is usually found in
:source:`psi4/src/core.cc`.
