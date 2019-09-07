.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2019 The Psi4 Developers.
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

.. index:: postg
.. _`sec:postg`:

Interface to postg (XDM dispersion correction)
==============================================

.. codeauthor:: Alberto Otero de la Roza
.. sectionauthor:: Alberto Otero de la Roza

*Module:* :ref:`Samples <apdx:testSuitepostg>`

The `exchange-hole dipole moment (XDM) model <http://schooner.chem.dal.ca/wiki/XDM>`_ 
[Becke:2007:154108]_ is a dispersion correction used to get accurate intermolecular
interaction energies in density-functional theory calculations. The
XDM energy is calculated post-SCF using an external program called 
`postg <http://schooner.chem.dal.ca/wiki/Postg>`_. This document
describe the installation of postg and how to run calculations using
XDM-corrected functionals in |PSIfour|.

Installation
~~~~~~~~~~~~

Postg must be compiled from the source code, available on
`github <https://github.com/aoterodelaroza/postg>`_. Compilation
of postg is simple. Usually you just need to enter the postg
directory, select a compiler in the Makefile, and run `make`. Detailed
installation instructions for postg can be found 
`here <https://github.com/aoterodelaroza/postg/blob/master/README>`_
To use any XDM-corrected functional in |PSIfour|, the program binary
(``postg``) must be in your :envvar:`PSIPATH` or :envvar:`PATH` (in
that order).

Running XDM calculations
~~~~~~~~~~~~~~~~~~~~~~~~

Calculations for selected XDM-corrected functionals known to |PSIfour|
can be performed using the corresponding keywords. These functionals
work in the same way as any other: you can calculate single-point
energies, optimize geometries, calculate vibrational frequencies,
etc. The list of XDM-enabled functionals includes: BLYP-XDM, PBE-XDM, 
PW86PBE-XDM, B3LYP-XDM, B3PW91-XDM, PBE0-XDM, CAM-B3LYP-XDM,
BHANDHLYP-XDM, TPSS-XDM, HSE06-XDM, BP86-XDM, and B86bPBE-XDM. The
damping function parameters for these functionals are automatically
chosen by the program [OterodelaRoza:2013:054103]_.

For functionals that are not on the list, XDM can also be applied
manually using the `dft_xdm_*` variables to set the damping function
coefficients and the volume token. The variables are:

* `dft_xdm_a1` is the a1 damping function coefficient.

* `dft_xdm_a2` is the a2 damping function coefficient (in angstrom).

* `dft_xdm_vol` is the volume token. This can be either a functional
  keyword known by postg or a real value between 0 and 1 indicating
  the percentage of exact exchange in the functional.

The use of these keywords immediately activates the use of the XDM
correction in the code.

Some examples:

* A B86bPBE-XDM single-point calculation ::

   energy('b86bpbe-xdm')

* A B3LYP-XDM geometry optimization ::

   optimize('b3lyp-xdm')

* A B3LYP-XDM single-point calculation using damping function parameters
  appropriate for a `small basis set <http://schooner.chem.dal.ca/wiki/XDM>`_ ::

   set basis 6-31+G*
   set dft_xdm_a1 0.4515
   set dft_xdm_a2 2.1357
   set dft_xdm_vol b3lyp
   energy('b3lyp')

