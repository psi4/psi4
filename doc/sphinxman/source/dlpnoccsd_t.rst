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

.. index::
   single: DLPNO-CCSD(T)

.. _`sec:dlpnoccsd_t`:

DLPNO-CCSD(T): Domain-Based Local Pair Natural Orbital CCSD(T)
==============================================================

.. codeauthor:: Andy Jiang 
.. sectionauthor:: Andy Jiang

*Module:* :ref:`Keywords <apdx:dlpno>`, :ref:`PSI Variables <apdx:dlpno_psivar>`, :source:`DLPNOCCSD(T) <psi4/src/psi4/dlpno>`

Introduction
------------

The gold-standard CCSD(T) method in quantum chemistry provides chemically reliable 
answers in regards to relative energies (< 1 kcal/mol deviation) with respect to
experimental results or a Full CI computation.

An example input file for a DLPNO-CCSD computation is::

   molecule h2o {
   0 1
   O
   H 1 1.0
   H 1 1.0 2 104.5
   symmetry c1
   }
   
   set basis cc-pvdz
   set scf_type df
   set freeze_core True
   set pno_convergence tight
   
   energy('dlpno-ccsd')

An example input file for a DLPNO-CCSD(T) computation is::

   molecule h2o {
   0 1
   O
   H 1 1.0
   H 1 1.0 2 104.5
   symmetry c1
   }
   
   set basis cc-pvdz
   set scf_type df
   set freeze_core True
   set pno_convergence tight
   
   energy('dlpno-ccsd(t)') # dlpno-ccsd(t0) for the semicanonical (T0) computation

Theory
------

Recommendations
---------------

