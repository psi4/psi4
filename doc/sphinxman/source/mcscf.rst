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

.. index:: MCSCF

.. index::
   pair: MCSCF; theory

.. index::
   pair: CI; multi-configurational self-consistent-field

.. _`sec:mcscf`:

MCSCF: Multi-Configurational Self-Consistent-Field
==================================================

.. codeauthor:: Daniel G. A. Smith, C. David Sherrill, and Matthew L. Leininger
.. sectionauthor:: Daniel G. A. Smith and C. David Sherrill

*Module:* :ref:`Keywords <apdx:detci>`, :ref:`PSI Variables <apdx:detci_psivar>`, :source:`DETCI <psi4/src/psi4/detci>`

As the cost of Full CI scales exponentially with respect to the number of
active orbitals it is often advantageous to neglect orbitals that do not
exhibit strong correlation. These orbitals are variationally optimized
simultaneously with the CI coefficients and known as Multi-Configurational
Self-Consistent Field (MCSCF). The most commonly used MCSCF procedure is the
complete-active-space self-consistent-field (CASSCF) approach [Roos:1980]_,
which includes all possible determinants (with the proper symmetry) that can be
formed by distributing a set of active electrons among a set of active
orbitals. The MCSCF module performs CASSCF optimization of molecular orbitals
via a two-step procedure in which the CI coefficients and orbitals are
optimized in an alternating manner. The program uses a fairly simple
approximate orbital Hessian [Chaban:1997:88]_ and a Newton-Raphson update,
accelerated by Pulay's DIIS procedure [Pulay:1980]_. We have also implemented
the RASSCF method [Malmqvist:1990:RASSCF]_, which is another kind of MCSCF
which is typically less complete (and less expensive) than CASSCF.

Inactive orbitals in the MCSCF may be specified by the
|globals__restricted_docc| and |globals__restricted_uocc| keywords. These
orbitals will remain doubly-occupied or doubly-unoccupied, respectively, in the
MCSCF wavefunction.  However, the form of these orbitals will be optimized in
the MCSCF procedure.  It is also possible to literally freeze inactive orbitals
in their original (SCF) form using the |globals__frozen_docc| and
|globals__frozen_uocc| keywords.  This is not normally what one wishes to do in
an MCSCF computation (*e.g.*, it complicates the computation of gradients), but
it can make the computations faster and is helpful in some circumstances where
unphysical mixing of inactive and active occupied orbitals might occur.
Presently, it is not possible to mix the use of restricted and frozen orbitals
in |PSIfour|.

An illustrative CASSCF example is as follows::

    molecule {
    O
    H 1 1.00
    H 1 1.00 2 103.1
    }
    
    set {
        basis           6-31G**
        restricted_docc [1, 0, 0, 0]
        active          [3, 0, 1, 2]
    }
    energy('casscf')

This input will compute the CASSCF energy of water where the 1s Oxygen orbital
and several virtual orbitals are not included in the CI expansion, but are
still optimized. The following is a full list of spaces within the various MCSCF
types.

.. _`table:mcscf_spaces`:

.. table:: Orbital spaces for MCSCF computations

    +----------------------------+----------------------------+
    | RASSCF                     | CASSCF                     |
    +============================+============================+
    | |globals__frozen_uocc|     | |globals__frozen_uocc|     |
    +----------------------------+----------------------------+
    | |globals__restricted_uocc| | |globals__restricted_uocc| |
    +----------------------------+----------------------------+
    | |globals__ras4|            | |globals__active|          |
    +----------------------------+                            +
    | |globals__ras3|            |                            |
    +----------------------------+                            +
    | |globals__ras2|            |                            |
    +----------------------------+                            +
    | |globals__ras1|            |                            |
    +----------------------------+----------------------------+
    | |globals__restricted_docc| | |globals__restricted_docc| |
    +----------------------------+----------------------------+
    | |globals__frozen_docc|     | |globals__frozen_docc|     |
    +----------------------------+----------------------------+

Basic MCSCF Keywords
~~~~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/detci__mcscf_e_convergence.rst
.. include:: autodir_options_c/detci__mcscf_r_convergence.rst
.. include:: autodir_options_c/detci__mcscf_type.rst
.. include:: autodir_options_c/detci__mcscf_algorithm.rst
.. include:: autodir_options_c/detci__mcscf_maxiter.rst
.. include:: autodir_options_c/detci__mcscf_rotate.rst
.. include:: autodir_options_c/detci__mcscf_diis_start.rst


