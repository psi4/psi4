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

.. index::
   single: harmonic vibrational analysis, vibrational analysis, thermochemical analysis

.. _`sec:thermo`:

Vibrational and Thermochemical Analysis
=======================================

.. codeauthor:: Rollin A. King
.. comment.. sectionauthor:: Rollin A. King and Lori A. Burns

*Module:* :ref:`Keywords <apdx:thermo>`, :ref:`PSI Variables <apdx:thermo_psivar>`, :source:`THERMO <psi4/src/psi4/thermo>`


.. caution:: It is important to know that |PSIfour|, like any other
   quantum chemistry program, does *not* compute the usual enthalpies,
   entropies, or Gibbs free energies *of formation* provided by most
   reference books.  Instead, quantum chemistry programs compute "absolute"
   thermodynamic properties relative to infinitely separated nuclei and
   electrons, not "formation" values relative to elements in their standard
   states.  If you are computing thermodynamic differences, like a reaction
   enthalpy computed as the enthalpy of the products minus the enthalpy
   of the reactants, then these "absolute" enthalpies are perfectly valid
   and usable.  However, they cannot be mixed and matched with enthalpies of
   formation from reference books, since the zero of energy is not the same.
   Additionally, the "thermal energies" reported in kcal/mol are the 
   finite-temperature *corrections* to the electronic total energy, and 
   not the overall thermal energies themselves.  If in doubt, use the
   reported Total Energies in Hartree/particle.

Keywords
^^^^^^^^

.. include:: autodir_options_c/thermo__t.rst
.. include:: autodir_options_c/thermo__p.rst
.. include:: autodir_options_c/thermo__rotational_symmetry_number.rst

.. _`sec:thermoExamples`:

Examples
^^^^^^^^


.. index::
   pair: vibrational analysis; output

Output
^^^^^^

The full list of keywords for thermo is provided in Appendix :ref:`apdx:thermo`.

Information on the Psithon function that drives frequency analyses is provided
at :py:func:`~psi4.frequency`.

