.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2025 The Psi4 Developers.
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
   single: harmonic vibrational analysis, vibrational analysis, thermochemical analysis

.. _`sec:thermo`:

Vibrational and Thermochemical Analysis
=======================================

.. codeauthor:: Rollin A. King and Lori A. Burns
.. sectionauthor:: Lori A. Burns

*Module:* :ref:`Keywords <apdx:thermo>`, :ref:`PSI Variables <apdx:thermo_psivar>`, :source:`THERMO <psi4/driver/qcdb/vib.py>`


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

A thermochemical analysis is performed after any full (not just specific
symmetry subgroups). If the wavefunction is retained, it may be reused
at a different temperature, pressure, rotational symmetry number, or
isotopic substitution through the function :py:func:`psi4.driver.qcdb.vib.thermo`
as is shown in :srcsample:`freq-isotope2`.

A few summary psivars are set: "ZPVE", "THERMAL ENERGY CORRECTION",
"ENTHALPY CORRECTION", "GIBBS FREE ENERGY CORRECTION", "ZERO K
ENTHALPHY", "THERMAL ENERGY", "ENTHALPY", "GIBBS FREE ENERGY".
But additionally, every valid combination of {S, Cv, Cp, ZPE, E, H, G}
with {elec, trans, rot, vib, corr, tot} (e.g., vibrational entropy,
S_vib, and enthalpy correction, H_corr) is returned by dictionary
from the ``thermo`` function. See :source:`tests/pytests/test_vibanalysis.py`
(near the end) for an example.


.. index::
   pair: vibrational analysis; output

Output
^^^^^^

The full list of keywords for thermo is provided in Appendix :ref:`apdx:thermo`.

Information on the Psithon function that drives frequency analyses is provided
at :py:func:`~psi4.driver.frequency`.

