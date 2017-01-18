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
   triple: setting; keywords; frequency()
   pair: vibrational analysis; function call
   pair: hessian; function call
   see: freq(); frequency();
   see: frequencies(); frequency();

.. _`sec:freq()`:

Harmonic Vibrational Analysis |w---w| :py:func:`~psi4.frequency` and :py:func:`~psi4.hessian`
=============================================================================================

* :ref:`Psi4 Native Hessian Methods <table:freq_gen>`

For further discussion of vibrational and thermochemical analysis,
see Sec. :ref:`sec:thermo`.

:py:func:`~psi4.frequency` is the only command most users will ever
need to access directly to perform frequency calculations. Behind
the scenes, :py:func:`~psi4.frequency` is a light wrapper over
:py:func:`~psi4.hessian` that computes the Hessian then adds a
thermochemical analysis.

.. autofunction:: psi4.frequency(name [, molecule, return_wfn, func, mode, dertype, irrep])

.. autofunction:: psi4.hessian(name [, molecule, return_wfn, func, dertype, irrep])

