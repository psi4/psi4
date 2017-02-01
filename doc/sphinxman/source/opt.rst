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
   triple: setting; keywords; optimize()
   pair: geometry optimization; function call
   pair: gradient; function call
   see: opt(); optimize()

.. _`sec:opt()`:

Geometry Optimization |w---w| :py:func:`~psi4.optimize` and :py:func:`~psi4.gradient`
=====================================================================================

* :ref:`Psi4 Native Gradient Methods <table:grad_gen>`
* :ref:`Psi4 Native DFT Gradient Methods (excepting double-hybrids) <table:grad_gen>`
* :ref:`CFOUR Interfaced Gradient Methods <table:grad_cfour>`

For further discussion of geometry optimization, see
Sec. :ref:`sec:optking`.

:py:func:`~psi4.optimize` is the only command most users will ever
need to access directly to perform geometry optimizations. Behind
the scenes, :py:func:`~psi4.optimize` is a wrapper that repeatedly
calls :py:func:`~psi4.gradient` that computes the gradient then adds a
call to the :ref:`geometry optimization module <sec:optking>`.

.. autofunction:: psi4.optimize(name [, molecule, return_wfn, func, mode, dertype, hessian_with])

.. autofunction:: psi4.gradient(name [, molecule, return_wfn, func, dertype])

