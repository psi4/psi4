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

Notes on Options
================

.. note:: The options referred to in the :ref:`sec:methods` section below
   and indexed in :ref:`apdx:options_c_module` are placed in ``set`` blocks as
   described in :ref:`sec:jobControl`, not as arguments to a Python function
   (like ``energy()``).

.. note:: All |PSIfour| keyword names and values are insensitive to case, both
   those that are placed in ``set`` blocks and as Python function arguments.
   The few exceptions are documented for the :py:func:`~psi4.driver.wrapper_database.database` function,
   where case structure must match the database file.

.. _`op_c_bool`:
.. _`op_c_boolean`:
.. note:: Boolean options can be specified by ``yes``, ``on``, ``true``, or ``1``
    for affirmative and ``no``, ``off``, ``false``, or ``0`` for negative,
    all insensitive to case.

.. _`op_c_conv`:
.. note:: Certain convergence and tolerance keywords, of type *double* (real numbers),
   may be specified using either a real number or an integer; and integer *X* is then
   treated as the number of converged decimal digits required. For example, to request
   as energy converged to :math:`10^{-6} E_h`, the user may set the ``e_convergence``
   keyword to ``0.000001``, ``1.0e-6``, or ``6``.

.. _`sec:psivarnotes`:

Notes on PSI Variables
======================

.. note:: Starting in 1.6, there are three standard ways to access an excited state
   property. We give examples below, but the method name and property name may change.

   * ``method ROOT 0 -> ROOT m property`` to get root ``m``.
   * ``method ROOT 0 -> ROOT m property - h TRANSITION`` to get root m and
     independently specify that the total transition symmetry is ``h``, e.g., A2.
   * ``method ROOT 0 (h) -> ROOT m (i) property`` to get the transition
     between two roots, specifying the symmetry of both states and the index of the target
     roots among states of their own symmetry.

   For example, to target the second excited-state, which is also the lowest energy state
   of its irrep, the first two calls will take m = 2, while the last takes m = 0.
   Methods that use this interface are: TD-fctl.
   Note that numberings are associated with the calculation much more strongly than 
   with the molecular system. Changing the number of roots sought, the symmetry 
   subspace or the symmetry apportionment of roots under which the computation is run, 
   or the excited state method are all likely to scramble root numberings.

