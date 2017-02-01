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

Notes on Options
================

.. comment warning:: Python naming practices of file_that_includes_function.function_name()
   are followed below. In psi4 input files, it is only necessary to call the
   function name alone. That is, use ``energy('scf')``, not ``driver.energy('scf')``.

.. note:: The Python options referred to in the :ref:`sec:psithonFunc` section below
   are placed as arguments to a Python
   function (like ``energy()``), not in ``set`` blocks or commands.
.. comment and indexed in :ref:`apdx:options_py`

.. note:: All |PSIfour| keyword names and values are insensitive to case, both
   those that are placed in ``set`` blocks and as Python function arguments.
   The one exception is documented for the *subset* option in the :py:func:`~wrapper_database.database`
   function, where case structure must match the database file.

.. _`op_py_bool`:

.. _`op_py_boolean`:

.. note:: Boolean options can be specified by ``yes``, ``on``, ``true``, or ``1``
    for affirmative and ``no``, ``off``, ``false``, or ``0`` for negative,
    all insensitive to case.

.. _`op_py_dertype`:

.. note:: The derivative level type for :py:func:`~psi4.optimize` and :py:func:`~psi4.frequency` functions can be specified by ``energy``, ``none``, or ``0`` for 0th derivative, ``gradient``, ``first``, or ``1`` for 1st derivative, and ``hessian``, ``second``, or ``2`` for 2nd derivative. For finite difference, as opposed to analytic, derivatives, the :term:`POINTS <POINTS (FINDIF)>` option can be increased to ``5`` for greater accuracy at increased cost.

.. _`op_py_function`:

.. note:: Function option for the Psithon function called by the current function;
   the default is usually :py:func:`~psi4.energy`. See Sec. :ref:`sec:intercalls`
   for a fuller description. Note that the value of the keyword is a Python object
   and so is not wrapped in quotes like a string.

.. _`op_py_molecule`:

.. note:: The molecule to be acted upon by the current function; the default is the
   "active" molecule, which is the nearest preceeding molecule declared in a
   ``molecule mymol {...}`` block or in an ``activate(mymol)`` statement. Note
   that the value of this keyword (``mymol`` in the example) is a Python object
   and so is not wrapped in quotes like a string. Technically, this is a
   :py:class:`~psi4.core.Molecule` object.

