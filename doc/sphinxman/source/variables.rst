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

.. _`sec:qcvars_functions`:

QCVariables: setting and accessing
==================================

.. _`table:qcvars_access`:

.. table:: API and conventions to access and program with QCVariables.
   :align: left

   +-------------------------------+-----------------------------------------+---------------------------------------------------+--------------------------------+----------------------------------------------+
   | Operation                     | Py-side Globals [#e1]_ [#e5]_           | Py-side Wavefunction [#e1]_ [#e5]_                | C-side Globals [#e3]_ [#e4]_   | C-side Wavefunction [#e5]_                   |
   +===============================+=========================================+===================================================+================================+==============================================+
   | access var                    | :py:func:`~psi4.core.variable()`        | :py:func:`~psi4.core.Wavefunction.variable()`     | ``val = P::e.globals["KEY"];`` | ``Wavefunction::scalar_variable("kEy")``     |
   |                               |                                         |                                                   | ``val = P::e.arrays["KEY"];``  | ``Wavefunction::array_variable("kEy")``      |
   +-------------------------------+-----------------------------------------+---------------------------------------------------+--------------------------------+----------------------------------------------+
   | access all vars               | :py:func:`~psi4.core.variables()`       | :py:func:`~psi4.core.Wavefunction.variables()`    | ``P::e.globals;``              | ``Wavefunction::scalar_variables()``         |
   |                               |                                         |                                                   | ``P::e.arrays;``               | ``Wavefunction::array_variables()``          |
   +-------------------------------+-----------------------------------------+---------------------------------------------------+--------------------------------+----------------------------------------------+
   | is var set?                   | :py:func:`~psi4.core.has_variable()`    | :py:func:`~psi4.core.Wavefunction.has_variable()` | ``P::e.globals.count("KEY");`` | ``Wavefunction::has_scalar_variable("kEy")`` |
   |                               |                                         |                                                   | ``P::e.arrays.count("KEY");``  | ``Wavefunction::has_array_variable("kEy")``  |
   +-------------------------------+-----------------------------------------+---------------------------------------------------+--------------------------------+----------------------------------------------+
   | set a var value [#e2]_        | :py:func:`~psi4.core.set_variable()`    | :py:func:`~psi4.core.Wavefunction.set_variable()` | ``P::e.globals["KEY"] = val;`` | ``Wavefunction::set_scalar_variable("kEy")`` |
   |                               |                                         |                                                   | ``P::e.arrays["KEY"] = val;``  | ``Wavefunction::set_array_variable("kEy")``  |
   +-------------------------------+-----------------------------------------+---------------------------------------------------+--------------------------------+----------------------------------------------+
   | unset var                     | :py:func:`~psi4.core.del_variable()`    | :py:func:`~psi4.core.Wavefunction.del_variable()` | ``P::e.globals.erase("KEY")``  | ``Wavefunction::del_scalar_variable("kEy")`` |
   |                               |                                         |                                                   | ``P::e.arrays.erase("KEY")``   | ``Wavefunction::del_array_variable("kEy")``  |
   +-------------------------------+-----------------------------------------+---------------------------------------------------+--------------------------------+----------------------------------------------+
   | unset all vars                | :py:func:`~psi4.core.clean_variables()` |                                                   |                                |                                              |
   +-------------------------------+-----------------------------------------+---------------------------------------------------+--------------------------------+----------------------------------------------+
   | register var with docs [#e6]_ | ``core.set_variable("GIBBS", G) # P::e THERMO`` or                                          | ``P::e.globals["CC T1 DIAG"] = val;`` or                                      |
   |                               | ``# wfn.set_variable("ROOT n") # P::e ADC``                                                 | ``// P::e.globals["CCSD ROOT n"];``                                           |
   +-------------------------------+-----------------------------------------+---------------------------------------------------+--------------------------------+----------------------------------------------+

.. [#e1] Works on float or array variables
.. [#e2] Py-side, value can be float or :py:class:`psi4.core.Matrix` or ``np.ndarray``. C-side, value can be float or Matrix, and the appropriate function must be used.
.. [#e3] ``P::e`` stands for ``Process.environment``, which is where all the "globals" live in |PSIfour| .
.. [#e4] ``P::e`` access is direct, not through functions, so all variable names must be ALL CAPS.
.. [#e5] Note that QCVariables are stored all uppercase. Getter/setter functions will do the conversion, so you can use arbitrary case (indicated by ``"kEy"``.
.. [#e6] In-code way works for explicit keys; in-comment way needed with programmatic keys or those with variables (the lowercase "n" here). Py-side needs the module specified explicitly, while c-side gets extracted from the src directory.

