.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2022 The Psi4 Developers.
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
   triple: setting; keywords; frequency()
   pair: vibrational analysis; function call
   pair: hessian; function call
   see: freq(); frequency();
   see: frequencies(); frequency();

.. _`sec:freq()`:

Harmonic Vibrational Analysis and Visualization of Normal Modes |w---w| :py:func:`~psi4.driver.frequency` and :py:func:`~psi4.driver.hessian`
=============================================================================================================================================

* :ref:`Psi4 Native Hessian Methods <table:freq_gen>`

For further discussion of vibrational and thermochemical analysis,
see Sec. :ref:`sec:thermo`.

:py:func:`~psi4.driver.frequency` is the only command most users will ever
need to access directly to perform frequency calculations. Behind
the scenes, :py:func:`~psi4.driver.frequency` is a light wrapper over
:py:func:`~psi4.driver.hessian` that computes the Hessian then adds a
thermochemical analysis.

.. autofunction:: psi4.frequency(name [, molecule, return_wfn, func, mode, dertype, irrep])
   :noindex:

.. autofunction:: psi4.hessian(name [, molecule, return_wfn, func, dertype, irrep])
   :noindex:

It's handy to collect the wavefunction after a frequency
calculation through ``e, wfn = psi4.frequency(...,
return_wfn=True)`` as the frequencies can be accessed through
:py:func:`psi4.core.Wavefunction.frequencies()`, the Hessian through
:py:func:`psi4.core.Wavefunction.hessian()`, and much other computation
info through ``psi4.core.Wavefunction.frequency_analysis``
(note no parentheses). Examples of using this data
structure can be found :srcsample:`fd-freq-gradient` and
:source:`tests/pytests/test_vibanalysis.py`. Formatted printing of vibrational
results is available through :py:func:`psi4.driver.qcdb.vib.print_vibs`.

.. _`table:frequency_analysis`:

.. table:: Results accessible through ``psi4.core.Wavefunction.frequency_analysis``

    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | key           | description (lbl & comment)                | units     | data (real/imaginary modes)                          |
    +===============+============================================+===========+======================================================+
    | omega         | frequency                                  | cm^-1     | ndarray(ndof) complex (real/imag)                    |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | q             | normal mode, normalized mass-weighted      | a0 u^1/2  | ndarray(ndof, ndof) float                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | w             | normal mode, un-mass-weighted              | a0        | ndarray(ndof, ndof) float                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | x             | normal mode, normalized un-mass-weighted   | a0        | ndarray(ndof, ndof) float                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | degeneracy    | degree of degeneracy                       |           | ndarray(ndof) int                                    |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | TRV           | translation/rotation/vibration             |           | ndarray(ndof) str 'TR' or 'V' or '-' for partial     |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | gamma         | irreducible representation                 |           | ndarray(ndof) str irrep or None if unclassifiable    |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | mu            | reduced mass                               | u         | ndarray(ndof) float (+/+)                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | k             | force constant                             | mDyne/A   | ndarray(ndof) float (+/-)                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | DQ0           | RMS deviation v=0                          | a0 u^1/2  | ndarray(ndof) float (+/0)                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | Qtp0          | Turning point v=0                          | a0 u^1/2  | ndarray(ndof) float (+/0)                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | Xtp0          | Turning point v=0                          | a0        | ndarray(ndof) float (+/0)                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+
    | theta_vib     | char temp                                  | K         | ndarray(ndof) float (+/0)                            |
    +---------------+--------------------------------------------+-----------+------------------------------------------------------+


Visualization of Normal Modes
-----------------------------

|PSIfour| has the ability to export a Molden file that stores information about
the harmonic frequencies and normal modes computed via :py:func:`~psi4.driver.frequency`.
This feature can be enabled by setting the option |findif__normal_modes_write| to true.
The filename of the Molden file ends in ``.molden_normal_modes``, and the prefix is
determined by |globals__writer_file_label| (if set), or else by the name of the
output file plus the name of the current molecule.
The normal coordinates saved in the Molden file are normalized and are not
mass weighted.

Molden Interface Keywords
~~~~~~~~~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/findif__normal_modes_write.rst

.. include:: autodir_options_c/globals__writer_file_label.rst


.. automodapi:: psi4.driver.qcdb.vib
   :headings: -~


API
---

.. autopydantic_model:: psi4.driver.driver_findif.FiniteDifferenceComputer

