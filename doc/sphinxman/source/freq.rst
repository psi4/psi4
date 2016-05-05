
.. include:: autodoc_abbr_options_c.rst

.. index::
   triple: setting; keywords; frequency()
   pair: vibrational analysis; function call
   pair: hessian; function call
   see: freq(); frequency();
   see: frequencies(); frequency();

.. _`sec:freq()`:

Harmonic Vibrational Analysis, :py:func:`~driver.frequency` and :py:func:`~driver.hessian`
==========================================================================================

* :ref:`Psi4 Native Hessian Methods <table:freq_gen>`

For further discussion of vibrational and thermochemical analysis,
see Sec. :ref:`sec:thermo`.

:py:func:`~driver.frequency` is the only command most users will ever
need to access directly to perform frequency calculations. Behind
the scenes, :py:func:`~driver.frequency` is a light wrapper over
:py:func:`~driver.hessian` that computes the Hessian then adds a
thermochemical analysis.

.. autofunction:: driver.frequency(name [, molecule, return_wfn, func, mode, dertype, irrep])

.. autofunction:: driver.hessian(name [, molecule, return_wfn, func, dertype, irrep])

