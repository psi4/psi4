
.. include:: autodoc_abbr_options_c.rst

.. index::
   triple: setting; keywords; optimize()
   pair: geometry optimization; function call
   pair: gradient; function call
   see: opt(); optimize()

.. _`sec:opt()`:

Geometry Optimization, :py:func:`~driver.optimize` and :py:func:`~driver.gradient`
==================================================================================

* :ref:`Psi4 Native Gradient Methods <table:grad_gen>`
* :ref:`Psi4 Native DFT Gradient Methods (excepting double-hybrids) <table:grad_gen>`
* :ref:`CFOUR Interfaced Gradient Methods <table:grad_cfour>`

For further discussion of geometry optimization, see
Sec. :ref:`sec:optking`.

:py:func:`~driver.optimize` is the only command most users will ever 
need to access directly to perform geometry optimizations. Behind
the scenes, :py:func:`~driver.optimize` is a wrapper that repeatedly
calls :py:func:`~driver.gradient` that computes the gradient then adds a
call to the :ref:`geometry projection module <sec:optking>`.

.. autofunction:: driver.optimize(name [, molecule, return_wfn, func, mode, dertype, hessian_with])

.. autofunction:: driver.gradient(name [, molecule, return_wfn, func, dertype])

