
Notes on Options
================

.. note:: The options referred to in the :ref:`sec:methods` section below
   and indexed in :ref:`apdx:options_c_module` are placed in ``set`` blocks as
   described in :ref:`sec:jobControl`, not as arguments to a Python function
   (like ``energy()``).

.. note:: All |PSIfour| keyword names and values are insensitive to case, both
   those that are placed in ``set`` blocks and as Python function arguments.
   The few exceptions are documented for the :py:func:`~driver.database` function,
   where case structure must match the database file.

.. _`op_c_bool`:
.. _`op_c_boolean`:
.. note:: Boolean options can be specified by ``yes``, ``on``, ``true``, or ``1``
    for affirmative and ``no``, ``off``, ``false``, or ``0`` for negative,
    all insensitive to case.

.. _`op_c_conv`:
.. note:: Certain convergence and tolerance keywords, of type *double* (real numbers),
   may be specified using either a real number of an integer; and integer *X* is then
   treated as the number of converged decimal digits required. For example, to request
   as energy converged to :math:`10^{-6} E_h`, the user may set the ``e_convergence``
   keyword to 0.000001, 1.0e-6, or 6.

