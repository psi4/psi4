
.. warning:: Python naming practices of file_that_includes_function.function_name()
    are followed below. In psi4 input files, it is only necessary to call the
    function name alone. That is, use ``energy('scf')``, not ``driver.energy('scf')``.

.. note:: The options documented below are placed as arguments in the command that
    calls the Python function, not in the ``set globals`` block or with any 
    other ``set`` command.

.. note:: Psithon keyword names and values are insensitive to case. The few
    exceptions are documented for the ``database()`` function, where case
    structure must match the database file.

.. _`bool`:
.. _`boolean`:
.. note:: Boolean arguments can be specified by ``yes``, ``on``, ``true``, or ``1``
    for affirmative and ``no``, ``off``, ``false``, or ``0`` for negative,
    all insensitive to case.

.. _`conv double`:
.. note:: Certain convergence and tolerance keywords, of type *double* (real numbers), may be specified using either a real numberor an integer; an integer *X* is then treated as the number of converged decimal digits required. For example, to request an energy converged to :math:`10^{-6} E_h`, the user may set the ``e_convergence`` keyword to 0.000001, 1.0e-6, or 6.

.. _`dertype string`:
.. note:: The derivative level type for :py:func:`driver.optimize` and :py:func:`driver.frequency`
    functions can be specified by ``energy``, ``none``, or ``0`` for 0th derivative,
    ``gradient``, ``first``, or ``1`` for 1st derivative, and ``hessian``,
    ``second``, or ``2`` for 2nd derivative.

