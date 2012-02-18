
.. note:: Boolean arguments can be specified by ``yes``, ``on``, ``true``, or ``1``
    for affirmative and ``no``, ``off``, ``false``, or ``0`` for negative,
    all insensitive to case.

.. note:: Similarly, the derivative level for ``optimize()`` and ``frequencies()``
    functions can be specified by ``energy``, ``none``, or ``0`` for 0th derivative,
    ``gradient``, ``first``, or ``1`` for 1st derivative, and ``hessian``,
    ``second``, or ``2`` for 2nd derivative.

.. note:: Psithon keyword names and values are insensitive to case. The few
    exceptions are documented for the ``database()`` function, where case
    structure must match the database file.

.. warning:: Python naming practices of file_that_includes_function.function_name()
    are followed below. In psi4 input files, it is only necessary to call the
    function name alone. That is, use ``energy('scf')``, not ``driver.energy('scf')``.


