
.. warning:: Python naming practices of file_that_includes_function.function_name()
    are followed below. In psi4 input files, it is only necessary to call the
    function name alone. That is, use ``energy('scf')``, not ``driver.energy('scf')``.

.. note:: The options documented below are placed as arguments in the command that
    calls the Python function, not in the ``set globals`` block or with any 
    other ``set`` command.

.. note:: Psithon keyword names and values are insensitive to case. The few
    exceptions are documented for the ``database()`` function, where case
    structure must match the database file.

.. note:: Boolean arguments can be specified by ``yes``, ``on``, ``true``, or ``1``
    for affirmative and ``no``, ``off``, ``false``, or ``0`` for negative,
    all insensitive to case.


