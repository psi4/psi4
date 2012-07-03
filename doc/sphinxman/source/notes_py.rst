
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
   The one exception is documented for the *subset* option in the :py:func:`~wrappers.database` 
   function, where case structure must match the database file.

.. _`op_py_bool`:

.. _`op_py_boolean`:

.. note:: Boolean options can be specified by ``yes``, ``on``, ``true``, or ``1``
    for affirmative and ``no``, ``off``, ``false``, or ``0`` for negative,
    all insensitive to case.

.. _`op_py_dertype`:

.. note:: The derivative level type for :py:func:`~driver.optimize` and :py:func:`~driver.frequency` functions can be specified by ``energy``, ``none``, or ``0`` for 0th derivative, ``gradient``, ``first``, or ``1`` for 1st derivative, and ``hessian``, ``second``, or ``2`` for 2nd derivative.

.. _`op_py_function`:

.. note:: Function option for the Psithon function called by the current function;
   the default is usually :py:func:`~driver.energy`. See Sec. :ref:`sec:intercalls`
   for a fuller description. Note that the value of the keyword is a Python object
   and so is not wrapped in quotes like a string.

.. _`op_py_molecule`:

.. note:: The molecule to be acted upon by the current function; the default is the
   nearest preceeding molecule declared in a ``molecule name {...}`` block. Note
   that the value of this keyword (``name`` in the example) is a Python object and
   so is not wrapped in quotes like a string.

