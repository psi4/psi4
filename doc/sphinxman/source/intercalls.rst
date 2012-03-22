
.. _`sec:intercalls`:

Function Intercalls
===================

For many of the PSI4 Python functions described above, it makes scientific
sense that they could be called in combination. For instance, one could
optimize all the reagents in a database or compute a
counterpoise-corrected interaction energy with an extrapolated method. The
table below outlines permitted intercalls between functions, showing that
db(opt(cbs(energy()))) is allowed, while db(cp(energy())) is not. This
table is not yet validated for calls with cp().

.. _`table:intercalls`:

.. table:: Permitted nesting of Psithon functions

    +-------------------------+-----+-----+-----+-----+--------+
    | Caller                  | Callee                         |
    +-------------------------+-----+-----+-----+-----+--------+
    |                         | cp  | db  | opt | cbs | energy |
    +=========================+=====+=====+=====+=====+========+
    | :ref:`sec:cp()`         |     | --- |  Y  |  Y  |   Y    |
    +-------------------------+-----+-----+-----+-----+--------+
    | :ref:`sec:db()`         | --- |     |  Y  |  Y  |   Y    |
    +-------------------------+-----+-----+-----+-----+--------+
    | :ref:`sec:opt()`        | --- | --- |     |  Y  |   Y    |
    +-------------------------+-----+-----+-----+-----+--------+
    | :ref:`sec:cbs()`        | --- | --- | --- |     |   Y    |
    +-------------------------+-----+-----+-----+-----+--------+
    | :ref:`sec:energy()`     | --- | --- | --- | --- |        |
    +-------------------------+-----+-----+-----+-----+--------+

- The command db(opt(cbs(energy()))) is actually expressed as ``db(...,
  db_func=opt, opt_func=cbs)``. The perhaps expected final argument of
  ``cbs_func=energy`` is not necessary since energy() is always the function
  called by default. Also, the outermost internal function call (``db_func``
  above can be called as just ``func``. Several examples of intercalls
  between Python functions can be found in sample input :srcsample:`pywrap_all`.

- All keyword arguments are passed along to each function traversed in the
  Python driver, so there should be no concern for separating them, grouping
  them, or designating them for a particular function when undertaking a
  nested calculation. Where the same keyword is used by multiple functions,
  prefixes are added, e.g., **db_mode** and **opt_mode**.

- Function intercalls should not be used in sow/reap mode.

