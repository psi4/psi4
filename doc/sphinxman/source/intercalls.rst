.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2017 The Psi4 Developers.
.. #
.. # The copyrights for code used from other parties are included in
.. # the corresponding files.
.. #
.. # This program is free software; you can redistribute it and/or modify
.. # it under the terms of the GNU General Public License as published by
.. # the Free Software Foundation; either version 2 of the License, or
.. # (at your option) any later version.
.. #
.. # This program is distributed in the hope that it will be useful,
.. # but WITHOUT ANY WARRANTY; without even the implied warranty of
.. # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. # GNU General Public License for more details.
.. #
.. # You should have received a copy of the GNU General Public License along
.. # with this program; if not, write to the Free Software Foundation, Inc.,
.. # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
.. #
.. # @END LICENSE
.. #

.. include:: autodoc_abbr_options_c.rst

.. _`sec:intercalls`:

Function Intercalls
===================

This topic is in transition. As of 1.0, the functionality provided by
``cp()`` and ``cbs()`` should not be accessed directly. Instead, use the
``cp`` kwarg and/or the extended method syntax like
``'mp3/aug-cc-pv[dt]z'`` to ``energy()``, ``opt()``, *etc.*, respectively.

For many of the |PSIfour| Python functions described above, it makes scientific
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
  between Python functions can be found in sample input :srcsample:`pywrap-all`.

- All keyword arguments are passed along to each function traversed in the
  Python driver, so there should be no concern for separating them, grouping
  them, or designating them for a particular function when undertaking a
  nested calculation. Where the same keyword is used by multiple functions,
  prefixes are added, *e.g.*, **db_mode** and **opt_mode**.

- Function intercalls should not be used in sow/reap mode.

