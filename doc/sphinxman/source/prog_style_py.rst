.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2025 The Psi4 Developers.
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

.. _`sec:style_python`:

Python Style
============




.. _`faq:ignoringadvice`:

How to Ignore the Bots
----------------------

Formatting and analysis bots are great because it takes more effort
to defy them than to accept their criticism. Nevertheless, for code
clarity, they can be honestly wrong, so we need a way to specifically
clear their findings.

* Py Formatting (yapf) ``# yapf: disable`` (single line or block)  ``# yapf: enable`` (resume)

* C++ Formatting (clang-format)  ``// clang-format off`` (single line or block) ``// clang-format on`` (resume)

* Py Dynamic Analysis (coverage.py)  ``# pragma: no cover``

* C++ Dynamic Analysis (gcov) https://stackoverflow.com/a/30078276 untested

* Py Static Analysis (lgtm) ``# lgtm[py/not-named-self]`` (click on the "?" to get the "Query ID")

* C++ Static Analysis (lgtm) ``// lgtm[cpp/wrong-type-format-argument]`` (click on the "?" to get the "Query ID")


