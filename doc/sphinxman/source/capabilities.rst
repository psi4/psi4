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

.. _`sec:mtd_capabilities`:

Capabilities and Alternate Implementations
==========================================

Depending on the reference (e.g., RHF, UHF, ROHF), the integral treatment
(conventional CONV, density-fitted DF, and Cholesky-decomposed CD),
active orbitals, and derivative level, computational methods are
sometimes assembled from implementations by multiple coders. Conversely,
some capabilities may be implemented multiple times. Capabilities,
modules, algorithm types, and defaults are detailed for many ground-state
methods at Table :ref:`Module Capabilities <table:managedmethods>`. Its
analogous summary table (with modules collapsed) is at :ref:`Summary
Capabilities <table:stdsuite>`. |PSIfour| transparently selects the
most efficient implementation, so one generally needn't consult this
table. However, to understand the details of what combinations are
accessible or what alternate implementations are available, read on.

Note that HF, DFT, and MPn (n<=3)
default to density-fitted integrals, while all higher methods default to
conventional integrals. Therefore, for a closed-shell molecule::

    # runs MP2 with default algorithm type ``DF`` with default implementation DFMP2
    energy('mp2')

    # runs MP2 with algorithm type ``CONV`` with default implementation OCC
    set mp2_type conv
    energy('mp2')

    # runs MP2 with default algorithm type ``DF`` with implementation OCC
    set qc_module occ
    energy('mp2')

    # runs MP2 with algorithm type ``CONV`` with implementation FNOCC
    set mp2_type conv
    set qc_module fnocc
    energy('mp2')


.. include:: autodoc_capabilities_details.rst

