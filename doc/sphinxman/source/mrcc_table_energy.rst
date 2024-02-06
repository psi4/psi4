.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2024 The Psi4 Developers.
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

.. _`table:energy_mrcc`:

   +-----------------------+--------------------------------------------------------------------------------------------+
   | :term:`QC_MODULE <QC_MODULE (GLOBALS)>`\ =MRCC                                                                     |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | name                  | calls method in Kallay's MRCC program :ref:`[manual] <sec:mrcc>`                           |
   +=======================+============================================================================================+
   | ccsd                  | CC through doubles :ref:`[details] <dd_ccsd>`                                              |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdt                 | CC through triples                                                                         |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdtq                | CC through quadruples                                                                      |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdtqp               | CC through quintuples                                                                      |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdtqph              | CC through sextuples                                                                       |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsd(t)               | CC through doubles with perturbative triples :ref:`[details] <dd_ccsd_prt_pr>`             |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdt(q)              | CC through triples with perturbative quadruples                                            |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdtq(p)             | CC through quadruples with pertubative quintuples                                          |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdtqp(h)            | CC through quintuples with pertubative sextuples                                           |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsd(t)_l             | CC through doubles with asymmetric perturbative triples :ref:`[details] <dd_accsd_prt_pr>` |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdt(q)_l            | CC through triples with asymmetric perturbative quadruples                                 |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdtq(p)_l           | CC through quadruples with asymmetric perturbative quintuples                              |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdtqp(h)_l          | CC through quintuples with asymmetric perturbative sextuples                               |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdt-1a              | CC through doubles with iterative triples (cheapest terms)                                 |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdtq-1a             | CC through triples with iterative quadruples (cheapest terms)                              |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdtqp-1a            | CC through quadruples with iterative quintuples (cheapest terms)                           |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdtqph-1a           | CC through quintuples with iterative sextuples (cheapest terms)                            |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdt-1b              | CC through doubles with iterative triples (cheaper terms)                                  |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdtq-1b             | CC through triples with iterative quadruples (cheaper terms)                               |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdtqp-1b            | CC through quadruples with iterative quintuples (cheaper terms)                            |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdtqph-1b           | CC through quintuples with iterative sextuples (cheaper terms)                             |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | cc2                   | approximate CC through doubles :ref:`[details] <dd_cc2>`                                   |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | cc3                   | approximate CC through triples :ref:`[details] <dd_cc3>`                                   |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | cc4                   | approximate CC through quadruples                                                          |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | cc5                   | approximate CC through quintuples                                                          |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | cc6                   | approximate CC through sextuples                                                           |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdt-3               | CC through doubles with iterative triples (all but the most expensive terms)               |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdtq-3              | CC through triples with iterative quadruples (all but the most expensive terms)            |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdtqp-3             | CC through quadruples with iterative quintuples (all but the most expensive terms)         |
   +-----------------------+--------------------------------------------------------------------------------------------+
   | ccsdtqph-3            | CC through quintuples with iterative sextuples (all but the most expensive terms)          |
   +-----------------------+--------------------------------------------------------------------------------------------+

