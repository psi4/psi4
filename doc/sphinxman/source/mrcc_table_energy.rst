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

.. _`table:energy_mrcc`:

    +-------------------------+---------------------------------------------------------------------------------------+
    | name                    | calls method in Kallay's MRCC program :ref:`[manual] <sec:mrcc>`                      |
    +=========================+=======================================================================================+
    | mrccsd                  | CC through doubles                                                                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdt                 | CC through triples                                                                    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtq                | CC through quadruples                                                                 |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqp               | CC through quintuples                                                                 |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqph              | CC through sextuples                                                                  |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsd(t)               | CC through doubles with perturbative triples                                          |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdt(q)              | CC through triples with perturbative quadruples                                       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtq(p)             | CC through quadruples with pertubative quintuples                                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqp(h)            | CC through quintuples with pertubative sextuples                                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsd(t)_l             |                                                                                       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdt(q)_l            |                                                                                       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtq(p)_l           |                                                                                       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqp(h)_l          |                                                                                       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdt-1a              | CC through doubles with iterative triples (cheapest terms)                            |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtq-1a             | CC through triples with iterative quadruples (cheapest terms)                         |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqp-1a            | CC through quadruples with iterative quintuples (cheapest terms)                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqph-1a           | CC through quintuples with iterative sextuples (cheapest terms)                       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdt-1b              | CC through doubles with iterative triples (cheaper terms)                             |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtq-1b             | CC through triples with iterative quadruples (cheaper terms)                          |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqp-1b            | CC through quadruples with iterative quintuples (cheaper terms)                       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqph-1b           | CC through quintuples with iterative sextuples (cheaper terms)                        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrcc2                   | approximate CC through doubles                                                        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrcc3                   | approximate CC through triples                                                        |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrcc4                   | approximate CC through quadruples                                                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrcc5                   | approximate CC through quintuples                                                     |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrcc6                   | approximate CC through sextuples                                                      |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdt-3               | CC through doubles with iterative triples (all but the most expensive terms)          |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtq-3              | CC through triples with iterative quadruples (all but the most expensive terms)       |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqp-3             | CC through quadruples with iterative quintuples (all but the most expensive terms)    |
    +-------------------------+---------------------------------------------------------------------------------------+
    | mrccsdtqph-3            | CC through quintuples with iterative sextuples (all but the most expensive terms)     |
    +-------------------------+---------------------------------------------------------------------------------------+

