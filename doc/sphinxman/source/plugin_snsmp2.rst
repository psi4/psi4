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

.. index:: SNS-MP2

Spin-Network-Scaled MP2 (SNS-MP2) by D. E. Shaw
===============================================

.. codeauthor:: D. E. Shaw
.. sectionauthor:: Shannon E. Houck

.. image:: https://img.shields.io/badge/home-sns_mp2-5077AB.svg
   :target: https://github.com/DEShawResearch/sns-mp2

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: https://github.com/DEShawResearch/sns-mp2/blob/master/README.md

This plugin is an implementation of the SNS-MP2 algorithm, originally 
developed by McGibbon et. al. [McGibbon:2017:161725]_, which uses neural networking to 
improve the accuracy of MP2 (ref:`dfmp2`) interaction energies for dimers. 
This |PSIfour| plugin allows the user to compute both energies and 
confidence intervals.

Installation
~~~~~~~~~~~~

**Source**

* Download the plugin from the GitHub repository:

  .. code-block:: bash
    >>> git clone https://github.com/DEShawResearch/sns-mp2

* Once dowloaded, the plugin can be installed:

  .. code-block:: bash
    >>> cd {top-level-sns-mp2-directory}
    >>> PSI4_PYTHON=$(head $(which psi4) -n 1 | sed -r 's/^.{2}//')
    >>> PSI4_PYTHON -m pip install .

Sample Input
~~~~~~~~~~~~

A sample input file, borrowed from D. E. Shaw's documentation, is shown below::

   # Sample SNS-MP2 calculation for two helium atoms

   molecule dimer {
   He 0 0 0
   --
   He 2 0 0
   }

   energy('sns-mp2')
    
Note that these two atoms are separated by double dashes, indicating that
the two should be treated as separate molecules. (See 
:ref:`sec:analysis-of-intermolecular-interactions` for more details on 
setting up dimer molecules.) This input file can be run in the usual fashion:

  .. code-block:: bash
    >>> psi4 input.dat


