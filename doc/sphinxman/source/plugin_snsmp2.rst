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

.. include:: autodoc_abbr_options_c.rst

.. index:: SNS-MP2

.. _`sec:snsmp2`:

Spin-Network-Scaled MP2 (SNS-MP2) by D. E. Shaw
===============================================

.. codeauthor:: D. E. Shaw Research
.. sectionauthor:: Shannon E. Houck

.. image:: https://img.shields.io/badge/home-sns--mp2-5077AB.svg
   :target: https://github.com/DEShawResearch/sns-mp2

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: https://github.com/DEShawResearch/sns-mp2/blob/master/README.md

This plugin is an implementation of the SNS-MP2 algorithm developed by McGibbon 
et. al. [McGibbon:2017:161725]_. The SNS-MP2 method uses neural networking to 
improve the accuracy of MP2 (:ref:`sec:dfmp2`) interaction energies for dimer molecules.
The plugin is distributed under the 2-clause BSD license.

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/snsmp2/badges/version.svg
     :target: https://anaconda.org/psi4/snsmp2

* snsmp2 is available as a conda package for Linux and macOS (and Windows, through the Ubuntu shell).

* If using the Psi4conda installer, snsmp2 has already been installed alongside.

* If using the |PSIfour| conda package, the snsmp2 conda package can
  be obtained through ``conda install snsmp2 -c psi4`` or ``conda install
  psi4-rt -c psi4``.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  snsmp2 can be obtained through ``conda install snsmp2 -c psi4``.
  Then, hint its location with :envvar:`PYTHONPATH`.

* To remove a conda installation, ``conda remove snsmp2``.

**Source**

* .. image:: https://img.shields.io/github/tag/DEShawResearch/sns-mp2.svg?maxAge=2592000
   :target: https://github.com/DEShawResearch/sns-mp2

* Download the plugin from the GitHub repository:

  .. code-block:: bash

     >>> git clone https://github.com/DEShawResearch/sns-mp2

* Once dowloaded, the plugin can be installed as outlined in the documentation:

  .. code-block:: bash

     >>> cd {top-level-sns-mp2-directory}
     >>> PSI4_PYTHON=$(head $(which psi4) -n 1 | sed -r 's/^.{2}//')
     >>> $PSI4_PYTHON -m pip install .

Sample Input
~~~~~~~~~~~~

A sample input file, adapted from the documentation, is shown below::

   # Sample SNS-MP2 calculation for two helium atoms

   molecule dimer {
   He 0 0 0
   --
   He 2 0 0
   }

   energy('sns-mp2')
    
Note that the two monomers are separated by double dashes, indicating that
they should be treated as separate molecules. (See 
:ref:`sec:tutorial_tu5` for more details on
setting up dimer molecules.) This input file can be run in the usual fashion:

.. code-block:: bash

   >>> psi4 input.dat


