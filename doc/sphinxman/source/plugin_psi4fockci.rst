.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2022 The Psi4 Developers.
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

.. index:: PSI4FOCKCI

.. _`sec:fockci`:

Psi4FockCI: A General Fock-Space CI For Spin-Flip And IP/EA
===========================================================

.. codeauthor:: Shannon E. Houck
.. sectionauthor:: Shannon E. Houck

.. image:: https://img.shields.io/badge/home-psi4fockci-5077AB.svg
   :target: https://github.com/shannonhouck/psi4fockci

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: https://shannonhouck.github.io/psi4fockci/build/index.html

This plugin is an implementation of the RAS-nSF-IP/EA approach detailed 
in the paper by Houck et. al. [Houck:2019:2278]_. 
This approach handles systems with both spin and spatial degeneracies 
by combining the spin-flip (SF) [Krylov:2001:522]_
and ionization potential/electron affinity (IP/EA) [Nooijen:1995:3629]_ 
approaches. 

The Psi4FockCI plugin allows one to perform spin-flip (SF), ionization 
potential (IP), and electron affinity (EA) calculations, as well as 
combined RAS-SF-IP/EA calculations, through the DETCI (:ref:`sec:ci`) module.

Installation
~~~~~~~~~~~~

* Download the plugin from the GitHub repository:

  .. code-block:: bash

     >>> git clone https://github.com/shannonhouck/psi4fockci.git

* Once downloaded, the plugin can be installed as follows:

  .. code-block:: bash

     >>> cd {top-level-psi4fockci-directory}
     >>> python -m pip install .

Sample Input
~~~~~~~~~~~~

To run a CAS-nSF-IP/EA calculation, start with a molecule with the 
correct charge and multiplicity for the reference state (i.e. some 
state well-represented by a single reference). Then, run an energy 
calculation, passing in the charge and multiplicity of the 
desired state; the number of spin-flips and IP/EA will be automatically 
determined based on this input. If additional excitations outside of the 
RAS II space are desired, one can set the level of external 
excitations using the ``conf_space`` keyword.

A sample input file for a RAS(h)-2SF-IP calculation is shown below:

.. code-block:: python

    molecule {
    0 7
    N 0.0 0.0 0.0
    N 0.0 0.0 1.3
    symmetry c1
    }

    set {
      basis cc-pVDZ
    }

    energy('psi4fockci', new_charge=1, new_multiplicity=1, conf_space="h")
    
Note that for calculations involving IP/EA, inclusion of hole (for IP) and 
particle (for EA) excitations is strongly recommended. Additional keywords 
can be found in the documentation.

This input file can be run with Psi4:

.. code-block:: bash

   >>> psi4 input.dat


