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

.. index:: BSSE

.. index::
   triple: setting; keywords; cp
   triple: setting; keywords; vmfc
   single: counterpoise correction
   triple: setting; keywords; mbe
   triple: setting; keywords; ssfc
   single: QCManyBody

.. _`sec:cp()`:

Basis Set Superposition Corrections
===================================

.. codeauthor:: Daniel G. A. Smith

.. autofunction:: psi4.driver.driver_nbody.nbody(func, method_string [, molecule, bsse_type, max_nbody, ptype, return_total_data, supersystem_ie_only])


The nbody function computes counterpoise-corrected (CP), non-CP (noCP), and Valiron-Mayer Function Counterpoise (VMFC) interaction energies for complexes composed of arbitrary numbers of monomers.


**Examples :** ::

    set {
      basis def2-svp
    }

    # Counterpoise corrected CCSD(T) energies for the Helium dimer
    molecule mol {
      He
      --
      He 1 3
    }
    # Calculate interaction energies only (skips monomers in monomer basis):
    energy('CCSD(T)', bsse_type='cp')
    # Calculate interaction and total energies, return interaction energies:
    energy('CCSD(T)', bsse_type=['cp','nocp'])
    # Calculate and return counterpoise-corrected gradient
    # Useful for e.g. CP-corrected geometry optimization
    gradient('CCSD(T)', bsse_type='cp', return_total_data=True)


    # noCP, VMFC, and CP energy for a helium cluster, limited at 3 bodies
    molecule mol {
      He 0 0 0
      --
      He 0 0 4
      --
      He 0 4 0
      --
      He 4 0 0
    }

    # Returns the nocp energy as its first in the list
    energy('CCSD(T)', bsse_type=['nocp', 'cp', 'vmfc'], max_nbody=3)
    # Calculate cp geometry optimization skipping the MBE intermediate levels
    optimize("ccsd(t)/cc-pv[dt]z", bsse_type="cp", supersystem_ie_only=True)

API
---

.. autoclass:: psi4.driver.driver_nbody.BsseEnum
   :members:
   :undoc-members:

.. autopydantic_model:: psi4.driver.driver_nbody.ManyBodyComputer
   :members:
   :undoc-members:
   :inherited-members: BaseModel, ProtoModel

.. autopydantic_model:: qcmanybody.ManyBodyComputer
   :members:
   :undoc-members:

