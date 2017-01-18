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

.. index:: BSSE

.. index::
   triple: setting; keywords; cp
   triple: setting; keywords; vmfc
   single: counterpoise correction

.. _`sec:cp()`:

Basis Set Superposition Corrections
===================================

.. codeauthor:: Daniel G. A. Smith

.. autofunction:: psi4.driver.driver_nbody.nbody_gufunc(func, method_string [, molecule, bsse_type, max_nbody, ptype, return_total_data])


The nbody function computes counterpoise-corrected (CP), non-CP (noCP), and Valiron-Mayer Function Counterpoise (VMFC) interaction energies for complexes composed of arbitrary numbers of monomers.


**Examples :** ::

    # Counterpoise corrected CCSD(T) energy for the Helium dimer
    molecule mol {
      He
      -- 
      He 1 3
    }

    energy('CCSD(T)', bsse_type='cp')   

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

