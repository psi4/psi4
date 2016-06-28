
.. include:: autodoc_abbr_options_c.rst

.. index::
   triple: setting; keywords; cp
   see: counterpoise_correct

.. _`sec:cp()`:

Counterpoise Correct
====================

.. codeauthor:: Daniel G. A. Smith

.. autofunction:: driver_nbody._nbody_gufunc(func, method_string [, molecule, bsse_type, max_nbody, ptype, return_total_data])


The nbody function computes counterpoise-corrected (CP), non-CP (noCP), and Valiron-Mayer Function Couterpoise (VMFC) interaction energies for complexes composed of arbitrary numbers of monomers.


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

