General Notes About the Code

1. While the primary purpose of the code is transformations between the AO and MO bases, 
   some codes (mrcc, detci, and especially cc) request that tasks involving frozen core
   orbitals be either done by libtrans or converted into effective quantities free of frozen
   core orbitals - that way, they don't need to worry about frozen core orbitals at all, and
   the calculation retains all the simplicity of one where core electrons simply don't exist.
   In particular, libtrans:integraltransform_sort_so_tei.cc has the following responsibilities:
   * Computing all energy contributions involving frozen core orbitals only and putting that
     result into frozen_core_energy_. Used to sanity-check the HF energy.
   * Constructing the "frozen-core operator", which is the core hamiltonian for non-frozen orbitals
     plus the Couloumb and exchange contributions terms arising from the electric field of
     the frozen core orbitals. Think of it as halfway between the core Hamiltonian and the Fock operator.
