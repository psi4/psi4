
.. autofunction:: wrappers.complete_basis_set(name [, scf_basis, scf_scheme, corl_wfn, corl_basis, corl_scheme, delta_wfn, delta_wfn_lesser, delta_basis, delta_scheme, delta2_wfn, delta2_wfn_lesser, delta2_basis, delta2_scheme])

Output
======

At the beginning of a cbs() job is printed a listing of the individual
energy calculations which will be performed. The output snippet below is
from the example job [2] above. It shows first each model chemistry needed
to compute the aggregate model chemistry requested through cbs(). Then,
since, for example, an ``energy('ccsd(t)')`` yields CCSD(T), CCSD, MP2,
and SCF energy values, the wrapper condenses this task into the second
list of minimum number of calculations which will actually be run. ::

    Naive listing of computations required.
            scf / aug-cc-pvqz              for  SCF TOTAL ENERGY
            mp2 / aug-cc-pvtz              for  MP2 CORRELATION ENERGY
            mp2 / aug-cc-pvqz              for  MP2 CORRELATION ENERGY
        ccsd(t) / aug-cc-pvdz              for  CCSD(T) CORRELATION ENERGY
        ccsd(t) / aug-cc-pvtz              for  CCSD(T) CORRELATION ENERGY
            mp2 / aug-cc-pvdz              for  MP2 CORRELATION ENERGY
            mp2 / aug-cc-pvtz              for  MP2 CORRELATION ENERGY

    Enlightened listing of computations required.
            mp2 / aug-cc-pvqz              for  MP2 CORRELATION ENERGY
        ccsd(t) / aug-cc-pvdz              for  CCSD(T) CORRELATION ENERGY
        ccsd(t) / aug-cc-pvtz              for  CCSD(T) CORRELATION ENERGY

At the end of a cbs() job is printed a summary section like the one below. First,
in the components section, are listed the results for each model chemistry available, whether
required for the cbs job (*) or not. Next, in the stages section, are listed the results for
each extrapolation. The energies of this section must be dotted with the weightings in column Wt
to get the total cbs energy. Finally, in the CBS section, are listed the results for each stage
of the cbs procedure. The stage energies of this section sum outright to the total cbs energy. ::

    ==> Components <==
    
    ----------------------------------------------------------------------------------
                   Method / Basis            Rqd   Energy [H]   Variable
    ----------------------------------------------------------------------------------
                      scf / aug-cc-pvqz        *  -1.11916375   SCF TOTAL ENERGY
                      mp2 / aug-cc-pvqz        *  -0.03407997   MP2 CORRELATION ENERGY
                      scf / aug-cc-pvdz           -1.11662884   SCF TOTAL ENERGY
                      mp2 / aug-cc-pvdz        *  -0.02881480   MP2 CORRELATION ENERGY
                  ccsd(t) / aug-cc-pvdz        *  -0.03893812   CCSD(T) CORRELATION ENERGY
                     ccsd / aug-cc-pvdz           -0.03893812   CCSD CORRELATION ENERGY
                      scf / aug-cc-pvtz           -1.11881134   SCF TOTAL ENERGY
                      mp2 / aug-cc-pvtz        *  -0.03288936   MP2 CORRELATION ENERGY
                  ccsd(t) / aug-cc-pvtz        *  -0.04201004   CCSD(T) CORRELATION ENERGY
                     ccsd / aug-cc-pvtz           -0.04201004   CCSD CORRELATION ENERGY
    ----------------------------------------------------------------------------------
    
    ==> Stages <==
    
    ----------------------------------------------------------------------------------
     Stage         Method / Basis             Wt   Energy [H]   Scheme
    ----------------------------------------------------------------------------------
       scf            scf / aug-cc-pvqz        1  -1.11916375   highest_1
      corl            mp2 / aug-cc-pv[tq]z     1  -0.03494879   corl_xtpl_helgaker_2
     delta        ccsd(t) / aug-cc-pv[dt]z     1  -0.04330347   corl_xtpl_helgaker_2
     delta            mp2 / aug-cc-pv[dt]z    -1  -0.03460497   corl_xtpl_helgaker_2
    ----------------------------------------------------------------------------------
    
    ==> CBS <==
    
    ----------------------------------------------------------------------------------
     Stage         Method / Basis                  Energy [H]   Scheme
    ----------------------------------------------------------------------------------
       scf            scf / aug-cc-pvqz           -1.11916375   highest_1
      corl            mp2 / aug-cc-pv[tq]z        -0.03494879   corl_xtpl_helgaker_2
     delta  ccsd(t) - mp2 / aug-cc-pv[dt]z        -0.00869851   corl_xtpl_helgaker_2
     total            CBS                         -1.16281105
    ----------------------------------------------------------------------------------


.. _sec_cbs_xtpl:

Extrapolation Schemes
=====================

.. autofunction:: wrappers.highest_1

.. autofunction:: wrappers.scf_xtpl_helgaker_2

.. autofunction:: wrappers.scf_xtpl_helgaker_3

.. autofunction:: wrappers.corl_xtpl_helgaker_2

