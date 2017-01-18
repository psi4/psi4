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

.. index::
   triple: setting; keywords; cbs()
   see: complete_basis_set(); cbs()
   single: basis set; delta correction

.. _`sec:cbs()`:

Complete Basis Set
==================

.. toctree::
   :hidden:

   cbs_eqn

.. codeauthor:: Lori A. Burns and Daniel G. A. Smith
.. sectionauthor:: Lori A. Burns

The :py:func:`psi4.cbs` function described below is
powerful but complicated, requiring many options. For most common
calculations, a shorthand can be accessed directly though
:py:func:`psi4.energy`, :py:func:`psi4.gradient`, *etc.* For example,
a MP2 single-point DT extrapolation can be accessed through the first item
below more conveniently than the equivalent second item.

* ``energy('mp2/cc-pv[dt]z')``

* ``energy(cbs, corl_wfn='mp2', corl_basis='cc-pv[dt]z')``

A CCSD(T) DT coupled-cluster correction atop a TQ MP2 extrapolation
geometry optimization can be accessed through the first item below more
conveniently than the equivalent second item.

* ``optimize('mp2/cc-pv[tq]z + D:ccsd(t)/cc-pvdz')``

* ``optimize(cbs, corl_wfn='mp2', corl_basis='cc-pv[tq]z', delta_wfn='ccsd(t)', delta_basis='cc-pvdz')``

Many examples can be found at :srcsample:`cbs-xtpl-energy`,
:srcsample:`cbs-xtpl-gradient`, :srcsample:`cbs-xtpl-opt`,
:srcsample:`cbs-xtpl-freq`, :srcsample:`cbs-xtpl-func`,
:srcsample:`cbs-xtpl-wrapper`.


.. autofunction:: psi4.cbs(name [, scf_basis, scf_scheme, corl_wfn, corl_basis, corl_scheme, delta_wfn, delta_wfn_lesser, delta_basis, delta_scheme, delta2_wfn, delta2_wfn_lesser, delta2_basis, delta2_scheme, delta3_wfn, delta3_wfn_lesser, delta3_basis, delta3_scheme, delta4_wfn, delta4_wfn_lesser, delta4_basis, delta4_scheme, delta5_wfn, delta5_wfn_lesser, delta5_basis, delta5_scheme])

.. note:: Presently (May 2016), only two of the five delta possibilities are active. Also, temporarily extrapolations are performed on differences of target and scf total energies, rather than on correlation energies directly. This doesn't affect the extrapolated values of the particular formulas defined here (though it does affect the betas, which are commented out), but it is sloppy and temporary and could affect any user-defined corl extrapolations.

.. index::
   pair: cbs(); output

Output
^^^^^^

At the beginning of a cbs() job is printed a listing of the individual
energy calculations which will be performed. The output snippet below is
from the example job [2] above. It shows first each model chemistry needed
to compute the aggregate model chemistry requested through cbs(). Then,
since, for example, an ``energy('ccsd(t)')`` yields CCSD(T), CCSD, MP2,
and SCF energy values, the wrapper condenses this task list into the second
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

.. _`sec:cbs_xtpl`:

.. index::
   single: cbs(); extrapolation schemes
   single: extrapolation schemes
   single: basis set; extrapolation

Extrapolation Schemes
^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: psi4.driver.driver_cbs.xtpl_highest_1

.. autofunction:: psi4.driver.driver_cbs.scf_xtpl_helgaker_2

.. autofunction:: psi4.driver.driver_cbs.scf_xtpl_helgaker_3

.. autofunction:: psi4.driver.driver_cbs.corl_xtpl_helgaker_2

Aliases
^^^^^^^

When a particular composite method or its functional form is going to be
reused often, it is convenient to define an alias to it. A convenient
place for such Python code to reside is in :source:`psi4/driver/aliases.py`
(source location) or ``psi4/lib/psi4/driver/aliases.py`` (installed
location). No recompilation is necessary after defining an alias. Some
existing examples are below.

.. autofunction:: psi4.driver.aliases.sherrill_gold_standard

.. autofunction:: psi4.driver.aliases.allen_focal_point

