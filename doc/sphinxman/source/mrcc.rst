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

.. index:: MRCC
.. _`sec:mrcc`:

Interface to MRCC by M. K\ |a_acute|\ llay
==========================================

.. codeauthor:: Justin M. Turney and Andrew C. Simmonett
.. sectionauthor:: Justin M. Turney


*Module:* :ref:`Keywords <apdx:mrcc>`, :ref:`PSI Variables <apdx:mrcc_psivar>`, :source:`MRCC <src/bin/mrcc>`, :ref:`Samples <apdx:testSuitedftd3>`

|PSIfour| contains code to interface to the MRCC program of M. K\ |a_acute|\ llay
and J. Gauss.  The license and source code of the MRCC program must be
obtained from Mih\ |a_acute|\ ly K\ |a_acute|\ llay (`http://www.mrcc.hu/ <http://www.mrcc.hu/>`_).

Installation
~~~~~~~~~~~~

Follow the instructions provided with the source to build the MRCC programs.
To be used by |PSIfour|, ensure that the program binary (``dmrcc``) can be
found in your :envvar:`PATH`. If |PSIfour| is unable to execute the binary, an
error will be reported.

Running MRCC
~~~~~~~~~~~~
MRCC can be invoked in similar fashion as other theories provided in |PSIfour|.
For example, if you want to obtain the CCSDT energy for water with cc-pVDZ using
MRCC simply provide the following::

   molecule h2o {
        O
        H 1 1.0
        H 1 1.0 2 104.5
   }
   set {
        basis cc-pVDZ
   }
   energy('mrccsdt')

``'mrccsdt'`` in the call to :py:func:`~psi4.energy` instructs |PSIfour| to first
perform an RHF calculation and then call MRCC to compute the CCSDT energy.
For a CCSDT(Q) energy, simply use ``'mrccsdt(q)'`` in the call to
:py:func:`~psi4.energy`. MRCC can be used to perform geometry optimization and
frequency calculations for electronic ground states only.

At this time, |PSIfour| is only able to automatically generate the proper
input file for MRCC for the methods listed in table below.
To utilize any method described in the table, you must prefix
the method name with ``MR``. For other methods, you will be required to
use the MRCC keywords described in Appendix :ref:`apdx:mrcc`.
Perturbative methods (``ccsd(t)``, ``ccsdtqp(h)_l``, etc.)
are available with |scf__reference| ROHF in versions of MRCC published
at least after July 1, 2014.

When using ROHF-CCSDT(Q), MRCC will compute and report two variants:
CCSDT(Q)/A and CCSDT(Q)/B.[Kallay:2008:144101]_ Psi4 will save both energies but will use
the CCSDT(Q)/B as the CCSDT(Q) energy. CCSDT(Q)/B has been found to be
more robust by Martin.[Martin:2014:785]_

.. include:: mrcc_table_energy.rst

Frozen-core approximation is also supported in the MRCC interface.
To optimize CH\ :sub:`4` with CCSDT freezing the 1\ *s* on carbon, run::

   molecule H2O {
       O
       H 1 r
       H 1 r 2 104.5
   
       r = 1.0
   }
   
   set {
       basis cc-pVDZ
       freeze_core true
   }
   
   optimize('mrccsdt')

Interface Details
~~~~~~~~~~~~~~~~~

.. _`table:mrcc__mrcc_method`:

.. table:: MRCC methods 

    +---------------------+--------------+-------------------------------------------------------------+
    | |mrcc__mrcc_method| | Method       | Description                                                 | 
    +=====================+==============+=============================================================+ 
    | 1                   | CC           |                                                             |
    +---------------------+--------------+-------------------------------------------------------------+
    | 2                   | CC(n-1)[n]   |                                                             |
    +---------------------+--------------+-------------------------------------------------------------+
    | 3                   | CC(n-1)(n)   | (CC(n-1)[n] energy is also calculated)                      | 
    +---------------------+--------------+-------------------------------------------------------------+
    | 4                   | CC(n-1)(n)_L | (CC(n-1)[n] and CC(n-1)(n) energies are also calculated)    | 
    +---------------------+--------------+-------------------------------------------------------------+
    | 5                   | CC(n)-1a     |                                                             |
    +---------------------+--------------+-------------------------------------------------------------+
    | 6                   | CC(n)-1b     |                                                             |
    +---------------------+--------------+-------------------------------------------------------------+
    | 7                   | CCn          |                                                             |
    +---------------------+--------------+-------------------------------------------------------------+
    | 8                   | CC(n)-3      |                                                             |
    +---------------------+--------------+-------------------------------------------------------------+

