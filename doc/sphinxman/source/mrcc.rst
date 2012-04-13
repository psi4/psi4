
.. include:: autodoc_abbr_options_c.rst

.. index:: MRCC
.. _`sec:mrcc`:

Interface to MRCC by M. K\ |a_acute|\ llay
==========================================

.. codeauthor:: Justin M. Turney, Andrew C. Simmonett
.. sectionauthor:: Justin M. Turney

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

``'mrccsdt'`` in the call to :py:func:`~driver.energy` instructs |PSIfour| to first
perform an RHF calculation and then call MRCC to compute the CCSDT energy.
For a CCSDT(Q) energy, simply use ``'mrccsdt(q)'`` in the call to
:py:func:`~driver.energy`. MRCC can be used to perform geometry optimization and
frequency calculations for electronic ground states only.

At this time, |PSIfour| is only able to automatically generate the proper
input file for MRCC for the methods listed in table below.
To utilize any method described in the table, you must prefix
the method name with ``MR``. For other methods, you will be required to
use the MRCC keywords described in Appendix :ref:`apdx:mrcc`.

.. _`table:mrccauto`:

.. table:: Methods available in automatic interface with MRCC

   +------------+-------------------+--------------------+
   |CCSD        | CCSD(T) [1]_      | CCSD(T)\_L [1]_    |
   +------------+-------------------+--------------------+
   |CCSDT       | CCSDT(Q) [1]_     | CCSDT(Q)\_L [1]_   |
   +------------+-------------------+--------------------+
   |CCSDTQ      | CCSDTQ(P) [1]_    | CCSDTQ(P)\_L [1]_  |
   +------------+-------------------+--------------------+
   |CCSDTQP     | CCSDTQP(H) [1]_   | CCSDTQP(H)\_L [1]_ |
   +------------+-------------------+--------------------+
   |CCSDTQPH    |                   |                    |
   +------------+-------------------+--------------------+
   |CCSDT-1a    | CCSDT-1b          | CCSDT-3            |
   +------------+-------------------+--------------------+
   |CCSDTQ-1a   | CCSDTQ-1b         | CCSDTQ-3           |
   +------------+-------------------+--------------------+
   |CCSDTQP-1a  | CCSDTQP-1b        | CCSDTQP-3          |
   +------------+-------------------+--------------------+
   |CCSDTQPH-1a | CCSDTQPH-1b       | CCSDTQPH-3         |
   +------------+-------------------+--------------------+
   |CC2         |                   |                    |
   +------------+-------------------+--------------------+
   |CC3         |                   |                    |
   +------------+-------------------+--------------------+
   |CC4         |                   |                    |
   +------------+-------------------+--------------------+
   |CC5         |                   |                    |
   +------------+-------------------+--------------------+
   |CC6         |                   |                    |
   +------------+-------------------+--------------------+

.. [1] Pertubative methods not available with ROHF reference.

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


