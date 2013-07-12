
.. include:: autodoc_abbr_options_c.rst

.. index:: LIBEFP, EFP
.. _`sec:libefp`:

Interface to LIBEFP by I. Kaliman
=================================

.. codeauthor:: A. Eugene DePrince, Andrew C. Simmonett, Rollin A. King, and Lori A. Burns
.. sectionauthor:: Lori A. Burns

*Module:* :ref:`Keywords <apdx:efp>`, :ref:`PSI Variables <apdx:efp_psivar>`, :source:`LIBEFP <src/lib/libefp_solver>`

|PSIfour| contains code to interface to the LIBEFP library developed
in L. Slipchenko's group by I. Kaliman.  LIBEFP at version 0.9.9
is distributed with |Psifour| and requires no additional licence,
downloads, or configuration.  More information about the LIBEFP project
is available at `http://www.libefp.org/ <http://www.libefp.org/>`_
and source is hosted at `https://github.com/libefp/libefp
<https://github.com/libefp/libefp>`_.

Installation
~~~~~~~~~~~~

No installation is required as LIBEFP is packaged with |Psifour|.

EFP Fragments
~~~~~~~~~~~~~

LIBEFP comes with a couple dozen ready-to-use fragments (water, benzene,
common solvents, etc.)  See LINK for available fragments and LINK for
how to use them in a |PSIfour| input file.

.. note:: As of Q-Chem 4.0.1, their built-in fragment library is *not*
   equivalent to that distributed with LIBEFP. Although many of the same
   molecules are present and should perform similarly in computations,
   exact matches of fragment geometries or efp energies should not be
   expected. See files in test case directories :source:`qc-efpefp-sp1
   <tests/libefp/qc-efpefp-sp1>` and :source:`qc-scfefp-sp1
   <tests/libefp/qc-scfefp-sp1>` for equivalent Q-Chem and |PSIfour|
   efp input files.

Creating new efp fragments requires the `GAMESS
<http://www.msg.ameslab.gov/gamess/>`_ quantum chemistry
package. Instructions on building new fragments are `here
<https://github.com/libefp/libefp#how-to-create-custom-efp-fragment-types>`_.

Once your new fragments are ready, make them assessible to |PSIfour|
by ensuring that the directory in which they are located is in the
environment variable :envvar:`PSI_PATH`.  If |PSIfour| is unable to
find the fragment, an error will be reported.

Running EFP 
~~~~~~~~~~~~
EFP can be invoked in similar fashion as other theories provided in |PSIfour|.
For example, if you want to obtain the EFP interaction energy for benzene and two waters,
simply provide the following::

   molecule {
     efp c6h6_l  0.0 0.0 0.0   0.0 0.0 0.0
     --
     efp h2o_l   4.0 0.0 0.0   0.0 0.0 0.0
     --
     efp h2o_l  -4.0 0.0 0.0   0.0 0.0 0.0
   }
   
   energy('efp')

This computation involves purely EFP/EFP fragment interactions and is
performed entirely by the LIBEFP library.  |PSIfour| can also handle
mixed systems of quantum mechanical (QM) and EFP fragments through
the native :ref:`SCF <sec:scf>` code augmented by calls to the LIBEFP
library. For example, turning one of the waters in the example above into a QM fragment, can be done below.

   molecule {
     efp c6h6_l  0.0 0.0 0.0   0.0 0.0 0.0
     --
     O  4.0   0.0   0.0
     H  4.7   0.7   0.0
     H  3.3  -0.7   0.0 
     --
     efp h2o_l  -4.0 0.0 0.0   0.0 0.0 0.0
   }
   
   energy('scf')

Anytime an EFP fragment is present in the active molecule, the SCF energy will include EFP contributions.

ACK, I don't think it's even possible to discard the EFP fragments
once molecule's been defined. There's no deactivate. If one wrote a
mol.qmonly or mol.efponly function, I don't see how to recall the EFP
instance back. At present, to run a plain scf, then an efpscf, you'd
have to define a molecule with only qm fragments, run energy('scf'),
then define the molecule again with qm and efp fragments. yikes.

``'mrccsdt'`` in the call to :py:func:`~driver.energy` instructs |PSIfour| to first
perform an RHF calculation and then call MRCC to compute the CCSDT energy.
For a CCSDT(Q) energy, simply use ``'mrccsdt(q)'`` in the call to
:py:func:`~driver.energy`. MRCC can be used to perform geometry optimization and
frequency calculations for electronic ground states only.

At this time, |PSIfour| is only able to perform pure-efp single-points and geometry optimizations and mixed qm/efp SCF single-points.

.. _`table:libefpauto`:

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


Molecule Specification
~~~~~~~~~~~~~~~~~~~~~~


.. index:: EFP; library fragments

.. _`sec:availableEFPFragments`:

Fragment Library
~~~~~~~~~~~~~~~~

Below are documented the EFP fragments available from the LIBEFP library.
These are the molecules accessible with ``_l`` in ``molecule {...}``
blocks.

----

.. comment This toctree directive only here to suppress warning at build time.
   include line below is doing the work.

.. toctree::
   :hidden:

   autodoc_available_fraglib

.. include:: autodoc_available_fraglib.rst


