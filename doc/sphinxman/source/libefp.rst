
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


.. index:: EFP; library fragments
   pair: EFP; adding new

.. _`sec:findingEFPFragments`:

EFP Fragments
~~~~~~~~~~~~~

LIBEFP comes with a couple dozen ready-to-use fragments (water, benzene,
common solvents, etc.) listed :ref:`here <sec:availableEFPFragments>` with source :source:`lib/libfrag`.
Any of these may be used directly in a |PSIfour| input files as described
:ref:`here <sec:usingEFPFragments>`.

.. note:: The built-in fragment library distributed with Q-Chem (as of version 4.0.1) is *not*
   equivalent to that distributed with LIBEFP. Although many of the same
   molecules are present and should perform similarly in computations,
   exact matches of fragment geometries or efp energies should not be
   expected. See files in test case directories :source:`qc-efpefp-sp1
   <tests/libefp/qc-efpefp-sp1>` and :source:`qc-scfefp-sp1
   <tests/libefp/qc-scfefp-sp1>` for equivalent Q-Chem and |PSIfour|
   EFP input files.

Creating new efp fragments requires the `GAMESS
<http://www.msg.ameslab.gov/gamess/>`_ quantum chemistry package.
Instructions on building new fragments are `here
<https://github.com/libefp/libefp#how-to-create-custom-efp-fragment-types>`_.
Once your new fragment is ready, make it assessible to |PSIfour| by
including the directory in which the ``.efp`` file is located in the colon
separated environment variable :envvar:`PSIPATH`. Fragments are searched
for first in the library, next in the paths of :envvar:`PSIPATH`, and
finally in the current directory. If |PSIfour| is unable to find the
fragment, an error will be reported.

.. note:: When constructing new fragment files, the name of the name of the
   file should be lowercase and have extension ``.efp``. The molecule name
   within the file, e.g., ``$NH3`` must correspond to the name of the
   fragment file.

.. index:: molecule; EFP
   single: EFP; molecule specification

.. _`sec:usingEFPFragments`:

Molecule Specification
~~~~~~~~~~~~~~~~~~~~~~

EFP fragment geometries are specified alongside the quantum mechanical
(QM) molecule and make use of the ``--`` fragment separation scheme
described :ref:`here <sec:fragments>`. Each EFP fragment has its own
fragment section that includes the label ``efp``, the name of the file
*fragname* from which EFP parameters are to be read, and the position
specification for the fragment in one of two ways, XYZABC or POINTS. For
XYZABC, the fragment specification is all on one line. ``efp`` and
*fragname* are followed by two sets of three numbers: the coordinates of
the center of mass of the fragment and the three Euler angles that specify
orientation about the center of mass. (Ilya- aren't there different
conventions for Euler angles?) ::

    efp  nh3  0.0 0.0 5.0  5 2 8

For POINTS, the fragment specification is four lines, the first of which
is ``efp`` and *fragname*. The next lines are the coordinates
(without element labels) of the first three atoms in the fragment. Note
that EFP fragment geometries are rigid, so the first atom will be placed
exactly where specified by the first point, the second atom will be placed
along the vector between the first and second points, and the third atom
will be placed in the plane formed by the three points. ::

    efp ch3oh
    1.275    -2.447    -4.673
    0.709    -3.191    -3.592
    2.213    -1.978    -4.343

:ref:`Just as for QM <sec:moleculeKeywords>`, the center of mass
coordinates in the XYZABC format and all coordinates in the POINTS format are
taken to be in Angstroms by default or in Bohr if ``units au`` is present.
Charge and multiplicity specifications are encoded in the fragment file
and so are not read from input.

Any combination of EFP and QM fragments can be placed in a molecule; even
the oddity below is legitimate. Note that reorientation is automatically
turned off when EFP fragments are present (``no_com`` and ``no_reorient``
are implied). ::

    molecule qmefp {
      efp nh3 0.0 0.0 5.0 5 2 8
      --
      C  0.0 0.0 0.0
      O  0.0 1.5 0.0
      O  0.0 -1.5 0.0
      --
      efp h2o 5.0 0.0 0.0 5 2 8
      --
      He  -3.0 4.0 4.0
      He  -4.0 5.0 4.0
      --
      efp ch3oh
      1.275    -2.447    -4.673
      0.709    -3.191    -3.592
      2.213    -1.978    -4.343
    }


Running EFP 
~~~~~~~~~~~~
EFP can be invoked in similar fashion as other theories provided in |PSIfour|.
For example, if you want to obtain the EFP interaction energy for benzene and two waters,
simply provide the following::

   molecule {
     efp c6h6  0.0 0.0 0.0   0.0 0.0 0.0
     --
     efp h2o   4.0 0.0 0.0   0.0 0.0 0.0
     --
     efp h2o  -4.0 0.0 0.0   0.0 0.0 0.0
   }
   
   energy('efp')

This computation involves purely EFP/EFP fragment interactions and is
performed entirely by the LIBEFP library.  |PSIfour| can also handle mixed
systems of quantum mechanical (QM) and EFP fragments through the native
:ref:`SCF <sec:scf>` code augmented by calls to the LIBEFP library. For
example, turning one of the waters in the example above into a QM
fragment, can be done below. ::

   molecule {
     efp c6h6  0.0 0.0 0.0   0.0 0.0 0.0
     --
     O  4.0   0.0   0.0
     H  4.7   0.7   0.0
     H  3.3  -0.7   0.0 
     --
     efp h2o  -4.0 0.0 0.0   0.0 0.0 0.0
   }
   
   energy('scf')

Anytime an EFP fragment is present in the active molecule, the SCF energy
will include EFP contributions.

ACK, I don't think it's even possible to discard the EFP fragments
once molecule's been defined. There's no deactivate. If one wrote a
mol.qmonly or mol.efponly function, I don't see how to recall the EFP
instance back. At present, to run a plain scf, then an efpscf, you'd
have to define a molecule with only qm fragments, run energy('scf'),
then define the molecule again with qm and efp fragments. yikes.

At this time, |PSIfour| is only able to perform pure-efp single-points and
geometry optimizations and mixed qm/efp SCF single-points.

.. _`table:libefpauto`:

    .. _`table:fnocc_methods`:

    +-------------------------+-------------------------------------------------------------+
    | name                    | calls method                                                |
    +=========================+=============================================================+
    | efp                     | EFP-only interaction energy                                 |
    +-------------------------+-------------------------------------------------------------+


.. index:: EFP; library fragments

.. _`sec:availableEFPFragments`:

Fragment Library
~~~~~~~~~~~~~~~~

Below are documented the EFP fragments available from the LIBEFP library.
These systems are accessible in ``molecule {...}`` blocks without
additional configuration.

----

.. comment This toctree directive only here to suppress warning at build time.
   include line below is doing the work.

.. toctree::
   :hidden:

   autodoc_available_fraglib

.. include:: autodoc_available_fraglib.rst


