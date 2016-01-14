.. include:: autodoc_abbr_options_c.rst

.. index:: PCMSolver, PCM

.. _`sec:pcmsolver`:

Interface to PCMSolver
======================

.. codeauthor:: Roberto Di Remigio, T. Daniel Crawford, Andrew C. Simmonett
.. sectionauthor:: Roberto Di Remigio

*Module:* :ref:`Keywords <apdx:pcm>`, :ref:`PSI Variables <apdx:pcm_psivar>`, :source:`PCMSolver <src/lib/libpsipcm>`

|PSIfour| contains code to interface to the PCMSolver library developed
by R. Di Remigio and L. Frediani.
The version 1.1.0 of the PCMSolver library is distributed with |Psifour|
and requires no additional licence, downloads, or configuration.
The library is documented at `http://pcmsolver.readthedocs.org/
<http://pcmsolver.readthedocs.org/en/latest/>`_, while the source code is hosted at
`https://github.com/PCMSolver/pcmsolver/ <https://github.com/PCMSolver/pcmsolver>`_
The library allows for calculations in solution with the polarizable continuum model (PCM),
a continuum solvation model.
Compilation of the library and its interface to |Psifour| can be *disabled* by passing the
``--pcmsolver=off`` to the ``setup`` script or ``-DENABLE_PCMSOLVER=OFF`` directly to CMake.

.. index:: PCM; Using PCM

.. _`sec:usingPCM`:

Using the polarizable continuum model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The inclusion of a PCM description of the solvent into your calculation
is achieved by setting ``pcm true`` in your input file.
|Psifour| understands the additional option ``pcm_scf_type`` with possible values ``total``
(the default) or ``separate``.
The latter forces the separate handling of nuclear and electronic electrostatic potentials and
polarization charges. It is mainly useful for debugging.

.. note:: At present PCM can only be used for energy calculations with SCF wavefunctions.
   Moreover, the PCMSolver library **cannot** exploit molecular point group symmetry.

The PCM model and molecular cavity are specified in a ``pcm`` section that has
to be explicitly typed in by the user. This additional section follows a syntax
that is slightly different from that of |Psifour| and is fully documented
`here <http://pcmsolver.readthedocs.org/en/latest/users/input.html>`_

A typical input for a Hartree-Fock calculation with PCM would look like the following: ::

    molecule NH3 {
    symmetry c1
    N     -0.0000000001    -0.1040380466      0.0000000000
    H     -0.9015844116     0.4818470201     -1.5615900098
    H     -0.9015844116     0.4818470201      1.5615900098
    H      1.8031688251     0.4818470204      0.0000000000
    units bohr
    no_reorient
    no_com
    }

    set {
      basis STO-3G
      scf_type pk
      pcm true
      pcm_scf_type total
    }

    pcm = {
       Units = Angstrom
       Medium {
       SolverType = IEFPCM
       Solvent = Water
       }

       Cavity {
       RadiiSet = UFF
       Type = GePol
       Scaling = False
       Area = 0.3
       Mode = Implicit
       }
    }

More examples can be found in the directories with PCM tests
:source:`pcm_scf <samples/pcm_scf/input.dat>`, :source:`pcm_dft <samples/pcm_dft/input.dat>`
and :source:`pcm_dipole <samples/pcm_dipole/input.dat>`
