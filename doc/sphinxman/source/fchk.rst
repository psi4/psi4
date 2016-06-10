
.. include:: autodoc_abbr_options_c.rst

.. index:: 
   FCHK
   Gaussian Formatted Checkpoint

.. _`sec:fchk`:

Interface to programs through FCHK files, :py:func:`~driver.fchk()`
===================================================================

.. codeauthor:: Andrew C. Simmonett
.. sectionauthor:: Andrew C. Simmonett

Many post-processing tools can read information from `Gaussian's formatted
checkpoint (FCHK) files <http://www.gaussian.com/g_tech/g_ur/f_formchk.htm>`_.
To allow interoperability with such tools, |PSIfour| includes a utility to
generate FCHK files.  Wavefunction information, such as orbitals, densities,
orbital energies and basis set information is currently supported, but geometry
optimization and vibrational frequency information are not supported at this
time.  To generate a FCHK file, simply store the wavefunction from the energy
calculation, and use it to create an FCHK writer::

    energy, wfn = energy('scf', return_wfn=True)
    fchk_writer = psi4.FCHKWriter(wfn)
    fchk_writer.write('output.fchk')

The file will be written to the name passed to the FCHK writer's *write()*
method.  Note that for MP2 and CCSD methods, the energy can be computed without
the expensive steps required to compute the density, so energy calls for these
methods will return a wavefunction that has the Hartree-Fock density.  If a
density is required for these methods, the user should instead request a
gradient computation, to ensure that the density is updated appropriately::

    grad, wfn = gradient('mp2', return_wfn=True)
    fchk_writer = psi4.FCHKWriter(wfn)
    fchk_writer.write('output.fchk')

.. autofunction:: driver.fchk(wfn, filename)

