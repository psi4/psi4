
.. include:: autodoc_abbr_options_c.rst

.. index:: 
   DMA
   GDMA
   Distributed Multipole Analysis

.. _`sec:gdma`:

Interface to GDMA Distributed Multipole Analysis by A. J. Stone, :py:func:`~driver.gdma`
========================================================================================

.. codeauthor:: Anthony J. Stone, Andrew C. Simmonett
.. sectionauthor:: Andrew C. Simmonett

*Module:* :ref:`Keywords <apdx:gdma>`, :ref:`PSI Variables <apdx:gdma_psivar>`, :source:`PCMSolver <src/lib/libgdma>`

Input
~~~

The distributed multipole analysis (DMA) technique, developed by Anthony J.
Stone and implemented by him into the `GDMA package
<http://www-stone.ch.cam.ac.uk/programs.html>`_, is available in |PSIfour|.
The current implementation simply embeds Stone's GDMA code into the main
executable, and generates the appropriate input files automatically.  The
program takes as input a data file, and a Gaussian formatted checkpoint (see
Section :ref:`FCHK <sec:fchk>`) file.  The simplest usage of the GDMA code is
demonstrated below, along with a listing of the options supported; these
options correspond to the options described in the 
:download:`GDMA manual <gdma-2.2.06.pdf>`.

If more advanced usage is desired, which is not is permitted by the options
listed below, the user may provide their own data file containing keywords to
control the GDMA code.  Simply place the data file in the directory |PSIfour|
is called from, and provide the file name as the datafile argument to the
:py:func:`~driver.gdma` routine.  For example, if GDMA data file is called
*control.dma*, the GDMA code is called as follows::

    grad, wfn = gradient('mp2', return_wfn=True)
    gdma(wfn, datafile='control.dma')

An FCHK file will be generated for the GDMA code to read; this file will have
the prefix given by |globals__writer_file_label| (if set), or else by the name
of the output file plus the name of the current molecule, and the suffix will
be '.fchk'.  This FCHK file name should be passed to the 'File' keyword in the
DGMA data file, to ensure that the GDMA code reads the correct wavefunction
information.

After running, two matrices of results can be accessed::

    dma_results = get_array_variable('DMA DISTRIBUTED MULTIPOLES')
    tot_results = get_array_variable('DMA TOTAL MULTIPOLES')

The first contains distributed multipoles, in units given by
|gdma__gdma_multipole_units|, with the row index corresponding to the site and
the column index referencing the multipole component.  Both indices are zero
based, and the :math:`Q^l_m` components of the multipoles are ordered as
:math:`Q^0_0, Q^1_0, Q^1_{1c}, Q^1_{1s}, Q^2_0, Q^2_{1c}, Q^2_{1s}, Q^2_{2c},
Q^2_{2s}, \ldots`  The second matrix returned has a single row, whose columns
are the total multipoles, translated to |gdma__gdma_origin|, and summed.


.. autofunction:: driver.gdma(wfn)

Options
~~~~~~~

.. include:: autodir_options_c/gdma__gdma_limit.rst
.. include:: autodir_options_c/gdma__gdma_origin.rst
.. include:: autodir_options_c/gdma__gdma_multipole_units.rst
.. include:: autodir_options_c/gdma__gdma_radius.rst
.. include:: autodir_options_c/gdma__gdma_switch.rst
