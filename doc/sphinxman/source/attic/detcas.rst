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

.. _`sec:casscf`:

Complete-Active-Space Self-Consistent-Field (CASSCF)
====================================================

Multi-configurational self-consistent-field (MCSCF) 
is a general method for obtaining qualitatively correct
wavefunctions for highly strained molecules, diradicals, or bond
breaking reactions.  The most commonly used MCSCF procedure
is the complete-active-space self-consistent-field (CASSCF)
approach [Roos:1980]_, which includes all possible determinants
(with the proper symmetry) that can be formed by distributing 
a set of active electrons among a set of active orbitals.
The detcasman module performs
CASSCF optimization of molecular orbitals via a two-step
procedure in which the CI wavefunction is computed using
detci, and the orbital rotation step is computed using
detcas.  The detcas program is fairly simple
and uses an approximate orbital Hessian [Chaban:1997:88]_
and a Newton-Raphson update,
accelerated by Pulay's DIIS procedure [Pulay:1980]_.
We have also implemented a prototype version of the RASSCF method
[Malmqvist:1990:RASSCF]_, which is another kind of MCSCF which 
is typically less complete (and less expensive) than CASSCF.
However, orbital convergence for RASSCF can be difficult in our
current implementation.

Inactive orbitals in the MCSCF may be specified by the 
|detci__restricted_docc| and |detci__restricted_uocc| keywords. These
orbitals will remain doubly-occupied or doubly-unoccupied, respectively,
in the MCSCF wavefunction.  However, the form of these orbitals will
be optimized in the MCSCF procedure.  It is also possible to 
literally freeze inactive orbitals in their original (SCF) form
using the |globals__frozen_docc| and |globals__frozen_uocc| keywords.
This is not normally what one wishes to do in an MCSCF computation
(*e.g.*, it complicates the computation of gradients),
but it can make the computations faster and is helpful in some
circumstances where unphysical mixing of inactive and active
occupied orbitals might occur.  Presently, it is not possible
to mix the use of restricted and frozen orbitals in |PSIfour|.

The division of the molecular orbitals into various subspaces such as RAS
spaces, or frozen vs active orbitals, etc, needs to be clear not only to
the detci program, but also at least to the transformation program
(and in the case of MCSCF, to other programs as well).  Thus,
orbital subspace keywords such as |detci__ras1|, |detci__ras2|, |detci__ras3|, |detci__ras4|,
|globals__frozen_docc|, |globals__frozen_uocc|, |detci__active|, *etc.*, need to be
in the global ``set {...}`` block section of the input file so they may
be accessed by other modules.

The ability to perform state-averaged 
[Docken:1972:4928]_ [Ruedenberg:1979:1069]_
CASSCF or RASSCF computations has been added.  This is accomplished using the 
|detci__avg_states| keyword.

See the :srcsample:`casscf-sp` and :srcsample:`casscf-sa-sp` examples in the 
samples directory and the example below.

Basic Keywords
--------------

WFN = string
This may be ``casscf`` or ``rasscf``.

REFERENCE = string
Any of the references allowed by detci should work (*i.e.*, not
{\tt uhf}), but there should be no reason not to use {\tt rhf}.

DERTYPE = string
At present, only energies ({\tt none}) are supported; future
releases will implement gradients ({\tt first}).

CONVERGENCE = integer
Convergence desired on the orbital gradient.  Convergence is achieved when
the RMS of the error in the orbital gradient is less than 10**(-n).  The 
default is 4 for energy calculations and 7 for gradients.  Note that
this is a different convergence criterion than for the \PSIdetci\
program itself.  These can be differentiated, if changed by the user,
by placing the {\tt CONVERGENCE} keywords within separate sections of
input, such as {\tt detcas: ( convergence = x )}.

ENERGY\_CONVERGENCE = integer
Convergence desired on the total MCSCF energy.  The default is 7.
\item[RESTRICTED\_DOCC = (integer array)]\mbox{}\\
Should be in {\tt psi:()} or {\tt default:()} sections of input.
The number of lowest energy doubly occupied orbitals in each irreducible
representation from which there will be no excitations.  
These orbitals are optimized in the MCSCF.
The Cotton ordering of the irredicible representations is used.
The default is the zero vector.

RESTRICTED\_UOCC = (integer array)
Should be in {\tt psi:()} or {\tt default:()} sections of input.
The number of highest energy unoccupied orbitals in each irreducible
representation into which there will be no excitations.
These orbitals are optimized in the MCSCF.
The default is the zero vector.

FROZEN\_DOCC = (integer array)]\mbox{}\\
Should be in {\tt psi:()} or {\tt default:()} sections of input.
The number of lowest energy doubly occupied orbitals in each irreducible
representation from which there will be no excitations.  
These orbitals are literally frozen and are not optimized in the MCSCF;
usually one wishes to use {\tt RESTRICTED\_DOCC} instead.
The current version of the program does not allow both
{\tt RESTRICTED\_DOCC} and {\tt FROZEN\_DOCC}.
Should be in {\tt psi:()} or {\tt default:()} sections of input.
The Cotton ordering of the irredicible representations is used.
The default is the zero vector.

FROZEN\_UOCC = (integer array)]\mbox{}\\
Should be in {\tt psi:()} or {\tt default:()} sections of input.
The number of highest energy unoccupied orbitals in each irreducible
representation into which there will be no excitations.
These orbitals are literally frozen and are not optimized in the MCSCF;
usually one wishes to use {\tt RESTRICTED\_UOCC} instead.
The current version of the program does not allow both
{\tt RESTRICTED\_UOCC} and {\tt FROZEN\_UOCC}.
Should be in {\tt psi:()} or {\tt default:()} sections of input.
The default is the zero vector.

NCASITER = integer]\mbox{}\\
Maximum number of iterations to optimize the orbitals.  This option
should be specified in the DEFAULT section of input, because
it needs to be visible to the control program PSI.  Defaults to 20.

AVERAGE\_STATES = (integer array)]\mbox{}\\
This gives a list of what states to average for the orbital 
optimization.  States are numbered starting from 1.

PRINT = integer]\mbox{}\\
This option determines the verbosity of the output.  A value of 1 or
2 specifies minimal printing, a value of 3 specifies verbose printing.
Values of 4 or 5 are used for debugging.  Do not use level 5 unless
the test case is very small (e.g. STO H\ :sub:`2`\ O CISD).

Examples
--------

Example of a CASSCF single-point calculation for H\ :sub:`2`\ O using
a valence active space 3a\ :sub:`1` 1b\ :sub:`1` 2b\ :sub:`2`. ::

    % 6-31G** H2O Test CASSCF Energy Point
                                                                                    
    psi: (
      label = "6-31G** CASSCF H2O"
      jobtype = sp
      wfn = casscf
      reference = rhf
      restricted_docc = (1 0 0 0)
      active          = (3 0 1 2)
      basis = "6-31G**"
      zmat = (
        o
        h 1 1.00
        h 1 1.00 2 103.1
      )
    )

