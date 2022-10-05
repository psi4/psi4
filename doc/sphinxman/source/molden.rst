.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2022 The Psi4 Developers.
.. #
.. # The copyrights for code used from other parties are included in
.. # the corresponding files.
.. #
.. # This file is part of Psi4.
.. #
.. # Psi4 is free software; you can redistribute it and/or modify
.. # it under the terms of the GNU Lesser General Public License as published by
.. # the Free Software Foundation, version 3.
.. #
.. # Psi4 is distributed in the hope that it will be useful,
.. # but WITHOUT ANY WARRANTY; without even the implied warranty of
.. # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. # GNU Lesser General Public License for more details.
.. #
.. # You should have received a copy of the GNU Lesser General Public License along
.. # with Psi4; if not, write to the Free Software Foundation, Inc.,
.. # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
.. #
.. # @END LICENSE
.. #

.. include:: autodoc_abbr_options_c.rst

.. index:: 
   Molden
   WebMO
   visualization

.. _`sec:molden`:

Interface to Molden |w---w| :py:func:`~psi4.driver.molden`
==========================================================

.. codeauthor:: Justin M. Turney
.. sectionauthor:: C. David Sherrill

|PSIfour| contains an interface to the Molden program.  Molden is a 
visualization program for electronic structure developed by Gijs Schaftenaar
at the University of of Nijmegen, Netherlands.  It is available at 
https://www3.cmbi.umcn.nl/molden/ . Molden can
plot atomic orbitals, densities, electrostatic potentials (ESPs), etc.
|PSIfour| can create a file containing
atomic coordinates, basis set, and SCF orbital coefficients in the 
so-called Molden format.  This file is
written by the SCF module (see Section :ref:`SCF <sec:scf>`) 
if the user sets the |scf__molden_write| keyword to true.  This Molden file is 
also used to pass information between |PSIfour| and WebMO, if |PSIfour| 
computations are invoked using the WebMO GUI.  The filename of the 
Molden file ends in ".molden", and the prefix is determined by 
|globals__writer_file_label| (if set), or else by the name of the output
file plus the name of the current molecule. If |globals__molden_with_virtual|
is set to false, the unoccupied orbitals are not written to the Molden
file.

.. autofunction:: psi4.molden(wfn, filename)
   :noindex:

Options
~~~~~~~

.. include:: autodir_options_c/scf__molden_write.rst
.. include:: autodir_options_c/globals__writer_file_label.rst
.. include:: autodir_options_c/globals__molden_with_virtual.rst

