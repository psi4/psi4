
.. include:: autodoc_abbr_options_c.rst

.. index:: 
   Molden
   WebMO
   visualization

.. _`sec:molden`:

Interface to Molden 
==========================================

.. codeauthor:: Justin M. Turney
.. sectionauthor:: C. David Sherrill

|PSIfour| contains an interface to the Molden program.  Molden is a 
visualization program for electronic structure developed by Gijs Schaftenaar
at the University of of Nijmegen, Netherlands.  It is available at 
`http://www.cmbi.ru.nl/molden/ <http://www.cmbi.ru.nl/molden/>`_.  Molden can
plot atomic orbitals, densities, electrostatic potentials (ESP's), etc.  
|PSIfour| can create a file containing
atomic coordinates, basis set, and SCF orbital coefficients in the 
so-called Molden format.  This file is
written by the SCF module (see Section :ref:`SCF <sec:scf>`) 
if the user provides a 
name for the file via the |scf__molden_file| keyword.  This Molden file is 
also used to pass information between |PSIfour| and WebMO, if |PSIfour| 
computations are invoked using the WebMO GUI (|PSIfour| input files created 
by WebMO automatically set |scf__molden_file| to `output.molden`.


