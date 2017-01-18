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

.. _`sec:psithonFunc`:

=========================================
Psithon Functions: Invoking a Calculation
=========================================

To allow arbitrarily complex computations to be performed, |PSIfour| is built
upon the Python interpreter, with modifications termed Psithon. Sec. 
:ref:`sec:psithonInput` describes the non-standard Python associated with
clean molecule, basis, and option specification in the |PSIfour| input file.
This documentation addresses the pure Python side- what functions allow
the efficient compiled code to be run, what functions post-process and
interact with that output, and how the ordinary (or ambitious) user can
extent |PSIfours| functionality.

.. toctree::
   :maxdepth: 2

   notes_py
   energy
   prop
   nbody
   opt
   freq
   db
   cbs
   diatomic
   intercalls
   sowreap
   cubeprop

