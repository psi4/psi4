
.. _`sec:psithonFunc`:

=========================================
Psithon Functions: Invoking a Calculation
=========================================

To allow arbitrarily complex computations to be performed, PSI4 is built
upon the Python interpreter, with modifications termed Psithon. Sec. 
:ref:`sec:psithonInput` describes the non-standard Python associated with
clean molecule, basis, and option specification in the PSI4 input file.
This documentation addresses the pure Python side- what functions allow
the efficient compiled code to be run, what functions post-process and
interact with that output, and how the ordinary (or ambitious) user can
extent PSI4's functionality.

.. toctree::
   :maxdepth: 2

   notes_py
   energy
   cp
   opt
   freq
   db
   cbs
   intercalls
   sowreap

