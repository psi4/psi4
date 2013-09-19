
.. include:: autodoc_abbr_options_c.rst

.. index:: Cfour
.. _`sec:cfour`:

Interface to CFOUR by J. Stanton and J. Gauss
=============================================

.. codeauthor:: Lori A. Burns
.. sectionauthor:: Lori A. Burns

*Module:* :ref:`Keywords <apdx:cfour>`, :ref:`PSI Variables <apdx:cfour_psivar>`, :source:`CFOUR <src/bin/cfour>`

|PSIfour| contains code to interface to the Cfour quantum chemistry suite of
John F. Stanton (U. Texas, Austin) and J\ |u_dots|\ rgen Gauss (U. Mainz),
which is available after a license agreement from 
`http://www.cfour.de/ <http://www.cfour.de/>`_.

Installation
~~~~~~~~~~~~

Follow the instructions provided with the Cfour download to install the
executable or to build the source. To by used by |PSIfour|, the program
binary (``xcfour``) must be found in your :envvar:`PATH` or
:envvar:`PSIPATH`.  The ``GENBAS`` file containing basis sets in Cfour
format is not necessary for this interface, but if you prefer to access
basis sets the "Cfour way", ``GENBAS``, too, must be in :envvar:`PATH` or
:envvar:`PSIPATH`. If |PSIfour| is unable to execute the binary, an error
will be reported.

NOTE: Need to check in a GENBAS so tests can run?

Cfour for |PSIfour| Users
~~~~~~~~~~~~~~~~~~~~~~~~~

* Set memory as usual

* Set molecule as usual

* Set basis set as usual (Cfour only cares about orbital basis, no fitting
  bases)

* For the type of computation intended, find appropriate options at
  :ref:`Keywords <apdx:cfour>` which contains the same information as the
  `proper CFOUR options list
  <http://slater.chemie.uni-mainz.de/cfour/index.php?n=Main.ListOfKeywordsInAlphabeticalOrder>`_
  plus notes on keyword relevance when run through |PSIfour|.  Information
  at the `CFOUR manual
  <http://slater.chemie.uni-mainz.de/cfour/index.php?n=Main.Manual>`_ may
  also be useful, as may the many samples at :source:`samples/cfour`.

* Generally, the p4c4 interface will handle best practices for path of
  execution: ``vcc``/``ecc``, derivative type, *etc.* Whereas the user is
  still responsible for setting convergence, frozen core, guess, diis,
  *etc.*

* Set Cfour keywords just like |PSIfour| keywords. The names of keywords
  are unchanged beyond a prepended "cfour\_". (Though be aware that common
  abbreviations like CALC and REF must be fully spelled out when used in
  |PSIfour|.)

* Set the task as usual, indicating Cfour as the intended code by
  prepending "c4-" to the method argument. So ``energy('scf')`` becomes
  ``energy('c4-scf')`` and ``optimize('ccsd(t)')`` becomes
  ``optimize('c4-ccsd(t)')``.




|PSIfour| for Cfour Users
~~~~~~~~~~~~~~~~~~~~~~~~~

friendly abbreviations CALC, ANHARM, CONV don't work if they're not the full keyword name

Output
~~~~~~
