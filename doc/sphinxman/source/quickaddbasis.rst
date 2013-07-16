
.. include:: autodoc_abbr_options_c.rst

.. index::
   pair: basis set; adding new

.. note:: No recompile of the PSI program is necessary for changes made to
    files in ``$PSIDATADIR``, including those described below.

.. _`sec:basisUserDefined`:

User-Defined Basis Sets
=======================

There are three routes by which a basis set in G94 format can be introduced to |PSIfours| notice.


.. rubric:: (1) Install new basis set file into |PSIfour| basis library.

Copy the basis set definitions for all elements into a blank file. Exclamation points denote comments.
As the first line of the file, add the word ``spherical`` or ``cartesian`` to indicate
whether the basis set will run in (5D/7F) or (6D/10F). ::

   cartesian
   ****
   H     0
   S   3   1.00
         3.42525091             0.15432897
         0.62391373             0.53532814
         0.16885540             0.44463454
   ****
   O     0
   S   3   1.00
       130.7093200              0.15432897
        23.8088610              0.53532814
         6.4436083              0.44463454
   SP   3   1.00
         5.0331513             -0.09996723             0.15591627
         1.1695961              0.39951283             0.60768372
         0.3803890              0.70011547             0.39195739
   ****

Name the file with the name of the basis set and a ``.gbs`` extension,
after applying the following transformations.

* All letters lowercase
* Replace all ``*`` with ``s``
* Replace all ``+`` with ``p``
* Replace all ``(`` ``)`` ``,`` with ``_`` (underscores replace parentheses and commas)

For example, basis 6-31++G** is stored in :source:`lib/basis/6-31ppgss.gbs`, 
and cc-pV(D+d)Z is stored in :source:`lib/basis/cc-pv_dpd_z.gbs`.
Only one basis set may be specified per file.
Copy the new basis set file into :source:`lib/basis`.
Request the new basis set in an input file in the usual manner. ::

   set basis new_basis_name


.. rubric:: (2) Use new basis set file in arbitrary location.

Prepare a basis set file exactly as above. Append the directory
containing the basis set file to the environment variable
:envvar:`PSIPATH`.

Request the new basis set in an input file in the usual manner. ::

   set basis new_basis_name

.. rubric:: (3) Include new basis set in input file.

Construct for a basis set a section like the one below that includes
``[basis name]``, |globals__puream| value, and element basis set
specifications. Hash signs denote comments.  This format is exactly like
the stand-alone basis file except for the addition of the basis name in
brackets. ::

   [ sto-3g ]
   cartesian
   ****
   H     0
   S   3   1.00
         3.42525091             0.15432897
         0.62391373             0.53532814
         0.16885540             0.44463454
   ****
   O     0
   S   3   1.00
       130.7093200              0.15432897
        23.8088610              0.53532814
         6.4436083              0.44463454
   SP   3   1.00
         5.0331513             -0.09996723             0.15591627
         1.1695961              0.39951283             0.60768372
         0.3803890              0.70011547             0.39195739
   ****

Copy the section into a |PSIfour| input file and surround it with the
command ``basis {...}``, as shown below.  Multiple basis sets can be
specified by adding additional sections within the surrounding brackets.
Use ``assign`` statements to actually request the basis set. (See
:srcsample:`mints2` for an example.) ::

   basis {

   # assign basset to all atoms and addl to hydrogens
   assign basset
   assign H addl

   # basis set section like in snippet above goes here
   [basset]
   ...

   # additional basis set sections follow
   [addl]
   ...
   }

