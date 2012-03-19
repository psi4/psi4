
.. _`sec:psimod`:

==============================
PsiMod: Linking C++ and Python
==============================

.. literalinclude:: man_PsiMod.rst 

..  * NOTES (LAB 2-28-2012)
    * How to get PsiMod.py documentation into Sphinx
    * (1) In psi4/src/bin/psi4/python.cc , uncomment the following lines
          s = strdup("help(PsiMod)");
          PyRun_SimpleString(s);
    * (2) cd into objdir and make psi4
    * (3) cd into a testcase directory
          cd tests/tu1-h2o-energy/
          make clean
          make > $PSIDATADIR/../doc/sphinxman/source/man_PsiMod.rst 
    * (4) Edit man_PsiMod.rst to remove junk at end of file.
    * (5) Rebuild Sphinx docs.
    * (6) LaTeX build can't handle including man_PsiMod.rst so comment out by default.
    * (7) Rebuild psi4 with lines in (1) commented.

.. toctree::
   
   man_PsiMod

