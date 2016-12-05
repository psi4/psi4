
=====================
Compiling and Testing
=====================


.. _`faq:buildquick`:

How to build and install |PSIfour|, the compact version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section outlines the main steps of configuring, compiling, and
installing |PSIfour|. More detail is given :ref:`here
<faq:builddetailed>`. ::

    >>> cd {top-level-psi4-dir}
    >>> cmake -H. -Bobjdir [your configuration options]
    >>> cd objdir
    >>> make -j`getconf _NPROCESSORS_ONLN`
    >>> make install


.. _`faq:builddetailed`:

How to build, test, and install |PSIfour|, in detail
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**1. Plan Directories**

   Get ahold of the |PSIfour| codebase, and navigate to the top level source
   directory, hereafter :samp:`{top-level-psi4-dir}`.
   
   * :ref:`How to obtain Psi4: start with find-the-code quiz, end in {top-level-psi4-dir} <faq:obtainpsi4>`
   
   ::
   
    >>> cd {top-level-psi4-dir}
   
   Choose a compilation directory, hereafter :samp:`{objdir}`
   
   * :ref:`How to compile elsewhere than {top-level-psi4-dir/objdir} <faq:setupobjdir>`
   * :ref:`How to choose the compilation directory, {objdir} <faq:chooseobjdir>`
   
   Choose an installation directory, hereafter :samp:`{prefix}`
   
   * :ref:`How to install elsewhere than "/usr/local/psi4" <faq:setupprefix>`

**2. Plan Configuration**

   Examine the strict and optional software requirements to make sure the
   target computer has all the necessary dependencies installed.

   * :ref:`What are the tools and dependencies strictly required for building Psi4 <faq:coredepend>`
   * :ref:`What are the add-on capabilities for Psi4 and what are their dependencies <faq:addondepend>`

   Prepare any necessary or desired configuration options for ``cmake``,
   hereafter ``[your configuration options]``

   * :ref:`How to see what build configuration options are available <faq:setuphelp>`
   * :ref:`How to configure Psi4 for invoke CMake <faq:cmakeviasetup>`

**3. Configure**

   Run CMake with planned options and directories, as below. It reports on
   software found or unfound as it scans the computer, then (upon success)
   creates :samp:`{objdir}` ready for compilation.

   ::

    >>> cmake -H. -B{objdir} -DCMAKE_INSTALL_PREFIX={prefix} [your configuration options]

**4. Compile**

   Compile the code (optional ``-j`` triggers parallel compilation).

   ::

    >>> cd {objdir}
    >>> make -j`getconf _NPROCESSORS_ONLN`

**5. Test**

   Optionally, use [CTest](9_CMake) to test the build.

   * :ref:`How to run a minute's worth of tests <faq:minutetests>`
   * :ref:`How to run a subset of tests <faq:subsettests>`
   * :ref:`How to see CTest testing errors <faq:testsoutput>`

   ::

   >>> ctest -j`getconf _NPROCESSORS_ONLN`

**6. Install**

   If tests pass, install the code.

   ::

   >>> make install

**7. Configure Runtime**

   To run Psi4 after installation, you need to configure a few variables:

   * :ref:`How to run installed Psi4 <faq:runfromprefix>`
   * :ref:`How to run Psi4 from the compilation directory <faq:runfromobjdir>`





.. _`faq:recompile`:

How to update and rebuild |PSIfour|
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Obtain code updates as appropriate from LINKTOVARMODIES. Move into
:samp:`{objdir}` and reissue ``make``, whereupon CMake may reconfigure but
will only rebuild objects and libraries depending on changed files.  It is
scarcely ever necessary for the user to reinvoke ``cmake`` to update
:samp:`{objdir}`.


.. _`faq:minutetests`:

How to run a minute's worth of tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When you want to do a very minimal test of the build and have
[CTest](9_CMake) installed, the following command can be useful. ::

    >>> ctest -L mini -j`getconf _NPROCESSORS_ONLN`


.. _`faq:subsettests`:

How to run a subset of tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[CTest](9_CMake) allows flexibly partitioned running of the test suite. In
the examples below, *testname* are regex of [test
names](../blob/master/tests), and *testlabel* are regex of labels (*e.g.*,
``cc``, ``mints``, ``libefp``).

* Run tests in parallel with ``-j`` flag. For maximum parallelism: :samp:`make -j\`getconf _NPROCESSORS_ONLN\`\ `
* Run about a third of the tests in 10--20 minutes, the so-called *quicktests*: ``ctest -L quicktests``
* Run tests matching by name: ``ctest -R testname``
* Run tests excluding those by name: ``ctest -E testname``
* Run tests matching by label: ``ctest -L testlabel``
* Run tests excluding those by label: ``ctest -LE testlabel``


.. _`faq:testsoutput`:

How to see ``ctest`` testing errors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

 >>> ctest
 Test project /your/path/2/psi4/build/directory/tests
     Start 248: tu1-h2o-energy
 1/2 Test #248: tu1-h2o-energy ...................   Passed    1.73 sec
      Start  6: cc1
 2/2  Test  #6: cc1 ..............................***Failed    0.07 sec
 ...

When ``ctest`` reports that some (or all) tests have failed, look in your
build directory for file
:samp:`{objdir}/tests/Testing/Temporary/LastTest.log`. It may have a
``.tmp`` extension, depending on whether the last test was interrupted and
a few other factors. Either way, this file should contain CMake's testing
output, as well as everything that was printed to the screen.

