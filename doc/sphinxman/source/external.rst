
.. include:: autodoc_abbr_options_c.rst

======================================
Installation and Runtime Configuration
======================================

Obtaining |PSIfour|
===================

The latest version of the |PSIfour| program package may be
obtained at `www.psicode.org <http://www.psicode.org>`_. The
package is available as a binary for Linux (:ref:`Installing
from Binary <sec:conda>`) or as source code (zipped
archive or git repository from `www.github.com/psi4/psi4
<http://www.github.com/psi4/psi4>`_).


Installing from Binary
======================

.. toctree::
   :maxdepth: 2

   conda


.. index:: prerequisites, compiling, installing
.. _`sec:installFile`:

Compiling and Installing from Source
====================================

Detailed directions on 
`obtaining <https://github.com/psi4/psi4public/wiki/1_Obtaining>`_, 
`prerequisites <https://github.com/psi4/psi4public/wiki/2_Planning#-what-are-the-tools-and-dependencies-strictly-required-for-building-psi4>`_,
`building and installing <https://github.com/psi4/psi4public/wiki/3_Building>`_,
and `FAQ <https://github.com/psi4/psi4public/wiki/0_FAQ>`_
are maintained on the `GitHub Wiki <https://github.com/psi4/psi4public/wiki>`_. 
If uncertain, `start here <https://github.com/psi4/psi4public/wiki/1_Obtaining#quiz>`_.


.. index:: scratch files, restart
.. _`sec:Scratch`:

Scratch Files and Elementary Restart
====================================

One very important part of user configuration at the end of the
installation process
is to tell |PSIfour| where to write its temporary
("scratch") files.  Electronic structure packages like |PSIfour| can
create rather large temporary disk files.  It is very important to 
ensure that |PSIfour| is writing its temporary files to a disk drive
physically attached to the computer running the computation.  If it
is not, it will significantly slow down the program and the network.
By default, |PSIfour| will write temporary files to ``/tmp``, but this
directory is often not large enough for typical computations.  Therefore,
you need to (a) make sure there is a sufficiently large directory on a
locally attached disk drive (100GB–1TB or more, depending on the size of
the molecules to be studied) and (b) tell |PSIfour| the path to this
directory. Scratch file location can be specified through the 
:envvar:`PSI_SCRATCH` environment variable or through the |psirc| file
(see section :ref:`sec:psirc`). Most of the time, :envvar:`PSI_SCRATCH`
is preferred, and it overrides any existing |psirc| setting. You can set up 
:envvar:`PSI_SCRATCH` by issuing the following commands in a terminal,
or including them in the appropriate ``rc`` file. For C shell (``~/.tcshrc`` file): ::

    setenv PSI_SCRATCH /scratch/user

For Bash (``~/.bashrc`` file): ::

    export PSI_SCRATCH=/scratch/user

|PSIfour| has a number of utilities that manage
input and output (I/O) of quantities to and from the hard disk.  Most
quantities, such as molecular integrals, are intermediates that are not of
interest to the user and can be deleted after the computation finishes, but
pertinent details of computations are also written to a checkpoint file and
might be useful in subsequent computations.  All files are written to the
designated scratch numbered by :ref:`content <apdx:psiFiles>` and labeled
with the process id, then are deleted at the end of the computation,
unless otherwise instructed by the user.

A Python callable handle to the |PSIfour| I/O management routines is available,
and is called ``psi4_io``.  To instruct the I/O manager to send all files to
another location, say ``/scratch/user``, add the following command to your input
file: ::

    psi4_io.set_default_path('/scratch/user')

For batch jobs running through a queue, it might be more convenient to use an
environmental variable (in this case ``$MYSCRATCH``) to set the scratch directory;
the following code will do that::

    import os
    scratch_dir = os.environ.get('MYSCRATCH')
    if scratch_dir:
        psi4_io.set_default_path(scratch_dir + '/')

Individual files can be sent to specific locations.  For example, file 32 is
the checkpoint file that the user might want to retain in the working directory
(*i.e.*, where |PSIfour| was launched from) for restart purposes.  This is
accomplished by the commands below::

    psi4_io.set_specific_path(32, './')
    psi4_io.set_specific_retention(32, True)

which is equivalent to ::

    psi4_io.set_specific_path(PSIF_CHKPT, './')
    psi4_io.set_specific_retention(PSIF_CHKPT, True)

A guide to the contents of individual scratch files may be found at :ref:`apdx:psiFiles`.
To circumvent difficulties with running multiple jobs in the same scratch, the
process ID (PID) of the |PSIfour| instance is incorporated into the full file
name; therefore, it is safe to use the same scratch directory for calculations
running simultaneously. This also means that if the user wants |PSIfour| to use
information from a previous file, like molecular orbitals, he needs to provide the
name of the file. This can be done through the ``restart_file`` option ::

  energy('scf',restart_file='./psi.PID.name.filenumber')

where by default, PID is the process number, name the name of the molecule,
and filenumber is listed in :ref:`content <apdx:psiFiles>`. Only the filenumber
is necessary for the driver to appropriately rename the file for the next |PSIfour|
job, and if none is found it defaults to 32, a checkpoint file. If two or more files
are to be read, they need to be provided as a Python list ::

  energy('scf',restart_file=['./file1.filenumber','./file2.filenumber'])

Note that the ``restart_file`` options is only available for energy procedures as of now.

Executing |PSIfour| with the :option:`psi4 -m` (for
messy) flag will prevent files being deleted at the end of the run::

    psi4 -m


.. index:: psirc, psi4rc
.. _`sec:psirc`:

|psirc| File
============

If using the environment variable :envvar:`PSI_SCRATCH` is inconvenient,
or if some ``psi4_io`` commands must be present in all input files,
the |psirc| resource file can be used (example :source:`samples/example_psi4rc_file`). 

All the commands mentioned in section :ref:`sec:Scratch` can be used in this file,
namely: ::

    psi4_io.set_default_path('/scratch/user')

to set up the scratch path, ::

    import os
    scratch_dir = os.environ.get('MYSCRATCH')
    if scratch_dir:
        psi4_io.set_default_path(scratch_dir + '/')

to set up the scratch path from a variable ``$MYSCRATCH``, ::

    psi4_io.set_specific_path(32, './')
    psi4_io.set_specific_retention(32, True)

which is equivalent to ::

    psi4_io.set_specific_path(PSIF_CHKPT, './')
    psi4_io.set_specific_retention(PSIF_CHKPT, True)

to set up a specific path for the checkpoint file and instruct |PSIfour| not to delete it.

The Python interpreter will execute the contents of the
|psirc| file in the current user's home area (if present) before performing any
tasks in the input file. As a consequence, the commands in the input files supersede
any instructions in the |psirc| file. During
excecution, the |psirc| defaults will be loaded in first, but then the commands
in the input file will be executed.  

The |psirc| file can also be used to define constants that are accessible
in input files or to place any Python statements that should be executed
with every |PSIfour| instance.

.. index:: parallel operation, threading
.. _`sec:threading`:

Threading
=========

Most new modules in |PSIfour| are designed to run efficiently on SMP architectures
via application of several thread models. The de facto standard for |PSIfour|
involves using threaded BLAS/LAPACK (particularly Intel's excellent MKL package)
for most tensor-like operations, OpenMP for more general operations, and Boost
Threads for some special-case operations. Note: Using OpenMP alone is a really
bad idea. The developers make little to no effort to explicitly parallelize
operations which are already easily threaded by MKL or other threaded BLAS. Less
than 20% of the threaded code in |PSIfour| uses OpenMP, the rest is handled by
parallel DGEMM and other library routines. From this point forward, it is
assumed that you have compiled |PSIfour| with OpenMP and MKL (Note that it is
possible to use g++ or another compiler and yet still link against MKL).

Control of threading in |PSIfour| can be accomplished at a variety of levels,
ranging from global environment variables to direct control of thread count in
the input file, to even directives specific to each model. This hierarchy is
explained below. Note that each deeper level trumps all previous levels.

.. rubric:: (1) OpenMP/MKL Environment Variables

The easiest/least visible way to thread |PSIfour| is to set the standard OpenMP/MKL
environment variables :envvar:`OMP_NUM_THREADS` and :envvar:`MKL_NUM_THREADS`. 
For instance, in tcsh::

    setenv OMP_NUM_THREADS 4
    setenv MKL_NUM_THREADS 4

|PSIfour| then detects these value via the API routines in ``<omp.h>`` and
``<mkl.h>``, and runs all applicable code with 4 threads. These environment
variables are typically defined in a ``.tcshrc`` or ``.bashrc``.

.. rubric:: (2) The -n Command Line Flag

To change the number of threads at runtime, the :option:`psi4 -n` flag may be used. An
example is::

    psi4 -i input.dat -o output.dat -n 4

which will run on four threads.

.. rubric:: (3) Setting Thread Numbers in an Input

For more explicit control, the Process::environment class in |PSIfour| can
override the number of threads set by environment variables. This functionality
is accessed via the :py:func:`~p4util.util.set_num_threads` Psithon function, which controls
both MKL and OpenMP thread numbers. The number of threads may be changed
multiple times in a |PSIfour| input file. An example input for this feature is::

    # A bit small-ish, but you get the idea
    molecule h2o {
    0 1
    O
    H 1 1.0
    H 1 1.0 2 90.0
    }
    
    set scf {
    basis cc-pvdz
    scf_type df
    }

    # Run from 1 to 4 threads, for instance, to record timings
    for nthread in range(1,5):
        set_num_threads(nthread)
        energy('scf')

.. rubric:: (4) Method-Specific Control

Even more control is possible in certain circumstances. For instance, the
threaded generation of AO density-fitted integrals involves a memory requirement
proportional to the number of threads. This requirement may exceed the total
memory of a small-memory node if all threads are involved in the generation of
these integrals. For general DF algorithms, the user may specify::

    set MODULE_NAME df_ints_num_threads n

to explicitly control the number of threads used for integral formation. Setting
this variable to 0 (the default) uses the number of threads specified by the
:py:func:`~p4util.util.set_num_threads` Psithon method or the default environmental variables.

.. index:: PBS queueing system, threading
.. _`sec:PBS`:

PBS job file
============

To run a |PSIfour| job on a PBS queueing system, you need to properly set up
all necessary variables in the PBS job file. Below is a minimal example of
a PBS job file for a threaded job, and a short explanation for each section. ::

    #!/bin/tcsh
    #PBS -j oe
    #PBS -l pmem=2120mb
    #PBS -N jobname
    #PBS -V
    
    
    setenv OMP_NUM_THREADS 4
    setenv MKL_NUM_THREADS 4
    cd $PBS_O_WORKDIR
    setenv myscratch /scratch/user/psi4.$PBS_JOBID
    
    foreach i (`sort $PBS_NODEFILE | uniq`)
        echo "Creating scratch directory " $myscratch " on " $i
        ssh $i rm -rf $myscratch
        ssh $i mkdir -p $myscratch
    end
    
    unsetenv PSI4DATADIR
    unsetenv PSIDATADIR
    setenv PSI_SCRATCH $myscratch
    if ! ( $?PSIPATH ) setenv PSIPATH ""
    setenv PSIPATH /path/to/external/modules:${PSIPATH}
    setenv PSIPATH /path/to/python/modules:${PSIPATH}
    /psi/install/directory/bin/psi4 -i input.in -o input.out
    
    foreach i (`sort $PBS_NODEFILE | uniq`)
        echo "Removing scratch directory " $myscratch " on " $i
        ssh $i rm -rf $myscratch
    end

The top section features PBS-specific commands. These depend on the 
specific characteristics of your PBS queuing system but they may include: ::

    #!/bin/tcsh
    #PBS -j oe 
    #PBS -l pmem=2120mb
    #PBS -N jobname
    #PBS -V
    
The ``PBS -j oe`` option instructs PBS to write any output or error message
from the queuing system in dedicated files. ``PBS -l pmem=2120mb`` requests 
2120 MB of memory for each thread on the node. The total memory requested for 
the job by PBS should generally be slightly greater than what indicated 
in the input file (see :ref:`memory setting <sec:memory>`).

In the next section, we define :envvar:`OMP_NUM_THREADS` and :envvar:`MKL_NUM_THREADS`
to use 4 threads for OpenMP parallelization and in threaded BLAS (see section :ref:`sec:threading`). ::

    setenv OMP_NUM_THREADS 4
    setenv MKL_NUM_THREADS 4

Then, we move to the working directory using PBS variable ``$PBS_O_WORKDIR`` and 
we create scratch directories on every node, using the ``$PBS_NODEFILE`` which 
points to a file containing a list of the nodes attributed to the job. ::

    cd $PBS_O_WORKDIR
    setenv myscratch /scratch/user/psi4.$PBS_JOBID
    
    foreach i (`sort $PBS_NODEFILE | uniq`)
        echo "Creating scratch directory " $myscratch " on " $i
        ssh $i rm -rf $myscratch
        ssh $i mkdir -p $myscratch
    end

The next section is *very important* as it sets the environment variables needed
by |PSIfour|: ::

    unsetenv PSIDATADIR
    setenv PSI_SCRATCH $myscratch
    if ! ( $?PSIPATH ) setenv PSIPATH ""
    setenv PSIPATH /path/to/external/modules:${PSIPATH}
    setenv PSIPATH /path/to/python/modules:${PSIPATH}

:envvar:`PSIDATADIR` does *not* need to be set if |PSIfour| has been *properly installed*.
In the present example we unset it to make sure it does not interfere with the location
of the installed directory. :envvar:`PSIPATH` is needed only if you are using external modules or 
plugins in |PSIfour| and should point to the directories where they can be found. In the
present example, we make sure the variable is set with ``if ! ( $?PSIPATH ) setenv PSIPATH ""``
before adding more paths to it. Finally, :envvar:`PSI_SCRATCH` should point to a fast, 
local disk for temporary file storage. The next step is then to actually run the computation: ::

    /psi/install/directory/bin/psi4 -i input.in -o input.out

And then to clean up the scratch directories previously created: ::

    foreach i (`sort $PBS_NODEFILE | uniq`)
        echo "Removing scratch directory " $myscratch " on " $i
        ssh $i rm -rf $myscratch
    end

Note again that the specific commands for your PBS system may differ. Refer
to your system administrator.

.. _`sec:commandLineOptions`:

Command Line Options
====================

|PSIfour| can be invoked with no command line arguments, as it takes as input
by default the file "input.dat" and directs output by default to "output.dat".
The set of three commands below are completely equivalent, while the fourth is,
perhaps, the most common usage. ::

   >>> psi4
   >>> psi4 -i input.dat -o output.dat
   >>> psi4 input.dat output.dat

   >>> psi4 descriptive_filename.in descriptive_filename.out

Command-line arguments to |PSIfour| can be accessed through :option:`psi4 --help`.

.. program:: psi4

.. option:: -a, --append

   Append results to output file. Default: Truncate first

.. option:: -d, --debug

   Flush the outfile at every fprintf. Default: true iff ``--with-debug``

.. option:: -h, --help

   Display the command-line options and usage information.

.. option:: -i <filename>, --input <filename>

   Input file name. Default: input.dat

.. option:: -l <name>, --psidatadir <name>

   Mainly for use by developers, this overrides the value of
   :envvar:`PSIDATADIR` and specifies the path to the Psi data
   library (psi4/share) 

.. option:: -m, --messy

   Leave temporary files after the run is completed.

.. option:: -n <threads>, --nthread <threads>

   Number of threads to use (overrides :envvar:`OMP_NUM_THREADS`)

.. option:: -o <filename>, --output <filename>

   Output file name. Use ``stdout`` as <filename> to redirect 
   to the screen. Default: when the input filename is "input.dat",
   then the output filename defaults to "output.dat".  Otherwise, the
   output filename defaults to the the input filename (subtracting
   any ".in" or ".dat" suffix) plus ".out"

.. option:: -p <prefix>, --prefix <prefix>

   Prefix for psi files. Default: psi

.. option:: -s <name>, --scratch <name>

   This overrides the value of :envvar:`PSI_SCRATCH` and provides
   a path to the location of scratch files

.. option:: --new-plugin <name>

   Creates a new directory <name> with files for writing a
   new plugin. An additional argument specifies a template
   to use, for example: ``--new-plugin name +mointegrals``.
   See Sec. :ref:`sec:plugins` for available templates.

.. option:: -v, --verbose

   Print a lot of information, including the Psithon translation of the input file

.. option:: -V, --version

   Print version information. ::

     >>> psi4 --version
     0.4.262

.. option:: -w, --wipe

   Clean out scratch area.


.. _`sec:environmentVariables`:

Environment Variables
=====================

These environment variables will influence |PSIfours| behavior.

.. envvar:: MKL_NUM_THREADS

   Number of threads to use by operations with Intel threaded BLAS libraries.

.. envvar:: OMP_NESTED

   Do access nested DGEMM in OpenMP sections in DFMP2 for multi-socket
   platforms. This is very low-level access to OpenMP functions for
   experienced programmers. Users should leave this variable unset or set
   to ``False``.

.. envvar:: OMP_NUM_THREADS

   Number of threads to use by modules with OpenMP threading.

.. envvar:: PATH

   Path for interfaced executables. 

   .. note:: Configuring |PSIfour| through :envvar:`PSIPATH` is preferred
      to modifying this environment variable.

   To run K\ |a_acute|\ llay's MRCC program 
   (see :ref:`MRCC <sec:mrcc>`), the ``dmrcc`` executable must be in :envvar:`PATH`.
   Likewise to run Grimme's dftd3 program (see :ref:`dftd3 <sec:dftd3>`), the 
   ``dftd3`` executable must be in :envvar:`PATH`.

.. envvar:: PSI_SCRATCH

   Directory where scratch files are written. Overrides settings in |psirc|.
   It is very important to ensure that |PSIfour| is writing its scratch files 
   to a disk drive physically attached to the computer running the computation. 
   If it is not, it will significantly slow down the program and the network. 

   Modify :envvar:`PSI_SCRATCH` through normal Linux shell commands before invoking ``psi4`` ::

      # csh, tcsh
      >>> setenv PSI_SCRATCH /scratch/user

      # bash
      >>> export PSI_SCRATCH=/scratch/user

   You can also include the above commands in the respective ``rc`` file, i.e.
   ``~/.tcshrc`` for csh and tcsh or ``~/.bashrc`` for Bash.

.. envvar:: PSIPATH

   Path in which |PSIfour| looks for user extensions to the built-in
   libraries. Specifically, directories containing 
   :ref:`user basis sets <sec:basisUserDefined>`,
   :ref:`EFP fragments <sec:findingEFPFragments>`,
   :ref:`databases <sec:createDatabase>`, 
   :ref:`plugins <sec:plugins>`, and 
   interfaced executables (
   ``dmrcc`` for :ref:`MRCC <sec:mrcc>` and 
   ``dftd3`` for :ref:`DFTD3 <sec:dftd3>`
   ) should be placed in this colon-separated list.

   |PSIfour| is designed so that user extensions that are findable through
   :envvar:`PSIPATH` can be used in input files entirely like their
   built-in counterparts, without additional tagging as non-standard.

   The typical search path is first the built-in libraries, next each
   :envvar:`PSIPATH` directory in order, and finally the execution
   directory (I won't swear everything tacks on the execution directory).

   Path in which the Python interpreter looks for modules to import. For 
   |PSIfour|, these are generally plugins (see :ref:`sec:plugins`) or databases.

   Modify :envvar:`PSIPATH` through normal Linux shell commands before invoking ``psi4`` ::

      # csh, tcsh
      >>> setenv PSIPATH /home/user/psiadditions:/home/user/gbs

      # bash
      >>> export PSIPATH=/home/user/psiadditions:/home/user/gbs

.. envvar:: PYTHONPATH

   Path in which the Python interpreter looks for modules to import. For 
   |PSIfour|, these are generally plugins (see :ref:`sec:plugins`) or databases.

   .. note:: Configuring |PSIfour| through :envvar:`PSIPATH` is preferred
      to modifying this environment variable.

   Modification of :envvar:`PYTHONPATH` can be done in three ways, equivalently.

   * Normal Linux shell commands. ::

        # csh/tcsh
        setenv PYTHONPATH /home/user/psiadditions:$PYTHONPATH
        # sh/bash
        PYTHONPATH=/home/user/psiadditions:$PYTHONPATH; export PYTHONPATH

   * Place the path in the |psirc| file so that it is available for 
     every |PSIfour| instance. ::

        sys.path.insert(0, '/home/user/psiadditions')

   * Place the path in the input file, either absolute or relative. ::

        sys.path.insert(0, '../../psiadditions')
        sys.path.insert(0, '/home/user/psiadditions')

.. envvar:: PSIDATADIR

   Path in which the |PSIfour| executable looks for its non-compiled
   dependencies (*i.e.*, Python driver, basis sets, databases, *etc.*).
   Not used when running from an installed (``make install``) executable
   or when running from a conda binary,
   so this variable is relevant primarily to developers running the
   executable directly from the compilation directory. Value should be set
   to directory containing driver, basis, *etc.* directories, generally
   ``psi4/share``.

