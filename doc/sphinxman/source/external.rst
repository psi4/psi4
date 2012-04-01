
.. _`sec:config`:

==================================================
Configuration: Preparing |PSIfour|\ 's Environment
==================================================

.. index:: scratch files
.. index:: psi4rc
.. _`sec:psirc`:

Scratch Files and the |psirc| File
==================================

One very important part of user configuration at the end of the
installation process is to tell |PSIfour| where to write its temporary
("scratch") files.  Electronic structure packages like |PSIfour| can
create rather large temporary disk files.  It is very important to 
ensure that |PSIfour| is writing its temporary files to a disk drive
phsyically attached to the computer running the computation.  If it
is not, it will significantly slow down the program and the network.
By default, PSI4 will write temporary files to ``/tmp``, but this
directory is often not large enough for typical computations.  Therefore,
you need to (a) make sure there is a sufficiently large directory on a
locally attached disk drive (100GB–1TB or more, depending on the size of
the molecules to be studied) and (b) tell |PSIfour| the path to this
directory.  The |PSIfour| installation instructions explain how to set up a
resource file, |psirc| (example :source:`samples/example_psi4rc_file`),
for each user providing this information.

For convenience, the Python interpreter will execute the contents of the
|psirc| file in the current user's home area (if present) before performing any
tasks in the input file. The primary use of the |psirc| file is to control the
handling of scratch files.  |PSIfour| has a number of utilities that manage
input and output (I/O) of quantities to and from the hard disk.  Most
quantities, such as molecular integrals, are intermediates that are not of
interest to the user and can be deleted after the computation finishes, but
pertinent details of computations are also written to a checkpoint file and
might be useful in subsequent computations.  All files are sequentially
numbered and are written to ``/tmp``, then deleted at the end of the computation,
unless otherwise instructed by the user.

A Python callable handle to the |PSIfour| I/O management routines is available,
and is called ``psi4_io``.  To instruct the I/O manager to send all files to
another location, say ``/scratch/user``, add the following command to the |psirc|
file (note the trailing "/")::

    psi4_io.set_default_path('/scratch/user/')

For batch jobs running through a queue, it might be more convenient to use an
environmental variable (in this case ``$MYSCRATCH``) to set the scratch directory;
the following code will do that::

    scratch_dir = os.environ.get('MYSCRATCH')
    if scratch_dir:
        psi4_io.set_default_path(scratch_dir + '/')

Individual files can be send to specific locations.  For example, file 32 is
the checkpoint file that the user might want to retain in the working directory
(*i.e.*, where |PSIfour| was launched from) for restart purposes.  This is
accomplished by the commands below::

    psi4_io.set_specific_path(32, './')
    psi4_io.set_specific_retention(32, True)

To circumvent difficulties with running multiple jobs in the same scratch, the
process ID (PID) of the |PSIfour| instance is incorporated into the full file
name; therefore, it is safe to use the same scratch directory for calculations
running simultaneously.

To override any of these defaults for selected jobs, simply place the
appropriate commands from the snippets above in the input file itself.  During
excecution, the |psirc| defaults will be loaded in first, but then the commands
in the input file will be executed.  Executing |PSIfour| with the :option:`psi4 -m` (for
messy) flag will prevent files being deleted at the end of the run::

    psi4 -m

Alternately, the scratch directory can be set through the environment
variable :envvar:`PSI_SCRATCH` (overrides |psirc| settings).

The |psirc| file can also be used to define constants that are accessible
in input files or to place any Python statements that should be executed
with every |PSIfour| instance.

.. _`sec:commandLineOptions`:

Command Line Options
====================

|PSIfour| can be invoked with no command line arguments, as it takes as input
by default the file "input.dat" and directs output by default to "output.dat".
The set of three commands below are completely equivalent, while the fourth is,
perhaps, the most common usage. ::

   psi4
   psi4 -i input.dat -o output.dat
   psi4 input.dat output.dat
   
   psi4 more_descriptive_filename.in more_descriptive_filename.out

Command-line arguments to |PSIfour| can be accessed through ``psi4 --help``.

.. program:: psi4

.. option:: -a, --append

   Append results to output file. Default: Truncate first

.. option:: -h, --help

   Display the command-line options and usage information.

.. option:: -i <filename>, --input <filename>

   Input file name. Default: input.dat

.. option:: -o <filename>, --output <filename>

   Output file name. Use ``stdout`` as <filename> to redirect 
   to the screen. Default: output.dat

.. option:: -m, --messy

   Leave temporary files after the run is completed.

.. option:: -n <threads>, --nthread <threads>

   Number of threads to use (overrides :envvar:`OMP_NUM_THREADS`)

.. option:: --new-plugin <name>

   Creates a new directory <name> with files for writing a
   new plugin. An additional argument specifies a template
   to use, for example: ``--new-plugin name +mointegrals``.
   See Sec. :ref:`sec:plugins` for available templates.

.. option:: -p <prefix>, --prefix <prefix>

   Prefix for psi files. Default: psi

.. option:: -v, --verbose

   Print a lot of information

.. option:: -d, --debug

   Flush the outfile at every fprintf. Default: true iff ``--with-debug``

.. option:: -V, --version

   Print version information.

.. option:: -w, --wipe

   Clean out scratch area.


.. _`sec:environmentVariables`:

Environment Variables
=====================

These environment variables will influence |PSIfour|\ 's behavior.

.. envvar:: OMP_NUM_THREADS

   Number of threads to use by modules with OpenMP threading.

.. envvar:: MKL_NUM_THREADS

   Number of threads to use by operations with Intel threaded BLAS libraries.

.. envvar:: PSI_SCRATCH

   Directory where scratch files are written. Overrides settings in |psirc|.

