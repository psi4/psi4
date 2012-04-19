
.. include:: autodoc_abbr_options_c.rst

.. _`sec:mp2`:

Second-order M\ |o_slash|\ ller--Plesset Theory: MP2 and MP2-R12 methods
========================================================================

Second-order M\ |o_slash|\ ller--Plesset theory is one of the most basic
wavefunction approaches which includes electron correlation
directly.
Due to its simplicity, the MP2 method is often the best
level one can afford for a larger molecular system.
At the other end of the spectrum, the MP2-R12 method
of Kutzelnigg, Klopper, and co-workers is a promising
approach to computing MP2 energies in the complete
basis set limit for smaller systems. |PSIfour| is
one of the very few publicly available programs to
feature a robust implementation of the MP2-R12 method.

|PSIfour| is capable of computing closed-shell
MP2 and MP2-R12/A energies using integral-direct techniques and a
multithreaded algorithm, which lends itself perfectly for execution 
on symmetric multiprocessor (SMP) machines. |PSIfour| is also
capable of computing RHF, UHF, and ROHF (using semicanonical orbitals)
MP2 energies and one-particle density matrices, and RHF MP2 analytic 
gradients.  Occupied and virtual orbitals can be frozen during the 
energy calculation, but not for the calculation of the 
one-particle density matrix or the analytic gradient.

.. table:: Summary of MP2 and MP2-R12 capabilities in |PSIfour|

   +-----------+-----------+---------------+--------------------------+----------+
   | Reference | Method    | Energy (conv) | Energy (integral-direct) | Gradient |
   +===========+===========+===============+==========================+==========+
   | RHF       | MP2       | Y             | Y                        | Y        |
   +-----------+-----------+---------------+--------------------------+----------+
   | UHF       | MP2       | Y             | ---                      | ---      |
   +-----------+-----------+---------------+--------------------------+----------+
   | ROHF      | MP2       | Y             | ---                      | ---      |
   +-----------+-----------+---------------+--------------------------+----------+
   | RHF       | MP2-R12/A | ---           | Y                        | ---      |
   +-----------+-----------+---------------+--------------------------+----------+

Basic Keywords
--------------

To compute a ground-state MP2 or MP2-R12 energy at a fixed geometry,
the following keywords are common:

\item[WFN = string]\mbox{}\\
Acceptable values are {\tt mp2} for MP2, {\tt mp2r12} [for MP2-R12/A]
There is no default.  

\item[REFERENCE = string]\mbox{}\\
The only acceptable value are {\tt rhf, uhf, and rohf}.
There is no default.

\item[JOBTYPE = string]\mbox{}\\
Acceptable values are {\tt sp} and {\tt opt}.  There is no default.

\item[MEMORY = (real MB)]\mbox{}\\
Specified the amount of core memory to be used, in MB.  Defaults to 256.
Other units (e.g., KB or GB) are also allowed.

.. comment include:: autodir_options_c/detci__reference.rst

\item[DIRECT = boolean]\mbox{}\\
Specifies whether to use the conventional ({\tt false}) or
integral-direct ({\tt true}) algorithm. Default is {\tt false}.

\item[NUM\_THREADS = integer]\mbox{}\\
Specified the number of threads to be used in the integral-direct
computation (only valid if {\tt DIRECT} is set to {\tt true}).
Default is 1.

\item[FREEZE\_CORE = boolean]\mbox{}\\
Specifies whether core orbitals (which are determined automatically) are to
be excluded from the correlated calculations.  Default is {\tt false}.

\item[PRINT = integer]\mbox{}\\
The desired print level for detailed output.  Setting this to 2 is a good
idea for larger calculations so that the progress of the calculation may be
easily followed.  Defaults to 0.

\item[OPDM = boolean]\mbox{}\\
If {\tt true}, calculate the one-particle density matrix.  The default is false.

\item[OPDM\_WRITE = boolean]\mbox{}\\
If {\tt true}, write the one-particle density matrix to disk.

\item[OPDM\_PRINT = boolean]\mbox{}\\
If {\tt true}, print the one-particle density matrix to the output file.

Using the MP2-R12 method
------------------------

Although this manual is not a how-to on running
quantum chemistry applications, the MP2-R12 method is
a rather non-standard tool, hence a few comments on its
use are appropriate.

* The version of the MP2-R12 method implemented in |PSIfour|
  is a so-called single-basis MP2-R12 method
  in standard approximation A. This means that a basis set
  rather complete in Hartree--Fock (or one-particle) sense
  is absolutely mandatory for meaningful computations with the MP2-R12
  method. The user is strongly urged to read literature on
  linear R12 methods before using |PSIfour| to compute MP2-R12
  energies.

* More robust, two-basis versions
  of the MP2-R12 method, also known as the auxiliary basis
  MP2-R12 method, have been implemented
  in a publicly available Massively Parallel Quantum Chemistry (MPQC)
  package (see \url{http://aros.ca.sandia.gov/~cljanss/mpqc/}).
  The two-basis version of the MP2-R12 method is a theoretically more
  sound approach, and thus should be preferred to the single-basis method.
  In some situations, however, it may make sense to use
  the single-basis method.

Larger Calculations
-------------------

Here are a few recommendations for carrying out extended integral-direct MP2 and
MP2-R12 calculations with |PSIfour|: 

* While the integral-direct MP2 algorithm doesn't need any
  significant disk storage,
  the integral-direct algorithm for the MP2-R12 energy
  stores the transformed integrals to disk, hence very large
  computations will require a lot of disk space. In general
  the storage requirement is :math:`16 o^2N^2` bytes, where :math:`o`
  is the number of occupied orbitals, and :math:`N` is the size of the basis.

* If there is not enough memory to perform the computation in one pass,
  the program will do multiple passes through the entire set of integrals,
  hence your computation will run that many times longer.
  In such case, find the machine with the most memory and processors available.

* On SMP machines, set the {\tt NUM\_THREADS} to the number of
  processors available for the job, or, if all processors are allocated for
  your job, set {\tt NUM\_THREADS} to {\em twice} the number of processors
  you have. Modern operating systems schedulers are usually very efficient
  at handling multithreaded programs, so the overhead of thread context
  switching is not significant, but using more threads may lead to better
  load balancing, and lower execution times. For example, on a 32-processor
  IBM eServer p690 we found that the optimal number of threads was 128.
  For the optimal performance, do a few runs with different number of threads
  and see which number works best.
  Avoid excessively large
  number of threads, as this descreases the net amount of memory available to
  the computation and thus may increase the number of passes. 

* Set the {\tt MEMORY} keyword to the 90% of the available physical
  memory, at most. There is a small amount of overhead associated with the
  integral-direct algorithms that is not accounted for by the internal memory
  handling routines.

* The implementation of the integral-direct MP2-R12 (and MP2) method
  in |PSIfour| can run efficiently on SMP, or shared-memory, machines,
  by utilizing multiple processors via multithreaded approach.
  However, it cannot utilize distributed memory machines,
  such as commodity (PC) clusters and massively parallel machines,
  to their full potential, since one computation can only take advantage
  of one node of such machine at a time.
  In such environments, the aformentioned MPQC implementation of
  the MP2-R12 method should be preferred
  (see \url{http://aros.ca.sandia.gov/~cljanss/mpqc/}).

