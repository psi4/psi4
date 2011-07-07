///  \file mainpage.h
///  \brief MADNESS 

///  \defgroup configuration MADNESS installation and configuration
///  \defgroup libraries MADNESS libraries

///          \defgroup parallel_runtime Parallel programming environment
///          \ingroup libraries

///                  \defgroup world Distributed computing environment (World and its relations)
///                  \ingroup parallel_runtime

///                          \defgroup rmi Remote method invocation
///                          \ingroup parallel_runtime

///                          \defgroup mpi Interfaces from World to MPI
///                          \ingroup parallel_runtime

///                          \defgroup worldobj Globally addressable objects (WorldObject)
///                          \ingroup parallel_runtime

///                          \defgroup worlddc Distributed containers (WorldContainer)
///                          \ingroup parallel_runtime

///                          \defgroup futures Futures
///                          \ingroup parallel_runtime

///                          \defgroup serialization Serialization
///                          \ingroup parallel_runtime

///                          \defgroup hashing Hashing
///                          \ingroup parallel_runtime

///                          \defgroup threading Multi-threading 
///                          \ingroup parallel_runtime

///                                  \defgroup mutexes Mutexes
///                                  \ingroup threading

///                                  \defgroup atomics Atomic operations
///                                  \ingroup threading

///                                  \defgroup threads Threads
///                                  \ingroup threading

///                                  \defgroup threadpool Thread pool
///                                  \ingroup threading

///                                  \defgroup taskq Task queue
///                                  \ingroup threading

///                                  \defgroup conhash Concurrent hash table
///                                  \ingroup threading

///         \defgroup mra Multiresolution analaysis
///         \ingroup libraries

///                 \defgroup funcplot Function plotting routines
///                 \ingroup mra

///                 \defgroup mrabcext Exterior boundary conditions
///                 \ingroup mra

///                 \defgroup mrabcint Preliminary support for interior boundary conditions
///                 \ingroup mra

///         \defgroup tensor Tensors or multidimension arrays
///         \ingroup libraries

///         \defgroup linalg Linear algebra (interface to LAPACK)
///         \ingroup libraries

///         \defgroup solvers Iterative solvers for linear/non-linear equations and optimizers
///         \ingroup libraries

///         \defgroup misc Miscellany
///         \ingroup libraries

/// \defgroup applications MADNESS applications

///          \defgroup examples Examples 
///          \ingroup applications

/*!

\mainpage

These pages serve as the main programmer's reference manual for MADNESS and 
were automatically generated from the source using Doxygen.

A good place to start is the <a href="modules.html">modules page</a>
that provides access to documentation for libraries, examples, and 
applications.  The examples are intended to meet the needs of those 
starting to develop new applications.

Additional technical documentation can be found in the MADNESS \c doc directory
  - The getting started guide introduces use of the numerical functionality
  - The implementation notes document the internal numerical workings 
    as a resource for both users and MADNESS developers.
  - The parallel runtime and API documents document the operation and 
    interface for the parallel programming environment.

The <a href="http://code.google.com/p/m-a-d-n-e-s-s/">project home</a>
contains information on configuring, building, etc.

\section status Status and supported platforms

MADNESS is roughly in beta release.  Most things work well enough for
several applications to be in production mode, but on the other hand there
are still lots of "rough" spots and the documentation in particular does not
yet meet our expectations.

Currently, we intend that MADNESS work correctly and efficiently
on the following platforms
  - Cray XT, 
  - IBM BG/P,
  - Linux workstations with x86 32-bit and 64-bit multi-core processors, and 
    clusters thereof
  - Apple Macintosh (Intel processors) with recent versions of OS X.

To build on any of these you will need a recent software stack.

\subsection used_by Software used by MADNESS

Our deep gratitude to these other projects whose software we employ
  - <a href="http://sourceforge.net/projects/tinyxml">TinyXML</a>
  - <a href="http://www.librow.com/articles/article-10/appendix-a-2">CFFT</a>
  - <a href="http://muparser.sourceforge.net">muParser</a>
  - <a href="http://www.tddft.org/programs/octopus/wiki/index.php/Libxc">libxc</a>
  - <a href="http://www.doxygen.org">Doxygen</a>
  - <a href="http://www.graphviz.org">Graphviz</a>
  - <a href="http://www.mcs.anl.gov/research/projects/mpich2">MPICH</a>
  - <a href="http://www.open-mpi.org">Open-MPI</a>
  - <a href="http://www.netlib.org/lapack">LAPACK</a>
  - <a href="http://gcc.gnu.org">GCC</a>
  - <a href="http://code.google.com/p/googletest">Google test</a>
  - <a href="http://code.google.com/p/google-perftools">Google performance tools</a>

\section license License

\verbatim
  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
\endverbatim


*/
