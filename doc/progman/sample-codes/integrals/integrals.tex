As noted previously, we want to start from a code that's not too tightly
integrated with the \PSIfour\ code itself, so we begin with a \file{Makefile}
that will allow us to write a standalone code that includes all requisite \PSI\
libraries.  We're going to write a small sample code that generates integrals,
which involves two just two source files.  We begin by defining a
\file{Makefile} that will include all of the \PSIfour\ libraries and header
files, so that we can take full advantage of the wide range of features
implemented without having to worry about the details of their implementation.

\includesource{sample-codes/integrals/Makefile}{make}

Only a few lines of this makefile need to be modified to utilize it for other
programming projects; we'll concentrate on them.  On the second line, we define
the name of the executable to be generated, in this example we opt for the
unimaginative title of \module{integrals}.  Line 4 provides the list of source files
that the project comprises; these will be detailed below.  The top source
directory for the \PSIfour\ installation and the top object directory (where
\PSIfour\ was compiled) should be provided on lines 6 and 8, respectively.
Lines 10 and 11 describe the flags needed to link in the \module{BLAS} and
\module{LAPACK} libraries and might need a combination of ``\file{-L
folder\_name}'' and ``\file{-l library\_name}'', depending on your system's
setup.  Finally, the compiler and flags are detailed on lines 12--17.  It's a
good idea to use the flags described on line 16 for development; they speed up
code compilation and provide lots of information for standard debugging tools.
As noted in the \file{Makefile} itself, nothing below line 17 should require
modification for any other \PSIfour\ project.

The \PSIfour\ driver program provides a lot of functionality that we forgo in
writing a standalone code; this is instead emulated in the {\tt main.cc} file,
shown below.

\includesource{sample-codes/integrals/main.cc}{C++}

All modules in \PSIfour\ must have the argument list and return type shown on
line 13.  The possible return types, defined by an enumeratable constant are
documented in \file{psi4-dec.h}, which lives in \$PSI4/include.  Notice that
all of the code must live in it's own namespace within the \module{psi}
namespace, in this case it's in the \module{psi::integrals} namespace.  Without
this nesting, functions belonging to different parts of the code, but having
the same name, would cause conflicts.  The \module{read\_options} function is
responsible for setting up the \module{Options} object, which contains the list
of user-provided options.  Lines 25--32 are important - these provide the list
of keywords expected by the code, their types, and their default values (if
any).  This part of the code will be inserted into the \PSIfour\ driver when
the module is ready for merging with the \PSIfour\ distribution; this process
will be detailed later in the chapter.  Notice the special format of the
comments on lines 27 and 30.  These are still valid \module{C++} comments, but
the extra hyphens inside are essential in this context.  Whenever adding any
options for any module, you must comment them as shown - this will ensure that
the keywords are automatically inserted into the \PSIfour\ users' manual.  The
\module{main} function does a little setting up of the \PSI\ input and output
environments, before calling the module code we're developing (on line 53) and
shutting down the \PSIfour\ I/O systems.

The module we're developing is in the following source file.

\includesource{sample-codes/integrals/integrals.cc}{C++}

Given the extensive documentation within the code, we'll not describe this file
line-by-line; however, some points warrant elaboration.  Notice that the entire
module is encapsulated in the \module{psi::integrals} namespace (lines 6 and
92).  This simple example has only one function body, which lives in a single
source file - if more functions and/or source files were added, these too would
have to live in the \module{psi::integrals} namespace.  On lines 29 and 31 of
\file{main.cc} we told the parser which keywords to expect, and provided
default values in case the user omited them from the input.  This makes
retrieving these options very clean and simple ({\it c.f.} lines 11 and 12 of
\file{integrals.cc}).  Each \PSIfour\ module will have to initialize its own
local \module{PSIO} and \module{Chkpt} objects to perform I/O and to retrieve
information from previously run modules.  Notice that these objects are created
within smart pointers (see section XXX for more information) so that they are
automatically deleted when they go out of scope, thus reducing the burden on
the programmer.  Likewise, the basis sets, matrices and integral objects are
allocated using smart pointers.

The code described above can be built by simply typing ``make'' on the command
line.  To run this code, you must first run the \module{input} module to read
in the basis set information.  A \PSI\ input for this code should look some
thing like the following:

\includeinput{sample-codes/integrals/input.dat}
