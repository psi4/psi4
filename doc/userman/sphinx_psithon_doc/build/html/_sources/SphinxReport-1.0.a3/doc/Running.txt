.. _Running:

=================
 Running a build
=================

Building a restructured text document is an iterative process that involves
usually the following steps:

0. Setting up a new :mod:`SphinxReport` project using :ref:`sphinxreport-quickstart`.
1. Writing some text.
2. Adding the python code for :term:`Tracker` as a data source.
3. Testing the :term:`Tracker` and customizing the plot using :ref:`sphinxreport-test`
4. Inserting the :report: directive into the text.
5. Optionally, remove existing cached data of this :term:`Tracker` using :ref:`sphinxreport-clean`.
6. Building the text using :ref:`sphinxreport-build` or the :ref:`Makefile`.

Steps 1-6 are repeated until the document is finished.

Command line utilities
======================

This page explains the various utilities that come with :mod:`SphinxReport`.

.. _Makefile:

Makefile
--------

If a project has been set up with the :file:`sphinxreport-quickstart` utility,
then the build process can be controlled from a Makefile. In the 
:term:`source directory` type::

   make help

for a list of all commands available.

.. _sphinxreport-quickstart:

sphinxreport-quickstart
-----------------------

The :file:`sphinxreport-quickstart` utility sets up a :mod:`sphinxreport`
project. It is called as::

   $ sphinxreport-build -d new_project

**-d/--destination** directory
   setup project using *directory* as the :term:`source directory`.. The default 
   is ".", the current directory.

sphinx-build
------------

At its simplest, :mod:`SphinxReport` is a :mod:`Sphinx` extension
and all images are simply built using the usual :mod:`Sphinx` build.
See the `Sphinx documentation <http://sphinx.pocoo.org/intro.html#running-a-build>`
on how to run a build.

.. _sphinxreport-build:

sphinxreport-build
------------------

Rendering many images and extracting data takes time. The :file:`sphinxreport-build`
utility can speed up this process by running several rendering processes in parallel.
Note that :ref:`Caching` needs to be enabled for this to work. It also takes care of 
building the :ref:`Gallery`. It is invoked as a prefix to the :file:`sphinx-build`
command, for example::
   
   $ sphinxreport-build --num-jobs=4 sphinx-build -b html -d _build/doctrees   . _build/html

Options for sphinxreport-build are:

**-a/--num-jobs** *jobs*
   Number of parallel jobs to execute. The default value is 2.

.. _sphinxreport-clean:

sphinxreport-clean
------------------

The :file:`sphinxreport-clean` utility removes files from a previous built. It is called as::

   $ sphinxreport-clean [*target* [[*tracker] ...]

Where *target* can be one of 

**clean**
   Remove the latest rendered documents, but leaves cached data.

**distclean**
   Remove all build information including cached data.

**<tracker>**
   The name of a tracker. All images, cached data and text elements based
   on this tracker are removed so that they will be re-build during the 
   next build. Multiple trackers can be named on the command line.

.. _sphinxreport-test:

sphinxreport-test
-----------------

The :file:`sphinxreport-test` utility presents previews of graphs and tables. It
can also generate templater restructured text for cutting and pasting into a 
document. It is called as

   $ sphinxreport-test tracker renderer

**-t/--tracker** tracker
   Tracker to use.

**-a/--tracks** tracks
   Tracks to display as a comma-separated list.

**-s/--slices** slices
   Slices to display as a comma-separated list.

**-r/--renderer** renderer
   Renderer to use.

**-o/--option** option
   Options for the renderer. Supply as key=value pairs (without spaces). Several **-o** options can
   be supplied on the command line.

**--no-print**
   Do not print an rst text template corresponding to the displayed plots.

**--no-plot**
   Do not plot.

If no command line arguments are given all Trackers are build in parallel.

.. _sphinxreport-gallery:

sphinxreport-gallery
--------------------

The :file:`sphinxreport-gallery` utility examines the build directory for images
and constructs a gallery. It should be called from the :term:`source directory`.

   $ sphinxreport-gallery

Calling :file:`sphinxreport-gallery` is usually not necessary if :file:`sphinxreport-build`
is used.

.. _Debugging:

Debugging
=========

Information and debugging messages from to the ``report`` directive are
written to the file :file:`sphinxreport.log` in the current directory.

.. _Caching:

Caching
=======

Extracting data from a database potentially takes much time if a lot of processing
is involved or the data set is large. To speed up the writing process :mod:`SphinxReport`
is able to cache function calls to a :term:`Tracker` if the configuration variable
``cachedir`` is set, for example to::

   cachedir=os.path.abspath("_cache")

Enabling caching will speed up the build process considerably, in particular as
:ref:`sphinxreport-build` can make use of parallel data gathering and plotting.
Unfortunately currently there is no :ref:`Dependency` checking for cached data.
Thus, changes in the code of a :term:`Tracker` or changes in the data will not
result in an automatic update of the cache. The best solution is to manually 
delete the cached data using the command :ref:`sphinxreport-clean`.

.. _Dependency:

Dependency checking
===================

:mod:`Sphinx` implements dependency checking such that existing documents are only rebuilt
if the underlying sources have changed. The same dependency checking is still available in 
:mod:`SphinxReport`, however currently there is no dependency checking between the data
source and an existing image. As long as an image or table is present on the file system, it
will not be re-rendered even if the document or the underlying data has changed. To force
re-rendering, use the command :ref:`sphinxreport-clean`.

.. _Gallery:

Gallery
=======

:mod:`SphinxReport` builds a gallery of all plots created similar to the 
`matplotlib gallery <matplotlib.sourceforge.net/gallery.html>`_. The gallery
can be built manually with :file:`sphinxreport-gallery`, but is also built
automatically by :file:`sphinxreport-build`.
