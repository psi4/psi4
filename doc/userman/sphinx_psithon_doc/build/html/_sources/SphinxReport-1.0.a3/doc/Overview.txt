.. _Overwiew:

********
Overview
********

The :mod:`SphinxReport` module is an extension for :mod:`sphinx`
that provides facilities for data retrieval and data rendering
within reStructured text. 

.. _Features:

Features
********

:mod:`SphinxReport` is a report generator that is implemented as an extension
to :mod:`Sphinx`. It is easy to use and powerful enough to give all the flexibility 
needed during the development of computational pipelines and robustness during the
production use of a pipeline.

It is intended for developers and computational scientists with ``python`` scripting experience.

Briefly, the :mod:`SphinxReport` is a report generator that is implemented as an extension
to :mod:`Sphinx`. It

* uses simple markup in the form of restructured text
* supports both automated and narrative analysis
* keeps code and annotation together
* takes care of indexing, formatting to produce something pretty
* produces static html and pdf
* provides the power of python within text

Usage
*****

The :mod:`SphinxReport` extension has three parts:

1. ``Renderers`` are provided by the :mod:`SphinxReport` extension and are used to display data
   as graphs or tables. They can be incorporated into restructured text files using the ``:report::`` directive.

2. ``Trackers`` are python classes and are provided by the user and collect data. The data can be obtained from
   data sources like flat files, SQL tables and more. Trackers provide data by :term:`track` and :term:`slice`. 
   Tracks are principal collections of data, for example data measurements taken from different species. Slices 
   are subsections of the data. For example, a :class:`Tracker` might provide *weight* measurements of different species
   according to the slice *gender*.

3. :file:``sphinxreport-build`` is used to build a document. 

The following minimal example illustrates how Renderers and Trackers work together. A ``:report:``
directive like::

   .. report:: BarData
      :render: bars

will insert a barplot (:class:`RendererInterleavedBars`) at the current
location. The Renderer will obtain the data from a python class or function *BarData* in the file
:file:`Trackers.py` that should be somewhere in the python search path (see :ref:`Configuration`).
The function *BarData* might look like this::

   @returnLabeledData
   def BarData(): return [("bar1", 10), ("bar2", 20)]

The document is built using the usual :mod:`Sphinx` process::

   sphinx-build -b html -d _build/doctrees   . _build/html

See the :ref:`Tutorials` for a more complete introduction on how to use the extension. See
:ref:`Running` on more advanced building methods.



.. _Background:

Background
**********

Scientific datasets these days are large and are usually processed by
computational pipelines creating a wealth of derived data, very often 
stored in a database. With computational power always increasing, 
the bottleneck is usually the subsequent analysis. 

Especially during code development and in the early exploratory stages, the data 
are sliced and plotted in multiple ways to find problems and understand the data. 
At the same time, the plots and tables are embedded into text with comments and 
notes that should later result in a publication. As bugs are fixed and the data 
are understood better, the plots and tables need to be frequently updated. Statically
copying and pasting images into a document becomes tedious quickly.

The interactive analysis is later followed by re-runs of the pipeline
on different data sets or with different parameters. Again the data is sliced
and plotted, this time to confirm the successful completion of the pipeline
and to compare results to those of previous runs. This is a mostly automatic
task, in which diagnostic plots are created to provide a high-level view
of the results. There is also an interactive component, where plots are 
selected to highlight unexpected deviations that are the bread-and-butter of science.

We found no tool that easily bridges the divide of interactive analysis and
automatic updating. On one end of the spectrum is office software with macros
or embedded images linked to physical files. Writing in office software is easy, 
there is drag & drop and the result is very close to the desired product: a
publishable manuscript. However, with complicated analyses the macros become 
unwieldy. Images on the hard-disc separate the code to create the images from 
the document and there is always the danger of links being broken. Taking a live
document and applying it to a new dataset is difficult.

At the other end of the spectrum are full-fledged content management systems
that provide dynamic access to the data. These have a steep learning curve and
require a lot of work to build and maintain. Some design is necessary beforehand
to prevent uncontrolled growth. Unfortunately this is usually at odds with
our experience how computational pipelines in science develop. Such effort is 
usually only justifyable for large pipelines, big projects and big teams.

Somewhere in the middle of the spectrum are report generators. These create 
static documents, but are designed to be run often and on different datasets. 
These are powerful, but often have a steep learning curve. We also found them
lacking in plotting capabilities. 

We thought the combination of :mod:``Sphinx`` and :mod:``matplotlib``
and ideal combination and extended the ``matplotlib`` ``:plot:`` directive
to interactively collect data. We are heavily indebted to these two
projects.

.. seealso::

   Sphinx: 
      http://sphinx.pocoo.org

   Matplotlib:
      http://matplotlib.sourceforge.net

   Python:
      http://www.python.org

   A restructured text quick reference: 
      http://docutils.sourceforge.net/docs/user/rst/quickref.html






