.. _Tutorial1:

***********************
Tutorial 1: First steps
***********************

This first tutorial shows step-by-step how a new plot
can be added to a restructured text document.

=============================
Creating the skeleton project
=============================

First, we create a skeleton report in the directory tutorial::

   sphinxreport-quickstart -d tutorial

This will create a skeleton report layout in the directory :file:`tutorial`.
Enter this directory and build the report as html::

   make html

This will create the skeleton report. View it by opening the file 
:file:`tutorial/_build/html/index.html` in a web browser.

==============
Adding a graph
==============

In this tutorial, will insert a simple bar plot into a restructured text document.
There are three steps involved:

1. Define a data source
2. Test the data source
3. Insert ``report`` directive in document

Adding a data source
====================

Create the file :file:`Tutorial1.py` in the :file:`python` subdirectory and add 
the following code::

   from SphinxReport.DataTypes import *

   @returnLabeledData
   def MyDataFunction():
      return [ ("header1", 10), ("header2", 20) ]

In the first line the module :mod:`DataTypes` is imported which
defines decorators used for ensuring type consistency between
data sources and ``Renderers``. This is all that is required in this
simple example. 

In the above example, the function :meth:`MyDataFunction` 
returns ``DataTypes.LabeledData``, which is simply a list of tuples 
of ``(key,value)`` pairs. The type of data returned requires on the renderer. 
For example, a histogram requires a list or tuple of values. 

Testing the data source
=======================

The utility :file:`sphinxreport-test` can be used to check if a
data source works. To test your data source, type::

   sphinxreport-test -t MyDataFunction -r bars

on the command prompt in the root directory. This should create a 
bar plot (``-r bars`` or ``--renderer=bars`` ) of your data source
(``-t MyDataFunction`` or ``--tracker=MyDataFunction``).

:file:`sphinxreport-test` will also produce the restructured text
required to render this graph within a restructured text document.
This utility is very useful for fine-tuning the appearance
of a plot before inserting it into the main document.

Creating a restructured text document
=====================================

Create the following rst file:`Tutorial1.rst`::

    ==========
    Tutorial 1
    ==========

    My first bar blot:

    .. report:: Tutorial1.MyDataFunction
       :render: bars

       My first bar plot.

The the ``report`` directive is used to insert the graph into 
the restructured text document. Behind the scenes, :mod:`sphinx` will call 
the :mod:`SphinxReport` extension and request a barplot. The :mod:`SphinxReport` in 
turn will look for a data source :meth:``MyDataFunction`` in the :file:``Tutorial1.py`` 
that should be somewhere in your :env:``PYTHONPATH`. The default location for
these is in the :file:``python`` subdirectory under the main installation
directory. The content of the ``report`` directive is the figure or table caption.

Add a link to the contents section in the :file:`index.rst` and rebuild the sources::

    make html

There should now be a tutorial1 section in your document 
with a barplot. See :ref:`Tutorial1Demo` how it should look
like.

The next Tutorial (:ref:`Tutorial2`) will cover more complex
data sources and plots.














