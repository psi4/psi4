.. _Tutorial2:

==========================
 Tutorial 2: Using tracks
==========================

The previous tutorial created a simple bar-plot of a single data source. More useful plots
show several data sources, called ``tracks`` in a single plot. This tutorial will show how 
to do this.

**********************
Converting to functors
**********************

Instead of using a python function as in the first tutorial, we will define a functor.
The functor provides additional methods that allows the renderer to query for available
``tracks``.

Adding a data source
********************

Create the file :file:`Tutorial2.py` in the :file:`python` subdirectory and add 
the following code::

   from SphinxReport.DataTypes import *
   from SphinxReport.Tracker import *

   class MyDataOneTrack(Tracker):
      """My one-tracked data."""

      def getTracks( self ):
          return ["all",]

      @returnLabeledData
      def __call__(self, track, slice = None ):
          return [ ("header1", 10), ("header2", 20) ]

Apart from :mod:`DataTypes` the module :mod:`Trackers` is imported
and the data source ``MyData`` is derived from it::
   
   class MyData(Tracker):

The docstring::

      """My one tracked data."""

will later reappear in the caption of the table or figure. The method :meth:`getTracks` returns
the tracks provided by this tracker::

      def getTracks( self ):
      	  return ["all",]

In this example, there is only one track called ``all``.

Finally, the method :meth:`__call__` provides the data::

      @returnLabeledData
      def __call__(self, track, slice = None ):
          return [ ("header1", 10), ("header2", 20) ]

Testing the data source
***********************

Testing the current implementation::

   sphinxreport-test -t MyDataOneTrack -r bars

will show a familiar plot - the functor is equivalent
to the single funtion case.

******************
Adding more tracks
******************

The functor now permits adding more flexibility to the data
source.

Adding a data source
********************

Add the following code to :file:`Tutorial1.py`::

   class MyDataTwoTracks(Tracker):
      """My one-tracked data."""

      def getTracks( self ):
          return ["track1","track2"]

      @returnLabeledData
      def __call__(self, track, slice = None ):
      	  if track == "track1":
	     return [ ("header1", 10), ("header2", 20) ]
      	  elif track == "track2":
	     return [ ("header1", 20), ("header2", 10) ]

Testing the data source
***********************

Testing the current implementation::

   sphinxreport-test -t MyDataTwoTracks -r bars

will now show two bars side-by-side. Try out::

   sphinxreport-test -t MyDataTwoTracks -r stacked-bars

Creating a restructured text document
*************************************

Create the following :file:`Tutorial2.rst` (and add it to :file:`index.rst`)::

    ==========
    Tutorial 2
    ==========

    My new bar plots:

    .. report:: Tutorial2.MyDataOneTrack
       :render: bars

       My first bar plot - this time as a functor

    .. report:: Tutorial2.MyDataTwoTracks
       :render: bars

       My new bar plot - two tracks

    .. report:: Tutorial2.MyDataTwoTracks
       :render: stacked-bars

       My new bar plot - same data, different renderer

Note that the same data can appear several times in the same document
with different renderers. See :ref:`Tutorial2Demo` to check 
how the result should look like.
