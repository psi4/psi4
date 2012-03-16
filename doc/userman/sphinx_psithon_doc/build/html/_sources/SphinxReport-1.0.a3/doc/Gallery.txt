*******
Gallery
*******

The gallery lists examples of all plots and available renderers.

Labeled data
************

These Renderers require :class:`DataTypes.LabeledData`. Labeled data is
consists of a list of labels and a list of data, for
example ``[ ("column1", "column2" ), ( 10, 20) ]``.

Table
=====

The :class:`Renderer.RendererTable` class inserts labeled data into
a table.

.. report:: Trackers.LabeledDataExample
   :render: table

   Table caption.

Stacked Barplot
===============

The :class:`Renderer.RendererStackedBars` class presents labeled data
as stacked bars.

.. report:: Trackers.LabeledDataExample
   :render: stacked-bars

   Figure caption.

Interleaved Barplot
===================

The :class:`Renderer.RendererInterleavedBars` class presents labeled data
as interleaved bars. Both *interleaved-bars* and *bars* can be used.

.. report:: Trackers.LabeledDataExample
   :render: interleaved-bars

   Figure caption.

Single column data
******************

Renderers that accept a single column of data of type :class:`DataTypes.SingleColumnData` 
per :term:`track` and :term:`slice`, for example ``[2,3,1,3,4]``.

Histogram
=========

The :class:`Renderer.RendererHistogram` class computes a histogram
of data and inserts it as a table.

.. report:: Trackers.SingleColumnDataExample
   :render: histogram
   :bins: arange(0,10)

   Table caption.

HistogramPlot
=============

The :class:`Renderer.RendererHistogramPlot` class computes a histogram
of data and inserts a plot.

.. report:: Trackers.SingleColumnDataExample
   :render: histogram-plot
   :bins: arange(0,10)

   Figure caption.

Stats
=====

The :class:`Renderer.RendererStats` class computes summary
statistics and displays them in a table.

.. report:: Trackers.SingleColumnDataExample
   :render: stats

   Table caption.

Boxplot
=======

The :class:`Renderer.RendererBoxplot` class computes boxplots.

.. report:: Trackers.SingleColumnDataExample
   :render: box-plot

   Figure caption.

Multi-column data
*****************

Renderers that accept multiple columns of data of type :class:`DataTypes.SingleColumnData` 
per :term:`track` and :term:`slice`. The data is asscociated with column headers, for example 
``[ ('column1', 'column2'), ( ( 1,2,3), (4,5,6) )]``.

Pairwise statistics table
=========================

Compute correlation statistics between all columns.

.. report:: Trackers.MultipleColumnDataExample
   :render: pairwise-stats

   Figure caption.

Pairwise statistics plot
========================

Plot correlation coefficients between all columns.

.. report:: Trackers.MultipleColumnDataExample
   :render: pairwise-stats-plot

   Figure caption.

Pairwise scatter plot
========================

A scatter plot.

.. report:: Trackers.MultipleColumnDataExample
   :render: scatter-plot

   Figure caption.


Grouped table
=============

A grouped table.

.. report:: Trackers.MultipleColumnsExample
   :render: grouped-table

   Figure caption.
