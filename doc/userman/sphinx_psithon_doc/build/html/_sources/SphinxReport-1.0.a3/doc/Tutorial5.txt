.. _Tutorial5:

=============================
Tutorial 5: Customizing plots
=============================

Plots are customized by options supplied to
the ``:report:`` directive. Different renderers
permit different options. For a full list, see
the :ref:`Reference`.

Building up on the example of :ref:`Tutorial4` 
we will now 

   * make a single plot using the ``groupby`` option,
   * plot only lines using the ``as-lines`` option, and
   * add an x-axis label using the ``xtitle`` option.

First, check the plot's look on the command line::

   sphinxreport-test -t ExpressionLevelWithSlices -r histogram-plot -o range=0,100,4 -o groupby=all -o as-lines -o xtitle="expression level"

If you are satisfied, add the following text to your document::

   .. report:: Tutorial4.ExpressionLevelWithSlices
      :render: histogram-plot
      :range: 0,100,4
      :xtitle: expression level
      :groupby: all
      :as-lines:

      Expression level by experiment and function.

Note that the same data can often be presented in different ways. For example,
you might be interested in the stats of these functions::

   .. report:: Tutorial4.ExpressionLevelWithSlices
      :render: stats

      Expression level by experiment and function.

display a box plot::

   .. report:: Tutorial4.ExpressionLevelWithSlices
      :render: box-plot
      :ytitle: expression level
      :groupby: all

      Expression level by experiment and function.


or might want to present the histogram literally::

   .. report:: Tutorial4.ExpressionLevelWithSlices
      :render: histogram
      :range: 0,100,4

See :ref:`Tutorial5Demo` to check how the result should look like.
