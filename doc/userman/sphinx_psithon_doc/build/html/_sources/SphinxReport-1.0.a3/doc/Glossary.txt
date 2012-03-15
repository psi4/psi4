*********
Glossary
*********

.. glossary::
   :sorted:

   render
      The restructured text directive supplied by the SphinxSqlPlot extension.

   track
      A data set, for example species like "frog", "mouse", and "dog".

   slice
      A slice through a data set, for example gender like "male" and "female". 

   Tracker
      A python function or functor returning data for a track, see :class:`Tracker.Tracker`.

   SingleColumnData
      Return type of a :class:`Tracker.Tracker`. SingleColumnData is a single list or tuple of data,
      for example ``(1,2,3,4)``.
      
   MultipleColumnData
      Return type of a :class:`Tracker.Tracker`. MultipleColumnData is a list/tuple of lists/tuples.
      The first tuple/list contains column labels, while each subsequent tuple/list contains the values
      of a column. For example ``( ("column1","column2), ((1,2,3,4), (5,6,7,8)) )`` corresponds to the
      following data::
      
         column1 column2
         1	 5
         2	 6
	 3	 7
	 4	 8

       In contrast to :term:`LabeledData`, all the columns are required to have the same length.

   LabeledData
      Return type of a :class:`Tracker.Tracker`. MultipleColumnData is a list/tuple of lists/tuples.
      The first tuple/list contains column labels, while each subsequent tuple/list contains the values
      of a column. For example ``( ("data1", 1), ("data2", 2))`` corresponds to the
      following data::
         
         data1   data2 
         1	 2

   source directory
      The directory which, including its subdirectories, contains all source
      files for one Sphinx project.

   configuration directory
      The directory containing :file:`conf.py`.  By default, this is the same as
      the :term:`source directory`, but can be set differently with the **-c**
      command-line option.
