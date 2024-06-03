Tips 
#######################################

Silencing status updates
---------------------------------------

To silence status updates to standard out, set ``silent=True``::

   import xagg as xa
   xa.set_defaults(silent=True)

You can also silence a single operation using ``with``::

   # Calculate weightmaps without status updates
   with xa.set_defaults(silent=True):
      weightmap = xa.pixel_overlaps(ds,gdf)

Note that the `silent=True` option in individual functions will be
slowly deprecated over the next few versions in favor of using one of
the two options above.

Saving weights file 
---------------------------------------
If calculating weights from rasters is taking a substantial amount of time (e.g., your raster is very high resolution), you can save the calculated weights using::

   # Create weightmap
   weightmap = xa.pixel_overlaps(ds,gdf,silent=True)

   # Save weightmap
   weightmap.to_file('weights')

   # Read weightmap
   weightmap = xa.read_wm('weights')

   # Continue as usual... 
   aggregated = xa.aggregate(ds,weightmap)

Note that ``weightmap.to_file(fn)`` creates and populates a separate _directory_ named ``fn`` to be able to store all the relevant components of the :py:class:`weightmap` class, including shapefiles with the geometry of the input polygons, the dataframe with the pixel overlap data, the source grid, and any additional weight grids.

This feature is still slightly experimental, so please let us know your experiences! 

Speed up overlap calculation
---------------------------------------
At the expense of increased memory usage, processing may be sped up using an alternate calculation method (``impl='dot_product'``) :: 

   import xagg as xa
   # Set dot_product as default backend implementation 
   xa.set_defaults(impl='dot_product')

Note that the `impl=dot_product` option in individual functions will be
slowly deprecated over the next few versions in favor of using the option
setting approach above (or using ``with xa.set_defaults(impl='dot_product'):`` ).

Create diagnostic figure to inspect raster/polygon overlaps 
------------------------------------------------------------
Once you have created a :py:class:`weightmap`, the following code will create a diagnostic figure, showing a particular polygon (or groups of polygons) + the grid cells that overlap it, colored by the relative overlap of each grid cell with the polygon (NB: this currently only works if :py:meth:`xa.pixel_overlaps` is run with ``subset_to_bbox=False``, or using :py:meth:`xa.subset_find` as detailed in `Detailed Code Run <./notebooks/full_run.ipynb>`_)::

   # Querying polygon by column of the polygon `gdf`
   weightmap.diag_fig({'name':'Alaska'},ds)

   # Plotting the first polygon in the polygon `gdf`
   weightmap.diag_fig(0,ds)




