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
:py:class:`xagg` has a few backend algorithms for aggregation. By default, :py:class:`xagg` uses the ``'for_loop'`` implementation, which minimizes memory use, but can be very slow for certain large datasets. Alternatively, you can use: 

- ``impl = 'dot_product'``: faster, but at the expense of increased memory usage
- ``impl = 'numba'``: likely fastest for large datasets once compiled, uses the just-in-time compiler from `numba` in ``xa.aggregate()``. Requires :py:class:`numba` to be installed. 


Switching to these algorithms is best achieved through setting defaults: :: 

   import xagg as xa
   # Set dot_product as default backend implementation 
   xa.set_defaults(impl='dot_product')

   # Set numba as default backend implementation
   xa.set_defaults(impl='numba')

Alternatively, defaults can be set for an individual ``with`` block: ::
   with xa.set_defaults(impl='...'):
      wm = xa.pixel_overlaps(ds,gdf)
      agg = xa.aggregate(ds,wm)

Note that the `impl=...` option in individual functions will be
slowly deprecated over the next few versions.

Create diagnostic figure to inspect raster/polygon overlaps 
------------------------------------------------------------
Once you have created a :py:class:`weightmap`, the following code will create a diagnostic figure, showing a particular polygon (or groups of polygons) + the grid cells that overlap it, colored by the relative overlap of each grid cell with the polygon (NB: this currently only works if :py:meth:`xa.pixel_overlaps` is run with ``subset_to_bbox=False``, or using :py:meth:`xa.subset_find` as detailed in `Detailed Code Run <./notebooks/full_run.ipynb>`_)::

   # Querying polygon by column of the polygon `gdf`
   weightmap.diag_fig({'name':'Alaska'},ds)

   # Plotting the first polygon in the polygon `gdf`
   weightmap.diag_fig(0,ds)

IndexErrors in :py:meth:`xa.pixel_overlaps`
------------------------------------------------------------
If you're running into an `IndexError` when running :py:meth:`xa.pixel_overlaps` (e.g., `IndexError: too many indices for array: array is 1-dimensional, but 3 were indexed`), try reprojecting the input :py:meth:`geodataframe` to `'EPSG:4326'` before running :py:meth:`xa.pixel_overlaps`. See `here <https://github.com/ks905383/xagg/issues/80>`_ for more discussion on this issue.

Non-rectangular grids
------------------------------------------------------------
:py:mod:`xagg` unfortunately currently only works with rectangular grids. 



