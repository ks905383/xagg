Tips 
#######################################

Silencing status updates
---------------------------------------

To silence status updates to standard out, use the ``silent=True`` flag::

   import xagg as xa

   # Get overlap between pixels and polygons
   weightmap = xa.pixel_overlaps(ds,gdf,silent=True)

   # Aggregate data in [ds] onto polygons
   aggregated = xa.aggregate(ds,weightmap,silent=True)

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

Note that ``weightmap.to_file(fn)`` creates and populates a separate _directory_ named ``fn`` to be able to store all the relevant components of the ``weightmap`` class, including shapefiles with the geometry of the input polygons, the dataframe with the pixel overlap data, the source grid, and any additional weight grids.

This feature is still slightly experimental, so please let us know your experiences! 

Speed up overlap calculation
---------------------------------------
At the expense of increased memory issue, processing may be sped up using an alternate calculation method (``impl='dot_product'``) :: 

   # Get overlap between pixels and polygons
   weightmap = xa.pixel_overlaps(ds,gdf,impl='dot_product')

   # Aggregate data in [ds] onto polygons
   aggregated = xa.aggregate(ds,weightmap,impl='dot_product')

This feature is still slightly experimental, so please let us know your experiences! 




