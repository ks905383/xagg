Tips 
#######################################

Silencing status updates
---------------------------------------

To silence status updates to standard out, use the `silent=True` flag::

   import xagg as xa

   # Get overlap between pixels and polygons
   weightmap = xa.pixel_overlaps(ds,gdf,silent=True)

   # Aggregate data in [ds] onto polygons
   aggregated = xa.aggregate(ds,weightmap,silent=True)

