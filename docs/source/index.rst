.. xagg documentation master file, created by
   sphinx-quickstart on Thu Feb 15 16:33:56 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

xagg
=======================================

A package to aggregate gridded data in :py:mod:`xarray` to polygons in :py:mod:`geopandas` using area-weighting from the relative area overlaps between pixels and polygons.

The simplest code run, involving raster data in an :py:mod:`xarray` :py:class:`Dataset` ``ds`` and polygons in a :py:mod:`geopandas` :py:class:`GeoDataFrame` ``gpd``, is::

 import xagg as xa
    
 # Get overlap between pixels and polygons
 weightmap = xa.pixel_overlaps(ds,gdf)

 # Aggregate data in [ds] onto polygons
 aggregated = xa.aggregate(ds,weightmap)


:py:class:`aggregated` can then be turned into an :py:mod:`xarray` :py:class:`Dataset`, a :py:mod:`geopandas` :py:class:`GeoDataFrame`, or directly exported to a CSV (for use in, e.g., STATA), NetCDF, or Shapefile.

.. toctree::
   :caption: Contents

   intro
   installation
   notebooks/base_run
   notebooks/full_run.ipynb
   tips
   modules


Indices and tables
=================================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
