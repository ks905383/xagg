Intro
#######################################

Overview
=======================================
Science often happens on grids - gridded weather products, interpolated pollution data, night time lights, remote sensing all approximate the continuous real world for reasons of data resolution, processing time, or ease of calculation.

However, living things don't live on grids, and rarely play, act, or observe data on grids either. Instead, humans tend to work on the county, state, township, Bezirk, or city level; birds tend to fly along complex migratory corridors; and rain- and watersheds follow valleys and mountains. 

So, whenever we need to work with both gridded and geographic data products, we need ways of getting them to match up. We may be interested for example what the average temperature over a county is, or the average rainfall rate over a watershed. 

Enter :py:mod:`xagg`. 

:py:mod:`xagg` provides an easy-to-use (2 lines!), standardized way of aggregating raster data to polygons. All you need is some gridded data in an :class:`xarray.Dataset` or :class:`xarray.DataArray` and some polygon data in a :class:`geopandas.GeoDataFrame`. Both of these are easy to use for the purposes of :py:mod:`xagg` - for example, all you need to use a shapefile is to open it::

   import xarray as xr
   import geopandas as gpd
    
   # Gridded data file (netcdf/climate data)
   ds = xr.open_dataset('file.nc')

   # Shapefile
   gdf = gpd.open_dataset('file.shp')


:py:mod:`xagg` will then figure out the geographic grid (lat/lon) in ``ds``, create polygons for each pixel, and then generate intersects between every polygon in the shapefile and every pixel. For each polygon in the shapefile, the relative area of each covering pixel is calculated - so, for example, if a polygon (say, a US county) is the size and shape of a grid pixel, but is split halfway between two pixels, the weight for each pixel will be 0.5, and the value of the gridded variables on that polygon will just be the average of both. 

The two lines mentioned before? ::

   import xagg as xa

   # Get overlap between pixels and polygons
   weightmap = xa.pixel_overlaps(ds,gdf)

   # Aggregate data in [ds] onto polygons
   aggregated = xa.aggregate(ds,weightmap)

The :py:class:`aggregated` object can now be converted into an xarray dataset (using :py:meth:`aggregated.to_dataset()`) or a :py:mod:`geopandas` :py:class:`GeoDataFrame` (using :py:meth:`aggregated.to_dataframe()`), or directly exported to netcdf, csv, or shp files using :py:meth:`aggregated.to_csv()`/:py:meth:`.to_netcdf()`/:py:meth:`.to_shp()`


Researchers often need to weight your data by more than just its relative area overlap with a polygon (for example, do you want to weight pixels with more population more?). :py:mod:`xagg` has a built-in support for adding an additional weight grid (another :py:mod:`xarray` :py:class:`DataArray`) into :py:func:`pixel_overlaps`. 

Finally, :py:mod:`xagg` allows for direct exporting of the aggregated data in several commonly used data formats (please open issues if you'd like support for something else!):

- NetCDF 
- CSV for STATA, R
- Shapefile for QGIS, further spatial processing

Best of all, :py:mod:`xagg` is flexible. Multiple variables in your dataset? :py:mod:`xagg` will aggregate them all, as long as they have at least ``lat/lon`` dimensions. Fields in your shapefile that you'd like to keep? :py:mod:`xagg` keeps all fields (for example FIPS codes from county datasets) all the way through the final export. Weird dimension names? :py:mod:`xagg` is trained to recognize all versions of "lat", "Latitude", "Y", "nav_lat", "Latitude_1"... etc. that the author has run into over the years of working with climate data; and this list is easily expandable as a keyword argument if needed. 

How to support `xagg`
=======================================
The easiest way to support `xagg` is to star the `GitHub repository <https://github.com/ks905383/xagg>`_ and spread the word!

Please also consider citing `xagg` if you use it in your research. The preferred citation can be found at the "Cite this repository" button in the About section on the top right of the `GitHub repository <https://github.com/ks905383/xagg>`_. 

`xagg`, like much of open-source software, is a volunteer-run effort. It means a lot to the developers if you reach out and tell us that you're using our software, how it's helped you, and how it can be improved - it makes the long hours fixing bugs feel that much more worth it. (If you're feeling particularly generous, the lead developer would not say no to additional thanks through [contributions to his tea fund through Ko-Fi](ko-fi.com/ks905383) ;) ) 

Getting Help and Contributing
=======================================
If you have any questions about how to use `xagg`, please ask them in the `GitHub Discussions <https://github.com/ks905383/xagg/discussions>`_ forum!

If you spot a bug (`xagg` not working as advertised), please `open an issue <https://github.com/ks905383/xagg/issues>`_ if it hasn't yet been raised (or comment on an existing one if you see it listed already). To make sure the issue gets solved as quickly as possible: 
- Include a `minimally reproducible example <https://stackoverflow.com/help/minimal-reproducible-example>`_ that triggers the bug
- Include a copy of your environment (for example, the output of `conda list`) in which the bug occurred

If you'd like to go the extra mile and help us fix the bug, feel free to `contribute a pull request <https://github.com/ks905383/xagg/pulls>`_! We ask that any PR: 
- Follows a standard development workflow, like `this <https://docs.xarray.dev/en/stable/contributing.html#development-workflow>`_ one. 
- If fixing a bug, `includes unit tests <https://stackoverflow.com/questions/3258733/new-to-unit-testing-how-to-write-great-tests>`_ that fail when confronted with the original bug. GitHub Actions are set up to automatically run all tests in `xagg/tests/` upon a push.

If there's a feature that you'd like `xagg` to have, please start a Discussion in the [GitHub Discussions](https://github.com/ks905383/xagg/discussions) forum, or implement it yourself in a pull request.  

For more information on contributing in general, the `contribution guidelines <https://docs.xarray.dev/en/stable/contributing.html>`_ to the `xarray` package are a great starting point (not everything will be directly relevant to `xagg`, but much of this guide is generally relevant!). 


Use cases
=======================================

Climate econometrics
--------------------------------------
Many climate econometrics studies use societal data (mortality, crop yields, etc.) at a political or administrative level (for example, counties) but climate and weather data on grids. Oftentimes, further weighting by population or agricultural density is needed. 

Area-weighting of pixels onto polygons ensures that aggregating weather and climate data onto polygons occurs in a robust way. Consider a (somewhat contrived) example: an administrative region is in a relatively flat lowlands, but a pixel that slightly overlaps the polygon primarily covers a wholly different climate (mountainous, desert, etc.). Using a simple mask would weight that pixel the same, though its information is not necessarily relevant to the climate of the region. Population-weighting may not always be sufficient either; consider Los Angeles, which has multiple significantly different climates, all with high densities. 

:py:mod:`xagg` allows a simple population *and* area-averaging, in addition to export functions that will turn the aggregated data into output easily used in STATA or R for further calculations.
