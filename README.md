xagg
==============================
[![Build Status](https://github.com/ks905383/xagg/workflows/Tests/badge.svg)](https://github.com/ks905383/xagg/actions)
[![codecov](https://codecov.io/gh/ks905383/xagg/branch/main/graph/badge.svg)](https://codecov.io/gh/ks905383/xagg)
[![pypi](https://img.shields.io/pypi/v/xagg.svg)](https://pypi.org/project/xagg)
[![conda-forge](https://anaconda.org/conda-forge/xagg/badges/version.svg)](https://anaconda.org/conda-forge/xagg/)
[![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/xagg.svg)](https://anaconda.org/conda-forge/xagg)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.07239/status.svg)](https://doi.org/10.21105/joss.07239)
[![DOI](https://zenodo.org/badge/367399550.svg)](https://zenodo.org/doi/10.5281/zenodo.11476311)
[![Documentation Status](https://readthedocs.org/projects/xagg/badge/?version=latest)](https://xagg.readthedocs.io/en/latest/?badge=latest)



A package to aggregate gridded data in `xarray` to polygons in `geopandas` using area-weighting from the relative area overlaps between pixels and polygons.

## Installation 
The easiest way to install the latest version of `xagg` is using `conda` or `mamba`:  

```
conda install -c conda-forge xagg==0.3.3.1

# or
mamba install -c conda-forge xagg==0.3.3.1
```

We recommend installing `xagg` in a new environment whenever possible, to ensure all (sub)dependencies are correctly loaded. 

Alternatively, you can use `pip`, though not all optional dependencies are available through `pip`, meaning that certain features may not be available:
```
pip install xagg
```

## Documentation 

See the latest documentation at https://xagg.readthedocs.io/en/latest/index.html

## Intro 

Science often happens on grids - gridded weather products, interpolated pollution data, night time lights, remote sensing all approximate the continuous real world for reasons of data resolution, processing time, or ease of calculation.

However, living things don't live on grids, and rarely play, act, or observe data on grids either. Instead, humans tend to work on the county, state, township, Bezirk, or city level; birds tend to fly along complex migratory corridors; and rain- and watersheds follow valleys and mountains. 

So, whenever we need to work with both gridded and geographic data products, we need ways of getting them to match up. We may be interested for example what the average temperature over a county is, or the average rainfall rate over a watershed. 

Enter `xagg`. 

`xagg` provides an easy-to-use (2 lines!), standardized way of aggregating raster data to polygons. All you need is some gridded data in an `xarray` Dataset or DataArray and some polygon data in a `geopandas` GeoDataFrame. Both of these are easy to use for the purposes of `xagg` - for example, all you need to use a shapefile is to open it:

```
   import xarray as xr
   import geopandas as gpd
    
   # Gridded data file (netcdf/climate data)
   ds = xr.open_dataset('file.nc')

   # Shapefile
   gdf = gpd.open_dataset('file.shp')
```

`xagg` will then figure out the geographic grid (lat/lon) in `ds`, create polygons for each pixel, and then generate intersects between every polygon in the `GeoDataFrame` and every pixel. For each polygon in the `GeoDataFrame`, the relative area of each covering pixel is calculated - so, for example, if a polygon (say, a US county) is the size and shape of a grid pixel, but is split halfway between two pixels, the weight for each pixel will be 0.5, and the value of the gridded variables on that polygon will just be the average of both. 

Here is a sample code run, using the loaded files from above: 

```

   import xagg as xa

   # Get overlap between pixels and polygons
   weightmap = xa.pixel_overlaps(ds,gdf)

   # Aggregate data in [ds] onto polygons
   aggregated = xa.aggregate(ds,weightmap)

   # aggregated can now be converted into an xarray dataset (using aggregated.to_dataset()), 
   # or a geopandas geodataframe (using aggregated.to_geodataframe() or aggregated.to_dataframe()
   # for a pure pandas result), or directly exported to netcdf, csv, or shp files using
   # aggregated.to_csv()/.to_netcdf()/.to_shp()
```

Researchers often need to weight your data by more than just its relative area overlap with a polygon (for example, do you want to weight pixels with more population more?). `xagg` has a built-in support for adding an additional weight grid (another `xarray` DataArray) into `xagg.pixel_overlaps()`. 

Finally, `xagg` allows for direct exporting of the aggregated data in several commonly used data formats:

- NetCDF 
- CSV for STATA, R
- Shapefile for QGIS, further spatial processing

Best of all, `xagg` is flexible. Multiple variables in your dataset? `xagg` will aggregate them all, as long as they have at least `lat/lon` dimensions. Fields in your shapefile that you'd like to keep? `xagg` keeps all attributes/fields (for example FIPS codes from county datasets) all the way through the final export. Weird dimension names? `xagg` is trained to recognize all versions of "lat", "Latitude", "Y", "nav_lat", "Latitude_1"... etc. that the author has run into over the years of working with climate data; and this list is easily expandable as a keyword argument if needed. 

## How to support `xagg`
The easiest way to support `xagg` is to star the repository and spread the word!

Please also consider citing `xagg` if you use it in your research. The preferred citation can be found at the "Cite this repository" button in the About section on the top right of this page. It links to our [paper](https://doi.org/10.21105/joss.07239) in the Journal of Open Source Software (JOSS). 

`xagg`, like much of open-source software, is a volunteer-run effort. It means a lot to the developers if you reach out and tell us that you're using our software, how it's helped you, and how it can be improved - it makes the long hours fixing bugs feel that much more worth it. (If you're feeling particularly generous, the lead developer would not say no to additional thanks through [contributions to his tea fund through Ko-Fi](https://ko-fi.com/ks905383) ;) ) 

## Getting Help and Contributing
If you have any questions about how to use `xagg`, please ask them in the [GitHub Discussions](https://github.com/ks905383/xagg/discussions) forum!

If you spot a bug (`xagg` not working as advertised), please [open an issue](https://github.com/ks905383/xagg/issues) if it hasn't yet been raised (or comment on an existing one if you see it listed already). To make sure the issue gets solved as quickly as possible: 
- Include a [minimally reproducible example](https://stackoverflow.com/help/minimal-reproducible-example) that triggers the bug
- Include a copy of your environment (for example, the output of `conda list`) in which the bug occurred

If you'd like to go the extra mile and help us fix the bug, feel free to [contribute a pull request](https://github.com/ks905383/xagg/pulls)! We ask that any PR: 
- Follows a standard development workflow, like [this](https://docs.xarray.dev/en/stable/contributing.html#development-workflow) one. 
- If fixing a bug, [includes unit tests](https://stackoverflow.com/questions/3258733/new-to-unit-testing-how-to-write-great-tests) that fail when confronted with the original bug. GitHub Actions are set up to automatically run all tests in `xagg/tests/` upon a push.

If there's a feature that you'd like `xagg` to have, please start a Discussion in the [GitHub Discussions](https://github.com/ks905383/xagg/discussions) forum, or implement it yourself in a pull request.  

For more information on contributing in general, the [contribution guidelines](https://docs.xarray.dev/en/stable/contributing.html) to the `xarray` package are a great starting point (not everything will be directly relevant to `xagg`, but much of this guide is generally relevant!). 

## Use Cases

### Climate econometrics
Many climate econometrics studies use societal data (mortality, crop yields, etc.) at a political or administrative level (for example, counties) but climate and weather data on grids. Oftentimes, further weighting by population or agricultural density is needed. 

Area-weighting of pixels onto polygons ensures that aggregating weather and climate data onto polygons occurs in a robust way. Consider a (somewhat contrived) example: an administrative region is in a relatively flat lowlands, but a pixel that slightly overlaps the polygon primarily covers a wholly different climate (mountainous, desert, etc.). Using a simple mask would weight that pixel the same, though its information is not necessarily relevant to the climate of the region. Population-weighting may not always be sufficient either; consider Los Angeles, which has multiple significantly different climates, all with high densities. 

`xagg` allows a simple population *and* area-averaging, in addition to export functions that will turn the aggregated data into output easily used in STATA or R for further calculations. 

--------

<p><small>Project based on the <a target="_blank" href="https://github.com/jbusecke/cookiecutter-science-project">cookiecutter science project template</a>.</small></p>
