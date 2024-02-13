import warnings
import xarray as xr
import copy

from . core import (create_raster_polygons,get_pixel_overlaps)


def pixel_overlaps(ds,gdf_in,
                   weights=None,weights_target='ds',
                   subset_bbox = True,
                   impl='for_loop',
                   silent=False):
    """ Wrapper function for determining overlaps between grid and polygon
    
    For a geodataframe `gdf_in`, takes an `xarray` structure `ds` (Dataset or 
    DataArray) and for each polygon in `gdf_in` provides a list of pixels 
    given by the `ds` grid which overlap that polygon, in addition to their
    relative area of overlap with the polygon. 
    
    The output is then ready to be fed into :func:`xagg.core.aggregate`, which aggregates
    the variables in `ds` to the polygons in `gdf_in` using area- (and 
    optionally other) weights. 
    
    (NB: the wrapper uses ``subset_bbox = True`` in :func:`xagg.core.create_raster_polygons`)

    Parameters
    ---------------
    ds : :class:`xarray.Dataset`, :class:`xarray.DataArray`
      an xarray Dataset or DataArray containing at least
      grid variables ("lat"/"lon", though several other names
      are supported; see docs for :func:`xagg.aux.fix_ds`) and at least 
      one variable on that grid

    gdf_in : :class:`geopandas.GeoDataFrame`
      a geopandas GeoDataFrame containing polygons (and 
      any other fields, for example fields from shapefiles)

    weights : :class:`xarray.DataArray` or ``None``, optional, default = ``None``
      (by default, None) if additional weights are desired,
      (for example, weighting pixels by population in addition
      to by area overlap), `weights` is an :class:`xarray.DataArray`
      containing that information. It does *not* have to 
      be on the same grid as `ds` - grids will be homogonized
      (see below).

    weights_target : :class:`str`, optional
      if 'ds', then weights are regridded to the grid in [ds];
      if 'weights', then the `ds` variables are regridded to
      the grid in 'weights' (LATTER NOT SUPPORTED YET, raises
      a `NotImplementedError`)

    impl : :class:`str`, optional, default = ``'for_loop'``
      whether to use the ``'for_loop'`` or ``'dot_product'``
      methods for aggregating; the former uses less memory,
      the latter may be faster in certain circumstances

    silent : bool, default = `False`
        if True, then no status updates are printed to std out
   
    Returns
    ---------------
    wm_out : dict
      the output of :func:`xagg.core.get_pixel_overlaps`  which gives the mapping of pixels to polygon aggregation; to be
      input into :func:`xagg.core.aggregate`. 
    
    """
    # Create deep copy of gdf to ensure the input gdf doesn't 
    # get modified
    gdf_in = copy.deepcopy(gdf_in)

    # Turn into dataset if dataarray
    if type(ds)==xr.core.dataarray.DataArray:
      if ds.name is None:
        ds = ds.to_dataset(name='var')
      else:
        ds = ds.to_dataset()
    
    # Create a polygon for each pixel
    if not silent:
      print('creating polygons for each pixel...')
    if subset_bbox:
        pix_agg = create_raster_polygons(ds,subset_bbox=gdf_in,weights=weights,silent=silent)
        if not silent:
          if pix_agg['gdf_pixels'].empty:
            warnings.warning('Bounding box around polygon(s) from `gdf_in` includes no grid cells in `ds`...')
    else:
        pix_agg = create_raster_polygons(ds,subset_bbox=None,weights=weights,silent=silent)
    
    # Get overlaps between these pixel polygons and the gdf_in polygons
    if not silent:
      print('calculating overlaps between pixels and output polygons...')
    wm_out = get_pixel_overlaps(gdf_in,pix_agg,impl=impl)
    if not silent:
      print('success!')
    
    return wm_out


