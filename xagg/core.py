import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
import warnings
import re 
import os
try:
    import xesmf as xe
    _has_xesmf=True
except ImportError:
    _has_xesmf=False

from . auxfuncs import (find_rel_area,normalize,fix_ds,get_bnds,subset_find,list_or_first)
from . classes import (weightmap,aggregated)

class NoOverlapError(Exception):
    """ Exception for when there's no overlap between pixels and polygons """
    pass


def read_wm(path):
    """ Load temporary weightmap files from wm.to_file()

    Builds a weightmap out of saved weightmap component 
    files. Particularly useful if the weightmap took a 
    particularly long time to calculated (i.e., if the 
    grid is particularly high resolution). 

    Assumes the files are generated from `wm.to_file()`; 
    i.e., the files are all in a directory `name`: 

    `name/name.shp` : the geometry of the input polygons
    `name/name.agg` : the dataframe with the pixel overlap
                      data
    `name/name_lat.nc`, `name/name_lon.nc` : 
                      the source grid of the raster data
    `name/name_weights.nc` :
                      the additional weights grid, if used
                      (this file is optional; if no file
                      with this name is found, no weights 
                      are assumed, and 
                      `wm.weights='noweights'`)

    Parameters
    ---------------
    path : str 
        The directory in which the files are stored. They
        are assumed to follow the filename convention of
        sharing the name of the directory (i.e., the last
        part of this path.)

    Returns
    ---------------
    wm : :class:`xagg.weightmap`

    """
    # the last bit of the path is also the filename
    fn = re.split('/',path)[-1]
    

    ###### Load geometry
    geo = gpd.read_file(path+'/'+fn+'.shp')
    geo = geo['geometry']

    ####### Load agg 
    agg = pd.read_hdf(path+'/'+fn+'.h5', 'wm')

    ###### Load source grid
    source_grid = {k:xr.open_dataset(path+'/'+fn+'_'+k+'.nc').set_index({'loc':('lat','lon')})[k+'v'] for k in ['lat','lon']}
    # Rename, removing the v added for the multi-index issue in export
    source_grid = {k:v.rename(k) for k,v in source_grid.items()}

    ###### Load weights
    if os.path.exists(path+'/'+fn+'_weights.csv'):
        # Specifying column because it saves it as with 
        # a dummy index column that gets loaded in an 
        # unproductive way
        weights = pd.read_csv(path+'/'+fn+'_weights.csv')['weights'].astype(object)
        # ^^ Setting astype(object) to make sure integral weights
        # don't change the general type of the frame. This may
        # only affect the testing routines, but setting this here
        # to be explicit 
    else:
        weights = 'nowghts'

    ###### Combine into weightmap
    wm = weightmap(agg=agg,
                   geometry=geo,
                   source_grid=source_grid,
                   weights=weights)

    return wm


def process_weights(ds,weights=None,target='ds',silent=False):
    """ Process weights - including regridding
    
    If ``target == 'ds'``, regrid `weights` to `ds`. If ``target == 'weights'``,
    regrid `ds` to `weights`. 

    Parameters
    ---------------
    ds : :class:`xarray.Dataset`, :class:`xarray.DataArray`
        an :class:`xarray.Dataset`/:class:`xarray.DataArray` to regrid
    
    weights : :class:`xarray.DataArray`, optional, default = ``None``
        an :class:`xarray.DataArray` containing a weight (numeric) 
        at each location
    
    target : str, optional, default = ``'ds'``
        whether weights should be regridded to the `ds` grid (by
        default) or vice-versa (not yet supported, returns 
        `NotImplementedError`)

    silent : bool, default = `False`
        if True, then no status updates are printed to std out

    Returns
    ---------------
    ds : :class:`xarray.Dataset`, :class:`xarray.DataArrays`
        the input :class:`xarray.Dataset`/:class:`xarray.DataArray`, with a new variable 
        `weights` specifying weights for each pixel

    weights_info : dict 
        a dictionary storing information about the 
        weights regridding process, with the fields:

        - ``target``: showing which of the two grids was retained
        - ``ds_grid``: a dictionary with the grid ``{"lat":ds.lat,"lon",ds.lon}``
        - ``weights_grid``: a dictionary with the grid ``{"lat":weights.lat,"lon":weights.lon}``
    """
    
    if weights is None:
        # (for robustness against running this without an extra if statement
        # in a wrapper function)
        weights_info = 'nowghts'
    else:
        # Check types
        if type(weights) is not xr.core.dataarray.DataArray:
            raise TypeError('[weights] must be an xarray DataArray.')
        if type(ds) not in [xr.core.dataarray.DataArray,
                                xr.core.dataset.Dataset]:
            raise TypeError('[ds] must be an xarray structure (DataArray or Dataset)')
            
        # Stick weights into the same supported input format as ds
        weights = fix_ds(weights)
        
        # Set regridding info
        weights_info = {'target':target,
                        'ds_grid':{'lat':ds.lat,'lon':ds.lon},
                        'weights_grid':{'lat':weights.lat,'lon':weights.lon}}

        # Regrid, if necessary (do nothing if the grids match up to within
        # floating-point precision)
        if ((not ((ds.sizes['lat'] == weights.sizes['lat']) and (ds.sizes['lon'] == weights.sizes['lon']))) or 
            (not (np.allclose(ds.lat,weights.lat) and np.allclose(ds.lon,weights.lon)))):
            # Import xesmf here to allow the code to work without it (it 
            # often has dependency issues and isn't necessary for many 
            # features of xagg)
            if not _has_xesmf:
                raise ImportError('If the `weights` grid and the `ds` grid are different, '+
                                 '`xesmf` is needed for `xagg` to regrid them to match; however, '+
                                 '`xesmf` is not installed. Either install `xesmf` or '+
                                 'manually regrid them to match each other.')
            if target == 'ds':
                if not silent:
                    print('regridding weights to data grid...')
                # Create regridder to the [ds] coordinates
                rgrd = xe.Regridder(weights,ds,'bilinear')
                # Regrid [weights] to [ds] grids
                weights = rgrd(weights)

            elif target == 'weights':
                raise NotImplementedError('The '+target+' variable is not *yet* supported as a target for regridding. Please choose "ds" for now.')
                # This is because of lack of downstream capability right now... 
                if not silent:
                    print('regridding data to weights grid...')
                # Create regridder to the [weights] coordinates
                rgrd = xe.Regridder(ds,weights,'bilinear')
                # Regrid [ds] to [weights] grid
                ds = rgrd(ds)

            else:
                raise KeyError(target+' is not a supported target for regridding. Choose "weights" or "ds".')
        else:
            # Make sure the values are actually identical, not just "close", 
            # otherwise assigning may not work below 
            weights['lat'] = ds['lat'].values
            weights['lon'] = ds['lon'].values

        # Add weights to ds
        ds['weights'] = weights

        # Add warnings
        if np.isnan(ds['weights']).all():
            warnings.warn('All inputted `weights` are np.nan after regridding.')
        if (ds['weights'] == 0).all():
            warnings.warn('All inputted `weights` are 0 after regridding.')
            
    # Return
    return ds,weights_info


def make_multipoly(pts):
    ''' Split pixel overlapping the antimeridian into MultiPolygon with 
        each sub-Polygon in its own hemisphere

    NB: to be used in `create_raster_polygons()`
        
    '''
    pts = np.array(pts)
    # Get which pixels are east of the antimeridian
    neg = pts[:,0]<0

    # Get point order for split up pixel 
    pts = [ # west of antimeridian
            np.vstack([pts[~neg],np.array([[180,x[1]] for x in pts[~neg][::-1]])]),
            # east of antimeridian
           np.vstack([np.array([[-180,x[1]] for x in pts[neg][::-1]]),pts[neg]])]
    pts = [[tuple(x) for x in pt] for pt in pts]

    # Create multipolygon
    return MultiPolygon([Polygon(pt) for pt in pts])



def create_raster_polygons(ds,
                           mask=None,subset_bbox=None,
                           weights=None,weights_target='ds',
                           wrap_around_thresh=5,
                          silent=False):
    """ Create polygons for each pixel in a raster

    Note: 
    'lat_bnds' and 'lon_bnds' can be created through the
    :func:`xagg.aux.get_bnds` function if they are not already included
    in the input raster file. 
    
    Note:
    Currently this code only supports regular 
    rectangular grids (so where every pixel side is
    a straight line in lat/lon space). Future versions
    may include support for irregular grids. 
    
    Parameters
    ---------------
    ds : :class:`xarray.Dataset`
        an xarray dataset with the variables 
        'lat_bnds' and 'lon_bnds', which are both
        lat/lon x 2 arrays giving the min and 
        max values of lat and lon for each pixel
        given by lat/lon
    subset_bbox : :class:`geopandas.GeoDataFrame`, optional, default = ``None`` 
        if a :class:`geopandas.GeoDataFrame`
        is entered, the bounding
        box around the geometries in the `gdf` 
        are used to mask the grid, to reduce
        the number of pixel polygons created
          
    Returns
    ---------------
    pix_agg: dict
        a dictionary containing:

        - ``'gdf_pixels'`` 
            a :class:`geopandas.GeoDataFrame` containing a 
            'geometry' giving the pixel boundaries 
            for each 'lat' / 'lon' pair
        - ``'source_grid'`` 
            a dictionary containing the original
            lat and lon inputs under the keys
            "lat" and "lon" (just the 
            :class:`xarray.DataArray` of those variables in
            the input `ds`)
                      
    """
    
    # Standardize inputs (including lat/lon order)
    ds = fix_ds(ds)
    ds = get_bnds(ds,silent=silent)
    #breakpoint()
    # Subset by shapefile bounding box, if desired
    if subset_bbox is not None:
        if type(subset_bbox) == gpd.geodataframe.GeoDataFrame:
            # Using the biggest difference in lat/lon to make sure that the pixels are subset
            # in a way that the bounding box is fully filled out
            # bbox_thresh = np.max([ds.lat.diff('lat').max(),ds.lon.diff('lon').max()])+0.1 
            grid_dist = np.max([ds.lat.diff('lat').max(),ds.lon.diff('lon').max()]) # first get the max grid size
            bbox_thresh = grid_dist*2. # then set threshold to twice grid size, avoids huge subsets for high res grids
            ds = ds.sel(lon=slice(subset_bbox.total_bounds[0]-bbox_thresh,subset_bbox.total_bounds[2]+bbox_thresh),
                        lat=slice(subset_bbox.total_bounds[1]-bbox_thresh,subset_bbox.total_bounds[3]+bbox_thresh))
        else:
            warnings.warn('[subset_bbox] is not a geodataframe; no mask by polygon bounding box used.')
            
    # Process weights
    ds,winf = process_weights(ds,weights,target=weights_target)
            
    # Mask
    if mask is not None:
        raise NotImplementedError('Masking by grid not yet supported. Stay tuned...')
        
    # Create dataset which has a lat/lon bound value for each individual pixel, 
    # broadcasted out over each lat/lon pair
    (ds_bnds,) = (xr.broadcast(ds.isel({d:0 for d in [k for k in ds.sizes if k not in ['lat','lon','bnds']]}).
                              drop_vars([v for v in ds.keys() if v not in ['lat_bnds','lon_bnds']])))
    # Stack so it's just pixels and bounds
    ds_bnds = ds_bnds.stack(loc=('lat','lon'))
        
    # In order:
    # (lon0,lat0),(lon0,lat1),(lon1,lat1),(lon1,lat1), but as a single array; to be 
    # put in the right format for Polygon in the next step
    pix_poly_coords = np.transpose(np.vstack([ds_bnds.lon_bnds.isel(bnds=0).values,ds_bnds.lat_bnds.isel(bnds=0).values,
                                                ds_bnds.lon_bnds.isel(bnds=0).values,ds_bnds.lat_bnds.isel(bnds=1).values,
                                                ds_bnds.lon_bnds.isel(bnds=1).values,ds_bnds.lat_bnds.isel(bnds=1).values,
                                                ds_bnds.lon_bnds.isel(bnds=1).values,ds_bnds.lat_bnds.isel(bnds=0).values]))

    # Reshape so each location has a 4 x 2 (vertex vs coordinate) array, 
    # and convert each of those vertices to tuples. This means every element
    # of pix_poly_coords is the input to shapely.geometry.Polygon of one pixel
    pix_poly_coords = tuple(map(tuple,np.reshape(pix_poly_coords,(np.shape(pix_poly_coords)[0],4,2))))

    # Figure out if any pixels cross the antimeridian; we'll have to deal with 
    # those separately... Identify them by seeing which pixels have longitudes
    # that are within the `wrap_around_thresh` (by default 5 degs) of both 
    # +180 and -180 degrees
    cross_antimeridian_idxs = ((np.abs(np.array(pix_poly_coords)[:,:,0] - -180) < wrap_around_thresh).any(axis=1) & 
     (np.abs(np.array(pix_poly_coords)[:,:,0] - 180) < wrap_around_thresh).any(axis=1))

    # Create empty geodataframe
    gdf_pixels = gpd.GeoDataFrame(pd.DataFrame({v:[None]*ds_bnds.sizes['loc']
                               for v in ['lat','lon','geometry']}),
                              geometry='geometry')
    if weights is not None:
        # Stack weights so they are linearly indexed like the ds (and fill
        # NAs with 0s)
        weights = ds.weights.stack(loc=('lat','lon')).fillna(0)
        # Preallocate weights column
        gdf_pixels['weights'] = [None]*ds_bnds.sizes['loc']
    
    # Now populate with a polygon for every pixel
    poly_dict = {'poly_pts': pix_poly_coords}
    df_poly = pd.DataFrame(poly_dict, columns=['poly_pts'])
    df_poly['poly'] = df_poly.poly_pts.apply(lambda pts: Polygon(pts))
    # Make MultiPolygons for pixels crossing the antimeridian
    df_poly.loc[np.where(cross_antimeridian_idxs)[0],'poly'] = df_poly.poly_pts.iloc[np.where(cross_antimeridian_idxs)[0]].apply(lambda pts: make_multipoly(pts))
    # Set geometry 
    gdf_pixels['geometry']=df_poly['poly']
    # Add lat/lons of pixels for identification 
    gdf_pixels['lat']=ds_bnds.lat.values
    gdf_pixels['lon']=ds_bnds.lon.values  
    if weights is not None:
        gdf_pixels['weights'] = weights.values

    # Add a "pixel idx" to make indexing better later
    gdf_pixels['pix_idx'] = gdf_pixels.index.values
    
    # Add crs (normal lat/lon onto WGS84)
    gdf_pixels = gdf_pixels.set_crs("EPSG:4326")
    #gdf_pixels.crs = {'init':'EPSG:4326'}
    
    # Save the source grid for further reference
    source_grid = {'lat':ds_bnds.lat,'lon':ds_bnds.lon}
    
    pix_agg = {'gdf_pixels':gdf_pixels,'source_grid':source_grid}
    
    # Return the created geodataframe
    return pix_agg


def get_pixel_overlaps(gdf_in,pix_agg,impl='for_loop'):
    """ Get, for each polygon, the pixels that overlap and their area of overlap
    
    Finds, for each polygon in `gdf_in`, which pixels intersect it, and by how much. 
    
    Note: 
    Uses EASE-Grid 2.0 on the WGS84 datum to calculate relative areas
    (see https://nsidc.org/data/ease)
    
    Parameters
    ---------------
    gdf_in : :class:`geopandas.GeoDataFrane`
        a :class:`geopandas.GeoDataFrame` giving the polygons over which 
        the variables should be aggregated. Can be just a read
        shapefile (with the added column of "poly_idx", which 
        is just the index as a column).
    pix_agg : dict
        the output of :func:`xagg.core.create_raster_polygons`; a dict containing:

        - ``'gdf_pixels'``
            a :class:`geopandas.GeoDataFrame` giving for each row 
            the columns "lat" and "lon" (with coordinates) 
            and a polygon giving the boundary of the pixel 
            given by lat/lon
        - ``'source_grid'``
            ``[da.lat,da.lon]`` of the grid used to create
            the pixel polygons
    impl : str
        whether the output will be used for the dot-product aggregation 
        calculation (needs a slightly different format), either of: 
        - ``'for_loop'`` (default behavior)
        - ``'dot_product'`` (to set up for ``impl='dot_product'`` in 
            ``xagg.core.aggregate``)
    

    Returns
    ---------------
    wm_out: dict
        A dictionary containing: 
        
        - ``'agg'``
            a dataframe containing all the fields of `gdf_in` (except
            geometry) and the additional columns: 

            - `coords`:  the lat/lon coordiates of all pixels that overlap
             the polygon of that row
            - `pix_idxs`: the linear indices of those pixels within the 
             `gdf_pixels` grid
            - `rel_area`: the relative area of each of the overlaps between
             the pixels and the polygon (summing to 1 - e.g. 
             if the polygon is exactly the size and location of
             two pixels, their rel_areas would be 0.5 each)
        - ``'source_grid'``: 
            a dictionary with keys 'lat' and 'lon' giving the 
            original lat/lon grid whose overlaps with the polygons
            was calculated
        - ``'geometry'``: 
            just the polygons from `gdf_in` 


    """
    
    # Add an index for each polygon as a column to make indexing easier
    #if 'poly_idx' not in gdf_in.columns:
    #    gdf_in['poly_idx'] = gdf_in.index.values
    gdf_in['poly_idx'] = np.arange(0,len(gdf_in))
        
    # Match up CRSes
    pix_agg['gdf_pixels'] = pix_agg['gdf_pixels'].to_crs(gdf_in.crs)

    # Choose a common crs for both, just to minimize the chance
    # of geographic shenanigans
    # (using the EASE grid https://nsidc.org/data/ease)
    if np.all(gdf_in.total_bounds[[1,3]]>0):
        # If min/max lat are both in NH, use North grid
        #epsg_set = {'init':'EPSG:6931'} (change to below bc of depreciation of {'init':...} format in geopandas)
        epsg_set = 'EPSG:6931'
    elif np.all(gdf_in.total_bounds[[1,3]]<0):
        # If min/max lat are both in SH, use South grid
        #epsg_set = {'init':'EPSG:6932'}
        epsg_set = 'EPSG:6932'
    else:
        # Otherwise, use the global/temperate grid
        #epsg_set = {'init':'EPSG:6933'}
        epsg_set = 'EPSG:6933'
    
    # Get GeoDataFrame of the overlaps between every pixel and the polygons
    with warnings.catch_warnings():
        # Filter UserWarnings that flag when overlapping would result in 
        # lines, or point overlaps (but we only care about 3D overlaps, so
        # keep_geom_type=True is the right call, and we don't need the 
        # warning)
        warnings.filterwarnings('ignore',category=UserWarning)
        overlaps = gpd.overlay(gdf_in.to_crs(epsg_set),
                               pix_agg['gdf_pixels'].to_crs(epsg_set),
                               how='intersection')
    
    if overlaps.empty:
        raise NoOverlapError('No `ds` grid cells overlapped with any polygon in `gdf_in`. Check the input `ds` and `gdf_in`.')
    else:
        if impl=='dot_product':
            # Get relative area of each pixel
            overlaps = overlaps.groupby('poly_idx',group_keys=False).apply(find_rel_area)
            overlaps['lat'] = overlaps['lat'].astype(float)
            overlaps['lon'] = overlaps['lon'].astype(float)

            # Now, group by poly_idx (each polygon in the shapefile)
            ov_groups = overlaps.groupby('poly_idx')

            overlap_info = ov_groups.agg(list_or_first)

            overlap_info = overlap_info.rename(columns={'pix_idx': 'pix_idxs'})
        elif impl=='for_loop':
            # Now, group by poly_idx (each polygon in the shapefile)
            ov_groups = overlaps.groupby('poly_idx')
            # Calculate relative area of each overlap (so how much of the total 
            # area of each polygon is taken up by each pixel), the pixels 
            # corresponding to those areas, and the lat/lon coordinates of 
            # those pixels
            overlap_info = ov_groups.agg(rel_area=pd.NamedAgg(column='geometry',aggfunc=lambda ds: [normalize(ds.area)]),
                                           pix_idxs=pd.NamedAgg(column='pix_idx',aggfunc=lambda ds: [idx for idx in ds]),
                                           lat=pd.NamedAgg(column='lat',aggfunc=lambda ds: [x for x in ds]),
                                           lon=pd.NamedAgg(column='lon',aggfunc=lambda ds: [x for x in ds]))

    # Zip lat, lon columns into a list of (lat,lon) coordinates
    # (separate from above because as of 12/20, named aggs with 
    # multiple columns is still an open issue in the pandas github)
    overlap_info['coords'] = overlap_info.apply(lambda row: list(zip(row['lat'],row['lon'])),axis=1)
    overlap_info = overlap_info.drop(columns=['lat','lon'])

    # Reset index to make poly_idx a column for merging with gdf_in
    overlap_info = overlap_info.reset_index()

    if impl=='dot_product':
        # Merge in pixel overlaps to the input polygon geodataframe
        overlap_columns = ['pix_idxs', 'rel_area', 'coords', 'poly_idx']
        gdf_in = pd.merge(gdf_in, overlap_info[overlap_columns],'outer', on='poly_idx')

        # make the weight grid an xarray dataset for later dot product
        idx_cols = ['lat', 'lon', 'poly_idx']
        overlap_da = overlaps.set_index(idx_cols)['rel_area'].to_xarray()
        overlap_da = overlap_da.stack(loc=['lat', 'lon'])
        overlap_da = overlap_da.fillna(0)
        wm_out = weightmap(agg=gdf_in.drop('geometry', axis=1),
                   source_grid=pix_agg['source_grid'],
                   geometry=gdf_in.geometry,
                   overlap_da = overlap_da)
    elif impl=='for_loop':
        # Merge in pixel overlaps to the input polygon geodataframe
        gdf_in = pd.merge(gdf_in,overlap_info,'outer')

        # Drop 'geometry' eventually, just for size/clarity
        wm_out = weightmap(agg=gdf_in.drop('geometry',axis=1),
                           source_grid=pix_agg['source_grid'],
                           geometry=gdf_in.geometry)
    
    if 'weights' in pix_agg['gdf_pixels'].columns:
        wm_out.weights = pix_agg['gdf_pixels'].weights

    return wm_out


def aggregate(ds,wm,impl='for_loop',silent=False):
    """ Aggregate raster variable(s) to polygon(s)
    
    Aggregates (N-D) raster variables in `ds` to the polygons
    in `gfd_out` - in other words, gives the weighted average
    of the values in [ds] based on each pixel's relative area
    overlap with the polygons. 
    
    The values will be additionally weighted if a `weight` was
    inputted into :func:`xagg.core.create_raster_polygons`
    
    The code checks whether the input lat/lon grid in `ds` is 
    equivalent to the linearly indexed grid in `wm`, or 
    if it can be cropped to that grid. 
    
    
    Parameters
    ---------------
    ds : :class:`xarray.Dataset`
        an :class:`xarray.Dataset` containing one or more
        variables with dimensions lat, lon (and possibly
        more). The dataset's geographic grid has to 
        include the lat/lon coordinates used in 
        determining the pixel overlaps in 
        :func:`xagg.core.get_pixel_overlaps` (and saved in 
        ``wm['source_grid']``)
               
    wm : :class:`xagg.classes.weightmap`
        the output to :func:`xagg.core.get_pixel_overlaps`; a 
        :class:`xagg.classes.weightmap` object containing 

        - ``['agg']``
            a dataframe, 
            with one row per polygon, and the columns `pix_idxs`
            and `rel_area`, giving the linear indices and 
            the relative area of each pixel over the polygon,
            respectively 
        - ``['source_grid']``  
            the lat/lon grid on which the aggregating parameters
            were calculated (and on which the linear indices 
            are based)

    impl : :class:str (def: ``'for_loop'``)
        which aggregation calculation method to use, either of: 

        - ``'for_loop'`` 
            default behavior, aggregation loops through
            all polygons in a for loop, requires less memory 
        - ``'dot_product'`` 
            aggregation is calculated using a dot product, 
            requires much more memory (due to broadcasting of
            variables) but may be faster in certain circumstances

    silent : bool, default = `False`
        if True, then no status updates are printed to std out
               
    Returns
    ---------------
    agg_out : :class:`xagg.classes.aggregated`
        an :class:`xagg.classes.aggregated` object with the aggregated variables 
    
    """
    # Make sure pixel_overlaps was correctly run if using dot product
    if (impl=='dot_product') and (wm.overlap_da is None):
        raise ValueError("no 'overlap_da' was found in the `wm` input - since you're using the dot product implementation, "+
                         "make sure to run `pixel_overlaps()` with `impl='dot_product'` to avoid this error.")

    # Turn into dataset if dataarray
    if type(ds)==xr.core.dataarray.DataArray:
      if ds.name is None:
        warnings.warn('An unnamed xr.DataArray was inputted instead of a xr.Dataset; the output variable will be "var"')
        ds = ds.to_dataset(name='var')
      else:
        ds = ds.to_dataset()


    # Run ds through fix_ds (to fix lat/lon names, lon coords)
    ds = fix_ds(ds)

    # Stack 
    ds = ds.stack(loc=('lat','lon'))
    
    # Adjust grid of [ds] if necessary to match 
    ds = subset_find(ds,wm.source_grid)
    
    # Set weights; or replace with ones if no additional weight information
    #if wm.weights != 'nowghts':
    if type(wm.weights) == pd.core.series.Series:
        weights = np.array([float(k) for k in wm.weights])
    else:
        if wm.weights != 'nowghts':
            warnings.warn('wm.weights is: \n '+print(wm.weights)+
                            ', \n which is not a supported weight vector (in a pandas series) '+
                            'or "nowghts" as a string. Assuming no weights are included...')
        if impl=='dot_product':
            weights = np.ones((len(wm.overlap_da['loc'])))
        elif impl=='for_loop':
            weights = np.ones((len(wm.source_grid['lat'])))
    
    if impl=='dot_product':
        data_dict = dict()
        for var in ds:
            # Process for every variable that has locational information, but isn't a 
            # bound variable
            if ('bnds' not in ds[var].sizes) & ('loc' in ds[var].sizes):
                if not silent:
                    print('aggregating '+var+'...')

                # select just data for which we have overlaps
                var_array = ds[var].sel(loc=wm.overlap_da['loc'])

                # multiply percent-overlaps by user-supplied weights 
                weights_and_overlaps = wm.overlap_da * weights

                # fill any nans in the gridded data to zero and change the corresponding
                # weights to zero so that they will not affect the aggregated value
                var_array_filled = var_array.fillna(0)
                weights_and_overlaps = weights_and_overlaps.where(var_array.notnull(), 0)

                # the weights now may add up to less than one because of nans or
                # because of the user-supplied weights. here we normalize the
                # weights so they add up to 1 and fill any nan's with 0 so they
                # won't be counted
                normed_weights = weights_and_overlaps/weights_and_overlaps.sum(dim = "loc", skipna=True)
                normed_weights = normed_weights.fillna(0)

                # finally we do the dot product to get the weighted averages
                aggregated_array = normed_weights.dot(var_array_filled)

                # if the original gridded values were all nan, make the final
                # aggregation nan
                if np.isnan(var_array.values).all():
                    aggregated_array = aggregated_array * np.nan

                data_dict[var] = aggregated_array

        ds_combined = xr.Dataset(data_dict)    
        df_combined = ds_combined.to_dataframe().reset_index()
        df_combined = df_combined.groupby('poly_idx').agg(list_or_first)

        wm.agg = pd.merge(wm.agg, df_combined, on='poly_idx')
        for var in ds:
            if ('bnds' not in ds[var].sizes) & ('loc' in ds[var].sizes):
                # convert to list of arrays - NOT SURE THIS IS THE RIGHT THING TO
                # DO, JUST TRYING TO MATCH ORIGINAL FORMAT
                wm.agg[var] = wm.agg[var].apply(np.array).apply(lambda x: [x])

        # Put in class format
        agg_out = aggregated(agg=wm.agg,source_grid=wm.source_grid,
                            geometry=wm.geometry,ds_in=ds_combined,weights=wm.weights)
    elif impl=='for_loop':
        for var in ds:
            # Process for every variable that has locational information, but isn't a 
            # bound variable
            if ('bnds' not in ds[var].sizes) & ('loc' in ds[var].sizes):
                if not silent:
                    print('aggregating '+var+'...')
                # Create the column for the relevant variable
                wm.agg[var] = None
                
                # Get weighted average of variable based on pixel overlap + other weights
                for poly_idx in wm.agg.poly_idx:
                    # Get average value of variable over the polygon; weighted by 
                    # how much of the pixel area is in the polygon, and by (optionally)
                    # a separate gridded weight. This if statement checks 
                    # - if the weights output for this polygon is just "np.nan", which 
                    #   indicates there's no overlap between the pixel grid and polygons
                    # - [this has been pushed down a few lines] or, if all the pixels in 
                    #   the grid have just nan values for this variable
                    # in both cases; the "aggregated variable" is just a vector of nans. 
                    if not np.isnan(wm.agg.iloc[poly_idx,:].pix_idxs).all():
                        # Get the dimensions of the variable that aren't "loc" (location)
                        other_dims = [k for k in np.atleast_1d(ds[var].sizes) if k != 'loc']
                        # Check first if any nans are "complete" (meaning that a pixel 
                        # either has values for each step, or nans for each step - if
                        # there are random nans within a pixel, throw a warning)
                        if not xr.Dataset.equals(np.isnan(ds[var].isel(loc=wm.agg.iloc[poly_idx,:].pix_idxs)).any(other_dims),
                                                  np.isnan(ds[var].isel(loc=wm.agg.iloc[poly_idx,:].pix_idxs)).all(other_dims)):
                            warnings.warn('One or more of the pixels in variable '+var+' have nans in them in the dimensions '+
                                           ', '.join([k for k in ds[var].sizes if k != 'loc'])+
                                           '. The code can currently only deal with pixels for which the '+
                                           '*entire* pixel has nan values in all dimensions, however there is currently no '+
                                           ' support for data in which pixels have only some nan values. The aggregation '+
                                           'calculation is likely incorrect.')

                        if not np.isnan(ds[var].isel(loc=wm.agg.iloc[poly_idx,:].pix_idxs)).all(): 
                            # Get relative areas for the pixels overlapping with this Polygon
                            tmp_areas = np.atleast_1d(np.squeeze(wm.agg.iloc[poly_idx,:].rel_area))
                            # Replace overlapping pixel areas with nans if the corresponding pixel
                            # is only composed of nans
                            tmp_areas[np.array(np.isnan(ds[var].isel(loc=wm.agg.iloc[poly_idx,:].pix_idxs)).all(other_dims).values)] = np.nan
                            # Calculate the normalized area+weight of each pixel (taking into account
                            # nans)
                            normed_areaweights = normalize(tmp_areas*weights[wm.agg.iloc[poly_idx,:].pix_idxs],drop_na=True)

                            # Take the weighted average of all the pixel values to calculate the 
                            # aggregated value for the shapefile
                            wm.agg.loc[poly_idx,var] = [[((ds[var].isel(loc=wm.agg.iloc[poly_idx,:].pix_idxs)*
                                                             normed_areaweights).
                                                            sum('loc')).values]]

                        else:
                            wm.agg.loc[poly_idx,var] = [[(ds[var].isel(loc=0)*np.nan).values]]

                    else:
                        wm.agg.loc[poly_idx,var] = [[(ds[var].isel(loc=0)*np.nan).values]]

        # Put in class format
        agg_out = aggregated(agg=wm.agg,source_grid=wm.source_grid,
        					 geometry=wm.geometry,ds_in=ds,weights=wm.weights)

    # Return
    if not silent:
        print('all variables aggregated to polygons!')
    return agg_out

