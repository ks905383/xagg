import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
import warnings
import os
import re

def normalize(a,drop_na = False):
    """ Normalizes the vector `a`

    The vector `a` is divided by its sum. If all non-`np.nan` elements of `a` are 0, 
    then all `np.nan`s are returned.  

    Parameters
    ---------------
    a : array_like
      A vector to be normalized

    drop_na : bool, optional, default = ``False``
      If drop_na = True, and there are nans in the vector a, then the normalization
      is calculated using only the non-nan locations in a, and the vector is returned
      with the nans in their original location. 
      In other words, `np.nansum(normalize(a),drop_na=True) == 1.0`

      If drop_na = False, and nans are present in a, then normalize just returns a 
      vector of np.nan the same length of a.
    
    Returns
    ---------------
    a : vector
      `a`, but normalized. 
  
    """
    
    if (drop_na) and (np.any(np.isnan(a))):
        a2 = a[~np.isnan(a)]
        if np.all(a2.sum()==0):
            # Return nans if the vector is only 0s
            # (/ 0 error)
            return a*np.nan

        else:
            a2 = a2/a2.sum()
            a[~np.isnan(a)] = a2

            return a

    elif (np.all(~np.isnan(a))) and (not np.all(a.sum()==0)):
        return a/a.sum()
    else:
        return a*np.nan

def find_rel_area(df):
    """ 
    Find the relative area of each row in a geodataframe
    """
    df['rel_area'] = df.area/df.area.sum()
    return df


def list_or_first(ser):
    lis = list(ser)
    # only the columns associated with the pixels should have multiple values;
    # for all other columns (those associated with the polygons), it should be
    # safe to return just the first item
    if all(x == lis[0] for x in lis) and ser.name not in ['pix_idx', 'coords', 'rel_area', 'lat', 'lon']:
        return lis[0]
    else:
        return lis


def fix_ds(ds,var_cipher = {'latitude':{'latitude':'lat','longitude':'lon'},
                            'Latitude':{'Latitude':'lat','Longitude':'lon'},
                            'Lat':{'Lat':'lat','Lon':'lon'},
                            'latitude_1':{'latitude_1':'lat','longitude_1':'lon'},
                            'nav_lat':{'nav_lat':'lat','nav_lon':'lon'},
                            'Y':{'Y':'lat','X':'lon'},
                            'y':{'y':'lat','x':'lon'}},
          chg_bnds = True):
    """ Puts the input `ds` into a format compatible with the rest of the package
    
    1) grid variables are renamed "lat" and "lon"
    2) the lon dimension is made -180:180 to be consistent with most geographic
       data (as well as any lon_bnds variable if ``chg_bnds=True`` (by default))
    3) the dataset is sorted in ascending order in both lat and lon


    NOTE: there probably should be a safeguard in case "y" and "x" are multiindex 
    dimension names instead of lat/lon names... maybe a warning for now... (TO DO)
    
    
    Parameters
    ---------------
    ds : :class:`xarray.Dataset`
      an input :class:`xarray.Dataset`, which may or may not need adjustment to be 
      compatible with this package
    var_cipher : dict, optional
      a dict of dicts for renaming lat/lon variables to "lat"/"lon". The 
      form is ``{search_str:{lat_name:'lat',lon_name:'lon'},...}``, 
      the code looks for `search_str` in the dimensions of the `ds`; 
      based on that, it renames `lat_name` to 'lat' and `lon_name` to 'lon'.
      Common names for these variables ('latitude', 'Latitude', 'Lat',
      'latitude_1','nav_lat','Y') are included out of the box.
    chg_bnds  : bool, optional, default = ``True``
      if ``True``, the names of variables with "_bnd" in their names
      are assumed to be dimension bound variables, and are changed as well if
      the rest of their name matches 'o' (for lon) or 'a' (for lat). 
      ## DOES THIS WORK FOR "X" and "Y"?
                  
    Returns
    ---------------
    ds : :class:`xarray.Dataset`
      a dataset with lat/lon variables in the format necessary for this package to function
    
    
    """
    
    # List of variables that represent bounds
    if type(ds) == xr.core.dataset.Dataset:
        bnd_vars = [k for k in list(ds) if 'bnds' in k]
        # Ignore time bounds, since not needed for xagg
        bnd_vars = [k for k in bnd_vars if 'time' not in k]
    elif type(ds) == xr.core.dataarray.DataArray:
        bnd_vars = []
    else:
        raise TypeError('[ds] needs to be an xarray structure (Dataset or DataArray).')
    
    # Fix lat/lon variable names (sizes instead of dims to be compatible with both ds, da...)
    if 'lat' not in ds.sizes:
        test_dims = [k for k in var_cipher if k in ds.sizes]
        if len(test_dims) == 0:
            raise NameError('No valid lat/lon variables found in the dataset.')
        else:
            # should there be a elif len()>1? If there are multiple things found?
            # Could be an issue in terms of ones where x/y are the coordinates of a 
            # non-rectangular grid with variables lat, lon, which is problematic,
            # or just duplicate dimension names, which is weird. 
            ds = ds.rename(var_cipher[test_dims[0]])
        
        # Now same for the bounds variable, if they exist
        if chg_bnds:
            if len(bnd_vars)>0:
                try:
                    ds = ds.rename({(key+'_bnds'):(value+'_bnds') for (key,value) in var_cipher[test_dims[0]].items()})
                except ValueError:
                    try: 
                        warnings.warn('Assuming dim '+[k for k in bnd_vars if 'o' in k][0]+' is longitude bounds and '+
                                      ' dim '+[k for k in bnd_vars if 'a' in k][0] + ' is latitude bounds.')
                        ds = ds.rename({[k for k in bnd_vars if 'o' in k][0]:'lon_bnds',
                                        [k for k in bnd_vars if 'a' in k][0]:'lat_bnds'})
                    except:
                        warnings.warn('Could not identify which of the following bounds '+
                                      'variables corresponds to lat/lon grid: '+', '.join(bnd_vars)+
                                     '; no bound variables renamed.')
                        

    # Switch longitude to -180:180
    if ds.lon.max()>180:
        ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
        
    # Do the same for a longitude bound, if present (the check for dataset is to 
    # avoid an error for .keys() for a dataarray; it seems like a good assumption 
    # that you're not running a dataarray of just lon_bnds through this). 
    if (type(ds) == xr.core.dataset.Dataset):
        # Three if statements because of what I believe to be a jupyter error
        # (where all three statements are evaluated instead of one at a time)
        if ('lon_bnds' in ds):
            if (ds.lon_bnds.max()>180):
                ds['lon_bnds'] = (ds['lon_bnds'] + 180) % 360 - 180

    # Sort by lon; this should be robust (and necessary to avoid fix_ds(fix_ds(ds)) 
    # from failing)
    # PERHAPS THERE SHOULD BE A TEST HERE - like `if ds.lon ~= ds.lon.sortby(ds.lon)`
    # Because currently this blows up memory use... 
    if not np.allclose(ds.lon.values,np.sort(ds.lon)):
        ds = ds.sortby(ds.lon)
    
    # Sort by latitude as well (to make latitude consistent)
    if not np.allclose(ds.lat.values,np.sort(ds.lat)):
        ds = ds.sortby(ds.lat)
        
    # Return fixed ds
    return ds


def get_bnds(ds,wrap_around_thresh='dynamic',
             break_window_width=3,
             break_thresh_x=2,
             silent=False):

    """ Builds vectors of lat/lon bounds if not present in `ds`
    
    Assumes a regular rectangular grid - so each lat/lon bound
    is 0.5*(gap between pixels) over to the next pixel.
    
    
    Parameters
    ---------------
    ds : :class:`xarray.Dataset`
      an xarray dataset that may or may not 
      contain variables "lat_bnds" and 
      "lon_bnds"

    silent : bool, by default False
      if True, then suppresses standard output
    
    wrap_around_thresh : numeric or str, optional, default = ``'dynamic'``
      the minimum distance between the last grid cell 
      in longitude and the 'edges' of the coordinate 
      system for which the pixels are 'wrapped around'. 
      By default "dynamic", which sets this to twice
      the median difference between longitude values. 

      In either case, if the lowest and highest `lon`
      values have opposite sign, and are both within 
      `wrap_around_thresh` of -180 or 180, then the 
      grid is assumed to be "wrapping around", and pixels
      will line up around (or across, depending on the 
      grid) the anti-meridian.

    break_thresh_x : numeric, default = 2
      If the difference between consecutive coordinate
      values is `break_thresh_x` times the surrounding
      coordinate steps (determined by `break_window_width`), 
      then a break in the grid is assumed and the corresponding 
      grid cell is assumed to be as wide as the preceding 
      one, instead of reaching all the way across the "break". 
    
    break_window_width : numeric, default = 3
      Window used to determine how anomalous a difference
      between coordinate is when determining the location
      of breaks in the grid. 
          
    Returns
    ---------------
    ds : :class:`xarray.Dataset`
      the same dataset as inputted, unchanged if "lat/lon_bnds" 
      already existed, or with new variables "lat_bnds" and "lon_bnds"
      if not.
    """
    #----------- Setup -----------
    if (type(wrap_around_thresh) == str) and (wrap_around_thresh != 'dynamic'):
        raise ValueError('`wrap_around_thresh` must either be numeric or the string "dynamic"; instead, it is '+str(wrap_around_thresh)+'.')
    
    if ds.lon.max()>180:
        raise ValueError('Longitude seems to be in the 0:360 format.'+
                         ' -180:180 format required (e.g., run `xa.fix_ds(ds)` before inputting.')
        # Future versions should be able to work with 0:360 as well...
        # honestly, it *may* already work by just changing edges['lon']
        # to [0,360], but it's not tested yet. 
        
    if ('lat' not in ds) or ('lon' not in ds):
        raise KeyError('"lat"/"lon" not found in [ds]. Make sure the '+
                       'geographic dimensions follow this naming convention (e.g., run `xa.fix_ds(ds)` before inputting.')
    
    if ('lat_bnds' in ds) and ('lon_bnds' in ds):
        # `xa.fix_ds()` should rename bounds to `lat/lon_bnds`
        # If bounds present, do nothing
        return ds
    else:
        if not silent:
            print('lat/lon bounds not found in dataset; they will be created.')
        # Build lat / lon bound 
        for var in ['lat','lon']:
            # Get coordinate spacing
            diffs = ds[var].diff(var)

            if wrap_around_thresh == 'dynamic':
                wat_tmp = diffs.median().values*2
            else:
                wat_tmp = wrap_around_thresh
            
            # Get edge of grid coordinates 
            # (this does double duty by checking for antimeridian coordinates
            # since `xa.fix_ds()` sorts in ascending order by lon)
            edge_coords = ds[var][[0,-1]]

            #---------- Wrapping ----------
            # Figure out if grid wraps around / is continuous
            if var == 'lon':
                wrap_flag = (# 1) are they different sign
                             (np.prod(np.sign(edge_coords)) == -1) and
                             # 2) are they both within wrap_around_thresh of 180
                             np.all(np.abs(np.abs(edge_coords) - 180) < wat_tmp))
            else:
                wrap_flag = False
            
            if wrap_flag:
                # Wrap around lons
                ec = edge_coords.copy(deep=True).values
            
                # Change edge_lons to 0:360, to get difference
                # across antimeridian
                ec[ec<0] = 360+ec[ec<0]
            
                diffs = xr.concat([xr.DataArray(data = [np.abs(np.diff(ec))[0]],
                                     coords = {var:([var],[ds[var][0].values])}),
                                   diffs],
                                  dim=var)
            else:
                # Just do second diff as the first diff
                diffs = xr.concat([xr.DataArray(data = [diffs[0].values],
                                     coords = {var:([var],[ds[var][0].values])}),
                                   diffs],
                                  dim=var)
            
            #---------- Breaks ----------
            # Now, have to identify possible breaks in the coordinate system
            # i.e., [-179,-178,-177,178,179,180] would have a break between 
            # -177 and 178, which would create a pixel width of 355 degrees
            # Identify those breaks, and replace their pixel width with the 
            # preceding pixel width
            
            # Build a rolling filter that gets the average of surrounding 
            # coordinate steps for each coordinate (but ignores the coordinate's
            # step itself)
            break_check = xr.DataArray(np.ones(break_window_width),dims=['window'])
            break_check[int(np.floor(break_window_width/(break_window_width-1)))] = 0
            break_check = break_check/break_check.sum()
            
            # Figure out if any coordinate steps are more than > break_thresh_x
            # larger than surrounding coordinate steps
            breaks = (diffs / 
                      (diffs.rolling({var:break_window_width},center=True).
                       construct('window').dot(break_check))) > break_thresh_x
            
            # Make those coordinate steps the preceding coordinate step instead
            nbreaks = breaks.sum().values
            if (not silent) and nbreaks>0:
                print('Found '+str(nbreaks)+' break(s)/jump(s) in '+var+" (i.e., grid doesn't cover whole planet); replacing grid cell width(s) with preceding grid cell width(s) at breakpoint.")
            diffs[np.where(breaks)[0]] = diffs[np.where(breaks)[0]-1].values

            # Somewhat hack-y solution to deal with the fact that if only one lon 
            # value is found east of the anti-meridian, but the grid does wrap 
            # around, the break is not currently correctly captured
            if wrap_flag and ((ds[var]>0).sum(var) == 1):
                diffs[-1] = np.abs(np.diff(ec))[0]
            
            #---------- Create bounds ----------
            # Now, create bounds using those diffs
            bnds_tmp = xr.concat([ds[var]-0.5*diffs,
                                  ds[var]+0.5*diffs],dim='bnds').transpose(var,'bnds')
            
            #---------- Clean up ----------
            # And now, fix lingering issues
            # 1. lon bounds > 180 / < -180 (which occurs when the grid crosses the antimeridian
            #    in the step above)
            # 2. lat bounds not at the right edges
            # 3. making sure the lon bounds line up right around the antimeridian, if relevant
            
            # Change lon_bounds > 180 to lon_bounds > -180
            if var == 'lon':
                bnds_tmp = bnds_tmp.where(bnds_tmp<=180,-360 + bnds_tmp.where(bnds_tmp>180))
                bnds_tmp = bnds_tmp.where(bnds_tmp>=-180,360 + bnds_tmp.where(bnds_tmp<-180))
            
            # Curtail lat bands to -90, 90
            if var == 'lat':
                bnds_tmp = bnds_tmp.where(bnds_tmp<=90,90)
                bnds_tmp = bnds_tmp.where(bnds_tmp>=-90,-90)
            
            # If wrapping, make sure the bounds line up
            if wrap_flag: 
                if ((not bnds_tmp[0,0] == bnds_tmp[-1,-1]) and  
                    (not ((bnds_tmp[0,0] == -180) and (bnds_tmp[-1,-1] == 180)))):
                    bnds_tmp[-1,-1] = bnds_tmp[0,0]
            # Add to ds
            ds[var+'_bnds'] = bnds_tmp
            del bnds_tmp

    # This is the cf_xarray standard, to have them as coordinates
    # Currently only keeping them as variables instead (which is 
    # CMIP standard I think?)
    #ds = ds.set_coords([var+'_bnds' for var in ['lat','lon']])

     # Return
    return ds    


def subset_find(ds0,ds1,silent=False):
    """ Finds the grid of `ds1` in `ds0`, and subsets `ds0` to the grid in `ds1`
    
    Parameters
    ---------------
    ds0 : :class:`xarray.Dataset`
      an xarray Dataset to be subset based on the grid of ds1; must 
      contain grid variables "lat" or "lon" (could add a fix_ds call)
    ds1 : :class:`xarray.Dataset`, :class:`xarray.DataArray`
      either an `xarray` structrue (Dataset, DataArray) with "lat" "lon"
      variables, or a dictionary with DataArrays ['lat'] and ['lon'].
      IMPORTANT: `ds1` HAS TO BE BROADCAST - i.e. one value of lat, lon 
      each coordinate, with lat and lon vectors of equal length. This can 
      be done e.g. using ``ds1.stack(loc=('lat','lon'))``. 

    silent : bool, by default False
      if True, then suppresses standard output
        
    Returns
    ---------------
    ds0 : :class:`xarray.Dataset`
      The input `ds0`, subset to the locations in `ds1`. 
    
    """
    
    if 'loc' not in ds0.sizes:
        ds0 = ds0.stack(loc = ('lat','lon'))
        was_stacked = True
    else:
        was_stacked = False
    #if 'loc' not in ds1.sizes:
    #    ds1 = ds1.stack(loc = ('lat','lon'))
    
    # Need a test to make sure the grid is the same. So maybe the gdf_out class 
    # has the lat/lon grid included - and then we can skip the lat/lon column
    # and just keep the pix_idxs
    if (not np.allclose(len(ds0.lat),len(ds1['lat'])) or (not np.allclose(len(ds0.lon),len(ds1['lon']))) or
         (not np.allclose(ds1['lat'],ds0.lat)) or (not np.allclose(ds1['lon'],ds0.lon))):
        if not silent:
            print('adjusting grid... (this may happen because only a subset of pixels '+
              'were used for aggregation for efficiency - i.e. [subset_bbox=True] in '+
             'xa.pixel_overlaps())') #(this also happens because ds and ds_bnds above was already subset)
        # Zip up lat,lon pairs to allow comparison
        latlons = list(zip(ds0.lat.values,ds0.lon.values))
        latlons0 = list(zip(ds1['lat'].values,ds1['lon'].values))
        
        # Find indices of the used grid for aggregation in the input grid
        loc_idxs = [latlons.index(i) for i in latlons0]
        
        if np.allclose(len(loc_idxs),len(latlons0)):
            if not silent:
                print('grid adjustment successful')
            # Subset by those indices
            ds0 = ds0.isel(loc=loc_idxs)
        else:
            raise ValueError('Was not able to match grids!')
        
    if was_stacked:
        # THIS MAY NOT WORK IN ALL CASES
        ds0 = ds0.unstack()
        
    return ds0

