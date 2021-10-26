import xarray as xr
import numpy as np
import pandas as pd
import warnings


def normalize(a,drop_na = False):
    """ Normalizes the vector `a`

    The vector `a` is divided by its sum. 

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
    
    if (drop_na) & (np.any(np.isnan(a))):
        a2 = a[~np.isnan(a)]
        a2 = a2/a2.sum()
        a[~np.isnan(a)] = a2

        return a

    elif (np.all(~np.isnan(a))) & (a.sum()>0):
        return a/a.sum()
    else:
        return a*np.nan

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
        bnd_vars = [k for k in list(ds.keys()) if 'bnds' in k]
    elif type(ds) == xr.core.dataarray.DataArray:
        bnd_vars = []
    else:
        raise TypeError('[ds] needs to be an xarray structure (Dataset or DataArray).')
    
    # Fix lat/lon variable names (sizes instead of dims to be compatible with both ds, da...)
    if 'lat' not in ds.sizes.keys():
        test_dims = [k for k in var_cipher.keys() if k in ds.sizes.keys()]
        if len(test_dims) == 0:
            raise NameError('No valid lat/lon variables found in the dataset.')
        else:
            # shoudl there be a elif len()>1? If there are multiple things found?
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
        if ('lon_bnds' in ds.keys()):
            if (ds.lon_bnds.max()>180):
                ds['lon_bnds'] = (ds['lon_bnds'] + 180) % 360 - 180

    # Sort by lon; this should be robust (and necessary to avoid fix_ds(fix_ds(ds)) 
    # from failing)
    ds = ds.sortby(ds.lon)
    
    # Sort by latitude as well (to make latitude consistent)
    ds = ds.sortby(ds.lat)
        
    # Return fixed ds
    return ds


def get_bnds(ds,
             edges={'lat':[-90,90],'lon':[-180,180]},
             wrap_around_thresh=5):
    """ Builds vectors of lat/lon bounds if not present in `ds`
    
    Assumes a regular rectangular grid - so each lat/lon bound
    is 0.5*(gap between pixels) over to the next pixel.
    
    
    Parameters
    ---------------
    ds : :class:`xarray.Dataset`
      an xarray dataset that may or may not 
      contain variables "lat_bnds" and 
      "lon_bnds"

    wrap_around_thresh : numeric, optional, default = ``5`` (degrees)
      the minimum distance between the last 
      pixel edge and the 'edges' of the 
      coordinate system for which the pixels
      are 'wrapped around'. For example, given
      'lon' edges of [-180,180] and a 
      wrap_around_thresh of 5 (default), if 
      the calculated edges of pixels match
      the edge on one side, but not the other
      (i.e. -180 and 179.4) and this gap 
      (180-179.4) is less than 5, the -180 
      edge is changed to 179.4 to allow the pixel
      to 'wrap around' the edge of the coordinate
      system.
          
    Returns
    ---------------
    ds : :class:`xarray.Dataset`
      the same dataset as inputted, unchanged if "lat/lon_bnds" 
      already existed, or with new variables "lat_bnds" and "lon_bnds"
      if not.
    """
    if ds.lon.max()>180:
        raise ValueError('Longitude seems to be in the 0:360 format.'+
                         ' -180:180 format required.')
        # Future versions should be able to work with 0:360 as well...
        # honestly, it *may* already work by just changing edges['lon']
        # to [0,360], but it's not tested yet. 
        
    if ('lat' not in ds.keys()) | ('lon' not in ds.keys()):
        raise KeyError('"lat"/"lon" not found in [ds]. Make sure the '+
                       'geographic dimensions follow this naming convention.')
    
    
    if 'lat_bnds' in ds.keys():
        return ds
    else:
        print('lat/lon bounds not found in dataset; they will be created.')
        # Build lat / lon bound 
        for var in ['lat','lon']:
            bnds_tmp = xr.DataArray(data=np.zeros((ds.dims[var],2))*np.nan,
                                    dims=[var,'bnds'],
                                    coords=[ds[var],np.arange(0,2)])

            # Assign all non-edge bounds as just half of the distance from the center
            # of each pixel to the center of the next pixel
            bnds_tmp[1:,:] = xr.concat([ds[var]-0.5*ds[var].diff(var),
                                          ds[var]+0.5*ds[var].diff(var)],dim='bnd').transpose(var,'bnd')

            # Fill in last missing band before edge cases (the inner band of the 
            # first pixel, which is just equal to the next edge)
            bnds_tmp[0,1] = bnds_tmp[1,0]
            #print(bnds_tmp)
            # Now deal with edge cases; basically either use the diff from the last
            # interval between pixels, or max out at 90. 
            if ds[var].diff(var)[0]>0:
                bnds_tmp[0,0] = np.max([edges[var][0],ds[var][0].values-0.5*(ds[var][1]-ds[var][0]).values])
                bnds_tmp[-1,1] = np.min([edges[var][1],ds[var][-1].values+0.5*(ds[var][-1]-ds[var][-2]).values])
            else:
                bnds_tmp[0,0] = np.min([edges[var][1],ds[var][0].values+0.5*(ds[var][1]-ds[var][0]).values])
                bnds_tmp[-1,1] = np.max([edges[var][0],ds[var][-1].values-0.5*(ds[var][-1]-ds[var][-2]).values])
            #print(bnds_tmp)

            # Fix crossing-over-360 issues in the lon
            if var == 'lon':
                # To be robust to partial grids; setting the rolled over edges equal to
                # each other if one of the edges is the -180/180 and the other one isn't, 
                # but 'close enough' (within 5 degrees + warning)
                if (bnds_tmp[0,0] in edges[var]) & (bnds_tmp[-1,1] not in edges[var]):
                    # Make sure that the other edge is within the wrap-around threshold
                    # (to avoid wrapping around if a grid is only -180:45 for example)
                    if np.min(np.abs(bnds_tmp[-1,1].values-edges[var])) <= wrap_around_thresh:
                        if np.min(np.abs(bnds_tmp[-1,1].values-edges[var])) > ds[var].diff(var).max():
                            warnings.warn('Wrapping around '+[var]+' value of '+bnds_tmp[-1,1].values+', '+
                                          'because it is closer to a coordinate edge ('+
                                          ', '.join([str(n) for n in edges[var]])+') than the '+
                                          '[wrap_around_thresh] ('+str(wrap_around_thresh)+'); '+
                                          'however, it is farther away from that edge than the '+
                                          'maximum pixel width in the '+var+' direction. If this is '+
                                          'intended, no further action is necessary. Otherwise, reduce '+
                                          'the [wrap_around_thresh].')
                        bnds_tmp[0,0] = bnds_tmp[-1,1]
                elif (bnds_tmp[0,0] not in edges[var]) & (bnds_tmp[-1,1] in edges[var]):
                    if np.min(np.abs(bnds_tmp[0,0].values-edges[var])) <= wrap_around_thresh:
                        if np.min(np.abs(bnds_tmp[0,0].values-edges[var])) > ds[var].diff(var).max():
                            warnings.warn('Wrapping around '+[var]+' value of '+bnds_tmp[0,0].values+', '+
                                          'because it is closer to a coordinate edge ('+
                                          ', '.join([str(n) for n in edges[var]])+') than the '+
                                          '[wrap_around_thresh] ('+str(wrap_around_thresh)+'); '+
                                          'however, it is farther away from that edge than the '+
                                          'maximum pixel width in the '+var+' direction. If this is '+
                                          'intended, no further action is necessary. Otherwise, adjust '+
                                          'the [wrap_around_thresh].')
                    bnds_tmp[-1,1] = bnds_tmp[0,0]
            # Add to ds
            ds[var+'_bnds'] = bnds_tmp
            del bnds_tmp
        
    # Return
    return ds    



def subset_find(ds0,ds1):
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
        
    Returns
    ---------------
    ds0 : :class:`xarray.Dataset`
      The input `ds0`, subset to the locations in `ds1`. 
    
    """
    
    if 'loc' not in ds0.sizes.keys():
        ds0 = ds0.stack(loc = ('lat','lon'))
        was_stacked = True
    else:
        was_stacked = False
    #if 'loc' not in ds1.sizes.keys():
    #    ds1 = ds1.stack(loc = ('lat','lon'))
    
    # Need a test to make sure the grid is the same. So maybe the gdf_out class 
    # has the lat/lon grid included - and then we can skip the lat/lon column
    # and just keep the pix_idxs
    if (not np.allclose(len(ds0.lat),len(ds1['lat'])) or (not np.allclose(len(ds0.lon),len(ds1['lon']))) or
         (not np.allclose(ds1['lat'],ds0.lat)) or (not np.allclose(ds1['lon'],ds0.lon))):
        print('adjusting grid... (this may happen because only a subset of pixels '+
              'were used for aggregation for efficiency - i.e. [subset_bbox=True] in '+
             'xa.pixel_overlaps())') #(this also happens because ds and ds_bnds above was already subset)
        # Zip up lat,lon pairs to allow comparison
        latlons = list(zip(ds0.lat.values,ds0.lon.values))
        latlons0 = list(zip(ds1['lat'].values,ds1['lon'].values))
        
        # Find indices of the used grid for aggregation in the input grid
        loc_idxs = [latlons.index(i) for i in latlons0]
        
        if np.allclose(len(loc_idxs),len(latlons0)):
            print('grid adjustment successful')
            # Subset by those indices
            ds0 = ds0.isel(loc=loc_idxs)
        else:
            raise ValueError('Was not able to match grids!')
        
    if was_stacked:
        # THIS MAY NOT WORK IN ALL CASES
        ds0 = ds0.unstack()
        
    return ds0

