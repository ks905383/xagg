import xarray as xr
import pandas as pd
import geopandas as gpd
import numpy as np
import os
import re
import warnings
import shutil

from . auxfuncs import (fix_ds,get_bnds,subset_find)


def export_weightmap(wm_obj,fn,overwrite=False):
    """ Save a copy of the weightmap, to avoid recalculating it
    """
    ###### NEED TO SETUP HD5 DEPENDENCY? THROUGH OPTIONAL DEPENDENCY MAYBE? 
    warnings.warn('export_weightmap() is still an experimental feature. use with care.')


    if (not overwrite) and (os.path.exists(fn)):
        raise FileExistsError(fn+'/ already exists. Either change the [fn] or set overwrite=True.')

    try: 
        ###### Save geometry
        # This also creates a folder, into which we can add the rest
        # of the files
        wm_obj.geometry.to_file(fn)

        ####### Save agg 
        # Turn into a dataframe from a geodataframe; only dataframes
        # can be converted to HD5 files (and there's no geographic 
        # information here that's needed; can all be saved in the 
        # geographic save below)
        df_out = pd.DataFrame(wm_obj.agg)
        with warnings.catch_warnings():
            # Catch performance warning about pickling
            warnings.filterwarnings('ignore')
            df_out.to_hdf(fn+'/'+re.split('/',fn)[-1]+'.h5','wm')

        ###### Save source grid
        for k in wm_obj.source_grid:
            # Adding v to variable (i.e., 'latv', 'lonv') needed 
            # because of the multi-index being used here 
            source_grid_tmp = wm_obj.source_grid[k].reset_index('loc').to_dataset(name=k+'v')

            # Replicate dataset, to avoid xarray bug stemming from 
            # Pandas MultiIndex testing... 
            source_grid_tmp = source_grid_tmp.reset_coords()
            source_grid_tmp = xr.Dataset({k:([d for d in v.sizes],v.values) for k,v in source_grid_tmp.items()})
            source_grid_tmp = source_grid_tmp.set_coords([d for d in source_grid_tmp if d not in ['loc', k+'v']])

            # Save
            source_grid_tmp.to_netcdf(fn+'/'+re.split('/',fn)[-1]+'_'+k+'.nc')

        ###### Save weights
        # (If weights == 'nowghts', then no file to load)
        if (type(wm_obj.weights) is not str) or (wm_obj.weights != 'nowghts'):
            # Setting astype(object) to make sure integral weights
            # don't change the general type of the frame. This may
            # only affect the testing routines, but setting this here
            # to be explicit 
            wm_obj.weights.astype(object).to_csv(fn+'/'+re.split('/',fn)[-1]+'_weights.csv')
    except RuntimeError as error:
        print(error)
        # Remove files that have already been generated, to make 
        # sure no flawed files are floating around. 
        shutil.rmtree(fn)
        raise Exception('Error while saving. The directory '+fn+' has been removed.')


def prep_for_nc(agg_obj,loc_dim='poly_idx'):
    """ Preps aggregated data for output as a netcdf
    
    Concretely, aggregated data is placed in a new xarray dataset 
    with dimensions of location (the different polygons in gdf_out)
    and any other dimension(s) in the original input raster data. 
    All fields from the input polygons are kept as variables with 
    dimension of location.
    
    Parameters
    ---------------
    agg_obj : :class:`xagg.classes.aggregated`

    loc_dim : :py:class:`str`, optional
        the name of the location dimension; by definition 
        'poly_idx'. Values of that dimension are currently
        only an integer index (with further information given
        by the field variables). Future versions may allow, 
        if loc_dim is set to the name of a field in the input
        polygons, to replace the dimension with the values of
        that field. 
               
    Returns:
    ---------------
    ds: :class:`xarray.Dataset`
        an :class:`xarray.Dataset` containing the aggregated variables in addition
        to the original fields from the location polygons. Dimensions are
        a location dimension (counting down the polygons - this is the 
        dimension of all the field contents) and any other non-location
        dimensions contained in the variables before being aggregated
    """
    # To output as netcdf, first put data back into an xarray dataset
    
    # Create xarray dataset with the aggregation polygons (poly_idx) 
    # there. 
    ds_out = xr.Dataset(coords={'poly_idx':(['poly_idx'],agg_obj.agg.poly_idx.values)})
    
    # Add other polygon attributes
    for var in [c for c in agg_obj.agg.columns if c not in ['poly_idx','rel_area','pix_idxs','coords']]:
        if var not in agg_obj.ds_in.var():
            # For auxiliary variables (from the shapefile), just copy them wholesale into the dataset
            ds_out[var] = xr.DataArray(data=agg_obj.agg[var],coords=[agg_obj.agg.poly_idx],dims=['poly_idx'])
        else:
            # For data variables (from the input grid), create empty array
            ds_out[var] = xr.DataArray(data=np.zeros((len(agg_obj.agg),
                                                         *[agg_obj.ds_in[var].sizes[k] for k in agg_obj.ds_in[var].sizes.keys() if k not in ['lat','lon','loc']]))*np.nan,
                                       dims=['poly_idx',*[k for k in agg_obj.ds_in[var].sizes.keys() if k not in ['lat','lon','loc']]],
                                       coords=[[k for k in agg_obj.agg.poly_idx],*[agg_obj.ds_in[var][k].values for k in agg_obj.ds_in[var].sizes.keys() if k not in ['lat','lon','loc']]])
        
            # Now insert aggregated values 
            for poly_idx in agg_obj.agg.poly_idx:
                ds_out[var].loc[{'poly_idx':poly_idx}] = np.squeeze(agg_obj.agg.loc[poly_idx,var])
    
    # Add non-geographic coordinates for the variables to be aggregated
    for crd in [k for k in agg_obj.ds_in.sizes.keys() if (k not in ['lat','lon','loc','bnds'])]:
        ds_out[crd] = xr.DataArray(dims=[crd],
                                   data=agg_obj.ds_in[crd].values,
                                   coords=[agg_obj.ds_in[crd].values])
        
    # Remove "name" variable created in get_pixel_overlaps
    if 'name' in ds_out:
        if ((ds_out['name'].dims == ('poly_idx',)) and 
            (len(np.unique(ds_out['name'])) == 1) and
            np.unique(ds_out['name'])[0] in ds_out):
            ds_out = ds_out.drop_vars(['name'])

    # Rename poly_idx if desired
    if loc_dim != 'poly_idx':
        ds_out = ds_out.rename({'poly_idx':loc_dim})
        
    # Return ds_out
    return ds_out


def prep_for_csv(agg_obj,add_geom=False):
    """ Preps aggregated data for output as a csv
    
    Concretely, aggregated data is placed in a new pandas dataframe
    and expanded **wide** - each aggregated variable is placed in new 
    columns; one column per coordinate in each dimension that isn't
    the location (poolygon). So, for example, a lat x lon x time
    variable "tas", aggregated to location x time, would be reshaped 
    long to columns "tas0", "tas1", "tas2",... for timestep 0, 1, etc.

    For data **long**, use ``agg_obj.to_dataframe().to_csv()`` instead. 
    
    Note: 
    Currently no support for variables with more than one extra dimension
    beyond their location dimensions. Potential options: a multi-index
    column name, so [var]0-0, [var]0-1, etc...
    
    Parameters
    ----------------
    agg_obj : :class:`xagg.classes.aggregated`
        the output from :func:`aggregate`

    add_geom : bool, default `False`
        if `True`, adds the geometry from the original shapefile/
        geodataframe back into the prepped output; this is for 
        the `.to_geodataframe()` conversion option. 

               
    Returns
    ----------------
    df : :class:`pandas.DataFrame`
        a pandas dataframe containing all the fields from the original 
        location polygons + columns containing the values of the aggregated
        variables at each location. This can then easily be exported as a 
        csv directly (using ``df.to_csv()``) or to shapefiles by first turning into
        a geodataframe. 
    
    """
    # Test to make sure there's only one non-location dimension
    var_dims = {var:[d for d in agg_obj.ds_in[var].sizes if d != 'loc'] for var in agg_obj.ds_in}
    var_ndims = {var:len(dims) for var,dims in var_dims.items()}
    if np.max([n for d,n in var_ndims.items()]) > 1:
        raise NotImplementedError('The `agg` object has variables with more than 1 non-location dimension; '+
                                  'agg.to_csv() and agg.to_geodataframe() return wide arrays, but the code can not yet create wide arrays spanning data with multiple `wide` dimensions. '+
                                  'Try agg.to_dataframe() instead. (the offending variables with their non-location dimensions are '+
                                  str({var:dims for var,dims in var_dims.items() if var_ndims[var]>1})+')')

    # For output into csv, work with existing geopandas data frame
    csv_out = agg_obj.agg.drop(columns=['rel_area','pix_idxs','coords','poly_idx'])
    
    # Now expand the aggregated variable into multiple columns
    for var in [c for c in agg_obj.agg.columns if ((c not in ['poly_idx','rel_area','pix_idxs','coords']) & (c in agg_obj.ds_in.var()))]:
        # NOT YET IMPLEMENTED: dynamic column naming - so if it recognizes 
        # it as a date, then instead of doing var0, var1, var2,... it does
        # varYYYYMMDD etc.
        # These are the coordinates of the variable in the original raster
        dimsteps = [agg_obj.ds_in[var][d].values for d in agg_obj.ds_in[var].sizes.keys() if d not in ['lat','lon','loc']]
        # ALSO SHOULD check to see if the variables are multi-D - if they are
        # there are two options: 
        # - a multi-index column title (var0-0, var0-1)
        # - or an error saying csv output is not supported for this
    
        # Reshape the variable wide and name the columns [var]0, [var]1,...
        if len(dimsteps) == 0:
            # (in this case, all it does is move from one list per row to 
            # one value per row)
            #expanded_var = (pd.DataFrame(pd.DataFrame(csv_out[var].to_list())[0].to_list(),
            #                         columns=[var]))
            expanded_var = pd.DataFrame(csv_out[var].apply(np.squeeze).to_list(),
                                        columns=[var])
        else:
            #expanded_var = (pd.DataFrame(pd.DataFrame(csv_out[var].to_list())[0].to_list(),
            #                         columns=[var+str(idx) for idx in np.arange(0,len(csv_out[var][0][0]))]))
            expanded_var = pd.DataFrame(csv_out[var].apply(np.squeeze).to_list(),
                                        columns = [var+str(idx) for idx in range(len(csv_out[var].apply(np.squeeze))+1)])
        # Append to existing series
        csv_out = pd.concat([csv_out.drop(columns=(var)),
                             expanded_var],
                            axis=1)
        del expanded_var

    if add_geom:
        # Return the geometry from the original geopandas.GeoDataFrame
        csv_out['geometry'] = agg_obj.geometry
        csv_out = csv_out.set_geometry('geometry')
        
    # Return 
    return csv_out


def output_data(agg_obj,output_format,output_fn,loc_dim='poly_idx',silent=False):
    """ Wrapper for `prep_for_*` functions
    
    
    Parameters
    ---------------
    agg_obj : :class:`xagg.classes.aggregated` 
        object to be exported   

    output_format : :py:class:`str`
        'netcdf', 'csv', or 'shp'
    output_fn: :py:class:`str`
        the output filename   
    loc_dim : :py:class:`str`, optional. default = ``'poly_idx'``
        the name of the
        dimension with location indices; used
        by :func:`xagg.export.prep_for_nc` for nc, 
        csv output
                     
    Returns
    ---------------
    the variable that gets saved, so depending on the `output_format`:
        -  "netcdf": the :class:`xarray.Dataset` on which ``.to_netcdf``
           was called
        -  "csv": the :class:`pandas.Dataframe` on which ``.to_csv`` 
           was called (uses the `xarray` `ds.to_dataframe()` functionality)
        - "shp": the :class:`geopandas.GeoDataDrame` on which ``.to_file`` was called
     
    """

    if output_format == 'netcdf':

        ds_out = prep_for_nc(agg_obj,loc_dim=loc_dim)

        # Save as netcdf
        if not output_fn.endswith('.nc'):
            output_fn = output_fn+'.nc'
        ds_out.to_netcdf(output_fn)
        if not silent:
            print(output_fn+' saved!')

        # Return
        return ds_out

    elif output_format == 'csv':

        csv_out = prep_for_nc(agg_obj,loc_dim=loc_dim)
        csv_out = csv_out.to_dataframe()

        # Save as csv
        if not output_fn.endswith('.csv'):
            output_fn = output_fn+'.csv'
        csv_out.to_csv(output_fn)
        if not silent:
            print(output_fn+' saved!')

        # Return 
        return csv_out

    elif output_format == 'shp':
        # This uses the same processing as the to_csv option above, since 
        # shapefiles function similarly (each field can only have one value
        # for each polygon, etc.)
        csv_out = prep_for_csv(agg_obj)

        # gdf_in.geometry should be replaced with gdf_out['geometry']
        # which should be kept when gdf_out is created... 
        # Turn back into GeoDataFrame
        shp_out = gpd.GeoDataFrame(csv_out,geometry=agg_obj.geometry)

        # Export that GeoDataFrame to shapefile
        if not output_fn.endswith('.shp'):
            output_fn = output_fn+'.shp'
        shp_out.to_file(output_fn)
        if not silent:
            print(output_fn+' saved!')

        # Return
        return shp_out

    else: 
        raise KeyError(output_format+' is not a supported output format.')