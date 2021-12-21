from .export import (prep_for_nc,prep_for_csv,output_data)
import warnings

try:
    import cartopy 
    import cartopy.crs as ccrs
    import cmocean
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    no_plotting = False
except ImportError:
    no_plotting = True

# POSSIBLE CHANGE: I'm not quite sure how python deals with memory
# in this case - but [aggregated], which contains [ds] as one of its
# fields (for keeping the right variables / shapes upon export) may
# result in a bloated variable (especially if [ds] is particularly
# large). A fix could be to instead of [ds] replace it with just 
# what it needs; which is the sizes of all the variables, and the
# coordinates. However, in that case, prep_for_nc() and prep_for_csv()
# would have to be modified.

class weightmap(object):
    """ Class for mapping from pixels to polgyons, output from :func:`xagg.wrappers.pixel_overlaps`
    
    """
    def __init__(self,agg,source_grid,geometry,weights='nowghts'):
        self.agg = agg
        self.source_grid = source_grid
        self.geometry = geometry
        self.weights = weights
        
    def diag_fig(self,poly_idx):
        """ (NOT YET IMPLEMENTED) a diagnostic figure of the aggregation
        """
        raise NotImplementedError('diagnostic figure not yet implemented.')
        # Eventually, will print a diagnostic figure showing 
        # - the polygon defined by "poly_idx"
        # - the pixels overlapping polygon
        #    - shaded by the amount of overlapping area 
        #    - or the total weight (area+weight)


class aggregated(object):
    """ Class for aggregated data, output from :func:`xagg.core.aggregate`
    """
    def __init__(self,agg,source_grid,geometry,ds_in,weights='nowghts'):
        self.agg = agg
        self.source_grid = source_grid
        self.geometry = geometry
        self.ds_in = ds_in
        self.weights = weights
    
    # Conversion functions
    def to_dataset(self,loc_dim='poly_idx'):
        """ Convert to xarray dataset.
        """
        ds_out = prep_for_nc(self,loc_dim=loc_dim)
        return ds_out

    def to_geodataframe(self):
        """ Convert to geopandas geodataframe.
        """
        df_out = prep_for_csv(self,add_geom=True)
        return df_out

    def to_dataframe(self,loc_dim='poly_idx'):
        """ Convert to pandas dataframe.
        """
        df_out = prep_for_nc(self,loc_dim=loc_dim)
        df_out = df_out.to_dataframe()
        return df_out

    # Export functions
    def to_netcdf(self,fn,loc_dim='poly_idx'):
        """ Save as netcdf
        """
        output_data(self,
                   output_format = 'netcdf',
                   output_fn = fn,
                   loc_dim = loc_dim)
        
    def to_csv(self,fn):
        """ Save as csv
        """
        output_data(self,
                   output_format = 'csv',
                   output_fn = fn)
        
    def to_shp(self,fn):
        """ Save as shapefile
        """
        output_data(self,
                   output_format = 'shp',
                   output_fn = fn)



