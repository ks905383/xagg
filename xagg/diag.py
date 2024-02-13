
import xarray as xr
import numpy as np 
import geopandas as gpd
import warnings
from cartopy import crs as ccrs
from matplotlib import pyplot as plt
import matplotlib as mpl
import shapely

from . core import create_raster_polygons


def diag_fig(wm,poly_id,pix_overlap_info,
	 cmap='magma_r',
	  max_title_depth = 5,
	   fig=None,ax=None):
	""" Create a diagnostic figure showing overlap between pixels and a given polygon

	Parameters
    ---------------
    wm : :class:`xagg.classes.weightmap`
        the output to :func:`xagg.core.get_pixel_overlaps`

    poly_id : int, list, or dict
      	which polygon to plot. If `int`, then just the polygon in that indexed row of `wm.agg`
      	is plotted. If `list` of `int`s, then those polygons are plotted. If dict, then all
      	matches in the geodataframe are plotted (e.g., `poly_id = {'adm1_code':'URY-8'}`). 

    pix_overlap_info : `xarray.classes.dataset.Dataset`, `xarray.classes.dataarray.DataArray`, 
    				   or the output to `xa.core.create_raster_polygons`
    	if `ds` or `da`: the original dataset used to calculate the `wm`; needs to be re-entered 
    					 here because `wm` does not keep the original pixel polygons 
    	Otherwise, put in the output to `xa.core.create_raster_polygons(ds)`

    cmap : str, by default 'magma_r'
    	colormap, must be the name of a matplotlib-recognized colormap

	max_title_depth : int, by default 5
		if only showing one polygon, then the plot title is `', '.join()`, called on the 
		first `max_title_depth` columns of `wm.agg` that aren't added by `xagg`

	fig : `mpl.figure.Figure` or None, by default None
		if not None, then this figure handle is used

	ax : `mpl.axis.Axis` or None, by default None
		if not None, then this axis handle is used
    
    Returns
    ---------------
    fig,ax : the matplotlib figure, axis handles

	"""

	#---------- Setup ----------
	if type(poly_id) == int:
	    poly_idx = [poly_id]
	elif type(poly_id) == dict:
	    poly_idx = list((wm.agg.loc[np.prod(np.array([wm.agg[k] == v for k,v in poly_id.items()]),axis=0).astype(bool),:].
	                index.values))
	elif type(poly_id) == list:
	    if not np.all([type(k) == int for k in poly_id]):
	        raise TypeError('If using list polygon ids, all list members must be integers corresponding to polygon idxs in `wm.agg`.')
	    poly_idx = poly_id

	# Get pixel polygons/overlaps, if necessary
	if ((type(pix_overlap_info) == xr.core.dataset.Dataset)
	     or (type(pix_overlap_info) == xr.core.dataarray.DataArray)):
	    pix_polys = create_raster_polygons(pix_overlap_info)
	else:
	    pix_polys = pix_overlap_info

	# Get which polygon to plot
	plot_poly = [wm.geometry[idx] for idx in poly_idx]
	plot_overlap = [wm.agg.iloc[idx,:] for idx in poly_idx]

	# Get pixels overlapping with the desired polygon
	pix_polys = [pix_polys['gdf_pixels'].loc[overlaps['pix_idxs']]
	             for overlaps in plot_overlap]

	# Get colors for relative overlap
	for idx in np.arange(0,len(pix_polys)):
	    pix_polys[idx]['rel_area'] = plot_overlap[idx]['rel_area'][0].values
	    pix_polys[idx]['color'] = [tuple(x) for x in plt.get_cmap(cmap)(pix_polys[idx]['rel_area'] / 
	                                                               pix_polys[idx]['rel_area'].max())]


	#---------- Plot ----------
	if fig == None:
	    fig = plt.figure()

	if ax == None:
	    # Create sample figure, centered on polygon centroid
	    ax = plt.subplot(projection=ccrs.PlateCarree(central_longitude=[k for k in 
	                                                                    gpd.GeoDataFrame(geometry=[poly.centroid for poly in plot_poly]).dissolve().centroid.values[0].coords][0][0]))

	# Plot polygon
	for idx in np.arange(0,len(plot_poly)):
	    if type(plot_poly[idx]) == shapely.geometry.polygon.Polygon:
	        plt.plot(*plot_poly[idx].exterior.xy,transform=ccrs.PlateCarree(),
	                  color='tab:green',linewidth=1)
	    else:
	         [plt.plot(*k.exterior.xy,transform=ccrs.PlateCarree(),
	                  color='tab:green',linewidth=1) for k in plot_poly[idx].geoms] #.geoms

	# Add reference coastlines
	ax.coastlines(color='grey')

	# Add pixels, transparent, but colored by overlap
	for idx in np.arange(0,len(plot_poly)):
	    pix_polys[idx].plot(color=pix_polys[idx].color,transform=ccrs.PlateCarree(),ax=ax,
	    					alpha=0.8)

	#------- Annotate -------
	# Title
	if len(plot_overlap) == 1:
	    ax.set_title('Poly #'+str(poly_idx[0])+': '+
	             '; '.join([str(plot_overlap[0][k]) for k in plot_overlap[0].index 
	                        if k not in ['poly_idx','rel_area','pix_idxs','coords']][0:max_title_depth]))
	else:
	    ax.set_title('Poly #s: '+', '.join([str(k) for k in poly_idx]))

	# Colorbar
	fig.subplots_adjust(right=0.825)
	cax = fig.add_axes([0.875, 0.15, 0.025, 0.7])
	if len(plot_overlap) == 1:
	    clabel = 'Overlap area (relative to largest pixel overlap)'
	else:
	    clabel = ('Overlap area (relative to largest pixel overlap,\n'+
	              'recalculated for each polygon; so they may not be\ndirectly comparable'+
	              ' across polygons)')
	    
	cb1 = mpl.colorbar.ColorbarBase(ax=cax, cmap=cmap,
	                                norm=mpl.colors.Normalize(vmin=0, vmax=1),
	                                orientation='vertical',
	                                label = clabel)


	# Gridlines
	gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=['x','y','bottom','left'],
	                  linewidth=1, color='gray', alpha=0.5, linestyle=':')
