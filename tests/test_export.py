import pytest
import pandas as pd
import numpy as np
import xarray as xr
import shutil
import geopandas as gpd
from geopandas import testing as gpdt
from unittest import TestCase
from shapely.geometry import Polygon

from xagg.core import (process_weights,create_raster_polygons,get_pixel_overlaps,aggregate,read_wm)
from xagg.wrappers import (pixel_overlaps)


##### to_dataset() tests #####
# build raster polygons from a simple 2x2 grid of lat/lon pixels
ds = xr.Dataset({'test':(['lon','lat','run'],np.array([[[0,1],[2,3]],[[0,1],[2,3]]])),
				 'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5]])),
				 'lon_bnds':(['lon','bnds'],np.array([[-0.5,0.5],[0.5,1.5]]))},
				coords={'lat':(['lat'],np.array([0,1])),
						'lon':(['lon'],np.array([0,1])),
						'run':(['run'],np.array([0,1])),
						'bnds':(['bnds'],np.array([0,1]))})

# Create polygon covering multiple pixels
gdf = {'name':['test'],
			'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")

# Get pixel overlaps
wm = pixel_overlaps(ds,gdf)

# Get aggregate
agg = aggregate(ds,wm)

def test_to_dataset(agg=agg):
	# Change to dataset
	ds_out = agg.to_dataset()

	# Build reference output dataset
	ds_ref = xr.Dataset({'name':(['poly_idx'],np.array(['test']).astype(object)),
						 'test':(['poly_idx','run'],np.array([[1.0,2.0]]))},
					coords={'poly_idx':(['poly_idx'],np.array([0])),
							'run':(['run'],np.array([0,1]))})

	# Assert equal within tolerance, again likely due to very slight 
	# variation off from actual 1.0, 2.0 due to crs
	xr.testing.assert_allclose(ds_out,ds_ref,atol=0.0001)

def test_to_dataframe(agg=agg):
	# Change to dataframe
	df_out = agg.to_dataframe()

	# Build reference output dataframe
	df_ref = pd.DataFrame({'poly_idx':[0,0],'run':[0,1],'name':['test','test'],'test':[0.9999,1.9999]})
	df_ref = df_ref.set_index(['poly_idx','run'])

	# Assert equal within tolerance, again likely due to very slight 
	# variation off from actual 1.0, 2.0 due to crs
	pd.testing.assert_frame_equal(df_out,df_ref,atol=0.0001)


##### pixel_overlaps() export tests #####

def test_pixel_overlaps_export_and_import(ds=ds):
	# Testing the .to_file() --> read_wm() workflow. 
	# Rather complex because of the many different components of wm. 
	# wm.agg in particular is a dataframe with lists in it (which is 
	# of course not great) which makes it hard to compare. This should
	# work. The main issue is this converts agg to a dataframe from a 
	# geodataframe; but I don't think that changes anything (it doesn't
	# have any geographic information anyways) 

	fn = 'wm_export_test'

	# Add a simple weights grid
	weights = xr.DataArray(data=np.array([[0.,1.],[2.,3.]]),
								dims=['lat','lon'],
								coords=[ds.lat,ds.lon])

	# Create polygon covering one pixel
	gdf = {'name':['test'],
				'geometry':[Polygon([(-0.5,-0.5),(-0.5,0.5),(0.5,0.5),(0.5,-0.5),(-0.5,-0.5)])]}
	gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")

	# Calculate weightmap
	wm_out = pixel_overlaps(ds,gdf,weights=weights)

	# Export weightmap
	wm_out.to_file(fn,overwrite=True)

	# Load weightmap
	wm_in = read_wm(fn)


	###### agg
	# This is just checking the parts of the thing that aren't lists. Unfortunately putting lists into
	# dataframes is a bad idea, but that's the best method I've come up with so far. 
	pd.testing.assert_frame_equal(wm_out.agg[[v for v in wm_out.agg if v not in ['rel_area','pix_idxs','coords']]],
								  wm_in.agg[[v for v in wm_in.agg if v not in ['rel_area','pix_idxs','coords']]])
	# Now for the columns with lists
	for v in ['rel_area','pix_idxs','coords']:
		for it in np.arange(0,wm_out.agg[v].shape[0]):
			np.testing.assert_array_equal(wm_out.agg[v].values[it],wm_in.agg[v].values[it])
	# Now for the column names to make sure all the columns are remained
	pd.testing.assert_index_equal(wm_out.agg.columns,wm_in.agg.columns)

	###### geometry
	gpdt.assert_geoseries_equal(wm_out.geometry,wm_in.geometry)

	###### source grids
	for k in wm_out.source_grid:
		xr.testing.assert_allclose(wm_out.source_grid[k],wm_in.source_grid[k])

	###### weights
	if (type(wm_out.weights) is str) and (wm_out.weights=='nowghts'):
		np.testing.assert_string_equal(wm_in.weights,wm_out.weights)
	else:
		# `read_wm()` reads in weights as objects (see notes in relevant
		# section of `read_wm()`... this shouldn't have too big 
		# of a consequence, but does make this test more complicated
		pd.testing.assert_series_equal(wm_in.weights,wm_out.weights.astype(object))

	##### clean 
	shutil.rmtree(fn)

