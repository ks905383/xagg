import pytest
import pandas as pd
import numpy as np
import xarray as xr
import shutil
import os
import geopandas as gpd
from geopandas import testing as gpdt
from unittest import TestCase
from unittest.mock import patch
from io import StringIO
from shapely.geometry import Polygon

from xagg.core import (process_weights,create_raster_polygons,get_pixel_overlaps,aggregate,read_wm)
from xagg.wrappers import pixel_overlaps
from xagg.export import prep_for_csv
from xagg.options import (set_options,get_options)


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
	ds_ref = xr.Dataset({'test':(['poly_idx','run'],np.array([[1.0,2.0]]))},
					coords={'poly_idx':(['poly_idx'],np.array([0])),
							'run':(['run'],np.array([0,1]))})

	# Assert equal within tolerance, again likely due to very slight 
	# variation off from actual 1.0, 2.0 due to crs
	xr.testing.assert_allclose(ds_out,ds_ref,atol=0.0001)

def test_to_dataset_renamelocdim(agg=agg):
	# Change to dataset, renaming `poly_idx`
	ds_out = agg.to_dataset(loc_dim = 'county_id')

	# Build reference output dataset
	ds_ref = xr.Dataset({'test':(['county_id','run'],np.array([[1.0,2.0]]))},
					coords={'county_id':(['county_id'],np.array([0])),
							'run':(['run'],np.array([0,1]))})

	# Assert equal within tolerance, again likely due to very slight 
	# variation off from actual 1.0, 2.0 due to crs
	xr.testing.assert_allclose(ds_out,ds_ref,atol=0.0001)

##### to_dataframe() tests #####
def test_to_dataframe(agg=agg):
	# Change to dataframe
	df_out = agg.to_dataframe()

	# Build reference output dataframe
	df_ref = pd.DataFrame({'poly_idx':[0,0],'run':[0,1],'test':[0.9999,1.9999]})
	df_ref = df_ref.set_index(['poly_idx','run'])

	# Assert equal within tolerance, again likely due to very slight 
	# variation off from actual 1.0, 2.0 due to crs
	pd.testing.assert_frame_equal(df_out,df_ref,atol=0.0001)

def test_to_dataframe_renamelocdim(agg=agg):
	# Change to dataframe, with a new name for the location dimension
	df_out = agg.to_dataframe(loc_dim='county_id')

	# Build reference output dataframe
	df_ref = pd.DataFrame({'county_id':[0,0],'run':[0,1],'test':[0.9999,1.9999]})
	df_ref = df_ref.set_index(['county_id','run'])

	# Assert equal within tolerance
	pd.testing.assert_frame_equal(df_out,df_ref,atol=0.0001)

##### to_csv() tests #####
def test_to_csv(agg=agg):
	# Change to dataframe
	df_out = agg.to_csv('test.csv')

	try:
		# Build reference output dataframe
		df_ref = pd.DataFrame({'poly_idx':[0,0],'run':[0,1],'test':[0.9999,1.9999]})
	
		# Load
		df_in = pd.read_csv('test.csv')

		# Assert equal within tolerance
		pd.testing.assert_frame_equal(df_in,df_ref,atol=0.0001)
	
	finally:
		# Clean
		os.remove('test.csv')

##### to_geodataframe() tests #####
def test_to_geodataframe(agg=agg):
	# Change to geodatafarme
	df_out = agg.to_geodataframe()

	# Build reference output geodataframe
	gdf_ref = gpd.GeoDataFrame({'name':['test'],
	                  				'test0':[0.9999629411369734],
	                  				'test1':[1.9999629411369735]},
	                 				geometry = [Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])],
	                 				crs = 'EPSG:4326')

	# Assert equal within tolerance, again likely due to very slight 
	# variation off from actual 1.0, 2.0 due to crs
	gpdt.assert_geodataframe_equal(df_out,gdf_ref,check_less_precise=True)

def test_prep_for_csv_multd():
	# Test to make sure .prep_for_csv() (for .to_geodataframe()
	# and .to_csv()) fails if you have a variable with more than one
	# non-location dimension
	
	# Have a 4-D variable (run and time in addition to geographic data)
	ds_extrad = xr.Dataset({'test':(['lon','lat','run','time'],np.random.rand(2,2,2,5)),
					 'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5]])),
					 'lon_bnds':(['lon','bnds'],np.array([[-0.5,0.5],[0.5,1.5]]))},
					coords={'lat':(['lat'],np.array([0,1])),
							'lon':(['lon'],np.array([0,1])),
							'run':(['run'],np.array([0,1])),
	                        'time':(['time'],pd.date_range('2001-01-01','2001-01-05')),
							'bnds':(['bnds'],np.array([0,1]))})

	# Create polygon covering multiple pixels
	gdf = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")

	# Get pixel overlaps
	wm_extrad = pixel_overlaps(ds_extrad,gdf)

	# Get aggregate
	agg_extrad = aggregate(ds_extrad,wm_extrad)

	with pytest.raises(NotImplementedError):
		# Make sure 
		prep_for_csv(agg_extrad)

##### to_shp() tests #####
def test_to_shp(agg=agg):
	# Export to shapefile
	agg.to_shp('test.shp')

	try:
		# Make reference geodataframe
		gdf_ref = gpd.GeoDataFrame({'name':['test'],
	                  				'test0':[0.9999629411369734],
	                  				'test1':[1.9999629411369735]},
	                 				geometry = [Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])],
	                 				crs = 'EPSG:4326')

		# Read file
		gdf_in = gpd.read_file('test.shp')

		# Assert equal within tolerance
		gpdt.assert_geodataframe_equal(gdf_in,gdf_ref,check_less_precise=True)
	finally:
		# Clean
		for suff in ['cpg','dbf','prj','shp','shx']:
			os.remove('test.'+suff)		


##### to_netcdf() tests #####
def test_to_netcdf(agg=agg):
	# Export to netcdf
	agg.to_netcdf('test.nc')

	try:
		# Make reference dataset
		ds_ref = xr.Dataset({'test':(('poly_idx','run'),[[0.9999629411369734,1.9999629411369735]])},
	               			coords = {'poly_idx':(('poly_idx'),[0]),
	                         'run':(('run'),[0,1])})

		# Load and Test
		with xr.open_dataset('test.nc') as ds_out:
			xr.testing.assert_allclose(ds_ref,ds_out)
		
	finally:
		# Remove test export file 
		os.remove('test.nc')


##### silent export tests #####
@patch('sys.stdout', new_callable=StringIO)
def test_silent_filesaves(mock_stdout):
	try:
		fns = {'nc':'test.nc',
				'csv':'test.csv',
				'shp':'test'}

		with set_options(silent=True):
			agg.to_netcdf(fns['nc'])
			agg.to_csv(fns['csv'])
			agg.to_shp(fns['shp'])
		# Check that nothing was printed
		assert mock_stdout.getvalue() == ''

	finally:
		for k in fns:
			# Clean
			if (k == 'shp') and os.path.exists(fns[k]+'.shp'):
				# Clean shapefile aux files
				for suff in ['cpg','dbf','prj','shp','shx']:
					os.remove(fns[k]+'.'+suff)
			elif os.path.exists(fns[k]):
				os.remove(fns[k])

##### pixel_overlaps() export tests #####
## Create weightmap to export
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

def test_export_wm_nooverwrite(wm_out=wm_out):
	# Test to make sure FileExistsError is thrown if the target 
	# directory already exists and overwrite=False
	fn = 'wm_export_test'

	try:
		# Create temporary directory 
		os.mkdir(fn)

		# Try to save with overwrite=False
		with pytest.raises(FileExistsError):
			wm_out.to_file(fn,overwrite=False)
	finally:
		if os.path.exists(fn):
			# Clean
			shutil.rmtree(fn)

def test_export_wm_standard(wm_out=wm_out):
	# Testing the .to_file() --> read_wm() workflow. 
	# Rather complex because of the many different components of wm. 
	# wm.agg in particular is a dataframe with lists in it (which is 
	# of course not great) which makes it hard to compare. This should
	# work. The main issue is this converts agg to a dataframe from a 
	# geodataframe; but I don't think that changes anything (it doesn't
	# have any geographic information anyways) 

	fn = 'wm_export_test'

	try:
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
	finally:
		##### clean 
		if os.path.exists(fn):
			shutil.rmtree(fn)

