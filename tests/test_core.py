import pytest
import pandas as pd
import numpy as np
import xarray as xr
import geopandas as gpd
from geopandas import testing as gpdt
from unittest import TestCase
from shapely.geometry import Polygon
import xesmf as xe

from xagg.core import (process_weights,create_raster_polygons,get_pixel_overlaps,aggregate)


##### process_weights() tests #####
def test_process_weights_null():
	# (by default, process_weights has weights==None, which should 
	# return the original ds, and weights_info = 'nowghts')
	ds = xr.Dataset(coords={'lat':(['lat'],np.array([0,1])),
							'lon':(['lon'],np.array([0,1])),
							})
	ds_t,weights_info = process_weights(ds)

	assert weights_info == 'nowghts'
	xr.testing.assert_allclose(ds,ds_t) # Check that the ds was untouched

def test_process_weights_basic():
	# Now, test with a weights array the same size as the ds (so, 
	# no regridding)
	ds = xr.Dataset(coords={'lat':(['lat'],np.array([0,1])),
							'lon':(['lon'],np.array([0,1]))})

	weights = xr.DataArray(data=np.array([[0,1],[2,3]]),
							dims=['lat','lon'],
							coords=[ds.lat,ds.lon])

	#rgrd = xe.Regridder(ds,weights,'bilinear')

	ds_t,weights_info = process_weights(ds,weights=weights)

	ds_compare = xr.Dataset({'weights':(('lat','lon'),np.array([[0,1],[2,3]]))},
							coords={'lat':(['lat'],np.array([0,1])),
									'lon':(['lon'],np.array([0,1])),
							})

	# Check if weights were correctly added to ds
	xr.testing.assert_allclose(ds_compare,ds_t)
	# (weights_info isn't currently used by anything)

def test_process_weights_regrid_weights():
	# Now, test with a weights array twice the resolution as the
	# ds, so it needs to be regridded
	ds = xr.Dataset(coords={'lat':(['lat'],np.array([0,1])),
							'lon':(['lon'],np.array([0,1])),
							})

	# Synthetic weights grid, with double the resolution, and shifted
	# by a half degree. Should regrid to the same weights grid as above
	weights = xr.DataArray(data=np.array([[-0.5,0.5,0.5,1.5],
										  [0.5,-0.5,1.5,0.5],
										  [1.5,2.5,2.5,3.5],
										  [2.5,1.5,3.5,2.5]]),
							dims=['lat','lon'],
							coords=[np.array([-0.25,0.25,0.75,1.25]),
									np.array([-0.25,0.25,0.75,1.25])])

	ds_t,weights_info = process_weights(ds,weights=weights)

	ds_compare = xr.Dataset({'weights':(('lat','lon'),np.array([[0,1],[2,3]]))},
							coords={'lat':(['lat'],np.array([0,1])),
									'lon':(['lon'],np.array([0,1])),
							})

	# Check if weights were correctly added to ds
	xr.testing.assert_allclose(ds_compare,ds_t)
	



##### create_raster_polygons() tests #####
def test_create_raster_polygons_basic():
	# build raster polygons from a simple 2x2 grid of lat/lon pixels
	ds = xr.Dataset({'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5]])),
					 'lon_bnds':(['lon','bnds'],np.array([[-0.5,0.5],[0.5,1.5]]))},
					coords={'lat':(['lat'],np.array([0,1])),
							'lon':(['lon'],np.array([0,1])),
							'bnds':(['bnds'],np.array([0,1]))})
	pix_agg = create_raster_polygons(ds)

	# Create comparison geodataframe (what it should come out to)
	poly_test = {'lat':[np.array(v) for v in [0.,0.,1.,1.]],'lon':[np.array(v) for v in [0.,1.,0.,1.]],
						'geometry':[Polygon([(-0.5,-0.5),(-0.5,0.5),(0.5,0.5),(0.5,-0.5),(-0.5,-0.5)]),
								   Polygon([(0.5, -0.5), (0.5, 0.5), (1.5, 0.5), (1.5, -0.5), (0.5, -0.5)]),
								   Polygon([(-0.5, 0.5), (-0.5, 1.5), (0.5, 1.5), (0.5, 0.5), (-0.5, 0.5)]),
								   Polygon([(0.5, 0.5), (0.5, 1.5), (1.5, 1.5), (1.5, 0.5), (0.5, 0.5)])],
						'pix_idx':[0,1,2,3]}
	poly_test = gpd.GeoDataFrame(poly_test,crs="EPSG:4326")

	# Use gpd assert (there's a check_less_precise option for "almost equal", 
	# but that doesn't seem necessary for this simple example)
	gpdt.assert_geodataframe_equal(poly_test,pix_agg['gdf_pixels'])


	# Test second bit of create_raster_polygons output (that the 
	# 'source_grid' guy gives you the stacked lat/lon coordinates)
	source_grid_test = {'lat':ds.stack(loc=('lat','lon')).lat,'lon':ds.stack(loc=('lat','lon')).lon}
	xr.testing.assert_allclose(source_grid_test['lat'],pix_agg['source_grid']['lat'])
	xr.testing.assert_allclose(source_grid_test['lon'],pix_agg['source_grid']['lon'])


def test_create_raster_polygons_with_weights():
	# build raster polygons from a simple 2x2 grid of lat/lon pixels
	ds = xr.Dataset({'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5]])),
					 'lon_bnds':(['lon','bnds'],np.array([[-0.5,0.5],[0.5,1.5]]))},
					coords={'lat':(['lat'],np.array([0,1])),
							'lon':(['lon'],np.array([0,1])),
							'bnds':(['bnds'],np.array([0,1]))})

	# Synethetic weights grid
	weights = xr.DataArray(data=np.array([[-0.5,0.5,0.5,1.5],
										  [0.5,-0.5,1.5,0.5],
										  [1.5,2.5,2.5,3.5],
										  [2.5,1.5,3.5,2.5]]),
							dims=['lat','lon'],
							coords=[np.array([-0.25,0.25,0.75,1.25]),
									np.array([-0.25,0.25,0.75,1.25])])

	pix_agg = create_raster_polygons(ds,weights=weights)

	compare_series = pd.Series(data=[np.array(v) for v in [0.,1.,2.,3.]],
								 			index=[0,1,2,3],
								 			name='weights')


	# There's an issue here in pd.testing.assert_series_equal...
	#pd.testing.assert_series_equal(pix_agg['gdf_pixels'].weights,
	#							   compare_series)
	#np.testing.assert_allclose(compare_series,pix_agg['gdf_pixels'].weights)
	#assert np.allclose(compare_series,pix_agg['gdf_pixels'].weights)
	assert np.allclose([float(v) for v in compare_series],[float(v) for v in pix_agg['gdf_pixels'].weights])


##### get_pixel_overlaps() tests #####
# build raster polygons from a simple 2x2 grid of lat/lon pixels
ds = xr.Dataset({'test':(['lon','lat'],np.array([[0,1],[2,3]])),
				 'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5]])),
				 'lon_bnds':(['lon','bnds'],np.array([[-0.5,0.5],[0.5,1.5]]))},
				coords={'lat':(['lat'],np.array([0,1])),
						'lon':(['lon'],np.array([0,1])),
						'bnds':(['bnds'],np.array([0,1]))})

# add a simple weights grid
weights = xr.DataArray(data=np.array([[0,1],[2,3]]),
							dims=['lat','lon'],
							coords=[ds.lat,ds.lon])

# calculate the pix_agg variable tested above, to be used in the 
# tests below
pix_agg = create_raster_polygons(ds,weights=weights)



# Test if the shapefile completely covers one pixel
def test_get_pixel_overlaps_one_pixel(pix_agg=pix_agg):
	# Create polygon covering one pixel
	gdf_test = {'name':['test'],
				'geometry':[Polygon([(-0.5,-0.5),(-0.5,0.5),(0.5,0.5),(0.5,-0.5),(-0.5,-0.5)])]}
	gdf_test = gpd.GeoDataFrame(gdf_test,crs="EPSG:4326")

	# Get pixel overlaps; store output as pandas dataframe (no
	# more geometry info relevant here)
	wm_out = get_pixel_overlaps(gdf_test,pix_agg)
	df0 = pd.DataFrame(wm_out.agg)

	# Define what the output should be
	df_compare = pd.DataFrame({'name':['test'],'poly_idx':0,
									'rel_area':[[[1.0]]],'pix_idxs':[[0]],
									'coords':[[(0,0)]]})


	# Since the elements of some of these data frame columns are lists, 
	# pd.testing.assert_frame_equal() fails on a ValueError ("the truth value
	# of a Series is ambiguous..."). This is the current way around it; but 
	# it's not very robust... Maybe this should result in a rethink of how
	# the geodataframe is organized in the future? 
	assert np.allclose([v for v in df0.rel_area],[v for v in df_compare.rel_area])
	assert np.allclose([v for v in df0.pix_idxs],[v for v in df_compare.pix_idxs])
	assert np.allclose([v for v in df0.coords],[v for v in df_compare.coords])


# Test if the shapefile is a fraction of one pixel
def test_get_pixel_overlaps_fraction_of_pixel(pix_agg=pix_agg):
	# Create polygon covering less than one pixel
	gdf_test = {'name':['test'],
				'geometry':[Polygon([(-0.5,-0.5),(-0.5,0),(0,0),(0,-0.5),(-0.5,-0.5)])]}
	gdf_test = gpd.GeoDataFrame(gdf_test,crs="EPSG:4326")

	# Get pixel overlaps; store output as pandas dataframe (no
	# more geometry info relevant here)
	wm_out = get_pixel_overlaps(gdf_test,pix_agg)
	df0 = pd.DataFrame(wm_out.agg)

	# Define what the output should be
	df_compare = pd.DataFrame({'name':['test'],'poly_idx':0,
									'rel_area':[[[1.0]]],'pix_idxs':[[0]],
									'coords':[[(0,0)]]})


	# Since the elements of some of these data frame columns are lists, 
	# pd.testing.assert_frame_equal() fails on a ValueError ("the truth value
	# of a Series is ambiguous..."). This is the current way around it; but 
	# it's not very robust... Maybe this should result in a rethink of how
	# the geodataframe is organized in the future? 
	assert np.allclose([v for v in df0.rel_area],[v for v in df_compare.rel_area])
	assert np.allclose([v for v in df0.pix_idxs],[v for v in df_compare.pix_idxs])
	assert np.allclose([v for v in df0.coords],[v for v in df_compare.coords])

# Test if the shapefile perfectly covers several complete pixels
def test_get_pixel_overlaps_multiple_pixels_complete(pix_agg=pix_agg):
	# Create polygon covering multiple pixels
	gdf_test = {'name':['test'],
				'geometry':[Polygon([(-0.5,-0.5),(-0.5,1.5),(1.5,1.5),(1.5,-0.5),(-0.5,-0.5)])]}
	gdf_test = gpd.GeoDataFrame(gdf_test,crs="EPSG:4326")

	# Get pixel overlaps; store output as pandas dataframe (no
	# more geometry info relevant here)
	wm_out = get_pixel_overlaps(gdf_test,pix_agg)
	df0 = pd.DataFrame(wm_out.agg)

	# Define what the output should be
	# (rel_area isn't 0.25, 0.25, 0.25, 0.25 because of slight 
	# latitude differences)
	df_compare = pd.DataFrame({'name':['test'],'poly_idx':0,
									'rel_area':[[[0.250019,0.250019,0.249981,0.249981]]],'pix_idxs':[[0,1,2,3]],
									'coords':[[(0,0),(0,1),(1,0),(1,1)]]})


	# Since the elements of some of these data frame columns are lists, 
	# pd.testing.assert_frame_equal() fails on a ValueError ("the truth value
	# of a Series is ambiguous..."). This is the current way around it; but 
	# it's not very robust... Maybe this should result in a rethink of how
	# the geodataframe is organized in the future? 
	assert np.allclose([v for v in df0.rel_area],[v for v in df_compare.rel_area])
	assert np.allclose([v for v in df0.pix_idxs],[v for v in df_compare.pix_idxs])
	assert np.allclose([v for v in df0.coords],[v for v in df_compare.coords])


# Test if the shapefile covers parts of several complete pixels
def test_get_pixel_overlaps_multiple_pixels_partial(pix_agg=pix_agg):
	# Create polygon covering multiple pixels
	gdf_test = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf_test = gpd.GeoDataFrame(gdf_test,crs="EPSG:4326")

	# Get pixel overlaps; store output as pandas dataframe (no
	# more geometry info relevant here)
	wm_out = get_pixel_overlaps(gdf_test,pix_agg)
	df0 = pd.DataFrame(wm_out.agg)

	# Define what the output should be
	# (rel_area isn't 0.25, 0.25, 0.25, 0.25 because of slight 
	# latitude differences)
	df_compare = pd.DataFrame({'name':['test'],'poly_idx':0,
									'rel_area':[[[0.250009,0.250009,0.249991,0.249991]]],'pix_idxs':[[0,1,2,3]],
									'coords':[[(0,0),(0,1),(1,0),(1,1)]]})


	# Since the elements of some of these data frame columns are lists, 
	# pd.testing.assert_frame_equal() fails on a ValueError ("the truth value
	# of a Series is ambiguous..."). This is the current way around it; but 
	# it's not very robust... Maybe this should result in a rethink of how
	# the geodataframe is organized in the future? 
	assert np.allclose([v for v in df0.rel_area],[v for v in df_compare.rel_area])
	assert np.allclose([v for v in df0.pix_idxs],[v for v in df_compare.pix_idxs])
	assert np.allclose([v for v in df0.coords],[v for v in df_compare.coords])

# Make sure the source_grid is passed without change
def test_get_pixel_overlaps_passthru_source_grid(pix_agg=pix_agg):
	# Create polygon covering multiple pixels
	gdf_test = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf_test = gpd.GeoDataFrame(gdf_test,crs="EPSG:4326")

	# Get pixel overlaps; store output as pandas dataframe (no
	# more geometry info relevant here)
	wm_out = get_pixel_overlaps(gdf_test,pix_agg)

	xr.testing.assert_equal(wm_out.source_grid['lat'],pix_agg['source_grid']['lat'])
	xr.testing.assert_equal(wm_out.source_grid['lon'],pix_agg['source_grid']['lon'])

def test_get_pixel_overlaps_passthru_weights(pix_agg=pix_agg):
	# Create polygon covering multiple pixels
	gdf_test = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf_test = gpd.GeoDataFrame(gdf_test,crs="EPSG:4326")

	# Get pixel overlaps; store output as pandas dataframe (no
	# more geometry info relevant here)
	wm_out = get_pixel_overlaps(gdf_test,pix_agg)
	pd.testing.assert_series_equal(wm_out.weights,pix_agg['gdf_pixels'].weights)

# The last remaining attribute of wm_out, wm_out.geometry, is currently only
# saved for reference and not used by anything; so I won't test its passing. 

# Should probably test multiple polygons just to be sure... 


###### get_pixel_overlaps() and aggregate() coupling tests #####
def test_get_pixel_overlaps_gdf_wpreexisting_index(pix_agg=pix_agg):
	# Test to make sure it works with pre-existing indices in the gdf
	# Create polygon covering multiple pixels
	gdf_test = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf_test = gpd.GeoDataFrame(gdf_test,crs="EPSG:4326",index=np.arange(10,11))

	# Get pixel overlaps
	wm_out = get_pixel_overlaps(gdf_test,pix_agg)

	# The index error for an incorrectly-indexed gdf is thrown in aggregate()
	agg = aggregate(ds,wm_out)

	# this assert uses 2.1666 because of the weighting that creates 
	# the pix_agg variable that this whole section has used. Doesn't really 
	# matter, since this is testing an index error that would've 
	# happened during aggregate() above. 
	assert np.allclose([v for v in agg.agg.test.values],2.1666,rtol=1e-4)


##### aggregate() tests #####


def test_aggregate_basic(ds=ds):
	# Create polygon covering multiple pixels
	gdf = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")

	# calculate the pix_agg variable tested above, to be used in the 
	# tests below
	pix_agg = create_raster_polygons(ds)


	# Get pixel overlaps
	wm = get_pixel_overlaps(gdf,pix_agg)

	# Get aggregate
	agg = aggregate(ds,wm)

	# This requires shifting rtol to 1e-4 for some reason, in that 
	# it's actually 1.499981, whereas multiplying out 
	# np.sum(agg.agg.rel_area[0]*np.array([0,1,2,3]))gives 1.499963... 
	# Possibly worth examining more closely later
	assert np.allclose([v for v in agg.agg.test.values],1.4999,rtol=1e-4)


def test_aggregate_with_weights(ds=ds):
	# Create polygon covering multiple pixels
	gdf = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")

	# add a simple weights grid (equator pixels have weight 1, 
	# 1 N pixels have weight 0)
	weights = xr.DataArray(data=np.array([[1,1],[0,0]]),
								dims=['lat','lon'],
								coords=[ds.lat,ds.lon])

	# calculate the pix_agg variable tested above, to be used in the 
	# tests below
	pix_agg = create_raster_polygons(ds,weights=weights)


	# Get pixel overlaps
	wm = get_pixel_overlaps(gdf,pix_agg)

	# Get aggregate
	agg = aggregate(ds,wm)

	# Since the "test" for the input ds has [0,2] for the two 
	# equatorial pixels, the average should just be 1.0
	assert np.allclose([v for v in agg.agg.test.values],1.0)



def test_aggregate_with_mismatched_grid():
	# This is to see if the subset_find call works

	ds = xr.Dataset({'test':(['lon','lat'],np.array([[30,40,50],[10,0,1],[20,2,3]])),
				 'lat_bnds':(['lat','bnds'],np.array([[-1.5,-0.5],[-0.5,0.5],[0.5,1.5]])),
				 'lon_bnds':(['lon','bnds'],np.array([[-1.5,-0.5],[-0.5,0.5],[0.5,1.5]]))},
				coords={'lat':(['lat'],np.array([-1,0,1])),
						'lon':(['lon'],np.array([-1,0,1])),
						'bnds':(['bnds'],np.array([0,1]))})


	# Create polygon covering multiple pixels
	gdf = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")

	# calculate the pix_agg variable tested above, to be used in the 
	# tests below
	pix_agg = create_raster_polygons(ds)


	# Get pixel overlaps
	wm = get_pixel_overlaps(gdf,pix_agg)

	# Get aggregate
	agg = aggregate(ds,wm)

	# On change in rtol, see note in test_aggregate_basic
	assert np.allclose([v for v in agg.agg.test.values],1.4999,rtol=1e-4)


# Should probably test multiple polygons just to be sure... 


### NEED A TEST FOR NAN BEHAVIOR
# 1) That if a region is covered by pixels that are entirely nans, that it spits out only nans
# 2) That if a region is covered *partially* by pixels that are nans, that those are ignored 
# in the aggregation calculation
def test_aggregate_with_all_nans():
	ds = xr.Dataset({'test':(['lon','lat'],np.array([[np.nan,np.nan],[np.nan,np.nan]])),
					 'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5]])),
					 'lon_bnds':(['lon','bnds'],np.array([[-0.5,0.5],[0.5,1.5]]))},
					coords={'lat':(['lat'],np.array([0,1])),
							'lon':(['lon'],np.array([0,1])),
							'bnds':(['bnds'],np.array([0,1]))})

	# get aggregation mapping
	pix_agg = create_raster_polygons(ds)

	# Create polygon covering multiple pixels
	gdf = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")


	# Get pixel overlaps
	wm = get_pixel_overlaps(gdf,pix_agg)

	# Get aggregate
	agg = aggregate(ds,wm)

	# Should only return nan 
	# (this is not a great assert - but agg.agg.test[0] comes out as [array(nan)], 
	# which... I'm not entirely sure how to reproduce. It quaks like a single nan,
	# but it's unclear to me how to get it to work)
	assert np.all([np.isnan(k) for k in agg.agg.test])

def test_aggregate_with_some_nans():
	ds = xr.Dataset({'test':(['lon','lat'],np.array([[np.nan,1],[2,np.nan]])),
					 'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5]])),
					 'lon_bnds':(['lon','bnds'],np.array([[-0.5,0.5],[0.5,1.5]]))},
					coords={'lat':(['lat'],np.array([0,1])),
							'lon':(['lon'],np.array([0,1])),
							'bnds':(['bnds'],np.array([0,1]))})

	# get aggregation mapping
	pix_agg = create_raster_polygons(ds)

	# Create polygon covering multiple pixels
	gdf = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")


	# Get pixel overlaps
	wm = get_pixel_overlaps(gdf,pix_agg)

	# Get aggregate
	agg = aggregate(ds,wm)

	# Should be 1.5; with one pixel valued 1, one pixel valued 2. 
	assert np.allclose([agg.agg.test[0]],1.5,rtol=1e-4)







