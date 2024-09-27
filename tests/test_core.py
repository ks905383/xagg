import pytest
import copy
import pandas as pd
import numpy as np
import xarray as xr
import geopandas as gpd
from geopandas import testing as gpdt
import itertools
from unittest import TestCase
from shapely.geometry import Polygon
from unittest.mock import patch
from io import StringIO
try:
    import xesmf as xe
    _has_xesmf=True
except ImportError:
	# To be able to test the rest with environments without xesmf
	_has_xesmf=False

from xagg.core import (process_weights,create_raster_polygons,get_pixel_overlaps,aggregate,NoOverlapError)
from xagg.options import set_options

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

	if _has_xesmf:
		ds_t,weights_info = process_weights(ds,weights=weights)

		ds_compare = xr.Dataset({'weights':(('lat','lon'),np.array([[0,1],[2,3]]))},
								coords={'lat':(['lat'],np.array([0,1])),
										'lon':(['lon'],np.array([0,1])),
								})

		# Check if weights were correctly added to ds
		xr.testing.assert_allclose(ds_compare,ds_t,atol=1e-4)
	else:
		# Should raise ImportError in the no-xesmf environment
		with pytest.raises(ImportError):
			ds_t,weights_info = process_weights(ds,weights=weights)
	
def test_process_weights_close_weights():
	# Make sure weights that are within `np.allclose` but not exactly
	# the same grid as the input ds are correctly allocated
	# Robustness against floating point differences in grids)
	ds = xr.Dataset(coords={'lat':(['lat'],np.array([0,1])),
							'lon':(['lon'],np.array([0,1]))})

	weights = xr.DataArray(data=np.array([[0,1],[2,3]]),
							dims=['lat','lon'],
							coords=[np.array([0,1])+np.random.rand(2)*1e-10,
									np.array([0,1])+np.random.rand(2)*1e-10])

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
					coords={'lat':(['lat'],np.array([0.,1.])),
							'lon':(['lon'],np.array([0.,1.])),
							'bnds':(['bnds'],np.array([0.,1.]))})
	pix_agg = create_raster_polygons(ds)

    # Create comparison geodataframe (what it should come out to)
	poly_test = {'lat':np.array([0.,0.,1.,1.]),'lon':np.array([0.,1.,0.,1.]),
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

	# Synthetic weights grid that requires regridding
	weights = xr.DataArray(data=np.array([[-0.5,0.5,0.5,1.5],
										  [0.5,-0.5,1.5,0.5],
										  [1.5,2.5,2.5,3.5],
										  [2.5,1.5,3.5,2.5]]),
							dims=['lat','lon'],
							coords=[np.array([-0.25,0.25,0.75,1.25]),
									np.array([-0.25,0.25,0.75,1.25])])

	if _has_xesmf:
		pix_agg = create_raster_polygons(ds,weights=weights)

		compare_series = pd.Series(data=[np.array(v) for v in [0.,1.,2.,3.]],
									 			index=[0,1,2,3],
									 			name='weights')


		# There's an issue here in pd.testing.assert_series_equal...
		#pd.testing.assert_series_equal(pix_agg['gdf_pixels'].weights,
		#							   compare_series)
		#np.testing.assert_allclose(compare_series,pix_agg['gdf_pixels'].weights)
		#assert np.allclose(compare_series,pix_agg['gdf_pixels'].weights)
		assert np.allclose([float(v) for v in compare_series],[float(v) for v in pix_agg['gdf_pixels'].weights],
			 				atol=1e-4)
	else:
		# Should raise ImportError in the no-xesmf environment
		with pytest.raises(ImportError):
			pix_agg = create_raster_polygons(ds,weights=weights)

def test_create_raster_polygons_at180():
	# Make sure raster polygons are correctly built at the 180/-180
	# edge (to not create one really long pixel)
	ds = xr.Dataset({'test':(('lat','lon'),np.random.rand(2,3))},
		coords={'lat':(['lat'],np.array([0,1])),
							'lon':(['lon'],np.array([179.2,-179.8, -178.8]))})

	# Create polygon far away in lon, but at the same latitude
	gdf_test = {'name':['test'],
				'geometry':[Polygon([(-0.5,-0.5),(-0.5,0.5),(0.5,0.5),(0.5,-0.5),(-0.5,-0.5)])]}
	gdf_test = gpd.GeoDataFrame(gdf_test,crs="EPSG:4326")

	# Create polygon grid from raster grid cells
	pix_agg = create_raster_polygons(ds,weights=None)

	# There shouldn't be any overlaps between this dataset and those
	# pixels, so we expect the NoOverlapError to happen
	with pytest.raises(NoOverlapError):
		wm_out = get_pixel_overlaps(gdf_test,pix_agg)



##### get_pixel_overlaps() tests #####
# build raster polygons from a simple 2x2x3 grid of lat/lon/time pixels
ds = xr.Dataset({'test':(['lon','lat','time'],np.reshape(np.arange(1,13),(2,2,3))),
				 'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5]])),
				 'lon_bnds':(['lon','bnds'],np.array([[-0.5,0.5],[0.5,1.5]]))},
				coords={'lat':(['lat'],np.array([0,1])),
						'lon':(['lon'],np.array([0,1])),
						'bnds':(['bnds'],np.array([0,1])),
						'time':(['time'],pd.date_range('2019-01-01','2019-01-03'))})

# add a simple weights grid
weights = xr.DataArray(data=np.array([[0,1],[2,3]]),
							dims=['lat','lon'],
							coords=[ds.lat,ds.lon])

# calculate the pix_agg variable tested above, to be used in the 
# tests below
pix_agg = create_raster_polygons(ds,weights=weights)



# Test if the shapefile completely covers one pixel
def test_get_pixel_overlaps_one_pixel(pix_agg=copy.deepcopy(pix_agg)):
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
def test_get_pixel_overlaps_fraction_of_pixel(pix_agg=copy.deepcopy(pix_agg)):
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
def test_get_pixel_overlaps_multiple_pixels_complete(pix_agg=copy.deepcopy(pix_agg)):
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
def test_get_pixel_overlaps_multiple_pixels_partial(pix_agg=copy.deepcopy(pix_agg)):
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
def test_get_pixel_overlaps_passthru_source_grid(pix_agg=copy.deepcopy(pix_agg)):
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
def test_get_pixel_overlaps_gdf_wpreexisting_index(pix_agg=copy.deepcopy(pix_agg),
	 											   ds=copy.deepcopy(ds)):
	# Test to make sure it works with pre-existing indices in the gdf
	# Create polygon covering multiple pixels
	gdf_test = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf_test = gpd.GeoDataFrame(gdf_test,crs="EPSG:4326",index=np.arange(10,11))

	# Get pixel overlaps
	wm_out = get_pixel_overlaps(gdf_test,pix_agg)

	# The index error for an incorrectly-indexed gdf is thrown in aggregate()
	agg = aggregate(ds,wm_out)

	print(agg.agg.test)
	print(agg.agg.test.values)
	pd.testing.assert_series_equal(agg.agg.test,
									pd.Series([[[[7.4999,8.4999,9.4999]]]],
										name='test'),atol=1e-4)


##### aggregate() tests #####
def test_aggregate_basic(ds=ds):
	# Create polygon covering multiple pixels
	gdf = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")

	# calculate the pix_agg variable tested above, to be used in the 
	# tests below
	pix_agg = create_raster_polygons(ds.copy())

	# Get pixel overlaps
	wm = get_pixel_overlaps(gdf,pix_agg)

	# Get aggregate
	agg = aggregate(ds,wm)

	# This requires shifting rtol to 1e-4 for some reason, in that 
	# it's actually 1.499981, whereas multiplying out 
	# np.sum(agg.agg.rel_area[0]*np.array([0,1,2,3]))gives 1.499963... 
	# Possibly worth examining more closely later
	pd.testing.assert_series_equal(agg.agg.test,
									pd.Series([[[[5.4999,6.4999,7.4999]]]],
										name='test'),atol=1e-4)

def test_aggregate_wdataarray(ds=ds):
	da = ds['test'].copy()
	# Create polygon covering multiple pixels
	gdf = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")

	# calculate the pix_agg variable tested above, to be used in the 
	# tests below
	pix_agg = create_raster_polygons(ds.copy())

	# Get pixel overlaps
	wm = get_pixel_overlaps(gdf,pix_agg)

	# Get aggregate
	agg = aggregate(da,wm)

	# This requires shifting rtol to 1e-4 for some reason, in that 
	# it's actually 1.499981, whereas multiplying out 
	# np.sum(agg.agg.rel_area[0]*np.array([0,1,2,3])) gives 1.499963... 
	# Possibly worth examining more closely later
	pd.testing.assert_series_equal(agg.agg.test,
									pd.Series([[[[5.4999,6.4999,7.4999]]]],
										name='test'),atol=1e-4)

def test_aggregate_basic_wdotproduct(ds=ds):
	# Create polygon covering multiple pixels
	gdf = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")

	# calculate the pix_agg variable tested above, to be used in the 
    # tests below
	pix_agg = create_raster_polygons(ds.copy())

    # Get pixel overlaps
	wm = get_pixel_overlaps(gdf,pix_agg,impl='dot_product')

    # Get aggregate
	agg = aggregate(ds,wm,impl='dot_product')

    # Same as above with for loop implementation
	pd.testing.assert_series_equal(agg.agg.test,
									pd.Series([[[[5.4999,6.4999,7.4999]]]],
										name='test'),atol=1e-4)
    

def test_aggregate_twopolys_wdotproduct(ds=ds):
    # Create multiple polygons, to double check, since dot product
    # implementation takes a slightly different approach to indices 
    # in the geodataframe
    gdf = {'name':['test1','test2'],
           'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)]),
                            Polygon([(-1,0),(-1,1),(0,1),(0,0),(-1,0)])]}
    gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")

    # calculate the pix_agg variable tested above, to be used in the 
    # tests below
    pix_agg = create_raster_polygons(ds.copy())

    # Get pixel overlaps
    wm = get_pixel_overlaps(gdf,pix_agg,impl='dot_product')

    # Get aggregate
    agg = aggregate(ds,wm,impl='dot_product')
    
    # This requires shifting rtol to 1e-3 for some reason, in that 
    # it's actually 1.499981, whereas multiplying out 
    # np.sum(agg.agg.rel_area[0]*np.array([0,1,2,3]))gives 1.499963... 
    # Possibly worth examining more closely later
    pd.testing.assert_series_equal(agg.agg.test,
    							 pd.Series([[[[5.4999,6.4999,7.4999]]],
           								    [[[2.4999,3.4999,4.4999]]]],
           								    name='test'),atol=1e-4)


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
	pix_agg = create_raster_polygons(ds.copy(),weights=weights)


	# Get pixel overlaps
	wm_for = get_pixel_overlaps(gdf,pix_agg,impl='for_loop')
	wm_dot = get_pixel_overlaps(gdf,pix_agg,impl='dot_product')

	# Get aggregate
	agg_for = aggregate(ds,wm_for,impl='for_loop')
	agg_dot = aggregate(ds,wm_dot,impl='dot_product')

	# Since the "test" for the input ds has [1,8] for the two 
	# equatorial pixels, the average should just be 4 for timestep 1
	pd.testing.assert_series_equal(agg_for.agg.test,pd.Series([[[[4,5,6]]]],name='test'))
	pd.testing.assert_series_equal(agg_dot.agg.test,pd.Series([[[[4,5,6]]]],name='test'))



def test_aggregate_with_mismatched_grid():
	# This is to see if the subset_find call works

	ds = xr.Dataset({'test':(['lon','lat','time'],np.reshape(np.arange(1,28),(3,3,3))),
				 'lat_bnds':(['lat','bnds'],np.array([[-1.5,-0.5],[-0.5,0.5],[0.5,1.5]])),
				 'lon_bnds':(['lon','bnds'],np.array([[-1.5,-0.5],[-0.5,0.5],[0.5,1.5]]))},
				coords={'lat':(['lat'],np.array([-1,0,1])),
						'lon':(['lon'],np.array([-1,0,1])),
						'bnds':(['bnds'],np.array([0,1])),
						'time':(['time'],pd.date_range('2019-01-01','2019-01-03'))})


	# Create polygon covering multiple pixels
	gdf = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")

	# calculate the pix_agg variable tested above, to be used in the 
	# tests below
	pix_agg = create_raster_polygons(ds.copy())


	# Get pixel overlaps
	wm_for = get_pixel_overlaps(gdf,pix_agg,impl='for_loop')
	wm_dot = get_pixel_overlaps(gdf,pix_agg,impl='dot_product')

	# Get aggregate
	agg_for = aggregate(ds,wm_for,impl='for_loop')
	agg_dot = aggregate(ds,wm_dot,impl='dot_product')

	pd.testing.assert_series_equal(agg_for.agg.test,pd.Series([[[[18.9999,19.9999,20.9999]]]],name='test'),atol=1e-4)
	pd.testing.assert_series_equal(agg_dot.agg.test,pd.Series([[[[18.9999,19.9999,20.9999]]]],name='test'),atol=1e-4)



# Should probably test multiple polygons just to be sure... 


### Higher-D tests
# Create a 4-D dataframe
ds = xr.Dataset({'test':(['lat','lon','time','plev'],np.reshape(np.arange(1,3*2*4*3+1),(2,3,4,3))),
				 'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5]])),
				 'lon_bnds':(['lon','bnds'],np.array([[-0.5,0.5],[0.5,1.5],[1.5,2.5]]))},
				coords={'lat':(['lat'],np.array([0,1])),
						'lon':(['lon'],np.array([0,1,2])),
						'bnds':(['bnds'],np.array([0,1])),
						'time':(['time'],pd.date_range('2019-01-01','2019-01-04')),
                        'plev':(['plev'],np.array([1000,950,900]))})

# Create multiple polygons
gdf = {'name':['test1','test2'],
       'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)]),
                   Polygon([(1,0),(1,1),(2,1),(2,0),(1,0)])]}
gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")

def test_aggregate_4d(ds=ds.copy(),gdf=gdf.copy()):
	# Test 4D aggregation, to make sure it's only aggregating
	# across lat / lon

	pix_agg = create_raster_polygons(ds)


	# Get pixel overlaps
	wm_for = get_pixel_overlaps(gdf,pix_agg,impl='for_loop')
	wm_dot = get_pixel_overlaps(gdf,pix_agg,impl='dot_product')

	# Get aggregate
	agg_for = aggregate(ds,wm_for,impl='for_loop')
	agg_dot = aggregate(ds,wm_dot,impl='dot_product')

	# 
	series_out = pd.Series([[[np.reshape(np.arange(24,36)+0.99933294,(4,3))]],
           					[[np.reshape(np.arange(36,48)+0.99933294,(4,3))]]],name='test')

	pd.testing.assert_series_equal(agg_for.agg.test,series_out)
	pd.testing.assert_series_equal(agg_dot.agg.test,series_out)

def test_aggregate_4d_altdimorders(ds=ds.copy(),gdf=gdf.copy()):
	# Testing to make sure that, no matter the order of dimensions, 
	# the aggregation is always over lat / lon in the right way
	for perm in [k for k in itertools.permutations(['lon','lat','time','plev'])]:
		ds = ds.transpose(*[*perm,'bnds'])
		# Test to make sure dimension order doesn't matter in aggregating
		pix_agg = create_raster_polygons(ds)


		# Get pixel overlaps
		wm_for = get_pixel_overlaps(gdf,pix_agg,impl='for_loop')
		wm_dot = get_pixel_overlaps(gdf,pix_agg,impl='dot_product')

		# Get aggregate
		agg_for = aggregate(ds,wm_for,impl='for_loop')
		agg_dot = aggregate(ds,wm_dot,impl='dot_product')

		# If plev / time (the non-geographic dimensions) are in 
		# different orders, then the output array will be in a different 
		# order. Will test separately in to_dataset() to makes sure the
		# dimensions remain correctly aggregated. 
		if np.where([k == 'plev' for k in perm])[0][0] < np.where([k == 'time' for k in perm])[0][0]:
			series_out = pd.Series([[[np.reshape(np.arange(24,36)+0.99933294,(4,3)).T]],
	           						[[np.reshape(np.arange(36,48)+0.99933294,(4,3)).T]]],name='test')
		else:
			series_out = pd.Series([[[np.reshape(np.arange(24,36)+0.99933294,(4,3))]],
	           						[[np.reshape(np.arange(36,48)+0.99933294,(4,3))]]],name='test')

		pd.testing.assert_series_equal(agg_for.agg.test,series_out)
		pd.testing.assert_series_equal(agg_dot.agg.test,series_out)


### NEED A TEST FOR NAN BEHAVIOR
# 1) That if a region is covered by pixels that are entirely nans, that it spits out only nans
# 2) That if a region is covered *partially* by pixels that are nans, that those are ignored 
# in the aggregation calculation
def test_aggregate_with_all_nans():
	ds = xr.Dataset({'test':(['lon','lat','time'],np.reshape(np.arange(1,13)*np.nan,(2,2,3))),
					 'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5]])),
					 'lon_bnds':(['lon','bnds'],np.array([[-0.5,0.5],[0.5,1.5]]))},
					coords={'lat':(['lat'],np.array([0,1])),
							'lon':(['lon'],np.array([0,1])),
							'bnds':(['bnds'],np.array([0,1])),
							'time':(['time'],pd.date_range('2019-01-01','2019-01-03'))})

	# get aggregation mapping
	pix_agg = create_raster_polygons(ds)

	# Create polygon covering multiple pixels
	gdf = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")


	# Get pixel overlaps
	wm_for = get_pixel_overlaps(gdf,pix_agg,impl='for_loop')
	wm_dot = get_pixel_overlaps(gdf,pix_agg,impl='dot_product')

	# Get aggregate
	agg_for = aggregate(ds,wm_for,impl='for_loop')
	agg_dot = aggregate(ds,wm_dot,impl='dot_product')

	# Should only return nan 
	pd.testing.assert_series_equal(agg_for.agg.test,pd.Series([[[[np.nan,np.nan,np.nan]]]],name='test'))
	pd.testing.assert_series_equal(agg_dot.agg.test,pd.Series([[[[np.nan,np.nan,np.nan]]]],name='test'))

def test_aggregate_with_some_gridnans():
	# This is a test for aggregating across grid nans - i.e., 
	# pixels that are entirely nan, which are dropped. 
	ds = xr.Dataset({'test':(['lon','lat','time'],np.reshape(np.arange(1,13),(2,2,3))),
				 'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5]])),
				 'lon_bnds':(['lon','bnds'],np.array([[-0.5,0.5],[0.5,1.5]]))},
				coords={'lat':(['lat'],np.array([0,1])),
						'lon':(['lon'],np.array([0,1])),
						'bnds':(['bnds'],np.array([0,1])),
						'time':(['time'],pd.date_range('2019-01-01','2019-01-03'))})
	# Make some pixels nan
	ds['test'] = ds['test'].where(ds.lat!=1)

	# get aggregation mapping
	pix_agg = create_raster_polygons(ds)

	# Create polygon covering multiple pixels
	gdf = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")


	# Get pixel overlaps
	wm_for = get_pixel_overlaps(gdf,pix_agg,impl='for_loop')
	wm_dot = get_pixel_overlaps(gdf,pix_agg,impl='dot_product')

	# Get aggregate
	agg_for = aggregate(ds,wm_for,impl='for_loop')
	agg_dot = aggregate(ds,wm_dot,impl='dot_product')

	# Should successfully aggregate
	pd.testing.assert_series_equal(agg_for.agg.test,pd.Series([[[[4,5,6]]]],name='test'))
	pd.testing.assert_series_equal(agg_dot.agg.test,pd.Series([[[[4,5,6]]]],name='test'))


def test_aggregate_with_partialnans():
	# This is a test to make sure the warning comes up when there's 
	# an inconsistent nan pixel - i.e., if a third dimension is time,
	# it's a pixel that's nan sometimes and not other times
	ds = xr.Dataset({'test':(['lon','lat','time'],np.reshape(np.arange(1,13),(2,2,3))),
				 'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5]])),
				 'lon_bnds':(['lon','bnds'],np.array([[-0.5,0.5],[0.5,1.5]]))},
				coords={'lat':(['lat'],np.array([0,1])),
						'lon':(['lon'],np.array([0,1])),
						'bnds':(['bnds'],np.array([0,1])),
						'time':(['time'],pd.date_range('2019-01-01','2019-01-03'))})
	ds['test'] = ds['test'].where(((ds.lon!=1) | (ds.lat!=1) | (ds.time!=pd.to_datetime('2019-01-01'))))

	# get aggregation mapping
	pix_agg = create_raster_polygons(ds)

	# Create polygon covering multiple pixels
	gdf = {'name':['test'],
				'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
	gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")


	# Get pixel overlaps
	wm_for = get_pixel_overlaps(gdf,pix_agg,impl='for_loop')
	wm_dot = get_pixel_overlaps(gdf,pix_agg,impl='dot_product')

	# Get aggregate
	with pytest.warns(UserWarning):
		agg_for = aggregate(ds,wm_for,impl='for_loop')
	with pytest.warns(UserWarning):
		agg_dot = aggregate(ds,wm_dot,impl='dot_product')

##### aggregate() silencing tests #####
# Create polygon covering multiple pixels
gdf = {'name':['test'],
			'geometry':[Polygon([(0,0),(0,1),(1,1),(1,0),(0,0)])]}
gdf = gpd.GeoDataFrame(gdf,crs="EPSG:4326")

# Get pixel overlaps
with set_options(silent=True):
	wm = get_pixel_overlaps(gdf,pix_agg)


@patch('sys.stdout', new_callable=StringIO)
def test_aggregate_silent_true(mock_stdout):
	with set_options(silent=True):
		# Get aggregate
		agg = aggregate(ds,wm)


	# Check that nothing was printed
	assert mock_stdout.getvalue() == ''

@patch('sys.stdout', new_callable=StringIO)
def test_aggregate_silent_false(mock_stdout):
	with set_options(silent=False):
		# Get aggregate
		agg = aggregate(ds,wm)


	# Check that nothing was printed
	assert mock_stdout.getvalue() != ''


