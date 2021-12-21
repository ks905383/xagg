import pytest
import pandas as pd
import numpy as np
import xarray as xr
import geopandas as gpd
import copy
from geopandas import testing as gpdt
from unittest import TestCase
from shapely.geometry import Polygon
import xesmf as xe

from xagg.wrappers import (pixel_overlaps)


##### pixel_overlaps() tests #####
def test_pixel_overlaps_nochange_gdf():
	# To make sure the input gdf isn't changed by any processing
	da = xr.DataArray(data=np.array([[0,1],[2,3]]),
					  coords={'lat':(['lat'],np.array([0,1])),
							'lon':(['lon'],np.array([0,1]))},
					  dims=['lon','lat'])

	# Create polygon covering one pixel
	gdf_test = {'name':['test'],
				'geometry':[Polygon([(-0.5,-0.5),(-0.5,0.5),(0.5,0.5),(0.5,-0.5),(-0.5,-0.5)])]}
	gdf_test = gpd.GeoDataFrame(gdf_test,crs="EPSG:4326")

	# Make deep copy to make sure the thing isn't changed
	gdf_input = copy.deepcopy(gdf_test)

	# Calculate pixel_overlaps through the wrapper function
	wm = pixel_overlaps(da,gdf_test)

	# See if gdf_test has changed
	gpdt.assert_geodataframe_equal(gdf_test,gdf_input)



def test_pixel_overlaps_dataarray():
	# To make sure that pixel_overlaps works with an unnamed dataarray
	da = xr.DataArray(data=np.array([[0,1],[2,3]]),
					  coords={'lat':(['lat'],np.array([0,1])),
							'lon':(['lon'],np.array([0,1]))},
					  dims=['lon','lat'])

	# Create polygon covering one pixel
	gdf_test = {'name':['test'],
				'geometry':[Polygon([(-0.5,-0.5),(-0.5,0.5),(0.5,0.5),(0.5,-0.5),(-0.5,-0.5)])]}
	gdf_test = gpd.GeoDataFrame(gdf_test,crs="EPSG:4326")

	# Calculate pixel_overlaps through the wrapper function, 
	# which should change the dataarray to a dataframe
	wm = pixel_overlaps(da,gdf_test)

	df0 = pd.DataFrame(wm.agg)

	# Define what the output should be
	df_compare = pd.DataFrame({'name':['test'],'poly_idx':0,
									'rel_area':[[[1.0]]],'pix_idxs':[[0]],
									'coords':[[(0,0)]]})


	assert np.allclose([v for v in df0.rel_area],[v for v in df_compare.rel_area])
	assert np.allclose([v for v in df0.pix_idxs],[v for v in df_compare.pix_idxs])
	assert np.allclose([v for v in df0.coords],[v for v in df_compare.coords])


def test_pixel_overlaps_dataarray_wname():
	# To make sure that pixel_overlaps works with a named dataarray
	da = xr.DataArray(data=np.array([[0,1],[2,3]]),
					  coords={'lat':(['lat'],np.array([0,1])),
							'lon':(['lon'],np.array([0,1]))},
					  dims=['lon','lat'],
					  name='tas')

	# Create polygon covering one pixel
	gdf_test = {'name':['test'],
				'geometry':[Polygon([(-0.5,-0.5),(-0.5,0.5),(0.5,0.5),(0.5,-0.5),(-0.5,-0.5)])]}
	gdf_test = gpd.GeoDataFrame(gdf_test,crs="EPSG:4326")

	# Calculate pixel_overlaps through the wrapper function, 
	# which should change the dataarray to a dataframe
	wm = pixel_overlaps(da,gdf_test)

	df0 = pd.DataFrame(wm.agg)

	# Define what the output should be
	df_compare = pd.DataFrame({'name':['test'],'poly_idx':0,
									'rel_area':[[[1.0]]],'pix_idxs':[[0]],
									'coords':[[(0,0)]]})


	assert np.allclose([v for v in df0.rel_area],[v for v in df_compare.rel_area])
	assert np.allclose([v for v in df0.pix_idxs],[v for v in df_compare.pix_idxs])
	assert np.allclose([v for v in df0.coords],[v for v in df_compare.coords])