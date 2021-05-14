import pytest
import numpy as np
import xarray as xr
from xagg.aux import (normalize,fix_ds,get_bnds,subset_find)

##### normalize() tests #####
def test_normalize():
	# xa.normalize() just divides each element of an array by the sum of the array
	# I'm pretty sure I only use this with numpy arrays - but double-check 
	assert np.allclose(normalize(np.array([1,1])),np.array([0.5,0.5]))


##### fix_ds tests #####
def test_fix_ds_null():
	# If lat, lon correctly named, lon as -180:180, bounds correctly named, sorted, that nothing happens
	ds = xr.Dataset({'test':(['lat','lon'],np.array([[0,1,2],[3,4,5],[6,7,8]])),
					 'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5],[1.5,2.5]])),
					 'lon_bnds':(['lon','bnds'],np.array([[-1.5,-0.5],[-0.5,0.5],[0.5,1.5]]))},
					coords={'lat':(['lat'],np.array([0,1,2])),
							'lon':(['lon'],np.array([-1,0,1]))})
	xr.testing.assert_allclose(ds,fix_ds(ds))

def test_fix_ds_identity():
	# Make sure fix_ds(fix_ds(ds)) == fix_ds(ds)
	ds = xr.Dataset({'test':(['lat','lon'],np.array([[0,1,2],[3,4,5],[6,7,8]]))},
					coords={'lat':(['lat'],np.array([2,1,0])),
							'lon':(['lon'],np.array([4,5,3]))})
	xr.testing.assert_allclose(fix_ds(ds),fix_ds(fix_ds(ds)))

def test_fix_ds_rename():
	# Test whether the renaming to "lat"/"lon" works
	ds = xr.Dataset(coords={'Latitude':(['Latitude'],np.array([0,1,2])),
							'Longitude':(['Longitude'],np.array([0,1,2]))})
	ds = fix_ds(ds)
	assert ('lat' in ds.dims.keys())
	assert ('lon' in ds.dims.keys())

def test_fix_ds_rename_bnds():
	# Test whether the bounds are renamed correctly
	# (currently only triggers when lat/lon are wrong... is this an issue?)
	ds = xr.Dataset({'latitude_bnds':(['latitude','bnds'],np.array([[-0.5,0.5],[0.5,1.5],[1.5,2.5]])),
					 'longitude_bnds':(['longitude','bnds'],np.array([[-0.5,0.5],[0.5,1.5],[1.5,2.5]]))},
					coords={'latitude':(['latitude'],np.array([0,1,2])),
							'longitude':(['longitude'],np.array([0,1,2])),
							'bnds':(['bnds'],np.array([0,1]))})
	ds = fix_ds(ds)
	assert ('lat_bnds' in ds.var())
	assert ('lon_bnds' in ds.var())

def test_fix_ds_lon_wraparound():
	# Test the longitude re-mapping from 0:360 to -180:180
	ds = xr.Dataset({'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5],[1.5,2.5]])),
					 'lon_bnds':(['lon','bnds'],np.array([[359.5,0.5],[178.5,179.5],[179.5,180.5],[180.5,181.5]]))},
					 coords={'lat':(['lat'],np.array([0,1,2])),
							'lon':(['lon'],np.array([0,179,180,181]))})
	ds = fix_ds(ds)

	# What that should be (also sorted!)
	ds_compare = xr.Dataset({'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5],[1.5,2.5]])),
					 		 'lon_bnds':(['lon','bnds'],np.array([[179.5,-179.5],[-179.5,-178.5],[-0.5,0.5],[178.5,179.5]]))},
					 		coords={'lat':(['lat'],np.array([0,1,2])),
									'lon':(['lon'],np.array([-180,-179,0,179]))})

	xr.testing.assert_allclose(ds,ds_compare)

def test_fix_ds_coord_sorting():
	# Test sorting behavior
	ds = xr.Dataset({'test':(['lat','lon'],np.array([[0,1,2],[3,4,5],[6,7,8]]))},
					coords={'lat':(['lat'],np.array([2,1,0])),
							'lon':(['lon'],np.array([4,5,3]))})
	ds = fix_ds(ds)

	# What that should be (also sorted!)
	ds_compare = xr.Dataset({'test':(['lat','lon'],np.array([[8,6,7],[5,3,4],[2,0,1]]))},
							coords={'lat':(['lat'],np.array([0,1,2])),
									'lon':(['lon'],np.array([3,4,5]))})

	xr.testing.assert_allclose(ds,ds_compare)



###### get_bnds tests #####
def test_get_bnds_null():
	# If there are bounds already, do nothing
	ds = xr.Dataset({'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5],[1.5,2.5]])),
					 'lon_bnds':(['lon','bnds'],np.array([[-0.5,0.5],[0.5,1.5],[1.5,2.5]]))},
					coords={'lat':(['lat'],np.array([0,1,2])),
							'lon':(['lon'],np.array([0,1,2])),
							'bnds':(['bnds'],np.array([0,1]))})

	xr.testing.assert_allclose(ds,get_bnds(ds))

def test_get_bnds_basic():
	# Basic attempt to get the bounds (far away from the wraparound)
	ds = xr.Dataset(coords={'lat':(['lat'],np.array([0,1,2])),
							'lon':(['lon'],np.array([0,1,2]))})
	ds = get_bnds(ds)

	ds_compare = xr.Dataset({'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5],[1.5,2.5]])),
					 		 'lon_bnds':(['lon','bnds'],np.array([[-0.5,0.5],[0.5,1.5],[1.5,2.5]]))},
							coords={'lat':(['lat'],np.array([0,1,2])),
									'lon':(['lon'],np.array([0,1,2])),
									'bnds':(['bnds'],np.array([0,1]))})

	xr.testing.assert_allclose(ds,ds_compare)
	

def test_get_bnds_wraparound():
	# Basic attempt to get the bounds (far away from the wraparound)
	ds = xr.Dataset(coords={'lat':(['lat'],np.arange(-89.4,89.7)),
							'lon':(['lon'],np.arange(-179.4,179.7))})
	ds = get_bnds(ds)

	lats = np.array(list(zip(np.arange(-89.9,89.91),np.arange(-88.9,90.9))))
	lats[-1,1] = 90.0
	lons = np.array(list(zip(np.arange(-179.9,179.91),np.arange(-178.9,180.9))))
	lons[-1,1] = -179.9
	ds_compare = xr.Dataset({'lat_bnds':(['lat','bnds'],lats),
					 		 'lon_bnds':(['lon','bnds'],lons)},
							coords={'lat':(['lat'],np.arange(-89.4,89.7)),
									'lon':(['lon'],np.arange(-179.4,179.7)),
									'bnds':(['bnds'],np.array([0,1]))})

	xr.testing.assert_allclose(ds,ds_compare)

#def test_get_bnds_badwraparound():
	# THERE ARE A LOT OF ISSUES LEFT WITH THIS - 360s showing up somehow, even though it's not in the "edges"... idk man
	# Try a test taht does some weird subset, like (-178:183). Also make sure that -89.9 gets turned to -90. 



###### subset_find tests #####
def test_subset_find_basic():
	ds0 = xr.Dataset({'test':(['lat','lon'],np.array([[0,1,2],[3,4,5],[6,7,8]]))},
					 coords={'lat':(['lat'],np.array([0,1,2])),
							'lon':(['lon'],np.array([-1,0,1]))})

	ds1 = xr.Dataset({'lat':(['lat'],np.array([0,1])),
					  'lon':(['lon'],np.array([-1,0]))})
	ds1 = ds1.stack(loc=('lat','lon'))

	ds_compare = xr.Dataset({'test':(['lat','lon'],np.array([[0,1],[3,4]]))},
							coords={'lat':(['lat'],np.array([0,1])),
									'lon':(['lon'],np.array([-1,0]))})

	xr.testing.assert_allclose(subset_find(ds0,ds1),ds_compare)




