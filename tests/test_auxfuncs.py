import pytest
import numpy as np
import xarray as xr
from xagg.auxfuncs import (normalize,fix_ds,get_bnds,subset_find)

##### normalize() tests #####
def test_normalize():
	# xa.normalize() just divides each element of an array by the sum of the array
	# I'm pretty sure I only use this with numpy arrays - but double-check 
	assert np.allclose(normalize(np.array([1,1])),np.array([0.5,0.5]))

def test_normalize_all0s():
	# Should return a vector of nans if all elements of the input vector are 0
	test_vec = np.array([0,0])

	norm_vec = normalize(test_vec)

	assert np.allclose(norm_vec,np.array([np.nan,np.nan]),
						equal_nan=True)

def test_normalize_dropnans():
	# Make sure nans are accurately dropped
	test_vec = np.array([2,4,np.nan,4])

	norm_vec = normalize(test_vec,drop_na=True)

	assert np.allclose(norm_vec,np.array([0.2,0.4,np.nan,0.4]),
					   equal_nan=True)


##### fix_ds() tests #####
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
	assert ('lat' in ds.sizes)
	assert ('lon' in ds.sizes)

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



###### get_bnds() tests #####
def test_get_bnds_null():
	# If there are bounds already, do nothing
	ds = xr.Dataset({'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5],[1.5,2.5]])),
					 'lon_bnds':(['lon','bnds'],np.array([[-0.5,0.5],[0.5,1.5],[1.5,2.5]]))},
					coords={'lat':(['lat'],np.array([0,1,2])),
							'lon':(['lon'],np.array([0,1,2]))})

	xr.testing.assert_allclose(ds,get_bnds(ds))

def test_get_bnds_basic():
	# Basic attempt to get the bounds (far away from the wraparound)
	ds = xr.Dataset(coords={'lat':(['lat'],np.array([0,1,2])),
							'lon':(['lon'],np.array([0,1,2]))})
	ds = get_bnds(ds)

	ds_compare = xr.Dataset({'lat_bnds':(['lat','bnds'],np.array([[-0.5,0.5],[0.5,1.5],[1.5,2.5]])),
					 		 'lon_bnds':(['lon','bnds'],np.array([[-0.5,0.5],[0.5,1.5],[1.5,2.5]]))},
							coords={'lat':(['lat'],np.array([0,1,2])),
									'lon':(['lon'],np.array([0,1,2]))})

	xr.testing.assert_allclose(ds,ds_compare)
	

def test_get_bnds_fullgrid():
	# -179.5:179.5, with bounds exactly at 180/-180
	ds = xr.Dataset(coords={'lat':(['lat'],np.arange(-89.5,89.51)), 
	                    	'lon':(['lon'],np.arange(-179.5,179.51))})
	ds = get_bnds(ds)

	lat_bnds = np.array(list(zip(np.arange(-90,89.91),np.arange(-89,90.1))))
	lon_bnds = np.array(list(zip(np.arange(-180,179.01),np.arange(-179,180.01))))
	ds_compare = xr.Dataset({'lat_bnds':(['lat','bnds'],lat_bnds),
	                         'lon_bnds':(['lon','bnds'],lon_bnds)},
	                        coords={'lat':(['lat'],np.arange(-89.5,89.51)),
	                                'lon':(['lon'],np.arange(-179.5,179.51))})

	xr.testing.assert_allclose(ds,ds_compare)

def test_get_bnds_truncatedlats():
	# Tests to make sure grids are truncated to 90/-90 when needed
	ds = xr.Dataset(coords={'lat':(['lat'],np.arange(-90,90.01)), 
	                    	'lon':(['lon'],np.arange(-180,179.01))})
	ds = get_bnds(ds)

	lat_bnds = np.array(list(zip(np.arange(-90.5,89.51),np.arange(-89.5,90.51))))
	lat_bnds[0,0] = -90; lat_bnds[-1,-1] = 90
	lon_bnds = np.array(list(zip(np.arange(-180.5,178.51),np.arange(-179.5,179.51))))
	lon_bnds[0,0] = 179.5
	ds_compare = xr.Dataset({'lat_bnds':(['lat','bnds'],lat_bnds),
	                         'lon_bnds':(['lon','bnds'],lon_bnds)},
	                        coords={'lat':(['lat'],np.arange(-90,90.01)),
	                                'lon':(['lon'],np.arange(-180,179.01))})

	xr.testing.assert_allclose(ds,ds_compare)

def test_get_bnds_partialgrid():
	# Partial grid, that should _not_ be wrapped around (farther
	# than the dynamic wrap_around limit = 2*median lon diff)
	ds = xr.Dataset(coords={'lat':(['lat'],np.arange(-89.5,89.51)), 
	                    	'lon':(['lon'],np.arange(-179.5,177.51))})
	ds = get_bnds(ds)

	lat_bnds = np.array(list(zip(np.arange(-90,89.91),np.arange(-89,90.1))))
	lon_bnds = np.array(list(zip(np.arange(-180,177.01),np.arange(-179,178.01))))
	ds_compare = xr.Dataset({'lat_bnds':(['lat','bnds'],lat_bnds),
	                         'lon_bnds':(['lon','bnds'],lon_bnds)},
	                        coords={'lat':(['lat'],np.arange(-89.5,89.51)),
	                                'lon':(['lon'],np.arange(-179.5,177.51))})

	xr.testing.assert_allclose(ds,ds_compare)

def test_get_bnds_offsetgrid():
	# Full planetary grid, but with pixels crossing the antimeridian
	ds = xr.Dataset(coords={'lat':(['lat'],np.arange(-89.4,89.7)), #  -180:180, offset
                        	'lon':(['lon'],np.arange(-179.4,179.7))})

	ds = get_bnds(ds)

	lat_bnds = np.array(list(zip(np.arange(-89.9,89.11),np.arange(-88.9,90.11))))
	lat_bnds[-1,-1] = 90 # Truncate at 90
	lon_bnds = np.array(list(zip(np.arange(-179.9,179.11),np.arange(-178.9,180.11))))
	lon_bnds[-1,-1] = -179.9 # Wrap around across antimeridian
	ds_compare = xr.Dataset({'lat_bnds':(['lat','bnds'],lat_bnds),
	                         'lon_bnds':(['lon','bnds'],lon_bnds)},
	                        coords={'lat':(['lat'],np.arange(-89.4,89.7)),
	                                'lon':(['lon'],np.arange(-179.4,179.7))})

	xr.testing.assert_allclose(ds,ds_compare)

def test_get_bnds_1pixelinEH():
	# Crossing the antimeridian, with just one pixel that is on the
	# other side of the antimeridian, EH version (this can be an 
	# issue with the code)
	ds = xr.Dataset(coords={'lat':(['lat'],np.array([0,1])),
	                        'lon':(['lon'],np.array([-179.8, -178.8,179.2]))})

	ds = get_bnds(ds)

	lat_bnds = np.array(list(zip(np.arange(-0.5,0.51),np.arange(0.5,1.51))))
	lon_bnds = np.array([[179.7,-179.3],
	                      [-179.3,-178.3],
	                      [178.7,179.7]])
	ds_compare = xr.Dataset({'lat_bnds':(['lat','bnds'],lat_bnds),
	                         'lon_bnds':(['lon','bnds'],lon_bnds)},
	                        coords={'lat':(['lat'],np.array([0,1])),
	                        'lon':(['lon'],np.array([-179.8, -178.8,179.2]))})

	xr.testing.assert_allclose(ds,ds_compare)



def test_get_bnds_1pixelinWH():
	# Crossing the antimeridian, with just one pixel that is on the
	# other side of the antimeridian, WH version (this can be an 
	# issue with the code)
	ds = xr.Dataset(coords={'lat':(['lat'],np.array([0,1])),
	                        'lon':(['lon'],np.array([-179.8, 178.2,179.2]))})

	ds = get_bnds(ds)

	lat_bnds = np.array(list(zip(np.arange(-0.5,0.51),np.arange(0.5,1.51))))
	lon_bnds = np.array(list(zip(np.arange(176.7,178.71),np.arange(177.7,179.71))))
	lon_bnds[0,0] = 179.7
	lon_bnds[0,1] = -179.3
	ds_compare = xr.Dataset({'lat_bnds':(['lat','bnds'],lat_bnds),
	                         'lon_bnds':(['lon','bnds'],lon_bnds)},
	                        coords={'lat':(['lat'],np.array([0,1])),
	                        'lon':(['lon'],np.array([-179.8, 178.2,179.2]))})

	xr.testing.assert_allclose(ds,ds_compare)




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




