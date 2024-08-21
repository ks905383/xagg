import pytest
import pandas as pd
import numpy as np
import xarray as xr
import geopandas as gpd
from matplotlib import pyplot as plt

from xagg.core import aggregate
from xagg.wrappers import pixel_overlaps

try:
    from cartopy import crs as ccrs
    from matplotlib import pyplot as pyplot
    import cmocean
    _has_plotpckgs=True
except ImportError:
	# To be able to test the rest with environments without xesmf
	_has_plotpckgs=False

##### diag_fig() tests #####

# Load sample data
ds = xr.open_dataset('data/climate_data/tas_Amon_CCSM4_rcp85_monthavg_20700101-20991231.nc')
gdf = gpd.read_file('data/geo_data/UScounties.shp')
# Subset manually to not go through slow subset_find in every test
ds = ds.sel(lat=slice(30,55),lon=slice(360-130,360-65))

# Calculate overlaps
wm = pixel_overlaps(ds,gdf)

if _has_plotpckgs:
	def test_diag_fig_noerror():
		# Test whether error occurs when calling diag_fig()
		# With some random location
		wm.diag_fig(50,ds)

	def test_diag_fig_isfig():
		# Test to make sure a figure, axis is returned
		fig,ax = wm.diag_fig(50,ds)

		assert isinstance(fig,plt.Figure)
		assert isinstance(ax,plt.Axes)

	def test_diag_fig_subsetbypolyid():
		# Test to make sure the right location is returned 
		# when using an integer poly index
		fig,ax = wm.diag_fig(50,ds)

		assert ax.get_title() == 'Poly #50: Sanders; Montana; 30; 089; 30089'

	def test_diag_fig_subsetbyname():
		# Test to make sure the right location is returned
		# when using a column dictionary
		fig,ax=wm.diag_fig({'NAME':'Los Angeles'},ds)

		assert ax.get_title() == 'Poly #2384: Los Angeles; California; 06; 037; 06037'
else:
	def test_diag_fig_noimport():
		# Should raise ImportError in the no-plot environment
		with pytest.raises(ImportError):
			wm.diag_fig(50,ds)
