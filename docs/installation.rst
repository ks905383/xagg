Installation
=======================================
The easiest and recommended way to install :py:mod:`xagg` is through ``conda`` or ``mamba``::

   # Mamba
   mamba install -c conda-forge xagg==0.3.2.0

   # Conda
   conda install -c conda-forge xagg==0.3.2.0


``xagg`` can also be installed through ``pip``::

   pip install xagg

though certain depencies may not be available. 

Optional dependencies
----------------------------------------
If using `weights` grids to add another weight layer to the aggregation (e.g., weighting raster data additionally by population density), we recommend installing :py:mod:`xesmf`, which is required for :py:mod:`xagg` to perform regridding if the weight and raster grids are not equal. :py:mod:`xesmf` must be [installed manually](https://xesmf.readthedocs.io/en/stable/installation.html), since its dependencies are not available through ``pip`` (and ``conda`` does not support installing optional dependencies). 


If wanting to create diagnostic figures using :py:meth:`weightmap.diag_fig()`, :py:mod:`matplotlib`, :py:mod:`cartopy`, and :py:mod:`cmocean` are additionally required. These can be installed using ``pip`` ::

   pip install xagg[plots]


