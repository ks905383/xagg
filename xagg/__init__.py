# Eventually restrict to just pixel_overlaps and aggregate; with 
# everything else happening behind the scenes (and the exporting 
# happening as methods to the classes that are exported from those
# two functions)
from .wrappers import pixel_overlaps
from .aux import (normalize,fix_ds,get_bnds,subset_find)
from .core import aggregate