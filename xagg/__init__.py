# Eventually restrict to just pixel_overlaps and aggregate; with 
# everything else happening behind the scenes (and the exporting 
# happening as methods to the classes that are exported from those
# two functions)
from .wrappers import pixel_overlaps
from .auxfuncs import (normalize,fix_ds,get_bnds,subset_find)
from .core import (aggregate,read_wm)