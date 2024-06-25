# Partially adapted from xarray's xr.core.options
# NB: In a future version, setting `impl` and `silent`
# in individual function calls should be deprecated 
# in favor of setting global defaults or using 
# with blocks. 
from typing import TypedDict

# Create class specifying the needed
# type of each option
class T_Options(TypedDict):
	silent : bool
	impl : str
	rgrd_alg : str
	nan_to_zero_regridding : bool


# Set default options. Defining it in
# the module makes it global to the module
# (as opposed to within a function, where
# it should only be local to the function)
OPTIONS: T_Options = {
	'silent' : False,
	'impl' : 'for_loop',
	'rgrd_alg' : 'conservative',
	'nan_to_zero_regridding' : True
}


# Options for the backend implementation
_IMPL_OPTIONS = frozenset(['for_loop','dot_product'])
# Options for regridding
_RGRD_OPTIONS = frozenset(['bilinear','conservative'])

# Each item of this dictionary is a test for whether a
# modification for the corresponding option was correctly
# set. I.e., "silent" can only be True or False, so the
# 'silent' dict option here tests for whether it's a bool. 
_VALIDATORS = {
	'silent': lambda value: isinstance(value, bool),
	'impl': _IMPL_OPTIONS.__contains__,
	'rgrd_alg': _RGRD_OPTIONS.__contains__,
	'nan_to_zero_regridding': lambda value: isinstance(value,bool)
}

# Define options class
class set_options: 
	""" Set options for xagg.

	Parameters
	----------
	silent : bool, by default ``False``
		If True, then status updates are suppressed 

	impl : str, by default ``"for_loop"``
		Sets backend algorithm, can be 

		* ``for_loop``: slower, but lower memory use
		* ``dot_product``: faster, but higher memory use

	"""
	
	def __init__(self, **kwargs):
		# Keep track of changed options, to be able to change
		# them back if used in a `with` block (see __exit__ below)
		self.old = {}

		for k, v in kwargs.items():
			# Check to make sure the option you're looking to change
			# is an option changeable by set_options()
			if k not in OPTIONS:
				raise ValueError(
					f"argument name {k!r} is not in the set of valid options {set(OPTIONS)!r}"
				)

			# Check to make sure the new value of the option is 
			# acceptable
			if k in _VALIDATORS and not _VALIDATORS[k](v):
				if k == "impl":
					expected = f"Expected one of {_IMPL_OPTIONS!r}"
				else:
					expected = ""
				raise ValueError(
					f"option {k!r} given an invalid value: {v!r}. " + expected
				)
			# Note original value of changed options, to reset
			# defaults in the __exit__ block below
			self.old[k] = OPTIONS[k]
		# Update OPTIONS 
		self._apply_update(kwargs)

	def _apply_update(self, options_dict):
		""" Update OPTIONS """
		OPTIONS.update(options_dict)

	# This allows use with "with set_options(...):"
	# See e.g., https://stackoverflow.com/questions/1984325/explaining-pythons-enter-and-exit
	def __enter__(self):
		return
	# Resets to the original OPTIONS at the end 
	# of a `with` block
	def __exit__(self, type, value, traceback):
		self._apply_update(self.old)


def get_options():
	"""
	Get module-wide options for xagg.

	See Also
	----------
	:py:meth:`set_options`

	"""
	# Returning a copy, to make sure no unintended changes to
	# the output of get_options() trickle down to anything
	# else that uses the variable
	return OPTIONS.copy()

	