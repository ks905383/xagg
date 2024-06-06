import pytest

from xagg.options import set_options,get_options


##### set_options() tests #####
def test_set_options_withtemp():
	# Test to make sure using set_options in a 
	# with block doesn't change the original 
	# options
	# (Hardcoded... would have to change if the 
	# settings default is changed...)
	with set_options(silent=True,impl='dot_product'):
		pass
	assert not get_options()['silent']
	assert get_options()['impl'] == 'for_loop'

def test_set_options():
	# Test changing options works
	with set_options(silent=True):
		assert get_options()['silent']
	with set_options(impl='dot_product'):
		assert get_options()['impl'] == 'dot_product'

def test_set_options_badoption():
	# Test error for unsupported option name
	with pytest.raises(ValueError):
		set_options(fake_option=True)

def test_set_options_badoptioninput():
	# Test error for unsupported option input
	# Testing with silent = not a bool 
	# (seems too much to test for every option...
	# but maybe required?)
	with pytest.raises(ValueError):
		set_options(silent='a')
	with pytest.raises(ValueError):
		set_options(impl='fake_option')