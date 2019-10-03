
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for quick assertions.
'''

def assert_type(x, stype, label):
	assert isinstance(x,stype), f"{label} should be {stype}, {type(x)} instead"

def assert_multiTypes(x, types, label):
	cond = any([isinstance(x, t) for t in types])
	assert cond, f"{label} should be one of {types}, {type(x)} instead"

def assert_nonNeg(x, label):
	assert x > 0, f"{label} should be greater than 0"

def assert_inInterv(x, vmin, vmax, label, leftClose = False, rightClose = True):
	if leftClose:
		if rightClose:
			assert x >= vmin and x <= vmax, f"expected {vmin}<={label}<={vmax}"
		else:
			assert x >= vmin and x < vmax, f"expected {vmin}<={label}<{vmax}"
	else:
		if rightClose:
			assert x > vmin and x <= vmax, f"expected {vmin}<{label}<={vmax}"
		else:
			assert x > vmin and x < vmax, f"expected {vmin}<{label}<{vmax}"
