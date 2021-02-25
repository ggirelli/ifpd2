"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import logging
import sys
from typing import Callable


def ert_type(x, stype, label):
    assert isinstance(x, stype), f"{label} should be {stype}, {type(x)} instead"


def ert_multiTypes(x, types, label):
    cond = any([isinstance(x, t) for t in types])
    assert cond, f"{label} should be one of {types}, {type(x)} instead"


def ert_nonNeg(x, label, include_zero=False):
    if not include_zero:
        assert x > 0, f"{label} should be greater than 0"
    else:
        assert x >= 0, f"{label} should be greater than or equal to 0"


def ert_inInterv(x, vmin, vmax, label, leftClose=False, rightClose=True):
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


def enable_rich_assert(fun: Callable) -> Callable:
    def wrapper(*args, **kwargs):
        try:
            return fun(*args, **kwargs)
        except AssertionError as e:
            logging.exception(e)
            sys.exit()

    return wrapper
