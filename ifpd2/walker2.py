"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import logging
import pandas as pd  # type: ignore
from typing import Tuple, Union


class GenomicWindowSet(object):
    _window_data: pd.DataFrame = pd.DataFrame(
        columns=[
            "set_id",
            "window_id",
            "start",
            "end",
            "focus_start",
            "focus_end",
        ]
    )
    _window_id: int = -1

    feature: str  # Chromosome
    region: Tuple[int, int]

    n_probes: int
    window_size: int
    window_shift: float

    # Either in nt (x>1) or fraction of window_size (0<x<=1)
    focus_size: Union[int, float]
    # Either in nt (x>1) or fraction of focus_size (0<x<=1)
    focus_step: Union[int, float]

    _growing: bool = False

    def __init__(self, region: Tuple[int, int]):
        super(GenomicWindowSet, self).__init__()
        self.region = region
        assert self.region[0] <= self.region[1]

    @property
    def growing(self) -> bool:
        return self._growing

    @property
    def reached_last_window(self) -> bool:
        return (self._window_id + 1) == self._window_data.shape[0]

    def init_windows_for_n_probes(self, n_probes: int) -> None:
        if 0 < self._window_data.shape[0]:
            logging.warning("cannot re-initalize windows.")
            return

    def init_windows_by_size(self, size: int, shift: float) -> None:
        if 0 < self._window_data.shape[0]:
            logging.warning("cannot re-initalize windows.")
            return
