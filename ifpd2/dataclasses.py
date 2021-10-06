"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from dataclasses import dataclass
from os.path import isdir, isfile
from typing import Optional, Tuple


def path_exists(path: str) -> bool:
    return not isdir(path) and not isfile(path)


@dataclass(frozen=True)
class Folder:
    path: str
    exists: bool = False

    def __post_init__(self):
        if self.exists:
            assert isdir(self.path)
        else:
            assert not path_exists(self.path)


@dataclass(frozen=True)
class File:
    path: str
    exists: bool = False

    def __post_init__(self):
        if self.exists:
            assert isfile(self.path)
        else:
            assert not path_exists(self.path)


@dataclass(frozen=True)
class PositiveInteger:
    n: int

    def __post_init__(self):
        assert self.n >= 1


@dataclass(frozen=True)
class NonNegativeFloat:
    n: float

    def __post_init__(self):
        assert self.n >= 0


@dataclass(frozen=True)
class NonNegativeInteger:
    n: int

    def __post_init__(self):
        assert self.n > 0


@dataclass(frozen=True)
class PositiveFloat:
    n: float
    limit: Optional[float] = None
    limit_included: bool = True

    def __post_init__(self):
        if self.limit is not None:
            if self.limit_included:
                assert 0 < self.n <= self.limit
            else:
                assert 0 < self.n < self.limit
        else:
            assert 0 < self.n


@dataclass(frozen=True)
class GenomicRegion:
    start: int
    end: int

    def __post_init__(self):
        assert self.start >= 0
        assert self.end >= self.start or self.end == -1

    def astuple(self) -> Tuple[int, int]:
        return (self.start, self.end)


@dataclass(frozen=True)
class NonNegativeIntInterval:
    start: int
    end: int

    def __post_init__(self):
        assert self._from >= 0
        assert self.end >= self.first

    def astuple(self) -> Tuple[int, int]:
        return (self.start, self.end)


@dataclass(frozen=True)
class QueryWindow:
    size: Optional[int]
    shift: Optional[float]

    def __post_init__(self):
        if self.size is not None:
            assert self.size >= 1
        assert 0 < self.shift <= 1

    def astuple(self) -> Tuple[Optional[int], Optional[float]]:
        return (self.size, self.shift)


@dataclass(frozen=True)
class QueryFocus:
    size: float
    step: float

    def __post_init__(self):
        assert self.size > 0
        assert self.step > 0


@dataclass(frozen=True)
class FreeEnergyInterval:
    start: float
    end: float

    def astuple(self) -> Tuple[float, float]:
        return (self.start, self.end)
