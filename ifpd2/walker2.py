"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from ifpd2.database import DataBase
from typing import Tuple


class GenomicRegion(object):
    """docstring for GenomicRegion"""

    __chrom: bytes
    __chromStart: int
    __chromEnd: int

    def __init__(self, chrom: bytes, chromStart: int, chromEnd: int, focus: float = 1):
        super(GenomicRegion, self).__init__()
        assert 0 != len(chrom), "chromosome cannot be empty"
        assert (
            chromStart >= 0
        ), f"start should be greater than or equal to 0: {chromStart}"
        assert (
            chromEnd > chromStart
        ), f"end should be greater than start: {chromStart}-{chromEnd}"
        assert (
            focus > 0 and focus <= 1
        ), f"focus expected to be in the (0,1] range: {focus}"
        self.__chrom = chrom
        self.__chromStart = chromStart
        self.__chromEnd = chromEnd
        self.__focusStart = int(
            (chromStart + chromEnd) / 2 - (chromEnd - chromStart) * focus
        )
        self.__focusEnd = int(
            (chromStart + chromEnd) / 2 + (chromEnd - chromStart) * focus
        )

    @property
    def chromosome(self) -> bytes:
        return self.__chrom

    @property
    def start(self) -> int:
        return self.__chromStart

    @property
    def end(self) -> int:
        return self.__chromEnd

    @property
    def focus(self) -> Tuple[int, int]:
        return (self.__focusStart, self.__focusEnd)


class Walker(object):
    """docstring for Walker"""

    def __init__(self):
        super(Walker, self).__init__()

    def walk_region(self, db: DataBase, region: GenomicRegion):
        db.buffer(region.chromosome, region.start, region.end)
