"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from ifpd2.database import DataBase, Record
import itertools
import logging
from typing import Iterator, List, Tuple


class GenomicRegion(object):
    """docstring for GenomicRegion"""

    __chrom: bytes
    __chromStart: int
    __chromEnd: int
    __focusStart: int
    __focusEnd: int
    __focusSize: int
    __focusStep: int

    def __init__(
        self,
        chrom: bytes,
        chromStart: int,
        chromEnd: int,
        focus: float = 1,
        focus_step: float = 1,
    ):
        super(GenomicRegion, self).__init__()
        self.__init_region(chrom, chromStart, chromEnd)
        self.__init_focus(focus, focus_step)

    def __init_region(self, chrom: bytes, chromStart: int, chromEnd: int) -> None:
        assert 0 != len(chrom), "chromosome cannot be empty"
        assert (
            chromStart >= 0
        ), f"start should be greater than or equal to 0: {chromStart}"
        assert (
            chromEnd > chromStart
        ), f"end should be greater than start: {chromStart}-{chromEnd}"
        self.__chrom = chrom
        self.__chromStart = chromStart
        self.__chromEnd = chromEnd

    def __init_focus(self, focus: float, focus_step: float) -> None:
        if focus > 1:
            assert focus <= self.__chromEnd - self.__chromStart
            self.__focusSize = int(focus)
        else:
            self.__focusSize = int((self.__chromEnd - self.__chromStart) * focus)
        self.__focusStart = int(
            (self.__chromStart + self.__chromEnd) / 2 - self.__focusSize / 2
        )
        self.__focusEnd = int(
            (self.__chromStart + self.__chromEnd) / 2 + self.__focusSize / 2
        )
        assert focus_step > 0
        if focus_step > 1:
            self.__focusStep = int(focus_step)
        else:
            self.__focusStep = int(self.__focusSize * focus_step)

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
    def region(self) -> Tuple[int, int]:
        return (self.__chromStart, self.__chromEnd)

    @property
    def focus(self) -> Tuple[int, int]:
        return (self.__focusStart, self.__focusEnd)

    @property
    def focus_step(self) -> int:
        return self.__focusStep

    def can_increase_focus(self) -> bool:
        return self.focus != self.region

    def increase_focus(self) -> None:
        if self.can_increase_focus():
            logging.warning("cannot increase the focus region any further")
        self.__focusStart = max(
            self.__chromStart, self.__focusStart - self.__focusStep // 2
        )
        self.__focusEnd = min(self.__chromEnd, self.__focusEnd + self.__focusStep // 2)


class Walker(object):
    """docstring for Walker"""

    def __init__(self):
        super(Walker, self).__init__()

    def walk_region(
        self, db: DataBase, region: GenomicRegion
    ) -> Iterator[List[Record]]:
        focus_start, focus_end = region.focus
        record_list = [r for r in db.buffer(region.chromosome, focus_start, focus_end)]
        yield record_list
        while region.can_increase_focus():
            region.increase_focus()
            new_focus_start, new_focus_end = region.focus
            record_list = list(
                itertools.chain(
                    *[
                        [
                            r
                            for r in db.buffer(
                                region.chromosome, new_focus_start, focus_start
                            )
                        ],
                        record_list,
                        [
                            r
                            for r in db.buffer(
                                region.chromosome, focus_end, new_focus_end
                            )
                        ],
                    ]
                )
            )
            yield record_list
            focus_start, focus_end = region.focus
