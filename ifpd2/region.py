"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from ifpd2.chromosome import ChromosomeData
import logging
from typing import List, Tuple


class GenomicRegion(object):
    """Details on genomic region with a central focus region."""

    __chrom: bytes
    __chromStart: int
    __chromEnd: int
    __focusStart: int
    __focusEnd: int
    __focusSize: int
    __focusStep: int

    def __init__(
        self,
        chrom_region: Tuple[bytes, int, int],
        focus_style: Tuple[float, float] = (1, 1),
    ):
        super(GenomicRegion, self).__init__()
        self.__init_region(*chrom_region)
        self.__init_focus(*focus_style)

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


class GenomicRegionBuilder(object):
    """docstring for GenomicRegionBuilder"""

    __chromosome: bytes
    __chromosome_size_nt: int
    __focus_style: Tuple[float, float]

    def __init__(
        self, chromosome_data: ChromosomeData, focus: float = 1, focus_step: float = 1
    ):
        super(GenomicRegionBuilder, self).__init__()
        self.__chromosome = chromosome_data.name
        self.__chromosome_size_nt = chromosome_data.size_nt
        assert focus > 0
        assert focus_step > 0
        self.__focus_style = (focus, focus_step)

    def __build_overlapping(self, size: int, step: int) -> List[List[GenomicRegion]]:
        assert step < size
        region_set_list: List[List[GenomicRegion]] = []
        for range_start in range(0, size, step):
            genomic_region_set: List[GenomicRegion] = []
            for start in range(range_start, self.__chromosome_size_nt, step):
                end = start + size
                if end < self.__chromosome_size_nt:
                    genomic_region_set.append(
                        GenomicRegion(
                            (self.__chromosome, start, end), self.__focus_style
                        )
                    )
            region_set_list.append(genomic_region_set)
        return region_set_list

    def __build_non_overlapping(
        self, size: int, step: int
    ) -> List[List[GenomicRegion]]:
        assert step == size
        genomic_region_set: List[GenomicRegion] = []
        for start in range(0, self.__chromosome_size_nt, step):
            end = start + size
            if end <= self.__chromosome_size_nt:
                genomic_region_set.append(
                    GenomicRegion((self.__chromosome, start, end), self.__focus_style)
                )
        return [genomic_region_set]

    def build_by_number(self, n: int) -> List[List[GenomicRegion]]:
        step: int = self.__chromosome_size_nt // n
        return self.build_by_size(step, step)

    def build_by_size(self, size: int, step_style: float) -> List[List[GenomicRegion]]:
        if step_style > 1:
            step = int(step_style)
        else:
            step = int(size * step_style)
        if step < size:
            return self.__build_overlapping(size, step)
        else:
            return self.__build_non_overlapping(size, step)
