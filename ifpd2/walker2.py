"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

from ifpd2.database import DataBase, Record
import itertools
import logging
import os
from typing import IO, Iterator, List, Tuple, Union


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


class ChromosomeWalker(object):
    """docstring for Walker"""

    __chromosome: bytes
    __db: DataBase
    __IH: IO

    def __init__(self, db: Union[DataBase, str], chromosome: bytes):
        super(ChromosomeWalker, self).__init__()
        if isinstance(db, str):
            self.__db = DataBase(db)
        else:
            self.__db = db
        assert chromosome in self.__db._chromosomes.keys()
        self.__IH = open(
            os.path.join(self.__db._root, f"{chromosome.decode()}.bin"), "rb"
        )
        self.__chromosome = chromosome

    @property
    def db(self) -> DataBase:
        return self.__db

    def __jump_to_bin(self, start_from_nt: int = 0) -> None:
        """Move buffer pointer to the first record of a bin.

        Keyword Arguments:
            start_from_nt {int} -- position in nucleotides (default: {0})
        """
        position_in_bytes = self.__db._chromosomes.get_chromosome(
            self.__chromosome
        ).index[start_from_nt]
        self.__IH.seek(position_in_bytes)

    def read_next_record(self) -> bytes:
        """Read the next record.

        The buffer pointer moves to the start of the following record.

        Returns:
            bytes -- record bytes
        """
        return self.__IH.read(self.__db.record_byte_size)

    def read_previous_record(self) -> bytes:
        """Read the previous record.

        The buffer pointer does not move.

        Returns:
            bytes -- record bytes
        """
        self.rewind()
        return self.read_next_record()

    def read_next_record_and_rewind(self) -> bytes:
        """Reads the next record.

        The buffer pointer does not move.

        Returns:
            bytes -- record bytes
        """
        record = self.read_next_record()
        self.rewind()
        return record

    def next_record_exists(self) -> bool:
        """Tries reading the next record.

        Returns:
            bool -- whether next record exists
        """
        if 0 == len(self.read_next_record_and_rewind()):
            return False
        return True

    def rewind(self) -> None:
        """Rewind of one record.

        The buffer pointer moves to the beginning of the previous record.
        """
        self.__IH.seek(max(self.__IH.tell() - self.__db.record_byte_size, 0))

    def skip(self) -> None:
        """Skip one record.

        The buffer pointer moves to the beginning of the following record.
        """
        self.__IH.seek(self.__IH.tell() + self.__db.record_byte_size)

    def fastforward(self, start_from_nt: int) -> None:
        """Jump to the first record at a given position.

        First the buffer pointer is moved to the beginning of the bin containing the
        start_from_nt position. Then, records are read an parsed until their start is
        greater than the start_from_nt value. Finally, the buffer pointer is moved to
        the beginning of the last parsed record.

        Arguments:
            IH {IO} -- database input handle
            chromosome {bytes} -- chromosome label
            start_from_nt {int} -- position to fastforward to (in nucleotides)
        """
        if 0 == start_from_nt:
            self.__IH.seek(0)
            return
        self.__jump_to_bin(start_from_nt)

        record_start = 0
        while record_start < start_from_nt:
            if self.next_record_exists():
                record_start = Record(self.read_next_record(), self.__db.dtype)["start"]
            else:
                logging.warning(
                    " ".join(
                        [
                            "the specified location",
                            f"({self.__chromosome.decode()}:{start_from_nt})",
                            "is outside the database.",
                        ]
                    )
                )
                return

    def buffer(self, start_from_nt: int = 0, end_at_nt: int = -1) -> Iterator[Record]:
        """Buffer a chromosome's records.

        Buffer records from a chromosome withing the specified region.
        To buffer the whole records, specify a [0, -1] region.

        Keyword Arguments:
            start_from_nt {int} -- region starting position (default: {0})
            end_at_nt {int} -- region end position (default: {-1})

        Yields:
            Iterator[Record] -- parsed record
        """
        self.fastforward(start_from_nt)
        record = self.read_next_record()
        while 0 != len(record):
            parsed_record = Record(record, self.__db.dtype)
            if parsed_record["start"] > end_at_nt and end_at_nt > 0:
                break
            yield parsed_record
            record = self.read_next_record()

    def walk_region(self, region: GenomicRegion) -> Iterator[List[Record]]:
        assert region.chromosome == self.__chromosome
        focus_start, focus_end = region.focus
        record_list = [r for r in self.buffer(focus_start, focus_end)]
        yield record_list
        while region.can_increase_focus():
            region.increase_focus()
            new_focus_start, new_focus_end = region.focus
            record_list = list(
                itertools.chain(
                    *[
                        [r for r in self.buffer(new_focus_start, focus_start)],
                        record_list,
                        [r for r in self.buffer(focus_end, new_focus_end)],
                    ]
                )
            )
            yield record_list
            focus_start, focus_end = region.focus
