"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
import copy
from ifpd2 import const
import logging
import numpy as np  # type: ignore
import os
import pickle
import pandas as pd  # type: ignore
from rich.progress import Progress, TaskID  # type: ignore
from typing import Any, Dict, IO, Iterator, List, Set, Tuple


def get_dtype_length(dtype) -> int:
    return sum([int(label.strip("><|SUuif")) for label in dtype.values()])


class ChromosomeIndex(object):
    """ChromosomeIndex"""

    _bin_size: int
    _index: Dict[int, Tuple[int, int]]

    def __init__(self, bin_size: int):
        super(ChromosomeIndex, self).__init__()
        assert bin_size >= 1
        self._bin_size = bin_size

    def __init_index(self, chrom_db: pd.DataFrame) -> None:
        self._index = {}
        chrom_size_nt = chrom_db["end"].values.max()
        for bin_id in range(0, (chrom_size_nt // self._bin_size) + 1):
            self._index[bin_id] = (np.inf, 0)

    def __populate_bins(
        self,
        chrom_db: pd.DataFrame,
        record_byte_size: int,
        track: Tuple[Progress, TaskID],
    ) -> None:
        current_position = -1
        for i in range(chrom_db.shape[0]):
            track[0].update(track[1], advance=1)
            position_in_nt = chrom_db["start"].values[i]
            assert position_in_nt > current_position
            current_position = position_in_nt

            position_in_bytes = record_byte_size * i
            binned_to = position_in_nt // self._bin_size

            bin_start, bin_end = list(self._index[binned_to])
            if bin_start > position_in_bytes:
                bin_start = position_in_bytes
            if bin_end < position_in_bytes:
                bin_end = position_in_bytes
            self._index[binned_to] = (bin_start, bin_end)

    def __fill_empty_bins(self) -> None:
        if not np.isfinite(self._index[0][0]):
            self._index[0] = (0, 0)
        for bin_id, (start, end) in self._index.items():
            if not np.isfinite(start):
                self._index[bin_id] = (
                    self._index[bin_id - 1][1],
                    self._index[bin_id - 1][1],
                )

    def build(
        self,
        chrom_db: pd.DataFrame,
        record_byte_size: int,
        track: Tuple[Progress, TaskID],
    ) -> None:
        for colname in ("chromosome", "start", "end"):
            assert colname in chrom_db.columns, f"missing '{colname}' column"
        chromosome_set: Set[bytes] = set(chrom_db["chromosome"].values)
        assert 1 == len(chromosome_set)

        self.__init_index(chrom_db)
        self.__populate_bins(chrom_db, record_byte_size, track)
        self.__fill_empty_bins()

    def __getitem__(self, position_in_nt: int) -> int:
        assert self._index is not None
        binned_to = position_in_nt // self._bin_size
        if binned_to not in self._index:
            return -1

        position_in_bytes = self._index[binned_to][0]
        if not np.isfinite(position_in_bytes):
            return -1
        return position_in_bytes


class ChromosomeData(object):
    """Contains information on a chromosome"""

    _name: str
    _size_nt: int
    _size_bytes: int
    _recordno: int
    _index: ChromosomeIndex
    _record_byte_size: int

    def __init__(
        self,
        chromosome_db: pd.DataFrame,
        dtype: Dict[str, str],
        index_bin_size: int,
        progress: Progress,
    ):
        super(ChromosomeData, self).__init__()

        assert "chromosome" in chromosome_db.columns
        selected_chrom = chromosome_db["chromosome"][0]
        assert 1 == len(set(chromosome_db["chromosome"].values))

        self._record_byte_size = get_dtype_length(dtype)
        assert self._record_byte_size > 0

        self._name = selected_chrom
        self._recordno = chromosome_db.shape[0]
        self._size_nt = chromosome_db["end"].values.max()
        self._size_bytes = chromosome_db.shape[0] * self._record_byte_size

        self._build_index(chromosome_db, index_bin_size, progress)

    @property
    def record_byte_size(self) -> int:
        return self._record_byte_size

    @property
    def size_nt(self):
        return self._size_nt

    @property
    def size_bytes(self):
        return self._size_bytes

    @property
    def recordno(self):
        return self._recordno

    @property
    def index(self):
        return copy.copy(self._index)

    def _build_index(
        self, chromosome_db: pd.DataFrame, index_bin_size: int, progress: Progress
    ) -> None:
        assert index_bin_size > 0
        indexing_track = progress.add_task(
            f"indexing {self._name}.bin",
            total=chromosome_db.shape[0],
            transient=True,
        )
        self._index = ChromosomeIndex(index_bin_size)
        self._index.build(
            chromosome_db, self._record_byte_size, (progress, indexing_track)
        )
        progress.remove_task(indexing_track)


class ChromosomeDict(object):
    """Wraps all chromosomes"""

    _index_bin_size: int
    _data: Dict[bytes, ChromosomeData]

    def __init__(self, index_bin_size: int = const.DEFAULT_DATABASE_INDEX_BIN_SIZE):
        super(ChromosomeDict, self).__init__()
        self._data = {}
        assert index_bin_size > 0
        self._index_bin_size = index_bin_size

    def __len__(self) -> int:
        return len(self._data)

    def keys(self) -> List[bytes]:
        return list(self._data.keys())

    def items(self) -> List[Tuple[bytes, ChromosomeData]]:
        return list(self._data.items())

    @property
    def sizes_nt(self) -> Dict[bytes, int]:
        return dict([(name, data.size_nt) for name, data in self._data.items()])

    @property
    def sizes_bytes(self) -> Dict[bytes, int]:
        return dict([(name, data.size_bytes) for name, data in self._data.items()])

    @property
    def recordnos(self) -> Dict[bytes, int]:
        return dict([(name, data.recordno) for name, data in self._data.items()])

    def get_chromosome(self, chromosome: bytes) -> ChromosomeData:
        return copy.copy(self._data[chromosome])

    def add_chromosome(
        self, chromosome_db: pd.DataFrame, dtype: Dict[str, str], progress: Progress
    ) -> None:
        self._data[chromosome_db["chromosome"][0]] = ChromosomeData(
            chromosome_db, dtype, self._index_bin_size, progress
        )


class Record(object):
    """DataBase Record"""

    _data: Dict[str, Any]

    def __init__(self, record: bytes, column_dtypes: Dict[str, str]):
        super(Record, self).__init__()
        self.__parse_bytes(record, column_dtypes)

    @property
    def data(self) -> Dict[str, Any]:
        return copy.copy(self._data)

    def __parse_bytes(self, record: bytes, column_dtypes: Dict[str, str]) -> None:
        assert len(record) == get_dtype_length(column_dtypes)

        self._data = {}
        current_location = 0
        for label in const.database_columns:
            dtype = column_dtypes[label]
            byte_size = get_dtype_length({label: dtype})

            self._data[label] = np.frombuffer(
                record[current_location : (current_location + byte_size)], dtype
            )[0]

            current_location += byte_size

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self._data, index=[0])

    def to_csv(self, sep: str = ",") -> str:
        csv_fields: List[str] = []
        for x in const.database_columns:
            data = self._data[x]
            if isinstance(data, bytes):
                csv_fields.append(str(data.decode()))
            else:
                csv_fields.append(str(data))
        return sep.join(csv_fields)

    @staticmethod
    def __norm_value_in_range(v, r):
        if 0 == r[1]:
            return np.nan
        return (v - r[0]) / (r[1] - r[0])

    def __update_score_by_nOT(self, F):
        off_target_no = self._data["off_target_no"]
        if off_target_no <= F[0]:
            self._data["score"] = 0
            return
        if off_target_no > F[1]:
            self._data["score"] = np.inf
            return
        return self.__norm_value_in_range(off_target_no, F)

    def __update_score_by_dG_Tm(self, Gs):
        ss_dG = self._data["ss_dG"]
        tm_dG = self._data["Tm_dG"]
        if ss_dG >= tm_dG * min(Gs):
            self._data["score"] = 0
            return
        if ss_dG < tm_dG * max(Gs):
            self._data["score"] = np.inf
            return
        return self.__norm_value_in_range(ss_dG, [tm_dG * f for f in Gs])

    def __update_score_by_dG_Gs(self, Gs):
        ss_dG = self._data["ss_dG"]
        if ss_dG >= Gs[0]:
            self._data["score"] = 0
            return
        if ss_dG < Gs[1]:
            self._data["score"] = np.inf
            return
        return self.__norm_value_in_range(ss_dG, Gs)

    def add_score(self, F, Gs) -> None:
        norm_nOT = self.__update_score_by_nOT(F)
        if norm_nOT is None:
            return
        if all([x >= 0 for x in Gs]):
            norm_ss_dG = self.__update_score_by_dG_Tm(Gs)
        else:
            norm_ss_dG = self.__update_score_by_dG_Gs(Gs)
        if norm_ss_dG is None:
            return
        self._data["score"] = np.mean([norm_nOT, norm_ss_dG])

    def __getitem__(self, key: str) -> Any:
        allowed_columns = ["score"]
        allowed_columns.extend(const.database_columns)
        if key not in allowed_columns:
            raise KeyError(f"unrecognized key '{key}'")
        else:
            return self.data[key]

    def __repr__(self) -> str:
        return str(self._data)


class DataBase(object):
    """Buffering and checking class for ifpd2 database."""

    _root: str
    _args: argparse.Namespace
    _chromosomes: ChromosomeDict
    _dtype: Dict[str, str]
    _record_byte_size: int

    def __init__(self, path: str):
        super(DataBase, self).__init__()
        assert os.path.isdir(path), f"cannot find database folder '{path}'"
        db_pickle_path = os.path.join(path, "db.pickle")
        assert os.path.isfile(db_pickle_path), f"'db.pickle is missing in '{path}'"

        with open(db_pickle_path, "rb") as IH:
            details = pickle.load(IH)
            self._chromosomes = details["chromosomes"]
            self._dtype = details["dtype"]
            self._args = details["args"]

        assert isinstance(self._chromosomes, ChromosomeDict)
        for chromosome in self.chromosome_list:
            chromosome_path = os.path.join(path, f"{chromosome.decode()}.bin")
            assert os.path.isfile(
                chromosome_path
            ), f"missing expected chromosome file: '{chromosome_path}'"

        self._record_byte_size = get_dtype_length(self._dtype)
        assert self._record_byte_size > 0

        self._root = path

    @property
    def path(self):
        return self._root

    @property
    def chromosome_list(self) -> List[bytes]:
        return self._chromosomes.keys()

    @property
    def chromosome_sizes_nt(self) -> Dict[bytes, int]:
        return self._chromosomes.sizes_nt

    @property
    def chromosome_sizes_bytes(self) -> Dict[bytes, int]:
        return self._chromosomes.sizes_bytes

    @property
    def chromosome_recordnos(self) -> Dict[bytes, int]:
        return self._chromosomes.recordnos

    def log_details(self) -> None:
        logging.info(f"Database name: {self._args.output}")
        logging.info(f"Sequence max length: {self._dtype['sequence'][2:]}")
        logging.info("")
        logging.info("[bold]## Input files[/bold]")
        logging.info(f"hush files: {self._args.hush}")
        logging.info(f"oligo-melting files: {self._args.melting}")
        logging.info(f"OligoArrayAux files: {self._args.secondary}")
        logging.info("")
        logging.info("[bold]## Chromosome details[/bold]")
        logging.info(f"Expecting {len(self._chromosomes)} chromosomes.")
        logging.info("Chromosomes:")
        for chromosome, data in self._chromosomes.items():
            empty_label = "".join([" " for c in chromosome.decode()])
            logging.info(f"\t{chromosome.decode()} => size : {data.size_nt} (nt)")
            logging.info(f"\t{chromosome.decode()} => size : {data.size_bytes} (bytes)")
            logging.info(f"\t{empty_label} => recordno : {data.recordno}")
        logging.info("")
        logging.info("[bold]Record details[/bold]")
        logging.info(f"Record size in bytes: {self._record_byte_size}")
        logging.info(f"Dtype: {self._dtype}")

    def __read_next_record(self, IH: IO) -> bytes:
        return IH.read(self._record_byte_size)

    def buffer(
        self, chromosome: bytes, start_from_nt: int = 0, end_at_nt: int = -1
    ) -> Iterator[Record]:
        assert chromosome in self._chromosomes.keys()
        with open(os.path.join(self._root, f"{chromosome.decode()}.bin"), "rb") as IH:
            if start_from_nt > 0:
                position_in_bytes = self._chromosomes.get_chromosome(chromosome).index[
                    start_from_nt
                ]
                if position_in_bytes > 0:
                    IH.seek(position_in_bytes)
            record = self.__read_next_record(IH)
            while 0 != len(record):
                parsed_record = Record(record, self._dtype)
                if parsed_record["start"] > end_at_nt and end_at_nt > 0:
                    break
                yield parsed_record
                record = self.__read_next_record(IH)
