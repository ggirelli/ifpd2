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
    """ChromosomeData"""

    __allowed_fields = ("recordno", "size_nt", "size_bytes")
    _data: Dict[bytes, Dict[str, Any]]
    _record_byte_size: int

    def __init__(
        self,
        chromosome_set: Set[bytes],
        dtype: Dict[str, str],
        index_bin_size: int = const.DEFAULT_DATABASE_INDEX_BIN_SIZE,
    ):
        super(ChromosomeData, self).__init__()

        self._record_byte_size = get_dtype_length(dtype)
        assert self._record_byte_size > 0

        self._data = {}
        for chromosome in set(chromosome_set):
            self._data[chromosome] = dict(index=ChromosomeIndex(index_bin_size))

    @property
    def record_byte_size(self) -> int:
        return self._record_byte_size

    def __len__(self) -> int:
        return len(self._data)

    def keys(self) -> List[bytes]:
        return list(self._data.keys())

    def items(self) -> List[Tuple[bytes, Dict[str, Any]]]:
        return list(self._data.items())

    @property
    def sizes_nt(self) -> Dict[bytes, int]:
        return dict(
            [(name, details["size_nt"]) for name, details in self._data.items()]
        )

    @property
    def sizes_bytes(self) -> Dict[bytes, int]:
        return dict(
            [(name, details["size_bytes"]) for name, details in self._data.items()]
        )

    @property
    def recordnos(self) -> Dict[bytes, int]:
        return dict(
            [(name, details["recordno"]) for name, details in self._data.items()]
        )

    def set(self, chromosome: bytes, key: str, value: int) -> None:
        if "index" == key:
            logging.warning(f"access precluded to '{key}' key")
        if chromosome not in self._data:
            raise KeyError(f"key '{chromosome.decode()}' not found")
        if key not in self.__allowed_fields:
            raise KeyError(f"key '{key}' not found")
        self._data[chromosome][key] = value

    def chromosome_details(self, chromosome: bytes) -> Dict[str, int]:
        return copy.copy(self._data[chromosome])

    def chromosome_index(self, chromosome: bytes) -> ChromosomeIndex:
        return copy.copy(self._data[chromosome]["index"])

    def populate(self, chrom_db: pd.DataFrame, track: Tuple[Progress, TaskID]) -> None:
        assert "chromosome" in chrom_db.columns
        selected_chrom = chrom_db["chromosome"][0]
        self.set(selected_chrom, "recordno", chrom_db.shape[0])
        self.set(selected_chrom, "size_nt", chrom_db["end"].values.max())
        self.set(
            selected_chrom,
            "size_bytes",
            chrom_db.shape[0] * self._record_byte_size,
        )
        self._data[selected_chrom]["index"].build(
            chrom_db, self._record_byte_size, track
        )

    def __check__(self) -> None:
        for chromosome, details in self._data.items():
            assert (
                "size_nt" in details
            ), f"missing nt size information for '{chromosome.decode()}'"
            assert (
                "size_bytes" in details
            ), f"missing byte size information for '{chromosome.decode()}'"
            assert (
                "recordno" in details
            ), f"missing size information for '{chromosome.decode()}'"


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
    _chromosomes: ChromosomeData
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

        assert isinstance(self._chromosomes, ChromosomeData)
        self._chromosomes.__check__()
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
        for chromosome, details in self._chromosomes.items():
            empty_label = "".join([" " for c in chromosome.decode()])
            logging.info(f"\t{chromosome.decode()} => size : {details['size_nt']} (nt)")
            logging.info(
                f"\t{chromosome.decode()} => size : {details['size_bytes']} (bytes)"
            )
            logging.info(f"\t{empty_label} => recordno : {details['recordno']}")
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
                position_in_bytes = self._chromosomes.chromosome_index(chromosome)[
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
