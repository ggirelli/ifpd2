"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import copy
from ifpd2 import const
import logging
import numpy as np  # type: ignore
import os
import pickle
import pandas as pd  # type: ignore
from tqdm import tqdm  # type: ignore
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
        for bin_id in range(0, chrom_size_nt, self._bin_size):
            self._index[bin_id] = (0, 0)

    def build(self, chrom_db: pd.DataFrame, record_byte_size: int) -> None:
        for colname in ("chromosome", "start", "end"):
            assert colname in chrom_db.columns, f"missing '{colname}' column"

        self.__init_index(chrom_db)

        chromosome_set: Set[bytes] = set(chrom_db["chromosome"].values)
        assert 1 == len(chromosome_set)
        chromosome = list(chromosome_set)[0].decode()

        current_position = -1
        for i in tqdm(
            range(chrom_db.shape[0]),
            desc=f"building '{chromosome}' index",
            leave=False,
        ):
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

    def __getitem__(self, position_in_nt: int) -> int:
        assert self._index is not None
        binned_to = position_in_nt // self._bin_size
        return self._index[binned_to][0]


class ChromosomeData(object):
    """ChromosomeData"""

    __allowed_fields = ("recordno", "size_nt", "size_bytes")
    _data: Dict[bytes, Tuple[Dict[str, Any], ChromosomeIndex]]
    _record_byte_size: int

    def __init__(
        self,
        chromosome_set: Set[bytes],
        dtype: Dict[str, str],
        index_bin_size: int = 100000,
    ):
        super(ChromosomeData, self).__init__()

        self._record_byte_size = get_dtype_length(dtype)
        assert self._record_byte_size > 0

        self._data = {}
        for chromosome in set(chromosome_set):
            self._data[chromosome] = ({}, ChromosomeIndex(index_bin_size))

    @property
    def record_byte_size(self) -> int:
        return self._record_byte_size

    def __len__(self) -> int:
        return len(self._data)

    def keys(self) -> List[bytes]:
        return list(self._data.keys())

    @property
    def sizes_nt(self) -> Dict[bytes, int]:
        return dict(
            [(name, details[0]["size_nt"]) for name, details in self._data.items()]
        )

    @property
    def sizes_bytes(self) -> Dict[bytes, int]:
        return dict(
            [(name, details[0]["size_bytes"]) for name, details in self._data.items()]
        )

    @property
    def recordnos(self) -> Dict[bytes, int]:
        return dict(
            [(name, details[0]["recordno"]) for name, details in self._data.items()]
        )

    def set(self, chromosome: bytes, key: str, value: int) -> None:
        if chromosome not in self._data:
            raise KeyError(f"key '{chromosome.decode()}' not found")
        if key not in self.__allowed_fields:
            raise KeyError(f"key '{key}' not found")
        self._data[chromosome][0][key] = value

    def get_details(self, chromosome: bytes) -> Dict[str, int]:
        return copy.copy(self._data[chromosome][0])

    def get_index(self, chromosome: bytes) -> ChromosomeIndex:
        return copy.copy(self._data[chromosome][1])

    def populate(self, chrom_db: pd.DataFrame) -> None:
        assert "chromosome" in chrom_db.columns
        selected_chrom = chrom_db["chromosome"][0]
        self.set(selected_chrom, "recordno", chrom_db.shape[0])
        self.set(selected_chrom, "size_nt", chrom_db["end"].values.max())
        self.set(
            selected_chrom,
            "size_bytes",
            chrom_db.shape[0] * self._record_byte_size,
        )
        self._data[selected_chrom][1].build(chrom_db, self._record_byte_size)


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

    def to_csv(self, sep: str = ",") -> str:
        csv_fields: List[str] = []
        for x in const.database_columns:
            data = self._data[x]
            if isinstance(data, bytes):
                csv_fields.append(str(data.decode()))
            else:
                csv_fields.append(str(data))
        return sep.join(csv_fields)

    def __getitem__(self, key: str) -> Any:
        if key not in const.database_columns:
            raise KeyError(f"unrecognized key '{key}'")
        else:
            return self.data[key]

    def __repr__(self) -> str:
        return str(self._data)


class DataBase(object):
    """Buffering and checking class for ifpd2 database."""

    _root: str
    _details: Dict[str, Any]
    _record_byte_size: int

    def __init__(self, path: str):
        super(DataBase, self).__init__()
        assert os.path.isdir(path), f"cannot find database folder '{path}'"
        db_pickle_path = os.path.join(path, "db.pickle")
        assert os.path.isfile(db_pickle_path), f"'db.pickle is missing in '{path}'"

        with open(db_pickle_path, "rb") as IH:
            self._details = pickle.load(IH)

        for chromosome, details in self._details["chromosomes"].items():
            chromosome_path = os.path.join(path, f"{chromosome.decode()}.bin")
            assert os.path.isfile(
                chromosome_path
            ), f"missing expected chromosome file: '{chromosome_path}'"
            assert "size" in details, f"missing size information for '{chromosome}'"
            assert "recordno" in details, f"missing size information for '{chromosome}'"

        self._record_byte_size = get_dtype_length(self._details["dtype"])
        assert self._record_byte_size > 0

        self._root = path

    @property
    def chromosome_list(self) -> List[bytes]:
        return self._details["chromosomes"].keys()

    @property
    def chromosome_sizes_nt(self) -> Dict[bytes, int]:
        return self._details["chromosomes"].sizes_nt()

    @property
    def chromosome_sizes_bytes(self) -> Dict[bytes, int]:
        return self._details["chromosomes"].sizes_bytes()

    @property
    def chromosome_recordnos(self) -> Dict[bytes, int]:
        return self._details["chromosomes"].recordnos()

    def log_details(self) -> None:
        logging.info(f"Database name: {self._details['args'].output}")
        logging.info(f"Sequence max length: {self._details['dtype']['sequence'][2:]}")
        logging.info("")
        logging.info("[bold]## Input files[/bold]")
        logging.info(f"hush files: {self._details['args'].hush}")
        logging.info(f"oligo-melting files: {self._details['args'].melting}")
        logging.info(f"OligoArrayAux files: {self._details['args'].secondary}")
        logging.info("")
        logging.info("[bold]## Chromosome details[/bold]")
        logging.info(f"Expecting {len(self._details['chromosomes'])} chromosomes.")
        logging.info("Chromosomes:")
        for chromosome, details in self._details["chromosomes"].items():
            empty_label = "".join([" " for c in chromosome.decode()])
            logging.info(f"\t{chromosome.decode()} => size : {details['size']}")
            logging.info(f"\t{empty_label} => recordno : {details['recordno']}")
        logging.info("")
        logging.info("[bold]Record details[/bold]")
        logging.info(f"Record size in bytes: {self._record_byte_size}")
        logging.info(f"Dtype: {self._details['dtype']}")

    def __read_next_record(self, IH: IO) -> bytes:
        return IH.read(self._record_byte_size)

    def walk_chromosome(self, chromosome: bytes) -> Iterator[Record]:
        assert chromosome in self._details["chromosomes"].keys()
        with open(os.path.join(self._root, f"{chromosome.decode()}.bin"), "rb") as IH:
            record = self.__read_next_record(IH)
            while 0 != len(record):
                yield Record(record, self._details["dtype"])
                record = self.__read_next_record(IH)
