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
from typing import Any, Dict, IO, List, Iterator


def get_dtype_length(dtype) -> int:
    return sum([int(label.strip("><|SUuif")) for label in dtype.values()])


class Record(object):
    """DataBase Record"""

    _data: Dict[str, Any]

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

    def __init__(self, record: bytes, column_dtypes: Dict[str, str]):
        super(Record, self).__init__()
        self.__parse_bytes(record, column_dtypes)

    def __get__(self, key) -> Any:
        if key not in const.database_columns:
            raise KeyError(f"unrecognized key '{key}'")
        else:
            return self._data[key]

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
        return list(self._details["chromosomes"])

    @property
    def chromosome_sizes(self) -> Dict[bytes, int]:
        return dict(
            [
                (name, details["size"])
                for name, details in self._details["chromosomes"].items()
            ]
        )

    @property
    def chromosome_recordnos(self) -> Dict[bytes, int]:
        return dict(
            [
                (name, details["recordno"])
                for name, details in self._details["chromosomes"].items()
            ]
        )

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
