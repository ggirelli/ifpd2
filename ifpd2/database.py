"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
import copy
from ifpd2 import const
from ifpd2.chromosome import ChromosomeData, ChromosomeDict
from ifpd2.io import get_dtype_length
import logging
import numpy as np  # type: ignore
import os
import pickle
import pandas as pd  # type: ignore
from typing import Any, Dict, List


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
        """Parse record bytes.

        Parses a record from its bytes based on column dtypes.

        Arguments:
            record {bytes} -- record bytes
            column_dtypes {Dict[str, str]} -- column dtypes
        """
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
        """
        Arguments:
            path {str} -- absolute path to database folder
        """
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

    @property
    def record_byte_size(self) -> int:
        return self._record_byte_size

    @property
    def dtype(self) -> Dict[str, str]:
        return copy.copy(self._dtype)

    def get_chromosome(self, chromosome: bytes) -> ChromosomeData:
        return self._chromosomes.get_chromosome(chromosome)

    def log_details(self) -> None:
        """Log database details."""
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
