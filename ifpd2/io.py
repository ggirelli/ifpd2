"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import copy
from ifpd2 import asserts as ass
from ifpd2 import const
import logging
import os
import pandas as pd  # type: ignore
from typing import Dict, List, Set, Tuple
from tqdm import tqdm  # type: ignore


def get_dtype_length(dtype: Dict[str, str]) -> int:
    """Calculate dtype record length in bytes.

    Arguments:
        dtype {Dict[str, str]} -- column dtype Dict[column_name: column_dtype]

    Returns:
        int -- length of record in bytes
    """
    return sum([int(label.strip("><|SUuif")) for label in dtype.values()])


def parse_hush(path: str) -> pd.DataFrame:
    """Parse HUSH output.

    Arguments:
        path {str} -- path to HUSH output file

    Returns:
        pd.DataFrame -- parsed HUSH data
    """
    assert os.path.isfile(path), f"cannot find file '{path}'"
    logging.info(f"parsing: '{path}'")
    sequence_lengths: Set[int] = set()
    parsed_lines: List[Tuple[str, str, int]] = []
    with open(path) as IH:
        header = ""
        for line in tqdm(IH, desc="Parsing hush output", leave=False):
            if line.startswith(">"):
                header = line.strip()
            else:
                line_split = line.strip().split(",")
                parsed_lines.append(
                    (header[1:], line_split[0], int(line_split[1].strip()))
                )
                sequence_lengths.add(len(line_split[0]))
                header = ""
    hush_df = pd.DataFrame(parsed_lines, columns=["name", "sequence", "off_target_no"])
    hush_df.set_index("name", inplace=True)
    ass.ert_in_dtype(hush_df["off_target_no"].values.max(), "u4")
    dtype = copy.copy(const.dtype_hush)
    dtype["sequence"] = f"|S{max(sequence_lengths)}"
    return hush_df.astype(dtype)


def parse_melting(path: str, sep: str = "\t", header: bool = True) -> pd.DataFrame:
    """Parse oligo-melting output.

    Arguments:
        path {str} -- path to oligo-melting output

    Keyword Arguments:
        sep {str} -- column separator (default: {"t"})
        header {bool} -- whether to expect a header line (default: {True})

    Returns:
        pd.DataFrame -- parsed oligo-melting data
    """
    assert os.path.isfile(path), f"cannot find file '{path}'"
    logging.info(f"parsing: '{path}'")
    expected_columns = copy.copy(const.dtype_melting)
    melting_df = pd.read_csv(path, sep=sep, header=None, skiprows=1 if header else 0)
    assert melting_df.shape[1] == len(expected_columns)
    melting_df.columns = list(expected_columns.keys())
    melting_df.set_index("name", inplace=True)
    expected_columns.pop("name", None)
    expected_columns[
        "sequence"
    ] = f'|S{max(set([len(seq) for seq in melting_df["sequence"].values]))}'
    return melting_df.astype(expected_columns)


def parse_secondary(path: str) -> pd.DataFrame:
    """Parse OligoArrayAux .ct output.

    Arguments:
        path {str} -- path to OligoArrayAux .ct output

    Returns:
        pd.DataFrame -- parsed OligoArrayAux .ct data
    """
    assert os.path.isfile(path), f"cannot find file '{path}'"
    logging.info(f"parsing: '{path}'")
    parsed_lines: List[Tuple[str, float]] = []
    with open(path, "r") as IH:
        for line in tqdm(IH, desc="Parsing OligoArrayAux output", leave=False):
            if "dG = " not in line:
                continue
            _, dG_string, name = line.strip().split("\t")
            parsed_lines.append((name, float(dG_string[5:])))
    secondary_df = pd.DataFrame(parsed_lines, columns=["name", "ss_dG"])
    secondary_df.set_index("name", inplace=True)
    ass.ert_in_dtype(secondary_df["ss_dG"].values.max(), "f4")
    return secondary_df.astype(const.dtype_secondary)
