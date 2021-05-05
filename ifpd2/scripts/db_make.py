"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import argparse
import copy
from ifpd2 import asserts as ass
from ifpd2.asserts import enable_rich_assert
from ifpd2 import const, database as db, io
from ifpd2.scripts import arguments as ap
import logging
import numpy as np  # type: ignore
import os
import pandas as pd  # type: ignore
import pickle
from rich.progress import Progress, track  # type: ignore
from typing import Callable, Dict, List, Optional, Set, Tuple


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        "make",
        description="""
Assembles files into a database. At least one of the following is required,
to retrieve sequence: hush output (-O), or oligo-melting output (-T).

The name of each record should be in the format: 'name pos=chrom:start-end'.

Accepted formats:
hush            single-line sequence fasta file, where a ", XXX" is appended after the
                  sequence. XXX is the number of detected off-targets.
oligo-melting   .tsv file with six (6) columns: name, dG, dH, dS, Tm, Seq.
                  NOTE: the first line of the file is skipped during parsing because
                  it is expected to be a header line.
OligoArrayAux   .ct (Connectivity Table) file. Details on format available here:
                  https://rna.urmc.rochester.edu/Text/File_Formats.html

Files are merged by: hush header, oligo-melting name, and
OligoArrayAux name (3rd field of 1st column).
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Assemble a database.",
    )

    parser.add_argument("output", type=str, help="Path to database folder (output).")

    indata = parser.add_argument_group("input arguments")
    indata.add_argument(
        "-O",
        "--hush",
        metavar="hush_path",
        nargs="+",
        type=str,
        help="Path to hush output. Format details above.",
    )
    indata.add_argument(
        "-T",
        "--melting",
        metavar="oligomelting_path",
        nargs="+",
        type=str,
        help="Path to oligo-melting output. Format details above.",
    )
    indata.add_argument(
        "-S",
        "--secondary",
        metavar="oligoarrayaux_path",
        nargs="+",
        type=str,
        help="Path to OligoArrayAux output. Format details above.",
    )

    advanced = parser.add_argument_group("advanced arguments")
    advanced.add_argument(
        "-p",
        "--prefix",
        metavar="chromPrefix",
        default="",
        type=str,
        help="Prefix to be added to chromosome labels. Default: ''.",
    )
    advanced.add_argument(
        "-b",
        "--binsize",
        metavar="indexBin",
        default=const.DEFAULT_DATABASE_INDEX_BIN_SIZE,
        type=int,
        help="Binning for the index.",
    )

    parser = ap.add_version_option(parser)
    parser.set_defaults(parse=parse_arguments, run=run)

    return parser


def init_log(args: argparse.Namespace) -> None:
    logging.info(f"output folder: '{args.output}'")
    logging.info(f"chromosome prefix: '{args.prefix}'")
    if args.hush is not None:
        logging.info(f"hush: {args.hush}")
    if args.melting is not None:
        logging.info(f"oligo-melting: {args.melting}")
    if args.secondary is not None:
        logging.info(f"OligoArrayAux: {args.secondary}")


def check_files_exist(path_list: Optional[List[str]], software_name: str) -> None:
    if path_list is None:
        return
    for path in path_list:
        assert os.path.isfile(path), f"{software_name} output not found: {path}"


@enable_rich_assert
def parse_arguments(args: argparse.Namespace) -> argparse.Namespace:
    assert not os.path.isdir(
        args.output
    ), f"database folder already exists: {args.output}."

    check_files_exist(args.hush, "hush")
    check_files_exist(args.melting, "oligo-melting")
    check_files_exist(args.secondary, "OligoArrayAux")

    assert 1 <= ((args.hush is not None) + (args.melting is not None)), " ".join(
        [
            "one of the following is required:",
            "FASTA input, hush output, or oligo-melting output",
        ]
    )

    return args


def parse_input(
    path_list: Optional[List[str]], parse_function: Callable, software_name: str
) -> pd.DataFrame:
    if path_list is None:
        return pd.DataFrame()
    logging.info(f"parsing {software_name} output")
    return pd.concat([parse_function(path) for path in path_list])


def populate_db(
    db: pd.DataFrame,
    path_list: Optional[List[str]],
    parse_function: Callable,
    software_name: str,
) -> pd.DataFrame:
    if path_list is None:
        return db
    parsed_db = parse_input(path_list, parse_function, software_name)
    logging.info(f"adding {software_name} output to database")
    return db.merge(
        parsed_db,
        how="outer",
        left_index=True,
        right_index=True,
    )


def reduce_sequence_columns(df: pd.DataFrame) -> pd.DataFrame:
    logging.info("discarding redundant sequence columns")
    seq_columns = [c for c in df.columns if "sequence" in c]
    if 1 == len(seq_columns):
        return df
    df.drop(seq_columns[1:], axis=1, inplace=True)
    df.rename(columns=dict([(seq_columns[0], "sequence")]), inplace=True)
    return df


def parse_sequences(df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, str]]:
    logging.info("adding sequence feature columns: length and GC-content")
    sequence_length_list: List[int] = []
    gc_content_list: List[float] = []

    for sequence in track(
        df["sequence"].values, description="calculating GC-content", transient=True
    ):
        sequence = sequence.upper()
        sequence_length_list.append(len(sequence))
        gc_content_list.append(
            (sequence.count(b"G") + sequence.count(b"C")) / len(sequence)
        )

    df["gc_content"] = gc_content_list
    dtype = copy.copy(const.dtype_sequence_features)
    dtype["sequence"] = f"|S{max(sequence_length_list)}"
    return (df.astype(dtype), dtype)


def parse_record_headers(
    db: pd.DataFrame, chromosome_prefix: str = ""
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    logging.info("adding header feature columns: name, chromosome, start, end")
    name_list: List[str] = []
    name_length_set: Set[int] = set()
    chromosome_list: List[str] = []
    chromosome_length_set: Set[int] = set()
    start_list: List[int] = []
    end_list: List[int] = []

    for record in track(
        db.itertuples(),
        total=db.shape[0],
        description="parsing record headers",
        transient=True,
    ):
        name, position = record.Index.split(" ")
        name_list.append(name)
        name_length_set.add(len(name))
        chromosome, extremes = position.split("=")[1].split(":")
        chromosome_list.append(f"{chromosome_prefix}{chromosome}")
        chromosome_length_set.add(len(f"{chromosome_prefix}{chromosome}"))
        start, end = [int(x) for x in extremes.split("-")]
        assert (end - start) == len(
            record.sequence
        ), f"{end - start} != {len(record.sequence)}"
        start_list.append(start)
        end_list.append(end)

    db["name"] = name_list
    db["chromosome"] = chromosome_list
    db["start"] = start_list
    db["end"] = end_list
    ass.ert_in_dtype(db["start"].values.max(), "u4")
    ass.ert_in_dtype(db["end"].values.max(), "u4")
    db.reset_index(drop=True, inplace=True)

    dtype = copy.copy(const.dtype_header_features)
    dtype["name"] = f"|S{max(name_length_set)}"
    dtype["chromosome"] = f"|S{max(chromosome_length_set)}"
    return (db.astype(dtype), dtype)


def write_database(
    dbdf: pd.DataFrame, dtype: Dict[str, str], args: argparse.Namespace
) -> None:
    with Progress() as progress:
        chromosome_set: Set[bytes] = set(dbdf["chromosome"].values)
        chromosome_data = db.ChromosomeDict(args.binsize)
        chromosome_task = progress.add_task(
            "exporting chromosome",
            total=len(chromosome_set),
            transient=True,
        )
        for selected_chrom in chromosome_set:
            chromosome_db = dbdf.loc[selected_chrom == dbdf["chromosome"], :]

            logging.info(f"sorting records for {selected_chrom.decode()}")
            chromosome_db.sort_values(
                by="start", axis=0, kind="mergesort", inplace=True
            )
            logging.info(f"building index for {selected_chrom.decode()}")
            chromosome_data.add_chromosome(chromosome_db, dtype, progress)

            with open(
                os.path.join(args.output, f"{selected_chrom.decode()}.bin"), "wb"
            ) as IH:
                writing_track = progress.add_task(
                    f"writing {selected_chrom.decode()}.bin",
                    total=chromosome_db.shape[0],
                    transient=True,
                )
                for record in chromosome_db.to_records(
                    index=False, column_dtypes=dtype
                ):
                    IH.write(record.tobytes())
                    progress.update(writing_track, advance=1)
            progress.update(chromosome_task, advance=1)

    logging.info("writing db.pickle")
    with open(os.path.join(args.output, "db.pickle"), "wb") as OH:
        args.parse = None
        args.run = None
        pickle.dump(dict(chromosomes=chromosome_data, dtype=dtype, args=args), OH)


@enable_rich_assert
def run(args: argparse.Namespace) -> None:
    init_log(args)
    os.mkdir(args.output)

    dbdf = pd.DataFrame(columns=["name"])
    dbdf.set_index("name", inplace=True)

    dbdf = populate_db(dbdf, args.hush, io.parse_hush, "hush")
    dbdf = populate_db(dbdf, args.melting, io.parse_melting, "oligo-melting")
    dbdf = populate_db(dbdf, args.secondary, io.parse_secondary, "OligoArrayAux")

    dbdf = reduce_sequence_columns(dbdf)
    dbdf, dtype_sequence = parse_sequences(dbdf)
    dbdf, dtype_header = parse_record_headers(dbdf, args.prefix)

    dtype = dict()
    dtype.update(const.dtype_melting)
    dtype.update(const.dtype_hush)
    dtype.update(const.dtype_secondary)
    dtype.update(dtype_sequence)
    dtype.update(dtype_header)

    for column in dtype.keys():
        if column not in dbdf.columns:
            dbdf["column"] = np.repeat(np.nan, dbdf.shape[0])
    dbdf = dbdf.loc[:, const.database_columns]

    write_database(dbdf, dtype, args)

    logging.info("Done. :thumbs_up: :smiley:")
