"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

__version__ = "1.0.0-alpha"

dtype_melting = dict(
    [
        ("name", "|S"),
        ("Tm_dG", "<f4"),
        ("Tm_dH", "<f4"),
        ("Tm_dS", "<f4"),
        ("Tm", "<f4"),
        ("sequence", "|S"),
    ]
)
dtype_hush = dict(
    [
        ("sequence", "|S"),
        ("off_target_no", "<u4"),
    ]
)
dtype_secondary = dict([("ss_dG", "<f4")])
dtype_sequence_features = dict(
    [
        ("sequence", "|S"),
        ("sequence_length", "<u4"),
        ("gc_content", "<f4"),
    ]
)
dtype_header_features = dict(
    [
        ("name", "|S"),
        ("chromosome", "|S"),
        ("start", "<u4"),
        ("end", "<u4"),
    ]
)

DEFAULT_DATABASE_INDEX_BIN_SIZE = 100000

database_columns = [
    "name",
    "chromosome",
    "start",
    "end",
    "sequence_length",
    "sequence",
    "gc_content",
    "off_target_no",
    "Tm_dG",
    "Tm_dH",
    "Tm_dS",
    "Tm",
    "ss_dG",
]
