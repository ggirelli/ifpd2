"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import subprocess
import pickle


def remove_test_db():
    subprocess.run(["rm", "-r", "test_db"])


def test_db_make():
    remove_test_db()
    script_output = subprocess.check_output(
        [
            "ifpd2",
            "db",
            "make",
            "-O",
            "test_data/test.hush.fa",
            "-T",
            "test_data/test.tm.tsv",
            "-S",
            "test_data/test.ss.ct",
            "-p",
            "chr",
            "test_db",
        ]
    )
    assert b"Error" not in script_output
    test_db = open("test_data/test_db/chr16.bin", "rb").read()
    new_db = open("test_db/chr16.bin", "rb").read()
    assert test_db == new_db
    test_pickle = pickle.load(open("test_data/test_db/db.pickle", "rb"))
    new_pickle = pickle.load(open("test_db/db.pickle", "rb"))
    assert test_pickle == new_pickle


def test_db_check():
    script_output = subprocess.check_output(
        [
            "ifpd2",
            "db",
            "check",
            "test_db",
        ]
    )
    assert b"Error" not in script_output, script_output


def test_db_dump():
    script_output = subprocess.check_output(
        [
            "ifpd2",
            "db",
            "dump",
            "test_db",
        ]
    )
    test_dump = open("test_data/test.dump.txt", "rb").read()
    assert script_output == test_dump


def test_db_info():
    script_output = subprocess.check_output(
        [
            "ifpd2",
            "db",
            "info",
            "test_db",
        ]
    )
    assert b"Error" not in script_output


def test_cleanup():
    remove_test_db()
