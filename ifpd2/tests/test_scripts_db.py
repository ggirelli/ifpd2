"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import subprocess
import pickle


def remove_test_db():
    subprocess.run(["rm", "-r", "test_db"], check=True)


def test_db_make():
    remove_test_db()
    script_output = subprocess.check_output(
        [
            "ifpd2",
            "db",
            "make",
            "test_db",
            "-O",
            "test_data/test.hush.fa",
            "-T",
            "test_data/test.tm.tsv",
            "-S",
            "test_data/test.ss.ct",
            "-p",
            "chr",
        ]
    )
    if b"Error" in script_output:
        raise AssertionError
    with open("test_data/test_db/chr16.bin", "rb").read() as test_db, open(
        "test_db/chr16.bin", "rb"
    ).read() as new_db:
        if test_db != new_db:
            raise AssertionError
        test_pickle = pickle.load(open("test_data/test_db/db.pickle", "rb"))
        new_pickle = pickle.load(open("test_db/db.pickle", "rb"))
        if test_pickle != new_pickle:
            raise AssertionError


def test_db_check():
    script_output = subprocess.check_output(
        [
            "ifpd2",
            "db",
            "check",
            "test_db",
        ]
    )
    if b"Error" in script_output:
        raise AssertionError(script_output)


def test_db_dump():
    script_output = subprocess.check_output(
        [
            "ifpd2",
            "db",
            "dump",
            "test_db",
        ]
    )
    with open("test_data/test.dump.txt", "rb").read() as test_dump:
        if script_output != test_dump:
            raise AssertionError


def test_db_info():
    script_output = subprocess.check_output(
        [
            "ifpd2",
            "db",
            "info",
            "test_db",
        ]
    )
    if b"Error" in script_output:
        raise AssertionError


def test_cleanup():
    remove_test_db()
