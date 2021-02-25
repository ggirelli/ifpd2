"""
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
"""

import configparser as cp
import os
import struct


class Database(object):
    """Class for ifpd2 database buffering."""

    def __init__(self, path):
        super(Database, self).__init__()
        assert os.path.isdir(path)
        self._path = path

    @property
    def path(self):
        return self._path

    def buffer(self, chrom):
        chrom_path = os.path.join(self.path, f"{chrom}.tsv")
        if not os.path.isfile(chrom_path):
            return
        with open(chrom_path, "r") as DBH:
            header = next(DBH)
            if not all([isinstance(s, str) for s in header.split("\t")]):
                DBH.seek(0)
            for line in DBH:
                yield line


class DatabaseBinary(Database):
    """Class for ifpd2 binary database buffering."""

    def __init__(self, path):
        super(DatabaseBinary, self).__init__(path)
        self.__check_integrity()

    @property
    def dtype(self):
        return self._dtype

    @property
    def n_expected_fields(self):
        return self._n_expected_fields

    @property
    def n_bytes(self):
        return self._n_bytes

    @property
    def N(self):
        return self._N

    @property
    def c(self):
        return self._c

    @property
    def k(self):
        return self._k

    def __check_integrity(self):
        assert_msg = f"Database at '{self.path}' does not pass integrity check."
        config_path = os.path.join(self.path, ".config")
        assert os.path.isfile(config_path), assert_msg
        config = cp.ConfigParser()
        config.read(config_path)
        assert "IFPD2DB" in config.keys(), assert_msg
        assert "namelen" in config["IFPD2DB"].keys(), assert_msg
        assert "chromlen" in config["IFPD2DB"].keys(), assert_msg
        assert "oligok" in config["IFPD2DB"].keys(), assert_msg
        self._N = config["IFPD2DB"].getint("namelen")
        self._c = config["IFPD2DB"].getint("chromlen")
        self._k = config["IFPD2DB"].getint("oligok")
        self._dtype = f"{self.N}s {self.c}s i i f f f f {self.k}s i i f f"
        self._n_expected_fields = len(self.dtype.split(" "))
        self._n_bytes = struct.calcsize(self.dtype)

    def buffer(self, chrom):
        chrom_path = os.path.join(self.path, f"{chrom}.bin")
        if not os.path.isfile(chrom_path):
            return
        with open(chrom_path, "rb") as DBH:
            bytepack = DBH.read(self.n_bytes)
            while 0 != len(bytepack):
                yield struct.unpack(self.dtype, bytepack)
                bytepack = DBH.read(self.n_bytes)
