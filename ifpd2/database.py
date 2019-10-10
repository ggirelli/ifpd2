
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for oligo database management.
'''

import configparser as cp
import logging
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
from pathlib import Path
import shutil
import struct
from tqdm import tqdm

from ifpd2.logging import Loggable
from ifpd2.asserts import *
from ifpd2.oligo import Oligo, OligoBinary, OligoGroup

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
		if not os.path.isfile(chrom_path): return
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
		assert 'IFPD2DB' in config.keys(), assert_msg
		assert 'namelen' in config['IFPD2DB'].keys(), assert_msg
		assert 'chromlen' in config['IFPD2DB'].keys(), assert_msg
		assert 'oligok' in config['IFPD2DB'].keys(), assert_msg
		self._N = config['IFPD2DB'].getint('namelen')
		self._c = config['IFPD2DB'].getint('chromlen')
		self._k = config['IFPD2DB'].getint('oligok')
		self._dtype = f"{self.N}s {self.c}s i i f f f f {self.k}s i i f f"
		self._n_expected_fields = len(self.dtype.split(" "))
		self._n_bytes = struct.calcsize(self.dtype)

	def buffer(self, chrom):
		chrom_path = os.path.join(self.path, f"{chrom}.bin")
		if not os.path.isfile(chrom_path): return
		with open(chrom_path, "rb") as DBH:
			bytepack = DBH.read(self.n_bytes)
			while 0 != len(bytepack):
				yield struct.unpack(self.dtype, bytepack)
				bytepack = DBH.read(self.n_bytes)

class GenomicWindowSet(object):
	"""Genomic window manager."""

	_window_sets = None

	C = "chr18"			# Chromosome
	S = int(3000000)	# Region start coordinate (included)
	E = int(3500000)	# Region end coordinate (excluded)

	X = 20				# Number of probes to design
	Ws = None			# Window size (used when X is not provided)
	Wh = 0.1			# Window shift (as a percentage of the window size)

	Rs = int(8000)		# Region focus size, either in nt (Rs > 1)
						#  or fraction of Ws (0<Rs<=1)
						#  When provided in nt, it is applied only if Rs<Ws
	Rt = int(1000)		# Region focus step, either in nt (Rt > 1)
						#  or fraction of Rs (0<Rt<=1),
						#  for focus region expansion
	_growing = False

	def __init__(self):
		super(GenomicWindowSet, self).__init__()

	@property
	def window_sets(self):
		return self._window_sets
	
	@window_sets.setter
	def window_sets(self, window_sets):
		self._window_sets = pd.DataFrame(np.vstack(window_sets))
		self._window_sets.columns = [
			"start", "mid", "end", "cfr_start", "cfr_end", "w", "s"]
		self._window_sets.sort_values("end", inplace = True)
		self._window_sets.reset_index(inplace = True, drop = True)

	@property
	def wid(self):
		return self._w

	@property
	def current_window(self):
		return self.window_sets.iloc[self.wid, :]

	@property
	def reached_last_window(self):
		return self._reached_last_window

	@property
	def growing(self):
		return self._growing

	def _assert(self):
		assert_type(self.C, str, "C")
		assert_type(self.S, int, "S")
		assert_nonNeg(self.S, "S", True)
		assert_type(self.E, int, "E")
		assert_nonNeg(self.E, "E", True)
		assert self.S <= self.E

		assert_multiTypes(self.X, [int, type(None)], "X")
		assert_multiTypes(self.Ws, [int, type(None)], "Ws")
		if isinstance(self.X, int):
			assert_type(self.Ws, type(None), "Ws")
			assert_nonNeg(self.X, "X")
		else:
			assert_type(self.Ws, int, "Ws")

		assert_multiTypes(self.Ws, [int, type(None)], "Ws")
		if isinstance(self.Ws, int):
			assert_type(self.X, type(None), "X")
			assert_nonNeg(self.Ws, "Ws")
		else:
			assert_type(self.X, int, "X")

		assert_type(self.Wh, float, "Wh")
		assert_inInterv(self.Wh, 0, 1, "Wh")

		assert_multiTypes(self.Rs, [int, float], "Rs")
		if isinstance(self.Rs, int):
			assert self.Rs > 1
		else:
			assert_inInterv(self.Rs, 0, 1, "Rs")

		assert_multiTypes(self.Rt, [int, float], "Rt")
		if isinstance(self.Rt, int):
			assert self.Rt > 1
		else:
			assert_inInterv(self.Rt, 0, 1, "Rt")

	def _init_windows(self):
		# Build windows and central focus regions (CFR)
		if not os.path.isdir(os.path.join(self.out_path, "window_sets")):
			os.mkdir(os.path.join(self.out_path, "window_sets"))

		if self.S == self.E:
			assert_msg = "During full-chromosome search, provide a window size."
			assert_msg += " I.e., it is not possible to design X probes."
			assert not isinstance(self.Ws, type(None)), assert_msg

		if isinstance(self.X, int):
			if 1 == self.X:
				self.Ws = self.E - self.S
			else:
				self.Ws = np.floor((self.E-self.S)/(self.X+1)).astype("i")

		self._w = 0 # Window ID
		self._reached_last_window = False

		if self.S == self.E:
			self._growing = True
			self.window_sets = self.__mk_first_window()
		else:
			self.window_sets = self.__mk_all_window_sets()

	def __mk_all_window_sets(self):
		# Prepare all window sets in a region of interest
		window_starts = np.floor(np.arange(self.S,self.E,self.Ws)).astype("i")
		if 0 != (self.E-self.S)%self.Ws:
			window_starts = window_starts[:-1]
		if 1 != len(window_starts):
			window_starts = window_starts[:-1]

		window_mids = (window_starts+self.Ws/2
			).reshape((window_starts.shape[0], 1))
		window_borders = np.transpose(np.vstack(
			[window_starts, window_starts+self.Ws]))

		nWindows = window_borders.shape[0]

		if isinstance(self.Rs, float):
			self.Rs = int(self.Rs*self.Ws)
		if isinstance(self.Rt, float):
			self.Rt = int(self.Rs*self.Rt)
		if self.Rs < self.Ws:
			central_regions = np.hstack([np.floor(window_mids-self.Rs/2),
				np.floor(window_mids+self.Rs/2)])
			self.__focus_on_center = True
		else:
			nans = np.array([np.nan for x in range(window_mids.shape[0])])
			central_regions = np.vstack([nans, nans]).transpose()
			self.__focus_on_center = False

		window_sets = []
		for i in range(len(np.arange(0, 1, self.Wh))):
			winID = np.array(range(nWindows)).reshape((nWindows, 1))
			setID = np.repeat(i, nWindows).reshape((nWindows, 1))
			window_sets.append(np.hstack([
				window_borders+i*self.Wh*self.Ws,
				window_mids+i*self.Wh*self.Ws,
				central_regions+i*self.Wh*self.Ws,
				winID, setID])[ :, (0, 2, 1, 3, 4, 5, 6)])

		return window_sets

	def __mk_first_window(self):
		# Initialize only the first window
		mid = self.S+self.Ws/2
		if self.Rs < self.Ws:
			cstart, cend = (np.floor(mid-self.Rs/2), np.floor(mid+self.Rs/2))
			self.__focus_on_center = True
		else:
			cstart, cend = (np.nan, np.nan)
			self.__focus_on_center = False
		return [[self.S, mid, self.S+self.Ws, cstart, cend, 0, 0]]

	def _add_window(self):
		# Add a window, assigning it to the first set with no overlaps
		new_window = self.window_sets.iloc[-1, :].copy()
		new_window.iloc[:5] = new_window.iloc[:5] + self.Wh*self.Ws
		
		gData = np.array([(n, g['end'].max())
			for n,g in self.window_sets.groupby("s")])
		non_overlapping_group_ids = (np.where(gData[:,1] <= new_window[0])[0])

		new_window_id = 0
		new_set_id = self.window_sets['s'].max()+1
		if 0 != len(non_overlapping_group_ids):
			new_set_id = gData[non_overlapping_group_ids].min()
			new_window_id = self.window_sets[
				self.window_sets['s'] == new_set_id]['w'].max()+1

		new_window.iloc[5] = new_window_id
		new_window.iloc[6] = new_set_id
		self.window_sets = [self.window_sets.values,
			new_window.values.reshape((1,7))]

	def go_to_next_window(self):
		if self._growing:
			self._add_window()
			self._w += 1
		else:
			if self.wid < self.window_sets.shape[0]-1:
				self._w += 1
			if self.wid >= self.window_sets.shape[0]-1:
				self._reached_last_window = True

	def export_window_set(self):
		self.window_sets.loc[
			self.window_sets['s'] == self.current_window['s'],:].to_csv(
			os.path.join(self.window_set_path, "windows.tsv"), "\t")

class Walker(GenomicWindowSet, Loggable):
	"""Walker walks through the oligos stored in an ifpd2 database,
	assigns them to windows based on user-defined parameters, and then builds
	probe candidates."""

	out_path = "."
	reuse = False
	_threads = 1

	__current_oligos = []
	__walk_results = {}

	def __init__(self, db_path, logger = logging.getLogger()):
		GenomicWindowSet.__init__(self)
		Loggable.__init__(self, logger)
		self.__db = Database(db_path)

	@property
	def threads(self):
		return self._threads
	
	@threads.setter
	def threads(self, threads):
		assert_type(threads, int, threads)
		threads = max(1, threads)
		threads = min(mp.cpu_count(), threads)
		self._threads = threads

	@property
	def current_oligos(self):
		return self.__current_oligos

	def remove_oligos_starting_before(self, pos):
		self.__current_oligos = [o for o in self.current_oligos
			if o.start >= pos]

	@property
	def window_set_path(self):
		return os.path.join(self.out_path, "window_sets",
			f"set_{int(self.current_window['s'])}")

	@property
	def window_path(self):
		return os.path.join(self.window_set_path,
			f"window_{int(self.current_window['w'])}")

	@property
	def window_tag(self):
		window = self.current_window
		return f"{int(window['s'])}.{int(window['w'])}"

	@property
	def window_range(self):
		window = self.current_window
		return f"[{int(window['start'])}:{int(window['end'])}]"

	@property
	def walk_results(self):
		return self.__walk_results

	@property
	def config(self):
		config = cp.ConfigParser()

		config['AIM'] = {
			'Region' : f"{self.C}:{self.S}-{self.E}"
		}
		if not isinstance(self.X, type(None)):
			config['AIM']['Probe(s) number'] = str(self.X)

		config['WINDOWS'] = {}
		if not isinstance(self.Ws, type(None)):
			config['WINDOWS']['Window size'] = str(self.Ws)
		config['WINDOWS'].update({
			'Window step' : str(self.Wh),
			'Focus region size' :  str(self.Rs),
			'Focus region step' : str(self.Rt)
		})
		return config

	def _assert(self):
		GenomicWindowSet._assert(self)
		assert os.path.isdir(self.out_path)

	def start(self, fparse, fimport, fprocess, fpost, *args, **kwargs):
		self._assert()
		self._init_windows()
		self.print_prologue()
		self.__walk(fparse, fimport, fprocess, fpost, *args, **kwargs)

	def print_prologue(self):

		s  = f"* Walker *\n\n"
		s += f"Threads: {self.threads}\n"
		s += f"Database: '{self.__db.path}'\n"
		s += f"Region of interest: {self.C}:{self.S}-{self.E}\n"

		if not isinstance(self.X, type(None)):
			nWindows = int(self.window_sets['w'].max()+1)
			s += f"Aim to scout {nWindows} windows.\n\n"

		s += f"Using a central focus region of {self.Rs} nt,"
		s += f" in windows of size {self.Ws} nt,\n"
		s += f"built with a shift of {self.Ws*self.Wh} nt ({self.Wh*100}%).\n"
		if not isinstance(self.X, type(None)):
			nSets = int(self.window_sets["s"].max()+1)
			s += f"Thus, a total of {nSets} window sets will be explored.\n"

		self.log.info(s)

	def __walk(self, fparse, fimport, fprocess, fpost, *args, **kwargs):
		if 1 == self.threads:
			fexec = self.process_window
		else:
			pool = mp.Pool(np.min([self.threads, mp.cpu_count()]))
			self.log.info(f"Prepared a pool of {self.threads} threads.")
			fexec = lambda *args, **kwargs: pool.apply_async(
				self.process_window_parallel, args, kwargs)

		self._preprocess_window(fimport)
		self._load_windows_until_next_to_do(fimport)
		
		if self.reached_last_window and self._window_done():
			self.log.info("All windows pre-processed. Skipped database walk.")
			return

		self.r = 0 # Record ID
		self.rw = 0 # Walk step counter

		exec_results = []

		walk_destination = self.E
		if walk_destination == self.S:
			walk_destination = np.inf

		DBHpb = tqdm(self.__db.buffer(self.C),
			desc = "Parsing records", leave = None)
		for line in DBHpb:
			oligo_start, oligo_end = [int(x)
				for x in line.strip().split("\t")[2:4]]

			if oligo_start >= self.current_window['start']:
				if oligo_start >= self.current_window['end']:
					DBHpb.clear()
					
					exec_results.append(fexec(self.current_oligos,
						self.current_window, fprocess, fpost,
						*args, opath = self.window_path,
						loggerName = self.log.name, **kwargs))

					if self.reached_last_window:
						self.log.info("Reached last window")
						break

					self.go_to_next_window()
					self._preprocess_window(fimport)
					self._load_windows_until_next_to_do(fimport)
					if self.reached_last_window and self._window_done():
						self.log.info("Reached last window and done")
						break

					if 0 != len(self.current_oligos):
						self.remove_oligos_starting_before(
							self.current_window['start'])

				if oligo_end > walk_destination:	# End reached
					self.log.info("Reached destination")
					break

				oligo = Oligo(line, self.r)
				fparse(oligo, *args, **kwargs)

				if not np.isnan(oligo.score):
					self.current_oligos.append(oligo)

				self.rw += 1

			self.r += 1
		self.log.info(f"Parsed {self.rw}/{self.r} records.")

		if 1 < self.threads:
			for promise in exec_results:
				s, w, results = promise.get()
				if 0 == len(results): continue
				if s in self.walk_results.keys():
					self.walk_results[s][w] = results
				else:
					self.walk_results[s] = {w : results}

		DBHpb.close()
		pool.close()

	def _window_done(self):
		s = int(self.current_window['s'])
		w = int(self.current_window['w'])
		if s in self.walk_results.keys():
			if w in self.walk_results[s].keys():
				return isinstance(self.walk_results[s][w], list)
		return False

	def _load_windows_until_next_to_do(self, fimport):
		while self._window_done() and not self.reached_last_window:
			self.go_to_next_window()
			self._preprocess_window(fimport)

	def _preprocess_window(self, fimport):
		# Preprocess current window
		# Triggers import if previously run and matching current window

		if not os.path.isdir(self.window_set_path):
			os.mkdir(self.window_set_path)
		self.export_window_set()

		if not os.path.isdir(self.window_path):
			os.mkdir(self.window_path)
		else:
			if not self.reuse:
				shutil.rmtree(self.window_path)
				os.mkdir(self.window_path)

			win_file_path = os.path.join(self.window_path, "window.tsv")
			win_done = os.path.isfile(os.path.join(self.window_path, ".done"))

			if os.path.isfile(win_file_path) and win_done:
				win = pd.read_csv(win_file_path, sep='\t',
					header=None, index_col=0)

				if (win.transpose().values == self.current_window.values).all():
					self.log.info("Re-using previous results for window " +
						f"{self.window_tag} {self.window_range}")


					sid = int(self.current_window['s'])
					if not sid in self.walk_results.keys():
						self.walk_results[sid] = {}
					wid = int(self.current_window['w'])
					self.walk_results[sid][wid] = fimport(self.window_path)
					return
				else:
					shutil.rmtree(self.window_path)
					os.mkdir(self.window_path)
			else:
				shutil.rmtree(self.window_path)
				os.mkdir(self.window_path)

		self.current_window.to_csv(os.path.join(self.window_path, "window.tsv"),
			sep = "\t", index = True)
		with open(os.path.join(self.window_path, "walker.config"), "w+") as CPH:
			self.config.write(CPH)

	@staticmethod
	def process_window_parallel(oligos, window, fprocess, fpost, *args,
		N = 1 , opath = None, loggerName = None, **kwargs):
		# Wrapper of process_window function, for parallelization
		window_tag = f"{int(window['s'])}.{int(window['w'])}"

		logFormatter = logging.Formatter(Loggable.defaultfmt,
			datefmt = Loggable.datefmt)
		logger = logging.getLogger(f"ifpd2-window-{window_tag}")
		logger.setLevel(logging.DEBUG)
		logger = Loggable(logger)
		logPath = "{0}/{1}.log".format(opath, "window")
		logger.addFileHandler(logPath)
		logger.log.info(f"This log is saved at '{logPath}'.")


		mainLogger = logging.getLogger(loggerName)
		mainLogger.info(f"Window {window_tag} sent to pool.")
		results = Walker.process_window(oligos, window,
			fprocess, fpost, *args, N = N, opath = opath,
			loggerName = f"ifpd2-window-{window_tag}", **kwargs)
		mainLogger.info(f"Processor pool returned window {window_tag}.")
		return results

	@staticmethod
	def process_window(oligos, window, fprocess, fpost, *args,
		N = 1, opath = None, loggerName = None, **kwargs):
		# Process oligos from window using fprocess. Then, post-process them
		# with fopost. Requires at least N oligos to proceeed. If opath is
		# specified, a ".done" file is touched upon successful postprocessing.
		logger = logging.getLogger(loggerName)
		kwargs['loggerName'] = loggerName

		if len(oligos) >= N:
			oGroup = OligoGroup(oligos, logger)
			logger.info(f"Retrieved {oGroup.data.shape[0]} oligos for" +
				f" window {int(window['s'])}.{int(window['w'])} " +
				f"[{int(window['start'])}:{int(window['end'])}]")
			results = fprocess(oGroup, window, *args, **kwargs)
		else:
			logger.warning(f"Window {int(window['s'])}.{int(window['w'])}" +
				" does not have enough oligos " +
				f"{len(oligos)}/{N}, skipped.")
			results = []
		
		status, results = fpost(results, opath, *args, **kwargs)
		if status and not isinstance(opath, type(None)):
			Path(os.path.join(opath, ".done")).touch()

		return (int(window['s']), int(window['w']), results)

class WalkerBinary(Walker):
	"""Walker for binary databases."""

	def __init__(self, db_path, logger = logging.getLogger()):
		super(WalkerBinary, self).__init__(db_path, logger)
		self.__db = DatabaseBinary(db_path)
	
	def start(self, fparse, fimport, fprocess, fpost, *args, **kwargs):
		self._assert()
		self._init_windows()
		self.print_prologue()
		self.__walk(fparse, fimport, fprocess, fpost, *args, **kwargs)

	def __walk(self, fparse, fimport, fprocess, fpost, *args, **kwargs):
		if 1 == self.threads:
			fexec = self.process_window
		else:
			pool = mp.Pool(np.min([self.threads, mp.cpu_count()]))
			self.log.info(f"Prepared a pool of {self.threads} threads.")
			fexec = lambda *args, **kwargs: pool.apply_async(
				self.process_window_parallel, args, kwargs)

		self._preprocess_window(fimport)
		self._load_windows_until_next_to_do(fimport)
		
		if self.reached_last_window and self._window_done():
			self.log.info("All windows pre-processed. Skipped database walk.")
			return

		self.r = 0 # Record ID
		self.rw = 0 # Walk step counter

		exec_results = []

		walk_destination = self.E
		if walk_destination == self.S:
			walk_destination = np.inf

		self.log.info("Starting to walk...")
		for line in tqdm(self.__db.buffer(self.C),
			desc = f"Parsing binary records", leave = None):
			oligo = OligoBinary(line, self.r)

			if oligo.start >= self.current_window['start']:
				if oligo.start >= self.current_window['end']:
					
					exec_results.append(fexec(self.current_oligos,
						self.current_window, fprocess, fpost,
						*args, opath = self.window_path,
						loggerName = self.log.name, **kwargs))

					if self.reached_last_window:
						self.log.info("Reached last window")
						break

					self.go_to_next_window()
					self._preprocess_window(fimport)
					self._load_windows_until_next_to_do(fimport)
					if self.reached_last_window and self._window_done():
						self.log.info("Reached last window and done")
						break

					if 0 != len(self.current_oligos):
						self.remove_oligos_starting_before(
							self.current_window['start'])

				if oligo.end > walk_destination:	# End reached
					self.log.info("Reached destination")
					break

				fparse(oligo, *args, **kwargs)

				if not np.isnan(oligo.score):
					self.current_oligos.append(oligo)

				self.rw += 1

			self.r += 1
		self.log.info(f"Parsed {self.rw}/{self.r} records.")

		if 1 < self.threads:
			for promise in exec_results:
				s, w, results = promise.get()
				if 0 == len(results): continue
				if s in self.walk_results.keys():
					self.walk_results[s][w] = results
				else:
					self.walk_results[s] = {w : results}

		pool.close()
