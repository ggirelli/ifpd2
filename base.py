
# ==============================================================================
# 
# 190912 - Gabriele Girelli
# Project: COOLFISH
# 
# Aim: implement probe design
#	
# ==============================================================================

# PARAMETERS ===================================================================

out_path = "/mnt/data/COOLFISH/ifpd2_out/test"
db_path = "/mnt/data/COOLFISH/mm10_chr18_selected_regions.tsv"
reuse = True

# DEPENDENCIES =================================================================

import logging
import numpy as np
import os
import pandas as pd
from pathlib import Path
import shutil
from tqdm import tqdm

# CLASSES ======================================================================

class Loggable(object):
	"""docstring for Loggable"""
	def __init__(self, logger = logging.getLogger(),
		formatter = logging.Formatter(
			'%(asctime)s %(levelname)s:\t%(message)s',
			datefmt = '%d/%m/%Y %I:%M:%S')):
		super(Loggable, self).__init__()
		self._logger = logger
		self._formatter = formatter

	@property
	def log(self):
		return self._logger

	def addFileHandler(self, path, mode="a+", level=logging.DEBUG):
		fileHandler = logging.FileHandler(path, "a+")
		fileHandler.setFormatter(self._formatter)
		fileHandler.setLevel(level)
		self.log.addHandler(fileHandler)

class OligoProbeBuilder(Loggable):
	"""OligoProbeBuilder contains methods to select non-overlapping sets of
	oligos (i.e., probes) from an oligo DataFrame in input."""

	N = int(96)			# Number of oligos per probe
	D = int(2)			# Minimum distance between consecutive oligos
	Tr = 10.0			# Melting temperature range half-width
	Ps = int(10000)		# Probe size threshold, in nt (Ps > 1)

	def __init__(self, logger = logging.getLogger):
		super(OligoProbeBuilder, self).__init__(logger)

	def _assert(self):

		assert_type(self.N, int, "N")
		assert_nonNeg(self.N, "N")

		assert_type(self.D, int, "D")
		assert_nonNeg(self.D, "D")

		assert_type(self.Tr, float, "Tr")
		assert_nonNeg(self.Tr, "Tr")

		assert_type(self.Ps, int, "Ps")
		assert self.Ps > 1

	def get_non_overlapping_probes(self, oData):
		# Gets all sets of N consecutive non-overlapping oligos with minimum
		# distance equal to D and melting temperature range of +-Tr
		assert_type(oData, pd.DataFrame, "oData")

		edges = []
		start_positions = oData["start"].values
		end_positions = oData["end"].values

		for i in range(oData.shape[0]):
			edges.append(np.logical_or(
				end_positions+self.D < oData["start"].values[i],
				start_positions-self.D >= oData["end"].values[i]
			))

		A = np.vstack(edges).astype('i')

		path_set = set()
		path_start_set = set()
		for i in range(A.shape[0]-1):
			if i in path_start_set:
				continue
			if 0 != len(A[i, i+1].nonzero()[0]):
				continue

			path = [i]
			k = i
			if 1 in A[k,(k+1):]:
				j = A[k,(k+1):].argmax()+k+1
				path.append(j)
			else:
				continue

			while (j+1) < A.shape[0]:
				k = j
				if 1 in A[k,(k+1):]:
					j = A[k,(k+1):].argmax()+k+1
					path.append(j)
				else:
					break

			if self.N > len(path):
				continue
			elif self.N == len(path):
				pData = oData.iloc[path, :]
				probe_size = pData['end'].values[0]-pData['start'].values[0]
				if self.Ps < probe_size:
					continue
				Tm_range = pData['Tm'].max()-pData['Tm'].min()
				if 2*self.Tr < Tm_range:
					continue
				path_set.add(tuple(path))
				path_start_set.add(path[0])
			else:
				for j in range(len(path) - self.N + 1):
					subpath = path[j:(j + self.N)]
					pData = oData.iloc[subpath, :]
					probe_size = pData['end'].values[0]-pData['start'].values[0]
					if self.Ps < probe_size:
						continue
					Tm_range = pData['Tm'].max()-pData['Tm'].min()
					if 2*self.Tr < Tm_range:
						continue
					path_set.add(tuple(subpath))
					path_start_set.add(subpath[0])

		probe_list = []
		if 0 != len(path_set):
			for path in list(path_set):
				probe_list.append(OligoProbe(oData.iloc[list(path), :]))

		return(probe_list)

class OligoWalker(OligoProbeBuilder):
	"""OligoWalker walks through the oligos stored in an ifpd2 database,
	assigns them to windows based on user-defined parameters, and then builds
	probe candidates."""

	db_path = None
	out_path = "."
	reuse = False

	C = "chr18"			# Chromosome
	S = int(3000000)	# Region start coordinate (included)
	E = int(4000000)	# Region end coordinate (excluded)
	X = 20				# Number of probes to design
	Ws = None			# Window size (used when X is not provided)
	Wh = 0.1			# Window shift (as a percentage of the window size)

	Rs = int(1000)		# Region focus size, either in nt (Rs > 1)
						#  or fraction of Ws (0<Rs<=1)
						#  When provided in nt, it is applied only if Rs<Ws
	Rt = int(1000)		# Region focus step, either in nt (Rt > 1)
						#  or fraction of Rs (0<Rt<=1),
						#  for focus region expansion

	F = (100, 1000)		# Threshold on number of off-targets (range)
	Gs = (0.0, 0.2)		# dG of SS either as kcal/mol (negative)
						#  or as fraction dG of hybrid (0<=Gs<=1)
						#  (range)
	Ot = .1				# Step for oligo score relaxation

	def __init__(self, db_path, logger = logging.getLogger):
		super(OligoWalker, self).__init__(logger)
		self.db_path = db_path

	def start(self):
		self._assert()
		self.__mk_windows()
		self.__print_prologue()
		self.__walk_db()
		pass

	def _assert(self):
		super(OligoWalker, self)._assert()

		assert os.path.isfile(self.db_path)
		assert os.path.isdir(self.out_path)

		assert_type(self.C, str, "C")
		assert_type(self.S, int, "S")
		assert_nonNeg(self.S, "S")
		assert_type(self.E, int, "E")
		assert_nonNeg(self.E, "E")
		assert self.S < self.E

		assert_multiTypes(self.X, [int, type(None)], "X")
		assert_type(self.Ws, type(None), "Ws")
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
	
		assert_type(self.F, tuple, "F")
		assert 2 == len(self.F)
		for i in range(2):
			assert_type(self.F[i], int, f"F[{i}]")
			assert_nonNeg(self.F[i], f"F[{i}]")
		assert self.F[1] >= self.F[0]

		assert_type(self.Gs, tuple, "Gs")
		assert 2 == len(self.Gs)
		for i in range(2):
			assert_type(self.Gs[i], float, f"Gs[{i}]")
			assert self.Gs[i] <= 1
		assert all(np.array(self.Gs) < 0) or all(np.array(self.Gs) >= 0)
		if self.Gs[0] >= 0:
			assert self.Gs[1] >= self.Gs[0]
		else:
			assert self.Gs[1] <= self.Gs[0]

		assert_type(self.Ot, float, "Ot")
		assert 0 < self.Ot and 1 >= self.Ot

	def __mk_windows(self):
		# Build windows and central focus regions (CFR)
		
		if isinstance(self.X, int):
			self.Ws = np.floor((self.E-self.S)/(self.X+1)).astype("i")

		window_starts = np.floor(np.arange(self.S,self.E,self.Ws)).astype("i")
		if 0 != (self.E-self.S)%self.Ws:
			window_starts = window_starts[:-1]

		window_mids = (window_starts[:-1]+self.Ws/2
			).reshape((window_starts.shape[0]-1, 1))
		window_borders = np.transpose(np.vstack(
			[window_starts[:-1], window_starts[:-1]+self.Ws]))

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

		self.window_sets = []
		for i in range(len(np.arange(0, 1, self.Wh))):
			winID = np.array(range(nWindows)).reshape((nWindows, 1))
			setID = np.repeat(i, nWindows).reshape((nWindows, 1))
			self.window_sets.append(np.hstack([
				window_borders+i*self.Wh*self.Ws,
				window_mids+i*self.Wh*self.Ws,
				central_regions+i*self.Wh*self.Ws,
				winID, setID])[ :, (0, 2, 1, 3, 4, 5, 6)])

		self.window_sets = pd.DataFrame(np.vstack(self.window_sets))
		self.window_sets.columns = ["start", "mid", "end",
			"cfr_start", "cfr_end", "w", "s"]
		self.window_sets = self.window_sets.sort_values("end")
		self.window_sets.index = range(self.window_sets.shape[0])

	def __print_prologue(self):
		nProbes = self.window_sets['w'].max().astype('i')+1
		nSets = self.window_sets["s"].max().astype('i')+1

		s  = f"\nDatabase: '{self.db_path}'\n"
		s += f"Region of interest: {self.C}:{self.S}-{self.E}\n"
		s += f"Aim to build {nProbes} probes, each with {self.N} oligos.\n\n"

		s += f"Using a central focus region of {self.Rs} nt,"
		s += f" in windows of size {self.Ws} nt,\n"
		s += f"built with a shift of {self.Ws*self.Wh} nt ({self.Wh*100}%).\n"
		s += f"Thus, a total of {nSets} window sets will be explored.\n\n"

		s += f"Off-target threshold range set at {self.F}.\n"
		if isinstance(self.Gs[0], int):
			s += f"Threshold on the delta free energy of the most stable"
			s += f" secondary structure set at range {self.Gs} kcal/mol.\n"
		else:
			s += f"Threshold on the delta free energy of the most stable"
			s += f" secondary structure\nset at range"
			s += f" {[t*100 for t in self.Gs]}% of the delta free energy of"
			s += f" hybridization.\n"

		s += f"\nMelting temperature range of {2*self.Tr} degC.\n"
		s += f"Minimum distance between consecutive oligos in a probe"
		s += f" set at {self.D} nt.\n"
		s += f"Probe size threshold set at {self.Ps} nt.\n"

		self.log.info(s)

	def __norm_value_in_range(self, v, r):
		if 0 == r[1]: return np.nan
		return (v - r[0])/(r[1] - r[0])

	def __norm_off_targets(self, oligo):
		nOT = oligo['nOT'].values
		if nOT <= self.F[0]: return 0
		if nOT > self.F[1]: return np.inf
		return self.__norm_value_in_range(nOT, self.F)

	def __norm_secondary_structure_dG(self, oligo):
		ss_dG = oligo['ss_dG'].values
		if isinstance(self.Gs[0], int):
			if ss_dG >= self.Gs[0]: return 0
			if ss_dG < self.Gs[1]: return np.inf
			return self.__norm_value_in_range(ss_dG, self.Gs)
		else:
			tm_dG = oligo['tm_dG'].values
			if ss_dG >= tm_dG*self.Gs[0]: return 0
			if ss_dG < tm_dG*self.Gs[1]: return np.inf
			return self.__norm_value_in_range(ss_dG, [tm_dG*f for f in self.Gs])

	def __calc_oligo_score(self, oligo):
		norm_nOT = self.__norm_off_targets(oligo)
		norm_ss_dG = self.__norm_secondary_structure_dG(oligo)
		return(norm_nOT * norm_ss_dG)

	def __walk_db(self):
		dtype = ['name', 'chrom', 'start', 'end', 'tm_dG',
			'dH', 'dS', 'Tm', 'seq', 'nOT', 'k', 'ss_dG']

		DBH = open(self.db_path, "r")
		next(DBH)

		self.current_oligos = []
		self.probe_candidates = []

		self.w = 0 # Window ID
		self.r = 0 # Record ID
		pbDBH = tqdm(DBH, leave = None, desc = "Parsing records")
		for oligo in pbDBH:
			oligo = oligo.strip().split("\t")
			for i in [2, 3, 9, 10]:
				oligo[i] = int(oligo[i])
			for i in [4, 5, 6, 7, 11]:
				oligo[i] = float(oligo[i])
			oligo = pd.DataFrame(oligo).transpose()
			oligo.columns = dtype

			if oligo['start'].values[0] >= self.S: # Skip all oligos before this
				end_of_current_window = self.window_sets['end'].values[self.w]
				if oligo['start'].values[0] >= end_of_current_window:
					pbDBH.clear()
					self.__select_from_window()
					if self.w == self.window_sets.shape[0]-1:
						break
					self.w += 1
				if oligo['end'].values[0] > self.E:	# End reached
					break
				oligo['score'] = self.__calc_oligo_score(oligo)
				if not np.isnan(oligo['score'].values[0]):
					self.current_oligos.append(oligo)
				self.r += 1
		self.log.info(f"Walked through {self.r} records.")

		DBH.close()

	def __import_window_output(self, win_path):
		oPath = os.path.join(win_path, "oligos.tsv")
		assert os.path.isfile(oPath)
		pPath = os.path.join(win_path, "probe_paths.tsv")
		assert os.path.isfile(pPath)
		
		oligos = pd.read_csv(oPath, "\t", index_col = 0)
		paths = pd.read_csv(pPath, "\t", index_col = 0)

		self.log.info(f"Loading {paths.shape[0]} probe candidates.")

		probe_list = []
		for p in paths.iloc[:, 0]:
			probe_list.append(OligoProbe(oligos.loc[
				[int(o) for o in p.split(",")]]))

		return(probe_list)

	def __select_from_window(self):
		window = self.window_sets.iloc[self.w,:]
		window_tag = f"{int(window['s'])}.{int(window['w'])}"
		window_range = f"[{int(window['start'])}:{int(window['end'])}]"

		set_path = os.path.join(self.out_path,f"set_{int(window['s'])}")
		if not os.path.isdir(set_path):
			os.mkdir(set_path)
			self.window_sets.loc[
				self.window_sets['s'] == window['s'],:].to_csv(
				os.path.join(set_path, "windows.tsv"), "\t")
		win_path = os.path.join(set_path,f"window_{int(window['w'])}")

		if not os.path.isdir(win_path):
			os.mkdir(win_path)
		else:
			if not self.reuse:
				shutil.rmtree(win_path)
				os.mkdir(win_path)
			win_file_path = os.path.join(win_path, "window.tsv")
			win_done = os.path.isfile(os.path.join(win_path, ".done"))
			if os.path.isfile(win_file_path) and win_done:
				win = pd.read_csv(win_file_path, sep='\t',
					header=None, index_col=0)
				if (win.transpose().values == window.values).all():
					self.log.info("Re-using previous results for window " +
						f"{window_tag} {window_range}")
					return(self.__import_window_output(win_path))
		
		window.to_csv(
			os.path.join(win_path, "window.tsv"), sep = "\t", index = True)

		self.addFileHandler(os.path.join(win_path, "window.log"))

		if len(self.current_oligos) >= self.N:
			oGroup = OligoGroup(self.current_oligos, self.log)
			self.log.info(f"Retrieved {oGroup.data.shape[0]} oligos for" +
				f" window {window_tag} {window_range}")
			probe_list = self.__build_probe_candidates(oGroup)
		else:
			self.log.warning(f"Window {window_tag} does not have enough" +
				f" oligos {len(self.current_oligos)}/{self.N}, skipped.")
		
		self.log.info(f"Built {len(probe_list)} oligo probe candidates")
		self.log.handlers = self.log.handlers[:-1]

		probe_df = pd.concat([op.featDF for op in probe_list],
			ignore_index = True)
		probe_df.to_csv(os.path.join(win_path, "probe_feat.tsv"), "\t")

		pd.concat([p.data for p in probe_list]).drop_duplicates().to_csv(
			os.path.join(win_path, "oligos.tsv"), "\t")

		probe_paths = []
		for pi in range(len(probe_list)):
			probe_paths.append([",".join([str(x)
				for x in probe_list[pi].data.index.tolist()])])
		probe_paths = pd.DataFrame(probe_paths)
		probe_paths.columns = ["cs_oligos"]
		probe_paths.to_csv(os.path.join(win_path, "probe_paths.tsv"), "\t")
		
		Path(os.path.join(win_path, ".done")).touch()

		if self.w+1 < self.window_sets.shape[0]:
			next_window_start = self.window_sets.iloc[self.w+1, :]['start']
			self.current_oligos = [o for o in self.current_oligos
				if o['start'].values[0] >= next_window_start]

		return(probe_list)

	def __build_probe_candidates(self, oGroup):
		# Applies oligo filters to the oligo group,
		# expands the focus group if needed, and build probe candidates
		
		if np.isnan(self.window_sets.loc[self.w, 'cfr_start']):
			oGroup.focus_all()
		else:
			oGroup.set_focus_window(self.window_sets.loc[self.w, 'cfr_start'],
				self.window_sets.loc[self.w, 'cfr_end'])
			oGroup.expand_focus_to_n_oligos(self.N)

		probe_list = self.__explore_filter(oGroup)
		if not np.isnan(self.window_sets.loc[self.w, 'cfr_start']):
			while 0 == len(probe_list):
				if not oGroup.expand_focus_by_step(self.Rt):
					break
				probe_list = self.__explore_filter(oGroup)
		else:
			self.log.warning("No CFR expansion, whole window already included.")
		
		return(probe_list)

	def __explore_filter(self, oGroup):
		# Explores the 0-to-1 score threshold range and stops as soon as one
		# probe candidate passes all user-defined thresholds
		
		nOligos_in_focus_window = oGroup.get_n_focused_oligos(True)

		score_thr = 0
		oGroup.apply_threshold(score_thr)
		nOligos_prev_score_thr = oGroup.get_n_focused_oligos(True)

		self.log.info(f"Set oligo score threshold at {score_thr:.3f}" +
			f" ({nOligos_prev_score_thr} oligos usable)...")

		probe_list = self.get_non_overlapping_probes(
			oGroup.get_focused_oligos(True))
		while 0 == len(probe_list):
			score_thr += self.Ot
			if score_thr > 1:
				break

			oGroup.apply_threshold(score_thr)
			nOligos = oGroup.get_n_focused_oligos()
			nOligosUsable = oGroup.get_n_focused_oligos(True)
			if nOligosUsable == nOligos_prev_score_thr:
				continue
			if 0 == nOligosUsable:
				continue

			self.log.info(f"Relaxed oligo score threshold to {score_thr:.3f}" +
				f" ({nOligosUsable} oligos usable)...")
			probe_list = self.get_non_overlapping_probes(
				oGroup.get_focused_oligos(True))

			if nOligosUsable == nOligos_in_focus_window:
				self.log.warning(
					"All oligos included. Score relaxation ineffective.")
				break
			nOligos_prev_score_thr = nOligos

		return(probe_list)

class OligoGroup(Loggable):
	"""Allows to select oligos from a group based on a "focus" window of
	interest. The window can be expanded to the closest oligo or to retain at
	least a given number of oligos."""

	_focus_window = None		# Left-close, right-open
	_oligos_in_focus_window = None
	_oligos_passing_score_filter = None

	def __init__(self, oligo_data, logger = logging.getLogger):
		super(OligoGroup, self).__init__(logger)
		self._data = pd.concat(oligo_data, ignore_index = True)
		self._oligos_passing_score_filter = self._data['score'].values <= 1

	@property
	def data(self):
		return self._data.copy()

	@property
	def focus_window(self):
		return self._focus_window
	
	@property
	def oligos_in_focus_window(self):
		return self._oligos_in_focus_window
	
	@property
	def oligos_passing_score_filter(self):
		return self._oligos_passing_score_filter
	
	@property
	def usable_oligos(self):
		return np.logical_and(
			self.oligos_in_focus_window,
			self.oligos_passing_score_filter)
	
	def get_focused_oligos(self, onlyUsable = False):
		if onlyUsable:
			return self._data.loc[self.usable_oligos]
		else:
			return self._data.loc[self.oligos_in_focus_window, :]

	def get_n_focused_oligos(self, onlyUsable = False):
		if isinstance(self.oligos_in_focus_window, type(None)):
			return 0
		return self.get_focused_oligos(onlyUsable).shape[0]

	def focus_all(self):
		# Use all oligos to define a focus region
		self.set_focus_window(self._data.loc[:,'start'].min(),
			self._data.loc[:,'start'].max()+1)

	def reset_focus_window(self):
		if isinstance(self.focus_window, type(None)):
			return
		(start, end) = self.focus_window
		start_condition = self._data['start'].values >= self.focus_window[0]
		end_condition = self._data['end'].values < self.focus_window[1]
		self._oligos_in_focus_window = np.logical_and(
			start_condition, end_condition)

	def set_focus_window(self, start, end, verbose = True):
		# Set a sub-window of interest to focus on
		self._focus_window = (start, end)
		start_condition = self._data['start'].values >= self.focus_window[0]
		end_condition = self._data['end'].values < self.focus_window[1]
		self._oligos_in_focus_window = np.logical_and(
			start_condition, end_condition)

		if verbose:
			nOligos = self.get_n_focused_oligos()
			nOligosUsable = self.get_n_focused_oligos(True)
			self.log.info(f"Set focus region to {self.focus_window}" +
				f" ({nOligos} oligos, {nOligosUsable} usable)")

	def expand_focus_to_n_oligos(self, n, verbose = True):
		# Expand the sub-window of interest to retrieve at least n oligos
		assert not isinstance(self.focus_window, type(None))

		if n <= self.get_n_focused_oligos():
			return

		for i in range(n-self.get_n_focused_oligos()):
			if not self.expand_focus_to_closest():
				return

		if verbose:
			nOligos = self.get_n_focused_oligos()
			nOligosUsable = self.get_n_focused_oligos(True)
			self.log.info(f"Expanded focus region to {self.focus_window}" +
				f" ({nOligos} oligos, {nOligosUsable} usable)")

	def expand_focus_by_step(self, step, verbose = True):
		# Expand the current focus window of a given step (in nt)
		assert 0 < step

		if self.focus_window[0] <= self._data['start'].min():
			if self.focus_window[1] >= self._data['end'].max():
				self.log.warning("Cannot expand the focus region any further " +
					"(all oligos already included)")
				return False

		new_focus_start, new_focus_end = self.focus_window
		new_focus_start -= step
		new_focus_end += step
		new_focus_start = np.max([new_focus_start, self._data['start'].min()])
		new_focus_end = np.min([new_focus_end, self._data['end'].max()])

		self.set_focus_window(new_focus_start, new_focus_end, verbose)
		return True

	def expand_focus_to_closest(self):
		# Expand the sub-window of interest to add the closest oligo
		# Return False if not possible (e.g., all oligos already included)
		if self.get_n_focused_oligos() == self._data.shape[0]:
			self.log.warning("Cannot expand the focus region any further " +
				"(all oligos already included)")
			return False

		earl = self._data['start'].values < self.focus_window[0]
		if 0 != earl.sum():
			max_start = self._data.loc[earl, 'start'].max()
			d_earl = self.focus_window[0] - max_start
		else:
			d_earl = np.inf

		late = self._data['end'].values >= self.focus_window[1]
		if 0 != late.sum():
			min_end = self._data.loc[late, 'end'].min()
			d_late = min_end - self.focus_window[1]
		else:
			d_late = np.inf

		if np.isinf(d_late):
			if np.isinf(d_earl):
				return False
			else:
				self.set_focus_window(max_start, self.focus_window[1], False)
		elif np.isinf(d_earl):
			self.set_focus_window(self.focus_window[0], min_end+1, False)
		else:
			if d_earl <= d_late:
				self.set_focus_window(max_start, self.focus_window[1], False)
			else:
				self.set_focus_window(self.focus_window[0], min_end+1, False)

		return True

	def apply_threshold(self, threshold):
		# Unfocuses oligos with score higher than the threshold
		assert threshold <= 1 and threshold >= 0
		self.reset_focus_window()
		self._oligos_passing_score_filter = self._data['score'] <= threshold

class OligoProbe(object):
	"""docstring for OligoProbe"""

	def __init__(self, oligo_data):
		super(OligoProbe, self).__init__()
		self.data = oligo_data

	@property
	def data(self):
		return self._data.copy()
	
	@data.setter
	def data(self, oligo_data):
		assert isinstance(oligo_data, pd.DataFrame)
		required_columns = ["start", "end", "Tm"]
		assert all([col in oligo_data.columns for col in required_columns])
		self._data = oligo_data
		self._range = (self._data['start'].min(), self._data['end'].max())
		self._size = self._range[1] - self._range[0]
		oDists = self._data['start'].values[1:] - self._data['end'].values[:-1]
		self._spread = np.std(oDists)
		self._d_range = (oDists.min(), oDists.max())
		self._tm_range = self._data['Tm'].max() - self._data['Tm'].min()

	@property
	def range(self):
		return self._range

	@property
	def tm_range(self):
		return self._tm_range
	
	@property
	def d_range(self):
		return self._d_range

	@property
	def spread(self):
		return self._spread
	
	@property
	def size(self):
		return self._size

	@property
	def featDF(self):
		df = pd.DataFrame([self.range[0], self.range[1], self._data.shape[0],
			self.size, self.spread, self.d_range[0], self.d_range[1],
			self.tm_range]).transpose()
		df.columns = ["start", "end", "size", "nOligos",
			"spread", "d_min", "d_max", "tm_range"]
		return df

	def __repr__(self):
		rep  = f"<OligoProbe[{self.range[0]}:{self.range[1]}"
		rep += f":{self.size}:{self.spread}]>"
		return rep

# FUNCTIONS ====================================================================

def assert_type(x, stype, label):
	assert isinstance(x,stype), f"{label} should be {stype}, {type(x)} instead"

def assert_multiTypes(x, types, label):
	cond = any([isinstance(x, t) for t in types])
	assert cond, f"{label} should be one of {types}, {type(x)} instead"

def assert_nonNeg(x, label):
	assert x > 0, f"{label} should be greater than 0"

def assert_inInterv(x, vmin, vmax, label, leftClose = False, rightClose = True):
	if leftClose:
		if rightClose:
			assert x >= vmin and x <= vmax, f"expected {vmin}<={label}<={vmax}"
		else:
			assert x >= vmin and x < vmax, f"expected {vmin}<={label}<{vmax}"
	else:
		if rightClose:
			assert x > vmin and x <= vmax, f"expected {vmin}<{label}<={vmax}"
		else:
			assert x > vmin and x < vmax, f"expected {vmin}<{label}<{vmax}"

# RUN ==========================================================================

if not reuse:
	assert not os.path.isfile(out_path)
	assert not os.path.isdir(out_path)
else:
	assert not os.path.isfile(out_path)
	if not os.path.isdir(out_path):
		os.mkdir(out_path)

logFormatter = logging.Formatter('%(asctime)s %(levelname)s:\t%(message)s',
	datefmt='%d/%m/%Y %I:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
consoleHandler.setLevel(logging.DEBUG)
logger.addHandler(consoleHandler)
logger = Loggable(logger)

logPath = "{0}/{1}.log".format(out_path, "test")
logger.addFileHandler(logPath)

logger.log.info(f"This log is saved at '{logPath}'.")

oWalker = OligoWalker(db_path, logger.log)
oWalker.out_path = out_path
oWalker.reuse = reuse
oWalker.start()

logging.shutdown()

# END ==========================================================================

################################################################################
