
'''
@author: Gabriele Girelli
@contact: gigi.ga90@gmail.com
@description: methods for oligo manipulation.
'''

import configparser as cp
import itertools
import logging
import multiprocessing as mp
import numpy as np
import os
import pandas as pd
from pathlib import Path
import shutil
import sys
from tqdm import tqdm

from ifpd2.logging import Loggable
from ifpd2.asserts import *

class Oligo(object):
	"""Oligo database line/record parser.
	Presents oligo values as properties."""

	__colnames = ['name', 'chrom', 'start', 'end', 'tm_dG',
		'dH', 'dS', 'Tm', 'seq', 'nOT', 'k', 'ss_dG', 'GC']

	def __init__(self, oligo, i):
		super(Oligo, self).__init__()
		assert_type(i, int, "oligo id")
		assert i >= 0
		oligo = oligo.strip().split("\t")
		for i in [2, 3, 9, 10]:
			oligo[i] = int(oligo[i])
		for i in [4, 5, 6, 7, 11, 12]:
			oligo[i] = float(oligo[i])
		self.data = pd.DataFrame(oligo,
			columns=[i], index=self.__colnames).transpose()

	@property
	def start(self):
		return self.data['start'].values[0]

	@property
	def end(self):
		return self.data['end'].values[0]

	@property
	def score(self):
		if "score" in self.data.columns:
			return self.data['score'].values[0]
		return np.nan
	
	@staticmethod
	def __norm_value_in_range(v, r):
		if 0 == r[1]: return np.nan
		return (v - r[0])/(r[1] - r[0])

	def add_score(self, F, Gs):
		nOT = self.data['nOT'].values
		if nOT <= F[0]:
			self.data['score'] = 0
			return
		if nOT > F[1]:
			self.data['score'] = np.inf
			return
		norm_nOT = self.__norm_value_in_range(nOT, F)

		ss_dG = self.data['ss_dG'].values
		if isinstance(Gs[0], int):
			if ss_dG >= Gs[0]: score = 0
			if ss_dG < Gs[1]: score = np.inf
			norm_ss_dG = self.__norm_value_in_range(ss_dG, Gs)
		else:
			tm_dG = self.data['tm_dG'].values
			if ss_dG >= tm_dG*Gs[0]: score = 0
			if ss_dG < tm_dG*Gs[1]: score = np.inf
			norm_ss_dG = self.__norm_value_in_range(ss_dG,
				[tm_dG*f for f in Gs])
		self.data['score'] = np.mean([norm_nOT, norm_ss_dG])

class OligoGroup(Loggable):
	"""Allows to select oligos from a group based on a "focus" window of
	interest. The window can be expanded to the closest oligo or to retain at
	least a given number of oligos."""

	_focus_window = None		# Left-close, right-open
	_oligos_in_focus_window = None
	_oligos_passing_score_filter = None

	def __init__(self, oligos, logger = logging.getLogger()):
		super(OligoGroup, self).__init__(logger)
		self._data = pd.concat([o.data for o in oligos], ignore_index = True)
		self._data = self._data.loc[self._data['score'] <= 1, :]
		self._oligos_passing_score_filter = self._data['score'].values <= 1

	@property
	def data(self):
		return self._data.copy()

	@property
	def focus_window(self):
		return self._focus_window

	@focus_window.setter
	def focus_window(self, focus_window):
		assert_type(focus_window, tuple, "focus window")
		assert 2 == len(focus_window)
		assert focus_window[1] > focus_window[0]
		self._focus_window = focus_window

	@property
	def focus_window_repr(self):
		return f"[{self.focus_window[0]}:{self.focus_window[1]})"
	
	@property
	def focus_window_size(self):
		return self.focus_window[1] - self.focus_window[0]

	@property
	def oligos_in_focus_window(self):
		return self._oligos_in_focus_window

	@property
	def oligos_passing_score_filter(self):
		return self._oligos_passing_score_filter
	
	@property
	def usable_oligos(self):
		if isinstance(self.oligos_in_focus_window, type(None)):
			return self.oligos_passing_score_filter
		return np.logical_and(
			self.oligos_in_focus_window,
			self.oligos_passing_score_filter)
	
	def get_focused_oligos(self, onlyUsable = False):
		if onlyUsable:
			return self._data.loc[self.usable_oligos]
		else:
			return self._data.loc[self.oligos_in_focus_window, :]

	def get_n_focused_oligos(self, onlyUsable = False):
		if 0 == self.usable_oligos.sum():
			return 0
		return self.get_focused_oligos(onlyUsable).shape[0]

	def focus_all(self):
		# Use all oligos to define a focus region
		self.set_focus_window(self._data.loc[:,'start'].min(),
			self._data.loc[:,'start'].max()+1)

	def set_focus_window(self, start, end, verbose = True):
		# Set a sub-window of interest to focus on
		self.focus_window = (start, end)

		start_condition = self._data['start'].values >= self.focus_window[0]
		end_condition = self._data['end'].values < self.focus_window[1]
		self._oligos_in_focus_window = np.logical_and(
			start_condition, end_condition)

		if verbose:
			nOligos = self.get_n_focused_oligos()
			nOligosUsable = self.get_n_focused_oligos(True)
			self.log.info(f"Set focus region to {self.focus_window_repr}" +
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
		new_focus_start -= step/2
		new_focus_end += step/2
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
		self._oligos_passing_score_filter = self._data['score'] <= threshold

	def reset_threshold(self):
		self.apply_threshold(1)

	def discard_focused_oligos_safeDist(self, safeDist):
		# Discard focused oligos that are not within a safe distance from the
		# CFR borders
		start = self.data.loc[self.oligos_in_focus_window, "start"].min()
		end = self.data.loc[self.oligos_in_focus_window, "end"].max()+1
		self.__discard_oligos_in_range(start+safeDist, end-safeDist)	

	def discard_focused_oligos_safeN(self, safeN, D):
		# Discard focused oligos that are neither the first nor the last safeN
		
		start = self.data.loc[self.oligos_in_focus_window,"start"
			].values[0] + len(self.data.loc[self.oligos_in_focus_window,"seq"
			].values[0]) + D
		end = self.data.loc[self.oligos_in_focus_window,"end"
			].values[-1] + len(self.data.loc[self.oligos_in_focus_window,"seq"
			].values[-1]) + D

		oData = self.get_focused_oligos()

		c = 1
		while c <= safeN:
			passing_oData = oData.loc[oData['start'] > start, :]
			if 0 == passing_oData.shape[0]:
				self.log.info("Not enough oligos, skipped discard step.")
				return
			start = passing_oData['start'].values[0] + len(
				passing_oData['seq'].values[0]) + D
			c += 1

		c = 1
		while c <= safeN:
			passing_oData = oData.loc[oData['end'] <= start]
			if 0 == passing_oData.shape[-1]:
				self.log.info("Not enough oligos, skipped discard step.")
				return
			end = passing_oData['end'].values[-1] - len(
				passing_oData['seq'].values[-1]) - D
			c += 1

		if not self.__discard_oligos_in_range(start, end):
			self.log.info(f"No oligos to discard in range [{start}:{end}).")

	def __discard_oligos_in_range(self, start, end):
		if start >= end: return False
		start_condition = self.data['end'] < start
		end_condition = self.data['start'] >= end
		keep_condition = np.logical_or(start_condition, end_condition)
		nDiscarded = self.oligos_in_focus_window.sum()-keep_condition.sum()
		self.__apply_keep_condition(keep_condition)
		self.log.info(f"Discarded {nDiscarded}" +
			f" oligos from the [{start}:{end}) range.")
		return True

	def __apply_keep_condition(self, keep_condition):
		self._data = self._data.loc[keep_condition, :]
		self._oligos_in_focus_window = self._oligos_in_focus_window[
			keep_condition]
		self._oligos_passing_score_filter = self._oligos_passing_score_filter[
			keep_condition]

class OligoProbe(object):
	"""Converts a DataFrame of oligo data into an OligoProbe."""

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
		self._d_mean = oDists.mean()
		self._tm_range = (self._data['Tm'].min(), self._data['Tm'].max())

	@property
	def n_oligos(self):
		return self._data.shape[0]

	@property
	def path(self):
		return self._data.index

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
	def d_mean(self):
		return self._d_mean

	@property
	def spread(self):
		return self._spread
	
	@property
	def size(self):
		return self._size

	@property
	def featDF(self):
		df = pd.DataFrame([self.range[0], self.range[1], self.n_oligos,
			self.size, self.spread,
			self.d_range[0], self.d_range[1], self.d_mean,
			np.diff(self.tm_range)[0]]).transpose()
		df.columns = ["start", "end", "nOligos", "size",
			"spread", "d_min", "d_max", "d_mean", "tm_range"]
		return df

	def __repr__(self):
		rep  = f"<OligoProbe[{self.range[0]}:{self.range[1]}"
		rep += f":{self.size}:{self.spread}]>"
		return rep

	def count_shared_oligos(self, probe):
		# Counts oligos shared with another probe
		# based on their paths
		return np.intersect1d(self.path, probe.path).shape[0]

	def export(self, path):
		assert not os.path.isfile(path)
		if os.path.isdir(path):
			shutil.rmtree(path)
		os.mkdir(path)

		self.featDF.to_csv(os.path.join(path, "probe.tsv"), "\t", index = False)
		self.data.to_csv(os.path.join(path, "oligos.tsv"), "\t", index = False)

		with open(os.path.join(path, "probe.bed"), "w+") as BH:
			for i in range(self.data.shape[0]):
				oligo = self.data.iloc[i]
				BH.write(f">{oligo['name']}:" +
					f"{oligo['chrom']}:{oligo['start']}-{oligo['end']}\n")
				BH.write(f"{oligo['seq']}\n")

class OligoPathBuilder(object):
	"""OligoPathBuilder contains methods to select non-overlapping sets of
	oligos (i.e., probes) from an oligo DataFrame in input."""

	N = int(48)			# Number of oligos per probe
	D = int(2)			# Minimum distance between consecutive oligos
	Tr = 10.0			# Melting temperature range half-width
	Ps = int(10000)		# Probe size threshold, in nt (Ps > 1)
	Ph = .1				# Maximum hole size in probe as fraction of probe size
	Po = .5				# Probe oligo intersection threshold for path reduction

	def __init__(self):
		super(OligoPathBuilder, self).__init__()

	def _assert(self):
		assert_type(self.N, int, "N")
		assert_nonNeg(self.N, "N")

		assert_type(self.D, int, "D")
		assert_nonNeg(self.D, "D")

		assert_type(self.Tr, float, "Tr")
		assert_nonNeg(self.Tr, "Tr")

		assert_type(self.Ps, int, "Ps")
		assert self.Ps > 1

		assert_type(self.Ph, float, "Ph")
		assert_inInterv(self.Ph, 0, 1, "Ph")

		assert_type(self.Po, float, "Po")
		assert_inInterv(self.Po, 0, 1, "Po")

	@staticmethod
	def get_non_overlapping_paths(oData, D):
		# Gets all paths of consecutive non-overlapping oligos with minimum
		# distance equal to D
		assert_type(oData, pd.DataFrame, "oData")

		edges = []
		start_positions = oData["start"].values
		end_positions = oData["end"].values

		for i in range(oData.shape[0]):
			edges.append(np.logical_or(
				end_positions+D < oData["start"].values[i],
				start_positions-D >= oData["end"].values[i]
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

			path_set.add(tuple(path))
			path_start_set.add(path[0])

		return path_set

	def filter_paths(self, path_set, oData):
		# Selects oligo paths based on length, melting temperature, size, and
		# presence of gaps.
		assert_type(oData, pd.DataFrame, "oData")

		exit_polls = {'P':0, 'N':0, 'S':0, 'H':0, 'T':0}

		sized_paths = set()
		for path in list(path_set):
			if path in sized_paths: continue
			if self.N > len(path):
				exit_polls['N'] = exit_polls['N'] + 1
				continue
			if self.N == len(path):
				sized_paths.add(path)
			else:
				for j in range(len(path) - self.N + 1):
					subpath = path[j:(j + self.N)]
					sized_paths.add(subpath)

		selected_paths = set()
		for path in list(sized_paths):
			if path in selected_paths: continue
			passed, comment = self.__path_passes(path, oData)
			exit_polls[comment] = exit_polls[comment] + 1
			if not passed: continue
			selected_paths.add(path)
		
		comment = "".join([f"{r}{c}" for (c,r) in exit_polls.items()])
		return (list(selected_paths), comment)

	def __path_passes(self, path, oData):
		if not isinstance(path, list): path = list(path)
		pData = oData.iloc[path, :]
		probe_size = pData['end'].values[-1]
		probe_size -= pData['start'].values[0]
		if self.Ps < probe_size:
			return (False, 'S')
		dData = pData['start'].values[1:] - pData['end'].values[:-1]
		max_hole_size = dData.max()
		if self.Ph * probe_size < max_hole_size:
			return (False, 'H')
		Tm_range = pData['Tm'].max()-pData['Tm'].min()
		if 2*self.Tr < Tm_range:
			return (False, 'T')
		return (True, 'P')

	@staticmethod
	def convert_paths_to_probes(path_list, oData):
		# Converts a list of paths into a list of probes
		assert_type(path_list, list, "path list")
		probe_list = []
		if 0 != len(path_list):
			for path in path_list:
				probe_list.append(OligoProbeBuilder.path2probe(path, oData))
		return probe_list

	@staticmethod
	def path2probe(path, oData):
		# Convert an oligo path into an OligoProbe
		return OligoProbe(oData.iloc[list(path), :])

class OligoProbeBuilder(OligoPathBuilder):
	"""docstring for OligoProbeBuilder"""

	k = 40
	F = (0, 100)		# Threshold on number of off-targets (range)
	Gs = (0.0, 0.5)		# dG of SS either as kcal/mol (negative)
						#  or as fraction dG of hybrid (0<=Gs<=1)
						#  (range)
	Ot = .1				# Step for oligo score relaxation

	def __init__(self):
		super(OligoProbeBuilder, self).__init__()

	@property
	def config(self):
		config = cp.ConfigParser()
		config['AIM'] = {
			'Oligo length (nt)' : self.k,
			'Oligo(s) number' : self.N
		}
		config['OLIGO FILTERS'] = {
			'Off-target threshold' : self.F,
			'Secondary structure dG threshold' : self.Gs,
			'Oligo score relaxation step' : self.Ot
		}
		config['PROBE FILTERS'] = {
			'Melting temperature range (degC)' : self.Tr,
			'Min. consecutive oligo distance (nt)' : self.D,
			'Probe size threshold' : self.Ps,
			'Maximum hole size' : self.Ph
		}
		return config

	def _assert(self):
		OligoPathBuilder._assert(self)

		assert_type(self.k, int, "k")
		assert_nonNeg(self.k, "k")
		
		assert (self.k+self.D)*self.N <= self.Ps
	
		assert_type(self.F, tuple, "F")
		assert 2 == len(self.F)
		for i in range(2):
			assert_type(self.F[i], int, f"F[{i}]")
			assert self.F[i] >= 0
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

	def get_prologue(self):
		s  = f"* OligoProbeBuilder *\n\n"
		s += f"Aim to build probes with {self.N} oligos each.\n"
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
		s += f"Reducing probes when oligo intersection fraction is equal to"
		s += f" or greater than {self.Po}.\n"

		return s

	def start(self, oGroup, window, cfr_step, logger):
		self._assert()
		return self.__build_probe_candidates(oGroup, window, cfr_step, logger)

	def __build_probe_candidates(self, oGroup, window, cfr_step, logger):
		# Applies oligo filters to the oligo group,
		# expands the focus group if needed, and build probe candidates
		
		if np.isnan(window['cfr_start']):
			oGroup.focus_all()
		else:
			oGroup.set_focus_window(window['cfr_start'], window['cfr_end'])
			oGroup.expand_focus_to_n_oligos(self.N)

		probe_list = self.__explore_filter(oGroup, logger)
		if not np.isnan(window['cfr_start']):
			while 0 == len(probe_list):
				oGroup.reset_threshold()
				if oGroup.focus_window_size >= self.Ps:
					oGroup.discard_focused_oligos_safeN(self.N-1, self.D)
				if not oGroup.expand_focus_by_step(cfr_step):
					break
				probe_list = self.__explore_filter(oGroup, logger)
		
		return probe_list

	def __explore_filter(self, oGroup, logger):
		# Explores the 0-to-max_score score threshold range and stops as soon as
		# one probe candidate passes all user-defined thresholds
		
		if oGroup.focus_window_size < self.Ps:
			max_score = 0
			logger.warning("Score relaxation deactivated when focus region" +
				" size is smaller than probe size threshold.")
		else:
			max_score = 1

		nOligos_in_focus_window = oGroup.get_n_focused_oligos(True)

		score_thr = 0
		oGroup.apply_threshold(score_thr)
		nOligos_prev_score_thr = oGroup.get_n_focused_oligos(True)
		logger.info(f"Set oligo score threshold at {score_thr:.3f}" +
			f" ({nOligos_prev_score_thr} oligos usable)...")

		if 0 == max_score and 0 == nOligos_prev_score_thr:
			return([])

		while 0 == nOligos_prev_score_thr and score_thr <= max_score-self.Ot:
			score_thr += self.Ot
			oGroup.apply_threshold(score_thr)
			nOligos_prev_score_thr = oGroup.get_n_focused_oligos(True)
			logger.info(f"Relaxed oligo score threshold to {score_thr:.3f}" +
				f" ({nOligos_prev_score_thr} oligos usable)...")

		probe_list = self.__get_non_overlapping_probes(
			oGroup.get_focused_oligos(True), logger)
		while 0 == len(probe_list):
			score_thr += self.Ot
			if score_thr > max_score:
				break

			oGroup.apply_threshold(score_thr)
			nOligos = oGroup.get_n_focused_oligos()
			nOligosUsable = oGroup.get_n_focused_oligos(True)
			if nOligosUsable == nOligos_prev_score_thr:
				continue
			if 0 == nOligosUsable:
				continue

			logger.info(f"Relaxed oligo score threshold to {score_thr:.3f}" +
				f" ({nOligosUsable} oligos usable)...")
			probe_list = self.__get_non_overlapping_probes(
				oGroup.get_focused_oligos(True), logger)

			if nOligosUsable == nOligos_in_focus_window:
				logger.warning(
					"All oligos included. Score relaxation ineffective.")
				break
			nOligos_prev_score_thr = nOligosUsable

		return(probe_list)

	def __get_non_overlapping_probes(self, oData, logger, verbosity = True):
		# Builds non overlapping oligo paths from an oligo data table, filters
		# them based on class attributes, and then converts the remaining ones
		# to OligoProbe objects.
		

		paths = self.get_non_overlapping_paths(oData, self.D)
		if 0 == len(paths):
			logger.warning("No oligo paths found.")
			return []

		pathMaxLen = np.max([len(p) for p in list(paths)])
		if verbosity:
			nPaths = len(paths)
			logger.info(f"Found {nPaths} sets with up to {pathMaxLen} " +
				f"non-overlapping oligos.")
		if pathMaxLen < self.N:
			logger.warning(f"Longest path is shorter than requested. " + 
				"Skipped.")
			return []

		paths, comment = self.filter_paths(paths, oData)
		if verbosity:
			logger.info(f"{len(paths)}/{nPaths} oligo paths remaining " +
				f"after filtering. ({comment})")

		return self.convert_paths_to_probes(paths, oData)

	def reduce_probe_list(self, probe_list):
		try:
			if 0 == len(probe_list): return []
			sorted_probes = sorted(probe_list, key=lambda p: p.range[0])

			selected_probes = []
			probe_ref = sorted_probes[0]
			for probe in sorted_probes[1:-1]:
				n_shared_oligos = probe_ref.count_shared_oligos(probe)

				if probe_ref.n_oligos == n_shared_oligos:
					raise Exception("Encountered probe duplicates!")

				if self.Po * self.N <= n_shared_oligos:
					probe_ref = self.select_probe_from_pair(probe_ref, probe)
				else:
					selected_probes.append(probe_ref)
					probe_ref = probe

			if not probe_ref in selected_probes:
				selected_probes.append(probe_ref)

			n_shared_oligos = probe_ref.count_shared_oligos(sorted_probes[-1])
			if self.Po * self.N > n_shared_oligos:
				selected_probes.append(sorted_probes[-1])

			return selected_probes
		except Exception as e:
			print(e)
			raise

	def select_probe_from_pair(self, probeA, probeB):
		if probeA.size < probeB.size:
			return probeA
		else:
			return probeB
		if np.diff(probeA.tm_range)[0] < np.diff(probeB.tm_range)[0]:
			return probeA
		else:
			return probeB
		if probeA.spread/probeA.d_mean < probeB.spread/probeB.d_mean:
			return probeA
		else:
			return probeB
		return probeA

	@staticmethod
	def import_probes(ipath):
		# Imports output from given directory path
		assert os.path.isdir(ipath)

		pPath = os.path.join(ipath, "probe_paths.tsv")
		if not os.path.isfile(pPath): return []

		oPath = os.path.join(ipath, "oligos.tsv")
		if not os.path.isfile(oPath): return []
		
		oligos = pd.read_csv(oPath, "\t", index_col = 0)
		paths = pd.read_csv(pPath, "\t", index_col = 0)

		probe_list = []
		for p in paths.iloc[:, 0]:
			probe_list.append(OligoProbe(oligos.loc[
				[int(o) for o in p.split(",")]]))

		return probe_list

	@staticmethod
	def export_probes(probe_list, opath):
		# Exports a list of probes to a directory
		assert os.path.isdir(opath)

		probe_df = pd.concat([p.featDF for p in probe_list], ignore_index=True)
		probe_df.sort_values("start", inplace = True)
		probe_df.to_csv(os.path.join(opath, "probe_feat.tsv"), "\t")

		pd.concat([p.data for p in probe_list]).drop_duplicates().to_csv(
			os.path.join(opath, "oligos.tsv"), "\t")

		probe_paths = []
		for pi in range(len(probe_list)):
			probe_paths.append([",".join([str(x)
				for x in probe_list[pi].data.index.tolist()])])
		probe_paths = pd.DataFrame(probe_paths)
		probe_paths.columns = ["cs_oligos"]
		probe_paths.to_csv(os.path.join(opath, "probe_paths.tsv"), "\t")

class OligoProbeSet(object):
	"""docstring for OligoProbeSet"""

	def __init__(self, probe_list):
		super(OligoProbeSet, self).__init__()
		self.probe_list = probe_list

	@property
	def probe_list(self):
		return self._probe_list

	@probe_list.setter
	def probe_list(self, probe_list):
		assert all([isinstance(p, OligoProbe) for p in probe_list])
		self._probe_list = sorted(probe_list, key = lambda p: p.range[0])
		self._probe_tm_ranges = [p.tm_range for p in self.probe_list]
		self._tm_range = (np.min([t[0] for t in self._probe_tm_ranges]),
			np.max([t[1] for t in self._probe_tm_ranges]))
		self._sizes = [p.size for p in self.probe_list]
		self._spreads = [p.spread for p in self.probe_list]
		self._probe_ranges = [p.range for p in self.probe_list]
		probe_starts = np.array([r[0] for r in self.probe_ranges])
		probe_ends = np.array([r[1] for r in self.probe_ranges])
		self._range = (probe_starts.min(), probe_ends.max())
		self._ds = probe_starts[1:] - probe_ends[:-1]
		self._d_mean = np.mean(self.ds)
		self._d_range = (np.min(self.ds), np.max(self.ds))
		tm = [p.data['Tm'].tolist() for p in self.probe_list]
		self._oligo_tm = list(itertools.chain(*tm))

	@property
	def range(self):
		return self._range

	@property
	def probe_ranges(self):
		return self._probe_ranges

	@property
	def ds(self):
		return self._ds

	@property
	def d_mean(self):
		return self._d_mean

	@property
	def d_range(self):
		return self._d_range
	
	@property
	def tm_range(self):
		return self._tm_range

	@property
	def probe_tm_ranges(self):
		return self._probe_tm_ranges
	
	@property
	def oligo_tm(self):
		return self._oligo_tm
	
	@property
	def sizes(self):
		return self._sizes

	@property
	def spreads(self):
		return self._spreads

	@property
	def featDF(self):
		df = pd.DataFrame([self.range[0], self.range[1],
			len(self.probe_list), np.diff(self.range)[0],
			np.mean(self.sizes), np.std(self.sizes),
			np.min(self.spreads), np.max(self.spreads),
			np.mean(self.spreads), np.std(self.spreads),
			self.d_range[0], self.d_range[1], self.d_mean, np.std(self.ds),
			np.diff(self.tm_range)[0],
			np.mean(self.oligo_tm), np.std(self.oligo_tm),
			1/3*(np.std(self.sizes)/np.mean(self.sizes) +
				np.std(self.spreads)/np.mean(self.spreads) +
				np.std(self.oligo_tm)/np.mean(self.oligo_tm))
		]).transpose()
		df.columns = ["start", "end",
			"nProbes", "size", "size_mean", "size_std",
			"spread_min", "spread_max", "spread_mean", "spread_std",
			"d_min", "d_max", "d_mean", "d_std",
			"tm_range", "tm_mean", "tm_std", "score"]
		return df

	def export(self, path):
		assert not os.path.isfile(path)
		if os.path.isdir(path):
			shutil.rmtree(path)
		os.mkdir(path)
		
		self.featDF.to_csv(os.path.join(path, "set.tsv"), "\t", index = False)

		pd.concat([p.featDF for p in self.probe_list]
			).reset_index(drop = True).to_csv(
			os.path.join(path, "probes.tsv"), "\t")

		with open(os.path.join(path, "set.bed"), "w+") as BH:
			for pi in range(len(self.probe_list)):
				probe = self.probe_list[pi]
				probe.export(os.path.join(path, f"probe_{pi}"))
				for i in range(probe.data.shape[0]):
					oligo = probe.data.iloc[i]
					BH.write(f">probe_{pi}:{oligo['name']}" +
						f":{oligo['chrom']}:{oligo['start']}-{oligo['end']}\n")
					BH.write(f"{probe.data.iloc[i]['seq']}\n")

class OligoProbeSetBuilder(Loggable):
	"""docstring for OligoProbeSetBuilder"""

	probe_set_list = []

	def __init__(self, out_path, logger = logging.getLogger()):
		Loggable.__init__(self, logger)
		assert not os.path.isfile(out_path)
		self.out_path = out_path
		if os.path.isdir(out_path):
			shutil.rmtree(out_path)
		os.mkdir(out_path)

	def build(self, probe_candidates):
		nProbes = []
		for (wSet, window_list) in probe_candidates.items():
			probe_set_list = [(p,) for p in list(window_list.values())[0]]
			#nProbes = [len(probes) for probes in window_list.values()]
			#nProbeSets = np.sum([n for n in nProbes if 0 != n])

			for w in list(window_list.values())[1:]:
				if 0 == len(w): continue

				current_probe_set_list = list(tuple(probe_set_list))
				probe_set_list = []

				for probe in w:
					for probe_set in current_probe_set_list:
						current_probe_set = list(probe_set)
						current_probe_set.append(probe)
						probe_set_list.append(current_probe_set)

			self.probe_set_list.extend(probe_set_list)
			self.log.info(f"Built {len(probe_set_list)} probe set candidates " +
				f"from window set #{wSet+1}")

		self.probe_set_list = [OligoProbeSet(ps)
			for ps in self.probe_set_list]
		self.log.info(f"Built {len(self.probe_set_list)} " +
			"probe set candidates in total.")

		pd.concat([ps.featDF for ps in self.probe_set_list]
			).reset_index(drop = True).to_csv(
			os.path.join(self.out_path, "probe_sets.tsv"), "\t")

	def export(self):
		self.log.info("Exporting probe sets...")
		for psi in tqdm(range(len(self.probe_set_list)), desc = "Probe set"):
			probe_set = self.probe_set_list[psi]
			probe_set.export(os.path.join(self.out_path, f"probe_set_{psi}"))

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
		self._window_sets.reset_index(inplace = True)

	@property
	def wid(self):
		return self._w

	@property
	def current_window(self):
		return self.window_sets.iloc[self.wid, :]

	@property
	def reached_last_window(self):
		return self._reached_last_window

	def _assert(self):
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

	def _init_windows(self):
		# Build windows and central focus regions (CFR)
		if not os.path.isdir(os.path.join(self.out_path, "window_sets")):
			os.mkdir(os.path.join(self.out_path, "window_sets"))

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

		window_sets = []
		for i in range(len(np.arange(0, 1, self.Wh))):
			winID = np.array(range(nWindows)).reshape((nWindows, 1))
			setID = np.repeat(i, nWindows).reshape((nWindows, 1))
			window_sets.append(np.hstack([
				window_borders+i*self.Wh*self.Ws,
				window_mids+i*self.Wh*self.Ws,
				central_regions+i*self.Wh*self.Ws,
				winID, setID])[ :, (0, 2, 1, 3, 4, 5, 6)])

		self._w = 0 # Window ID
		self._reached_last_window = False
		self.window_sets = window_sets

	def go_to_next_window(self):
		if self.wid < self.window_sets.shape[0]-1:
			self._w += 1
		if self.wid == self.window_sets.shape[0]-1:
			self._reached_last_window = True

class OligoWalker(GenomicWindowSet, Loggable):
	"""OligoWalker walks through the oligos stored in an ifpd2 database,
	assigns them to windows based on user-defined parameters, and then builds
	probe candidates."""

	db_path = None
	out_path = "."
	reuse = False
	_threads = 1

	__current_oligos = []
	__walk_results = {}

	def __init__(self, db_path, logger = logging.getLogger()):
		GenomicWindowSet.__init__(self)
		Loggable.__init__(self, logger)
		self.db_path = db_path

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
			'Region' : f"{self.C}:{self.S}-{self.E}",
			'Probe(s) number' : self.X,
		}
		config['WINDOWS'] = {
			'Window size' : self.Ws,
			'Window step' : self.Wh,
			'Focus region size' :  self.Rs,
			'Focus region step' : self.Rt
		}
		return config

	def _assert(self):
		GenomicWindowSet._assert(self)

		assert os.path.isfile(self.db_path)
		assert os.path.isdir(self.out_path)

	def start(self, fparse, fimport, fprocess, fpost, *args, **kwargs):
		self._assert()
		self._init_windows()
		self.print_prologue()
		self.__walk(fparse, fimport, fprocess, fpost, *args, **kwargs)

	def print_prologue(self):
		nWindows = self.window_sets['w'].max().astype('i')+1
		nSets = self.window_sets["s"].max().astype('i')+1

		s  = f"* OligoWalker *\n\n"
		s += f"Threads: {self.threads}\n"
		s += f"Database: '{self.db_path}'\n"
		s += f"Region of interest: {self.C}:{self.S}-{self.E}\n"
		s += f"Aim to scout {nWindows} windows.\n\n"

		s += f"Using a central focus region of {self.Rs} nt,"
		s += f" in windows of size {self.Ws} nt,\n"
		s += f"built with a shift of {self.Ws*self.Wh} nt ({self.Wh*100}%).\n"
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
		DBH = open(self.db_path, "r")
		next(DBH)

		self.__preprocess_window(fimport)
		self.__load_windows_until_next_to_do(fimport)
		
		if self.reached_last_window and self.__window_done():
			self.log.info("All windows pre-processed. Skipped database walk.")
			return

		self.r = 0 # Record ID
		self.rw = 0 # Walk step counter

		exec_results = []

		DBHpb = tqdm(DBH, leave = None, desc = "Parsing records")
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
						break

					self.go_to_next_window()
					self.__preprocess_window(fimport)
					self.__load_windows_until_next_to_do(fimport)
					if self.reached_last_window and self.__window_done():
						break

					if 0 != len(self.current_oligos):
						self.remove_oligos_starting_before(
							self.current_window['start'])

				if oligo_end > self.E:	# End reached
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
				if s in self.walk_results.keys():
					self.walk_results[s][w] = results
				else:
					self.walk_results[s] = {w : results}

		DBH.close()
		pool.close()

	def __window_done(self):
		s = int(self.current_window['s'])
		w = int(self.current_window['w'])
		if s in self.walk_results.keys():
			if w in self.walk_results[s].keys():
				return isinstance(self.walk_results[s][w], list)
		return False

	def __load_windows_until_next_to_do(self, fimport):
		while self.__window_done() and not self.reached_last_window:
			self.go_to_next_window()
			self.__preprocess_window(fimport)

	def __preprocess_window(self, fimport):
		# Preprocess current window
		# Triggers import if previously run and matching current window

		if not os.path.isdir(self.window_set_path):
			os.mkdir(self.window_set_path)
			self.window_sets.loc[
				self.window_sets['s'] == self.current_window['s'],:].to_csv(
				os.path.join(self.window_set_path, "windows.tsv"), "\t")

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
		results = OligoWalker.process_window(oligos, window,
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
			return
		
		status, results = fpost(results, opath, *args, **kwargs)
		if status and not isinstance(opath, type(None)):
			Path(os.path.join(opath, ".done")).touch()

		return (int(window['s']), int(window['w']), results)