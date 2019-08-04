"""
Class representing a collection of features on a sequence.
"""

import operator

import grimoire.toolbox as toolbox
from grimoire.feature import Feature, mRNA, Gene

class FeatureTableError(Exception):
	pass

class FeatureTable:
	"""A collection of features from the same parent DNA"""
	
	def __init__(self, features=None):
		"""
		Parameters
		----------
		+ features `list` of `Feature` objects (optional)
		"""
		
		if features:
			self.features = features
		else:
			self.features = []

	def add_feature(self, feature):
		"""
		Parameters
		----------
		+ feature - an object of type `Feature`
		"""
		
		self.features.append(feature)

	def build_genes(self):
		"""
		Returns a list of `Gene` objects. Currently, only `mRNA` objects
		are constructed from features of type 'exon' and 'CDS' but in
		the future, there ought to be `ncRNA` objects as well. See the
		`Gene`, `Transcript`, `mRNA`, and `ncRNA` objects.

	
		In order for the exon/CDS features to be grouped into transcripts and then into genes, the features must use the `Feature`
		id and pid attributes correctly.
		"""
	
		genes = {}
		mRNAs = {}
		parts = []
		for f in self.features:
			if f.type == 'gene':
				if f.id == None: raise FeatureError('genes need ids')
				genes[f.id] = Gene(f.dna, f.beg, f.end, f.strand, f.type,
					id=f.id)
			elif f.type == 'mRNA':
				if f.id == None: raise FeatureError('mRNAs need ids')
				if f.pid == None: raise FeatureError('mRNAs need pids')
				mRNAs[f.id] = mRNA(f.dna, f.beg, f.end, f.strand, f.type,
					id=f.id, pid=f.pid)
			elif f.type == 'exon' or f.type == 'CDS':
				if f.pid == None: raise FeatureError('need parent')
				parts.append(f)
	
		# add parts to mRNAs
		for f in parts:
			for pid in f.pid:
				if pid in mRNAs:
					mRNAs[pid].add_child(f)

		# add mRNAs to genes
		for txid in mRNAs:
			f = mRNAs[txid]
			if len(f.pid) != 1: raise FeatureError('mRNA pids != 1')
			pid = f.pid[0]
			if pid in genes:
				genes[pid].add_child(f)
	
		# perform sanity checks on all genes
		for gid in genes:
			genes[gid].validate()
	
		return list(genes.values())

	def compare(self, other):
		"""
		Compares this feature table with some other feature table,
		returning a dictionary of comparision values.
		"""
				
		# NT-level comparisons
		same, diff = 0, 0
		s, p = [], []
		for f in self.features:
			for i in range(f.beg, f.end + 1): s.append(f.type)
		for f in other.features:
			for i in range(f.beg, f.end + 1): p.append(f.type)
		for i in range(len(s)):
			if s[i] == p[i]: same += 1
			else:            diff += 1

		# Exact-level comparisons
		exact, inexact = 0, 0
		if diff == 0: exact = 1
		else:         inexact = 1

		# Feature-type-level comparisons
		matrix = {}
		for i in range(len(s)):
			if s[i] not in matrix: matrix[s[i]] = {}
			if p[i] not in matrix[s[i]]: matrix[s[i]][p[i]] = 0
			matrix[s[i]][p[i]] += 1

		stats = {
			'nt_same' : same,
			'nt_diff' : diff,
			'exact' : exact,
			'inexact' : inexact,
			'matrix' : matrix
		}
		
		return stats

