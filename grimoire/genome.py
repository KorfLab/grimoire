"""
Classes for handling genomes with their attached annotation.
"""

import re
import json
import sys
import operator
import gzip

import grimoire.io as io
from grimoire.sequence import DNA
from grimoire.feature import Feature, Gene, mRNA, ncRNA

def build_genes(dna):
	"""
	Returns gene objects from a `DNA` object using its feature table. `Gene`
	objects contain `Transcript` objects which may be of subclass `mRNA` or
	`ncRNA` depending on if they are coding. If the feature table has both
	exon and CDS objects, `build_genes()` will attempt to construct `mRNA`
	objects (otherwise `ncRNA`).
	
	In order for the exon/CDS features to be grouped into transcripts and
	then into genes, the GFF must properly specify ID and Parent_ID as in
	the GFF3 spec.
	
	Parameter
	---------
	+ dna `grimoire.DNA` DNA object with features
	"""
	
	genes = {}
	mRNAs = {}
	parts = []
	for f in dna.features:
		if f.type == 'gene':
			if f.id == None:
				raise FeatureError('genes need ids')
			genes[f.id] = Gene(dna, f.beg, f.end, f.strand, f.type, id=f.id)
		elif f.type == 'mRNA':
			if f.id == None:
				raise FeatureError('mRNAs need ids')
			if f.pid == None:
				raise FeatureError('mRNAs need pids')
			mRNAs[f.id] = mRNA(dna, f.beg, f.end, f.strand, f.type,
				id=f.id, pid=f.pid)
		elif f.type == 'exon' or f.type == 'CDS':
			if f.pid == None:
				raise FeatureError('exons and CDSs need parents')
			parts.append(f)
	
	# add parts to mRNAs
	for f in parts:
		for pid in f.pid:
			if pid in mRNAs:
				mRNAs[pid].add_child(f)

	# add mRNAs to genes
	for txid in mRNAs:
		f = mRNAs[txid]
		if len(f.pid) != 1:
			raise FeatureError('mRNAs should have exactly one parent id')
		pid = f.pid[0]
		if pid in genes:
			genes[pid].add_child(f)
	
	# perform sanity checks on all genes
	for gid in genes:
		genes[gid].validate()
	
	return list(genes.values())

class GenomeError(Exception):
	pass

class Reader:
	"""Class for iterating through DNA objects with attached feature tables."""

	def __init__(self, fasta=None, gff=None, check=False):
		"""
		Parameters
		----------
		+ fasta= `str`  path to fasta file (may be compressed)
		+ gff=  `str`  path to gff file (may be compressed)
		+ check= `bool` check that the alphabet conforms to IUPAC DNA
		"""

		self._fp = None
		self._gz = False
		self._gff = io.GFF_file(gff)
		self._check = check

		if re.search('\.gz$', fasta):
			self._fp = gzip.open(fasta)
			self._gz = True
		else:
			self._fp = open(fasta, 'r')
		self._lastline = ''
		self._done = False

	def __iter__(self):
		return self

	def __next__(self):
		return self.next()

	def next(self):
		"""
		Retrieves the next entry of the FASTA file as DNA with features bound.
		"""

		if self._done: raise StopIteration()
		header = None
		if self._lastline[0:1] == '>':
			header = self._lastline
		else:
			header = self._fp.readline()
			if self._gz: header = str(header, 'utf-8')

		m = re.search('>\s*(\S+)\s*(.*)', header)
		id = m[1]
		desc = m[2]
		seq = []

		while (True):
			line = self._fp.readline()
			if self._gz: line = str(line, 'utf-8')
			if line[0:1] == '>':
				self._lastline = line
				break
			if line == '':
				self._done = True
				self._fp.close()
				break

			line = line.replace(' ', '')
			seq.append(line.strip())
		
		dna = DNA(name=id, seq=''.join(seq))
		if self._check: dna.check_alphabet()

		# add features
		for g in self._gff.get(chrom=dna.name):
			attr = g.attr.rstrip()
			
			id = None
			im = re.search('ID=([^;]+)', attr)
			if im:
				id = im[1].rstrip()
			
			pid = None
			pm = re.search('Parent=([^;]+)', attr)			
			if pm:
				if ',' in pm[1]: pid = pm[1].split(',')
				else:            pid = pm[1]
			
			dna.features.append(Feature(dna, g.beg, g.end, g.strand,
				g.type, source=g.source, score=g.score, id=id, pid=pid))

		return dna

