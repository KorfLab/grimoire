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

class GenomeError(Exception):
	pass

class Reader:
	"""Class for iterating through DNA objects with attached feature tables."""

	def __init__(self, fasta=None, gff=None, check=False):
		"""
		Parameters
		----------
		+ fasta= `str`  path to fasta file (may be compressed)
		+ gff=   `str`  path to gff (or other file, may be compressed)
		+ check= `bool` check that the alphabet conforms to IUPAC DNA
		"""

		self._fp = None
		self._gz = False
		self._gff = io.GFF_file(gff)
		self._check = check

		if re.search(r'\.gz$', fasta):
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

		m = re.search(r'>\s*(\S+)\s*(.*)', header)
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
			
			dna.ftable.add_feature(Feature(dna, g.beg, g.end, g.strand,
				g.type, source=g.source, score=g.score, id=id, pid=pid))

		return dna

