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

def gff_to_feature(dna, gff):
	"""Converts a GFF object into a Feature object."""
	
	attr = gff.attr.rstrip()
	
	id = None
	im = re.search('ID=([^;]+)', attr)
	if im:
		id = im[1].rstrip()
	
	pid = None
	pm = re.search('Parent=([^;]+)', attr)			
	if pm:
		if ',' in pm[1]: pid = pm[1].split(',')
		else:            pid = pm[1]
			
	return Feature(dna, gff.beg, gff.end, gff.strand, gff.type,
		source=gff.source, score=gff.score, id=id, pid=pid)

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

		self._gff = io.GFF_file(gff)
		self._fasta = io.FASTA_stream(fasta)
		self._check = check

	def __iter__(self):
		return self

	def __next__(self):
		return self.next()

	def next(self):
		"""
		Retrieves the next entry of the FASTA file as DNA with features bound.
		"""
		try:
			entry = next(self._fasta)
			dna = DNA(name=entry.id, seq=entry.seq)
			if self._check: dna.check_alphabet()
			
			# add features
			for g in self._gff.get(chrom=dna.name):
				dna.ftable.add_feature(gff_to_feature(dna, g))
			
			return dna
			
		except StopIteration:
			raise StopIteration()

