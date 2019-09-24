"""
Classes for represnting biological sequences.
"""
import grimoire.toolbox as toolbox
from grimoire.feature import FeatureTable

class SequenceError(Exception):
	pass

class BioSequence:
	"""Parent base class for all biological sequences."""

	def fasta(self, wrap=80):
		"""
		Returns the object as a string in FASTA format.

		Parameters
		----------
		+ wrap= `int` number of characters per line
		"""

		s = '>'
		if self.name: s += self.name
		if self.desc: s += ' ' + self.desc
		s += '\n'
		for i in range(0, len(self.seq), wrap):
			s += self.seq[i:i+wrap] + '\n'
		return s

	def __str__(self):
		return self.fasta()

class DNA(BioSequence):
	"""Class for DNA sequences. IUPAC alphabet. Uppercase only."""

	canonical = ['A', 'C', 'G', 'T']
	extended = ['A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'W', 'S', 'B', 'D', 'H', 'V', 'N']

	def __init__(self, name=None, seq=None, desc=None):
		"""
		Parameters
		----------
		+ name=		`str` ideally a unique identifier
		+ desc=		`str` free text description of sequence
		+ seq=		`str` nucleotide sequence
		
		Attributes
		----------
		+ name
		+ desc
		+ seq
		+ ftable	`FeatureTable` object
		"""

		self.name = name
		self.seq = seq
		self.desc = desc
		self.ftable = FeatureTable(dna=self)

	def check_alphabet(self):
		"""Checks if sequence is in nucleotide alphabet. Throws `SequenceError`."""

		for i in range(len(self.seq)):
			nt = self.seq[i:i+1]
			if nt not in self.extended:
				raise SequenceError('letter not in DNA alphabet: ' + nt)

	def revcomp(self):
		"""Reverse-complements sequence and feature table."""
		self.check_alphabet()
		self.seq = toolbox.revcomp_str(self.seq)
		self.ftable._revcomp()

	def translate(self):
		"""Returns translated sequence as `str` protein sequence with no name or desc."""
		self.check_alphabet()
		pro = toolbox.translate_str(self.seq)
		return Protein(seq=pro)

class Protein(BioSequence):
	"""Class for protein sequences. Uppercase only. 20 aa + 'X' and '*'."""

	canonical = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
		'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	extended = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
		'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X', '*']

	def __init__(self, name=None, seq=None, desc=None, species=None):
		"""
		Parameters & Attributes
		-----------------------
		+ name= `str` ideally a unique identifier
		+ desc= `str` free text description
		+ seq=  `str` amino acid sequence
		"""

		self.name = name
		self.seq = seq
		self.desc = desc

	def check_alphabet(self):
		"""Check if sequence is in amino acid alphabet. Throws `SequenceError`."""
		for i in range(len(self.seq)):
			aa = self.seq[i:i+1]
			if aa not in self.extended:
				raise SequenceError('letter not in protein alphabet: ' + aa)

