"""
Classes for reading standard bioinformatics file formats.
"""

import re
import gzip
import math

class FASTA_error(Exception):
	pass

class FASTA_entry:
	"""Class representing a FASTA entry."""

	def __init__(self, id, desc, seq):
		"""
		Parameters & Attributes
		-----------------------
		+ id   `str` unique id
		+ desc `str` free text description, which may be None
		+ seq  `str` sequence of symbols
		"""

		self.id = id
		self.desc = desc
		self.seq = seq

	def string(self, wrap=80):
		"""
		Returns a string formatted as a fasta file.
		
		Parameters
		----------
		+ wrap=`int` number of characters per line
		"""
		
		s = '>'
		if self.id: s += self.id
		if self.desc: s += ' ' + self.desc
		s += '\n'
		for i in range(0, len(self.seq), wrap):
			s += self.seq[i:i+wrap] + '\n'
		return s

	def __str__(self):
		return self.string()

class FASTA_file:
	"""Class for reading a FASTA file with random access."""

	def __init__(self, filename):
		"""
		Parameters
		----------
		+ filename `str` path to file
		
		Attributes
		----------
		+ ids - list of identifiers
		"""

		self._filename = filename
		if re.search(r'\.gz$', filename):
			raise NotImplementedError('.gz files not supported in FASTA_file')
		self._offset = {} # indexes identifiers to file offsets
		self.ids = []
		self._fp = open(self._filename, 'r')
		while (True):
			line = self._fp.readline()
			if line == '': break
			if line[0:1] == '>':
				m = re.search(r'>\s*(\S+)', line)
				if m[1] in self._offset:
					raise FASTA_error('duplicate id: ' + m[1])
				self.ids.append(m[1])
				self._offset[m[1]] = self._fp.tell() - len(line)
		self._fp.close()

	def get(self, id):
		"""
		Returns a `FASTA_entry` given an id.
		
		Parameters
		----------
		+ id `str` identifier
		"""

		self._fp = open(self._filename, 'r')
		self._fp.seek(self._offset[id])
		header = self._fp.readline()
		m = re.search(r'>\s*(\S+)\s*(.*)', header)
		id = m[1]
		desc = m[2]
		seq = []
		while (True):
			line = self._fp.readline()
			if line[0:1] == '>': break
			if line == '': break
			line = line.replace(' ', '')
			seq.append(line.strip())
		self._fp.close()
		return FASTA_entry(id, desc, "".join(seq))

class FASTA_stream:
	"""Class for iterating through a FASTA file."""

	def __init__(self, filename=None, filepointer=None):
		"""
		Parameters
		----------
		+ filename    `str` name of a file
		+ filepointer `obj` bytes-like object
		
		Specify filename or filepointer, not both
		"""

		self._fp = None
		self._gz = False

		if filename != None:
			if re.search(r'\.gz$', filename):
				self._fp = gzip.open(filename)
				self._gz = True
			else:
				self._fp = open(filename, 'r')
		elif filepointer != None:
			self._fp = filepointer
		else:
			raise IOError('no file or filepointer given')
		self._lastline = ''
		self.done = False

	def __iter__(self):
		return self

	def __next__(self):
		return self.next()

	def next(self):
		"""
		Retrieves the next entry of the FASTA file/stream.
		"""

		if self.done: raise StopIteration()
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
				self.done = True
				self._fp.close()
				break

			line = line.replace(' ', '')
			seq.append(line.strip())

		return FASTA_entry(id, desc, "".join(seq))

class GFF_skip(Exception):
	pass

class GFF_error(Exception):
	pass

class GFF_entry:
	"""Represents a GFF entry (row)."""
	
	def __init__(self, line):
		"""
		Parameters
		----------
		+ line `str` a line from a GFF file
		
		Attributes
		----------
		+ chrom  `str`   chromosome
		+ source `str`   entity that created this feature (arbitrary)
		+ type   `str`   string describing feature type (hopefully from SO)
		+ beg    `int`   1-based coordinate
		+ end    `int`   1-based coordinate, always bigger than beg
		+ score  `float` number or '.' if not defined
		+ strand `str`   either '+' or '-', or '.' if not defined
		+ phase  `int`   0, 1, or 2, or '.' if not defined
		+ attr   `str`  structured string for extra information (e.g. grouping)
		"""

		if line[0:1] == '#':
			raise GFF_skip
		if line[0:1] == '\n':
			raise GFF_skip
			
		column = line.split('\t')
		if len(column) < 8:
			raise GFF_error('badly formatted gff')
			
		self.chrom = column[0]
		self.source = column[1]
		self.type = column[2]
		self.beg = int(column[3])
		self.end = int(column[4])
		self.score = column[5]
		self.strand = column[6]
		self.phase = column[7]
		if len(column) == 9:
			self.attr = column[8]

class GFF_file:
	"""Reading and searching GFF files (slurps all into memory)."""
	
	def __init__(self, filename):
		"""
		Parameters
		----------
		+ filename `str` path to GFF file, which may be gzipped
		
		Attributes
		----------
		+ chroms - list of chromosome names
		+ types - list of feature type names
		"""

		self._chroms = {}
		self._types = {}
		gz = False
		fp = None
		if re.search(r'\.gz$', filename):
			fp = gzip.open(filename)
			gz = True
		else:
			fp = open(filename, 'r')

		while (1):
			line = fp.readline()
			if gz: line = str(line, 'utf-8')
			if line == '': break
			
			gff = None
			try:
				gff = GFF_entry(line)
			except GFF_skip:
				continue
			except GFF_error:
				raise GFF_error('badly formatted gff')
			
			if gff.chrom not in self._chroms: self._chroms[gff.chrom] = []
			self._chroms[gff.chrom].append(gff)
			if gff.type not in self._types: self._types[gff.type] = []
			self._types[gff.type].append(gff)

		fp.close()
		self.chroms = list(self._chroms.keys())
		self.types = list(self._types.keys())

	def get(self, type=None, chrom=None, beg=None, end=None):
		"""
		Searches for GFF entries and returns them in a list of `GFF_entry`.

		Parameters
		----------
		+ type=  `str`  type of GFF entry (e.g. exon), all if not specified
		+ chrom= `str` chromosome (e.g. I), all if not specified
		+ beg=   `int`   1-based begin coordinate, all if not specified
		+ end=   `int`   1-based end coordinate, all if not specified
		"""

		type_search = {}
		if type:
			if type not in self._types:
				raise IOErrorboxError('type not defined: ' + type)
			type_search[type] = True
		else:
			for t in self._types: type_search[t] = True

		chrom_search = []
		if chrom:
			if chrom in self._chroms:
				chrom_search.append(chrom)
		else:
			chrom_search = self.chroms

		beg = -1 if not beg else beg
		end = math.inf if not end else end
		if beg > end:
			raise IOError('beg > end: ' + beg + '-' + end)

		found = []
		for c in chrom_search:
			for entry in self._chroms[c]:
				if entry.type in type_search:
					if entry.beg >= beg and entry.end <= end:
						found.append(entry)

		return found

class GFF_stream:
	"""Class for iterating through GFF records."""

	def __init__(self, filename=None, filepointer=None):
		"""
		Parameters
		----------
		+ filename    `str` name of a file
		+ filepointer `obj` bytes-like object
		
		Specify filename or filepointer, not both
		"""
		
		self._fp = None
		self._gz = False

		if filename != None:
			if re.search(r'\.gz$', filename):
				self._fp = gzip.open(filename)
				self._gz = True
			else:
				self._fp = open(filename, 'r')
		elif filepointer != None:
			self._fp = filepointer
		else:
			raise IOError('no file or filepointer given')

	def __iter__(self):
		return self

	def __next__(self):
		return self.next()

	def next(self):
		"""
		Returns the next `GFF+entry` in the GFF file/stream.
		"""

		line = self._fp.readline()
		if self._gz: line = str(line, 'utf-8')
		if line == '':
			self._fp.close()
			raise StopIteration()
		
		gff = None
		try:
			gff = GFF_entry(line)
		except GFF_skip:
			return self.next()
		except GFF_error:
			raise GFF_error('badly formatted gff')
		
		return gff

