"""Mostly IO functions for standard interchange formats"""

import math
import re
from functools import reduce


class ToolboxError(Exception):
	pass

def log(p):
	if p < 0: raise ValueError('p < 0')
	if p == 0: return -999
	else:      return math.log(p)

def sumlog(v1, v2):
	if v1 < v2: v1, v2 = v2, v1
	return math.log(1 + math.exp(v2 - v1)) + v1

def prod(iterable):
	return reduce(operator.mul, iterable, 1)

class GFF_entry:
	"""Class representing a GFF entry (row)"""

	def __init__(self, column):
		self.chrom = column[0]
		self.source = column[1]
		self.type = column[2]
		self.beg = int(column[3])
		self.end = int(column[4])
		self.score = column[5]
		self.strand = column[6]
		self.phase = column[7]
		self.attr = column[8]

class GFF_file:
	"""Class for reading and searching GFF files (slurps all into memory)."""

	def __init__(self, filename):
		self._chroms = {}
		self._types = {}
		with open(filename, 'r') as self.file:
			while (1):
				line = self.file.readline()
				if line == '': break
				if line[0:1] == '#': continue
				col = line.split('\t')
				chrom = col[0]
				type = col[2]
				entry = GFF_entry(col)
				if chrom not in self._chroms: self._chroms[chrom] = []
				self._chroms[chrom].append(entry)
				if type not in self._types: self._types[type] = []
				self._types[type].append(entry)
				self.chroms = list(self._chroms.keys())
				self.types = list(self._types.keys())

	def get(self, type=None, chrom=None, beg=None, end=None):
		type_search = {}
		if type:
			if type not in self._types:
				raise ToolboxError('type not defined: ' + type)
			type_search[type] = True
		else:
			for t in self.types: type_search[t] = True

		chrom_search = []
		if chrom:
			if chrom in self._chroms:
				chrom_search.append(chrom)
		else:
			chrom_search = self.chroms

		beg = 0 if not beg else beg
		end = 1e300 if not end else end
		if beg > end:
			raise ToolboxError('beg > end: ' + beg + '-' + end)


		found = []
		for c in chrom_search:
			for entry in self._chroms[c]:
				if entry.type in type_search:
					if entry.beg >= beg and entry.end <= end:
						found.append(entry)

		return found

class GFF_stream:
	"""Class for reading GFF records one at a time"""

	def __init__(self, filename=None, filepointer=None):
		self.fp = None
		if   filename    != None: self.fp = open(filename, 'r')
		elif filepointer != None: self.fp = filepointer
		else: raise ToolboxError('no file or filepointer given')

	def __iter__(self):
		return self

	def __next__(self):
		return self.next()

	def next(self):
		line = self.fp.readline()
		if line == '':
			self.fp.close()
			raise StopIteration()
		if line[0:1] == '#': return self.next()
		col = line.split('\t')
		chrom = col[0]
		type = col[2]
		return GFF_entry(col)

class FASTA_entry:
	"""Class representing a FASTA entry (id, desc, seq)"""

	def __init__(self, id, desc, seq):
		self.id = id
		self.desc = desc
		self.seq = seq

class FASTA_file:
	"""Class for reading a FASTA file with random acess"""

	def __init__(self, filename):
		self.filename = filename
		self.offset = {} # indexes identifiers to file offsets
		self.ids = []
		self.file = open(self.filename, 'r')
		while (True):
			line = self.file.readline()
			if line == '': break
			if line[0:1] == '>':
				m = re.search('>\s*(\S+)', line)
				if m[1] in self.offset:
					raise ToolboxError('duplicate id: ' + m[1])
				self.ids.append(m[1])
				self.offset[m[1]] = self.file.tell() - len(line)
		self.file.close()

	def get(self, id):
		self.file = open(self.filename, 'r')
		self.file.seek(self.offset[id])
		header = self.file.readline()
		m = re.search('>\s*(\S+)\s*(.*)', header)
		id = m[1]
		desc = m[2]
		seq = []
		while (True):
			line = self.file.readline()
			if line[0:1] == '>': break
			if line == '': break
			line = line.replace(' ', '')
			seq.append(line.strip())
		self.file.close()
		return FASTA_entry(id, desc, "".join(seq))

class FASTA_stream:
	"""Class for reading FASTA records in a stream"""

	def __init__(self, filename=None, filepointer=None):
		self.fp = None
		if   filename    != None: self.fp = open(filename, 'r')
		elif filepointer != None: self.fp = filepointer
		else: raise ToolboxError('no file or filepointer given')
		self.lastline = ''
		self.done = False

	def __iter__(self):
		return self

	def __next__(self):
		return self.next()

	def next(self):
		if self.done: raise StopIteration()
		header = None
		if self.lastline[0:1] == '>': header = self.lastline
		else:                         header = self.fp.readline()

		m = re.search('>\s*(\S+)\s*(.*)', header)
		id = m[1]
		desc = m[2]
		seq = []

		while (True):
			line = self.fp.readline()
			if line[0:1] == '>':
				self.lastline = line
				break
			if line == '':
				self.done = True
				self.fp.close()
				break

			line = line.replace(' ', '')
			seq.append(line.strip())

		return FASTA_entry(id, desc, "".join(seq))
