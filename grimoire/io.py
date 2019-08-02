import re
import gzip

class IOError(Exception):
	pass

class GFF_entry:
	"""Class representing a GFF entry (row)"""

	def __init__(self, column):
		"""
		Parameters
		----------
		column: list
			A list of columns in GFF file
		"""

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
		"""
		Parameters
		----------
		filename: str
			Path to the GFF file
		"""

		self._chroms = {}
		self._types = {}
		gz = False
		fp = None
		if re.search('\.gz$', filename):
			fp = gzip.open(filename)
			gz = True
		else:
			fp = open(filename, 'r')

		while (1):
			line = fp.readline()
			if gz: line = str(line, 'utf-8')
			if line == '': break
			if line[0:1] == '#': continue
			col = line.split('\t')
			if len(col) < 8: continue
			chrom = col[0]
			type = col[2]
			entry = GFF_entry(col)
			if chrom not in self._chroms: self._chroms[chrom] = []
			self._chroms[chrom].append(entry)
			if type not in self._types: self._types[type] = []
			self._types[type].append(entry)
			self.chroms = list(self._chroms.keys())
			self.types = list(self._types.keys())
		fp.close()

	def get(self, type=None, chrom=None, beg=None, end=None):
		"""
		Retrieves GFF entries with given parameters. If parameter isn't
		specified, all are returned

		Parameters
		----------
		type: str
			Type of GFF entry (e.g. exon, gene) (default is None)
		chrom: str
			Chromosome of interest (default is None)
		beg: int
			Beginning coordinate (default is None)
		end: int
			Ending coordinate (default is None)
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
	"""Class for iterating through GFF records"""

	def __init__(self, filename=None, filepointer=None):
		"""
		Use either a path to file or filepointer object.

		Parameters
		----------
		filename: str
			Path to GFF file

		filepointer:
			Filepointer
		"""

		self._fp = None
		self._gz = False

		if filename != None:
			if re.search('\.gz$', filename):
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
		Returns the next row in the GFF file/stream.
		"""

		line = self._fp.readline()
		if self._gz: line = str(line, 'utf-8')
		if line == '':
			self._fp.close()
			raise StopIteration()
		if line[0:1] == '#': return self.next()
		col = line.split('\t')
		if len(col) < 8: return self.next()
		chrom = col[0]
		type = col[2]
		return GFF_entry(col)

class FASTA_entry:
	"""Class representing a FASTA entry"""

	def __init__(self, id, desc, seq):
		"""
		Parameters
		----------
		id: str
			Identifier of current FASTA entry
		desc: str
			Description/Info of entry (often empty)
		seq: str
			Sequence of entry
		"""

		self.id = id
		self.desc = desc
		self.seq = seq

class FASTA_file:
	"""Class for reading a FASTA file with random access"""

	def __init__(self, filename):
		"""
		Reads in FASTA file. Raises errors on duplicate ids.

		Parameters
		----------
		filename: str
			Path to FASTA file
		"""

		self.filename = filename
		if re.search('\.gz$', filename):
			raise IOError('.gz files not supported in FASTA_file')
		self._offset = {} # indexes identifiers to file offsets
		self.ids = []
		self._fp = open(self.filename, 'r')
		while (True):
			line = self.file.readline()
			if line == '': break
			if line[0:1] == '>':
				m = re.search('>\s*(\S+)', line)
				if m[1] in self.offset:
					raise IOError('duplicate id: ' + m[1])
				self.ids.append(m[1])
				self._offset[m[1]] = self.fp.tell() - len(line)
		self.file.close()

	def get(self, id):
		"""
		Retrieves FASTA entry with given identifier.

		Parameters
		----------
		id: str
			Identifier name
		"""

		self._fp = open(self.filename, 'r')
		self._fp.seek(self._offset[id])
		header = self._fp.readline()
		m = re.search('>\s*(\S+)\s*(.*)', header)
		id = m[1]
		desc = m[2]
		seq = []
		while (True):
			line = self._fp.readline()
			if line[0:1] == '>': break
			if line == '': break
			line = line.replace(' ', '')
			seq.append(line.strip())
		self.file.close()
		return FASTA_entry(id, desc, "".join(seq))

class FASTA_stream:
	"""Class for iterating through a FASTA file"""

	def __init__(self, filename=None, filepointer=None):
		"""
		Use path to file or file pointer.

		Parameters
		----------
		filename: str
			Path to file
		filepointer:
			File pointer object
		"""

		self._fp = None
		self._gz = False

		if filename != None:
			if re.search('\.gz$', filename):
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
				self.done = True
				self._fp.close()
				break

			line = line.replace(' ', '')
			seq.append(line.strip())

		return FASTA_entry(id, desc, "".join(seq))

