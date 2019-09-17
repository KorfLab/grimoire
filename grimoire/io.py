"""
Classes for reading standard bioinformatics file formats.
"""

import re
import gzip
import math
import tempfile

class FASTA_error(Exception):
	pass

class FASTA_entry:
	"""Class representing a FASTA entry."""

	def __init__(self, id, desc, seq):
		"""
		Parameters & Attributes
		-----------------------
		+ id   `str` unique id
		+ desc `str` free text description, which may be `None`
		+ seq  `str` sequence of symbols
		"""

		self.id = id
		self.desc = desc
		self.seq = seq
		if not self.seq or self.seq == '\n':
			raise FASTA_error('empty string is not a valid sequence')

	def string(self, wrap=80):
		"""
		Returns a string formatted as a FASTA file.
		
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
		+ ids	`list` identifiers of entries in file
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
		if not self.ids: raise FASTA_error('file is not properly formatted FASTA')

	def get(self, id):
		"""
		Returns a `FASTA_entry` corresponding to given identifier.
		
		Parameters
		----------
		+ id `str` unique identifier
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
	"""Class for reading a FASTA file as a stream. Iterable allowing access to
	`FASTA_entry` objects representing contents of FASTA file."""

	def __init__(self, filename=None, filepointer=None):
		"""
		Parameters
		----------
		+ filename=    `str` name of a file
		+ filepointer= `obj` bytes-like object
		
		Specify filename or filepointer, not both.
		"""

		self._fp = None

		if filename != None:
			if re.search(r'\.gz$', filename):
				self._fp = gzip.open(filename)
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
		if self.done: raise StopIteration()
		header = None
		if self._lastline[0:1] == '>':
			header = self._lastline
		else:
			header = self._fp.readline()
			if isinstance(header, bytes): header = header.decode()

		m = re.search(r'>\s*(\S+)\s*(.*)', header)
		if m == None: raise FASTA_error('file is not properly formatted FASTA')
		id = m[1]
		desc = m[2]
		seq = []

		while (True):
			line = self._fp.readline()
			if isinstance(line, bytes): line = line.decode()
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
	"""Class representing a GFF entry (row)."""
	
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
		+ attr   `str`   structured string for extra information (e.g. grouping)
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
		else:
			self.attr = ''

def _from_GTF(file):
	
	fp = None
	if re.search(r'\.gz$', file): fp = gzip.open(file)
	else:                         fp = open(file, 'r')
		
	# stuff
	temp = tempfile.TemporaryFile()
	gffs = []
	stops = []

	# parse the file
	for line in fp:
		if isinstance(line, bytes): line = line.decode()
		if line[0:1]=='#': continue
		col = line.rstrip().split('\t')
		if len(col) < 9: continue
			
		if col[2] == 'gene':
			id = re.search('gene_id "(\S+)"', col[8])[1]
			col[8] = 'ID=' + id
			gffs.append(col)
		elif col[2] == 'transcript':
			pid = re.search('gene_id "(\S+)"', col[8])[1]
			tid = re.search('transcript_id "(\S+)"', col[8])[1]
			col[2] = 'mRNA'
			col[8] = 'ID=' + tid + ';Parent=' + pid
			gffs.append(col)
		elif col[2] == 'exon' or col[2] == 'CDS':
			pid = re.search('transcript_id "(\S+)"', col[8])[1]
			col[8] = 'Parent=' + pid
			gffs.append(col)
		elif col[2] == 'stop_codon':
			stops.append(col[3])
		
	for gff in gffs:
		if gff[2] == 'CDS':
			if gff[6] == '+':
				end = int(gff[4])
				if str(end + 1) in stops:
					gff[4] = str(end + 3)
			elif gff[6] == '-':
				beg = int(gff[3])
				if str(beg -3) in stops:
					gff[3] = str(beg -3)
	
	# return
	for gff in gffs: temp.write(str.encode('\t'.join(gff) + '\n'))
	temp.seek(0)
	return temp

def _from_BED12(file):

	fp = None
	if re.search(r'\.gz$', file): fp = gzip.open(file)
	else:                         fp = open(file, 'r')
		
	# stuff
	temp = tempfile.TemporaryFile()
	genes = {}
	
	# parse the file
	for line in fp:
		if isinstance(line, bytes): line = line.decode()
		col = line.split('\t')
		if len(col) < 12: continue
		
		chr_id = col[0]
		chr_beg = int(col[1]) + 1
		chr_end = int(col[2])
		txid = col[3]
		score = col[4]
		strand = col[5]
		cds_beg = int(col[6]) + 1
		cds_end = int(col[7])
		rgb = col[8]
		n = int(col[9])
		sizes = col[10].split(',')
		starts = col[11].split(',')
		gid = re.search('(\w+)\.\d+', txid)[1]
		attr = 'ID=' + txid + ';Parent=' + gid
		
		# create the gene feature if necessary
		if gid not in genes:
			temp.write(str.encode('\t'.join([chr_id, '.', 'gene', str(chr_beg),
				str(chr_end), score, strand, '.', 'ID='+gid]) + '\n'))
			genes[gid] = True

		# create the mRNA feature
		temp.write(str.encode('\t'.join([chr_id, '.', 'mRNA', str(chr_beg),
			str(chr_end), score, strand, attr]) + '\n'))

		# create the exons and CDS features
		for i in range(n):
			beg = chr_beg + int(starts[i])
			end = beg + int(sizes[i]) -1
			attr = 'Parent=' + txid
			temp.write(str.encode('\t'.join([chr_id, '.', 'exon', str(beg),
				str(end), score, strand, ',', attr]) + '\n'))
			if beg <= cds_end and end >= cds_beg:
				cb, ce = None, None
				if beg > cds_beg: cb = beg
				else:             cb = cds_beg
				if end > cds_end: ce = cds_end
				else:             ce = end
				temp.write(str.encode('\t'.join([chr_id, 'araport', 'CDS',
					str(cb), str(ce), score, strand, '.', attr]) + '\n'))
	
	# return the bytes object
	temp.seek(0)
	return temp

class GFF_file:
	"""Class for reading and searching a GFF (and other) annotation files."""
	
	def __init__(self, file=None):
		"""
		Parameters
		----------
		+ file=   `str` path to file, which may be compressed
		
		Attributes
		----------
		+ chroms	`list` chromosome names
		+ types		`list` feature type names

		File extensions
		---------------
		GFF_file reads GFF3 files natively. It also autoconverts files in GTF
		and BED12 into GFF3. All conversions are by file extension mappings.
		
		+ GFF3: *.gff3, *.gff3.gz
		+ GTF: *.gtf, *.gtf.gz
		+ BED12: *.bed, *.bed.gz
		"""
		
		# if it's not GFF3, turn it into GFF3 before proceeding
		fp = None
		if   re.search(r'\.gff3\.gz$', file): fp = gzip.open(file)
		elif re.search(r'\.gff3$',     file): fp = open(file, 'r')
		elif re.search(r'\.gff\.gz$',  file): fp = gzip.open(file)
		elif re.search(r'\.gff$',      file): fp = open(file, 'r')
		elif re.search(r'\.gtf\.gz$',  file): fp = _from_GTF(file)
		elif re.search(r'\.gtf$',      file): fp = _from_GTF(file)
		elif re.search(r'\.bed\.gz$',  file): fp = _from_BED12(file)
		elif re.search(r'\.bed$',      file): fp = _from_BED12(file)
		else: raise IOError('unknown format: ' + file)
				
		self._chroms = {}
		self._types = {}

		while (1):
			line = fp.readline()
			if isinstance(line, bytes): line = line.decode()
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
		Returns list of `GFF_entry` objects matching search criteria.

		Parameters
		----------
		+ type=  `str`  type of GFF entry (e.g. exon), all if not specified
		+ chrom= `str` 	chromosome (e.g. I), all if not specified
		+ beg=   `int`  1-based begin coordinate, all if not specified
		+ end=   `int`  1-based end coordinate, all if not specified
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
	"""Class for reading a GFF file as a stream. Iterable allowing access to
	`GFF_entry` objects representing contents of GFF file."""

	def __init__(self, filename=None, filepointer=None):
		"""
		Parameters
		----------
		+ filename=    `str` name of a file
		+ filepointer= `obj` bytes-like object
		
		Specify filename or filepointer, not both.
		"""
		
		self._fp = None

		if filename != None:
			if   re.search(r'\.gff3\.gz$', filename): filepointer = gzip.open(filename)
			elif re.search(r'\.gff3$',     filename): filepointer = open(filename, 'r')
			elif re.search(r'\.gff\.gz$',  filename):
				raise NotImplementedError('GFF_file accepts GFF3 formatted files only')
			elif re.search(r'\.gff$',      filename):
				raise NotImplementedError('GFF_file accepts GFF3 formatted files only')
			elif re.search(r'\.gtf\.gz$',  filename):
				raise NotImplementedError('GFF_file does not accept GTF files')
			elif re.search(r'\.gtf$',      filename):
				raise NotImplementedError('GFF_file does not accept GTF files')
			elif re.search(r'\.bed\.gz$',  filename):
				raise NotImplementedError('GFF_file does not accept BED files')
			elif re.search(r'\.bed$',      filename):
				raise NotImplementedError('GFF_file does not accept BED files')
			else:
				raise NotImplementedError('Unrecognized file type. Please use GFF3.')
			self._fp = filepointer
		else:
			raise IOError('no file or filepointer given')

	def __iter__(self):
		return self

	def __next__(self):
		line = self._fp.readline()
		if isinstance(line, bytes): line= line.decode()
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


