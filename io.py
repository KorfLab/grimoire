"""Module for IO of various standard formats"""

import re

class IOError(Exception):
	pass

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
	"""Class representing a GFF file"""
	
	def __init__(self, filename):
		self._chroms = {} 
		self._types = {}
		self.file = open(filename, 'r')
		while (1):
			line = self.file.readline()
			if line == '': break
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
				raise IOError('type not defined: ' + type)
			type_search[type] = True
		else:
			for t in self.types: type_search[t] = True

		chrom_search = []
		if chrom:
			if chrom not in self._chroms:
				raise IOError('chrom not defined: ' + chrom)
			chrom_search.append(chrom)
		else:
			chrom_search = self.chroms
		
		beg = 0 if not beg else beg
		end = 1e300 if not end else end
		if beg > end:
			raise IOError('beg > end: ' + beg + '-' + end)


		found = []
		for c in chrom_search:
			for entry in self._chroms[c]:
				if entry.type in type_search:
					if entry.beg >= beg and entry.end <= end:
						found.append(entry)
		
		return found


class FASTA_entry:
	"""Class representing a FASTA entry (header, seq)"""

	def __init__(self, id, desc, seq):
		self.id = id
		self.desc = desc
		self.seq = seq

class FASTA_file:
	"""Class representing a FASTA file"""
	
	def __init__(self, filename):
		self.offset = {} # indexes identifiers to file offsets
		self.ids = []
		self.file = open(filename, 'r')
		while (True):
			line = self.file.readline()
			if line == '': break
			if line[0:1] == '>':
				m = re.search('>\s*(\S+)', line)
				if m[1] in self.offset:
					raise IOError('duplicate id: ' + m[1])
				self.ids.append(m[1])
				self.offset[m[1]] = self.file.tell() - len(line)
	
	def get(self, id):
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
		return FASTA_entry(id, desc, "".join(seq))

class FASTA_stream:
	def __init__(self, filename=None, filepointer=None):
		self.fp = None
		if   filename    != None: self.fp = open(filename, 'r')
		elif filepointer != None: self.fp = filepointer
		else: raise IOError('no file or filepointer given')
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
				break

			line = line.replace(' ', '')
			seq.append(line.strip())

		return FASTA_entry(id, desc, "".join(seq))



"""
sys.exit(1)


## convert bed12 to GFF3

fp = open('araport.bed', 'r')
source = 'ARAPORT11'
# genes can't be created until all transcripts are created
genes = {}
while (1):
	line = fp.readline()
	if line == '': break
	col = line.split('\t')
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
	
	if gid not in genes:
		print('\t'.join([chr_id, source, 'gene', str(chr_beg), str(chr_end),
			score, strand, 'ID='+gid]))
		genes[gid] = True
	
	print('\t'.join([chr_id, source, 'mRNA', str(chr_beg), str(chr_end),
		score, strand, attr]))

	for i in range(n):
		beg = chr_beg + int(starts[i])
		end = beg + int(sizes[i]) -1
		attr = 'Parent=' + txid
		print('\t'.join([chr_id, source, 'exon', str(beg), str(end),
			score, strand, attr]))

sys.exit(1)


## build from GFF3
gen = genome.Genome(fasta='data/TAIR10_1.fasta', gff3='data/TAIR10_1.gff3')
for chr in gen.chromosomes:
	for gene in chr.genes:
		if gene.issues:
#			print(gene.id, 'has issues')
			for tx in gene.transcripts:
				if tx.issues:
#					print('\t', tx.id, 'has issues')
					for issue in tx.issues:
#						print('\t\t', issue)
						pass

"""

