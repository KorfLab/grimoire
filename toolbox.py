"""toolbox of miscellaneous sequence definitions, functions, and classes
"""

import os
import sys
import re
from operator import itemgetter

ALPHABET = {
	'nt' : ['A', 'C', 'G', 'T'],
	'aa' : ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
		'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'],
}

CODE = {
	'standard': {		
		'AAA' : 'K',	'AAC' : 'N',	'AAG' : 'K',	'AAT' : 'N',
		'AAR' : 'K',	'AAY' : 'N',	'ACA' : 'T',	'ACC' : 'T',
		'ACG' : 'T',	'ACT' : 'T',	'ACR' : 'T',	'ACY' : 'T',
		'ACK' : 'T',	'ACM' : 'T',	'ACW' : 'T',	'ACS' : 'T',
		'ACB' : 'T',	'ACD' : 'T',	'ACH' : 'T',	'ACV' : 'T',
		'ACN' : 'T',	'AGA' : 'R',	'AGC' : 'S',	'AGG' : 'R',
		'AGT' : 'S',	'AGR' : 'R',	'AGY' : 'S',	'ATA' : 'I',	
		'ATC' : 'I',	'ATG' : 'M',	'ATT' : 'I',	'ATY' : 'I',
		'ATM' : 'I',	'ATW' : 'I',	'ATH' : 'I',	'CAA' : 'Q',
		'CAC' : 'H',	'CAG' : 'Q',	'CAT' : 'H',	'CAR' : 'Q',
		'CAY' : 'H',	'CCA' : 'P',	'CCC' : 'P',	'CCG' : 'P',
		'CCT' : 'P',	'CCR' : 'P',	'CCY' : 'P',	'CCK' : 'P',
		'CCM' : 'P',	'CCW' : 'P',	'CCS' : 'P',	'CCB' : 'P',
		'CCD' : 'P',	'CCH' : 'P',	'CCV' : 'P',	'CCN' : 'P',
		'CGA' : 'R',	'CGC' : 'R',	'CGG' : 'R',	'CGT' : 'R',
		'CGR' : 'R',	'CGY' : 'R',	'CGK' : 'R',	'CGM' : 'R',
		'CGW' : 'R',	'CGS' : 'R',	'CGB' : 'R',	'CGD' : 'R',
		'CGH' : 'R',	'CGV' : 'R',	'CGN' : 'R',	'CTA' : 'L',
		'CTC' : 'L',	'CTG' : 'L',	'CTT' : 'L',	'CTR' : 'L',
		'CTY' : 'L',	'CTK' : 'L',	'CTM' : 'L',	'CTW' : 'L',
		'CTS' : 'L',	'CTB' : 'L',	'CTD' : 'L',	'CTH' : 'L',
		'CTV' : 'L',	'CTN' : 'L',	'GAA' : 'E',	'GAC' : 'D',
		'GAG' : 'E',	'GAT' : 'D',	'GAR' : 'E',	'GAY' : 'D',
		'GCA' : 'A',	'GCC' : 'A',	'GCG' : 'A',	'GCT' : 'A',
		'GCR' : 'A',	'GCY' : 'A',	'GCK' : 'A',	'GCM' : 'A',
		'GCW' : 'A',	'GCS' : 'A',	'GCB' : 'A',	'GCD' : 'A',
		'GCH' : 'A',	'GCV' : 'A',	'GCN' : 'A',	'GGA' : 'G',
		'GGC' : 'G',	'GGG' : 'G',	'GGT' : 'G',	'GGR' : 'G',
		'GGY' : 'G',	'GGK' : 'G',	'GGM' : 'G',	'GGW' : 'G',
		'GGS' : 'G',	'GGB' : 'G',	'GGD' : 'G',	'GGH' : 'G',
		'GGV' : 'G',	'GGN' : 'G',	'GTA' : 'V',	'GTC' : 'V',
		'GTG' : 'V',	'GTT' : 'V',	'GTR' : 'V',	'GTY' : 'V',
		'GTK' : 'V',	'GTM' : 'V',	'GTW' : 'V',	'GTS' : 'V',
		'GTB' : 'V',	'GTD' : 'V',	'GTH' : 'V',	'GTV' : 'V',
		'GTN' : 'V',	'TAA' : '*',	'TAC' : 'Y',	'TAG' : '*',
		'TAT' : 'Y',	'TAR' : '*',	'TAY' : 'Y',	'TCA' : 'S',
		'TCC' : 'S',	'TCG' : 'S',	'TCT' : 'S',	'TCR' : 'S',
		'TCY' : 'S',	'TCK' : 'S',	'TCM' : 'S',	'TCW' : 'S',
		'TCS' : 'S',	'TCB' : 'S',	'TCD' : 'S',	'TCH' : 'S',
		'TCV' : 'S',	'TCN' : 'S',	'TGA' : '*',	'TGC' : 'C',
		'TGG' : 'W',	'TGT' : 'C',	'TGY' : 'C',	'TTA' : 'L',
		'TTC' : 'F',	'TTG' : 'L',	'TTT' : 'F',	'TTR' : 'L',
		'TTY' : 'F',	'TRA' : '*',	'YTA' : 'L',	'YTG' : 'L',
		'YTR' : 'L',	'MGA' : 'R',	'MGG' : 'R',	'MGR' : 'R',
	}
}

def revcomp(seq):
	"""Returns a reverse-complement of the input sequence"""
	comp = str.maketrans('ACGTRYMKWSBDHV', 'TGCAYRKMWSVHDB')
	seq = seq.translate(comp)[::-1]
	return seq

def translate(seq, table='standard'):
	"""Returns an amino acid translation of the input DNA sequence"""
	pro = []
	for i in range(0, len(seq), 3):
		codon = seq[i:i+3]
		if codon in CODE[table]: pro.append(CODE[table][codon])
		else: pro.append('X')
	return "".join(pro)

def _kmers(alphabet, table, key, n, k, v):
	if (k == 0) :
		if key not in table:
			table[key] = v
			return

	for i in range(n):
		t = key + alphabet[i]
		_kmers(alphabet, table, t, n, k - 1, v)

def generate_kmers(alphabet='nt', k=1, pseudo=0):
	"""Creates a dictionary of all kmers of either nt or aa alphabet."""
	table = {}
	if (alphabet == 'nt') :
		_kmers(ALPHABET[alphabet], table, '', 4, k, pseudo)
	elif (alphabet == 'aa') :
		_kmers(ALPHABET[alphabet], table, '', 20, k, pseudo)
	return table


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
				sys.stderr.write('error: type not defined: ' + type + '\n')
				sys.exit(1)
			type_search[type] = True
		else:
			for t in self.types: type_search[t] = True

		chrom_search = []
		if chrom:
			if chrom not in self._chroms:
				sys.stderr.write('error: chrom not defined: ' + chrom + '\n')
				sys.exit(1)
			chrom_search.append(chrom)
		else:
			chrom_search = self.chroms
		
		beg = 0 if not beg else beg
		end = 1e300 if not end else end
		if beg > end:
			sys.stderr.write('error: beg > end: ' + beg + ' ' + end + '\n')
			sys.exit(1)

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
					sys.stderr.write('error: duplicate id: ' + m[1] + '\n')
					sys.exit(1)
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
		else: sys.exit(1)
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

		return Entry(id, desc, "".join(seq))



