"""Class for reading GFF files"""

import os
import sys
import re
from operator import itemgetter

class Entry:

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

class Gff:
	
	def __init__(self, filename):
		self._chroms = {} 
		self._types = {}
		self.file = open(filename, 'r')
		gff = self.file.readlines()
		for line in gff:
			col = line.split('\t')
			chrom = col[0]
			type = col[2]
			entry = Entry(col)
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

