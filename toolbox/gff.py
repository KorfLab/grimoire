"""Class for reading GFF files"""

import os
import sys
import re
from operator import itemgetter

class Gff:
	
	def __init__(self, filename):
		self._chroms = {} 
		self._types = {}
		self.file = open(filename, 'r')
		gff = self.file.readlines()
		for line in gff:	
			m = re.search('(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)', line)
			chrom = m[1]
			type = m[3]
			entry = {
				'chrom' : chrom,
				'type' : type,
				'beg' : m[4],
				'end' : m[5],
				'strand' : m[7],
				'phase' : m[8],
				'attr' : m[9]
			}
			if (chrom not in self._chroms): self._chroms[chrom] = []
			self._chroms[chrom].append(entry)
			if (type not in self._types): self._types[type] = []
			self._types[type].append(entry)
		self.chroms = list(self._chroms.keys())
		self.types = list(self._types.keys())
	
	def get(self, type=None, range=None):
		if (type and type not in self._types): return None
		if (range):
			m = re.search('(\w+):*(\d*)\D*(\d*)', range)
			if (not m): return None
			chrom = m[1]
			if (not chrom or chrom not in self._chroms): return None
			beg = m[2]
			end = m[3]
			range_list = []
			if (not beg or not end):
				range_list = self._chroms[chrom]
			else:
				chrom_sort = sorted(list(self._chroms[chrom]), key=itemgetter('beg'))
				for entry in chrom_sort:
					if (entry['beg'] >= beg and entry['beg'] < end):
						range_list.append(entry)
					elif (entry['end'] > beg and entry['end'] < end):
						range_list.append(entry)
			if (type):
				range_type = []
				for entry in range_list:
					if (entry['type'] == type): range_type.append(entry)
				return range_type
			else:
				return range_list
					
		if (type and not range):
			return self._types[type]
		
