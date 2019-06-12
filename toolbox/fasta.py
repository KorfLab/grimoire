"""Class for reading FASTA files"""

import os
import sys
import re

class Entry:

	def __init__(self, id, desc, seq):
		self.id = id
		self.desc = desc
		self.seq = seq

class Fasta:
	
	def __init__(self, filename):
		self.offset = {} # indexes identifiers to file offsets
		self.ids = []
		self.file = open(filename, 'r')
		while (True):
			line = self.file.readline()
			if line == '': break
			if line[0:1] == '>':
				m = re.search('>\s*(\S+)', line)
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
		return Entry(id, desc, "".join(seq))

