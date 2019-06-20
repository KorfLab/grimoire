"""Class for reading FASTA files"""

import os
import sys
import re

class Entry:

	def __init__(self, id, desc, seq):
		self.id = id
		self.desc = desc
		self.seq = seq

class FastaFile:
	
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
		return Entry(id, desc, "".join(seq))

class FastaStream:
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



