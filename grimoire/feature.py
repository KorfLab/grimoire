"""
Classes for sequence features
"""

import operator

import grimoire.toolbox as toolbox

class FeatureError(Exception):
	pass

class Feature:
	"""Class representing a sequence feature, which may have children"""

	def __init__(self, dna, beg, end, strand, type,
			phase='.', score='.', source='.', id=None, pid=None):
		"""
		Parameters
		----------
		+ dna     `DNA`   object of type `sequence.DNA`
		+ beg     `int`   1-based coordinate
		+ end     `int`   1-based coordinate
		+ strand  `str`   '+' or '-', or '.' for undefined
		+ type    `str`   token, hopefully SO-compliant
		+ phase=  `.`     may be {0, 1, 2} or '.' for undefined
		+ score=  `float` often unspecified as `.`
		+ source= `.`     creator of feature (e.g. WormBase)
		+ id=     `str`   unique identifier, optional
		+ pid=    `str`   a parent id or list of parent ids
		
		Attributes
		----------
		+ issues - a list of warnings and errors
		+ children - a list of sub-features
		+ validated - `bool` check if this has been validated
		"""

		self.dna = dna
		self.beg = beg
		self.end = end
		self.length = end - beg + 1
		self.strand = strand
		self.type = type
		self.phase = phase
		self.id = id
		self.pid = []
		if pid:
			if   isinstance(pid, list): self.pid = pid
			elif isinstance(pid, str):  self.pid.append(pid)
			else: raise FeatureError('pid must be string or list of strings')
		self.score = score
		self.source = source
		self.issues = {}
		self.children = []
		self.validated = False

		if self.dna == None:
			raise GenomeError('attempt to create feature unbound to DNA')
		self._validate()

	def _validate(self):
		if self.beg < 0: self.issues['beg<0'] = True
		if self.beg > self.end: self.issues['beg>end'] = True
		if self.end > len(self.dna.seq): self.issues['end>seq'] = True
		if self.children:
			for child in self.children:
				child.validate()
				if child.beg < self.beg:
					self.issues['child.beg<parent.beg'] = True
				if child.end > self.end:
					self.issues['child.end>parent.end'] = True
				if child.strand != self.strand:
					self.issues['mixed_strands'] = True
				if child.issues:
					self.issues['child_issues'] = True

	def validate(self):
		"""Check `Feature` instance for common errors."""

		if self.validated: return
		self._validate()
		self.validated = True

	def add_parent(self, pid):
		"""
		Add a parent identifier to a `Feature` object
		
		Parameters
		----------
		+ pid `str` parent identifier token (e.g. gene name)
		"""
		
		if id not in self.pid:
			self.pid.append(id)

	def add_child(self, child):
		"""
		Add a child to `Feature` object. Unsets validated flag.

		Parameters
		----------
		+ child `Feature` feature object
		"""

		self.validated = False
		if not self.id:
			raise GenomeError('parent feature requires ID')
		else:
			self.children.append(child)

	def seq_str(self):
		"""Returns the sequence of a `Feature` as a `str` in + strand."""

		seq = self.dna.seq[self.beg-1:self.end]
		if self.strand == '-': seq = toolbox.revcomp_str(seq)
		return seq

	def gff(self):
		"""Returns a string representation of a feature in GFF3 format."""
		if not self.validated: self.validate()

		attr = ''
		if self.id and self.pid:
			attr = 'ID=' + self.id + ';Parent=' + ','.join(self.pid)
		elif self.id:
			attr = 'ID=' + self.id
		elif self.pid:
			attr = 'Parent=' + ','.join(self.pid)

		string = '\t'.join([self.dna.name, self.source, self.type,
			str(self.beg), str(self.end), str(self.score),
			self.strand, self.phase, attr])
		if self.children:
			string += '\n'
			stuff = []
			for child in self.children:
				stuff.append(child.gff())
			string += '\n'.join(stuff) + '\n'
		return string

	def overlap(self, f2):
		"""
		Determines if two features overlap.

		Parameters
		----------
		+ f2 `Feature` the other feature (self implicitly f1)
		"""

		if self.dna.name == f2.dna.name:
			if self.beg >= f2.beg and self.beg <= f2.end: return True
			if self.end >= f2.beg and self.end <= f2.end: return True
			if self.beg <= f2.beg and self.end >= f2.end: return True
		return False
	
	def __str__(self):
		return self.gff()

class Transcript(Feature):
	"""Base class for transcripts (children of `Gene`s)."""

	def _check_overlaps(self, f, type):
		for i in range(1, len(f)):
			if f[i-1].end >= f[i].beg:
				self.issues['overlap_' + type] = True

	def _check_lengths(self, features, type):
		for f in features:
			if f.length < self.limit[type]['min']:
				self.issues['short_' + type] = True
			if f.length > self.limit[type]['max']:
				self.issues['long_' + type] = True

	def tx_str(self):
		""""Returns the transcript as a string with introns removed."""
		
		if not self.validated: self.validate()
		seq = []
		for exon in self.exons: seq.append(exon.seq_str())
		if self.strand == '-': seq.reverse()
		return ''.join(seq)

class mRNA(Transcript):
	"""Class for protein-coding mRNAs. Not for ncRNAs.
	
	Extends the base `Feature` class. mRNAs are created by first instantiating
	a parent `mRNA`. Then `add_child()` of type 'exon' _and_ 'CDS'.
	Introns and untranslated regions are inferred from the coordinates
	of the exon and CDS features. The mRNA must have an 'id' so that its
	children can reference it. You don't need to spcifiy 'pid' in the
	children, as this will be assigned automatically.
	
	Parameters
	----------
	See the `Feature` base class. Make sure you set id.
	
	Attributes (extended from the base class)
	----------
	+ exons - list of exons
	+ inrons - list of introns
	+ cdss - list of CDSs
	+ utr5s - list of 5' UTRs
	+ utr3s - list of 3' UTRs
	"""

	clade = 'std' #standard
	limit = {
		'exon':   {'min':20, 'max':10000},
		'cds':    {'min':1,  'max':10000},
		'utr5':   {'min':0,  'max':1000},
		'utr3':   {'min':0,  'max':1000},
		'intron': {'min':30, 'max':10000},
	}
	dons = {'GT':True}
	accs = {'AG':True}
	starts = {'ATG':True}
	stops = {'TAA':True, 'TGA':True, 'TAG':True}

	def set_rules(self, clade='std'):
		"""Set rules for mRNA by clade.
		
		+ dons - dictionary of donor sites (e.g. GT)
		+ accs - dictionary of acceptor sites (e.g. AG)
		+ starts - dictionary of start sites (e.g. ATG)
		+ stops - dictionary of stop sites (e.g. TAA, TGA, TAG)
		+ limit - 2d dictionary of length limits

		Parameters
		----------
		+ clade `str` name of clade ('std' or 'mammal')
		"""

		if clade == 'std':
			pass # the defaults are considered standard
		elif clade == 'mammal':
			self.limit['intron'][min] = 50
			self.limit['intron'][max] = 100000
		else:
			raise NotImplemented('clade not yet supported: ' + clade)
		self.clade = clade

	def validate(self):
		"""Checks mRNA for common errors, populating issues."""

		if self.validated: return
		if not self.id: raise FeatureError('mRNAs must have ids')
		self._validate()

		# mRNA properties (extended constructor)
		self.is_coding = True
		self.exons = []
		self.introns = []
		self.cdss = []
		self.utr5s = []
		self.utr3s = []

		for f in self.children:
			if   f.type == 'exon': self.exons.append(f)
			elif f.type == 'CDS': self.cdss.append(f)
			else: raise FeatureError('mRNA takes exon and CDS only')
			if f.pid == None:
				f.pid = self.id

		if len(self.cdss) == 0:
			self.issues['no_CDS'] = True

		self.exons.sort(key = operator.attrgetter('beg'))
		self.cdss.sort(key = operator.attrgetter('beg'))

		# create introns from exons
		for i in range(len(self.exons)-1):
			beg = self.exons[i].end +1
			end = self.exons[i+1].beg -1
			self.introns.append(
				Feature(self.dna, beg, end, self.strand, 'intron'))

		# create 5' and 3' UTRs from exons and CDSs
		if len(self.cdss) > 0:
			cds_beg = self.cdss[0].beg
			cds_end = self.cdss[-1].end
			for exon in self.exons:
				if exon.beg < cds_beg:
					ub = exon.beg
					ue = None
					if exon.end < cds_beg: ue = exon.end
					else:                  ue = cds_beg -1
					if exon.strand == '+':
						self.utr5s.append(Feature(self.dna, ub, ue, self.strand,
							'five_prime_UTR'))
					else:
						self.utr3s.append(Feature(self.dna, ub, ue, self.strand,
							'three_prime_UTR'))
				if exon.end > cds_end:
					ub = None
					ue = exon.end
					if exon.beg < cds_end: ub = cds_end + 1
					else:                  ub = exon.beg
					if exon.strand == '+':
						self.utr3s.append(Feature(self.dna, ub, ue, self.strand,
							'three_prime_UTR'))
					else:
						self.utr5s.append(Feature(self.dna, ub, ue, self.strand,
							'five_prime_UTR'))

		# check for overlapping features
		self._check_overlaps(self.exons, 'exon')
		self._check_overlaps(self.cdss, 'cds')
		self._check_overlaps(self.utr5s, 'utr5')
		self._check_overlaps(self.utr3s, 'utr3')
		self._check_overlaps(self.introns, 'intron')

		# check for unusual lengths
		self._check_lengths(self.exons, 'exon')
		self._check_lengths(self.cdss, 'cds')
		self._check_lengths(self.utr5s, 'utr5')
		self._check_lengths(self.utr3s, 'utr3')
		self._check_lengths(self.introns, 'intron')

		# canonical splicing
		for intron in self.introns:
			s = intron.seq_str()
			don = s[0:2]
			acc = s[-2:len(s)]
			if don not in self.dons: self.issues['donor'] = True
			if acc not in self.accs: self.issues['acceptor'] = True

		# translation checks
		self.validated = True # must be set now to use str methods
		cds = self.cds_str()
		pro = toolbox.translate_str(cds)
		start = cds[0:3]
		stop = cds[-3:len(cds)]
		if start not in self.starts: self.issues['start'] = True
		if stop not in self.stops: self.issues['stop'] = True
		for i in range(len(pro) - 1):
			if pro[i:i+1] == '*': self.issues['ptc'] = True

	def cds_str(self):
		if not self.validated: self.validate()
		seq = []
		for exon in self.cdss: seq.append(exon.seq_str())
		if self.strand == '-': seq.reverse()
		return ''.join(seq)

	def protein_str(self):
		if not self.validated: self.validate()
		return toolbox.translate_str(self.cds_str())

class ncRNA(Transcript):
	"""Class for non-coding RNAs.
	
	Extends the base `Feature` class for ncRNAs. Don't use this class
	for protein-coding mRNAs because it does not make sanity checks on
	protein-coding sequence.
	
	`ncRNAs`s are created by first instantiating this class and then using
	the `add_child()` method (from the parent `Feature` class) to add features
	of type 'exon'. All 'ncRNA`s must have an 'id'. Children will automatically
	reference this in their 'pid'.
	
	Parameters
	----------
	See the `Feature` base class, but make sure you set id.
	
	Attributes (extended from the base class)
	----------
	+ exons - list of exons
	+ inrons - list of introns
	"""

	clade = 'std' #standard
	limit = {
		'exon':   {'min':20, 'max':10000},
		'intron': {'min':30, 'max':10000},
	}
	dons = {'GT':True}
	accs = {'AG':True}

	def set_rules(self, clade='std'):
		"""Set rules for ncRNA by clade.
		
		+ dons - dictionary of donor sites (e.g. GT)
		+ accs - dictionary of acceptor sites (e.g. AG)
		+ limit - 2d dictionary of length limits

		Parameters
		----------
		+ clade `str` name of clade ('std' or 'mammal')
		"""

		if clade == 'std':
			pass # the defaults are considered standard
		elif clade == 'mammal':
			self.limit['intron'][min] = 50
			self.limit['intron'][max] = 100000
		else:
			raise NotImplemented('clade not yet supported: ' + clade)
		self.clade = clade

	def validate(self):
		"""Checks ncRNA for common errors, populating issues."""

		if self.validated: return
		if not self.id: raise FeatureError('ncRNAs must have ids')
		self._validate()

		# mRNA properties (extended constructor)
		self.is_coding = False
		self.exons = []
		self.introns = []

		for f in self.children:
			if   f.type == 'exon': self.exons.append(f)
			else: raise GenomeError('unknown type: ' + f.type)

		self.exons.sort(key = operator.attrgetter('beg'))

		# create introns from exons
		for i in range(len(self.exons)-1):
			beg = self.exons[i].end +1
			end = self.exons[i+1].beg -1
			self.introns.append(
				Feature(self.dna, beg, end, self.strand, 'intron'))

		# check for overlapping features
		self._check_overlaps(self.exons, 'exon')
		self._check_overlaps(self.introns, 'intron')

		# check for unusual lengths
		self._check_lengths(self.exons, 'exon')
		self._check_lengths(self.introns, 'intron')

		# canonical splicing
		for intron in self.introns:
			s = intron.seq_str()
			don = s[0:2]
			acc = s[-2:len(s)]
			if don not in self.dons: self.issues['donor'] = True
			if acc not in self.accs: self.issues['acceptor'] = True

class Gene(Feature):
	"""
	Class for genes, which have `Transcript` children.
	
	Extends the base `Feature` class. Genes are created by first instantiating
	a parent `Gene`. Then `add_child()` of some `Transcript` class such as
	`mRNA` or `ncRNA`. The gene must have an 'id' so that its children can
	reference it.
	
	Parameters
	----------
	See the `Feature` base class. Make sure to set id.
	"""

	def transcripts(self):
		"""Returns a list of transcripts. Just an alias for children."""

		if not self.validated:
			self.validate()
		return self.children

	def validate(self):
		"""Runs validators for self anc children."""
	
		if not self.id: raise FeatureError('genes must have ids')
		if self.pid: raise FeatureError('genes have no pids')
		self._validate()
		for f in self.children:
			if f.pid == None:
				f.pid = self.id
		self.validated = True
