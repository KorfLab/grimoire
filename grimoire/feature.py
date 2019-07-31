class Performance:
	"""Class for Performance evaluation"""

	def __init__(self, model):
		"""
		Parameters
		----------
		model: object
			Model used for performance evaluation
		"""

		self.model = model
		self.nt_same = 0
		self.nt_diff = 0
		self.full_same = 0
		self.full_diff = 0
		self.feature = {} # 2D table of [type][type] = count

	def compare(self, source=None, prediction=None):
		"""
		Compare source features with prediction at various levels.

		Parameters
		----------
		source: list[object]
			List of source Feature objects
		prediction: list[object]
			List of predicted Feature objects
		"""

		# NT-level comparisons
		same, diff = 0, 0
		s, p = [], []
		for f in source:
			for i in range(f.beg, f.end): s.append(f.type)
		for f in prediction:
			for i in range(f.beg, f.end): p.append(f.type)
		for i in range(len(s)):
			if s[i] == p[i]: same += 1
			else:            diff += 1
		self.nt_same += same
		self.nt_diff += diff

		# Complete-level comparisons
		if diff == 0: self.full_same += 1
		else:         self.full_diff += 1

		# Feature-type-level comparisons
		for i in range(len(s)):
			if s[i] not in self.feature: self.feature[s[i]] = {}
			if p[i] not in self.feature[s[i]]: self.feature[s[i]][p[i]] = 0
			self.feature[s[i]][p[i]] += 1

	def report(self):
		"""Create a report. Currently prints to STDOUT"""

		print('Exact:', self.full_same / (self.full_same + self.full_diff))
		print('Accuracy:', self.nt_same / (self.nt_same + self.nt_diff))
		for s1 in self.feature:
			print(s1)
			total = 0
			for s2 in self.feature[s1]: total += self.feature[s1][s2]
			for s2 in self.feature[s1]:
				print('',s2, self.feature[s1][s2] / total)

def overlap(f1, f2):
	"""
	Determine if two features overlap

	Parameters
	----------
	f1: object
		Feature 1
	f2: object
		Feature 2
	"""

	if f1.dna.name == f2.dna.name:
		if f1.beg >= f2.beg and f1.beg <= f2.end: return True
		if f1.end >= f2.beg and f1.end <= f2.end: return True
		if f1.beg <= f2.beg and f1.end >= f2.end: return True
	return False

class Feature:
	"""Class representing a sequence feature, which may have children"""

	def __init__(self, dna, beg, end, strand, type,
			id=None, phase='.', score='.', source='.', parent_id=None):
		"""
		Parameters
		----------
		dna: object
			DNA object
		beg: int
			Beginning position of feature (1-based)
		end: int
			Ending position of feature (always >= begin)
		strand: str
			Forward (+) or reverse (-) strand
		type: str
			Feature type (e.g. gene)
		id: str
			Feature ID  (default is None)
		score: float/str
			A floating point number (default is '.')
		phase: [0, 1, 2] or .
			Coding phase (for CDS features only)
		source: str
			Data source (default is '.')
		parent_id: str
			ID of parent (default is None)
		
		Notes
		-----
		features is a list of feature objects
		genes is a list of features constructed into genes
		"""

		self.dna = dna
		self.beg = beg
		self.end = end
		self.length = end - beg + 1
		self.strand = strand
		self.type = type
		self.phase = phase
		self.id = id
		self.parent_id = parent_id
		self.score = score
		self.source = source
		self.issues = {}
		self.children = []
		self.validated = False

		if self.dna == None:
			raise GenomeError('attempt to create feature unbound to DNA')
		self._validate()

	def _validate(self, cid='child', pid='parent'):
		if self.beg < 0: self.issues['beg<0'] = True
		if self.beg > self.end: self.issues['beg>end'] = True
		if self.end > len(self.dna.seq): self.issues['end>seq'] = True
		if self.children:
			for child in self.children:
				child.validate()
				if child.beg < self.beg:
					self.issues[cid + '.beg<' + pid + '.beg'] = True
				if child.end > self.end:
					self.issues[cid + '.end>' + pid + '.end'] = True
				if child.strand != self.strand:
					self.issues['mixed_strands'] = True
				if child.issues:
					self.issues[cid + '_issues'] = True

	def validate(self):
		"""Make sure feature is validated"""

		if self.validated: return
		self._validate()
		self.validated = True

	def add_child(self, child):
		"""
		Add child to Feature

		Parameters
		----------
		child: str
			Name of child feature
		"""

		self.validated = False
		if not self.id:
			raise GenomeError('parent feature requires ID')
		else:
			self.children.append(child)

	def seq_str(self):
		"""Return a sequence of feature"""

		seq = self.dna.seq[self.beg-1:self.end]
		if self.strand == '-': seq = sequence.revcomp_str(seq)
		return seq

	def gff(self):
		"""Return feature in GFF3 format"""

		attr = ''
		if self.id and self.parent_id:
			attr = 'ID=' + self.id + ';Parent=' + self.parent_id
		elif self.id:
			attr = 'ID=' + self.id
		elif self.parent_id:
			attr = 'Parent=' + self.parent_id

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

	def __str__(self):
		return self.gff()

class mRNA(Feature):
	"""Class for mRNA"""

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
		"""
		Set rules for mRNA by clade.
		Currently, rules include boundaries for the size of introns and
		canonical signals (start, stop, splices).

		Parameters
		----------
		clade: str
			Type of clade (default is 'std or standard') Currently, 'std' and
			'mammal' are supported.
		"""

		if clade == 'std':
			pass # the defaults are considered standard
		elif clade == 'mammal':
			self.limit['intron'][min] = 50
			self.limit['intron'][max] = 100000
		else:
			raise GenomeError('clade not yet supported: ' + clade)
		self.clade = clade

	def _check_overlaps(self, f, type):
		"""
		Check for overlap issues

		Parameters
		----------
		f: object
			Feature
		type: str
			Feature type
		"""

		for i in range(1, len(f)):
			if f[i-1].end >= f[i].beg:
				self.issues['overlap_' + type] = True

	def _check_lengths(self, features, type):
		"""
		Check for length issues

		Parameters
		----------
		features: list[object]
			A list of feature objects
		type: str
			Feature type
		"""

		for f in features:
			if f.length < self.limit[type]['min']:
				self.issues['short_' + type] = True
			if f.length > self.limit[type]['max']:
				self.issues['long_' + type] = True

	def validate(self):
		"""Validate Feature and make sure nothing breaks"""

		if self.validated: return
		if not self.id: raise GenomeError('mRNAs must have ids')
		if not self.parent_id: raise GenomeError('mRNAs must have parent_ids')
		self._validate()

		# mRNA properties (extended constructor)
		self.exons = []
		self.introns = []
		self.cdss = []
		self.utr5s = []
		self.utr3s = []

		for f in self.children:
			if   f.type == 'exon': self.exons.append(f)
			elif f.type == 'CDS': self.cdss.append(f)
			else: raise GenomeError('unknown type: ' + f.type)

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
		pro = sequence.translate_str(cds)
		start = cds[0:3]
		stop = cds[-3:len(cds)]
		if start not in self.starts: self.issues['start'] = True
		if stop not in self.stops: self.issues['stop'] = True
		for i in range(len(pro) - 1):
			if pro[i:i+1] == '*': self.issues['ptc'] = True

	def tx_str(self):
		if not self.validated: self.validate()
		seq = []
		for exon in self.exons: seq.append(exon.seq_str())
		if self.strand == '-': seq.reverse()
		return ''.join(seq)

	def cds_str(self):
		if not self.validated: self.validate()
		seq = []
		for exon in self.cdss: seq.append(exon.seq_str())
		if self.strand == '-': seq.reverse()
		return ''.join(seq)

	def protein_str(self):
		if not self.validated: self.validate()
		return sequence.translate_str(self.cds_str())

class ProteinCodingGene(Feature):
	"""Class for Protein Coding Genes"""

	def mRNAs(self):
		"""Protein Coding Gene is mRNA. Validates that mRNA has all components"""

		if not self.validated:
			self.validate()
		return self.children

	def validate(self):
		if not self.id: raise GenomeError('genes must have ids')
		if self.parent_id: raise GenomeError('genes have no parent_ids')
		self._validate(cid='mRNA', pid='gene')
		self.validated = True
