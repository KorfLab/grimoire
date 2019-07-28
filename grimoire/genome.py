"""
Genome

This modeule contains classes and functions used for the representation and
validation of genomes.

"""

import re
import json
import sys
import operator
import gzip

import grimoire.toolbox as toolbox
import grimoire.sequence as sequence

class GenomeError(Exception):
	pass

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

class Genome:
	"""Class representing a genome, which has chromosomes (DNA objects)"""

	def __init__(self, species=None, fasta=None, gff3=None,
			check_alphabet=False):
		"""
		Parameters
		----------
		species: str
		   Species of genome, not really used yet
		fasta: file
			Path to fasta file, may be compressed
		gff3: file
			Path to gff3 file, may be compressed
		check_alphabet: bool
			Check whether it is the correct alphabet (good idea)
		"""

		self.species = species
		self.chromosomes = []
		ff = toolbox.FASTA_stream(fasta)
		gf = toolbox.GFF_file(gff3)

		mRNA_parts = ['CDS', 'exon']
		for entry in ff:
			chrom = sequence.DNA(name=entry.id, seq=entry.seq)
			if check_alphabet: chrom.check_alphabet()

			# convert protein-coding gene-based GFF to Features
			genes = {}
			mRNAs = {}
			parts = []
			for g in gf.get(chrom=chrom.name):
				id, pid = None, None
				im = re.search('ID=([\w\.\:]+)', g.attr)
				pm = re.search('Parent=([\w\.:]+)', g.attr)
				if im: id = im[1]
				if pm: pid = pm[1]
				if g.type == 'gene':
					genes[id] = ProteinCodingGene(chrom, g.beg, g.end, g.strand, g.type, id=id, parent_id=pid)
				elif g.type == 'mRNA':
					mRNAs[id] = mRNA(chrom, g.beg, g.end, g.strand, g.type, id=id, parent_id=pid)
				elif g.type in mRNA_parts:
					parts.append(Feature(chrom, g.beg, g.end, g.strand, g.type, id=id, parent_id=pid))

			# add parts to mRNAs
			for f in parts:
				f.validate()
				if f.parent_id in mRNAs:
					mRNAs[f.parent_id].add_child(f)

			# add mRNAs to genes
			for txid in mRNAs:
				f = mRNAs[txid]
				if f.parent_id in genes:
					genes[f.parent_id].add_child(f)

			# add genes to chromosome
			for gid in genes:
				chrom.features.append(genes[gid])

			# add chromosome to genome
			self.chromosomes.append(chrom)

class Genomic:
	"""Class for iterating through DNA objects with attached feature tables"""

	def __init__(self, fasta=None, gff=None, check_alphabet=None):
		"""
		Use path to file or file pointer.

		Parameters
		----------
		fasta: str
			Path to fasta file (may be compressed)
		gff: str
			Path to gff file (may be compressed)
		"""

		self.fp = None
		self.gz = False
		self.gff = toolbox.GFF_file(gff)
		self.check = check_alphabet

		if re.search('\.gz$', fasta):
			self.fp = gzip.open(fasta)
			self.gz = True
		else:
			self.fp = open(fasta, 'r')
		self.lastline = ''
		self.done = False

	def __iter__(self):
		return self

	def __next__(self):
		return self.next()

	def next(self):
		"""
		Retrieves the next entry of the FASTA file as DNA with features
		"""

		if self.done: raise StopIteration()
		header = None
		if self.lastline[0:1] == '>':
			header = self.lastline
		else:
			header = self.fp.readline()
			if self.gz: header = str(header, 'utf-8')

		m = re.search('>\s*(\S+)\s*(.*)', header)
		id = m[1]
		desc = m[2]
		seq = []

		while (True):
			line = self.fp.readline()
			if self.gz: line = str(line, 'utf-8')
			if line[0:1] == '>':
				self.lastline = line
				break
			if line == '':
				self.done = True
				self.fp.close()
				break

			line = line.replace(' ', '')
			seq.append(line.strip())
		
		dna = sequence.DNA(name = id, seq=''.join(seq))
		if self.check: dna.check_alphabet()

		# add features
		for g in self.gff.get(chrom=dna.name):
			id, pid = None, None
			im = re.search('ID=([\w\.\:]+)', g.attr)
			pm = re.search('Parent=([\w\.:]+)', g.attr)
			if im: id = im[1]
			if pm: pid = pm[1]
			dna.features.append(Feature(dna, g.beg, g.end, g.strand,
				g.type, source=g.source, score=g.score,
				id=id, parent_id=pid))
				
		# add genes (intentionally redudant with above: separate addresses)
		mRNA_parts = ['CDS', 'exon']
		genes = {}
		mRNAs = {}
		parts = []
		for g in self.gff.get(chrom=dna.name):
			id, pid = None, None
			im = re.search('ID=([\w\.\:]+)', g.attr)
			pm = re.search('Parent=([\w\.:]+)', g.attr)
			if im: id = im[1]
			if pm: pid = pm[1]

			if g.type == 'gene':
				genes[id] = ProteinCodingGene(dna, g.beg, g.end, 
					g.strand, g.type, id=id, parent_id=pid)
			elif g.type == 'mRNA':
				mRNAs[id] = mRNA(dna, g.beg, g.end, g.strand, g.type,
					id=id, parent_id=pid)
			elif g.type in mRNA_parts:
				parts.append(Feature(dna, g.beg, g.end, g.strand, g.type,
					id=id, parent_id=pid))
		
			# add parts to mRNAs
			for f in parts:
				f.validate()
				if f.parent_id in mRNAs:
					mRNAs[f.parent_id].add_child(f)

			# add mRNAs to genes
			for txid in mRNAs:
				f = mRNAs[txid]
				if f.parent_id in genes:
					genes[f.parent_id].add_child(f)

			# add genes to dna
			for gid in genes:
				dna.genes.append(genes[gid])
			
		return dna