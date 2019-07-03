
import re
import json
import sys
import operator

import toolbox
import sequence

class GenomeError(Exception):
	pass


class Feature:
	"""Class representing a sequence feature, which may have children"""

	def __init__(self, dna, beg, end, strand, type,
			id=None, score='.', source='.', parent_id=None):
		self.dna = dna
		self.beg = beg
		self.end = end
		self.strand = strand
		self.type = type
		self.id = id
		self.parent_id = parent_id
		self.score = score
		self.source = source
		self.issues = []
		self.children = []
		
		if beg < 0: self.issues.append('beg < 0')
		if beg > end: self.issues.append('beg > end')
		if end > len(self.dna.seq): self.issues.append('end out of range:' + end)

	def add_child(self, child):
		if not self.id:
			raise GenomeError('parent feature without ID')
		else:
			self.children.append(child)
			# make sure children are inside parent?

	def seq_str(self):
		seq = self.dna.seq[self.beg-1:self.end]
		if self.strand == '-': seq = sequence.revcomp_str(seq)
		return seq

	def gff(self):
		attr = ''
		if self.id and self.parent_id:
			attr = 'ID=' + self.id + ';Parent=' + self.parent_id
		elif self.id:
			attr = 'ID=' + self.id
		elif self.parent_id:
			attr = 'Parent=' + self.parent_id
		
		string = '\t'.join([self.dna.name, self.source, self.type,
			str(self.beg), str(self.end), str(self.score),
			self.strand, '.', attr])
		if self.children:
			string += '\n'
			for child in self.children:
				string += child.gff() + '\n'
			string += '\n'
			
		return string

class mRNA:
	"""Class representing a mRNA, contains exons and such"""

	def __init__(self, id=None, features=None, rules='std'):
			
		## sanity checks for chrom, beg, end, strand
		self.issues = []
		dna = features[0].dna
		strand = features[0].strand
		for f in features:
			if f.issues: self.issues.append('feature issue')
			if f.dna != dna: self.issues.append('mixed chroms')
			if f.strand != strand: self.issues.append('mixed strands')
			if f.beg > f.end: self.issues.append('beg > end')
		self.dna = dna
		self.strand = strand
		
		## build mRNA from features
		self.exons = []
		self.cdss = []
		self.utr5s = []
		self.utr3s = []
		self.id = id
		for f in features:
			if f.type == 'exon': self.exons.append(
				Feature(f.dna, f.beg, f.end, f.strand, f.type))
			elif f.type == 'CDS': self.cdss.append(
				Feature(f.dna, f.beg, f.end, f.strand, f.type))
			elif f.type == 'five_prime_UTR': self.utr5s.append(
				Feature(f.dna, f.beg, f.end, f.strand, f.type))
			elif f.type == 'three_prime_UTR': self.utr3s.append(
				Feature(f.dna, f.beg, f.end, f.strand, f.type))
			else:
				raise GenomeError('unknown type in feature table')
				
		# sort lists
		self.exons.sort(key = operator.attrgetter('beg'))
		self.cdss.sort(key = operator.attrgetter('beg'))
		self.utr5s.sort(key = operator.attrgetter('beg'))
		self.utr3s.sort(key = operator.attrgetter('beg'))
		
		# check for overlapping features
		if self.overlaps(self.exons): self.issues.append('exon overlaps')
		if self.overlaps(self.cdss): self.issues.append('cds overlaps')
		if self.overlaps(self.utr5s): self.issues.append('utr5 overlaps')
		if self.overlaps(self.utr3s): self.issues.append('utr3 overlaps')

		# assignments
		self.beg = self.exons[0].beg
		self.end = self.exons[-1].end
		
		# create introns
		self.introns = []
		for i in range(len(self.exons)-1):
			beg = self.exons[i].end +1
			end = self.exons[i+1].beg -1
			self.introns.append(
				Feature(self.dna, beg, end, self.strand, 'intron'))

		# gene structure warnings
		min_intron = 30
		max_intron = 10000
		dons = {'GT':True}
		accs = {'AG':True}
		starts = {'ATG':True}
		stops = {'TAA':True, 'TGA':True, 'TAG':True}
		if rules == 'std':
			pass # don't override parameters above
		elif rules == 'mammal':
			min_intron = 50
			max_intron = 50000
		else:
			pass # there ought to be other rule sets...
			
	
		# canonical splicing and intron length
		for intron in self.introns:
			s = intron.seq_str()
			don = s[0:2]
			acc = s[-2:len(s)]
			if don not in dons: self.issues.append('don:' + don)
			if acc not in accs: self.issues.append('acc:' + acc)
			if len(s) < min_intron or len(s) > max_intron:
				self.issues.append('int_len:' + str(len(s)))
			
		# translation checks
		cds = self.cds_str()
		pro = sequence.translate_str(cds)
		start = cds[0:3]
		stop = cds[-3:len(cds)]
		if start not in starts: self.issues.append('start:' + start)
		if stop not in stops: self.issues.append('stop:' + stop)
		for i in range(len(pro) - 1):
			if pro[i:i+1] == '*': self.issues.append('ptc:' + str(i))
	
	def tx_str(self):
		seq = []
		for exon in self.exons: seq.append(exon.seq_str())
		if self.strand == '-': seq.reverse()
		return ''.join(seq)
	
	def cds_str(self):
		seq = []
		for exon in self.cdss: seq.append(exon.seq_str())
		if self.strand == '-': seq.reverse()
		return ''.join(seq)

	def protein_str(self):
		return sequence.translate_str(self.cds_str())

	def overlaps(self, list):
		if len(list) > 1:
			for i in range(len(list)- 1):
				if list[i].end >= list[i+1].beg: return True
		return False

class ncRNA:
	"""Class representing non-coding RNAs"""
	pass

class Gene:
	"""Class representing a gene, which contains mRNA children"""

	def __init__(self, id=None, transcripts=None):
		self.id = id
		self.beg = None
		self.end = None
		self.strand = None
		self.transcripts = transcripts
		self.issues = []
		
		self.beg = transcripts[0].beg
		self.end = transcripts[0].end
		self.strand = transcripts[0].strand
		for tx in transcripts:
			if tx.beg < self.beg: self.beg = tx.beg
			if tx.end > self.end: self.end = tx.end
			if tx.strand != self.strand: self.issues.append('mixed strands')
			if tx.issues: self.issues.append('tx issue')

class Genome:
	"""Class representing a genome, which has chromosomes"""

	def __init__(self, species=None, fasta=None, gff3=None, check_alphabet=False):
		self.species = species
		self.chromosomes = []
		ff = toolbox.FASTA_file(fasta)
		gf = toolbox.GFF_file(gff3)
		mRNA_parts = ['CDS', 'exon', 'five_prime_UTR', 'three_prime_UTR']
		
		for id in ff.ids:
			entry = ff.get(id)
			chrom = sequence.DNA(name=entry.id, seq=entry.seq)
			if check_alphabet: chrom.check_alphabet()

			# convert all GFF to features
			genes = {}
			mRNAs = {}
			parts = []
			for g in gf.get(chrom=chrom.name):
				id, pid = None, None
				im = re.search('ID=([\w\.]+)', g.attr)
				pm = re.search('Parent=([\w\.]+)', g.attr)
				if im: id = im[1]
				if pm: pid = pm[1]
				if g.type == 'gene':
					genes[id] = Feature(chrom, g.beg, g.end, g.strand, g.type, id=id, parent_id=pid)
				elif g.type == 'mRNA':
					mRNAs[id] = Feature(chrom, g.beg, g.end, g.strand, g.type, id=id, parent_id=pid)
				elif g.type in mRNA_parts:
					parts.append(Feature(chrom, g.beg, g.end, g.strand, g.type, id=id, parent_id=pid))
			
			# add parts to mRNAs
			for f in parts:
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
"""
			ft = {}	
			for part in gene_parts:
				for entry in gf.get(type=part, chrom=chrom.name):
					pid = re.search('Parent=([\w\.]+)', entry.attr)[1]
					if pid not in ft: ft[pid] = []
					ft[pid].append(Feature(chrom, entry.beg, entry.end,
						entry.strand, entry.type))
			gt = {}
			for entry in gf.get(type='mRNA', chrom=chrom.name):
				tid = re.search('ID=([\w\.]+)', entry.attr)[1]
				gid = re.search('Parent=([\w\.]+)', entry.attr)[1]
				if gid not in gt: gt[gid] = {}
				gt[gid][tid] = entry
				
			for gid in gt:
				txs = []
				for tid in gt[gid]:
					tx = mRNA(id=tid, features=ft[tid])
					txs.append(tx)
				chrom.features.append(Gene(id=gid, transcripts=txs))
			self.chromosomes.append(chrom)
"""

		