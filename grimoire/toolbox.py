"""
Miscellaneous functions for operating on numbers, strings, and such.
"""

import math
import re
import gzip
import operator
import random
import sys
from functools import reduce

class ToolboxError(Exception):
	pass

def prod(iterable):
	"""Returns product of the elements in an iterable."""
	return reduce(operator.mul, iterable) # replace with math.prod when  3.8 is standard

def log(p):
	"""Returns the value in log base e with a minimum value of -999."""
	if p < 0: raise ValueError('p < 0: ' + str(p))
	if p == 0: return -999
	else:      return math.log(p)

def sumlog(v1, v2):
	"""Returns the sum of two logspaced values in logspace."""
	if v1 < v2: v1, v2 = v2, v1
	return math.log(1 + math.exp(v2 - v1)) + v1

def _kmers(alphabet, table, key, n, k, v):
	if (k == 0) :
		if key not in table:
			table[key] = v
			return

	for i in range(n):
		t = key + alphabet[i]
		_kmers(alphabet, table, t, n, k - 1, v)

def generate_kmers(alphabet='nt', k=1, pseudo=0):
	"""Returns a dictionary of all k-mers of either nucleotide or amino acid alphabet.

	Parameters
	----------
	+ alphabet= `str` 'nt' or 'aa'
	+ k=        `int` length of k-mer
	+ pseudo=   `int` pseudocount
	"""

	table = {}
	if (alphabet == 'nt') :
		_kmers(['A', 'C', 'G', 'T'], table, '', 4, k, pseudo)
	elif (alphabet == 'aa') :
		_kmers(['A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'W', 'S', 'B', 'D', 'H', 'V', 'N'],
		table, '', 20, k, pseudo)
	return table

GCODE = {
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

def revcomp_str(seq):
	"""
	Returns the reverse complement of a `str` nucleotide sequence.

	Parameters
	----------
	+ seq `str` nucleotide sequence
	"""

	comp = str.maketrans('ACGTRYMKWSBDHV', 'TGCAYRKMWSVHDB')
	anti = seq.translate(comp)[::-1]
	return anti

def translate_str(seq, table='standard'):
	"""
	Returns translated nucleotide sequence as a `str` amino acid sequence.

	Parameters
	----------
	+ seq   `str` nucleotide sequence
	+ table= `str` flavor of translation table
	"""

	pro = []
	for i in range(0, len(seq), 3):
		codon = seq[i:i+3]
		if codon in GCODE[table]: pro.append(GCODE[table][codon])
		else: pro.append('X')
	return "".join(pro)

def longest_orf(seq):
	orfs = []
	for f in range(3):
		pro = translate_str(seq[f:])
		start = 0
		while start < len(pro):
			stop = 0
			if pro[start] == 'M':
				for i, s in enumerate(pro[start+1:]):
					if s == '*':
						stop = i + start +1
						break
			if stop != 0:
				orfs.append( (pro[start:stop], start*3 + f) )
			start += 1
	orfs.sort(key=lambda t: len(t[0]), reverse=True)
	
	if len(orfs) > 0: return orfs[0]
	else:             return (None, None)

def random_dna(length, a=0.25, c=0.25, g=0.25, t=0.25):
	"""
	Generates random nucleotide sequence.

	Parameters
	----------
	+ length `int` nucleotide sequence
	+ a= `float` probability of A
	+ c= `float` probability of C
	+ g= `float` probability of G
	+ t= `float` probability of T
	"""

	assert(math.isclose(a+c+g+t, 1.0))
	seq = ''
	for i in range(length):
		r = random.random()
		if   r < a:     seq += 'A'
		elif r < a+c:   seq += 'C'
		elif r < a+c+g: seq += 'G'
		else:           seq += 'T'
	return seq

def get_filepointer(filename):
	"""
	Returns a filepointer to a file based on file name (.gz or - for stdin).
	"""

	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	elif filename == '-':          fp = sys.stdin
	else:                          fp = open(filename)
	return fp

def read_fasta(filename):
	"""
	Simple fasta reader that returns name, seq for a filename.

	Parameters
	----------
	+ filename
	"""

	name = None
	seqs = []

	fp = get_filepointer(filename)

	while True:
		line = fp.readline()
		if line == '': break
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()
