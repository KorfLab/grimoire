"""
Miscellaneous functions for operating on numbers and strings.
"""

import math
import re
import gzip
import operator
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
