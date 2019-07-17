"""
Sequence

This script contains classes and functions used to make/interpret sequences.

The following functions are provided in HMM:
 	* generate_kmers
	* revcomp_str
	* translate_str

The following classes and methods are provided in HMM:
	* BioSequence
		* BioSequence.fasta
	* DNA
		* DNA.check_alphabet
		* DNA.revcomp
		* DNA.translate
	* Protein
		* Protein.check_alphabet
"""

class SequenceError(Exception):
	pass

##########################
### naked definitions  ###
##########################

def _kmers(alphabet, table, key, n, k, v):
	if (k == 0) :
		if key not in table:
			table[key] = v
			return

	for i in range(n):
		t = key + alphabet[i]
		_kmers(alphabet, table, t, n, k - 1, v)

def generate_kmers(alphabet='nt', k=1, pseudo=0):
	"""Creates a dictionary of all kmers of either nt
	or aa alphabet.

	Parameters
	----------
	alphabet: str
		The alphabet used. Either 'nt' or 'aa'
		(default is 'nt')
	k: int
		Any integer. Placeholder for the value in kmer
	Pseudo:
		Pseudocount value (default is 0)
	"""

	table = {}
	if (alphabet == 'nt') :
		_kmers(DNA.canonical, table, '', 4, k, pseudo)
	elif (alphabet == 'aa') :
		_kmers(Protein.canonical, table, '', 20, k, pseudo)
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
	Makes reverse compliment sequence

	Parameters
	----------
	seq: str
		Given sequence
	"""

	comp = str.maketrans('ACGTRYMKWSBDHV', 'TGCAYRKMWSVHDB')
	anti = seq.translate(comp)[::-1]
	return anti

def translate_str(seq, table='standard'):
	"""
	Translate sequence into protein

	Parameters
	----------
	seq: str
		Given sequence
	table: str
		Given translation table (default is 'standard')
	"""

	pro = []
	for i in range(0, len(seq), 3):
		codon = seq[i:i+3]
		if codon in GCODE[table]: pro.append(GCODE[table][codon])
		else: pro.append('X')
	return "".join(pro)

###############
### Classes ###
###############

class BioSequence:
	"""Generic parent class of biological sequences"""

	def fasta(self, wrap=80):
		"""
		Convert FASTA format to objects

		Parameters
		----------
		wrap: int
			Number of characters before wrapping/
			new line
		"""

		s = '>'
		if self.name: s += self.name
		if self.desc: s += ' ' + self.desc
		if self.species: s += '[' + self.species + ']'
		s += '\n'
		for i in range(0, len(self.seq), wrap):
			s += self.seq[i:i+wrap] + '\n'
		return s

	def __str__(self):
		return self.fasta()

class DNA(BioSequence):
	"""Class for DNA sequences. IUPAC alphabet. Uppercase only"""

	canonical = ['A', 'C', 'G', 'T']
	extended = ['A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'W', 'S', 'B', 'D', 'H', 'V', 'N']

	def __init__(self, name=None, seq=None, desc=None, species=None):
		"""
		Parameters
		----------
		name: str
			Name of Sequence (default is None)
		seq: str
			The actual sequence (default is None)
		desc: str
			Description of sequence (default is None)
		species: str
			Specie of sequence (default is None)
		"""

		self.name = name
		self.seq = seq
		self.desc = desc
		self.species = species
		self.features = []

	def check_alphabet(self):
		"""Check if sequence is in given alphabet"""

		for i in range(len(self.seq)):
			nt = self.seq[i:i+1]
			if nt not in self.extended:
				raise SequenceError('letter not in DNA alphabet: ' + nt)

	def revcomp(self):
		"""Return reverse compliment sequence"""
		anti = revcomp_str(self.seq)
		return DNA(seq=anti)

	def translate(self, table='standard'):
		"""Return translated protein sequence"""
		pro = translate_str(self.seq)
		return Protein(seq=pro)

class Protein(BioSequence):
	"""Class for protein sequences. Uppercase only. 20 aa + X and *"""

	canonical = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
		'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	extended = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
		'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X', '*']

	def __init__(self, name=None, seq=None, desc=None, species=None):
		"""
		Parameters
		----------
		name: str
			Name of Sequence (default is None)
		seq: str
			The actual sequence (default is None)
		desc: str
			Description of sequence (default is None)
		species: str
			Specie of sequence (default is None)
		"""

		self.name = name
		self.seq = seq
		self.desc = desc
		self.species = species

	def check_alphabet(self):
		"""Check if sequence is in given alphabet"""
		for i in range(len(self.seq)):
			aa = self.seq[i:i+1]
			if aa not in self.extended:
				raise SequenceError('letter not in protein alphabet: ' + aa)
