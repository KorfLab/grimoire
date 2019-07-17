"""
HMM(Hidden Markov Model)

This script contains classes and functions that used to create HMM.
It is assumed that the sequences inputted are in the list format.

The following functions are provided in HMM:
 	* emission_model
	* train_emission
	* train_emissions
	* train_cds
	* state_factory
	* null_state_factory
	* connect_all
	* connect2

The following classes and methods are provided in HMM:
	* State - Class for HMM States
	* HMM - Class for HMMs
		* HMM.read
		* HMM.write
		* HMM.convert2log
		* HMM.macro_labels
"""
#Documentation done in NumPy/SciPy format.

import json
import re
import gzip
import sys

import grimoire.sequence as sequence
import grimoire.toolbox as toolbox

class HMMError(Exception):
	pass

def emission_model(context=1, alphabet='nt'): #Why first context then alphabet?
	"""
	Makes an empty emission model with specified context and alphabet.

	Table is in dictionary format.

    Parameters
    ----------
    context: int
        The level of context, or the number of previous states taken into
		account when calculating the emission probabilities. Value must be a
		nonzero whole number (default is 1)
    alphabet: str
		The type of alphabet for model. Takes input nucleotides 'nt' or amino
		acids 'aa'. If other strings are inputted, an alphabet error is raised.
		(default is 'nt')
    """

	letters = None
	if alphabet == 'nt': letters = sequence.DNA.canonical
	elif alphabet == 'aa': letters = sequence.Protein.canonical
	else: raise HMMError('alphabet error: ', alphabet)

	if (context > 0):
		table = sequence.generate_kmers(alphabet=alphabet, k=context)
		for k in table:
			table[k] = {}
			for a in letters:
				table[k][a] = 0
	elif (context == 0) :
		table = sequence.generate_kmers(alphabet=alphabet, k=1, pseudo=0)
	else :
		raise HMMError('negative context')

	return table

def train_emission(seqs, context=0):
	"""
	Train the nucleotides emission model with customized context.

	This function outputs the emission frequencies as a dictionary.

    Parameters
    ----------
	seqs: list
		A list of all the sequences used in training emissions
    context: int
        The level of context, or the number of previous states taken into
		account when calculating the emission probabilities. Value must be a
		nonzero whole number (default is 1)
    """

	count = emission_model(context=context)
	freq = {}

	if context == 0:
		total = 0
		for seq in seqs:
			for i in range(len(seq)):
				nt = seq[i:i+1]
				if (nt in count):
					count[nt] += 1
					total += 1
		for nt in count: freq[nt] = round(count[nt] / total, 4)
	else:
		for seq in seqs:
			for i in range(len(seq) -context):
				pos = i + context
				ctx = seq[i:i+context]
				nt = seq[pos:pos+1]
				if ctx in count and nt in count[ctx]:
					count[ctx][nt] += 1
		for ctx in count:
			total = 0
			freq[ctx] = {}
			for nt in count[ctx]: total += count[ctx][nt]
			for nt in count[ctx]: freq[ctx][nt] = round(count[ctx][nt]/total, 4)

	return freq

def train_emissions(seqs, context=0):
	"""
	Train nucleotides emission models with customized context.

	This function outputs the emission frequencies as a list of dictionaries.
	If context = 0, outputs list with a single dictionary.
	If context > 0, outputs list with multiple dictionaries.

    Parameters
    ----------
	seqs: list
		A list of all the sequences used in training emissions
    context: int
        The level of context, or the number of previous states taken into
		account when calculating the emission probabilities. Value must be a
		nonzero whole number (default is 1)
    """

	counts = []
	freqs = []
	for i in range(len(seqs[0])):
		counts.append(emission_model(context=context))

	if (context == 0):
		for seq in seqs:
			for i in range(len(seq)):
				nt = seq[i:i+1]
				if (nt in counts[i]): counts[i][nt] += 1
		for count in counts:
			total = 0
			freq = {}
			for nt in count: total += count[nt]
			for nt in count: freq[nt] = round(count[nt] / total	, 4)
			freqs.append(freq)
	else:
		for i in range(context):
			for ctx in counts[i]:
				for nt in counts[i][ctx]:
					counts[i][ctx][nt] = 1 # they will all be 0.25
		for seq in seqs:
			for i in range(len(seq) - context):
				pos = i + context
				ctx = seq[i:i+context]
				nt = seq[i + context : i + context + 1]
				if (ctx in counts[pos] and nt in counts[pos][ctx]):
					counts[pos][ctx][nt] += 1
		for i in range(len(counts)):
			freqs.append({})
			for ctx in counts[i]:
				total = 0
				freqs[i][ctx] = {}
				for nt in counts[i][ctx]: total += counts[i][ctx][nt]
				for nt in counts[i][ctx]:
					if total != 0:
						freqs[i][ctx][nt] = round(counts[i][ctx][nt]/total, 4)
					else:
						freqs[i][ctx][nt] = 0
	return freqs

def train_cds(seqs, context=0):
	"""
	Train the coding sequence (CDS) emission model with customized context.

	This function outputs the emission frequencies as a list of dictionaries.
	If context = 0, outputs list with a single dictionary.
	If context > 0, outputs list with multiple dictionaries.

    Parameters
    ----------
	seqs: list
		A list of all the sequences used in training emissions
    context: int
        The level of context, or the number of previous states taken into
		account when calculating the emission probabilities. Value must be a
		nonzero whole number (default is 1)
    """

	counts = []
	freqs = []
	for i in range(3):
		counts.append(emission_model(context=context))

	if (context == 0):
		for seq in seqs:
			for i in range(0, len(seq) -3, 3):
				for j in range(3):
					nt = seq[i+j]
					if (nt in counts[j]) : counts[j][nt] += 1
		for count in counts:
			total = 0
			freq = {}
			for nt in count: total += count[nt]
			for nt in count: freq[nt] = round(count[nt] / total	, 4)
			freqs.append(freq)
	else:
		for seq in seqs:
			#print(seq[0:10])
			for i in range(0, len(seq) -context -3, 3):
				for j in range(3):
					pos = i + j + context
					ctx = seq[i+j:i+j+context]
					nt = seq[i + j + context]
					#print('>>', i, j, pos, ctx, nt)
					if (ctx in counts[j] and nt in counts[j][ctx]) :
						counts[j][ctx][nt] += 1
		for j in range(3):
			freqs.append({})
			for ctx in counts[j]:
				total = 0
				freqs[j][ctx] = {}
				for nt in counts[j][ctx]: total += counts[j][ctx][nt]
				for nt in counts[j][ctx]:
					if total != 0:
						freqs[j][ctx][nt] = round(counts[j][ctx][nt]/total, 4)
					else:
						freqs[j][ctx][nt] = 0

	return freqs

class State:
	"""Class for HMM States"""

	def __init__(self, name=None, context=None, emits=None, init=0, term=0, next=None):
		"""
		Parameters
		----------
		name: str
			Name of State (default is None)
		context: int
			The level of context, or the number of previous states taken into
			account when calculating the emission probabilities. Value must be a
			nonzero whole number (default is None)
		emits: float?
			Singular emission probability for state? (default is None)
		init:  float
			Probability of having the HMM start at this state (default is 0)
		term: float
			Probability of having the HMM end at this state (default is 0)
		next: dictionary
			Dictionary of all next states (default is None)
		"""

		self.name = name
		self.ctxt = context
		self.init = init
		self.term = term
		self.emit = emits
		self.next = {}

	@classmethod
	def from_json(cls, json_string):
		state = cls()
		state.__dict__ = json.loads(json_string)

	def to_json(self):
		return(json.dumps(self.__dict__, indent = 4))

def state_factory(prefix, emissions):
	"""
	Create a list of state objects.

	Parameters
	----------
	prefix: str
		State prefix (e.g. ACC, DON, GEN)
	emissions: list
		A list of emission probabilities
	"""

	state_list = []
	for i in range(len(emissions)):
		if type(emissions[i]) == type(dict()):
			if ('A' in emissions[i]):
				if type(emissions[i]['A']) == type(dict()):
					context = 1
				else:
					context = 0
			else:
				key = next(iter(emissions[i]))
				context = len(key)
		state = State(name=prefix + '-' + str(i),
			emits=emissions[i],
			context=context)
		state_list.append(state)
	return(state_list)

def null_state_factory(file=None, context=None):
	"""
	Create a list of null state objects.

	Parameters
	----------
	file: str
		File name
	context: int
        The level of context, or the number of previous states taken into
		account when calculating the emission probabilities. Value must be a
		nonzero whole number (default is None)
	"""

	fasta = toolbox.FASTA_stream(filename=file)
	count = emission_model(context=context)
	freq = {}

	if context == 0:
		total = 0
		for entry in fasta:
			seq = entry.seq
			for i in range(len(seq)):
				nt = seq[i:i+1]
				if (nt in count):
					count[nt] += 1
					total += 1
		for nt in count: freq[nt] = round(count[nt] / total, 4)
	else:
		for entry in fasta:
			seq = entry.seq
			for i in range(len(seq) -context):
				pos = i + context
				ctx = seq[i:i+context]
				nt = seq[pos:pos+1]
				if ctx in count and nt in count[ctx]:
					count[ctx][nt] += 1
		for ctx in count:
			total = 0
			freq[ctx] = {}
			for nt in count[ctx]: total += count[ctx][nt]
			for nt in count[ctx]: freq[ctx][nt] = round(count[ctx][nt]/total, 4)

	return State(name='null', context=context, emits=freq)


def connect_all (states) :
	"""
	Connects all of the state objects.
	Updates the dictionary of next state objects

	Parameters
	----------
	states: list
		A list of state objects
	"""

	for i in range(len(states)-1) :
		states[i].next[states[i+1].name] = 1

def connect2 (s1, s2, p) :
	"""
	Connects 2 state objects.
	Updates the dictionary of next state objects

	Parameters
	----------
	s1: object
		State 1
	s2: object
		State 2
	p: float
		transition probability
	"""

	s1.next[s2.name] = p

class HMMdecoder(json.JSONEncoder):
	def default(self, o):
		return o.__dict__

class HMM:
	"""Class for HMMs"""
	def __init__(self, name=None, log=False, states=None, null=None):
		"""
		Parameters
		----------
		name: str
		 	Name of HMM
		log: bool
			HMM is in Logspace
		states: list
			List of HMM states
		null: object
			State representing null model
		"""

		self.name = name
		self.log = log
		self.states = states
		self.null = null

	@classmethod
	def read(cls, filename):
		"""
		Read in HMM file

		Parameters
		----------
		filename: str
			Name of HMM file
		"""

		d = None
		if re.search('\.gz$', filename):
			with gzip.open(filename, mode='r') as fp:
				d = json.loads(fp.read())
		else:
			with open(filename, 'r') as fp:
				d = json.loads(fp.read())
		hmm = HMM()
		hmm.name = d['name']
		hmm.log = d['log']
		hmm.states = []
		for s in d['states']:
			st = State()
			st.name = s['name']
			st.init = s['init']
			st.term = s['term']
			st.ctxt = s['ctxt']
			st.emit = s['emit']
			st.next = s['next']
			hmm.states.append(st)
		hmm.null = State()
		hmm.null.name = d['null']['name']
		hmm.null.init = d['null']['init']
		hmm.null.term = d['null']['term']
		hmm.null.ctxt = d['null']['ctxt']
		hmm.null.emit = d['null']['emit']
		hmm.null.next = d['null']['next']
		return hmm

	def write(self, filename):
		"""
		Write in HMM file

		Parameters
		----------
		filename: str
			Name of HMM file
		"""

		if re.search('\.gz$', filename):
			with gzip.open(filename, mode='w') as fp:
				fp.write(json.dumps(self.__dict__, indent=4, cls=HMMdecoder))
		else:
			with open(filename, 'w+') as fp:
				fp.write(json.dumps(self.__dict__, indent=4, cls=HMMdecoder))

	def convert2log(self):
		"""
		Converts the HMM into logspace
		"""

		if self.log:
			raise HMMError('model already in log space')

		for state in self.states + [self.null]:
			state.init = toolbox.log(state.init)
			state.term = toolbox.log(state.term)
			if state.ctxt == 0:
				for nt in state.emit:
					state.emit[nt] = toolbox.log(state.emit[nt])
			else:
				for ctx in state.emit:
					for nt in state.emit[ctx]:
						state.emit[ctx][nt] = toolbox.log(state.emit[ctx][nt])
			for next in state.next:
				state.next[next] = toolbox.log(state.next[next])
		self.log = True

	def macro_labels(self):
		"""
		Returns a list of label names and strips off the probability attributes
		Output is not a true list.
		"""
		label = {}
		for state in self.states:
			macro = re.search('(\w+)', state.name)[1]
			label[macro] = True
		return list(label.keys())
