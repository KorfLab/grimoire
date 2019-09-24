"""
Classes for building and decoding HMMs and variants.
"""

import json
import re
import gzip
import sys
import math
import copy
import random

import grimoire.sequence as sequence
import grimoire.toolbox as toolbox
import grimoire.io as io
from grimoire.feature import Feature

class HMMError(Exception):
	pass

def emission_model(context=0, alphabet='nt'):
	"""
	Returns an empty emission table `dict` with specified context and alphabet.
	
	Parameters
	----------
	+ context=  `int` order of the Markov model (0 or greater)
	+ alphabet= `str` either 'nt' or 'aa' (but really 'aa'?)
	"""

	letters = None
	if alphabet == 'nt': letters = sequence.DNA.canonical
	elif alphabet == 'aa': letters = sequence.Protein.canonical
	else: raise HMMError('alphabet error: ', alphabet)

	if (context > 0):
		table = toolbox.generate_kmers(alphabet=alphabet, k=context)
		for k in table:
			table[k] = {}
			for a in letters:
				table[k][a] = 0
	elif (context == 0) :
		table = toolbox.generate_kmers(alphabet=alphabet, k=1, pseudo=0)
	else :
		raise HMMError('negative context')

	return table

def train_emission(seqs, context=0):
	"""
	Train an emission table for a single sequence (e.g. UTR).

	Parameters
	----------
	+ seqs    `list` list of sequences of type `str`
	+ context= `int` order of the Markov model (0 or greater)
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
	Train a list of emission tables. Usually for use with connected states.

	Parameters
	----------
	+ seqs     `list` list of sequences of type `str`
	+ context= `int` order of the Markov model (0 or greater)
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
	Train the coding sequence (CDS) emission model. Using a context of 2 or
	greater will ensure there are no in-frame stop codons (provided the
	sequences don't have in-frame stop codons).

	Parameters
	----------
	+ seqs     `list` list of sequences of type `str`
	+ context= `int`  order of the Markov model (0 or greater)
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
			for i in range(context, len(seq), 1):
				ctx = seq[i-context:i]
				nt = seq[i]
				frame = i % 3
				if (ctx in counts[frame] and nt in counts[frame][ctx]) :
						counts[frame][ctx][nt] += 1
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

def state_factory(prefix, emissions):
	"""
	Create a list of related state objects (e.g. splice donor site).

	Parameters
	----------
	+ prefix    `str`  name of the state (e.g. ACC, DON)
	+ emissions `list` of emission tables produced by `train_emissions()`
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

def null_state_factory(file=None, context=0):
	"""
	Create the null state object (usually from entire genome).

	Parameters
	----------
	+ file=    `str` path to fasta file, which may be gzipped
	+ context= `int` order of the Markov model (0 or greater)
	"""

	fasta = io.FASTA_stream(filename=file)
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
	Connects a list of state objects sequentially with probability 1. The first and
	last states are unconnected.

	Parameters
	----------
	+ states `list` state objects, probably all related (e.g. donors)
	"""

	for i in range(len(states)-1) :
		states[i].next[states[i+1].name] = 1

def connect2 (s1, s2, p) :
	"""
	Connects 2 state objects with a given probability.

	Parameters
	----------
	+ s1 `State` state 1
	+ s2 `State` state 2
	+ p  `float` transition probability
	"""

	s1.next[s2.name] = p

class State:
	"""Class representing HMM States and their associated transitions."""

	def __init__(self, name=None, context=None, emits=None, init=0, term=0, next=None):
		"""
		Parameters
		----------
		+ name=    `str`   name of the state, must be unique in the model
		+ context= `int`   order of the Markov emission model (0 or greater)
		+ emits=   `obj`   emission table/dict
		+ init=    `float` probability to start in this state
		+ term=    `float` probability to end in this state
		+ next=    `dict`  dictionary of connected states and probabilities
		"""

		self.name = name
		self.ctxt = context
		self.init = init
		self.term = term
		self.emit = emits
		self.next = {}

	@classmethod
	def from_json(cls, json_string):
		"""Generates `State` object from JSON-formatted string."""
		state = cls()
		state.__dict__ = json.loads(json_string)

	def to_json(self):
		"""Returns JSON-formatted `str` representing `State` object."""
		return(json.dumps(self.__dict__, indent = 4))

class HMMdecoder(json.JSONEncoder):
	"""Used for json-based io"""
	def default(self, o):
		return o.__dict__

class HMM:
	"""Class for HMMs, which is mostly a container for `State` objects."""
	
	def __init__(self, name=None, log=False, states=None, null=None):
		"""
		Parameters
		----------
		+ name=   `str`   name of the HMM (arbitrary)
		+ log=    `bool`  flag to indicate if the model is in logspace
		+ states= `list`  list of `State` objects
		+ null=   `State` the null state used for calculated scores
		"""

		self.name = name
		self.log = log
		self.states = states
		self.null = null

	@classmethod
	def read(cls, filename):
		"""
		File-based constructor for HMMs that are in the file system.

		Parameters
		----------
		+ filename `str` path to the HMM file, which may be gzipped
		"""

		d = None
		if re.search(r'\.gz$', filename):
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
		Write the HMM file to the file system.

		Parameters
		----------
		+ filename `str` path to the outputfile (.gz will compress)
		"""

		if re.search(r'\.gz$', filename):
			with gzip.open(filename, mode='w') as fp:
				fp.write(json.dumps(self.__dict__, indent=4, cls=HMMdecoder).encode())
		else:
			with open(filename, 'w+') as fp:
				fp.write(json.dumps(self.__dict__, indent=4, cls=HMMdecoder))

	def convert2log(self):
		"""
		Converts the HMM into logspace, which is critical for long sequences.
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
		Returns a list of label names stripped of the numeric identifier.
		For example [ACC-0, ACC-1] becomes ACC.
		"""
		label = {}
		for state in self.states:
			macro = re.search(r'(\w+)', state.name)[1]
			label[macro] = True
		return list(label.keys())

class DecodeError(Exception):
	pass

class Parse:
	"""Class representing a parse through a sequence."""
	
	def __init__(self, score=None, path=None, freq=None, decoder=None):
		"""
		Parameters
		----------
		+ score= `float` usually a log-odds ratio vs. null model parse
		+ path=  `list`  a list of state labels
		+ freq=  `float` frequency of taking this path if there are multiples
		+ decoder= `obj` decoder object that produced the path
		"""

		self.score = score
		self.path = path
		self.freq = freq
		self.decoder = decoder

	def features(self):
		"""Returns parse as a list of `Feature` objects."""

		dna = self.decoder.dna
		labels = self.decoder.model.macro_labels()

		mypath = []
		for name in self.path:
			found = False
			for label in labels:
				if re.search(label, name):
					mypath.append(label)
					found = True
					break
			if not found:
				raise DecodeError('name not matched to label: ' + name)

		features = []
		beg, end = 0, 0
		while (end < len(mypath)):
			for i in range(beg, len(mypath)):
				if mypath[beg] != mypath[i]:
					end = i-1
					features.append(Feature(dna, beg+1, end+1, '+',
						mypath[beg]))
					beg = end+1
					break
				elif i == len(mypath) -1:
					features.append(Feature(dna, beg+1, i+1, '+',
						mypath[beg]))
					end = len(mypath) +1
		return features

class HMM_NT_decoder:
	"""Base class for the HMM nucleotide decoders."""
	
	def __init__(self, dna=None, model=None):
		"""
		Parameters
		----------
		+ model= `HMM`
		+ dna=   `DNA`
		"""
		
		self.dna = dna
		self.model = model
		self.null_score = None

	def _state_map(self):
		"""
		Creates a dictionary to look up states by name (used internally)
		key: state name (str), value: state (object)
		"""

		smap = {}
		for state in self.model.states + [self.model.null]:
			smap[state.name] = state
		return smap

	def _transition_map(self):
		"""
		Reverse of the state.next dictionary (used internally)
		"""

		tmap = {}
		for state in self.model.states:
			tmap[state.name] = {}
		for state in self.model.states:
			for next in state.next:
				if next not in tmap:
					tmap[next] = {}
				tmap[next][state.name] = state.next[next]
		return tmap

	def _set_null_score(self):
		"""Generate/set the score for the null model. Must be set once."""

		self.null_score = 0 if self.model.log else 1
		for i in range(len(self.dna.seq)):
			ep = self.emission(self.model.null.name, i)
			if self.model.log: self.null_score += ep
			else:              self.null_score *= ep

	def emission(self, state_name, i):
		"""
		Returns the emission probability for any state at any position.

		Parameters
		----------
		+ state_name `str` name of the state
		+ i          `int` position in the sequence (1-based coordinates)
		"""

		nt = self.dna.seq[i:i+1]
		state = self.smap[state_name]

		if i < 0: raise DecodeError('position less than 0')

		if self.model.log:
			if i < state.ctxt: return math.log(0.25)
			if state.ctxt == 0:
				if nt in state.emit: return state.emit[nt]
				else: return math.log(0.25)
			else:
				ctx = self.dna.seq[i-state.ctxt:i]
				if nt in state.emit[ctx]: return state.emit[ctx][nt]
				else: return math.log(0.25)
		else:
			if i < state.ctxt: return 0.25
			if state.ctxt == 0:
				if nt in state.emit: return state.emit[nt]
				else: return 0.25
			else:
				ctx = self.dna.seq[i-state.ctxt:i]
				if nt in state.emit[ctx]: return state.emit[ctx][nt]
				else: return 0.25
		raise DecodeError('not possible')

	def score_path(self, path):
		"""
		Returns the score for a given path

		Parameters
		----------
		+ path `list` of state names
		"""
		
		# should this be state.next or tmap?
		if self.model.log:
			score = self.smap[path[0]].init + self.emission(path[0], 0)
			for i in range(1, len(path)):
				ep = self.emission(path[i], i)
				tp = self.tmap[path[i]][path[i-1]]
				score += ep + tp
			score += self.smap[path[-1]].term
		else:
			score = self.smap[path[0]].init * self.emission(path[0], 0)
			for i in range(1, len(path)):
				ep = self.emission(path[i], i)
				tp = self.tmap[path[i]][path[i-1]]
				score *= ep * tp
			score *= self.smap[path[-1]].term
		
		return score
		
	def _inspect(self, field, beg=None, end=None):
		# used internally for debugging
		if not beg: beg = 0
		if not end: end = len(self.matrix)

		# print numbers
		print('{:<6s}'.format(''), end='')
		for col in range(beg, end):
			print('{:<10d}'.format(col), end='')
		print()

		# print letters
		print('{:<6s}'.format(''), end='')
		for col in range(beg, end):
			nt = self.dna.seq[col:col+1]
			print('{:<10s}'.format(nt), end='')
		print()

		# print scores
		for state in self.model.states:
			print('{:<6s}'.format(state.name), end='')
			for i in range(beg, end):
				#if display == 'score':
				try:
					print('{:<10.3g}'.format(self.matrix[i][state.name][field]), end='')
				except TypeError:
					print('{:<10s}'.format('None'), end='')
			print()

class Viterbi(HMM_NT_decoder):
	"""Standard Viterbi decoding in probabilty or log space."""

	def __init__(self, model=None, dna=None):
		"""
		Parameters
		----------
		+ model= `HMM`
		+ dna=   `DNA`
		"""

		self.model = model
		self.dna = dna
		self.null_score = None
		self.max_score = None
		self.score = None
		self.path = None
		self.matrix = None
		self.tmap = self._transition_map()
		self.smap = self._state_map()
		self._set_null_score()
		
		if len(self.dna.seq) > 100 and not self.model.log:
			sys.stderr.write('underflowing? use logspace model')
		
		# initialize viterbi matrix
		v = []
		for i in range(len(self.dna.seq)):
			v.append({})
			for s in self.model.states:
				v[i][s.name] = {'score' : None, 'trace' : None}

		# set initial probabilities
		for s in self.model.states:
			if self.model.log:
				v[0][s.name]['score'] = s.init + self.emission(s.name, 0)
			else:
				v[0][s.name]['score'] = s.init * self.emission(s.name, 0)

		# fill
		for i in range(1, len(self.dna.seq)):
			for here in self.tmap:
				ep = self.emission(here, i)
				max = -math.inf
				trace = None
				for prev in self.tmap[here]:
					tp = self.tmap[here][prev]
					pp = v[i-1][prev]['score']
					p = None
					if self.model.log: p = ep + tp + pp
					else:              p = ep * tp * pp
					if p > max:
						max = p
						trace = prev
				v[i][here]['score'] = max
				v[i][here]['trace'] = trace
		
		# set terminal probabilities
		for s in self.model.states:
			if self.model.log: v[-1][s.name]['score'] += s.term
			else:              v[-1][s.name]['score'] *= s.term

		# find maximum score and maximum ending state (sid)
		self.max_score = -math.inf
		sid = None
		for s in self.model.states:
			if v[-1][s.name]['score'] > self.max_score:
				self.max_score = v[-1][s.name]['score']
				sid = s.name

		# trace back
		path = []
		for i in range(len(self.dna.seq) -1, -1, -1):
			path.append(sid)
			sid = v[i][sid]['trace']
		path.reverse()

		self.score = None
		if self.model.log:
			self.score = self.max_score -self.null_score
		else:
			self.score = self.max_score / self.null_score
		self.path = path
		self.matrix = v

	def generate_path(self):
		"""Generates a path using the Viterbi algorithm. Returns a `Parse` object
		representing the optimal trace-back as calculated by the Viterbi algorithm."""

		return(Parse(path=self.path, score=self.score, decoder=self))

class StochasticViterbi(HMM_NT_decoder):
	"""Viterbi supporting multiple sub-optimal trace-backs."""

	def __init__(self, model=None, dna=None, seed=None):
		"""
		Parameters
		----------
		+ model= `HMM`
		+ dna=   `DNA`
		+ seed=  `int` random seed for reproducibility
		"""

		self.model = model
		self.dna = dna
		self.seed = seed
		self.null_score = None
		self.max_score = None
		self.score = None
		self.path = None
		self.matrix = None

		if self.seed: random.seed(self.seed)
		self.tmap = self._transition_map()
		self.smap = self._state_map()
		self._set_null_score()

		# initialize matrix
		v = []
		for i in range(len(self.dna.seq)):
			v.append({})
			for s in self.model.states:
				v[i][s.name] = {'score' : None, 'traces' : []}

		# set initial probabilities
		for s in self.model.states:
			if self.model.log:
				v[0][s.name]['score'] = s.init + self.emission(s.name, 0)
			else:
				v[0][s.name]['score'] = s.init * self.emission(s.name, 0)

		# fill
		for i in range(1, len(self.dna.seq)):
			for here in self.tmap:
				ep = self.emission(here, i)
				fwd = None
				if self.model.log:	fwd = -math.inf
				else: 			    fwd = 0

				traces = []
				for prev in self.tmap[here]:
					tp = self.tmap[here][prev]
					pp = v[i-1][prev]['score']
					p = None
					if self.model.log:
						p = ep + tp + pp
						fwd = toolbox.sumlog(fwd, p)
					else:
						p = ep * tp * pp
						fwd += p
					traces.append({'name' : prev, 'score' : p})
				v[i][here]['score'] = fwd

				if self.model.log:
					rsum = -math.inf
					for trace in traces:
						rsum = toolbox.sumlog(rsum, trace['score'])
						v[i][here]['traces'].append({
							'name': trace['name'],
							'cdf' : math.exp(rsum - fwd)
						})
				else:
					if fwd == 0: continue # is this the right thing to do?
					rsum = 0
					for trace in traces:
						rsum += trace['score'] / fwd
						v[i][here]['traces'].append({
							'name' : trace['name'],
							'cdf' : rsum
						})

		# set terminal probabilities
		for s in self.model.states:
			if self.model.log: v[-1][s.name]['score'] += s.term
			else:              v[-1][s.name]['score'] *= s.term

		self.matrix = v

	def generate_paths(self, n):
		"""Generates a list of paths using the Stochastic Viterbi algorithm. Returns
		a `list` of `Parse` objects."""

		term_cdf = []
		if self.model.log:
			sum = -math.inf
			for s in self.model.states:
				sum = toolbox.sumlog(sum, self.matrix[-1][s.name]['score'])
			rsum = -math.inf
			for s in self.model.states:
				rsum = toolbox.sumlog(rsum, self.matrix[-1][s.name]['score'])
				term_cdf.append({'name':s.name, 'cdf':math.exp(rsum-sum)})
		else:
			sum = 0
			for s in self.model.states:
				sum += self.matrix[-1][s.name]['score']

			rsum = 0
			for s in self.model.states:
				rsum += self.matrix[-1][s.name]['score'] / sum
				term_cdf.append({'name':s.name, 'cdf':rsum})

		# traceback n times, storing unique paths in table
		trace_table = {}
		for _ in range(n):
			# select terminal state
			rand = random.random()
			sid = None
			for final in term_cdf:
				if rand < final['cdf']:
					sid = final['name']
					break

			# traceback
			path = []
			for i in range(len(self.matrix) -1, -1, -1):
				path.append(sid)
				rand = random.random()
				for trace in self.matrix[i][sid]['traces']:
					if rand < trace['cdf']:
						sid = trace['name']
						break
			path.reverse()

			# calculate Viterbi score for path
			if self.model.log:
				score = self.smap[path[0]].init + self.emission(path[0], 0)
				for i in range(1, len(path)):
					ep = self.emission(path[i], i)
					tp = self.tmap[path[i]][path[i-1]]
					score += ep + tp
				score += self.smap[path[-1]].term
			else:
				score = self.smap[path[0]].init * self.emission(path[0], 0)
				for i in range(1, len(path)):
					ep = self.emission(path[i], i)
					tp = self.tmap[path[i]][path[i-1]]
					score *= ep * tp
				score *= self.smap[path[-1]].term

			# store path
			sig = '-'.join(path)
			if sig not in trace_table:
				trace_table[sig]= {
					'path':path,
					'score':score,
					'count':1
				}
			else:
				trace_table[sig]['count'] += 1

		# return parses sorted by frequency
		parses = []
		pn = 1
		for sig in sorted(trace_table,
				key=lambda x:trace_table[x]['count'], reverse=True):
			score = None
			if self.model.log:
				score = trace_table[sig]['score'] - self.null_score
			else:
				score = trace_table[sig]['score'] / self.null_score
			parses.append(Parse(
				score=score,
				path=trace_table[sig]['path'],
				freq=trace_table[sig]['count'] / n,
				decoder=self))
			pn += 1

		return parses

class Transcoder(HMM_NT_decoder):
	"""Class for scoring a labeled sequence"""
	
	def __init__(self, model=None, dna=None):
		"""
		Parameters
		----------
		+ model `HMM` usually converted to logspace
		+ dna   `DNA` 
		"""
		
		self.model = model
		self.dna = dna
		self.smap = self._state_map()
		self.tmap = self._transition_map()
		self.label_count = {}
		for state in model.states:
			stub = re.search(r'(\w+)', state.name)[1]
			if stub not in self.label_count:
				self.label_count[stub] = 0
			self.label_count[stub] += 1
		self._set_null_score()
	
		cds_len = None # for tracking phase
		path = []
		for f in self.dna.ftable.features:
			if self.label_count[f.type] == 1:
				for i in range(f.beg -1, f.end):
					path.append(f.type)
			elif f.type == 'CDS':
				if cds_len == None:
					if   f.phase == '.': cds_len = 0
					elif f.phase == 0:   cds_len = 0
					elif f.phase == 1:   cds_len = 2
					elif f.phase == 2:   cds_len = 1
				for i in range(f.length):
					frame = cds_len % 3
					path.append(f.type + '-' + str(frame))
					cds_len += 1
			else:
				if f.length != self.label_count[f.type]:
					raise DecodeError('state length mismatch')
				for i in range(f.length):
					path.append(f.type + '-' + str(i))
		
		self.max_score = self.score_path(path)
		self.score = None
		if self.model.log:
			self.score = self.max_score - self.null_score
		else:
			self.score = self.max_score / self.null_score

class ForwardBackward(HMM_NT_decoder):
	"""Computes posterior probabilities using Forward-Backward algorithm."""

	def __init__(self, model=None, dna=None):
		"""
		Parameters
		----------
		+ model= `HMM`
		+ dna=   `DNA`
		"""

		self.model = model
		self.dna = dna
		self.tmap = self._transition_map()
		self.smap = self._state_map()
		self.matrix = None
		prod = toolbox.prod
		sumlog = toolbox.sumlog

		p = {}
		T =	len(self.dna.seq) - 1

		# Initial probabilities
		for k in self.model.states:
			ki = k.init
			ep = self.emission(k.name, 0)
			p[k.name] = [{} for _ in range(0,T+1)]
			p[k.name][0] = {}
			p[k.name][T] = {}
			if self.model.log:
				p[k.name][0]['f'] = ki + ep
				p[k.name][T]['b'] = 999
			else:
				p[k.name][0]['f'] = ki * ep
				p[k.name][T]['b'] = 1
		# Build forward/backward matrix
		for t in range(1,T+1):
			for k in self.model.states:
				if self.model.log:
					# Forward probability
					p[k.name][t]['f'] = -999
					for j in self.model.states:
						prob = p[j.name][t-1]['f'] + \
							(self.tmap[k.name][j.name] if j.name in self.tmap[k.name] else -999) \
							+ self.emission(k.name, t)
						p[k.name][t]['f'] = sumlog(p[k.name][t]['f'], prob)
					# Backward probability
					p[k.name][T-t]['b'] = -999
					for j in k.next:
						prob = p[j][T-t+1]['b'] + self.tmap[j][k.name] \
							+ self.emission(j, T-t+1)
						p[k.name][T-t]['b'] = sumlog(p[k.name][T-t]['b'], prob)
				else:
					# Forward probability
					p[k.name][t]['f'] = sum(p[j.name][t-1]['f'] \
						* (self.tmap[k.name][j.name] if j.name in self.tmap[k.name] else 0) \
						for j in self.model.states) * self.emission(k.name, t)
					# Backward probability
					p[k.name][T-t]['b'] = sum(p[j][T-t+1]['b'] \
						* self.tmap[j][k.name] \
						* self.emission(j, T-t+1) for j in k.next)

		# Compute posterior probabilities
		for t in range(0,T+1):
			for k in self.model.states:
				if self.model.log:
					p[k.name][t]['posterior'] = p[k.name][t]['f'] + p[k.name][t]['b']
				else:
					p[k.name][t]['posterior'] = p[k.name][t]['f'] * p[k.name][t]['b']

		self.matrix = p

	def posterior(self, state_name, i):
		"""
		Returns the posterior probability of the given state at the specified position
		within a sequence.

		Parameters
		----------
		+ state_name `str` name of the state
		+ i          `int` position in sequence (1-based)
		"""

		return self.matrix[state_name][i]['posterior']
