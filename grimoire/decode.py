"""
Decode

This modeule contains classes and functions used for decoding, as well as
generating probabilistic paths and tables for the Viterbi algorithm and
variants.

Currently, there are 2 classes of the Viterbi algorithm:
	* Viterbi: finds path with maximal probability
	* StochasticViterbi: finds multi paths
"""

import math
import sys
import json
import copy
import re
import random

import grimoire.hmm as hmm
import grimoire.genome as genome
import grimoire.toolbox as toolbox

class DecodeError(Exception):
	pass

class Parse:
	"""Class representing a parse through a sequence"""
	
	def __init__(self, score=None, path=None, freq=None, decoder=None):
		"""
		Parameters
		----------
		score: float
			Usually a log-odds ratio vs. null model
		path: list
		 	State path (a list of state labels)
		freq: float
			Frequency of path occurring (should be set if stochastic)
		decoder: object
			Decoder that produced the path
		"""

		self.score = score
		self.path = path
		self.freq = freq
		self.decoder = decoder

	def features(self):
		"""Compile into a list of genome Feature objects"""

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
					features.append(genome.Feature(dna, beg+1, end+1, '+',
						mypath[beg]))
					beg = end+1
					break
				elif i == len(mypath) -1:
					features.append(genome.Feature(dna, beg+1, i+1, '+',
						mypath[beg]))
					end = len(mypath) +1
		return features

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

class HMM_NT_decoder:
	"""Base class for the HMM Nucleotide Decoders"""

	def state_map(self):
		"""
		A dictionary to look up states by name (used internally)
		key: state name (str), value: state (object)
		"""

		smap = {}
		for state in self.model.states + [self.model.null]:
			smap[state.name] = state
		return smap

	def transition_map(self):
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

	def set_null_score(self):
		"""Generate/set the score for the null model."""

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
		state_name: str
			Name of state
		i: int
			Position of state
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
		path: list
			List of state names
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
	"""Standard Viterbi in probabilty or log space"""

	def __init__(self, model=None, dna=None):
		"""
		Parameters
		----------
		model: object
			HMM object
		dna: object
			DNA object
		"""

		self.model = model
		self.dna = dna
		self.null_score = None
		self.max_score = None
		self.score = None
		self.path = None
		self.matrix = None
		self.tmap = self.transition_map()
		self.smap = self.state_map()
		self.set_null_score()
		
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
		"""Generates a path using the Viterbi algorithm"""

		return(Parse(path=self.path, score=self.score, decoder=self))

class StochasticViterbi(HMM_NT_decoder):
	"""Viterbi supporting multiple sub-optimal trace-backs"""

	def __init__(self, model=None, dna=None, seed=None):
		"""
		Parameters
		----------
		model: object
			HMM object
		dna: object
			DNA object
		seed: int
			Random seed for reproducibility
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
		self.tmap = self.transition_map()
		self.smap = self.state_map()
		self.set_null_score()

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
		"""Generates a list of paths using the Stochastic Viterbi algorithm"""

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
		model: HMM
			An HMM, probably converted to log space
		"""
		
		self.model = model
		self.dna = dna
		self.smap = self.state_map()
		self.tmap = self.transition_map()
		self.label_count = {}
		for state in model.states:
			stub = re.search('(\w+)', state.name)[1]
			if stub not in self.label_count:
				self.label_count[stub] = 0
			self.label_count[stub] += 1
		self.set_null_score()
	
	def score(self):
		cds_len = None # for tracking phase
		path = []
		for f in self.dna.features:
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
		return self.score_path(path)

class ForwardBackward(HMM_NT_decoder):
	"""Compute posterior probabilities using Forward-Backward algorithm"""

	def __init__(self, model=None, dna=None):
		"""
		Parameters
		----------
		model: object
			HMM object
		dna: object
			DNA object
		"""

		self.model = model
		self.dna = dna
		self.tmap = self.transition_map()
		self.smap = self.state_map()
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
		Get the probability of being in @state at @time

		Parameters
		----------
		state_name: str
			State name of interest
		i: int
			Position in the sequence
		"""

		return self.matrix[state_name][i]['posterior']
