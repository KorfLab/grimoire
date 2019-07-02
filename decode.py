import math
import sys
import json
from operator import itemgetter
import copy
import re
import random

import hmm
import genome
import toolbox

class DecodeError(Exception):
	pass

def inspect_matrix(decoder, beg, end, field):
	matrix = decoder.matrix
	
	# print numbers
	print('{:<6s}'.format(''), end='')
	for col in range(beg, end):
		print('{:<10d}'.format(col), end='')
	print()
	
	# print letters
	print('{:<6s}'.format(''), end='')
	for col in range(beg, end):
		nt = decoder.dna.seq[col:col+1]
		print('{:<10s}'.format(nt), end='')
	print()
	
	# print scores
	for state in decoder.model.states:
		print('{:<6s}'.format(state.name), end='')
		for i in range(beg, end):
			#if display == 'score':
			try:
				print('{:<10.3g}'.format(matrix[i][state.name][field]), end='')
			except TypeError:
				print('{:<10s}'.format('None'), end='')
		print()

class Parse:
	def __init__(self, score=None, path=None):
		self.score = score
		self.path = path
	
	def features(self, labels=None, chrom=None):
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
					features.append(genome.Feature(chrom, beg+1, end+1, '+',
						mypath[beg]))
					beg = end+1
					break
				elif i == len(mypath) -1:
					features.append(genome.Feature(chrom, beg+1, i+1, '+',
						mypath[beg]))
					end = len(mypath) +1
		return features

class HMM_NT_decoder:

	def state_map(self):
		map = {}
		for state in self.model.states + [self.model.null]:
			map[state.name] = state
		return map

	def transition_map(self):
		map = {}
		for state in self.model.states:
			map[state.name] = {}
		for state in self.model.states:
			for next in state.next:
				if next not in map:
					map[next] = {}
				map[next][state.name] = state.next[next]
		return map
	
	def set_null_score(self):
		self.null_score = 0 if self.log else 1
		for i in range(len(self.dna.seq)):
			ep = self.emission(self.model.null.name, i)
			if self.log: self.null_score += ep
			else:        self.null_score *= ep

	def emission(self, state_name, i):
		nt = self.dna.seq[i:i+1]
		state = self.smap[state_name]
	
		if i < 0: raise DecodeError('position less than 0')
		
		if self.log:
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

class Viterbi(HMM_NT_decoder):
	"""Standard Viterbi in probabilty or log space"""
	
	def __init__(self, model=None, dna=None, log=False):
		self.model = model
		self.dna = dna
		self.log = log
		self.null_score = None
		self.max_score = None
		self.score = None
		self.path = None
		self.matrix = None
		if self.log and self.model.logspace == False:
			self.model.convert2log()
		self.tmap = self.transition_map()
		self.smap = self.state_map()
		self.set_null_score()

		# initialize viterbi matrix
		v = []
		for i in range(len(self.dna.seq)):
			v.append({})
			for s in self.model.states:
				v[i][s.name] = {'score' : None, 'trace' : None}
			
		# set initial probabilities
		for s in self.model.states:
			v[0][s.name]['score'] = s.init * self.emission(s.name, 0)
	
		# fill
		for i in range(1, len(self.dna.seq)):
			for here in self.tmap:
				ep = self.emission(here, i)
				max = -math.inf
				trace = None
				for prev in self.tmap[here]:
					tp = self.tmap[prev][here]
					pp = v[i-1][prev]['score']
					p = None
					if self.log: p = ep + tp + pp
					else:        p = ep * tp * pp
					if p > max:
						max = p
						trace = prev
				v[i][here]['score'] = max
				v[i][here]['trace'] = trace
	
		# set terminal probabilities
		for s in self.model.states:
			if self.log: v[-1][s.name]['score'] += s.term
			else:        v[-1][s.name]['score'] *= s.term

		# find maximum score and maximum ending state (sid)
		self.max_score = -math.inf
		sid = None
		for s in self.model.states:
			if v[-1][s.name]['score'] > self.max_score:
				self.max_score = v[-1][s.name]['score']
				sid = s.name
	
		# trace back
		path = []
		for i in range(len(self.dna.seq) - 1, -1, -1):
			path.append(sid)
			sid = v[i][sid]['trace']
		path.reverse()
		
		self.score = (self.max_score - self.null_score if self.log
			else self.max_score / self.null_score)
		self.path = path
		self.matrix = v
		
	def generate_path(self):
		return(Parse(path=self.path, score=self.score))

class StochasticViterbi(HMM_NT_decoder):
	"""Viterbi supporting multiple sub-optimal trace-backs"""
		
	def __init__(self, model=None, dna=None, log=True):
		self.model = model
		self.dna = dna
		self.log = log
		self.null_score = None
		self.max_score = None
		self.score = None
		self.path = None
		self.matrix = None
		if self.log and self.model.logspace == False:
			self.model.convert2log()
		self.tmap = self.transition_map()
		self.smap = self.state_map()
		self.set_null_score()
		if self.log: raise DecodeError('StochasticViterbi does not support log space')
		# initialize matrix
		v = []
		for i in range(len(self.dna.seq)):
			v.append({})
			for s in self.model.states:
				v[i][s.name] = {'score' : None, 'traces' : []}
			
		# set initial probabilities
		for s in self.model.states:
			v[0][s.name]['score'] = s.init * self.emission(s.name, 0)
	
		# fill
		for i in range(1, len(self.dna.seq)):
			for here in self.tmap:
				ep = self.emission(here, i)
				fwd = None
				if self.log:	fwd = -math.inf
				else: 			fwd = 0
				traces = []
				for prev in self.tmap[here]:
					tp = self.tmap[prev][here]
					pp = v[i-1][prev]['score']
					p = None
					if self.log:
						p = ep + tp + pp
						fwd = toolbox.sumlog(fwd, p)
					else:
						p = ep * tp * pp
						fwd += p
					traces.append({'name' : prev, 'score' : p})
				v[i][here]['score'] = fwd
				
				rsum = 0
				sum = []
				for trace in traces:
					rsum += trace['score'] / fwd
					v[i][here]['traces'].append({
						'name' : trace['name'],
						'cdf' : rsum
					})
		# terminate matrix
		if self.log:
			raise DecodeError('log space not supported')
		else:
			for s in self.model.states: v[-1][s.name]['score'] *= s.term
		
		self.matrix = v
		
	def generate_paths(self, n):
		if self.log: raise DecodeError('log space yet not supported in StochasticViterbi')
		# choose terminal state by weighted random guessing
		# make terminal ramp
		sum = 0
		term = []
		for s in self.model.states:
			sum +=  self.matrix[-1][s.name]['score']

		rsum = 0
		for s in self.model.states:
			rsum += self.matrix[-1][s.name]['score'] / sum
			term.append({'name' : s.name, 'score' : rsum})
				
		# traceback n times, creating parse object for each traceback
		parses = []

		for i in range(n):
			# select terminal state
			rand = random.random()
			sid = None
			for s in term:
				if rand < s['score']:
					sid = s['name']
					break
			# traceback
			path = []
			for j in range(len(self.matrix) - 1, -1, -1):
				path.append(sid)
				rand = random.random()
				for trace in self.matrix[j][sid]['traces']:
					if rand < trace['cdf']:
						sid = trace['name']
						break

			path.reverse()
			
			# calculate Viterbi score for path
			sid = path[0]
			pid = None
			score = self.state_map()[sid].init * self.emission(sid, 0)
			for pos, sid in enumerate(path):
				tp = self.transition_map()[path[pos]][sid]
				ep = self.emission(sid, pos)
				score *= tp * ep
				pid = sid				

			parses.append(Parse(path=path, score=score))
		
		parse_count = {}
		for parse in parses:
			if ''.join(parse.path) not in parse_count:
				parse_count[''.join(parse.path)] = 0
			parse_count[''.join(parse.path)] += 1
		parse_freq = {}
		for parse in parse_count:
			parse_freq[parse] = parse_count[parse] / n
		print(json.dumps(parse_freq, indent=4))

		return(parses)			
				
		
		
# 		path = []
# 		for i in range(len(self.dna.seq), 0, -1):
# 			path.append(sid)
# 			sid = v[i][sid]['trace']
# 		path.reverse()
# 		
# 		self.score = (self.max_score - self.null_score if self.log
# 			else self.max_score / self.null_score)
# 		self.path = path
# 		self.matrix = v

class ViterbiXD(HMM_NT_decoder):
	"""Viterbi supporting states with explicit durations"""

class Posterior(HMM_NT_decoder):
	"""Forward-Backward decoding"""
	
	def __init__(self, model=None, dna=None, log=False):
		self.model = model
		self.dna = dna
		self.log = log
		self.null_score = None
		self.max_score = None
		self.score = None
		self.path = None
		self.matrix = None
		if self.log and self.model.logspace == False:
			self.model.convert2log()
		self.tmap = self.transition_map()
		self.smap = self.state_map()
		self.set_null_score()

		# initialize viterbi matrix
		v = []
		for i in range(len(self.dna.seq)):
			v.append({})
			for s in self.model.states:
				v[i][s.name] = {'fwd' : None, 'bwd' : None}
			
		# set initial probabilities
		for s in self.model.states: v[0][s.name]['fwd'] = s.init
	
		# fill
		for i in range(1, len(self.dna.seq)):
			for here in self.tmap:
				ep = self.emission(here, i)
				
				if self.log: fwd = -math.inf
				else:        fwd = 0
				
				for prev in self.tmap[here]:
					tp = self.tmap[prev][here]
					pp = v[i-1][prev]['fwd']
					p = None
					if self.log:
						p = ep + tp + pp
						fwd = toolbox.sumlog(fwd, p)
					else:
						p = ep * tp * pp
						fwd += p
				v[i][here]['fwd'] = fwd
	
		# set terminal probabilities
		for s in self.model.states:
			if self.log: v[-1][s.name]['fwd'] += s.term
			else:        v[-1][s.name]['fwd'] *= s.term

		self.matrix = v
		
	def generate_path(self):
		return(Parse(path=self.path, score=self.score))
