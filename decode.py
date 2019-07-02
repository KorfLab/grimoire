import math
import sys
import json
from operator import itemgetter
import copy
import re

import hmm
import genome
import toolbox

class DecodeError(Exception):
	pass

def inspect_matrix(decoder, matrix, beg, end):
	
	# print numbers
	print('{:<6s}'.format(''), end='')
	for col in range(beg, end):
		print('{:<10d}'.format(col), end='')
	print()
	
	# print letters
	print('{:<6s}'.format(''), end='')
	for col in range(beg, end):
		nt = None
		if col == 0: nt = '_'
		else: nt = decoder.dna.seq[col-1:col]
		print('{:<10s}'.format(nt), end='')
	print()
	
	# print scores
	for state in decoder.model.states:
		print('{:<6s}'.format(state.name), end='')
		for i in range(beg, end):
			#if display == 'score':
			try:
				print('{:<10.3g}'.format(matrix[i][state.name]['score']), end='')
			except TypeError:
				print('{:<10s}'.format('None'), end='')
			#if display == 'trace':
			#	try:
			#		print('{:<10s}'.format(matrix[i][state.name]['trace']), end='')
			#	except TypeError:
			#		print('{:<10s}'.format('None'), end='')
		print()

def initialize_matrix(model, seq):
	tm = get_transitions(model)
	print(json.dumps(tm, indent=4))
	matrix = []
	for i in range(len(seq)):
		matrix.append({})
		for state in model.states:
			matrix[i][state.name] = {'score' : None, 'trace' : None, 'traces' : {}}
	emit = seq[0:1]
	for this in model.states:
		max_score = 0
		max_state = None
		cumulative = 0
		traces = {}
		for prev in model.states:
			if prev.name not in tm[this.name]: continue
			pp = prev.init
			tp = tm[this.name][prev.name]
			ep = None
			if this.ctxt == 0 and emit in this.emit:
				ep = this.emit[emit]
			else:
				ep = 0.25
			score = pp * tp * ep
			traces[prev.name] = score
			cumulative += score
			if score > max_score:
				max_score = score
				max_state = prev.name
		matrix[0][this.name]['score'] = max_score
		matrix[0][this.name]['trace'] = max_state
		for trace in traces:
			try:
				traces[trace] /= cumulative
			except ZeroDivisionError:
				traces[trace] = 0
		matrix[0][this.name]['traces'] = traces
	inspect_matrix(model, matrix, 0, len(seq))
	return(matrix)

def fill_matrix(model, seq):
	matrix = initialize_matrix(model, seq)
	tm = get_transitions(model)
	for pos in range(1, len(seq)):
		emit = seq[pos - 1 : pos]
		for this in model.states:
			max_score = 0
			max_state = None
			for prev in model.states:
				if prev.name not in tm[this.name]: continue
				tp = tm[this.name][prev.name]
				pp = matrix[pos - 1][prev.name]['score']
				ep = None
				if this.ctxt == 0:
					ep = this.emit[emit] if emit in this.emit else 0.25
				else:
					if pos - 1 > this.ctxt:
						ctx = seq[pos - this.ctxt - 1 : pos - 1]
						ep = this.emit[ctx][emit] if ctx in this.emit and emit in this.emit[ctx] else 0.25
					else:
						ep = 0.25
				score = pp * tp * ep
				if score > max_score:
					max_score = score
					max_state = prev.name
			matrix[pos][this.name]['score'] = max_score
			matrix[pos][this.name]['trace'] = max_state
	# terminate
	for this in model.states:
		max_score = 0
		max_state = None
		for final in model.states:
			score = None
			if this.name in tm[final.name]:
				score = tm[final.name][this.name] * final.term
				if score > max_score:
					max_score = score
					max_state = final.name
			else:
				score = 0
		matrix[-1][this.name]['score'] *= max_score
	return(matrix)

def traceback(matrix):
	max_score = 0
	max_state = None
	# locate maximum score in final column
	pos = len(matrix) - 1
	for state in matrix[pos]:
		if matrix[pos][state]['score'] > max_score:
			max_score = matrix[pos][state]['score']
			max_state = matrix[pos][state]['trace']
		print(pos, state, matrix[pos][state]['score'], matrix[pos][state]['trace'])
	prev = matrix[pos][max_state]['trace']
	# begin traceback from maximum scoring state in final column
	path = []
	while pos > -1:
		path.append(prev)
		max_state = prev
		prev = matrix[pos][max_state]['trace']
		pos -= 1
	path.reverse()
	return(path, max_score)
		
def decode(model=None, seq=None, null_state=None, inspect=None):
	# transition matrix
	tm = get_transitions(model)
	# fill in matrix
	viterbi = fill_matrix(model, seq)
	inspect_matrix(model, viterbi, 0, len(seq) - 1)

	# for verbose output
	if inspect:
		inspector(model, viterbi, 0, len(seq), display=inspect)
	# traceback
	path, score = traceback(viterbi)
	return(path, score)
		
def fill_stochastic(model, seq):
	matrix = initialize_matrix(model, seq)
	tm = get_transitions(model)
	for pos in range(1, len(seq)):
		emit = seq[pos - 1 : pos]
		for this in model.states:
			max_score = 0
			max_state = None
			traces = {}
			cumulative = 0
			for prev in model.states:
				tp = tm[this.name][prev.name]
				pp = matrix[pos - 1][prev.name]['score']
				ep = None
				if this.ctxt == 0:
					ep = this.emit[emit] if emit in this.emit else 0.25
				else:
					if pos - 1 > this.ctxt:
						ctx = seq[pos - this.ctxt - 1 : pos - 1]
						ep = this.emit[ctx][emit] if ctx in this.emit and emit in this.emit[ctx] else 0.25
					else:
						ep = 0.25
				score = pp * tp * ep
				traces[prev.name] = score
				cumulative += score
				if score > max_score:
					max_score = score
					max_state = prev.name
			matrix[pos][this.name]['score'] = max_score
			matrix[pos][this.name]['trace'] = max_state
			for trace in traces:
				try:
					traces[trace] /= cumulative
				except ZeroDivisionError:
					traces[trace] = 0
			matrix[pos][this.name]['traces'] = traces
	
	# terminate. do not adjust tracer probabilities, only score
	for this in model.states:
		max_score = 0
		for final in model.states:
			score = tm[final.name][this.name] * final.term
			if score > max_score:
				max_score = score
		matrix[-1][this.name]['score'] *= max_score
	
	return(matrix)

def stochastic_traceback(matrix, n):
	paths = []
	for i in range(n):
		max_score = 0
		max_state = None
		pos = len(matrix) - 1
		# locate maximum score in final column
		for state in matrix[pos]:
			if matrix[pos][state]['score'] > max_score:
				max_score = matrix[pos][state]['score']
				max_state = matrix[pos][state]['trace']
		path = []
		prev = matrix[pos][max_state]['traces']
		
		while pos > -1:
			choices = sorted(prev, key=prev.get)
			selector = random.random()
			for choice in choices:
 				if selector < prev[choice]:
 					path.append(max_state)
 					max_state = choice
 					prev = matrix[pos][max_state]['traces']
 					pos -= 1
 					break
		path.append(max_state)
		path.reverse()
		paths.append(path)
	return(paths)

def stochastic(model=None, seq=None, null_state=None, n=1):
	# add multiple traces to the matrix
	tm = get_transitions(model)
	viterbi = fill_stochastic(model, seq)
	# traceback
	paths = stochastic_traceback(viterbi, n)
	# store emission and context of each state for re-scoring paths
	et = {}
	ct = {}
	for state in model.states:
		et[state.name] = state.emit
		ct[state.name] = state.ctxt
	# re-score paths THIS IS CLEARLY WRONG. PROBABILITY OF A PATH IS NOT STABLE
	corrected_scores = []
	corrected_paths = []
	for path in paths:
		corrected_score = viterbi[0][path[0]]['score']
		for pos in range(1, len(seq)):
			emit = seq[pos - 1 : pos]
			tp = tm[path[pos]][path[pos - 1]]
			ep = None
			if ct[path[pos]] == 0:
				ep = et[path[pos]][emit] if emit in et[path[pos]] else 0.25
			else:
				if pos - 1 > ct[path[pos]]:
					ctx = seq[pos - ct[path[pos]] - 1 : pos - 1]
					ep = et[path[pos]][ctx][emit] if ctx in et[path[pos]] and emit in et[path[pos]][ctx] else 0.25
				else:
					ep = 0.25
			corrected_score *= tp * ep
		corrected_scores.append(corrected_score)
		corrected_paths.append(path[1:len(path)])
	return(corrected_paths, corrected_scores)

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
		for i in range(len(self.dna.seq)+1):
			v.append({})
			for s in self.model.states:
				v[i][s.name] = {'score' : None, 'trace' : None}
			
		# set initial probabilities
		for s in self.model.states: v[0][s.name]['score'] = s.init
	
		# fill
		for i in range(1, len(self.dna.seq)+1):
			for here in self.tmap:
				ep = self.emission(here, i-1)
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
		for i in range(len(self.dna.seq), 0, -1):
			path.append(sid)
			sid = v[i][sid]['trace']
		path.reverse()
		
		self.score = self.max_score - self.null_score if self.log else self.max_score / self.null_score
		self.path = path
		self.matrix = v
		
	def generate_path(self):
		return(Parse(path=self.path, score=self.score))

class StochasticViterbi(HMM_NT_decoder):
	"""Viterbi supporting multiple sub-optimal trace-backs"""
		
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

		if self.log:
			raise DecodeError('stochastic decoding in log space not yet supported')
		# initialize viterbi matrix
		v = []
		for i in range(len(self.dna.seq)+1):
			v.append({})
			for s in self.model.states:
				v[i][s.name] = {'score' : None, 'traces' : []}
			
		# set initial probabilities
		for s in self.model.states: v[0][s.name]['score'] = s.init
	
		# fill
		for i in range(1, len(self.dna.seq)+1):
			for here in self.tmap:
				ep = self.emission(here, i-1)
				fwd = 0
				traces = []
				for prev in self.tmap[here]:
					tp = self.tmap[prev][here]
					pp = v[i-1][prev]['score']
					p = ep * tp * pp
					fwd += p
					traces.append({'name' : prev, 'score' : p})
				v[i][here]['score'] = fwd
				
				rsum = 0
				sum = []
				for trace in traces:
					rsum += trace['score'] / fwd
					v[i][here]['traces'].append({'name' : trace['name'], 'sum' : rsum})

		self.matrix = v
		
	def generate_paths(self, n):
		# choose terminal state by weighted random guessing
		# make terminal ramp
		term_ramp = []
		for s in self.model.states:
			s.term
		# traceback n times, creating parse object for each traceback
		
		pass

class ViterbiXD(HMM_NT_decoder):
	"""Viterbi supporting states with explicit durations"""
