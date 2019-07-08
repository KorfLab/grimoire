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

class Parse:
	def __init__(self, score=None, path=None, freq=None, pid=None):
		self.score = score
		self.path = path
		self.freq = freq
		self.pid = pid

	def features(self, labels=None, dna=None):
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
						mypath[beg], parent_id=self.pid))
					beg = end+1
					break
				elif i == len(mypath) -1:
					features.append(genome.Feature(dna, beg+1, i+1, '+',
						mypath[beg], parent_id=self.pid))
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
		self.null_score = 0 if self.model.log else 1
		for i in range(len(self.dna.seq)):
			ep = self.emission(self.model.null.name, i)
			if self.model.log: self.null_score += ep
			else:              self.null_score *= ep

	def emission(self, state_name, i):
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

	def inspect(self, beg, end, field):

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
		return(Parse(path=self.path, score=self.score))

class StochasticViterbi(HMM_NT_decoder):
	"""Viterbi supporting multiple sub-optimal trace-backs"""

	def __init__(self, model=None, dna=None, seed=None):
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
				sum +=  self.matrix[-1][s.name]['score']
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
					tp = self.tmap[path[i-1]][path[i]]
					score += ep + tp
				score += self.smap[path[-1]].term
			else:
				score = self.smap[path[0]].init * self.emission(path[0], 0)
				for i in range(1, len(path)):
					ep = self.emission(path[i], i)
					tp = self.tmap[path[i-1]][path[i]]
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
				pid='parse-' + str(pn)))
			pn += 1

		return parses

class ForwardBackward(HMM_NT_decoder):
	"""Forward Backward algorithm generating a matrix of posterior probabilities"""

	def __init__(self, model=None, dna=None):
		self.model = model
		self.dna = dna
		self.tmap = self.transition_map()
		self.smap = self.state_map()

		f, sum_f = self._forward()
		b, sum_b = self._backward()

		self.p = [] # matrix of posterior probabilities
		for i in range(len(self.dna.seq)):
			dist = {}
			for j in self.model.states:
				if self.model.log:
					dist[j.name] = toolbox.sumlog(f[i][j.name], b[i][j.name]) - sum_f
				else:
					dist[j.name] = (f[i][j.name] * b[i][j.name]) / sum_f
			self.p.append(dist)

	def get_prob(self, state, t):
		"""Returns the probability of state at time t"""
		return self.p[t][state]

	def find_best_state(self, t):
		"""Finds the most likely state at time t"""
		return max(self.p[t].keys(), key=(lambda key: self.p[t][key]))

	def _forward(self):
		f = []
		if self.model.log:
			prev = {k.name: toolbox.sumlog(self.emission(k.name,0), k.init) \
				for k in self.model.states}
		else:
			prev = {k.name: self.emission(k.name,0) * k.init for k in self.model.states}

		for i in range(len(self.dna.seq)):
			curr = {}
			for j in self.model.states:
				if self.model.log:
					if i == 0:
						curr[j.name] = toolbox.sumlog(self.emission(j.name, i), j.init)
					else:
						curr[j.name] = toolbox.sumlog(self.emission(j.name, i),
							toolbox.prod(toolbox.sumlog(self.tmap[k.name][j.name],prev[k.name]) \
							for k in self.model.states))
				else:
					if i == 0:
						curr[j.name] = self.emission(j.name, i) * j.init
					else:
						curr[j.name] = self.emission(j.name, i) \
							* sum(self.tmap[k.name][j.name] * prev[k.name] \
							for k in self.model.states)
			f.append(curr)
			prev = curr
		if self.model.log:
			sum_f = toolbox.prod(curr[j.name] + j.term for j in self.model.states)
		else:
			sum_f = sum(curr[j.name] * j.term for j in self.model.states)
		return f, sum_f

	def _backward(self):
		b = [{} for _ in range(len(self.dna.seq))]
		prev = {k.name: k.term for k in self.model.states}
		for i in range(len(self.dna.seq), 0, -1):
			curr = {}
			for j in self.model.states:
				if i == len(self.dna.seq):
					curr[j.name] = j.term
				elif self.model.log:
					curr[j.name] = toolbox.prod(self.tmap[j.name][k.name] \
						+ self.emission(k.name, i) + prev[k.name] \
						for k in self.model.states)
				else:
					curr[j.name] = sum(self.tmap[j.name][k.name] \
					* self.emission(k.name	, i) * prev[k.name] \
					for k in self.model.states)
			b[i-1] = curr
			prev = curr
			if self.model.log:
				sum_b = toolbox.prod(k.init + self.emission(k.name,0) \
					+ curr[k.name] for k in self.model.states)
			else:
				sum_b = sum(k.init * self.emission(k.name,0) * curr[k.name] \
					for k in self.model.states)
		return b, sum_b
