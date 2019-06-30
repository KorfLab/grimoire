import math
import copy

import hmm
import toolbox




# def InputError(Exception):
# 	pass
# 	
# class Viterbi:
# 	def __init__(self, model):
# 		self.states = model.states
# 		self.transitions = {}
# 		for state in self.states:
# 			for next in state.next:
# 				if next not in self.transitions:
# 					self.transitions[next] = {}
# 				self.transitions[next][state.name] = state.next[next]
# 	def initialize(self, seq):
# 		matrix = []
# 		for i in range(len(seq)):
# 			matrix.append{}
# 			for state in self.states:
# 				matrix[i][state.name] = {'score' : None, 'trace' : None, 'traces' : {}}
# 		emit = seq[0:1]
# 		for this in self.states:
# 			max_score = 0
# 			max_state = None
# 			cumulative = 0
# 			traces = {}
# 			for prev in self.states:
# 				pp = prev.init
# 				tp = self.transitions[this.name][prev.name]
# 				ep = None
# 				if this.ctxt == 0 and emit in this.emit:
# 					ep = this.emit
# 				else:
# 					ep = 0.25
# 				score = pp * tp * ep
# 				traces[prev.name] = score
# 				cumulative += score
# 				if score > max_score:
# 					max_score = score
# 					max_state = prev.name
# 			matrix[0][this.name]['score'] = max_score
# 			matrix[0][this.name]['trace'] = max_state
# 			for trace in traces:
# 				try:
# 					traces[trace] /= cumulative
# 				except ZeroDivisionError:
# 					traces[trace] = 0
# 			matrix[0][this.name]['traces'] = traces
# 		self.matrix = matrix
# 	def compute(self, seq):
# 		
# 	def inspect(self, beg=0, end=len(self.matrix), display='score'):
# 		""""displays Viterbi matrix property 'score' or 'trace' from beg to end"""
# 		if display != 'score' and display != 'trace':
# 			raise InputError()
# 		print('{:<6s}'.format(''), end='')
# 		for i in range(beg, end):
# 			print('{:<10d}'.format(i), end='')
# 		for state in self.states:
# 			print('{:<6s}'.format(state.name), end='')
# 			for i in range(beg, end):
# 				if display == 'score':
# 					try:
# 						print('{:<10.3g}'.format(matrix[i][state.name]['score']), end='')
# 					except TypeError:
# 						print('{:<10s}'.format('None'), end='')
# 				else:
# 					try:
# 						print('{:<10s}'.format(matrix[i][state.name]['trace']), end='')
# 					except TypeError:
# 						print('{:<10s}'.format('None'), end='')
# 			print()
# 			
# 		
			
def safe_log(n):
	r = None
	try:
		r = math.log(n)
	except ValueError:
		r = -math.inf
	return(r)	

def inspect_matrix(model, matrix, beg, end, display='score'):
	print('{:<6s}'.format(''), end='')
	for col in range(beg, end):
		print('{:<10d}'.format(col), end='')
	print()
	for state in model.states:
		print('{:<6s}'.format(state.name), end='')
		for i in range(beg, end):
			if display == 'score':
				try:
					print('{:<10.3g}'.format(matrix[i][state.name]['score']), end='')
				except TypeError:
					print('{:<10s}'.format('None'), end='')
			if display == 'trace':
				try:
					print('{:<10s}'.format(matrix[i][state.name]['trace']), end='')
				except TypeError:
					print('{:<10s}'.format('None'), end='')
		print()

def log_space(model):
	hmm = hmm.HMM(name=model.name, states=model.states)
	for state in model.states:
		new_state = hmm.state.State(name=state.name, ctxt=state.ctxt,
									emit=hmm.state.emission_model(context=state.ctxt),
									init=safe_log(state.init),
									term=safe_log(state.init))
		
		if new_state.ctxt == 0:
			for nt in new_state.emit:
				new_state.emit[nt] = safe_log(state.emit[nt])
		else:
			for ctx in new_state.emit:
				for nt in new_state.emit[ctx]:
					new_state.emit[ctx][nt] = safe_log(state.emit[ctx][nt])
		for next in state.next:
			new_state.next[next] = safe_log(state.next[next])

def get_transitions(model):
	tm = {}
	for state in model.states:
		tm[state.name] = {}
	for state in model.states:
		for next in state.next:
			if next not in tm:
				tm[next] ={}
			tm[next][state.name] = state.next[next]
	return(tm)
		
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

def null_decode(model=None, seq=None, null_state=None):
	pass
		
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
	

def decodeP(model, fasta):
	print('decoding', model.name, fasta.seq)

def mylog(p):
	if p < 0: raise ValueError
	if p == 0: return -999
	else:      return math.log(p)

def convert2log(m1):
	m2 = copy.deepcopy(m1)
	for state in m2.states:
		state.init = mylog(state.init)
		state.term = mylog(state.term)
		if state.ctxt == 0:
			for nt in state.emit:
				state.emit[nt] = mylog(state.emit[nt])
		else:
			for ctx in state.emit:
				for nt in state.emit[ctx]:
					state.emit[ctx][nt] = mylog(state.emit[ctx][nt])
		for next in state.next:
			state.next[next] = mylog(state.next[next])
			
	return m2

hmm1 = hmm.HMM.read('toy.hmm')
hmm2 = convert2log(hmm1)

fasta = toolbox.FASTA_stream('toy.fasta')
for entry in fasta:
	decodeP(hmm1, entry)


