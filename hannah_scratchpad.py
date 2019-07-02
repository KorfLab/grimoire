import sys
import json

import toolbox
from toolbox import GFF_file as GFF
from toolbox import FASTA_file as FASTA
from hmm import HMM

gf = GFF('internal.gff')
ff = FASTA('data/TAIR10_1.fasta')

seqs = []

for id in ff.ids:
	seq = ff.get(id).seq
	for feature in gf.get(chrom='Chr1', type='model'):
		if feature.strand == '-': continue
		seqs.append(seq[feature.beg : feature.end])

hmm = HMM.read('internal.hmm')

for seq in seqs:
	print(seq)
	path, score = hmm.decode(model=hmm, seq=seq)
	print(path)
	print(score)
# 	paths, scores = viterbi.stochastic(model=hmm, seq=seq, n=100)
# 	for i in range(len(paths)):
# 		print(paths[i])
# 		print(scores[i])


#############################
# Stuff below was in decode #
#############################

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
