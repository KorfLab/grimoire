#!/usr/bin/python3
import sys

import toolbox.fasta
import toolbox.gff


def inspector(model, matrix, beg, end):
	print('{:<6s}'.format(''), end='')
	for col in range(beg, end):
		print('{:<10d}'.format(col), end='')
	print()
	for state in model.states:
		print('{:<6s}'.format(state.name), end='')
		for i in range(beg, end):
			print('{:<10.3g}'.format(matrix[i][state.name]['score']), end='')
		print()
			
def decode(model=None, seq=None, inspect=False):
	# transition matrix
	tm = {}
	for state in model.states:
		tm[state.name] = {}
	for s1 in model.states:
		for s2 in s1.next:
			tm[s2][s1.name] = s1.next[s2]
	# set up data structure
	viterbi = []
	for i in range(len(seq) + 1):
		viterbi.append({})
		for state in model.states:
			viterbi[i][state.name] = {'score' : None, 'trace' : None}
	# initialize matrix
	for state in model.states:
		viterbi[0][state.name]['score'] = state.init

	# fill in matrix
	for pos in range(1, len(seq) + 1):
		nt = seq[pos - 1 : pos]
		#print(pos)
		for this in model.states:
			max_score = 0
			max_state = None
			for prev in model.states:
				if this.name in tm[this.name]:
					tp = tm[this.name][prev.name]
					pp = viterbi[pos - 1][prev.name]['score']
					ep = None
					if this.ctxt == 0:
						ep = this.emit[nt] if nt in this.emit else 0.25
					else:
						if pos - 1 > this.ctxt:
							ctx = seq[pos - this.ctxt - 1: pos - 1]
							ep = this.emit[ctx][nt] if ctx in this.emit and nt in this.emit[ctx] else 0.25
						else:
							ep = 0.25
					score = tp * pp * ep
				if score > max_score:
					max_score = score
					max_state = prev.name
					#print(prev.name, this.name, max_score)
	
			viterbi[pos][this.name]['score'] = max_score
			viterbi[pos][this.name]['trace'] = max_state
	# terminate
	max_score = 0
	max_state = None
	for state in model.states:
		viterbi[-1][state.name]['score'] *= state.term
		if viterbi[-1][state.name]['score'] > max_score:
			max_score = viterbi[-1][state.name]['score']
			max_state = state.name
	# for verbose output
	if inspect:
		inspector(model, viterbi, 0, len(seq) + 1)	
	# traceback
	pos = len(seq)
	path = []
	prev = viterbi[pos][max_state]['trace']

	while pos > 0:
		path.append(prev)
		prev = viterbi[pos][max_state]['trace']
		max_state = prev
		pos -= 1
	path.reverse()
	return(path, max_score)
	
	
