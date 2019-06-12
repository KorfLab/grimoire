#!/usr/bin/python3

import argparse
import sys
import re
import json
import operator

import kmer

def revcomp(seq) :
	comp = str.maketrans('ACGTRYMKWSBDHV', 'TGCAYRKMWSVHDB')
	seq = seq.translate(comp)[::-1]
	return seq

def emission_model(context=1, alphabet='nt') :
	if (context > 0) :
		table = kmer.generate(alphabet=alphabet, k=context)
		for k in table :
			table[k] = {}
			for a in kmer.ALPHABET[alphabet] :
				table[k][a] = 0
	elif (context == 0) :
		table = kmer.generate(alphabet=alphabet, k=1, pseudo=0)
	else :
		sys.exit("negative context not allowed")
	
	return table

class State:
	def __init__(self, name=None, context=None, init=None, term=None, next={}):
		self.name = name
		self.ctxt = context
		self.init = init
		self.term = term
		self.emit = emission_model(context=context) if context is not None else {}
		self.next = next if next else {}
		
	@classmethod
	def from_json(cls, json_string):
		state = cls()
		state.__dict__ = json.loads(json_string)
		
	def to_json(self):
		return(json.dumps(self.__dict__, indent = 4))

def state_factory(stub, state_count=1, context=0):
	state_list = []
	for i in range(state_count) :
		state = State(name=stub + '-' + str(i), context=context)
		if (i < state_count - 1) :
			next = stub + '-' + str(i+1)
			state.next[next] = 1
		#	state.term = 0
		#if (i > 0) :
		#	state.init = 0
		state_list.append(state)
	return(state_list)

example_empty = State()

print(example_empty.to_json())

example_state_array = state_factory('example', state_count=2, context=1)

for state in example_state_array :
	print(state.to_json())
