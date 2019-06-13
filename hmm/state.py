#!/usr/bin/python3

import argparse
import sys
import re
import json
import operator

import kmer

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

def train_emissions(seqs, context=0) :
	counts = []
	freqs = []
	for i in range(len(seqs[0])):
		counts.append(emission_model(context=context))
	
	if (context == 0) :
		for seq in seqs:	
			for i in range(len(seq)) :
				nt = seq[i:i+1]
				if (nt in counts[i]) : counts[i][nt] += 1
		for count in counts:
			total = 0
			freq = {}
			for nt in count: total += count[nt]
			for nt in count: freq[nt] = count[nt] / total	
			freqs.append(freq)
	else :
		for i in range(context):
			for ctx in counts[i]:
				for nt in counts[i][ctx]:
					counts[i][ctx][nt] = 1
		for seq in seqs:
			for i in range(len(seq) - context):
				pos = i + context
				ctx = seq[i:i+context]
				nt = seq[i + context : i + context + 1]
				if (ctx in counts[pos] and nt in counts[pos][ctx]) :
					counts[pos][ctx][nt] += 1
		for i in range(len(counts)):
			freqs.append({})
			for ctx in counts[i]:
				total = 0
				freqs[i][ctx] = {}
				for nt in counts[i][ctx]: total += counts[i][ctx][nt]
				for nt in counts[i][ctx]:
					if total != 0:
						freqs[i][ctx][nt] = counts[i][ctx][nt] / total
					else:
						freqs[i][ctx][nt] = 0
	return freqs

def train_emission(seqs, context=0):
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
		for nt in count: freq[nt] = count[nt] / total
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
			for nt in count[ctx]: freq[ctx][nt] = count[ctx][nt] / total
				
	return freq

class State:
	def __init__(self, name=None, context=None, emits=None, init=0, term=0, next={}):
		self.name = name
		self.ctxt = context
		self.init = init
		self.term = term
		self.emit = emits
		self.next = next
		
	@classmethod
	def from_json(cls, json_string):
		state = cls()
		state.__dict__ = json.loads(json_string)
		
	def to_json(self):
		return(json.dumps(self.__dict__, indent = 4))

def state_factory(stub, emissions):
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
		state = State(name=stub + '-' + str(i),
			emits=emissions[i],
			context=context)
		state_list.append(state)
		# connect them?
	return(state_list)

