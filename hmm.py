
import json

import sequence
import toolbox

class HMMError(Exception):
	pass

def emission_model(context=1, alphabet='nt'):
	letters = None
	if alphabet == 'nt': letters = sequence.DNA.canonical
	elif alphabet == 'aa': letters = sequence.Protein.canonical
	else: raise HMMError('alphabet error: ', alphabet)
	
	if (context > 0):
		table = sequence.generate_kmers(alphabet=alphabet, k=context)
		for k in table:
			table[k] = {}
			for a in letters:
				table[k][a] = 0
	elif (context == 0) :
		table = sequence.generate_kmers(alphabet=alphabet, k=1, pseudo=0)
	else :
		raise HMMError('negative context')
	
	return table

def train_emissions(seqs, context=0):
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
						freqs[i][ctx][nt] = round(counts[i][ctx][nt]/total, 4)
					else:
						freqs[i][ctx][nt] = 0
	return freqs

def train_cds(seqs, context=0):
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
			#print(seq[0:10])
			for i in range(0, len(seq) -context -3, 3):
				for j in range(3):
					pos = i + j + context
					ctx = seq[i+j:i+j+context]
					nt = seq[i + j + context]
					#print('>>', i, j, pos, ctx, nt)
					if (ctx in counts[j] and nt in counts[j][ctx]) :
						counts[j][ctx][nt] += 1
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

class State:
	"""Class for HMM States"""
	
	def __init__(self, name=None, context=None, emits=None, init=0, term=0, next=None):
		self.name = name
		self.ctxt = context
		self.init = init
		self.term = term
		self.emit = emits
		self.next = {}
		
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
	return(state_list)


def connect_all (states) :
	for i in range(len(states)-1) :
		states[i].next[states[i+1].name] = 1

def connect2 (s1, s2, p) :
	s1.next[s2.name] = p

class HMMdecoder(json.JSONEncoder):
	def default(self, o):
		return o.__dict__

class HMM:
	"""Class for HMMs"""
	
	def __init__(self, name=None, logspace=False, states=None, null=None):
		self.name = name
		self.logspace = logspace
		self.states = states
		self.null = null
	
	@classmethod
	def read(cls, filename):
		fp = open(filename, 'r')
		d = json.loads(fp.read())
		hmm = HMM()
		hmm.name = d['name']
		hmm.logspace = d['logspace']
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
		fp = open(filename, 'w+')
		fp.write(json.dumps(self.__dict__, indent=4, cls=HMMdecoder))

	def convert2log(self):
		if self.logspace:
			raise HMMError('attempt to re-convert to logs, make a copy first')
		
		for state in self.states + [self.null]:
			state.init = toolbox.mylog(state.init)
			state.term = toolbox.mylog(state.term)
			if state.ctxt == 0:
				for nt in state.emit:
					state.emit[nt] = toolbox.mylog(state.emit[nt])
			else:
				for ctx in state.emit:
					for nt in state.emit[ctx]:
						state.emit[ctx][nt] = toolbox.mylog(state.emit[ctx][nt])
			for next in state.next:
				state.next[next] = toolbox.mylog(state.next[next])
		self.logspace = True