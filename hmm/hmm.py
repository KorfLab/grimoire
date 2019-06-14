import hmm.state

def connect_all (states) :
	for i in range(len(states)-1) :
		states[i].next[states[i+1].name] = 1

def connect2 (s1, s2, p) :
	s1.next[s2.name] = p

class Hmm:
	"""something"""
	
	def __init__(self, name, states):
		self.name = name
		self.state = states
	
	def dump(self):
		print('hello')

def create_exon_hmm(
		exon_seqs=None, exon_context=None,
		acc_seqs=None, acc_context=None,
		don_seqs=None, don_context=None,
		int_seqs=None, int_context=None):
	
	acc_emits = hmm.state.train_emissions(acc_seqs, context=acc_context)
	don_emits = hmm.state.train_emissions(don_seqs, context=don_context)
	exon_emits = hmm.state.train_emission(exon_seqs, context=exon_context)
	int_emits = hmm.state.train_emission(int_seqs, context=int_context)
	
	acc_states = hmm.state.state_factory('ACC', acc_emits)
	don_states = hmm.state.state_factory('DON', don_emits)
	exon_state = hmm.state.State(name='EXON', context=exon_context, emits=exon_emits)
	int_state = hmm.state.State(name='INT', context=int_context, emits=int_emits)

	connect_all(acc_states)
	connect2(acc_states[-1], exon_state, 1)
	connect2(exon_state, exon_state, 0.99)
	connect2(exon_state, don_states[0], 0.01)
	connect_all(don_states)
	connect2(don_states[-1], int_state, 1)
	connect2(int_state, int_state, 0.99)
	connect2(int_state, acc_states[0], 0.01)
	
	return Hmm('simple_exon', [[int_state] + acc_states + [exon_state] + don_states])
	