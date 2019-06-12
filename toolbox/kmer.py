"""kmer module for generating biological kmers
"""

ALPHABET = {
	'nt' : ['A', 'C', 'G', 'T'],
	'aa' : ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
		'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'],
}

def _kmers(alphabet, table, key, n, k, v):
	if (k == 0) :
		if key not in table:
			table[key] = v
			return

	for i in range(n):
		t = key + alphabet[i]
		_kmers(alphabet, table, t, n, k - 1, v)

def generate(alphabet='nt', k=1, pseudo=0):
	"""Creates dictionary of all kmers of either nt or aa alphabet."""
	table = {}
	if (alphabet == 'nt') :
		_kmers(ALPHABET[alphabet], table, '', 4, k, pseudo)
	elif (alphabet == 'aa') :
		_kmers(ALPHABET[alphabet], table, '', 20, k, pseudo)
	return table

