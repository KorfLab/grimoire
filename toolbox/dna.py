"""DNA functions"""

def revcomp(seq) :
	comp = str.maketrans('ACGTRYMKWSBDHV', 'TGCAYRKMWSVHDB')
	seq = seq.translate(comp)[::-1]
	return seq

