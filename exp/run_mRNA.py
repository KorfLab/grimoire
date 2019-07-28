import re
from multiprocessing import Pool

import grimoire.genome as genome
import grimoire.hmm as hmm
import grimoire.decode as decode
import grimoire.sequence as sequence

Model = None
Tcode = None

def get_predicted(dna, model):
	v = decode.Viterbi(dna=dna, model=model)
	p = v.generate_path()
	beg = None
	end = None
	for f in p.features():
		if f.type == 'ATG':
			beg = f.beg
		elif f.type == 'STOP':
			end = f.end
	cds = dna.seq[beg-1:end]
	return beg, end, v.score

def analyze(dna):

	# get the position of the 'true' start & stop
	b0 = None
	e0 = None
	for f in dna.features:
		if f.type == 'ATG':
			b0 = f.beg
		elif f.type == 'STOP':
			e0 = f.end
	t0 = decode.Transcoder(model=Model, dna=dna)
	s0 = t0.score

	# get the position of the predicted start and stop
	b1, e1, s1 = get_predicted(dna, Model)
	
	# create a difference flag
	flag = 'same'
	if b0 != b1 and e0 == e1:
		flag = 'start'
	elif e0 != e1:
		flag = 'stop'
	
	return ','.join([dna.name, flag, str(b0), str(e0), str(b1), str(e1),
		str(s0), str(s1)])


if __name__ == '__main__':
	Model = hmm.HMM.read('mRNA2.hmm')
	Model.convert2log()
	gen = genome.Genomic(fasta='test.fasta', gff='test.gff')
	print(','.join(['gene', 'diff', 'tbeg', 'tend', 'pbeg', 'pend',
		'tscore', 'pscore']))

#	for dna in gen:
#		analyze(dna)
	
	with Pool(6) as proc:
		outs = proc.map(analyze, gen)
		for out in outs:
			print(out)
			