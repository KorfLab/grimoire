#!/usr/bin/env python3

import sys
import json

import grimoire.toolbox as toolbox
import grimoire.hmm as hmm
import grimoire.decode as decode
from grimoire.sequence import DNA
from grimoire.genome import Genome, Feature

model = hmm.HMM.read('mRNA.hmm')
dna = DNA(seq='ACATGCCCTT', name='test')
v = decode.Viterbi(model=model, dna=dna)
v._inspect('score')
p = v.generate_path()
features = p.features()
for f in features:
	print(f)
