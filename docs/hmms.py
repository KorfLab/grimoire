#!/usr/bin/env python3

from os import system
import sys

models = ['acc', 'don', 'exon', 'intron', 'mRNA1', 'mRNA2']
forge = 'python3 bin/forge --fasta data/C.elegans.1percent.fasta.gz --gff3 data/C.elegans.1percent.gff3.gz'
graph = 'python3 bin/hmm_grapher'
fp = open('docs/hmms.html', 'w+')
fp.write('<html><head><title>HMMs</title><body><h1>HMMs</h1>')
for model in models:
	hmm = model + '.hmm'
	svg = model + '.svg'
	png = model + '.png'
	system(forge + ' --hmm ' + hmm + ' --model ' + model)
	system(graph + ' --hmm ' + hmm + ' --svg ' + svg)
	system('convert ' + svg + ' ' + png) # requires ImageMagick
	system('mv ' + png + ' docs')
	system('mv ' + svg + ' docs')
	system('rm ' + hmm)
	fp.write('<img src="' + svg + '"><hr>')
fp.write('</body></html>')