from grimoire.sequence import DNA
from grimoire.feature import Feature, FeatureTable

dna = DNA(seq='ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT', name='test')

ft = FeatureTable(dna=dna)
ft.add_feature(Feature(dna, 1, 10, '+', 'one'))
ft.add_feature(Feature(dna, 15, 20, '+', 'two'))
ft.add_feature(Feature(dna, 7, 17, '+', 'three'))
ft.add_feature(Feature(dna, 5, 16, '+', 'four'))

f = ft.fetch(0, 25)

for e in f:
	print(e.beg, e.end, e.strand, e.type)
