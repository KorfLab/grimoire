"""
Decode

This modeule contains classes and functions used for decoding, as well as
generating probabilistic paths and tables for the Viterbi algorithm and
variants.

Currently, there are 2 classes of the Viterbi algorithm:
	* Viterbi: finds path with maximal probability
	* StochasticViterbi: finds multi paths
"""



import grimoire.hmm as hmm
import grimoire.genome as genome
import grimoire.toolbox as toolbox
