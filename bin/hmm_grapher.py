from graphviz import Digraph
import grimoire.hmm as hmm
import argparse

## Command line w/ Argparse ##
parser = argparse.ArgumentParser(description='HMM Grapher')
parser.add_argument('--hmm', required=True, type=str,
	metavar='<path>', help='input HMM file name (%(type)s)')
parser.add_argument('--output', required=True, type=str,
	metavar='<path>', help='output graph name (%(type)s)')
parser.add_argument('--format', required=True, type=str,
	metavar='<path>', help='output graph format (%(type)s)')
arg = parser.parse_args()

'''Default Graph Attributes'''
graph = Digraph(format=arg.format)
graph.attr('node', shape='circle')

'''Read in File and Make Graph'''
hmm = hmm.HMM.read(arg.hmm)

for s1 in hmm.states:
	graph.node(s1.name)
	for s2 in s1.next:
		graph.edge(s1.name, s2, label=' '+str(s1.next[s2]))
graph.render(filename='test.gv', view=True)

'''Save Graph'''
with open(arg.output, 'w') as f:
    f.write(graph._repr_svg_())
