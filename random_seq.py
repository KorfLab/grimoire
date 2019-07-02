"""generate random sequences with given hmm

users need to inform the program
how many time could generate function run
to get desire sequence with given end state

when having empty emit, e.g.:
{
    "A": 0,
    "C": 0,
    "G": 0,
    "T": 0
}
random element in alphabet will be added

"""

import json
import random

class Random_Seq:
    def __init__(self, model, alphabet, length, starter="none", ender="none"):
        self.model = model
        if self.model["logspace"] == True:
            print("Error: logspace == True")
        self.alphabet = alphabet
        self.length = length
        self.start = {}
        self.end = {}
        for n in range(len(self.model["states"])):
            if self.model["states"][n]["init"] > 0:
                v = self.model["states"][n]["init"]
                self.start[self.model["states"][n]["name"]] = v
            if self.model["states"][n]["term"] > 0:
                v = self.model["states"][n]["term"]
                self.end[self.model["states"][n]["name"]] = v
        self.starter = starter
        self.ender = ender

    def _rand(self):
        x = random.randint(0, len(self.alphabet) - 1)
        return self.alphabet[x]

    def _find(self, name):
        for n in range(len(self.model["states"])):
            if self.model["states"][n]["name"] == name:
                return n

    def _get(self,dict):
        sum = 0
        x = random.uniform(0.00, 1.00)
        for s in dict:
            sum += dict[s]
            if sum >= x:
                return s
        raise Exception

    def _check(self,state_seq, cur_name, ctxt):
        for i in range(-(ctxt),0,1):
            if state_seq[i] != cur_name:
                return False
        return True

    def generate_null(self):
        seq = []
        while len(seq) < self.length:
            seq.append(self._get(self.model["null"]["emit"]))
        seq = ''.join(seq)
        return seq


    def generate(self, free_end = True, trial = 1):
        stat = []
        seq = []
        if self.starter != "none":
            start = self.starter
        else:
            start = self._get(self.start)
        cur = self._find(start)
        cur_name = start
        ctxt = self.model["states"][cur]["ctxt"]
        stat.append(cur_name)
        if ctxt == 0:
            seq.append(self._get(self.model[cur]["emit"]))
        else:
            seq.append(self._rand())

        while len(stat) < self.length:
            st = self._get(self.model["states"][cur]["next"])
            if st != cur_name:
                cur =  self._find(st)
                ctxt = self.model["states"][cur]["ctxt"]
            enough = self._check(stat, st, ctxt)
            stat.append(st)
            if ctxt == 0:
                seq.append(self._get(self.model["states"][cur]["emit"]))
            elif enough:
                temp_pre =[]
                for i in range(-(ctxt),0,1):
                    temp_pre.append(seq[i])
                prefix = ''.join(temp_pre)
                try:
                    seq.append(self._get(self.model["states"][cur]["emit"][prefix]))
                except:
                    seq.append(self._rand())
            else:
                seq.append(self._rand())
        seq = ''.join(seq)
        if free_end == True:
            return stat, seq
        else:
            if self.ender != "none":
                end = self.ender
            else:
                end = self._get(self.end)
            if stat[-1] != end:
                trial -= 1
                if trial == 0:
                    print("Error: More Trilas?")
                    return [], ''
                else:
                    stat, seq = self.generate(self, free_end, trial)
            return stat, seq


"""example code
"""
ALPHABET = {
    'nt' : ['A', 'C', 'G', 'T'],
    'aa' : ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
        'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'],
}

f = open('toy.hmm')
m =  json.load(f)
a = Random_Seq(model = m, alphabet = ALPHABET['nt'], length = 1000)
q = a.generate_null()
print("null:")
print(q)
q,b = a.generate()
print("free end:")
print(q,b)
q,b = a.generate(free_end = False, trial = 10)
print("no free end:")
print(q,b)
