import math
import collections
import itertools

import numpy as np
import pandas as pd

def generate_kmers(alphabet,k):
    '''

    '''
    for kmer in itertools.product(alphabet,repeat=k):
        yield ''.join(kmer)

def iterkmers(seq,k):
    '''

    '''

    for i in range(len(seq)-(k-1)):
        yield seq[i:i+k]

def kmer_frequencies(seq,k,relative=False):
    '''

    '''
    freqs = dict(collections.Counter((kmer for kmer in iterkmers(seq,k))))
    if relative:
        num_kmers = len(seq)//k
        relative_freqs = {}
        for kmer,count in freqs.items():
            relative_freqs[kmer] = count/num_kmers
        freqs = relative_freqs
    return freqs


class HMM:
    def __init__(self,transitions=None,emissions=None,initials=None,states=None):
        
        self.states = states
        self.transitions = transitions
        self.emissions = emissions
        self.initial_probs = initials
        
        self.k = 2
        self.OWNZERO = 10**(-30)
    
    def _transition(self,labeled_sequence,states):
        
        num_states = len(states)
        
        kmers = kmer_frequencies(labeled_sequence,self.k)
        all_possible = generate_kmers(states,self.k)
        
        state_space = {state:i for i,state in enumerate(states)}
        transition_matrix = np.zeros(((num_states),(num_states)),dtype=np.float)
        
        for kmer in all_possible:
            if kmer not in kmers:
                kmers[kmer] = self.OWNZERO
                
        for kmer,value in kmers.items():
            prob = value/sum((value for key,value in kmers.items() if kmer[0]==key[0]))
            base1,base2 = kmer
            transition_matrix[state_space[base1],state_space[base2]] = prob
        return transition_matrix
    
    def _initial(self,labeled_sequence):
        total=len(list(filter(lambda base:base.islower(),labeled_sequence)))
        initial_matrix = {
            key:value/total 
            for key,value in kmer_frequencies(labeled_sequence,1).items()
                         } 
        return initial_matrix
    
    def _emission(self,states):
        emissions = {s:[1 if a.upper() == s.upper() else self.OWNZERO for a in 'ATCG'] for s in states}
        emission_matrix = pd.DataFrame(emissions,index=['A','T','C','G'])
        return  emission_matrix
    
    def from_sequence(self,labeled_sequence,states):
        self.labeled_sequence = labeled_sequence
        self.states = states
        self.transitions = self._transition(self.labeled_sequence,self.states)
        self.emissions = self._emission(self.states)
        self.initial_probs = self._initial(self.labeled_sequence)
        return self

    def viterbi(self,observations,transitions=None,emissions=None,initial_probs=None,states=None):
        
        if transitions is None:
            transitions = self.transitions
        if emissions is None:
            emissions = self.emissions.to_numpy()
        if initial_probs is None:
            initial_probs = self.initial_probs
        if states is None:
            states = self.states

        num_obs = len(observations)
        num_states = len(self.states)

        V = np.zeros((num_states,num_obs))
        trace_back = np.zeros((num_states,num_obs))

        transition = np.log(transitions)
        emission = np.log(emissions)

        #initialize
        for i,st in enumerate(self.states):
            e = emission[observations[0],i]
            V[i,0] = e + math.log(initial_probs[st])

        #Fill V matrix
        for i in range(1,num_obs):
            for l,state_l in enumerate(states):
                result = V[:,i-1]+transition[:,l]

                V[l,i] = result.max()
                trace_back[l,i] = result.argmax()
                V[l,i] += emission[observations[i],l]



        #traceback
        output = np.zeros(num_obs).astype(int)
        output[num_obs-1]=V[:,num_obs-1].argmax()
        for t in range(num_obs-1,0,-1):
            output[t-1]=trace_back[output[t],t]
        trace_back_result = ''.join(self.states[i] for i in output) 
        return trace_back_result


