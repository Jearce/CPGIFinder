import argparse 
import json
import sys

import numpy as np
import pandas as pd
from Bio import SeqIO


import SeqHandler
import HMM

DEFAULT_FTYPE='gb'
STATES = ['A','T','G','C','a','t','g','c']

def parse_args():
    '''

    '''

    description = '''Reads one or more FASTA files and predicts start and end
                     locations of CpG islands.'''

    parser = argparse.ArgumentParser(
    prog='CPGIFinder',
    usage='%(prog)s FASTA_FILE',
    description=description)

    parser.add_argument(
            'fasta_files',
            nargs='*',
            metavar='FASTA_FILE',
            type=str,
            help='Input FASTA file')
    parser.add_argument(
            '-train',
            nargs='*',
            metavar='TRAINING_FILE',
            type=str,
            default=DEFAULT_FTYPE,
            help='''Input training file. Genbank or FASTA.
            If FASTA then another file containing start and end locations is required.
            ''')
    parser.add_argument(
            '-ftype',
            metavar='FILE_TYPE',
            type=str,
            help='File type, fasta or gb. Default is gb.')

    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit()

    return  parser.parse_args()


def check(seq):
    '''

    '''
    pos = []
    count = 0
    for i in range(len(seq)):
        if seq[i].isupper() and seq[i-1].islower():
            start = i
        elif seq[i].isupper() and seq[i+1].islower():
            count +=1
            end = i
            
            print('{}.'.format(count),(start,end+1))


def main():
    args = parse_args()
    observation = SeqHandler.convert_seq(args.fasta_files[0]) #TODO:return list of seqs

    tp = np.loadtxt('probability_matricies/transition_probabilites.csv',delimiter=',')
    ep = pd.read_csv('probability_matricies/emission_probabilites.csv',index_col=0)
    with open('probability_matricies/initial_probabilites.json','r') as json_file:
        ip = json.load(json_file)


    cpgi_finder = HMM.HMM(transitions=tp,emissions=ep,initials=ip,states=STATES)
    result = cpgi_finder.viterbi(observation)
    check(result)


if __name__ == '__main__':
    main()
