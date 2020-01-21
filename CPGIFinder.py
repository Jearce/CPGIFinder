#!/usr/bin/env python
import argparse 
import json
import sys
import os

import numpy as np
import pandas as pd

import SeqHandler
import HMM

program_dir = os.path.dirname(os.path.realpath(__file__))

TRANSITION_MATRIX = program_dir+'/probability_matricies/transition_probabilites.csv'
EMISSION_MATRIX = program_dir+'/probability_matricies/emission_probabilites.csv'
INITIAL_MATRIX = program_dir+'/probability_matricies/initial_probabilites.json'
STATES = ['A','T','G','C','a','t','g','c']

def parse_args():
    '''Parse command line options.

    Returns options object with command line values as attributes.
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

    if len(sys.argv)==1:
        parser.print_usage(sys.stderr)
        sys.exit()

    return  parser.parse_args()

def main():
    args = parse_args()

    if args.fasta_files:

        #load matricies 
        tp = np.loadtxt(TRANSITION_MATRIX,delimiter=',')
        ep = pd.read_csv(EMISSION_MATRIX,index_col=0)
        with open(INITIAL_MATRIX,'r') as json_file:
            ip = json.load(json_file)

        #initialize HMM
        cpgi_finder = HMM.HMM(transitions=tp,emissions=ep,initials=ip,states=STATES)

        #predict islands
        results = []
        for fasta_file in args.fasta_files:

            try:
                seq_id,observation = SeqHandler.convert_seq(fasta_file) 
            except ValueError as err:
                print('Incorrect file format:',err)
            else:
                result = cpgi_finder.viterbi(observation)
                results.append((seq_id,result))

        if results:
            SeqHandler.writeout(results)

if __name__ == '__main__':
    main()

