'''
Used to handle any sequence data.
'''
import sys

from Bio import SeqIO



def convert_seq(fastafile):
    '''Reads FASTA file and encode sequence for viterbi algorithm.

    Arguements:
     fastafile: FASTA file containing query sequence

    Results:
     Tuple with FASTA ID and encoded sequence
    '''
    encoding = {'A':0,'T':1,'C':2,'G':3}
    record = SeqIO.read(fastafile,'fasta')
    return record.id,[encoding[base] for base in record.seq]



def writeout(results):
    '''Writes start and end locations of predicited CpG islands.

    Arguments:
     results: list of tuples containing (id,result)

    Result:
     Tab-delimited output of sequence id, start, and end locations
    '''

    sys.stdout.write('{}\t{}\t{}\n'.format('#SequenceID','Start','End'))
    for id,seq in results:

        length = len(seq)
        #CpG islands are labeled with upper case chars
        for i in range(length):
            if (i == 0 and seq[i].isupper()) or (seq[i].isupper() and seq[i-1].islower()):
                start = i
            elif (i+1 == length and seq[i].isupper()) or (seq[i].isupper() and seq[i+1].islower()):
                end = i

                sys.stdout.write('{}\t{}\t{}\n'.format(id,start+1,end+1))
