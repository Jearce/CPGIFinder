from Bio import SeqIO

def convert_seq(fastafile):
    print(fastafile)
    '''

    '''
    observations = {'A':0,'T':1,'C':2,'G':3}
    record = SeqIO.read(fastafile,'fasta')
    return [observations[base] for base in record.seq]

