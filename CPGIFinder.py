import argparse 

from Bio import SeqIO

def parse_args():

    description = '''Reads one or more FASTA files and predicts start and end
    locations of CpG islands.'''

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
            'fasta_files',
            nargs='*',
            metavar='FASTA_File',
            type=str,
            help='Input FASTA file')
    parser.add_argument(
            '-train',
            nargs='*',
            metavar='TRAINING_FILE',
            type=str,
            help='''Input training file. Genbank or FASTA, if FASTA then another file containing start and end locations is required.
            ''')
    parser.add_argument(
            '-ftype',
            metavar='FILE_TYPE',
            type=str,
            help='File type, fasta or gb. Default is gb.')

    return  parser.parse_args()
