import argparse 

from Bio import SeqIO

DEFAULT_FTYPE='gb'

def parse_args():

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

    return  parser.parse_args()

def main():
    args = parse_args()

if __name__ == '__main__':
    main()
