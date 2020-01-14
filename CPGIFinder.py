import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-in',nargs='*',metavar='FASTA_FILE',type=str,help='Input FASTA file')
args = parser.parse_args()
