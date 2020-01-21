# CPGIFinder
CPGIFinder predicts CpG islands in a given genomic sequence using Hidden MarvokModel as described by [Durbin et al](https://www.semanticscholar.org/paper/Biological-Sequence-Analysis%3A-Probabilistic-Models-Durbin-Eddy/571f5bbecd3a083a2bb6844f59a3f8cea237252e)
## Prerequisites
Python 3 is required.

Packages: 
* numpy
* pandas
* biopython

## Installing

Cloning this repository.
> $ git clone https://github.com/Jearce/CPGIFinder.git

Move into respository directory.
> $ cd CPGIFinder

## Behaviour
CPGIFinder accepts one or more FASTA files on the command line. For each FASTA file CpG islands are predicted using viterbi algorithm. The FASTA ID's,start, and end locations of the predicted CpG islands are written to `stdout` in tab-delimited format.

Before running CPGIFinder first make executable:
> $ chmod +x CPGIFinder.py

## Help
To display usage information on the command line use -h or --help argument:

```
$ CPGIFinder.py -h

usage: CPGIFinder.py FASTA_FILE

Reads one or more FASTA files and predicts start and end locations of CpG
islands.

positional arguments:
  FASTA_FILE  Input FASTA file

optional arguments:
  -h, --help  show this help message and exit
```

## Reading FASTA files on command line
There are two sequences in the test_data directory that will be used in the examples below.

The example below illustrates using CPGIFinder on one FASTA file:

```
$ CPGIFinder.py test_data/human.fasta 
#SequenceID	Start	End
AF093117.1	29347	29477
AF093117.1	29631	30366
AF093117.1	46772	46902
AF093117.1	48165	48348
AF093117.1	49096	49350
AF093117.1	67899	68205
AF093117.1	78824	80327
AF093117.1	93317	94103
AF093117.1	94201	95120
```

The example below illustrates using CPGIFinder on two FASTA files.

```
$ CPGIFinder.py test_data/human.fasta test_data/zebrafish.fasta 
#SequenceID	Start	End
AF093117.1	29347	29477
AF093117.1	29631	30366
AF093117.1	46772	46902
AF093117.1	48165	48348
AF093117.1	49096	49350
AF093117.1	67899	68205
AF093117.1	78824	80327
AF093117.1	93317	94103
AF093117.1	94201	95120
AL672171.8	1457	1607
AL672171.8	50521	50737
AL672171.8	51058	51167
AL672171.8	68463	68618
```
## Authors
* Jessie Arce

## License
CPGIFinder is licensed under the GNU General Public License v3.0

## Acknowledgments
Hidden Markov Model implementation was based on [mchmm](https://github.com/maximtrp/mchmm) and [cpg-island-prediction-HMM](https://github.com/devanshdalal/cpg-island-prediction-HMM). Reading [Durbin et al](https://www.semanticscholar.org/paper/Biological-Sequence-Analysis%3A-Probabilistic-Models-Durbin-Eddy/571f5bbecd3a083a2bb6844f59a3f8cea237252e) inspired me to create this program and test what I understood. 



