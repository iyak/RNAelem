# usage: python _run-bertrbp-h0.py --positive <p.fa> --negative <n.fa> > train.tsv
#        python _run-bertrbp-h0.py --no-info  <t.fa>                   > test.tsv


from argparse import ArgumentParser
from pathlib import Path
import re
import sys


def parse_fasta(fasta_path, is_positive):
    with open(fasta_path) as fi:
        label = 0
        while True:
            head = fi.readline().strip(" \n")
            seq  = fi.readline().strip(" \n")
            if (not head) or (not seq):
                break
            seq = re.sub("U", "T", seq)
            kmer = [seq[x: x + 3] for x in range(len(seq) + 1 - 3)]
            kmers = " ".join(kmer)
            if   True  == is_positive: label = 1
            elif False == is_positive: label = 0
            else                     : label = 1 - label
            sys.stdout.write("%s\t%d\n" % (kmers, label))


if "__main__" == __name__:
    ap = ArgumentParser()
    ap.add_argument("--positive", metavar="FASTA")
    ap.add_argument("--negative", metavar="FASTA")
    ap.add_argument("--no-info",  metavar="FASTA")
    arg = ap.parse_args()

    sys.stdout.write("sequence\tlabel\n")
    if arg.positive:
        parse_fasta(Path(arg.positive), True)
    if arg.negative:
        parse_fasta(Path(arg.negative), False)
    if arg.no_info:
        parse_fasta(Path(arg.no_info), None)
