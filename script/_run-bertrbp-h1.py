from argparse import ArgumentParser
import sys
import numpy as np


if "__main__" == __name__:
    ap = ArgumentParser()
    ap.add_argument("--test", metavar="FASTA", required=True)
    ap.add_argument("--pred", metavar="NPY", required=True)
    arg = ap.parse_args()


    pred = np.load(arg.pred)
    with open(arg.test) as fi:

        n = 0
        while True:
            head = fi.readline().strip(" \n")[1:]
            seq  = fi.readline().strip(" \n")
            if (not head) or (not seq):
                break

            for i, base in enumerate(seq):
                sys.stdout.write("%s\t%d\t%d\t%f\t%f\n" % (
                    head,
                    i,
                    0,
                    pred[n],
                    np.nan,
                    ))

            n += 1

