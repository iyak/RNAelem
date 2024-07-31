#! /usr/bin/env python
# usage: python _run-graphprot-h0.py <GraphProt.profile> <GraphProt.predictions> <test.fa>

import sys
with open(sys.argv[1]) as prof, open(sys.argv[2]) as pred, open(sys.argv[3]) as fa:
    while True:
        fa0=fa.readline().strip()
        fa1=fa.readline().strip()
        if (not fa0) or (not fa1): break
        pos=[]
        for fld in fa0[1:].strip().split(";"):
            if fld.startswith("motif-site"):
                pos=[list(map(int, r.split("-"))) for r
                        in fld.split(":")[1].split(",")]
        L = len(fa1)
        seq_wise=pred.readline().strip().split()[2]
        for fld in [prof.readline().strip().split() for i in range(L)]:
            # The next conditional break deals with a bug of GraphProt;
            # there is no output in GraphProt.profile for the last base of the last sequence.
            if not fld: break
            ox=1 if any(be[0]<=int(fld[1])<be[1] for be in pos) else 0
            sys.stdout.write('\t'.join(map(str,[fa0[1:],fld[1],ox,seq_wise,fld[2]]))+'\n')
