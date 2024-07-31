#!/bin/env python
#usage: run-rnacontext-h0.py <RNAcontext/outut (dir)> <t.fa>

import sys,os,numpy as np,math

dname=sys.argv[1]
faname=sys.argv[2]
res_fname="{dname}/test_res_10.txt".format(**locals())
inf,nan=float("inf"),float("nan")

with open(res_fname) as res, open(faname) as fa:
    while True:
        fa0=fa.readline().strip()
        fa1=fa.readline().strip()
        if (not fa0) or (not fa1): break
        pos=[]
        for fld in fa0[1:].strip().split(";"):
            if fld.startswith("motif-site"):
                pos=[list(map(int, r.split("-"))) for r
                        in fld.split(":")[1].split(",")]
        L=len(fa1)
        seq_wise=res.readline().strip().split()[1]
        for i in range(L):
            ox=1 if any(be[0]<=i<be[1] for be in pos) else 0
            sys.stdout.write('\t'.join(map(str,[
                fa0[1:],i,ox,seq_wise,nan
                ]))+'\n')
