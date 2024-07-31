#!/bin/env python3
# usage: python _run-rnapromo-h0.py <rnapromo_out> <test.fa>

import sys
inf = float("inf")

ifname = sys.argv[1]
faname = sys.argv[2]

def get_correct():
    correct = {}
    with open(faname) as fi:
        b,e=-1,-1
        for l in fi:
            if l.startswith(">"):
                l=l[1:].strip()
                try:
                    flds=dict(fld.split(':') for fld in l.split(';'))
                    ID=l
                    b,e=list(map(int,flds["motif-site"].split('-')))
                except:
                    ID=l
                    b,e=-1,-1
            else: # sequence
                correct[ID]=[b,e,len(l.strip())]
    return correct

def get_hits():
    hits = {}
    with open(ifname) as fi:
        for l in fi:
            if not l.strip(): continue
            fld = l[1:].strip().split('\t')
            ID=fld[0]
            hits[ID] = []
            b,e,score=0,0,0
            try:
                b,e = int(fld[2]), int(fld[3])+1
                score = float(fld[1])
            except:
                continue
            hits[ID] += [[b,e,score]]
    return hits

def calc_ox(correct, hits):
    for ID in correct.keys():
        short_ID=ID[:42]
        try: hits[short_ID]
        except KeyError: hits[short_ID] = []
        b,e,l = correct[ID]
        scores=[max([-inf]+[h[2] for h in hits[short_ID] if h[0]<=i<h[1]]) for i in range(l)]
        for i in range(l):
            ox = 1 if b<=i<e else 0
            rank = -1
            sys.stdout.write('\t'.join(map(str,[ID,i,ox,max(scores),scores[i]]))+'\n') # only best match

def main():
    correct = get_correct()
    hits = get_hits()
    calc_ox(correct, hits)

main()

