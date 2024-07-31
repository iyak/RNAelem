#!/bin/env /usr/local/package/python/3.10.5/bin/python3
# usage: python _run-cmfinder-h0.py <cmsearch_output> <test.fa>

import sys
from util import *

ifname = sys.argv[1]
faname = sys.argv[2]


def get_correct():
    correct = {}
    with open(faname) as fi:
        for l in fi:
            if l.startswith(">"):
                l=l[1:].strip()
                ID=l
                try:
                    flds=dict(fld.split(':') for fld in l.split(';'))
                    b,e=list(map(int,flds["motif-site"].split('-')))
                except:
                    b,e=-1,-1
            else: # sequence
                correct[ID]=[b,e,len(l.strip())]
    return correct


def parse_wuss(wuss):
    res = ""
    for c in wuss:
        if c in "<({["  : res += "L"
        elif c in ">)}]": res += "R"
        elif c in "_"   : res += "H"
        elif c in "."   : res += "G"
        elif c in ":"   : res += "O"
        elif c in "-"   : res += "I"
        elif c in ","   : res += "M"
        else            : res += "X"
    return res


def get_hits():
    hits = {}
    with open(ifname) as fi:
        while True:
            l = fi.readline()
            if not l:
                break
            l = l.strip("\n")
            if l.startswith("sequence:"):
                #ID = l[10:].strip().split(";")[0].split(':')[1]
                #hits[ID] = []
                ID=l[10:].strip()
                hits[ID]=[]
            elif l.startswith("hit"):
                b,e = map(int, l.strip().split()[3:5])
                e += 1 # pythonify
                # if b < e: # forward strand
                #     bit = float(l.strip().split()[5])
                #     hits[ID] += [[b,e,bit]]
                if e < b:
                    b,e = e,b

                bit = float(l.strip().split()[5])
                ll = fi.readline().strip(" \n")
                _  = fi.readline()
                _  = fi.readline()
                sq = fi.readline().split()[1]
                i_dels = [i for i,c in enumerate(sq) if c == "-"]
                ll = [lli for i,lli in enumerate(ll) if i not in i_dels]

                rss_type = parse_wuss(ll)
                hits[ID] += [[b,e,bit,rss_type]]
    return hits


def calc_ox(gt, hits):
    for IDc in gt.keys():
        for k in hits.keys():
            if k in IDc:
                IDh=k
                break
        else:
            IDh=IDc
            hits[IDh]=[]
        b = gt[IDc]["b"]
        e = gt[IDc]["e"]
        l = len(gt[IDc]["seq"])

        # Find max bit motif
        max_bit, max_imotif = 0, 0
        rss = ""
        for i, (_, _, bit, _) in enumerate(hits[IDh]):
            if max_bit < bit:
                max_bit, max_imotif = bit, i
        max_motif = hits[IDh][max_imotif] if hits[IDh] else None
        rss_pred = "%s%s%s" % (
                "X" * (max_motif[0]),
                max_motif[3], # rss_type
                "X" * (l - max_motif[1])) if hits[IDh] else "X" * (l)

        # max_bit = max(h[2] for h in [[-1,-1,0.,""]] + hits[IDh])
        for i in range(l):
            ox = 1 if b<=i<e else 0
            rank = -1
            bit = max([0.]+[h[2] for h in hits[IDh] if h[0]<=i<h[1]])
            sys.stdout.write("%s\t%d\t%d\t%f\t%f\t%s\n" % (IDc,i,ox,max_bit,bit,rss_pred[i]))


def main():
    gt = get_gt_from_fa(faname)
    hits = get_hits()
    calc_ox(gt, hits)


main()
