#!/home/iyak/.pyenv/shims/python
# usage: ./analyse_dat.py <raw output>
# description: post-data process for RNAelem to evaluate the performance to infer the existence of motifs in each input sequence

#$ -cwd
#$ -e /home/iyak/.ugeerr
#$ -o /home/iyak/.ugeout
#$ -S /home/iyak/.pyenv/shims/python
##$ -m es -M h.m@tuta.io
#$ -V
#$ -l s_vmem=0.5G,mem_req=0.5G
#$ -r no
#$ -pe def_slot 1

import sys,numpy as np
inf=float("inf")
nan=float("nan")
for fname in sys.argv:
    with open(fname) as fi:
        while True:
            try:
                begin=list(map(np.exp,eval(fi.readline().strip().split(": ")[1])))[:-1]
                end=list(map(np.exp,eval(fi.readline().strip().split(": ")[1])))[:-1]
                inner=list(map(np.exp,eval(fi.readline().strip().split(": ")[1])))[:-1]
                motif_region=map(int,fi.readline().strip().split(": ")[1].split("-"))
                exist_prob=float(fi.readline().strip().split(": ")[1])
                annot=dict(fld.split(':') for fld in fi.readline().strip()[5:].split(';'))
                seq=fi.readline().strip().split(": ")[1]
                rss=fi.readline().rstrip("\n").split(": ")[1]
                motif=fi.readline().rstrip("\n").split(": ")[1]
            except: break
            try:
                b,e=map(int,annot["motif-site"].split('-'))
                positive=1
            except:
                positive=0
            sys.stdout.write("\t".join(map(str,[positive,0,exist_prob]))+"\n")
