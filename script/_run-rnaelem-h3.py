#!/home/iyak/.pyenv/shims/python
# usage: ./analyse_dat.py <raw output>
# description: post-data process for RNAelem to evaluate the performance to infer the motif position

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
nan,inf=float("nan"),float("inf")
for fname in sys.argv:
    with open(fname) as fi:
        while True:
            try:
                _=list(map(np.exp,eval(fi.readline().strip().split(": ")[1])))[:-1]
                _=list(map(np.exp,eval(fi.readline().strip().split(": ")[1])))[:-1]
                inner=list(map(np.exp,eval(fi.readline().strip().split(": ")[1])))[:-1]
                motif_region=map(int,fi.readline().strip().split(": ")[1].split("-"))
                _=float(fi.readline().strip().split(": ")[1])
                annot=dict(fld.split(':') for fld in fi.readline().strip()[5:].split(';'))
                seq=fi.readline().strip().split(": ")[1]
                _=fi.readline().rstrip("\n").split(": ")[1]
                _if=fi.readline().rstrip("\n").split(": ")[1]
            except IndexError: break
            #inner_cond=[p/exist_prob for p in inner]
            b,e=-1,-1
            try: b,e=map(int,annot["decoy-site"].split('-'))
            except: pass
            TP,FP,TN,FN=0,0,len(seq)-(e-b),e-b
            for j,i in enumerate(np.argsort(inner)[::-1]):
                x=1 if b<=i<e else 0
                sys.stdout.write("\t".join(map(str,[x,j,inner[i]]))+"\n")
