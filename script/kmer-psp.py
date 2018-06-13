#!/usr/bin/env python
#usage:python RNAelem_fa2fq.py <positive.fa> <negative.fa>
import sys,numpy as np,itertools,random,re
from scipy.stats import fisher_exact as fish
fanameP=sys.argv[1]
fanameN=sys.argv[2] if 2<len(sys.argv) else None
kmin,kmax,thresh,rich,poor,base=3,10,5e-2,[],[],10
def gen_kmer(k):
    for kmer in itertools.product(*(["ACGUT"]*k)):
        yield ''.join(kmer)
def parse_fa(faname):
    with open(faname) as fi:
        while True:
            try:
                ann=fi.readline().strip()
                seq=fi.readline().strip()
            except:break
            if not all([ann,seq]):break
            yield ann,seq
def count_kmer(pname,nname,k):
    nT,nF,nP,nN=0,0,dict(),dict()
    for _,seq in parse_fa(pname):
        nT+=1
        kmers=set(seq[i:i+k] for i in range(len(seq)-k))
        for kmer in kmers:
            if -1!=seq.find(kmer):
                if not kmer in nP:nP[kmer]=1
                else:nP[kmer]+=1
    for _,seq in parse_fa(nname):
        nF+=1
        kmers=set(seq[i:i+k] for i in range(len(seq)-k))
        for kmer in kmers:
            if -1!=seq.find(kmer):
                if not kmer in nN:nN[kmer]=1
                else:nN[kmer]+=1
    return nT,nF,nP,nN
def calc_pval(nT,nF,nP,nN,thresh):
    min_pval,rich,poor=1.,list(),list()
    for kmer in nP.keys():
        if not kmer in nN:continue
        p=fish([[nP[kmer],nN[kmer]],
            [nT-nP[kmer],nF-nN[kmer]]])[1]
        if nN[kmer]<nP[kmer]:
            if p<=min_pval:min_pval=p
            if p<thresh:
                sys.stderr.write("+%s\t%f\n"%(kmer,p))
                rich+=[(kmer,p)]
        else:
            if p<thresh:
                sys.stderr.write("-%s\t%f\n"%(kmer,p))
                poor+=[(kmer,p)]
    return rich,poor
def output_fq(fanameP,rich,poor,base):
    for ann,seq in parse_fa(fanameP):
        qlt=np.fromiter([base for _ in range(len(seq))],
                dtype=np.float32)
        for sub in rich:
            for m in re.finditer(sub[0],seq):
                qlt[m.start():m.start()+len(sub[0])]+=1
        for sub in poor:
            for m in re.finditer(sub[0],seq):
                qlt[m.start():m.start()+len(sub[0])]-=1
        sys.stdout.write('\n'.join([
            '@'+ann[1:],seq,'+',
            ''.join([chr(max(min(ord('!')+int(round(q))
                ,ord('~')),ord('!')))for q in qlt])+'!',
            ])+'\n')
if None==fanameN:
    output_fq(fanameP,{},{},base)
    exit()
#choose best k
min_k,min_pval=-1,1.
for k in range(kmin,kmax+1):
    nT,nF,nP,nN=count_kmer(fanameP,fanameN,k)
    rich,poor=calc_pval(nT,nF,nP,nN,thresh)
    if len(rich)==0: continue
    p=min(rich,key=lambda x:x[1])[1]
    if p<min_pval:min_k,min_pval=k,p
k=min_k
sys.stderr.write("k:%d\n"%k)
#calc inner prob
nT,nF,nP,nN=count_kmer(fanameP,fanameN,k)
rich,poor=calc_pval(nT,nF,nP,nN,thresh)
output_fq(fanameP,rich,poor,base)
