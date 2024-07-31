#usage: python _run-rnacontext-h3.py <model_res_10.txt> <nsites>
import sys,re,numpy as np
from mymo import util
mname=util.fullpath(sys.argv[1])
nsites=int(sys.argv[2])
table={}
def sumexp(v):return sum(np.exp(vi) for vi in v)
with open(mname) as m:
    while True:
        if "Base Parameters"==m.readline().strip():
            while True:
                l=m.readline().strip()
                if not l: break
                alph=l.split()[0].replace('U','T')
                table[alph]=list(map(float,l.replace('-'," -").split()[1:]))
            break
        if not m.readline():raise RuntimeError("no Base Parameters")
if not all(len(table[alph])==len(table[t]) for t in table):
    raise RuntimeError("not same len",[len(t) for t in table])
with open(mname+".tamo",'w') as fo:
    fo.write("Log-odds matrix for Motif 0 a (%d)\n"%nsites)
    fo.write('#\t\t'+'\t'.join(map(str,range(len(table[alph]))))+'\n')
    for k in "ACTG":fo.write('#'+k+'\t'+'\t'.join(map(str,table[k]))+'\n')
    fo.write("Source: a a Mono-syllabic 0\n")
_,e,_=util.po("tamo2meme %s.tamo"%mname,stdout="%s.meme.dna"%mname)
sys.stderr.write(e)
if not e:util.rm("%s.tamo"%mname)

with open("%s.meme.dna"%mname) as fi,open("%s.meme"%mname,'w') as fo:
    while True:
        l=fi.readline()
        if not l:break
        if l.startswith("ALPHABET="):
            fo.write("ALPHABET= ACGU\n")
            continue
        elif l.startswith("Background letter frequencies"):
            fo.write(l)
            l=fi.readline()
            ll=l.split()
            ll[6]='U'
            fo.write(' '.join(ll))
        else:fo.write(l)
util.rm("%s.meme.dna"%mname)
