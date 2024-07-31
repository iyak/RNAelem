#!/home/iyak/.pyenv/shims/python
# usage: ./analyse_dat.py <raw>
# description: post-data process for RNAelem to evaluate the performance to infer the motif position

#$ -cwd
#$ -e /home/iyak/.ugeerr
#$ -o /home/iyak/.ugeout
#$ -S /home/iyak/.pyenv/shims/python
##$ -m es -M h.m@tuta.io
#$ -V
#$ -l s_vmem=2G,mem_req=2G
#$ -r no
#$ -pe def_slot 1

import sys,numpy as np
nan=float("nan")
inf=float("inf")
raw=sys.argv[1]
### Util
def cry(*x):
    stderr.write(" ".join(map(str,x))+'\n')
def die(*x):
    cry(*x)
    exit()
def _run(cmd,stdin=None):
    cry("sh:",dedent(cmd))
    if not None==stdin:cry(dedent(stdin))
    p=Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE,stdin=PIPE)
    o,e=p.communicate(stdin.encode("utf-8") if not None==stdin else None)
    return o.decode("utf-8"),e.decode("utf-8")
def run(cmd,stdin=None):
    o,e=_run(cmd,stdin)
    stdout.write(o)
    stderr.write(e)
    return o,e
def nlines(fname):
    count=0
    for _ in open(fname):count+=1
    return count
def n_fa(faname):
    count=0
    for l in open(faname):count+=1 if '>'==l[0] else 0
    return count
def chunk(x,n):
    j,r=0,[]
    for w in [len(x)//n+int(i<len(x)%n) for i in range(n)]:
        r+=[x[j:j+w]]
        j+=w
    return r
def parse_fasta(fasta_name):
    fh=open(fasta_name)
    faiter=(x[1] for x in groupby(fh,lambda l:'>'==l[0]))
    for header in faiter:
        header=header.__next__()[1:].strip()
        seq="".join(s.strip() for s in faiter.__next__())
        yield header,seq
def parse_RNAelem_raw(fname):
    with open(fname) as fi:
        while True:
            _ID=fi.readline()
            if not _ID: break
            ID=_ID.strip().split(': ',1)[1][1:]
            S=eval(fi.readline().split(': ',1)[1])
            E=eval(fi.readline().split(': ',1)[1])
            I=eval(fi.readline().split(': ',1)[1])
            psi=eval(fi.readline().split(': ',1)[1])
            region=list(map(int, fi.readline().strip().split(': ',1)[1].split('-')))
            exist=float(fi.readline().split(': ',1)[1])
            seq=fi.readline().strip('\n').split(': ',1)[1]#seq
            rss=fi.readline().strip('\n').split(': ',1)[1]#rss
            mot=fi.readline().strip('\n').split(': ',1)[1]#aligned motif
            yield ID,S,E,I,psi,region,exist,seq,rss,mot

ox_p,s_p,ox_s,s_s=[],[],[],[]
for ID,_,_,I,_,(b,e),exist,seq,rss,mot in parse_RNAelem_raw(raw):
    pos=[]
    rss_pred = "".join("X" if not b<=i<e else "G" if "*"==m else r for i,(r,m) in enumerate(zip(rss,mot)))
    for fld in ID.strip().split(";"):
        if fld.startswith("motif-site"):
            pos=[list(map(int, r.split("-"))) for r
                    in fld.split(":")[1].split(",")]
    for i in range(len(seq)):
        ox=1 if any(be[0]<=i<be[1] for be in pos) else 0
        sys.stdout.write('\t'.join(map(str,[ID,i,ox,exist,I[i],rss_pred[i]]))+'\n')
