#!/usr/bin/env python
# usage draw_motif.py <elem_out/pattern-*> <eps> <svg>
import numpy as np
from subprocess import Popen,PIPE
from sys import stderr,stdout,argv,platform
from re import compile
inf,nan=float("inf"),float("nan")
dname=argv[1]
eps=argv[2]
svg=argv[3]
def cry(*x):
    stderr.write(' '.join(map(str,x))+'\n')
def die(*x):
    cry(*x)
    exit()
def _run(cmd,stdin=None):
    cry("sh:",cmd)
    if not None==stdin:cry(stdin)
    p=Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE,stdin=PIPE)
    o,e=p.communicate(input=None if None==stdin else bytes(stdin,encoding='utf-8'))
    return o.decode("utf-8"),e.decode("utf-8")
def run(cmd,stdin=None):
    o,e=_run(cmd,stdin)
    stdout.write(o)
    stderr.write(e)
    return o,e
def is_suitable_font(ttf):
    if "Italic" in ttf or "italic" in ttf:
        return False
    if "Bold" in ttf or "bold" in ttf:
        return False
    return True
def get_gothic_ttf():
    if platform.startswith("win") or platform.startswith("cygwin"):
        print("for windows, I don't get path for fonts.")
        pass
    elif platform.startswith("darwin"):
        ttfs,e=_run("find /Library/Fonts -name *.ttf")
        stderr.write(e)
        for ttf in ttfs.strip().split("\n"):
            if is_suitable_font(ttf):
                return ttf
        else:
            cry("could not find .ttf file.")
            cry("please set manually")
    else:
        ttfs,e=_run("find /usr/share/fonts -name *.ttf")
        stderr.write(e)
        for ttf in ttfs.strip().split("\n"):
            if is_suitable_font(ttf):
                return ttf
        else:
            cry("could not find font file.")
            cry("please set manually")
def parse_raw(raw):
    with open(raw) as fi:
        while True:
            try:yield dict(fi.readline().strip('\n').split(": ",1) for _ in range(10))
            except:return
# build pwm
fld=dict()
#fld=dict(l.strip().split(": ",1) for l in open(dname+"/log"))
with open(dname+"/log") as f:
    for l in f:
        try:
            k,v=l.strip('\n').split(": ",1)
            fld[k]=v
        except:
            cry("wrong format:",l)
N=eval(fld["E[N]"])
p=fld["motif pattern"]
z=0.
L_gap=[0 for _ in range(p.count('*'))]
L_gap_div=0.
pos_gap=[i+1 for i,pi in enumerate(p) if '*'==pi]
for fld in parse_raw(dname+"/train.raw"):
    try:
        z+=eval(fld["exist prob"])
        y=eval(fld["psihat"])
        for i,pg in enumerate(pos_gap):
            L_gap[i]+=z*sum(1 for yi in y if yi==pg)
        L_gap_div+=z
    except:pass
L_gap=[l/L_gap_div for l in L_gap]
pwm,stack,loop=[],[],[]
k,l=1,0
for j in range(len(p)):
    if '('==p[j]:
        pwm+=[[]]
        loop+=[1]
        stack+=[j]
    elif ')'==p[j]:
        i=stack.pop()
        pwm[i]=[N[k][4],N[k][0],N[k][1]+N[k][2],N[k][3]+N[k][5]]
        loop[i]=sum(N[k])/z
        pwm+=[[N[k][5],N[k][1],N[k][0]+N[k][3],N[k][2]+N[k][4]]]
        loop+=[sum(N[k])/z]
        k+=1
    elif '.'==p[j] or '_'==p[j]:
        pwm+=[N[k]]
        loop+=[sum(N[k])/z]
        k+=1
    elif '*'==p[j]:
        pwm+=[[1,1,1,1]]
        loop+=[L_gap[l]]
        l+=1
    else:raise die("unknown pattern:",p[j])
# get seq
if '_' in p:# no RSS
    run("RNAelem-plot >%s"%eps)
else:
    seq,rss="",""
    for v,l,pi in zip(pwm,loop,p):
        s=sum(v)
        w=[vi/s for vi in v if 0<vi]
        e=sum(-wi*np.log2(wi) for wi in w)
        if e<1.5:
            order=sorted(range(len(w)),key=lambda k:-w[k])
            if w[order[1]]*2<w[order[0]]:
                seq+="ACGU"[order[0]]*round(l)
                rss+=pi*round(l)
            else:
                seq+="AMRWMCSYRSGKWYKU"[order[0]*4+order[1]]*round(l)
                rss+=pi*round(l)
        else:
            seq+=' '*round(l)
            rss+=pi*round(l)
    run("RNAelem-plot '%s' '%s' >%s"%(seq,rss,eps))
# gen seq logo
stack,logo_seq,logo_color,logo_val,logo_meta=[],[],[],[],[]
k=1
if not p.startswith("*"):
    p = "*" + p
    loop += [1] + loop
if not p.endswith("*"):
    p = p + "*"
    loop += loop + [1]
for j in range(len(p)):
    if'('==p[j]:
        stack+=[j]
        logo_seq+=["[CG,GC,AU,UA,GU,UG]"]
        logo_color+=["[[,lightgray],[,lightgray],[,lightgray],[,lightgray],[,lightgray],[,lightgray]]"]
        logo_val+=['[]'] #fill at closing pair
        logo_meta+=['(']
    elif')'==p[j]:
        i=stack.pop()
        logo_seq+=["[CG,GC,AU,UA,GU,UG]"]
        logo_color+=["[[lightgray,],[lightgray,],[lightgray,],[lightgray,],[lightgray,],[lightgray,]]"]
        logo_val[i]='['+','.join(map(str,N[k]))+']'
        logo_val+=['['+','.join(map(str,N[k]))+']']
        logo_meta+=[')']
        k+=1
    elif'.'==p[j] or '_'==p[j]:
        logo_seq+=["[A,C,G,U]"]
        logo_color+=["[,,,]"]
        logo_val+=['['+','.join(map(str,N[k]))+']']
        logo_meta+=[' . ']
        k+=1
    elif'_'==p[j]:
        logo_seq+=["[A,C,G,U]"]
        logo_color+=["[,,,]"]
        logo_val+=['['+','.join(map(str,N[k]))+']']
        logo_meta+=['']
        k+=1
    elif'*'==p[j]:
        logo_seq+=["[A,C,G,U]"]
        logo_color+=["[,,,]"]
        logo_val+=["["+",".join(map(str,N[0]))+"]"]
        logo_meta+=['*']
    else:raise die("unknown pattern:",p[j])
stdin="""
seq:%s
val:%s
color:%s
meta:%s
"""%('['+','.join(s for s in logo_seq)+']',
    '['+','.join(s for s in logo_val)+']',
    '['+','.join(s for s in logo_color)+']',
    '['+','.join(s+str(round(l)) if 1<round(l) else s for s,l in zip(logo_meta,loop))+']')
run("RNAelem-logo %s >%s"%(get_gothic_ttf(),svg),stdin=stdin)
