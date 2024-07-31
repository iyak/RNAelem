# usage ./gen-preset
import itertools
import sys

minlen = int(sys.argv[1])#7
maxlen = int(sys.argv[1])#7
maxnnode = 99
maxsibl = 99

def enum_recu(dfs):
    if maxnnode<len(dfs)-1 or maxlen<2*len(dfs): return
    print_dfs(dfs)
    parent = len(dfs)-1
    while -1 != parent:
        enum_recu(dfs + [parent])
        parent = dfs[parent]
    dfs.pop()

def tree_to_brackets(rdfs, parent):
    brackets = []
    for child in rdfs[parent]:
        brackets += ['('] + tree_to_brackets(rdfs, child) + [')']
    return brackets

def print_dfs(dfs):
    nbp = (len(dfs)-1)*2
    maxnb = maxlen-nbp
    rdfs = [[] for _ in dfs]
    for i,di in enumerate(dfs[1:]): rdfs[di].append(i+1)
    brackets = tree_to_brackets(rdfs, 0)
    for nb in range(maxnb+1):
        for c in itertools.combinations(range(nbp+nb), nbp):
            rss = ['.' for _ in range(nbp+nb)]
            for i,ci in enumerate(c):
                rss[ci] = brackets[i]
            filter_print(rss,dfs,rdfs)

def filter_print(rss,dfs,rdfs):
    rss = ''.join(rss)#.replace('()','(*)').replace('(..)','(.*.)')
    if '()' in rss: return
    if '(.)' in rss: return
    if '(..)' in rss: return
    if ')(' in rss: return
    elif len(rss.replace('*','')) < minlen: return
    elif maxlen < len(rss.replace('*','')): return
    elif maxsibl < max(len(x) for x in rdfs): return
    else:
        print(rss)

def enum():
    enum_recu([-1])

enum()
