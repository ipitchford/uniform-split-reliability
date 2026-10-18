"""Independent ordinary-integer sanity checks; explicit checks survive Python -O."""
from itertools import combinations
from math import comb
from pathlib import Path
import json,time,sys
ROOT=Path(__file__).resolve().parent
def require(ok,why):
    if not ok:raise RuntimeError(why)
def connected_coefficients(n,edges):
    polys={0:[0]};ec={0:0}
    for mask in range(1,1<<n):
        ec[mask]=sum(bool(mask>>a&1 and mask>>b&1) for a,b in edges)
        p=[comb(ec[mask],j) for j in range(ec[mask]+1)];r=mask&-mask;sub=(mask-1)&mask
        while sub:
            if sub&r:
                other=mask^sub
                for i,a in enumerate(polys[sub]):
                    for j in range(ec[other]+1):p[i+j]-=a*comb(ec[other],j)
            sub=(sub-1)&mask
        require(all(x>=0 for x in p),'negative connected count');polys[mask]=p
    return polys[(1<<n)-1]
def cone(n,s):return connected_coefficients(n,list(combinations(range(n-1),2))+[(n-1,v) for v in range(s)])
def det(a):
    a=[r[:] for r in a];last=1
    for k in range(len(a)-1):
        require(a[k][k]!=0,'zero pivot');pivot=a[k][k]
        for i in range(k+1,len(a)):
            for j in range(k+1,len(a)):
                v=a[i][j]*pivot-a[i][k]*a[k][j];require(v%last==0,'noninteger Bareiss');a[i][j]=v//last
        last=pivot
    return a[-1][-1]
def trees(n,edges):
    a=[[0]*n for _ in range(n)]
    for u,v in edges:
        a[u][u]+=1;a[v][v]+=1;a[u][v]-=1;a[v][u]-=1
    return det([r[:-1] for r in a[:-1]])
def components(n,edges):
    parent=list(range(n))
    def find(x):
        while parent[x]!=x:x=parent[x]
        return x
    for u,v in edges:parent[find(u)]=find(v)
    return [find(i) for i in range(n)]
def split_brute(n,edges,s,t):
    out=[0]*(len(edges)+1)
    for mask in range(1<<len(edges)):
        c=components(n,[e for i,e in enumerate(edges) if mask>>i&1])
        if len(set(c))==2 and c[s]!=c[t]:out[mask.bit_count()]+=1
    return out
start=time.monotonic();report={'optimized':not __debug__,'submitted_code_imported':False,'zero_star':[],'excluded':[]}
for n in range(7,31):
    q=n-2;core=list(combinations(range(q),2));w=q
    distinct=core+[(w,0),(w,1)];same=core+[(w,0),(w,0)];competitor=[e for e in distinct if e!=(2,3)]+[(w,2)]
    a,b,c=trees(q+1,distinct),trees(q+1,same),trees(q+1,competitor)
    require(a==2*(q+1)*q**(q-3),'distinct attachments formula')
    require(b==2*q**(q-2),'same attachment formula')
    require(c-a==(q-4)*(q+1)**2*q**(q-5) and c>a>b,'zero star margin')
    report['zero_star'].append({'n':n,'distinct':a,'same':b,'competitor':c,'margin':c-a})
for N in range(7,11):
    low=cone(N,2)
    for eps in (0,1):
        h=N-2-eps;d=h-2;high=cone(N,h);dif=[0]*(len(high)+1)
        for i in range(1,d+2):
            for j,v in enumerate(low):dif[i+j]+=comb(d+1,i)*v
        for j,v in enumerate(high):dif[j+1]-=v
        report['excluded'].append({'N':N,'epsilon':eps,'d':d,'holds':all(x>=0 for x in dif),'difference':dif})
require([(x['N'],x['epsilon'],x['d']) for x in report['excluded'] if x['holds']]==[(7,0,3),(7,1,2),(8,1,3),(9,1,4),(10,1,5)],'five exceptions')
loopchecks=0;bad_without_loop=0
for n in range(3,6):
    all_edges=list(combinations(range(n),2))
    # Complete, path and cycle cover adjacent and nonadjacent terminal identification.
    for edges in [all_edges,[(i,i+1) for i in range(n-1)],[(i,(i+1)%n) for i in range(n)]]:
        original=connected_coefficients(n,edges)
        for s,t in combinations(range(n),2):
            vertices=[v for v in range(n) if v!=t];mapping={v:i for i,v in enumerate(vertices)};mapping[t]=mapping[s]
            quotient=[(mapping[u],mapping[v]) for u,v in edges if mapping[u]!=mapping[v]]
            c=connected_coefficients(n-1,quotient);hasloop=(s,t) in edges or (t,s) in edges
            padded=[(c[i] if i<len(c) else 0)+(c[i-1] if hasloop and 0<=i-1<len(c) else 0) for i in range(len(edges)+1)]
            exact=split_brute(n,edges,s,t)
            require([padded[i]-original[i] for i in range(len(original))]==exact,'loop-padded identity')
            if hasloop and [(c[i] if i<len(c) else 0)-original[i] for i in range(len(original))]!=exact:bad_without_loop+=1
            loopchecks+=1
require(bad_without_loop>0,'loop omission negative control ineffective')
report.update(loop_checks=loopchecks,loop_omission_rejections=bad_without_loop,status='passed',elapsed_seconds=time.monotonic()-start)
print(json.dumps(report,indent=2))
