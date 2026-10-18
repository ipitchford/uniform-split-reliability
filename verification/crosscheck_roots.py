"""Compare submitted root-set arithmetic with independent edge-subset enumeration."""
import sys,json,time
from pathlib import Path
from itertools import combinations
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from cone_star2 import check
def require(ok,why):
    if not ok:raise RuntimeError(why)
started=time.monotonic();rows=[]
for M in range(1,7):
    edges=list(combinations(range(M),2));counts=[[0]*(len(edges)+1) for _ in range(M+1)]
    for mask in range(1<<len(edges)):
        parts=[1<<i for i in range(M)]
        for j,(u,v) in enumerate(edges):
            if mask>>j&1:
                a=next(p for p in parts if p>>u&1);b=next(p for p in parts if p>>v&1)
                if a!=b:parts.remove(a);parts.remove(b);parts.append(a|b)
        for a in range(1,M+1):
            if all(p&((1<<a)-1) for p in parts):counts[a][mask.bit_count()]+=1
    cache={};_,cone=check(M+1,0,cache);gamma,_,bits=next(iter(cache.values()));digitmask=(1<<bits)-1
    for a in range(1,M+1):require(counts[a]==[(gamma[M,a]>>(bits*j))&digitmask for j in range(len(edges)+1)],f'rootset mismatch {M} {a}')
    require(not any(cone(0)),'zero-degree cone must be disconnected')
    rows.append({'M':M,'root_sizes':M,'edge_subsets':1<<len(edges)})
print(json.dumps({'status':'passed','optimized':not __debug__,'rows':rows,'root_pairs':sum(x['root_sizes'] for x in rows),'edge_subsets':sum(x['edge_subsets'] for x in rows),'elapsed_seconds':time.monotonic()-started},indent=2))
