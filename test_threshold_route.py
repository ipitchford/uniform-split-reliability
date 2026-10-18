"""For n in a range, enumerate threshold graphs T' on w = n-1-delta vertices with delta+eps edges
(eps=0: k=n-1; eps=1: k=n), embed in K_{n-1} (vertices V minus s), and test
   S_delta(p) * Rel(K_{n-1} - T')  >=  Rel(K_{n-1} - K_{1,1+eps})   for p in (0,1),
where S_delta = 1+(1-p)+...+(1-p)^{delta-1}.  Threshold graphs are generated as all graphs with the
given parameters and filtered by the P4/C4/2K2-free characterisation."""
import sys, subprocess, numpy as np, networkx as nx
from multigraph_lib import conn_poly_multigraph, simple_mult
def is_threshold(F):
    # threshold iff no induced P4, C4, 2K2
    nodes = list(F.nodes())
    for a in nodes:
        for b in nodes:
            if a>=b: continue
            for c in nodes:
                if c in (a,b): continue
                for d in nodes:
                    if d in (a,b,c) or d<c: continue
                    S = F.subgraph([a,b,c,d])
                    e = S.number_of_edges()
                    if e == 2 and max(dict(S.degree()).values()) == 1: return False      # 2K2
                    if e == 3 and sorted(dict(S.degree()).values()) == [1,1,2,2]: return False   # P4
                    if e == 4 and sorted(dict(S.degree()).values()) == [2,2,2,2]: return False   # C4
    return True
def rel_curve(q, missing, ps):
    adj=[0]*q
    for u in range(q):
        for v in range(u+1,q):
            if (u,v) in missing or (v,u) in missing: continue
            adj[u]|=1<<v; adj[v]|=1<<u
    c = conn_poly_multigraph(q, simple_mult(q, adj)); m=len(c)-1
    return sum(c[i]*ps**i*(1-ps)**(m-i) for i in range(m+1))
ps = np.linspace(0.002, 0.998, 300)
for n in range(7, int(sys.argv[1])+1):
    q = n-1
    for eps in (0,1):
        target = rel_curve(q, {(0,1)} if eps==0 else {(0,1),(0,2)}, ps)
        worst = (9, None)
        for delta in range(2, n-2):
            w = n-1-delta; E = delta+eps
            if w < 2 or E > w*(w-1)//2: continue
            g6s = subprocess.run(['./bin/geng','-q',str(w),'%d:%d'%(E,E)], capture_output=True, text=True).stdout.split()
            for g6 in g6s:
                T = nx.from_graph6_bytes(g6.encode())
                if not is_threshold(T): continue
                missing = {(u,v) for u,v in T.edges()}   # T' on vertices 0..w-1 of K_q
                R = rel_curve(q, missing, ps)
                S = sum((1-ps)**j for j in range(delta))
                ratio = (S*R/target).min()
                if ratio < worst[0]: worst = (ratio, (delta, g6, sorted(T.edges())))
        print("n=%d eps=%d  worst ratio %.4f at %s" % (n, eps, worst[0], worst[1]))
