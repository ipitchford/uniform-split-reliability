"""Generate all connected graphs on n vertices with at most kmax missing edges (as complements of
graphs with <= kmax edges), write them in graph6 to a file."""
import sys, subprocess, networkx as nx
n, kmax, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
g6s = subprocess.run(['./bin/geng', '-q', str(n), '0:%d' % kmax], capture_output=True, text=True).stdout.split()
cnt = 0
with open(out, 'w') as f:
    for g6 in g6s:
        F = nx.from_graph6_bytes(g6.encode())
        H = nx.complement(F)
        if not nx.is_connected(H): continue
        f.write(nx.to_graph6_bytes(H, header=False).decode().strip() + "\n"); cnt += 1
print("wrote", cnt, "connected graphs on", n, "vertices with <=", kmax, "missing edges")
