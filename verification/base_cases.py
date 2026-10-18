"""Fresh C replay at exactly the two boundary edge counts, n=7..10.
Corpus completeness is inherited from supplied corpus and reviewer regeneration;
this run validates counts and coefficients, not a new graph enumeration.
"""
from pathlib import Path
import json,subprocess,time,hashlib,argparse
parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path,required=True);parser.add_argument('--binary',type=Path,required=True);args=parser.parse_args()
ROOT=args.output_dir.resolve();ROOT.mkdir(parents=True,exist_ok=True);(ROOT/'logs').mkdir(exist_ok=True)
SRC=Path(__file__).resolve().parents[1]
def require(ok,why):
    if not ok:raise RuntimeError(why)
def graph6(line):
    n=ord(line[0])-63;bits=[(ord(c)-63)>>s&1 for c in line[1:] for s in range(5,-1,-1)];edges=[];i=0
    require(1<=n<=12 and len(bits)>=n*(n-1)//2,'graph6 encoding')
    for v in range(1,n):
        for u in range(v):
            if bits[i]:edges.append((u,v))
            i+=1
    return n,edges
def balloon(n,edges,k):
    adj=[set() for _ in range(n)]
    for u,v in edges:adj[u].add(v);adj[v].add(u)
    for leaf in range(n):
        if len(adj[leaf])!=1:continue
        core=[v for v in range(n) if v!=leaf];missing=[(u,v) for u in core for v in core if u<v and v not in adj[u]];e=len(missing)
        if e!=k-(n-2):continue
        if e==0 or any(sum(v in edge for edge in missing)==e for v in core):return True
    return False
start=time.monotonic();result={'cases':[],'corpora':[],'fresh_compilation':True}
counts={(7,0):40,(7,1):64,(8,0):114,(8,1):220,(9,0):344,(9,1):770,(10,0):1102,(10,1):2768}
for n in range(7,11):
    source=SRC/('g%d.g6'%n if n<=8 else 'g%ddense.g6'%n);data=source.read_bytes();selected=[]
    for line in data.decode().splitlines():
        N,edges=graph6(line);require(N==n,'wrong corpus order')
        if len(edges) in ((n-1)*(n-2)//2,(n-1)*(n-2)//2-1):selected.append(line)
    inputfile=ROOT/f'base-n{n}.g6';inputfile.write_text('\n'.join(selected)+'\n')
    proc=subprocess.run([str(args.binary.resolve()),'--conn'],input=inputfile.read_bytes(),capture_output=True)
    (ROOT/f'logs/base-n{n}.stderr').write_bytes(proc.stderr);(ROOT/f'base-n{n}.txt').write_bytes(proc.stdout)
    require(proc.returncode==0,'C process failure');rows={}
    for line in proc.stdout.decode().splitlines():
        parts=line.split();m=int(parts[1]);v=list(map(int,parts[2:]));require(len(v)==m+1,'coefficient length');rows.setdefault(m,[]).append((parts[0],v))
    require(sum(map(len,rows.values()))==len(selected),'dropped graph')
    for eps in (0,1):
        m=(n-1)*(n-2)//2-eps;r=rows[m];require(len(r)==counts[n,eps],'corpus boundary count')
        lower=[min(v[j] for _,v in r) for j in range(m+1)];wins=[g for g,v in r if v==lower]
        require(len(wins)==2+eps,'minimizer count')
        require(all(balloon(n,graph6(g)[1],n-1+eps) for g in wins),'minimizer not balloon')
        result['cases'].append({'n':n,'epsilon':eps,'graphs':len(r),'minimizers':wins,'coefficientwise_minimum_attained':True,'all_minimizers_balloon':True})
    result['corpora'].append({'file':source.name,'sha256':hashlib.sha256(data).hexdigest(),'selected_graphs':len(selected)})
result.update(status='passed',graphs=sum(x['graphs'] for x in result['cases']),elapsed_seconds=time.monotonic()-start)
require(result['graphs']==5422,'total count')
print(json.dumps(result,indent=2))
