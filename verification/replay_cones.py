"""Fail-closed wrapper of supplied algorithm; not an independent implementation."""
import sys,json,time,math,hashlib,platform,argparse
from pathlib import Path
SOURCE=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--output-dir',type=Path,required=True)
parser.add_argument('--max-N',type=int,default=50)
parser.add_argument('--inject-failure',action='store_true')
parser.add_argument('--inject-coverage-gap',action='store_true')
args=parser.parse_args()
if not 7<=args.max_N<=50:parser.error('max-N must be 7..50')
ROOT=args.output_dir.resolve();ROOT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(SOURCE))
from cone_star2 import check
def require(ok,why):
    if not ok: raise RuntimeError(why)
started=time.monotonic(); report={'algorithm':'cone_star2.py','sha256':hashlib.sha256((SOURCE/'cone_star2.py').read_bytes()).hexdigest(),'python':platform.python_version(),'optimized':not __debug__,'cases':[]}
for N in range(7,args.max_N+1):
    cache={}
    for eps in (0,1):
        out,cone=check(N,eps,cache);need=N-5-eps
        if '--inject-failure' in sys.argv and N==7 and eps==0:out[1]=False
        if '--inject-coverage-gap' in sys.argv and N==7 and eps==0:out.pop(1)
        require(set(out)==set(range(1,need+1)),f'Coverage mismatch {N} {eps}')
        require(all(out.values()),f'Required inequality failed {N} {eps}')
        h=N-2-eps;d=h-2;high,low=cone(h),cone(2)
        differences=[sum(math.comb(d+1,i)*low[j-i] for i in range(1,d+2) if 0<=j-i<len(low))-(high[j-1] if j else 0) for j in range(len(high)+1)]
        report['cases'].append({'N':N,'epsilon':eps,'required_triples':need,'scalar_inequalities':need*(math.comb(N-1,2)+N-eps),'required_passed':True,'first_excluded_d':d,'excluded_holds':all(v>=0 for v in differences),'first_negative':next(([j,v] for j,v in enumerate(differences) if v<0),None)})
        print(N,eps,'PASS',flush=True)
    (ROOT/('cone-progress-optimized.json' if not __debug__ else 'cone-progress-normal.json')).write_text(json.dumps(report,indent=2)+'\n')
report.update(status='passed',triples=sum(x['required_triples'] for x in report['cases']),scalar_inequalities=sum(x['scalar_inequalities'] for x in report['cases']),elapsed_seconds=time.monotonic()-started)
expected_triples=sum(N-5-e for N in range(7,args.max_N+1) for e in (0,1))
expected_scalars=sum((N-5-e)*(math.comb(N-1,2)+N-e) for N in range(7,args.max_N+1) for e in (0,1))
require(report['triples']==expected_triples and report['scalar_inequalities']==expected_scalars,'Coverage totals')
exceptions=[(x['N'],x['epsilon'],x['first_excluded_d']) for x in report['cases'] if x['excluded_holds']]
require(exceptions==[v for v in [(7,0,3),(7,1,2),(8,1,3),(9,1,4),(10,1,5)] if v[0]<=args.max_N],'Excluded cases differ from review')
report['excluded_exceptions']=exceptions
(ROOT/('cone-optimized.json' if not __debug__ else 'cone-normal.json')).write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ('status','triples','scalar_inequalities','elapsed_seconds')}))
