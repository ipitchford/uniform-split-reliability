"""Exact finite sanity checks of interval proof algebra; not a universal proof."""
from fractions import Fraction as F
from math import comb

def require(condition, message):
    if not condition:
        raise RuntimeError(message)

require(6*F(1,2)**4+15*F(1,2)**7+10*F(1,2)**8 == F(17,32), "M6 cut sum")
require(F(17,32)<F(4,7), "M6 margin")
require(F(11,32)<F(4,7), "M7+ margin")
cases=0
for den in range(2,21):
    for num in range(1,den):
        q=F(num,den)
        for s in range(3,15):
            for d in range(1,15):
                A=sum(q**i for i in range(d+1))
                D=1-(1-q**(s+d))/(A*(1-q**s))
                product=q*(1-q**d)/(1-q**(d+1))*(1-q**(s-1))/(1-q**s)
                require(D==product, "factor identity")
                require(D>=q/(1+q+q*q), "uniform margin")
                cases+=1
for N in range(7,101):
    for eps in (0,1):
        h=N-2-eps
        for d in range(1,N-4-eps):
            s=h-d
            R=F((d+1)*s,h)*F(N-1,N)**d
            require(R>1,"tree ratio")
            require(d+comb(N-1,2)+s-N+1==comb(N-1,2)-1-eps,"exponent identity")
print(f"PASS: {cases} rational factor/margin identities; cut constants; target tree ratios N=7..100.")
print("These finite checks support algebra only. The universal proof is in paper/split_reliability.md.")
