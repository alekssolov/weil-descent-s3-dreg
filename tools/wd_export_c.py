"""Export a descent system to the wdclose input format."""
import random, sys
sys.path.insert(0,'.')
from wd_generic_dims import instance
def export(eqs, N, D, path):
    with open(path,'w') as f:
        f.write(f"{N} {D}\n{len(eqs)}\n")
        for quad,lin,const in eqs:
            toks=[format((1<<i)|(1<<j),'x') for (i,j) in quad]+[format(1<<i,'x') for i in lin]+(["0"] if const else [])
            f.write(" ".join(toks)+"\n")
if __name__=="__main__":
    N,n,D=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]); seed=int(sys.argv[4]) if len(sys.argv)>4 else 7000+100*N+n
    eqs,NN=instance(n,N//2,N-N//2,random.Random(seed)); export(eqs,NN,D,f"c_N{N}_n{n}.txt"); print(f"written c_N{N}_n{n}.txt (seed {seed})")
