"""Strict small-capacity endpoint packing; see push027_audit_strict_proof.tex."""
import itertools
from fractions import Fraction as F

def labels_and_pairs(base_labels,K,h,factor):
 lw=lambda c:sum((i+h)*n for i,n in enumerate(c))
 labels=[(c,u,cap) for c,u in base_labels for cap in range(factor*(K-lw(c)-u)+1)]
 pairs=[(label,j) for label in labels for j in range(1,K+1) if (label[0][j-h-1]>0 if j>h else factor*(j-1)<label[2])]
 return labels,pairs

def packing_values(K,h,factor,label,L,U):
 """L,U map class j to exact endpoint envelope costs at fixed time/load."""
 c,u,cap=label;R0=K-sum((i+h)*n for i,n in enumerate(c))-u
 Js=list(range(1,min(h,(cap-1)//factor+1)+1)) if cap else []
 if cap==factor*R0:
  V=[0]*R0
  for w in range(1,R0):
   opts=[V[w-1]]
   for j in Js:
    if j<=w:opts.append(U[j]+V[w-j])
    if 1<j<=w+1:opts.append(L[j]+V[w-j+1])
   V[w]=max(opts)
  def ad(d):return sum((n-d[i])*L[i+h+1]+d[i]*U[i+h+1] for i,n in enumerate(c) if n)
  vals=[]
  for d in itertools.product(*[range(n+1) for n in c]):
   used=sum(d)
   if 0<used<=R0:vals.append(ad(d)+V[R0-used])
  base=ad((0,)*len(c));vals.append(base+V[R0-1])
  vals += [base+U[j]+V[R0-j] for j in Js]
  return vals
 V=[F(0)]*(cap+1)
 types=[((j-1)*factor+k,F(factor-k,factor)*L[j]+F(k,factor)*U[j],k) for j in Js for k in range(factor+1) if 0<(j-1)*factor+k<=cap]
 for w in range(1,cap+1):V[w]=max([V[w-1]]+[v+V[w-size] for size,v,k in types if size<=w and not (w==cap and k==0)])
 big=[i+h+1 for i,n in enumerate(c) for _ in range(n)]
 vals=[]
 for offsets in itertools.product(range(factor+1),repeat=len(big)):
  used=sum(offsets)
  if used<=factor*R0:
   ad=sum(F(factor-k,factor)*L[j]+F(k,factor)*U[j] for j,k in zip(big,offsets))
   vals.append(ad+V[min(cap,factor*R0-used)])
 return vals

def minorant_cases(K,h,factor,full,obs,j):
 lw=lambda c:sum((i+h)*n for i,n in enumerate(c))
 R0=K-lw(full);unseen=lw(full)-lw(obs)
 caps=range(factor*R0+1) if j>h else range(factor*(j-1)+1,factor*R0+1)
 for cap in caps:
  label=(obs,unseen,cap)
  low_size=F(j-1,K);high_size=F(j,K) if j>h else min(F(j,K),F(cap,factor*K))
  samples={low_size,high_size};kink=F(cap-1,factor*K)
  if low_size<kink<high_size:samples.add(kink)
  for size in sorted(samples):
   x=K*size-(j-1)
   low=max(F(0),kink) if j>h else max(kink,size)
   yield label,size,x,low

def quadratic_minorant_cases(K,h,factor,full,obs,j):
 """Yield (label,x,S Bernstein degree-one coeff,S² degree-two coeff)."""
 rows=list(minorant_cases(K,h,factor,full,obs,j))
 for label,size,x,low in rows:yield label,x,low,low*low
 if j<=h:
  for left,right in zip(rows,rows[1:]):
   if left[0]==right[0]:
    yield left[0],(left[2]+right[2])/2,(left[3]+right[3])/2,left[3]*right[3]
