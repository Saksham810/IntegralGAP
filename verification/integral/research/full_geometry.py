"""Exact spatial verification for size-affine interior densities."""
from fractions import Fraction as F

def crossing(K,j,z,s,a0,a1,gm,gp,atom=None):
 ans=F(0)
 for r,gr in enumerate(gm):
  if r<K-j:gr=(j-s)*gr+(s-j+1)*gp[r]
  ans+=gr*max(F(0),min(r+1,z,K-s)-max(r,z-s,F(0)))
 if s>0 and (z<=s if atom is None else atom):ans+=(j-s)*a0+(s-j+1)*a1
 return ans

def spatial_checks(K,j,a0,a1,gm,gp,lo,hi):
 def H(l,x,y):return (1-y)*((1-x)*lo[l]+x*hi[l])+y*((1-x)*lo[l+1]+x*hi[l+1])
 def inc(x):return x if j==1 else F(1)
 half=F(1,2)
 for l in range(K):
  pts=[(F(0),F(0)),(F(0),F(1)),(F(1),F(0)),(F(1),F(1))]
  seg=[(pts[0],pts[3])]
  if l==K-j:
   mid=(half,half);pts.append(mid)
   seg=[(pts[0],mid),(mid,pts[3]),(pts[1],mid),(mid,pts[2])]
  seg += [(pts[0],pts[2]),(pts[1],pts[3])]
  for x,y in pts:
   yield H(l,x,y)-crossing(K,j,l+y,j-1+x,a0,a1,gm,gp),inc(x)
  for (x0,y0),(x1,y1) in seg:
   xm,ym=(x0+x1)/2,(y0+y1)/2
   atom=(l+ym<=j-1+xm)
   v0=H(l,x0,y0)-crossing(K,j,l+y0,j-1+x0,a0,a1,gm,gp,atom)
   v1=H(l,x1,y1)-crossing(K,j,l+y1,j-1+x1,a0,a1,gm,gp,atom)
   vm=H(l,xm,ym)-crossing(K,j,l+ym,j-1+xm,a0,a1,gm,gp,atom)
   yield 2*vm-(v0+v1)/2,(inc(x0)+inc(x1))/2
