"""Exact spatial inequalities for truncated-density, affine-atom rules.
Only fractions and integer arithmetic; independent of the optimizer.
"""
from fractions import Fraction as F

def crossing(K,j,z,s,a0,a1,g):
    ans=F(0)
    for r,gr in enumerate(g):
        ans+=gr*max(F(0),min(r+1,z,K-s)-max(r,z-s,F(0)))
    if s>0 and z<=s:
        ans+=(j-s)*a0+(s-j+1)*a1
    return ans

def spatial_checks(K,j,a0,a1,g,lo,hi):
    """Yield (envelope minus bound, increase under a uniform repair).
    For j=1 only upper-size controls are increased, preserving H(z,0)=0.
    Other classes increase all envelope controls equally.
    """
    def H(l,x,y):return (1-y)*((1-x)*lo[l]+x*hi[l])+y*((1-x)*lo[l+1]+x*hi[l+1])
    def inc(x):return x if j==1 else F(1)
    half=F(1,2)
    for l in range(K):
        pts=[(F(0),F(0)),(F(0),F(1)),(F(1),F(0)),(F(1),F(1))]
        segments=[((F(0),F(0)),(F(1),F(1)))]
        if l==K-j:
            mid=(half,half);pts.append(mid)
            segments=[(pts[0],mid),(mid,pts[3]),(pts[1],mid),(mid,pts[2])]
        subdivision=2
        segments=[(tuple(p0[k]+(p1[k]-p0[k])*F(i,subdivision) for k in [0,1]),tuple(p0[k]+(p1[k]-p0[k])*F(i+1,subdivision) for k in [0,1])) for p0,p1 in segments for i in range(subdivision)]
        pts=sorted(set(pts+[p for segment in segments for p in segment]))
        for x,y in pts:
            yield H(l,x,y)-crossing(K,j,l+y,j-1+x,a0,a1,g),inc(x)
        for (x0,y0),(x1,y1) in segments:
            h0,h1=H(l,x0,y0),H(l,x1,y1)
            hm=H(l,(x0+x1)/2,(y0+y1)/2)
            b0=crossing(K,j,l+y0,j-1+x0,a0,a1,g)
            b1=crossing(K,j,l+y1,j-1+x1,a0,a1,g)
            yield 2*hm-(h0+h1+b0+b1)/2,(inc(x0)+inc(x1))/2


def availability_controls(av,v):
    """Bernstein coefficients after subdividing the cubic difference at 1/2.
    Integer weights need no division when testing nonnegativity.
    """
    d=[av[0]-v[0], av[1]-(v[0]+2*v[1])/F(3),
       av[2]-(2*v[1]+v[2])/F(3), av[3]-v[2]]
    for w in [(1,0,0,0),(1,1,0,0),(1,2,1,0),(1,3,3,1),
              (0,1,2,1),(0,0,1,1),(0,0,0,1)]:
        yield sum(wi*di for wi,di in zip(w,d))
