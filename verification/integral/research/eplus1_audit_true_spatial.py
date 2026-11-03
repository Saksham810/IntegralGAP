"""Exact actual-boundary audit, retaining one-sided traces at atom jumps.

On each size/load grid square, the density crossing is piecewise affine.
Its only interior break lines are y=x, and (when l=K-j) x+y=1.
The atom break z=s is y=x when l=j-1; in adjacent cells it only
meets a corner. Therefore the listed open segments contain no break.
Their true H-B quadratic is recovered from three INTERIOR samples;
closed endpoints are checked separately. This avoids interpolating an atom
present only at an endpoint into the open segment.
"""
import gzip,json,sys,time,hashlib,os
from fractions import Fraction as F
from pathlib import Path

def freeze(x):return tuple(map(freeze,x)) if isinstance(x,list) else x

def minimum_quadratic(d0,d1,d2):
 ans=min(d0,d2);den=d0-2*d1+d2
 if den>0:
  t=(d0-d1)/den
  if 0<t<1:ans=min(ans,(1-t)**2*d0+2*t*(1-t)*d1+t*t*d2)
 return ans

def audit(path):
 start=time.time();data=json.load(gzip.open(path,'rt'));K=data['summary']['K'];Q=data['summary']['denominator']
 a={n:{freeze(k):F(v) for k,v in data[n]} for n in ('A','AL','G','L','U')}
 count=interior_checks=jump_endpoints=0;worst=None
 for key,aa in a['A'].items():
  j=key[-1]
  def B(z,s):
   ans=sum(a['G'][key+(r,)]*max(F(0),min(r+1,z,K-s)-max(r,z-s,0)) for r in range(K-j+1))
   if s>0 and z<=s:ans+=(j-s)*a['AL'][key]+(s-j+1)*aa
   return ans
  def H(l,x,y):return (1-y)*((1-x)*a['L'][key+(l,)]+x*a['U'][key+(l,)])+y*((1-x)*a['L'][key+(l+1,)]+x*a['U'][key+(l+1,)])
  for l in range(K):
   pts=[(F(0),F(0)),(F(0),F(1)),(F(1),F(0)),(F(1),F(1))]
   segments=[(pts[0],pts[3])]
   if l==K-j:
    m=(F(1,2),F(1,2));pts.append(m);segments=[(pts[0],m),(m,pts[3]),(pts[1],m),(m,pts[2])]
   for x,y in pts:
    v=H(l,x,y)-B(l+y,j-1+x)
    assert v>=0,('closed corner',key,l,x,y,v)
    count+=1
   for p,q in segments:
    x0,y0=p;x1,y1=q
    def gap(t):
     x=x0+(x1-x0)*t;y=y0+(y1-y0)*t
     return H(l,x,y)-B(l+y,j-1+x)
    f1,f2,f3=(gap(t) for t in (F(1,4),F(1,2),F(3,4)))
    d0=3*f1-3*f2+f3;d2=f1-3*f2+3*f3;d1=5*f2-2*f1-2*f3
    t=F(1,3)
    assert gap(t)==(1-t)**2*d0+2*t*(1-t)*d1+t*t*d2,('open trace identity',key,l,p,q)
    interior_checks+=1
    e0,e2=gap(F(0)),gap(F(1))
    assert e0>=0 and e2>=0,('closed segment endpoint',key,l,p,q,e0,e2)
    if d0!=e0:jump_endpoints+=1
    if d2!=e2:jump_endpoints+=1
    v=minimum_quadratic(d0,d1,d2)
    assert v>=0,('one-sided quadratic',key,l,p,q,v)
    if worst is None or v<worst:worst=v
    count+=1
 result={'status':'PASS','certificate':path,'certificate_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),'exact_corner_and_quadratic_checks':count,'additional_interior_polynomial_identity_checks':interior_checks,'jump_endpoint_occurrences':jump_endpoints,'minimum_margin_certificate_units':str(worst),'minimum_margin_unscaled':str(worst/Q),'method':'Closed boundary values plus one-sided open-segment quadratics reconstructed from three interior samples; exact interior identity crosscheck.','elapsed':round(time.time()-start,3)}
 print(json.dumps(result),flush=True)
 return result

if __name__=='__main__':
 path=sys.argv[1];out=Path(sys.argv[2]) if len(sys.argv)>2 else Path(path.replace('_cert.json.gz','_spatial_audit.json'))
 result=audit(path);tmp=out.with_suffix(out.suffix+'.tmp')
 with tmp.open('w') as f:json.dump(result,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(tmp,out)
