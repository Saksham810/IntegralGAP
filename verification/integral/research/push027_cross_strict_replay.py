"""Independent exact replay for observed-small-mass reserve certificates.

Uses only the standard library and existing exact geometry/integration lemmas.
No optimizer, floating-point calculations, or repair is used for acceptance.
"""
import gzip,itertools,json,math,sys,time
from fractions import Fraction as F
from decimal import Decimal,localcontext
from exact_geometry import availability_controls
from push027_cross_geometry import spatial_checks
from push027_strict_geometry import labels_and_pairs,packing_values,minorant_cases,quadratic_minorant_cases
from replay_exact import counts,freeze,exp_bounds,integral_weights,require,signed_interval_sum

def replay(path):
 started=time.time()
 with gzip.open(path,'rt') as f:data=json.load(f)
 s=data['summary'];require(s.get('schema')=='finecap-moment-cross2-strict-v1','schema')
 K,h,M,Q=s['K'],s['h'],s['M'],s['denominator'];T=F(s['T']);D=3
 factor=s['factor'];tight_B=s['tight_B'];require(s.get('variance_bound')=='exact' and tight_B in (0,1),'variance model');require(type(factor) is int and factor>=1,'capacity factor')
 require(0<h<K and M>0 and Q>0 and 1<=T<102,'parameters')
 require(s.get('target_degree',3)==3,'degree')
 availability=s.get('availability','subdiv');require(availability in ('linear','subdiv'),'availability')
 cs=counts(K,h)
 def lw(c):return sum((i+h)*n for i,n in enumerate(c))
 basepairs=[(c,j) for c in cs for j in range(1,K+1) if (c[j-h-1]>0 if j>h else lw(c)+j-1<=K-1)]
 bases=sorted({(obs,lw(full)-lw(obs)) for full in cs for obs in itertools.product(*[range(n+1) for n in full])})
 labels,pairs=labels_and_pairs(bases,K,h,factor)
 def arr(name,keys,signed=False):
  rows=data[name];a={freeze(k):v for k,v in rows}
  require(len(a)==len(rows) and set(a)==set(keys),name+': indices')
  require(all(type(v) is int and (signed or v>=0) for v in a.values()),name+': entries')
  return a
 P=arr('P',[(b,l) for b in range(D*M+1) for l in range(K+1)])
 ak=[(b,e,label,j) for b in range(M) for e in range(D) for label,j in pairs]
 A=arr('A',ak);AL=arr('AL',ak)
 G=arr('G',[k+(r,) for k in ak for r in range(K-k[-1]+1)])
 hk=[k+(l,) for k in ak for l in range(K+1)]
 L=arr('L',hk);U=arr('U',hk)
 require(all(v<=Q for v in P.values()),'tail range')
 require(all(P[b,l]>=P[b,l+1] for b in range(D*M+1) for l in range(K)),'load monotonicity')
 require(all(P[b,l]<=P[b+1,l] for b in range(D*M) for l in range(K+1)),'time monotonicity')
 require(all(AL[k]==0 for k in ak if k[-1]==1),'zero-size atom')
 require(all(L[k]==0 for k in hk if k[-2]==1),'zero-size crossing')
 availability_checks=0
 for b in range(M):
  for label,j in pairs:
   for r in range(-2,K-j+1):
    av=[Q-P[D*b+k,0] if r<0 else P[D*b+k,r]-P[D*b+k,r+1] for k in range(D+1)]
    v=[AL[b,e,label,j] if r==-2 else A[b,e,label,j] if r==-1 else G[b,e,label,j,r] for e in range(D)]
    if availability=='subdiv':
     require(all(vv>=0 for vv in availability_controls(av,v)),'availability subdivision')
     availability_checks+=7
    else:
     require(v[0]<=av[0] and v[-1]<=av[-1],'availability endpoints')
     require(all(k*v[k-1]+(D-k)*v[k]<=D*av[k] for k in range(1,D)),'availability controls')
     availability_checks+=4
 crossing_checks=0
 for key in ak:
  b,e,label,j=key
  gg=[G[key+(r,)] for r in range(K-j+1)]
  ll=[L[key+(l,)] for l in range(K+1)];uu=[U[key+(l,)] for l in range(K+1)]
  for margin,_ in spatial_checks(K,j,AL[key],A[key],gg,ll,uu):
   require(margin>=0,'spatial crossing');crossing_checks+=1
 R=[F(v) for v in data['correction_slopes']]
 require(len(R)==M and min(R)>=0,'reserve slopes')
 Hstar=T/M*sum(R);alpha=1/(1+Hstar)
 require(F(s['correction_Hstar'])==Hstar and F(s['alpha'])==alpha,'reserve normalization')
 packing_checks=0
 for b in range(M):
  for e in range(D):
   for label in labels:
    for l in range(K+1):
     ll={j:L[b,e,label,j,l] for lab,j in pairs if lab==label}
     uu={j:U[b,e,label,j,l] for lab,j in pairs if lab==label}
     vals=packing_values(K,h,factor,label,ll,uu)
     allowance=F(D*M,T*Q)*(P[D*b+e+1,l]-P[D*b+e,l])+R[b]
     require(F(max(vals),Q)<=allowance,'packing');packing_checks+=len(vals)
 witnesskeys=[]
 for full,j in basepairs:
  m=list(full)
  if j>h:m[j-h-1]-=1
  for sub in itertools.product(*[range(n+1) for n in m]):
   obs=list(sub)
   if j>h:obs[j-h-1]+=1
   obs=tuple(obs)
   witnesskeys += [(full,obs,j,b,e) for b in range(M) for e in range(D)]
 minor=arr('minor',[k+(ep,) for k in witnesskeys for ep in [0,1]],True)
 slope=arr('slope',witnesskeys,True)
 curve=arr('curve',witnesskeys,True);require(all(v<=0 for v in curve.values()),'positive curve');delta=F(h,K)
 def Jmax(C):
  k=C//delta;return k*delta*delta+(C-k*delta)**2
 require(all(v<=0 for v in slope.values()),'positive slope')
 minorant_checks=0
 for key in witnesskeys:
  full,obs,j,b,e=key;R0=K-lw(full);unseen=lw(full)-lw(obs);v=slope[key]
  for label,x,low,low2 in quadratic_minorant_cases(K,h,factor,full,obs,j):
   q=(1-x)*AL[b,e,label,j]+x*A[b,e,label,j]+sum(G[b,e,label,j,r] for r in range(K-j))+(1-x)*G[b,e,label,j,K-j]
   u=(1-x)*minor[key+(0,)]+x*minor[key+(1,)]
   require(u+v*low+curve[key]*low2<=q,'minorant');minorant_checks+=1
 maxcompanions=(K-1)//h
 weights={(d,n,b,e):integral_weights(d,n,b,e,T,M,D) for d in range(maxcompanions+3) for n in range(d+1) for b in range(M) for e in range(D)}
 coefficients=[];worst=None
 for full,j in basepairs:
  m=list(full)
  if j>h:m[j-h-1]-=1
  d=sum(m);Bmax=F(K-lw(full),K)
  for ep in [0,1]:
   size=F(j-1+ep,K);lo=hi=F(0)
   for sub in itertools.product(*[range(n+1) for n in m]):
    mult=math.prod(math.comb(n,t) for n,t in zip(m,sub));n=sum(sub);obs=list(sub)
    if j>h:obs[j-h-1]+=1
    obs=tuple(obs)
    for b in range(M):
     for e in range(D):
      key=(full,obs,j,b,e);u=F(minor[key+(ep,)],Q);v=F(slope[key],Q);cc=F(curve[key],Q)
      wp=weights[d+1,n+1,b,e];w=weights[d,n,b,e];wpp=weights[d+2,n+2,b,e]
      if j>h:
       bt=Bmax-(F(ep,K) if tight_B else 0)
       terms=[(u,w),(v*bt+cc*Jmax(bt),wp),(cc*(bt*bt-Jmax(bt)),wpp)]
      else:
       rest=Bmax-size
       terms=[(u+v*size+cc*size*size,w),(v*rest+cc*(2*size*rest+Jmax(rest)),wp),(cc*(rest*rest-Jmax(rest)),wpp)]
      lv,hv=signed_interval_sum(terms);lo+=mult*lv;hi+=mult*hv
   val=alpha*lo;coefficients.append(val)
   if worst is None or val<worst[0]:worst=(val,full,j,ep)
 gamma=min(coefficients);claimed=F(s['limiting_coefficient'])
 require(0<=claimed<=gamma,'service guarantee')
 require(claimed<=1-exp_bounds(-T)[1],'zero-size guarantee')
 C=max(F(D*(D-1)*M*M,T*T*Q)*abs(P[D*b+v+2,l]-2*P[D*b+v+1,l]+P[D*b+v,l]) for b in range(M) for v in range(D-1) for l in range(K+1))
 reported_C=F(s['second_derivative_bound']);require(reported_C>=C,'curvature')
 H=max([F(1)]+[F(abs(minor[k+(ep,)])+abs(slope[k])+abs(curve[k]),Q) for k in witnesskeys for ep in [0,1]])
 eu=exp_bounds(T)[1];E=eu.numerator//eu.denominator+1;d=maxcompanions
 finite=H*(4*(d+2)**2+(20*E+20)*M+20)+F(E,2)*T*reported_C
 require(F(s['finite_loss_constant'])>=finite,'finite correction')
 with localcontext() as ctx:
  ctx.prec=35;decimal=str(Decimal(claimed.numerator)/Decimal(claimed.denominator))
 result={'status':'PASS','certificate':path,'counts':len(cs),'labels':len(labels),'rules':len(pairs),'availability_checks':availability_checks,'crossing_checks':crossing_checks,'packing_checks':packing_checks,'minorant_checks':minorant_checks,'service_coefficients':len(coefficients),'limiting_lower_bound':str(claimed),'decimal':decimal,'worst_pair':[worst[1],worst[2],worst[3]],'witness_bound':str(H),'finite_loss_constant':s['finite_loss_constant'],'elapsed':round(time.time()-started,3)}
 print(json.dumps(result),flush=True)
 return result

if __name__=='__main__':
 for path in sys.argv[1:]:replay(path)
