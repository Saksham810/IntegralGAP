"""Replay an unseen-large-reserve certificate using only Python's standard library.

Usage: python -S push_reserve_cert/replay_reserved.py path/to/certificate.json.gz
No numerical optimizer, NumPy, search data, or rational repair is used.
"""
import gzip, itertools, json, math, sys, time
from fractions import Fraction as F
from exact_geometry import spatial_checks, availability_controls
from full_geometry import spatial_checks as full_checks
from decimal import Decimal, localcontext

def counts(K,h):
 out=[]
 def rec(j,r,c):
  if j>K: out.append(tuple(c)); return
  for n in range(r//(j-1)+1): rec(j+1,r-(j-1)*n,c+[n])
 rec(h+1,K-1,[])
 return out

def freeze(x): return tuple(freeze(y) for y in x) if isinstance(x,list) else x
def exp_bounds(x,terms=100):
 if x<0:
  lo,hi=exp_bounds(-x,terms); return 1/hi,1/lo
 term=F(1); lo=term
 for j in range(1,terms+1): term=term*x/j; lo+=term
 tail=term*x/(terms+1)/(1-x/(terms+2))
 hi=lo+tail;grid=10**80
 return F(lo.numerator*grid//lo.denominator,grid),F((hi.numerator*grid+hi.denominator-1)//hi.denominator,grid)
def signed_interval_sum(terms):
 lo=hi=F(0)
 for v,(l,r) in terms:
  lo+=v*(l if v>=0 else r); hi+=v*(r if v>=0 else l)
 return lo,hi
def power(i,n): return i[0]**n,i[1]**n
def integral_weights(d,n,b,e,T,M,target_degree=2):
 delta=T/M
 l=exp_bounds(-T+b*delta); r=exp_bounds(-T+(b+1)*delta)
 terms=[]
 for j in range(d-n+1):
  m=n+j+1; den=delta*m*m
  J0=signed_interval_sum([(F(1,m),power(r,m)),(F(-1,m),power(l,m))])
  J1=signed_interval_sum([((m*delta-1)/den,power(r,m)),(1/den,power(l,m))])
  J2=signed_interval_sum([((delta*delta/m-2*delta/(m*m)+F(2,m**3))/(delta*delta),power(r,m)),(-F(2,m**3)/(delta*delta),power(l,m))])
  Js=[J0,J1,J2];rdeg=target_degree-1
  q=signed_interval_sum([(math.comb(rdeg,e)*(-1)**v*math.comb(rdeg-e,v),Js[e+v]) for v in range(rdeg-e+1)])
  terms.append((F((-1)**j*math.comb(d-n,j),1),q))
 lo,hi=signed_interval_sum(terms)
 grid=10**30
 return max(F(0),F(lo.numerator*grid//lo.denominator,grid)),F((hi.numerator*grid+hi.denominator-1)//hi.denominator,grid)

def require(condition,message):
 if not condition: raise ValueError(message)

def read_array(data,name,keys):
 rows=data[name];a={freeze(k):v for k,v in rows}
 require(len(a)==len(rows) and set(a)==set(keys),name+': missing, extra, or repeated index')
 require(all(type(v) is int and v>=0 for v in a.values()),name+': entries must be nonnegative integers')
 return a

def replay(path):
 started=time.time()
 with gzip.open(path,'rt') as f:data=json.load(f)
 s=data['summary'];require(s.get('schema')=='exact-size-reserve-v1','wrong certificate schema');K,h,M=s['K'],s['h'],s['M'];T=F(s['T']);Q=s['denominator']
 full=s.get('full_density',False);D=s.get('target_degree',2);mode=s['mode'];availability=s.get('availability','linear')
 require(0<h<K and M>0 and Q>0 and 0<T<102,'invalid basic parameters')
 require(D in (2,3) and mode in ('linear','conic'),'unsupported polynomial degree or crossing rule')
 require(availability in ('linear','conic','subdiv') and (availability in ('linear','subdiv') or D==2),'unsupported availability rule')
 cs=counts(K,h)
 def lower_weight(c):return sum((i+h)*n for i,n in enumerate(c))
 basepairs=[(c,j) for c in cs for j in range(1,K+1)
            if (c[j-h-1]>0 if j>h else lower_weight(c)+j-1<=K-1)]
 labels=sorted({(observed,lower_weight(full)-lower_weight(observed))
                for full in cs
                for observed in itertools.product(*[range(n+1) for n in full])})
 require(all(0<=u and lower_weight(c)+u<=K-1 for c,u in labels),'invalid regenerated reserve label')
 require(s.get('labels')==len(labels),'incorrect label count')
 pairs=[(label,j) for label in labels for j in range(1,K+1)
        if (label[0][j-h-1]>0 if j>h else lower_weight(label[0])+label[1]+j-1<=K-1)]
 P=read_array(data,'P',[(b,l) for b in range(D*M+1) for l in range(K+1)])
 ak=[(b,e,c,j) for b in range(M) for e in range(D) for c,j in pairs]
 A=read_array(data,'A',ak);AL=read_array(data,'AL',ak)
 G=read_array(data,'G',[k+(r,) for k in ak for r in range(K-k[-1]+1)])
 GU=read_array(data,'GU',[k+(r,) for k in ak for r in range(K-k[-1])]) if full else {}
 hk=[k+(l,) for k in ak for l in range(K+1)]
 L=read_array(data,'L',hk);U=read_array(data,'U',hk)
 require(all(v<=Q for v in P.values()),'target tail exceeds one')
 require(all(P[b,l]>=P[b,l+1] for b in range(D*M+1) for l in range(K)),'target fails monotonicity in load')
 require(all(P[b,l]<=P[b+1,l] for b in range(D*M) for l in range(K+1)),'target fails monotonicity in time')
 require(all(AL[k]==0 for k in ak if k[-1]==1),'tiny zero-size endpoint accepts from zero')
 require(all(L[k]==0 for k in hk if k[-2]==1),'tiny class has nonzero zero-size crossing bound')
 availability_checks=0
 for b in range(M):
  for c,j in pairs:
   for r in range(-2,2*(K-j)+1 if full else K-j+1):
    rr=r if r<=K-j else r-(K-j+1)
    av=[Q-P[D*b+k,0] if r<0 else P[D*b+k,rr]-P[D*b+k,rr+1] for k in range(D+1)]
    v=[AL[b,e,c,j] if r==-2 else A[b,e,c,j] if r==-1 else G[b,e,c,j,r] if r<=K-j else GU[b,e,c,j,rr] for e in range(D)]
    require(v[0]<=av[0] and v[-1]<=av[-1],'endpoint availability fails')
    if availability=='subdiv':
     require(D==3 and all(vv>=0 for vv in availability_controls(av,v)),'subdivided availability fails')
    elif availability=='linear':
     require(all(k*v[k-1]+(D-k)*v[k]<=D*av[k] for k in range(1,D)),'control availability fails')
    else:
     a=av[0]-v[0];cc=av[2]-v[1];dd=sum(v)-2*av[1]
     require(dd<=0 or dd*dd<=4*a*cc,'exact quadratic availability fails')
    availability_checks+=(7 if availability=='subdiv' else D+1)
 crossing_checks=0
 for b,e,c,j in ak:
  key=(b,e,c,j)
  gg=[G[key+(r,)] for r in range(K-j+1)]
  ll=[L[key+(l,)] for l in range(K+1)]
  uu=[U[key+(l,)] for l in range(K+1)]
  checks=full_checks(K,j,AL[key],A[key],gg,[GU[key+(r,)] for r in range(K-j)],ll,uu) if full else spatial_checks(K,j,AL[key],A[key],gg,ll,uu)
  for margin,_ in checks:
   require(margin>=0,'spatial crossing inequality fails')
   crossing_checks+=1
 R=[F(v) for v in data['correction_slopes']]
 require(len(R)==M and min(R)>=0,'invalid correction slopes')
 Hstar=T/M*sum(R);alpha=1/(1+Hstar)
 require(Hstar==F(s['correction_Hstar']) and alpha==F(s['alpha']),'incorrect normalization correction')
 packing_checks=0
 for b in range(M):
  for e in range(D):
   for label in labels:
    c,u=label
    R0=K-lower_weight(c)-u;Js=range(1,min(h,R0)+1)
    for l in range(K+1):
     V=[0]*R0
     for w in range(1,R0):
      opts=[V[w-1]]
      for j in Js:
       if j<=w:opts.append(U[b,e,label,j,l]+V[w-j])
       if 1<j<=w+1:opts.append(L[b,e,label,j,l]+V[w-j+1])
      V[w]=max(opts)
     def ad(d):return sum((n-d[i])*L[b,e,label,i+h+1,l]+d[i]*U[b,e,label,i+h+1,l] for i,n in enumerate(c) if n)
     sums=[]
     for d in itertools.product(*[range(n+1) for n in c]):
      used=sum(d)
      if 0<used<=R0:sums.append(ad(d)+V[R0-used])
     base=ad((0,)*len(c));sums.append(base+V[R0-1])
     sums += [base+U[b,e,label,j,l]+V[R0-j] for j in Js]
     allowance=F(D*M,T*Q)*(P[D*b+e+1,l]-P[D*b+e,l])+R[b]
     require(F(max(sums),Q)<=allowance,'endpoint packing bound fails')
     packing_checks+=len(sums)
 weights={(d,n,b,e):integral_weights(d,n,b,e,T,M,D)[0] for d in range((K-1)//h+1)
          for n in range(d+1) for b in range(M) for e in range(D)}
 coeffs=[];pair_set=set(pairs)
 for c,j in basepairs:
  for endpoint in [0,1]:
   m=list(c)
   if j>h:m[j-h-1]-=1
   d=sum(m);value=F(0)
   for sub in itertools.product(*[range(n+1) for n in m]):
    mult=math.prod(math.comb(n,t) for n,t in zip(m,sub));n=sum(sub);obs=list(sub)
    if j>h:obs[j-h-1]+=1
    obs=tuple(obs);unseen=lower_weight(c)-lower_weight(obs);obs=(obs,unseen)
    require((obs,j) in pair_set,'missing coefficient reserve rule')
    for b in range(M):
     for e in range(D):
      q=F((AL if endpoint==0 else A)[b,e,obs,j]+sum((GU if full and endpoint==1 else G)[b,e,obs,j,r] for r in range(K-j+1-endpoint)),Q)
      value+=mult*weights[d,n,b,e]*q
   coeffs.append(alpha*value)
 gamma=min(coeffs);claimed=F(s['limiting_coefficient'])
 require(0<=claimed<=gamma,'claimed limiting coefficient exceeds verified bound')
 require(claimed<=1-exp_bounds(-T)[1],'zero-size coefficient is inadequate')
 C=max(F(D*(D-1)*M*M,T*T*Q)*abs(P[D*b+v+2,l]-2*P[D*b+v+1,l]+P[D*b+v,l]) for b in range(M) for v in range(D-1) for l in range(K+1))
 reported_C=F(s['second_derivative_bound'])
 require(reported_C>=C,'second derivative bound is insufficient')
 eu=exp_bounds(T)[1];E=eu.numerator//eu.denominator+1;d=(K-1)//h;S=4 if availability=='subdiv' else 1 if D==2 else 3
 require(T>=1,"finite bound assumes T>=1")
 finite=d*(d+2)+(20*E+20)*M+10+F(E,2)*T*reported_C
 require(F(s['finite_loss_constant'])>=finite,'finite-horizon constant is insufficient')
 with localcontext() as ctx:
  ctx.prec=30;decimal=str(Decimal(claimed.numerator)/Decimal(claimed.denominator))
 result={'status':'PASS','certificate':path,'counts':len(cs),'labels':len(labels),'menu_rules':len(pairs),'value_coefficients':2*len(basepairs),
         'availability_checks':availability_checks,'crossing_checks':crossing_checks,
         'packing_checks':packing_checks,'limiting_lower_bound':str(claimed),
         'decimal_for_display':decimal,'finite_loss_constant':s['finite_loss_constant'],
         'elapsed':round(time.time()-started,3)}
 print(json.dumps(result))
 return result

if __name__=='__main__':
 if len(sys.argv)<2: raise SystemExit('Usage: python -S replay_reserved.py certificate.json.gz [...]')
 for filename in sys.argv[1:]:replay(filename)
