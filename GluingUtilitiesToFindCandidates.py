# cd into gluing_utilities/src/

'''
Attach("fast_gluing.m");
P<x> := PolynomialRing(QQ);
E := EllipticCurve(x^3 - x^2 - 191986*x + 32389771,P!1);
C := HyperellipticCurveOfGenus(2,[-x^5 + x^4 + x - 1,x^3+1]);
X3s := AllGeometricGluingsCCEfficient(E,C,QQ,3);
X3s;
'''

load("compatible_curves.py")
# load("gluing_wrapper.py")
# cd "~/lmfdb"
from lmfdb import db

R.<x> = PolynomialRing(QQ)
ell = 3

#g2c = list(db.g2c_curves.search({'non_maximal_primes':{'$contains': [ell]}, 'bad_primes':{'not $contains': [ell]}}))
#g2c = list(db.g2c_curves.search({'non_maximal_primes':{'$contains': [ell]}, 'cond':{'$mod ell':1}}))
g2c = list(db.g2c_curves.search({'non_maximal_primes':{'$contains': [ell]}}))
len(g2c)
g2c = [x for x in g2c if not ell in x['bad_primes']]
len(g2c)
g2c = [x for x in g2c if not x['tamagawa_product'] % ell != 0]
len(g2c)

candidate_C = []
candidate_C_cond = []
candidate_C_ellrk = []
for C in g2c:
  rkC = C['analytic_rank']
  torsC = C['torsion_subgroup']
  # print(torsC)
  # print([int(ord) for ord in torsC[1:-1].split(",")])
  elltorsCrk = len([int(ord) for ord in torsC[1:-1].split(",") if len(ord) > 0 and int(ord) % ell == 0])
  n = rkC + elltorsCrk
#  if n > 0:
  if n > 0 and elltorsCrk == 0: #TODO Change n to psi-version of weak Mordell-Weil rank. don't set elltorsCrk to 0 necessarily.
    candidate_C.append(C['eqn'])
    candidate_C_cond.append(C['cond'])
    candidate_C_ellrk.append(n)

len(candidate_C)
len(candidate_C_cond)
len(candidate_C_ellrk)

# magma.attach("fast_gluing.m")

'''
goodC = []
goodC_cond = []
goodC_ellrk = []
goodE = []
goodE_cond = []
goodE_rk = []
'''

fil = open('curvesinfo.txt','a')
# foundone = False
for i in range(19,len(candidate_C)):
  print(i)
  eqn = candidate_C[i]
  eqn_split = eqn.split("]")
  fcoeffs = eqn_split[0][2:]
  if len(fcoeffs) == 0:
    f = R(0)
  else:
    f = R([int(x) for x in fcoeffs.split(",")])
  hcoeffs = eqn_split[1][2:]
  if len(hcoeffs) == 0:
    h = R(0)
  else:
    h = R([int(x) for x in hcoeffs.split(",")])

  C = HyperellipticCurve(f,h)
  N = candidate_C_cond[i]
  L = compatible_curves(C,N,ell)
  L_cond = [x.conductor() for x in L[0]]
  L_tamagawaprod = [x.tamagawa_product() for x in L[0]]
  goodinds = [j for j in range(len(L[0])) if L_cond[j] % ell != 0 and L_tamagawaprod[j] % ell != 0]
  print(goodinds)
  if len(goodinds) > 0:
    candidate_E_conds = [L_cond[j] for j in goodinds]
    candidate_E = [L[0][j] for j in goodinds]
    candidate_E_ainvs = [L[0][j].ainvs() for j in goodinds]
    candidate_E_ranks = [x.analytic_rank() for x in candidate_E]
    print(C, N, candidate_C_ellrk[i])
    print(candidate_E_conds, candidate_E, candidate_E_ranks)
    fil.write(':'.join((str(f), str(h), str(N), str(candidate_C_ellrk[i]), str(candidate_E_conds), str(candidate_E_ainvs), str(candidate_E_ranks), '\n')))
#    for k in range(len(candidate_E)):
#      glue = magma.function_call("AllGeometricGluingsCCEfficient",[candidate_E[k],C,Rationals(),ell])
#      if len(glue) > 0:
#        print(candidate_E[k])
#        foundone = True
#        break
#    else:
#      continue
#    break

fil.close()


'''
    goodC.append(C)
    goodC_cond.append(N)
    goodC_ellrk.append(candidate_C_ellrk[i])
    goodE.append(candidate_E)
    goodE_cond.append(candidate_E_conds)
    goodE_rk.append(candidate_E_ranks)

len(goodC)
'''

'''
InterfaceError: connection already closed
sage: print(i)
18
sage: print(C)
Hyperelliptic Curve over Rational Field defined by y^2 + (x^3 + x^2 + 1)*y = -5*x^5 + 9*x^4 - 9*x^2 - 3*x + 6
sage: print(N)
51076
'''
