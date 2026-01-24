from lmfdb import db

R.<x> = PolynomialRing(QQ)

L = list(db.g2c_curves.search({'analytic_sha': {'$ne': 1},'is_simple_geom': False}))
len(L)
Es = [list(db.g2c_endomorphisms.search({'label':x['label']}, projection=['spl_facs_labels'])) for x in L]
len(Es)
{len(x[0]['spl_facs_labels']) for x in Es if len(x[0]) == 1}

Es_sha = []
goodinds = []
for i in range(len(Es)):
  if len(Es[i][0]) == 1:
    lbls = Es[i][0]['spl_facs_labels']
    if len(lbls) == 2:
      E1 = list(db.ec_curvedata.search({'lmfdb_label':lbls[0]},projection=['sha']))
      if len(E1) == 0:
        continue
      E2 = list(db.ec_curvedata.search({'lmfdb_label':lbls[1]},projection=['sha']))
      if len(E1) == 0:
        continue
      sha1 = E1[0]['sha']
      sha2 = E2[0]['sha']
      if (sha1 != 1 and sha2 == 1) or (sha1 == 1 and sha2 != 1):
        Es_sha.append([sha1, sha2])
        print(L[i]['label'])
        print(L[i]['analytic_sha'], {sha1,sha2})
        goodinds.append(i)

for i in range(len(goodinds)):
  ii = goodinds[i]
  print(L[ii]['label'], L[ii]['analytic_sha'], Es_sha[i])
'''
8730.a.235710.1 4 [4, 1]
10800.c.691200.1 16 [4, 1]
14400.e.432000.1 4 [1, 4]
15680.b.250880.1 16 [1, 4]
72000.b.72000.1 16 [1, 4]
75660.a.605280.1 4 [4, 1]
83200.d.832000.1 4 [1, 4]
95160.a.190320.1 4 [4, 1]
114240.a.114240.1 32 [1, 4]
177660.a.355320.1 32 [1, 4]
184320.b.552960.1 4 [1, 4]
187200.a.748800.1 16 [1, 4]
190320.a.190320.1 16 [4, 1]
196800.a.393600.1 4 [1, 4]
330624.b.330624.1 36 [9, 1]
439280.a.439280.1 16 [1, 4]
460362.a.460362.1 4 [1, 4]
483552.a.483552.1 16 [4, 1]
540800.a.540800.1 16 [4, 1]
540800.b.540800.1 16 [1, 4]
547344.a.547344.1 36 [9, 1]
549250.a.549250.1 16 [4, 1]
702720.a.702720.1 4 [4, 1]
723520.b.723520.1 32 [1, 4]
735488.c.735488.1 4 [4, 1]
735488.d.735488.1 16 [1, 4]
839680.a.839680.1 4 [4, 1]
'''

'''
Well... the map on Sha induced from 0 -> E1 -> J -> E2 -> 0, i.e.,
Sha(E1) -> Sha(J) -> Sha(E2) is not exact even in the middle.
eventhough H1(E1) -> H1(J) -> H1(E2) is exact in the middle.
'''