from lmfdb import db

R.<x> = PolynomialRing(QQ)

ell = 5
fil = open('../data/ellipticcurvesinfo5.txt','a')

query = {
    'sha_primes': { '$contains': [ell] },
    'bad_primes' :  { '$not': { '$contains': ell } },
    'rank' : 0,
    'torsion_primes': { '$not': { '$contains': ell } }
}
ec = list(db.ec_curvedata.search(query))
len(ec)

for x in ec:
  lbl = x['lmfdb_label']
  isolbl = x['lmfdb_iso']
  E_tamagawa = list(db.ec_mwbsd.search({'lmfdb_label': lbl}, projection=['tamagawa_product']))[0]['tamagawa_product']
  if E_tamagawa%ell != 0:
    aps = list(db.ec_classdata.search({'lmfdb_iso': isolbl}, projection=['aplist']))[0]['aplist']
    fil.write(':'.join((str(x['ainvs']), str(x['conductor']), lbl, str(aps), '\n')))

fil.close()
