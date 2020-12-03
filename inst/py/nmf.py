#!/usr/bin/env python3
import multiprocess
import numpy as np
import nimfa
import itertools
from itertools import repeat

'''Python interface to brunet NMF algorithm

Copyright @ Shixiang Wang <w_shixiang@163.com>
'''

def NMF(V, rank):
  '''Run NMF for a target matrix.
  V: the target matrix to decompose.'''
  nmf = nimfa.Nmf(
    V, 
    seed="random_c", 
    rank=rank,
    n_run=1,
    max_iter=2000,
    update='divergence',
    objective='div')
  nmf_fit = nmf()
  W = nmf_fit.basis()
  K = W.shape[1]
  print('Stop at iteration #: %d' % nmf_fit.summary()["n_iter"])
  return {"W":W, "H":nmf_fit.coef(), "K":K, "KLD":nmf_fit.distance(metric='kl')}

def NMF2(i, V, rank):
  print("Starting NMF run: %d." %i)
  return NMF(V, rank)

def MultiNMF(V, rank, nrun, cores, return_best = True):
  '''Run NMF for a target matrix multiple times'''
  available_cores = multiprocess.cpu_count()
  if available_cores < cores:
    print("Only %d cores available but %d cores requested, reset it to available number." %(available_cores, cores))
    cores = available_cores
  xs = range(1, nrun + 1)
  with multiprocess.Pool(processes=cores) as p:
    v = p.starmap(NMF2, zip(xs, repeat(V), repeat(rank)))

  if return_best:
    print("Compare KLD and return the best.")
    l = len(v)
    kld = v[0]["KLD"]
    idx = 0
    for i in range(1, l):
      k = v[i]["KLD"]
      print("KLD to compare: %f, KLD to be compared: %f" %(k, kld))
      if kld > k:
        print("Current KLD value updated to %f" %k)
        kld = k
        idx = i
    print("Return the run #%d" %(idx+1))
    v = v[idx]

  return v

def NMF3(i, V, rank):
  print("Starting NMF run: %d." %i)
  return {str(i): NMF(V, rank)}

def MultiVNMF(V_tuple, rank, nrun, cores):
  '''Run NMF for a target matrix multiple times'''
  available_cores = multiprocess.cpu_count()
  if available_cores < cores:
    print("Only %d cores available but %d cores requested, reset it to available number." %(available_cores, cores))
    cores = available_cores
  xs = range(1, len(V_tuple) * nrun + 1)
  v = list()
  with multiprocess.Pool(processes=cores) as p:
    v.append(
      p.starmap(
        NMF3, 
        zip(
          xs,
          list(itertools.chain.from_iterable(itertools.repeat(i, nrun) for i in V_tuple)),
          repeat(rank)
          )))
    
  return v
  
# if __name__ == "__main__":
#   V = np.array([[1, 2, 3], [4, 5, 6], [6, 7, 8]])
#   MultiNMF(V, 2, 2, 2)
#   V1 = np.array([[1, 2, 3], [4, 5, 6], [6, 7, 8]])
#   V2 = np.array([[1, 2, 8], [100, 5, 6], [6, 7, 8]])
#   x = MultiVNMF((V1, V2), 2, 2, 2)
#   print(x)
