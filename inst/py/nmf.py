import multiprocessing
import numpy as np
import nimfa
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

def MultiNMF(V, rank, nrun, cores):
  '''Run NMF for a target matrix multiple times'''
  available_cores = multiprocessing.cpu_count()
  if available_cores < cores:
    print("Only %d cores available but %d cores requested, reset it to available number." %(available_cores, cores))
    cores = available_cores
  xs = range(nrun)
  v = list()
  with multiprocessing.Pool(processes=cores) as p:
    v.append(p.starmap(NMF2, zip(xs, repeat(V), repeat(rank))))
    
  return tuple(v)
  
