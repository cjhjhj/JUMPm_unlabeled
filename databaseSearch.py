import sys, os, re, sqlite3, pickle, pandas as pd, numpy as np, multiprocessing as mp, timeit
import runMetFrag


# Pre-calculated (fully-aligned) features
features = pd.read_pickle("IROA_c18_target_full_features.pickle")

# Main part (multi-thread)
if __name__ == '__main__':
    start = timeit.default_timer()
    pool = mp.Pool(mp.cpu_count())
    res = pool.map(runMetFrag.run, [row.to_dict() for idx, row in features.iloc[0:100].iterrows()])
    pool.close()
    print(timeit.default_timer() - start)
    print()
