# Ensure the 'scripts' directory exists
import os
os.makedirs('scripts', exist_ok=True)

# Use a Python multiline string to define the script content
script_content = '''
"""stats.py
Permutation test, Bayesian credible intervals, bootstrap stability
"""
import numpy as np
from scipy.stats import beta
from sklearn.metrics import adjusted_rand_score

def permutation_test(observed_counts, reference_counts, n_perm=10000, seed=42):
    rng = np.random.default_rng(seed)
    obs = np.array(observed_counts)
    ref = np.array(reference_counts)
    obs_prop = obs / obs.sum()
    ref_prop = ref / ref.sum()
    diffs = []
    for _ in range(n_perm):
        perm = rng.permutation(ref_prop)
        diffs.append(np.abs(perm - ref_prop).max())
    p_val = np.mean(np.array(diffs) >= np.abs(obs_prop - ref_prop).max())
    return p_val

def credible_interval(k, n, alpha=0.05):
    lower = beta.ppf(alpha/2, k+1, n-k+1)
    upper = beta.ppf(1-alpha/2, k+1, n-k+1)
    return lower, upper

def bootstrap_cluster_stability(adata, n_bootstraps=20, sample_frac=0.8, leiden_resolution=1.0, seed=42):
    rng = np.random.default_rng(seed)
    base_labels = adata.obs['leiden'].astype(str)
    aris = []
    for i in range(n_bootstraps):
        idx = rng.choice(adata.n_obs, size=int(adata.n_obs * sample_frac), replace=False)
        sub = adata[idx,:].copy()
        import scanpy as sc
        sc.tl.pca(sub, svd_solver='arpack', n_comps=min(30, sub.n_vars))
        sc.pp.neighbors(sub, n_neighbors=10, n_pcs=20)
        sc.tl.leiden(sub, resolution=leiden_resolution)
        common = sub.obs_names
        ari = adjusted_rand_score(base_labels.loc[common], sub.obs['leiden'])
        aris.append(ari)
    return aris
'''

# Write the content to the file using Python's file I/O
with open('scripts/stats.py', 'w') as f:
    f.write(script_content)

print("scripts/stats.py created successfully.")
