#!/usr/bin/env python
"""03. Consensus NMF subtype discovery in RB+ samples.

Pick k by cophenetic coefficient + dispersion; output per-sample labels.
"""
import argparse, os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, cophenet

ap = argparse.ArgumentParser()
ap.add_argument("--counts", required=True); ap.add_argument("--meta", required=True)
ap.add_argument("--out_labels", required=True); ap.add_argument("--out_coph", required=True)
ap.add_argument("--out_consensus", required=True)
ap.add_argument("--k_min", type=int, default=2); ap.add_argument("--k_max", type=int, default=6)
ap.add_argument("--n_boot", type=int, default=50); ap.add_argument("--top_var", type=int, default=2000)
args = ap.parse_args()

os.makedirs(os.path.dirname(args.out_labels), exist_ok=True)

counts = pd.read_csv(args.counts, sep="\t", index_col=0)
meta   = pd.read_csv(args.meta,   sep="\t")
rb_pos = meta.loc[meta["Infection_Group"] == "RB(+)", "Sample"].tolist()
rb_pos = [s for s in rb_pos if s in counts.columns]

X = counts[rb_pos].astype(float)
lib = X.sum()
lcpm = np.log2(X.div(lib, axis=1) * 1e6 + 1)
top = lcpm.var(axis=1).nlargest(args.top_var).index
M = lcpm.loc[top, :].T.values                   # samples × genes
M = M - M.min() + 1e-6                          # non-negative

ks = list(range(args.k_min, args.k_max + 1))
results = {}
n = M.shape[0]

for k in ks:
    C = np.zeros((n, n))
    for b in range(args.n_boot):
        rng = np.random.RandomState(b)
        idx = rng.choice(M.shape[1], int(0.8 * M.shape[1]), replace=False)
        nmf = NMF(n_components=k, init="nndsvd", max_iter=500, random_state=b)
        W = nmf.fit_transform(M[:, idx])
        lab = np.argmax(W, axis=1)
        C += (lab[:, None] == lab[None, :]).astype(float)
    C /= args.n_boot
    d = 1 - C
    np.fill_diagonal(d, 0)
    Z = linkage(squareform(d, checks=False), method="average")
    coph, _ = cophenet(Z, squareform(d, checks=False))
    # dispersion (Brunet)
    disp = np.sum(4 * (C - 0.5) ** 2) / (n * n)
    results[k] = {"consensus": C, "coph": coph, "dispersion": disp, "linkage": Z}
    print(f"k={k}  coph={coph:.3f}  disp={disp:.3f}")

best_k = max(ks, key=lambda k: results[k]["coph"])
labels = fcluster(results[best_k]["linkage"], t=best_k, criterion="maxclust")
pd.DataFrame({"Sample": rb_pos, "Subtype": [f"S{int(x)}" for x in labels],
              "best_k": best_k}).to_csv(args.out_labels, index=False)

# Cophenetic/dispersion plot
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(ks, [results[k]["coph"] for k in ks], "o-"); ax[0].set(title="Cophenetic", xlabel="k")
ax[1].plot(ks, [results[k]["dispersion"] for k in ks], "o-"); ax[1].set(title="Dispersion", xlabel="k")
ax[0].axvline(best_k, ls="--", color="red"); ax[1].axvline(best_k, ls="--", color="red")
fig.tight_layout(); fig.savefig(args.out_coph, dpi=200); plt.close(fig)

# Consensus heatmap for best_k
order = np.argsort(labels)
C = results[best_k]["consensus"][np.ix_(order, order)]
fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(C, cmap="viridis", aspect="auto", vmin=0, vmax=1)
ax.set(title=f"Consensus matrix, k={best_k}", xticks=[], yticks=[])
fig.colorbar(im, ax=ax, label="Co-cluster fraction"); fig.tight_layout()
fig.savefig(args.out_consensus, dpi=200); plt.close(fig)
print(f"Best k = {best_k}")
