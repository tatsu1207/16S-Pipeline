"""Generate a synthetic BIOM table with known differentially abundant ASVs.

16 samples (8 Control, 8 Treatment) x 300 ASVs. Taxonomy and sequences are
borrowed from real ASVs (the tutorial dataset) so labels and plots look
realistic; counts are simulated. 20 ASVs are spiked: 10 increased and 10
decreased in Treatment by 3-10 fold in true (absolute) abundance.

Outputs (next to this script):
  da_demo.biom           counts + taxonomy + sequences + sample metadata
  da_demo_metadata.tsv   sample metadata for the DA page's metadata upload
  da_demo_truth.tsv      the spiked ASVs and their true fold change

Usage (from the repo root):
  conda run -n microbiome_16S python test_samples/da_demo/make_da_demo.py [source.biom]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from biom import Table, load_table

OUT = Path(__file__).parent
SOURCE = sys.argv[1] if len(sys.argv) > 1 else "data/datasets/1/asv_table.biom"
N_ASV, N_PER_GROUP, N_SPIKED = 300, 8, 20

rng = np.random.default_rng(16)

# Real ASV IDs, taxonomy and sequences, most abundant first
src = load_table(SOURCE)
order = np.argsort(-np.asarray(src.sum(axis="observation")))
obs_ids = [src.ids(axis="observation")[i] for i in order[:N_ASV]]
obs_md = []
for oid in obs_ids:
    md = src.metadata(oid, axis="observation") or {}
    obs_md.append({"taxonomy": list(md.get("taxonomy", [])),
                   "sequence": md.get("sequence", "")})

# Baseline community: skewed abundances
base = rng.lognormal(mean=0, sigma=2.0, size=N_ASV)
base = np.sort(base)[::-1]
base /= base.sum()

# Spike 20 ASVs among moderately abundant ones (ranks 5-120), alternating up/down
spiked = rng.choice(np.arange(5, 120), size=N_SPIKED, replace=False)
fold = np.ones(N_ASV)
truth = []
for k, idx in enumerate(spiked):
    fc = rng.uniform(3, 10)
    direction = "up" if k % 2 == 0 else "down"
    fold[idx] = fc if direction == "up" else 1 / fc
    truth.append({"feature": obs_ids[idx], "direction": direction,
                  "true_fold_change": round(fold[idx], 3),
                  "true_log2fc": round(float(np.log2(fold[idx])), 3),
                  "taxonomy": "; ".join(t for t in obs_md[idx]["taxonomy"] if t)})

# Chance of being absent from a sample rises with rarity (rank 0 = most abundant)
p_absent = 0.7 * (np.arange(N_ASV) / N_ASV) ** 3

samples, groups, counts = [], [], []
for group in ("Control", "Treatment"):
    for i in range(1, N_PER_GROUP + 1):
        sid = f"DA_{group[0]}{i:02d}"
        # True abundances: baseline (x fold in Treatment) with biological noise
        mu = base * (fold if group == "Treatment" else 1.0)
        mu = mu * rng.lognormal(0, 0.5, size=N_ASV)
        # Rare taxa are often absent from a sample; abundant ones almost never
        mu[rng.random(N_ASV) < p_absent] = 0
        props = mu / mu.sum()
        depth = int(rng.integers(15_000, 45_000))
        counts.append(rng.multinomial(depth, props))
        samples.append(sid)
        groups.append(group)

data = np.array(counts).T
sample_md = [{"group": g, "depth": int(data[:, j].sum())} for j, g in enumerate(groups)]
table = Table(data, obs_ids, samples, observation_metadata=obs_md,
              sample_metadata=sample_md, table_id="da_demo")
with __import__("h5py").File(OUT / "da_demo.biom", "w") as f:
    table.to_hdf5(f, "16S-Pipeline DA demo (synthetic)")

pd.DataFrame({"sample-id": samples, "group": groups}).to_csv(
    OUT / "da_demo_metadata.tsv", sep="\t", index=False)
pd.DataFrame(truth).sort_values(["direction", "true_log2fc"]).to_csv(
    OUT / "da_demo_truth.tsv", sep="\t", index=False)

print(f"{len(samples)} samples x {len(obs_ids)} ASVs, "
      f"{int((data > 0).mean() * 100)}% non-zero, depth {data.sum(0).min()}-{data.sum(0).max()}")
