"""Shared helpers for the mapquik ONT-R10 vs HiFi evaluation.

Per (strain, tech): subsample reads -> minimap2 truth -> mapquik run -> compare.
Large files live in scratch/r10/, small results (tsv) in experiments/r10/results/.
"""
import gzip
import json
import random
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
MAPQUIK = REPO / "target/release/mapquik"
DATA = Path("/home/thomas/Projects/mindtherest/data")
SCRATCH = REPO / "scratch/r10"
RESULTS = HERE / "results"

STRAINS = ["FAM14197", "FAM20838", "FAM20851", "FAM20853", "FAM20856",
           "FAM20862", "FAM20866", "FAM22000", "FAM29080"]
TECHS = ["pb", "ont"]
MM2_PRESET = {"pb": "map-hifi", "ont": "map-ont"}
N_READS = 20000
SEED = 1
THREADS = 16
MIN_TRUTH_MAPQ = 30

ID_BINS = [0, .90, .95, .97, .98, .99, .995, 1.01]
ID_LABELS = ["<90", "90-95", "95-97", "97-98", "98-99", "99-99.5", ">=99.5"]
LEN_BINS = [0, 5000, 10000, 20000, np.inf]
LEN_LABELS = ["<5k", "5-10k", "10-20k", ">=20k"]

# Mapquik parameter name -> CLI flag. hpc=False -> --nohpc.
PARAM_FLAGS = {"k": "--k", "l": "--l", "density": "--density", "chain": "--chain",
               "seed": "--seed", "gap_diff": "--gap-diff"}
DEFAULTS = {"k": 5, "l": 31, "density": 0.01, "chain": 4, "seed": 11, "gap_diff": 2000, "hpc": True}


def workdir(strain, tech):
    d = SCRATCH / "work" / strain / tech
    d.mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------- timing
def run_timed(cmd, stdout, stderr, workdir_):
    """Run cmd under /usr/bin/time -v; return dict(wall_s, max_rss_mb)."""
    tfile = workdir_ / "time.tmp"
    with open(stdout, "w") as so, open(stderr, "w") as se:
        subprocess.run(["/usr/bin/time", "-v", "-o", str(tfile)] + [str(c) for c in cmd],
                       stdout=so, stderr=se, check=True)
    txt = tfile.read_text()
    wall = re.search(r"Elapsed \(wall clock\) time.*: (.*)", txt).group(1).split(":")
    wall_s = sum(float(x) * 60 ** i for i, x in enumerate(reversed(wall)))
    rss = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", txt).group(1))
    return {"wall_s": wall_s, "max_rss_mb": rss / 1024}


# ---------------------------------------------------------------- inputs
def reference(strain):
    """Symlink the PacBio assembly as ref.fasta (mapquik needs a .fasta extension)."""
    d = SCRATCH / "work" / strain
    d.mkdir(parents=True, exist_ok=True)
    ref = d / "ref.fasta"
    if not ref.exists():
        ref.symlink_to(DATA / "pb" / strain / "assembly.fasta")
    return ref


def subsample(strain, tech, n=N_READS, seed=SEED):
    """Reservoir-sample n reads (fixed seed). Writes reads.fastq.gz and reads.tsv (name, length)."""
    wd = workdir(strain, tech)
    out, tsv = wd / "reads.fastq.gz", wd / "reads.tsv"
    if out.exists() and tsv.exists():
        return out
    rng = random.Random(seed)
    res = []  # reservoir of (index, record)
    with gzip.open(DATA / tech / strain / "reads.fastq.gz", "rt") as fh:
        i = -1
        while True:
            rec = [fh.readline() for _ in range(4)]
            if not rec[0]:
                break
            i += 1
            if len(res) < n:
                res.append((i, rec))
            else:
                j = rng.randint(0, i)
                if j < n:
                    res[j] = (i, rec)
    res.sort(key=lambda x: x[0])
    names, lens = [], []
    with gzip.open(out, "wt", compresslevel=3) as fo:
        for _, rec in res:
            fo.writelines(rec)
            names.append(rec[0][1:].split()[0])
            lens.append(len(rec[1].strip()))
    pd.DataFrame({"name": names, "length": lens}).to_csv(tsv, sep="\t", index=False)
    return out


# ---------------------------------------------------------------- minimap2 truth
def truth(strain, tech):
    """Run minimap2 once, return per-read primary alignment table (cached as truth.tsv)."""
    wd = workdir(strain, tech)
    tsv, info = wd / "truth.tsv", wd / "truth.json"
    if tsv.exists() and info.exists():
        return pd.read_csv(tsv, sep="\t"), json.loads(info.read_text())
    reads = subsample(strain, tech)
    paf = wd / "minimap2.paf"
    t = run_timed(["minimap2", "-x", MM2_PRESET[tech], "-c", "--secondary=no", "-t", THREADS,
                   reference(strain), reads], paf, wd / "minimap2.log", wd)
    rows = []
    with open(paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            tags = {x[:2]: x[5:] for x in f[12:]}
            if tags.get("tp") != "P":
                continue
            rows.append((f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], int(f[6]), int(f[7]),
                         int(f[8]), int(f[9]), int(f[10]), int(f[11]), float(tags["de"])))
    df = pd.DataFrame(rows, columns=["name", "qlen", "qs", "qe", "strand", "tname", "tlen", "ts", "te",
                                     "nmatch", "alen", "mapq", "de"])
    # Split reads (e.g. across the origin of a circular contig) give several tp:P lines.
    # `best` marks the longest query span; it defines MAPQ/identity of the read. Other pieces
    # are only used to accept a mapquik hit as correct.
    df = df.assign(span=df.qe - df.qs).sort_values(["name", "span"], ascending=[True, False])
    df["best"] = ~df.name.duplicated()
    df["identity"] = 1 - df.de
    df["identity_blast"] = df.nmatch / df.alen
    df = df.drop(columns="span")
    df.to_csv(tsv, sep="\t", index=False)
    t["n_alignments"] = len(df)
    info.write_text(json.dumps(t))
    return df, t


# ---------------------------------------------------------------- mapquik
def setting_args(params):
    p = {**DEFAULTS, **params}
    args = []
    for key, flag in PARAM_FLAGS.items():
        if p[key] != DEFAULTS[key]:
            args += [flag, str(p[key])]
    if not p["hpc"]:
        args.append("--nohpc")
    return args


def run_mapquik(strain, tech, name, params, force=False):
    """Run mapquik with the given parameter overrides; return (paf path, json path, timing dict)."""
    wd = workdir(strain, tech)
    prefix, js = wd / f"mq.{name}", wd / f"mq.{name}.json"
    paf, tinfo = Path(f"{prefix}.paf"), wd / f"mq.{name}.time.json"
    if not force and paf.exists() and js.exists() and tinfo.exists():
        return paf, js, json.loads(tinfo.read_text())
    cmd = [MAPQUIK, "--reference", reference(strain), "-p", prefix, "--threads", THREADS,
           "--json", js] + setting_args(params) + [subsample(strain, tech)]
    t = run_timed(cmd, wd / f"mq.{name}.stdout", wd / f"mq.{name}.stderr", wd)
    tinfo.write_text(json.dumps(t))
    return paf, js, t


def read_mapquik_paf(paf):
    cols = ["name", "qlen", "qs", "qe", "strand", "tname", "tlen", "ts", "te", "score", "tlen2", "mq_mapq"]
    df = pd.read_csv(paf, sep="\t", header=None, names=cols, usecols=range(12))
    ndup = int(df.name.duplicated().sum())
    df = df.sort_values("mq_mapq", ascending=False).drop_duplicates("name")
    return df, ndup


# ---------------------------------------------------------------- comparison
def union_len(iv):
    iv = sorted(iv)
    tot, end = 0, -1
    for s, e in iv:
        s = max(s, end)
        if e > s:
            tot += e - s
            end = e
    return tot


def counts(sub):
    """Metrics for a set of truth-set reads (columns: cls, mq_mapq)."""
    n = len(sub)
    c, w = (sub.cls == "correct"), (sub.cls == "wrong")
    c60, w60 = c & (sub.mq_mapq == 60), w & (sub.mq_mapq == 60)
    div = lambda a, b: float(a / b) if b else float("nan")
    return {"n_truth": n, "correct": int(c.sum()), "wrong": int(w.sum()), "unmapped": int((sub.cls == "unmapped").sum()),
            "correct60": int(c60.sum()), "wrong60": int(w60.sum()),
            "sens": div(c.sum(), n), "wrong_rate": div(w.sum(), c.sum() + w.sum()),
            "sens60": div(c60.sum(), n), "wrong_rate60": div(w60.sum(), c60.sum() + w60.sum())}


def evaluate(strain, tech, name, params, truth_df, mq_paf, mq_json, timing, mm2_timing):
    reads = pd.read_csv(workdir(strain, tech) / "reads.tsv", sep="\t")
    mq, ndup = read_mapquik_paf(mq_paf)
    best = truth_df[truth_df.best]
    tr = best[["name", "mapq", "identity", "identity_blast"]].rename(columns={"mapq": "t_mapq"})
    mq = mq[["name", "tname", "strand", "ts", "te", "mq_mapq"]]
    df = reads.merge(tr, on="name", how="left").merge(mq, on="name", how="left")
    mapped = df.tname.notna()
    # a mapquik hit agrees if it matches ANY minimap2 primary piece of the read
    m = mq.merge(truth_df[["name", "tname", "strand", "ts", "te"]], on="name", suffixes=("", "_t"))
    ov = np.minimum(m.te, m.te_t) - np.maximum(m.ts, m.ts_t)
    short = np.minimum(m.te - m.ts, m.te_t - m.ts_t)
    ok = (m.tname == m.tname_t) & (m.strand == m.strand_t) & (ov >= 0.5 * short)
    agree = df.name.isin(set(m.name[ok]))
    in_truth = df.t_mapq >= MIN_TRUTH_MAPQ
    df["cls"] = np.where(~in_truth, "", np.where(~mapped, "unmapped", np.where(agree, "correct", "wrong")))
    t = df[in_truth].copy()
    t["id_bin"] = pd.cut(t.identity, ID_BINS, labels=ID_LABELS, right=False)
    t["len_bin"] = pd.cut(t.length, LEN_BINS, labels=LEN_LABELS, right=False)

    strata = [("all", "all", counts(t))]
    for col, lab in (("id_bin", ID_LABELS), ("len_bin", LEN_LABELS)):
        for b in lab:
            strata.append((col, b, counts(t[t[col] == b])))
    bins = pd.DataFrame([{"strain": strain, "tech": tech, "setting": name, "stratum_type": a, "stratum": b, **m}
                         for a, b, m in strata])

    # mapquik mappings for reads outside the truth set
    nt = df[~in_truth]
    unm = nt[nt.t_mapq.isna()]
    low = nt[nt.t_mapq.notna()]
    extra = {"n_mm2_unmapped": len(unm), "mq_mapped_mm2_unmapped": int(unm.tname.notna().sum()),
             "n_mm2_lowmapq": len(low), "mq_mapped_mm2_lowmapq": int(low.tname.notna().sum()),
             "mq_agrees_mm2_lowmapq": int(agree[low.index].sum())}

    # contig coverage: mapquik JSON vs minimap2 primary alignments (all MAPQ)
    mqj = json.loads(Path(mq_json).read_text())
    cov_rows = []
    for contig, c in mqj["contig_coverage"].items():
        sub = truth_df[(truth_df.tname == contig) & truth_df.best]
        L = c["length_bp"]
        cov_rows.append({"strain": strain, "tech": tech, "setting": name, "contig": contig, "length": L,
                         "mq_breadth": c["breadth_fraction"], "mq_depth": c["average_depth"], "mq_reads": c["mapped_reads"],
                         "mm2_breadth": union_len(list(zip(sub.ts, sub.te))) / L,
                         "mm2_depth": float((sub.te - sub.ts).sum() / L), "mm2_reads": len(sub)})
    cov = pd.DataFrame(cov_rows)

    allm = counts(t)
    run = {"strain": strain, "tech": tech, "setting": name, **{k: v for k, v in {**DEFAULTS, **params}.items()},
           "n_reads": len(reads), "n_mq_mapped": int(mapped.sum()), "n_mq_dup_lines": ndup,
           **allm, **extra,
           "mq_wall_s": timing["wall_s"], "mq_rss_mb": timing["max_rss_mb"],
           "mm2_wall_s": mm2_timing["wall_s"], "mm2_rss_mb": mm2_timing["max_rss_mb"],
           "mm2_mapped": int(df.t_mapq.notna().sum()),
           "median_identity": float(df.identity.median()), "mean_read_len": float(reads.length.mean())}
    return run, bins, cov


def run_and_evaluate(strain, tech, name, params, truth_cache=None):
    """Full per-run pipeline with results cached under scratch/r10/results/<strain>/<tech>/<name>.*"""
    out = SCRATCH / "results" / strain / tech
    out.mkdir(parents=True, exist_ok=True)
    f = out / f"{name}.json"
    if f.exists():
        return json.loads(f.read_text())["run"]
    tdf, mm2_t = truth_cache or truth(strain, tech)
    paf, js, t = run_mapquik(strain, tech, name, params)
    run, bins, cov = evaluate(strain, tech, name, params, tdf, paf, js, t, mm2_t)
    f.write_text(json.dumps({"run": run, "bins": bins.to_dict("records"), "cov": cov.to_dict("records")}))
    return run


def load_results(strains=None, techs=None, settings=None):
    """Concatenate cached results -> (runs, bins, cov) DataFrames."""
    runs, bins, cov = [], [], []
    for f in sorted((SCRATCH / "results").glob("*/*/*.json")):
        strain, tech = f.parts[-3], f.parts[-2]
        if (strains and strain not in strains) or (techs and tech not in techs):
            continue
        d = json.loads(f.read_text())
        if settings and d["run"]["setting"] not in settings:
            continue
        runs.append(d["run"])
        bins += d["bins"]
        cov += d["cov"]
    return pd.DataFrame(runs), pd.DataFrame(bins), pd.DataFrame(cov)
