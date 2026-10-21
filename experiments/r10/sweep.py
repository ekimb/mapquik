"""ONT parameter sweep (pooled over the sweep strains, ONT reads only; truth is reused).

Stage 1: one-factor-at-a-time around mapquik defaults (all six parameters).
Stage 2: full grid over k x l x density x hpc.
Stage 2b: extension of the grid to higher densities / smaller k, l (the stage-2 optimum was on the grid edge).
Stage 3: one-factor-at-a-time of chain / seed / gap_diff around the best stage-2 setting.
Selection: maximise sensitivity subject to wrong-rate <= max(0.5%, wrong-rate of HiFi defaults);
because sensitivity saturates, see select() for the plateau rule used to break the near-ties.

    pixi run python sweep.py --strains FAM14197 FAM22000
"""
import argparse
import itertools
import subprocess

import pandas as pd

import r10lib as r

KS, LS, DS = [2, 3, 4, 5, 6], [13, 15, 17, 19, 21, 25, 31], [0.005, 0.01, 0.02, 0.03, 0.05]


def name_of(p):
    p = {**r.DEFAULTS, **p}
    return (f"sw_k{p['k']}_l{p['l']}_d{p['density']}_hpc{int(p['hpc'])}"
            f"_c{p['chain']}_s{p['seed']}_g{p['gap_diff']}")


def run_settings(strains, settings):
    caches = {s: r.truth(s, "ont") for s in strains}
    for i, p in enumerate(settings):
        for s in strains:
            try:
                r.run_and_evaluate(s, "ont", name_of(p), p, caches[s])
            except subprocess.CalledProcessError as e:
                print(f"FAILED {s} {name_of(p)}: {e}", flush=True)
        print(f"  [{i + 1}/{len(settings)}] {name_of(p)}", flush=True)


def pooled(strains):
    runs, _, _ = r.load_results(strains=strains, techs=["ont"])
    runs = runs[runs.setting.str.startswith("sw_")]
    g = runs.groupby("setting")
    out = g[["n_truth", "correct", "wrong", "unmapped", "correct60", "wrong60", "mq_wall_s"]].sum()
    out["n_strains"] = g.size()
    out = out[out.n_strains == len(strains)].copy()
    out["sens"] = out.correct / out.n_truth
    out["wrong_rate"] = out.wrong / (out.correct + out.wrong)
    out["sens60"] = out.correct60 / out.n_truth
    out["wrong_rate60"] = out.wrong60 / (out.correct60 + out.wrong60)
    first = g[["k", "l", "density", "hpc", "chain", "seed", "gap_diff"]].first()
    return first.join(out).reset_index()


def select(table, max_wrong):
    """Plateau rule: sensitivity saturates (>99.9%) for many settings, so rather than picking an
    arbitrary extreme point take the fastest setting that is within 0.02 pp of the best sensitivity
    and within 0.5 pp of the best MAPQ-60 sensitivity (among settings meeting the wrong-rate limit)."""
    ok = table[table.wrong_rate <= max_wrong].copy()
    ok["plateau"] = (ok.sens >= ok.sens.max() - 0.0002) & (ok.sens60 >= ok.sens60.max() - 0.005)
    return ok.sort_values(["plateau", "mq_wall_s"], ascending=[False, True]).reset_index(drop=True)


def best_row(strains, max_wrong):
    return select(pooled(strains), max_wrong).iloc[0]


def ofat(center, grid):
    out = []
    for key, vals in grid.items():
        out += [{**center, key: v} for v in vals]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--strains", nargs="+", default=["FAM14197", "FAM22000"])
    a = ap.parse_args()
    D = {k: v for k, v in r.DEFAULTS.items()}
    # HiFi-default wrong-rate (pooled) defines the allowed error budget
    runs, _, _ = r.load_results(strains=a.strains, techs=["pb"], settings=["default"])
    max_wrong = max(0.005, runs.wrong.sum() / (runs.correct.sum() + runs.wrong.sum())) if len(runs) else 0.005
    print("max wrong rate:", max_wrong)

    print("stage 1: one-factor-at-a-time")
    run_settings(a.strains, [D] + ofat(D, {"k": KS, "l": LS, "density": DS, "hpc": [True, False],
                                           "chain": [2, 3, 6, 8], "seed": [3, 5, 8, 15], "gap_diff": [500, 1000, 5000]}))
    print("stage 2: grid k x l x d x hpc")
    run_settings(a.strains, [dict(D, k=k, l=l, density=d, hpc=h)
                             for k, l, d, h in itertools.product(KS, LS, DS, [True, False])])
    print("stage 2b: grid extension")
    run_settings(a.strains, [dict(D, k=k, l=l, density=d, hpc=h)
                             for k, l, d, h in itertools.product([1, 2, 3, 5], [11, 13, 17, 31], [0.07, 0.1, 0.2], [True, False])])
    b = best_row(a.strains, max_wrong)
    center = {k: b[k].item() if hasattr(b[k], "item") else b[k] for k in ["k", "l", "density", "hpc", "chain", "seed", "gap_diff"]}
    print("stage 3: chain/seed/gap around", center)
    run_settings(a.strains, ofat(center, {"chain": [2, 3, 4, 6, 8], "seed": [3, 5, 8, 11, 15], "gap_diff": [500, 1000, 2000, 5000]}))
    table = pooled(a.strains)
    r.RESULTS.mkdir(exist_ok=True)
    table.sort_values("sens", ascending=False).to_csv(r.RESULTS / "sweep.tsv", sep="\t", index=False)
    sel = select(table, max_wrong)
    pd.set_option("display.width", 250)
    print(sel.head(15)[["setting", "plateau", "sens", "wrong_rate", "sens60", "mq_wall_s"]].to_string())


if __name__ == "__main__":
    main()
