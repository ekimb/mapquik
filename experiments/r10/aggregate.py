"""Collect cached per-run results into tidy TSVs in results/ (runs, bins, coverage, identity)."""
import pandas as pd

import r10lib as r


def identity_table(strains):
    rows = []
    for strain in strains:
        for tech in r.TECHS:
            f = r.workdir(strain, tech) / "truth.tsv"
            if not f.exists():
                continue
            t = pd.read_csv(f, sep="\t")
            t = t[t.best]
            reads = pd.read_csv(r.workdir(strain, tech) / "reads.tsv", sep="\t")
            q = t.identity.quantile([.01, .05, .25, .5, .75]).values
            rows.append({"strain": strain, "tech": tech, "n_reads": len(reads), "n_mm2_mapped": len(t),
                         "median_len": reads.length.median(), "n50_len": n50(reads.length),
                         "id_median": q[3], "id_mean": t.identity.mean(), "id_q01": q[0], "id_q05": q[1], "id_q25": q[2],
                         "id_q75": q[4], "blast_id_median": t.identity_blast.median(),
                         "frac_ge97": (t.identity >= .97).mean(), "frac_ge99": (t.identity >= .99).mean(),
                         "frac_lt95": (t.identity < .95).mean()})
    return pd.DataFrame(rows)


def n50(x):
    x = sorted(x, reverse=True)
    tot, acc = sum(x), 0
    for v in x:
        acc += v
        if acc >= tot / 2:
            return v


def main():
    r.RESULTS.mkdir(exist_ok=True)
    runs, bins, cov = r.load_results()
    sw = runs.setting.str.startswith("sw_")
    runs[sw].to_csv(r.RESULTS / "sweep_runs.tsv", sep="\t", index=False)  # per strain; pooled table is sweep.tsv
    runs, bins, cov = runs[~sw], bins[~bins.setting.str.startswith("sw_")], cov[~cov.setting.str.startswith("sw_")]
    runs.to_csv(r.RESULTS / "runs.tsv", sep="\t", index=False)
    bins.to_csv(r.RESULTS / "bins.tsv", sep="\t", index=False)
    cov.to_csv(r.RESULTS / "coverage.tsv", sep="\t", index=False)
    ident = identity_table(sorted(runs.strain.unique()))
    ident.to_csv(r.RESULTS / "identity.tsv", sep="\t", index=False)
    print(f"{len(runs)} runs, {len(bins)} bin rows, {len(cov)} coverage rows, {len(ident)} identity rows -> {r.RESULTS}")


if __name__ == "__main__":
    main()
