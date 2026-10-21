"""Markdown tables from results/*.tsv (no mappings are run). `python report.py > results/tables.md`."""
import pandas as pd

import r10lib as r

R = r.RESULTS
DEV = ["FAM14197", "FAM22000"]  # strains used to choose best_ont
HELD = [s for s in r.STRAINS if s not in DEV]
GROUPS = {"all 9": r.STRAINS, "tuning (2)": DEV, "held-out (7)": HELD}
CONFIGS = [("pb", "default"), ("pb", "best_ont"), ("ont", "default"), ("ont", "best_ont")]
TECH = {"pb": "HiFi", "ont": "ONT"}


def md(df, floatfmt=None):
    floatfmt = floatfmt or {}
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, row in df.iterrows():
        cells = [floatfmt.get(c, "{:.4g}").format(v) if isinstance(v, float) else str(v) for c, v in row.items()]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def pct(x, d=2):
    return f"{100 * x:.{d}f}%" if pd.notna(x) else "-"


def rates(g):
    """Pooled rates from summed counts."""
    n = g.n_truth
    return {"n_truth": int(n), "sens": g.correct / n, "sens60": g.correct60 / n,
            "wrong_rate": g.wrong / (g.correct + g.wrong), "wrong_rate60": g.wrong60 / (g.correct60 + g.wrong60)}


def pooled_bins(bins, kind, strains, tech, setting):
    b = bins[(bins.stratum_type == kind) & bins.strain.isin(strains) & (bins.tech == tech) & (bins.setting == setting)]
    s = b.groupby("stratum")[["n_truth", "correct", "wrong", "unmapped", "correct60", "wrong60"]].sum()
    labs = {"id_bin": r.ID_LABELS, "len_bin": r.LEN_LABELS, "all": ["all"]}[kind]
    out = s.reindex(labs)
    out["sens"] = out.correct / out.n_truth
    out["sens60"] = out.correct60 / out.n_truth
    out["wrong_rate"] = out.wrong / (out.correct + out.wrong)
    return out


def sec_identity(ident):
    rows = []
    for tech in r.TECHS:
        t = pd.concat([pd.read_csv(r.workdir(s, tech) / "truth.tsv", sep="\t").query("best") for s in r.STRAINS])
        reads = pd.concat([pd.read_csv(r.workdir(s, tech) / "reads.tsv", sep="\t") for s in r.STRAINS])
        rows.append({"tech": TECH[tech], "reads": len(reads), "minimap2-mapped": len(t), "median len": f"{reads.length.median():,.0f}",
                     "median id (1-de)": pct(t.identity.median()), "median id (BLAST)": pct(t.identity_blast.median()),
                     "q05": pct(t.identity.quantile(.05)), "q01": pct(t.identity.quantile(.01)),
                     ">=97%": pct((t.identity >= .97).mean()), ">=99%": pct((t.identity >= .99).mean()), "<95%": pct((t.identity < .95).mean())})
    out = ["### Read identity vs the HiFi assembly, pooled over 9 strains\n", md(pd.DataFrame(rows)), "\n",
           "### Read identity per strain\n"]
    t = ident[["strain", "tech", "n_reads", "median_len", "n50_len", "id_median", "id_q05", "id_q01", "blast_id_median", "frac_ge97", "frac_ge99", "frac_lt95"]].copy()
    t["tech"] = t.tech.map(TECH)
    for c in t.columns[5:]:
        t[c] = t[c].map(pct)
    out.append(md(t))
    return "\n".join(out)


def sec_overall(runs):
    out = ["### Pooled results (tech x setting x strain group)\n",
           "Counts are summed over strains. `mapquik s`, `minimap2 s` and RSS are per-strain means; `speedup` = sum(minimap2 s) / sum(mapquik s).\n"]
    rows = []
    for gname, strains in GROUPS.items():
        for tech, setting in CONFIGS:
            g = runs[runs.strain.isin(strains) & (runs.tech == tech) & (runs.setting == setting)]
            s = g[["n_truth", "correct", "wrong", "unmapped", "correct60", "wrong60"]].sum()
            rows.append({"strains": gname, "tech": TECH[tech], "setting": setting, "n_truth": int(s.n_truth), **{
                "sens": pct(s.correct / s.n_truth, 3), "sens60": pct(s.correct60 / s.n_truth, 3),
                "wrong-rate": pct(s.wrong / (s.correct + s.wrong), 3), "wrong60": pct(s.wrong60 / (s.correct60 + s.wrong60), 3),
                "unmapped": int(s.unmapped), "mapquik s": f"{g.mq_wall_s.mean():.2f}", "minimap2 s": f"{g.mm2_wall_s.mean():.2f}",
                "speedup": f"{g.mm2_wall_s.sum() / g.mq_wall_s.sum():.1f}x", "mapquik RSS MB": f"{g.mq_rss_mb.mean():.0f}",
                "minimap2 RSS MB": f"{g.mm2_rss_mb.mean():.0f}"}})
    out.append(md(pd.DataFrame(rows)))
    return "\n".join(out)


def sec_per_strain(runs):
    out = ["### Per strain\n"]
    t = runs.sort_values(["strain", "tech", "setting"])[["strain", "tech", "setting", "n_reads", "n_truth", "sens", "wrong_rate", "sens60", "wrong_rate60",
                                                          "unmapped", "mq_wall_s", "mm2_wall_s"]].copy()
    t["tech"] = t.tech.map(TECH)
    for c in ["sens", "sens60"]:
        t[c] = t[c].map(lambda x: pct(x, 3))
    for c in ["wrong_rate", "wrong_rate60"]:
        t[c] = t[c].map(lambda x: pct(x, 3))
    t["mq_wall_s"] = t.mq_wall_s.map("{:.2f}".format)
    t["mm2_wall_s"] = t.mm2_wall_s.map("{:.2f}".format)
    t = t.rename(columns={"mq_wall_s": "mapquik s", "mm2_wall_s": "minimap2 s"})
    out.append(md(t))
    return "\n".join(out)


def sec_bins(bins, kind, strains, title):
    labs = r.ID_LABELS if kind == "id_bin" else r.LEN_LABELS
    out = [f"### {title}\n", "Cells: sensitivity / MAPQ-60 sensitivity. Wrong-rate is below 0.1% in every cell with n>=100 (see `results/bins.tsv`).\n"]
    cols = {(t, s): pooled_bins(bins, kind, strains, t, s) for t, s in CONFIGS}
    rows = [{"bin": "n (HiFi truth reads)", **{}}]
    hdr = ["bin", "n ONT", "n HiFi"] + [f"{TECH[t]} {s}" for t, s in CONFIGS]
    rows = []
    for lab in labs:
        row = {"bin": lab, "n ONT": int(cols[("ont", "default")].loc[lab, "n_truth"]), "n HiFi": int(cols[("pb", "default")].loc[lab, "n_truth"])}
        for t, s in CONFIGS:
            x = cols[(t, s)].loc[lab]
            row[f"{TECH[t]} {s}"] = f"{pct(x.sens, 1)} / {pct(x.sens60, 1)}" if x.n_truth > 0 else "-"
        rows.append(row)
    out.append(md(pd.DataFrame(rows)[hdr]))
    return "\n".join(out)


def sec_coverage(cov):
    out = ["### Contig coverage: mapquik JSON vs minimap2 primary alignments\n",
           "Depth ratio = mapquik average depth / minimap2 average depth per contig (mapquik reports one alignment per read and does not extend to read ends).\n"]
    c = cov.assign(ratio=cov.mq_depth / cov.mm2_depth)
    rows = []
    for tech, setting in CONFIGS:
        g = c[(c.tech == tech) & (c.setting == setting)]
        rows.append({"tech": TECH[tech], "setting": setting, "contigs": len(g), "min breadth (mapquik)": f"{g.mq_breadth.min():.4f}",
                     "min breadth (minimap2)": f"{g.mm2_breadth.min():.4f}", "median depth ratio": f"{g.ratio.median():.3f}",
                     "min depth ratio": f"{g.ratio.min():.3f}", "contigs with ratio<0.9": int((g.ratio < .9).sum())})
    out.append(md(pd.DataFrame(rows)))
    return "\n".join(out)


def sec_sweep():
    import sweep
    sw = pd.read_csv(R / "sweep.tsv", sep="\t")
    fm = {"sens": "{:.5f}", "sens60": "{:.5f}", "wrong_rate": "{:.5f}", "wrong_rate60": "{:.5f}", "mq_wall_s": "{:.2f}"}
    cols = ["k", "l", "density", "hpc", "chain", "seed", "gap_diff", "sens", "wrong_rate", "sens60", "wrong_rate60", "mq_wall_s"]
    d = sw[(sw.k == 5) & (sw.l == 31) & (sw.density == .01) & sw.hpc & (sw.chain == 4) & (sw.seed == 11) & (sw.gap_diff == 2000)]
    n_diff = sum((sw[k] != v).astype(int) for k, v in r.DEFAULTS.items())
    chosen = sw[(sw.k == 2) & (sw.l == 21) & (sw.density == .03) & sw.hpc & (sw.chain == 4) & (sw.seed == 11) & (sw.gap_diff == 2000)]
    return "\n".join([
        "### Sweep (ONT, FAM14197+FAM22000 pooled): defaults and chosen setting\n", md(pd.concat([d, chosen])[cols], fm), "",
        "### Sweep: one factor at a time around the defaults\n",
        md(sw[n_diff == 1].sort_values(["k", "l", "density", "hpc", "chain", "seed", "gap_diff"])[cols], fm), "",
        "### Sweep: top 10 by plateau rule (fastest among near-best)\n", md(sweep.select(sw, 0.005).head(10)[cols], fm), "",
        "### Sweep: top 5 by raw sensitivity / sens60\n", md(sw.sort_values(["sens", "sens60"], ascending=False).head(5)[cols], fm)])


def main():
    ident = pd.read_csv(R / "identity.tsv", sep="\t")
    runs = pd.read_csv(R / "runs.tsv", sep="\t")
    bins = pd.read_csv(R / "bins.tsv", sep="\t")
    cov = pd.read_csv(R / "coverage.tsv", sep="\t")
    parts = [sec_identity(ident), sec_overall(runs),
             sec_bins(bins, "id_bin", r.STRAINS, "Pooled by read identity (all 9 strains)"),
             sec_bins(bins, "id_bin", HELD, "Pooled by read identity (7 held-out strains)"),
             sec_bins(bins, "len_bin", r.STRAINS, "Pooled by read length (all 9 strains)"),
             sec_bins(bins, "len_bin", HELD, "Pooled by read length (7 held-out strains)"),
             sec_per_strain(runs), sec_coverage(cov), sec_sweep()]
    print("\n\n".join(parts))


if __name__ == "__main__":
    main()
