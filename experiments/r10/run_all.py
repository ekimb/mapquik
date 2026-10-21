"""Run the evaluation: for each strain and tech, subsample, run minimap2 truth, run mapquik with
each named setting (settings.py), compare. Results are cached, so re-running is cheap.

    pixi run all                                   # all 9 strains, pb+ont, default + best_ont
    pixi run python run_all.py --strains FAM14197 FAM22000 --settings default
"""
import argparse

import r10lib as r
from settings import SETTINGS


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--strains", nargs="+", default=r.STRAINS)
    ap.add_argument("--techs", nargs="+", default=r.TECHS)
    ap.add_argument("--settings", nargs="+", default=list(SETTINGS))
    a = ap.parse_args()
    for strain in a.strains:
        for tech in a.techs:
            tc = r.truth(strain, tech)
            for name in a.settings:
                run = r.run_and_evaluate(strain, tech, name, SETTINGS[name], tc)
                print(f"{strain} {tech} {name}: sens={run['sens']:.4f} wrong={run['wrong_rate']:.5f} "
                      f"sens60={run['sens60']:.4f} mq={run['mq_wall_s']:.2f}s mm2={run['mm2_wall_s']:.2f}s", flush=True)


if __name__ == "__main__":
    main()
