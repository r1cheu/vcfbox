# Rebuild the F1 pedigree from a `vcfbox parentage` report.
#
# usage: uv run --with pandas,numpy python scripts/f1_pedigree.py \
#            f1.parentage.tsv data/f0.02/parents.len.bcf OUT_DIR [--twin 0.003]
#
# For each sample, the first rule that holds: depth < MIN_DEPTH drops; a
# paper parent in DISCARDED_LINES drops; the paper pair best up to twins keeps
# the paper names (fit_flag marks a poor fit); a paper parent outside the
# panel drops; ONE_PARENT_ABSENT / MISFIT / hom_mm above q99 drop; a parent
# sample (paper pair A x A) with a detected cross drops; lod_vs_second below
# MIN_LOD drops unless the runner-up differs from the best only by twins;
# otherwise the detected pair (the lone GS3 as father), with a twin of a paper
# parent put back to its paper name, kept only as a GS1/GS2 x GS3 cross or,
# for a parent sample, as the detected line B x B.
# Writes OUT_DIR/new_pedigree.tsv (sample, father, mother, expected_father,
# expected_mother),
# OUT_DIR/new_pedigree_annotated.tsv (every F1 with its features and reason),
# OUT_DIR/panel_distance.tsv.gz (pairwise hom-site discordance of the panel),
# OUT_DIR/calibration.txt.
import argparse
import os
import subprocess

import numpy as np
import pandas as pd

MIN_DEPTH = 0.2
DEPTH_BINS = [0.2, 1.0, 3.0, np.inf]
# GS1 and GS2 are male-sterile lines: a cross is a GS1/GS2 mother by a GS3
# father.
MOTHER_GROUPS = {"GS1", "GS2"}
FATHER_GROUP = "GS3"
# lod_vs_second (natural log) below which a sample is AMBIGUOUS, as
# `vcfbox parentage --min-lod`.
MIN_LOD = 50
# Lines discarded by the breeders (segregating); a sample naming one on paper
# or in its final pair drops.
DISCARDED_LINES = {"GS3-110"}


def panel_distance(bcf, cache):
    # Fraction of sites, among those homozygous in both lines, where the two
    # lines carry different alleles; every 10th panel site.
    if os.path.exists(cache):
        return pd.read_csv(cache, sep="\t", index_col=0)
    names = subprocess.run(
        ["bcftools", "query", "-l", bcf], check=True, capture_output=True, text=True
    ).stdout.split()
    query = subprocess.Popen(
        ["bcftools", "query", "-f", "[%GT ]\n", bcf], stdout=subprocess.PIPE, text=True
    )
    rows = []
    for i, line in enumerate(query.stdout):
        if i % 10 == 0:
            rows.append([int(gt[0]) + int(gt[2]) for gt in line.split()])
    if query.wait() != 0:
        raise RuntimeError(f"bcftools query failed on {bcf}")
    g = np.array(rows, dtype=np.int8)
    hom0 = (g == 0).astype(np.float32)
    hom2 = (g == 2).astype(np.float32)
    disc = hom0.T @ hom2
    disc = disc + disc.T
    both_hom = (hom0 + hom2).T @ (hom0 + hom2)
    dist = pd.DataFrame(disc / both_hom, index=names, columns=names)
    dist.to_csv(cache, sep="\t", float_format="%.6f")
    return dist


def group(name):
    return name.split("-")[0]


def valid_cross(mother, father):
    return group(mother) in MOTHER_GROUPS and group(father) == FATHER_GROUP


def orient(paper_mother, paper_father, a, b, is_twin):
    # Returns the detected pair as (mother, father): a lone GS3 is the father;
    # otherwise (a pair that is dropped as an invalid cross anyway) a detected
    # parent equal or twin to a paper parent keeps that side, else (a, b).
    if (group(a) == FATHER_GROUP) != (group(b) == FATHER_GROUP):
        return (b, a) if group(a) == FATHER_GROUP else (a, b)
    same = lambda p, q: p == q or is_twin(p, q)
    straight = same(paper_mother, a) or same(paper_father, b)
    crossed = same(paper_mother, b) or same(paper_father, a)
    return (b, a) if crossed and not straight else (a, b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("report")
    ap.add_argument("panel")
    ap.add_argument("out_dir")
    ap.add_argument(
        "--twin",
        type=float,
        default=0.003,
        help="panel distance below which two lines are twins (not told apart"
        " for AMBIGUOUS, listed in name_ambiguous)",
    )
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    dist = panel_distance(
        args.panel, os.path.join(args.out_dir, "panel_distance.tsv.gz")
    )
    report = pd.read_csv(args.report, sep="\t")
    # F1s (dropping a detected selfing) and parent samples (paper A x A).
    paper_self = report.expected_a == report.expected_b
    report["kind"] = np.where(paper_self, "parent", "F1")
    samples = report[
        report.expected_a.notna() & (paper_self | (report.best_a != report.best_b))
    ].copy()
    samples["depth_bin"] = pd.cut(samples.depth_mean, DEPTH_BINS, right=False)
    f1 = samples[samples.kind == "F1"]

    known_good = f1[
        (f1.expected_status == "MATCH")
        & (f1.verdict == "PASS")
        & (f1.depth_mean >= MIN_DEPTH)
    ]
    q99 = known_good.groupby("depth_bin", observed=False).hom_mismatch_rate.quantile(
        0.99
    )
    samples["hom_mm_q99"] = samples.depth_bin.map(q99).astype(float)

    def d(a, b):
        return float(dist.loc[a, b])

    def is_twin(a, b):
        return d(a, b) < args.twin

    def twins(name):
        row = dist.loc[name]
        return sorted(n for n in row.index if n != name and is_twin(name, n))

    def same_pair(a, b, p, q):
        # (a, b) is (p, q) up to twins, in either order.
        same = lambda x, y: x == y or is_twin(x, y)
        return (same(a, p) and same(b, q)) or (same(a, q) and same(b, p))

    def twin_of_paper(r):
        return r.expected_status != "ABSENT" and same_pair(
            r.best_a, r.best_b, r.expected_a, r.expected_b
        )

    def twin_second(r):
        # The runner-up pair differs from the best only by twins.
        return isinstance(r.second_a, str) and same_pair(
            r.best_a, r.best_b, r.second_a, r.second_b
        )

    records = []
    for r in samples.itertuples():
        parent = r.kind == "parent"
        # vcfbox never calls a selfing best pair AMBIGUOUS, so the lod rule is
        # applied here for both kinds. A best pair ahead only of its twin pair
        # is not ambiguous.
        ambiguous = r.lod_vs_second < MIN_LOD and not twin_second(r)
        rec = {
            "sample": r.sample,
            "kind": r.kind,
            "paper_mother": r.expected_a,
            "paper_father": r.expected_b,
            "best_a": r.best_a,
            "best_b": r.best_b,
            "status": r.expected_status,
            "verdict": r.verdict,
            "depth": r.depth_mean,
            "hom_mm": r.hom_mismatch_rate,
            "hom_mm_q99": r.hom_mm_q99,
            "expected_rank": r.expected_rank,
            "expected_lod_vs_best": r.expected_lod_vs_best,
        }
        mother = father = None
        if r.depth_mean < MIN_DEPTH:
            reason = "drop_low_depth"
        elif {r.expected_a, r.expected_b} & DISCARDED_LINES:
            reason = "drop_discarded_line"
        elif twin_of_paper(r):
            # MATCH (paper pair ranked first) is the exact case of this.
            mother, father = r.expected_a, r.expected_b
            reason = "match" if r.expected_status == "MATCH" else "twin_paper"
            if r.verdict in ("ONE_PARENT_ABSENT", "MISFIT"):
                rec["fit_flag"] = r.verdict
            elif ambiguous:
                rec["fit_flag"] = "AMBIGUOUS"
            elif r.hom_mismatch_rate > r.hom_mm_q99:
                rec["fit_flag"] = "HOM_MM_ABOVE_Q99"
        elif r.expected_status == "ABSENT":
            reason = "drop_absent"
        elif r.verdict in ("ONE_PARENT_ABSENT", "MISFIT"):
            reason = f"drop_verdict_{r.verdict}"
        elif r.hom_mismatch_rate > r.hom_mm_q99:
            reason = "drop_hom_mm_above_q99"
        elif parent and r.best_a != r.best_b:
            reason = "drop_detected_cross"
        elif ambiguous:
            reason = "drop_verdict_AMBIGUOUS"
        else:
            # The best pair is not the paper pair up to twins, so at least one
            # side is a new (detected) parent.
            det_mother, det_father = (
                (r.best_a, r.best_a)
                if parent
                else orient(r.expected_a, r.expected_b, r.best_a, r.best_b, is_twin)
            )
            sides = []
            for paper, det in ((r.expected_a, det_mother), (r.expected_b, det_father)):
                dd = d(paper, det)
                if paper == det:
                    sides.append((paper, "same", dd))
                elif is_twin(paper, det):
                    sides.append((paper, "twin_paper", dd))
                else:
                    sides.append((det, "detected", dd))
            (mother, m_src, m_d), (father, f_src, f_d) = sides
            rec.update(
                det_mother=det_mother,
                det_father=det_father,
                mother_source=m_src,
                father_source=f_src,
                d_mother=m_d,
                d_father=f_d,
            )
            if not parent and not valid_cross(mother, father):
                reason = "drop_invalid_cross"
                mother = father = None
            else:
                reason = "detected"
                amb = dict.fromkeys(
                    f"{n} or {'/'.join(twins(n))}"
                    for n, src in ((mother, m_src), (father, f_src))
                    if src == "detected" and twins(n)
                )
                rec["name_ambiguous"] = "; ".join(amb)
        if {mother, father} & DISCARDED_LINES:
            reason, mother, father = "drop_discarded_line", None, None
        rec.update(reason=reason, mother=mother, father=father)
        records.append(rec)

    out = pd.DataFrame(records)
    kept = out[out.mother.notna()]

    # Every pedigree pair must be either the paper pair, the detected pair, or
    # the detected pair with twins put back to their paper names.
    for r in kept.itertuples():
        pair = {r.mother, r.father}
        if r.kind == "parent":
            if r.mother != r.father:
                raise AssertionError(f"{r.sample}: parent sample as a cross")
        elif r.mother == r.father:
            raise AssertionError(f"{r.sample}: selfing in pedigree")
        elif not valid_cross(r.mother, r.father):
            raise AssertionError(
                f"{r.sample}: {r.mother} x {r.father} is not a GS1/GS2 x GS3 cross"
            )
        if r.reason == "detected":
            for new, det in ((r.mother, r.det_mother), (r.father, r.det_father)):
                if new != det and not is_twin(new, det):
                    raise AssertionError(
                        f"{r.sample}: {new} is neither detected {det} nor its twin"
                    )
        elif pair != {r.paper_mother, r.paper_father}:
            raise AssertionError(f"{r.sample}: {r.reason} row differs from paper")

    kept = kept.rename(
        columns={"paper_father": "expected_father", "paper_mother": "expected_mother"}
    )
    kept[["sample", "father", "mother", "expected_father", "expected_mother"]].to_csv(
        os.path.join(args.out_dir, "new_pedigree.tsv"), sep="\t", index=False
    )
    out.to_csv(
        os.path.join(args.out_dir, "new_pedigree_annotated.tsv"),
        sep="\t",
        index=False,
        float_format="%.6g",
    )
    with open(os.path.join(args.out_dir, "calibration.txt"), "w") as fh:
        fh.write(f"twin threshold: panel distance < {args.twin}\n")
        fh.write(
            f"known-good rows (MATCH, PASS, depth >= {MIN_DEPTH}): {len(known_good)}\n"
        )
        fh.write("hom_mismatch_rate 99th percentile by depth bin:\n")
        fh.write(q99.to_string() + "\n")
    print(open(os.path.join(args.out_dir, "calibration.txt")).read())
    print(out.groupby(["kind", "reason"]).size().to_string())


if __name__ == "__main__":
    main()
