# Rebuild the F1 pedigree from a `vcfbox parentage` report.
#
# usage: uv run --with pandas,numpy python scripts/f1_pedigree.py \
#            PREFIX.tsv PANEL.bcf OUT_DIR [--twin 0.003]
# PREFIX.topk.tsv is read next to PREFIX.tsv.
#
# "The pair" is the best pair (topk rank 1), whatever model vcfbox calls.
# For each sample, the first rule that holds: depth < MIN_DEPTH drops; a paper
# parent in DISCARDED_LINES drops; the paper pair equal to the pair up to
# twins keeps the paper names (fit_flag marks a single / null call, hom_mm
# above the cut or an ambiguous call); a paper parent outside the panel drops;
# a single / null call or hom_mm above the cut drops; a parent sample (paper
# pair A x A) with a detected cross drops; an ambiguous call drops unless it
# is a pair contest whose runner-up differs from the pair only by twins;
# otherwise the pair (the lone GS3 as father), with a twin of a paper parent
# put back to its paper name, kept only as a GS1/GS2 x GS3 cross or, for a
# parent sample, as the detected line B x B.
# Writes OUT_DIR/new_pedigree.tsv (sample, father, mother, expected_father,
# expected_mother),
# OUT_DIR/new_pedigree_annotated.tsv (every F1 with its features and reason),
# OUT_DIR/panel_distance.tsv.gz (pairwise hom-site discordance of the panel),
# OUT_DIR/calibration.txt.
# The cut is the HOM_MM_QUANTILE of hom_mm in the sample's depth bin over the
# F1s whose paper pair ranks first and is called pair.
import argparse
import os
import subprocess

import numpy as np
import pandas as pd

MIN_DEPTH = 0.2
HOM_MM_QUANTILE = 0.90
DEPTH_BINS = [0.2, 1.0, 3.0, np.inf]
# GS1 and GS2 are male-sterile lines: a cross is a GS1/GS2 mother by a GS3
# father.
MOTHER_GROUPS = {"GS1", "GS2"}
FATHER_GROUP = "GS3"
# Lines discarded by the breeders (segregating); a sample naming one on paper
# or in its final pair drops.
DISCARDED_LINES = {"GS3-110"}


def panel_distance(bcf, cache):
    # Fraction of sites, among those homozygous in both lines, where the two
    # lines carry different alleles; every 10th panel site. A missing call is
    # not homozygous.
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
            rows.append(
                [-1 if "." in gt else int(gt[0]) + int(gt[2]) for gt in line.split()]
            )
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
        " for an ambiguous pair contest, listed in name_ambiguous)",
    )
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    dist = panel_distance(
        args.panel, os.path.join(args.out_dir, "panel_distance.tsv.gz")
    )
    report = pd.read_csv(args.report, sep="\t", dtype={"expected_pair_rank": str})
    if not args.report.endswith(".tsv"):
        raise ValueError(f"{args.report}: expected PREFIX.tsv")
    topk = pd.read_csv(f"{args.report[: -len('.tsv')]}.topk.tsv", sep="\t")
    top = topk[topk["rank"] == 1].set_index("sample")
    second = topk[topk["rank"] == 2].set_index("sample")
    report["pair_a"] = report["sample"].map(top.a)
    report["pair_b"] = report["sample"].map(top.b)
    report["second_a"] = report["sample"].map(second.a)
    report["second_b"] = report["sample"].map(second.b)
    # F1s (dropping a detected selfing) and parent samples (paper A x A).
    paper_self = report.expected_a == report.expected_b
    report["kind"] = np.where(paper_self, "parent", "F1")
    samples = report[
        report.expected_a.notna() & (paper_self | (report.pair_a != report.pair_b))
    ].copy()
    samples["depth_bin"] = pd.cut(samples.coverage, DEPTH_BINS, right=False)
    f1 = samples[samples.kind == "F1"]

    known_good = f1[
        (f1.expected_pair_rank == "1")
        & (f1.best == "pair")
        & (f1.coverage >= MIN_DEPTH)
    ]
    cut = known_good.groupby(
        "depth_bin", observed=False
    ).pair_hom_mismatch_rate.quantile(HOM_MM_QUANTILE)
    samples["hom_mm_cut"] = samples.depth_bin.map(cut).astype(float)

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
        return r.expected_pair_rank != "ABSENT" and same_pair(
            r.pair_a, r.pair_b, r.expected_a, r.expected_b
        )

    def twin_second(r):
        # The runner-up pair differs from the best only by twins.
        return isinstance(r.second_a, str) and same_pair(
            r.pair_a, r.pair_b, r.second_a, r.second_b
        )

    def misfit(r):
        # The reason the pair does not explain the sample, or None.
        if r.best in ("single", "null"):
            return r.best.upper()
        if r.pair_hom_mismatch_rate > r.hom_mm_cut:
            return "HOM_MM_ABOVE_CUT"
        return None

    records = []
    for r in samples.itertuples():
        parent = r.kind == "parent"
        # A pair contest against a twin of the best pair is not ambiguous.
        ambiguous = r.best == "ambiguous" and not (
            r.contest == "pair" and twin_second(r)
        )
        rec = {
            "sample": r.sample,
            "kind": r.kind,
            "paper_mother": r.expected_a,
            "paper_father": r.expected_b,
            "best": r.best,
            "pair_a": r.pair_a,
            "pair_b": r.pair_b,
            "contest": r.contest,
            "depth": r.coverage,
            "unmapped_rate": r.unmapped_rate,
            "hom_mm": r.pair_hom_mismatch_rate,
            "hom_mm_cut": r.hom_mm_cut,
            "expected_rank": r.expected_pair_rank,
            "single_lod_vs_pair": r.single - r.best_pair,
            "expected_lod_vs_best": r.expected - r.best_pair,
        }
        mother = father = None
        flag = misfit(r)
        if r.coverage < MIN_DEPTH:
            reason = "drop_low_depth"
        elif {r.expected_a, r.expected_b} & DISCARDED_LINES:
            reason = "drop_discarded_line"
        elif twin_of_paper(r):
            mother, father = r.expected_a, r.expected_b
            reason = "match" if r.expected_pair_rank == "1" else "twin_paper"
            if flag is not None:
                rec["fit_flag"] = flag
            elif ambiguous:
                rec["fit_flag"] = "AMBIGUOUS"
        elif r.expected_pair_rank == "ABSENT":
            reason = "drop_absent"
        elif flag is not None:
            reason = f"drop_{flag.lower()}"
        elif parent and r.pair_a != r.pair_b:
            reason = "drop_detected_cross"
        elif ambiguous:
            reason = "drop_ambiguous"
        else:
            # The best pair is not the paper pair up to twins, so at least one
            # side is a new (detected) parent.
            det_mother, det_father = (
                (r.pair_a, r.pair_a)
                if parent
                else orient(r.expected_a, r.expected_b, r.pair_a, r.pair_b, is_twin)
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
            f"known-good rows (expected pair rank 1, pair call,"
            f" depth >= {MIN_DEPTH}): {len(known_good)}\n"
        )
        fh.write(f"hom_mismatch_rate cut (quantile {HOM_MM_QUANTILE}) by depth bin:\n")
        fh.write(cut.to_string() + "\n")
    print(open(os.path.join(args.out_dir, "calibration.txt")).read())
    print(out.groupby(["kind", "reason"]).size().to_string())


if __name__ == "__main__":
    main()
