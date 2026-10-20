# List label problems found in a `vcfbox parentage` report, one row per
# sample and issue.
#
# usage: uv run --with pandas python scripts/label_issues.py \
#            f1.parentage.tsv f1_pedigree/panel_distance.tsv.gz OUT.tsv
#
# Issues:
#   tube_shift         in a run of >= 3 samples, each detected pair is the
#                      record of the previous (N-1) or next (N+1) sample
#   neighbour_swap     samples N and N+1 carry each other's records
#   line_substitution  an F1 whose paper line is replaced by the line that
#                      replaces it in most of its F1s (LINE_SUBSTITUTIONS)
#   parent_mislabel    a parent sample (paper A x A) detected as B x B
# Only samples with depth >= MIN_DEPTH and a PASS / SELFING verdict count.
import sys

import pandas as pd

MIN_DEPTH = 0.2
# Paper line -> line found instead in most of its F1s; the reverse does not
# happen, so these are not panel label swaps.
LINE_SUBSTITUTIONS = {"GS3-110": "GS3-081-2", "GS3-063-3": "GS3-064",
                      "GS3-174-2": "GS3-147", "GS3-186-3": "GS3-147",
                      "GS3-187": "GS3-186-3", "GS3-026": "GS3-117"}
MIN_RUN = 3


def pair_text(pair):
    return " x ".join(sorted(pair)) if pair else "NA"


def main():
    report_path, distance_path, out_path = sys.argv[1:4]
    r = pd.read_csv(report_path, sep="\t")
    r["num"] = pd.to_numeric(r["sample"].str[2:], errors="coerce")
    r["best"] = [frozenset((a, b)) for a, b in zip(r.best_a, r.best_b)]
    r["paper"] = [frozenset((a, b)) if isinstance(a, str) else None
                  for a, b in zip(r.expected_a, r.expected_b)]
    paper = dict(zip(r.num, r.paper))
    dist = pd.read_csv(distance_path, sep="\t", index_col=0)
    ok = r[(r.depth_mean >= MIN_DEPTH) & r.verdict.isin(["PASS", "SELFING"])]
    ok = ok.set_index("num", drop=False)

    rows = []

    def add(x, issue, detail):
        rows.append({"sample": x["sample"], "issue": issue,
                     "paper_pair": pair_text(x.paper),
                     "detected_pair": pair_text(x.best),
                     "verdict": x.verdict, "depth": x.depth_mean,
                     "detail": detail})

    def neighbour_hits(k):
        return sorted(n for n, x in ok.iterrows()
                      if x.best != x.paper and paper.get(n + k) == x.best)

    explained = set()
    for k in (-1, 1):
        hits = neighbour_hits(k)
        run = []
        for n in hits + [None]:
            if run and (n is None or n - run[-1] > 2):
                if len(run) >= MIN_RUN:
                    for m in run:
                        add(ok.loc[m], "tube_shift",
                            f"run SH{run[0]:05d}-SH{run[-1]:05d}: carries "
                            f"the record of SH{m + k:05d}")
                        explained.add(m)
                run = []
            if n is not None:
                run.append(n)

    for n in neighbour_hits(1):
        if n + 1 in ok.index and n + 1 in neighbour_hits(-1):
            for m, o in ((n, n + 1), (n + 1, n)):
                add(ok.loc[m], "neighbour_swap",
                    f"carries the record of SH{o:05d}")
                explained.add(m)

    for n, x in ok.iterrows():
        if n in explained or x.paper is None or len(x.paper) != 2:
            continue
        for line, sub in LINE_SUBSTITUTIONS.items():
            if line in x.paper and sub in x.best \
                    and x.paper - {line} <= x.best:
                add(x, "line_substitution",
                    f"{line} found as {sub} (panel distance "
                    f"{dist.loc[line, sub]:.4f})")

    for n, x in ok.iterrows():
        if n in explained or x.paper is None or len(x.paper) != 1 \
                or x.expected_status == "ABSENT":
            continue
        if len(x.best) == 1 and x.best != x.paper:
            (a,) = x.paper
            (b,) = x.best
            d = dist.loc[a, b]
            add(x, "parent_mislabel",
                f"parent sample of {a} is line {b} (panel distance {d:.4f})")

    out = pd.DataFrame(rows).sort_values(["issue", "sample"])
    out.to_csv(out_path, sep="\t", index=False, float_format="%.3f")
    print(out.issue.value_counts().to_string())


if __name__ == "__main__":
    main()
