# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

```sh
pixi run build      # produces .build/vcfbox and .build/tests/vcfbox_tests
pixi run test       # builds, then runs the Catch2 suite
```

htslib/CLI11/Catch2/cmake/ninja/compiler come from `.pixi/envs/default` — always go through pixi. `data/` is gitignored and holds real fixtures for manual checks only; tests never read it.

## Tests

- Layout mirrors `src/`: `tests/parentage/test_panel.cpp` tests `src/parentage/panel.*`. Shared helpers live in `tests/` (`file_fixture.h` for a temp dir, `vcf_fixture.h` / `sam_fixture.h` for literal VCF / SAM text; htslib reads plain `.vcf` / `.sam`, so no binary fixtures).
- One behavior per `TEST_CASE`/`SECTION`, named after the action and expected outcome, laid out as given/when/then.
- Cases of the same behavior (e.g. different boundaries of one rule) are `SECTION`s of one `TEST_CASE`, not separate `TEST_CASE`s.
- Test through the public API and assert on results, not on calls made. A test changes only when behavior changes. Use real objects; no mocks.
- Hermetic: no network, sleeping, or shared mutable state; I/O only through fixtures in the test tree.
- No logic: no loops, branches, or computed expectations. Write inputs and expected values as literals; use `GENERATE` instead of loops. Generation that needs a loop belongs in a fixture helper.
- DAMP over DRY. Shared setup holds only values the test does not care about.
- Diagnosable failures: matchers for scalars (`ContainsSubstring`, `WithinAbs`), `CAPTURE` before compound comparisons.
- Flaky tests are fixed or deleted, never retried.

## CLI

Single binary with subcommands:

- `combine` — implemented; merges paired samples into pseudo-homozygous samples named `A~B`.
- `parentage` — being rewritten; `docs/parentage_design.md` is the design (purpose, data structures, plan), `docs/parentage_review.md` is the review of the old module. Implemented and tested so far: panel loading (`src/parentage/{bitrows.h,sites,panel}`), pileup + counts files (`src/parentage/{counts,pileup,counts_io}`) and scoring (`src/parentage/score`: single-parent model first, integer pair kernel over the `--candidates` best single parents, float stats of the top pairs); `parentage.cpp` runs them per sample on a `--threads` worker pool and prints scores to stderr until the report stage exists.

## Naming

- Types (class/struct/enum): `PascalCase`
- Functions/variables/files: `snake_case`
- Private instance members: trailing underscore (`member_`)
- `const`/`constexpr` objects, including static members: `snake_case`
- Macros: `UPPER_SNAKE_CASE`

## Design

- Each abstraction layer must encapsulate exactly one independent axis of change. If every implementation of a layer looks the same, collapse that layer; if one layer has to carry two axes of change, split it.
- Use `class` when a type independently owns runtime invariants: keep its representation private and enforce validity through construction and controlled mutation. Use `struct` only for transparent aggregates, compile-time traits, or internal algorithm state whose invariants are explicitly owned elsewhere.
- Templates and concepts constrain types and compile-time values; they do not replace validation of runtime numeric invariants.

### Construction APIs

- A constructor accepts the final representation or retained dependencies and validates every invariant determined by that stored state. Delegating constructors only supply defaults or normalize arguments — no separate validation, I/O, or policy.
- Use a namespace-level free factory when construction involves anything not retained by the result: source objects from another abstraction, combining independent owners, relationships that need validating but aren't stored, or a deduced/selected return type. If only a factory can establish validity, make the constructor non-public instead of duplicating checks.
- Naming: no generic static `T::make(...)`. `make_*` free functions for in-memory construction, `open_*`/`read_*` for I/O, type-local statics only for named presets like `defaults()`.

## Code Style

- Only comment when the logic is not obvious; names and signatures should carry the intent
- A comment states what the thing is or does, plus the contract facts the signature cannot express: what a returned value or error means, what is mutated, when it fails. No derivations or rationale; cite a source instead of explaining it.
- Comments use `//` only, never `/* */`
- Use `noexcept` only when the function contract is unconditionally non-throwing
- `[[nodiscard]]` only where discarding the result is almost certainly a bug: factories and `open_*`/`read_*`/`make_*`, computations whose return value is the sole product, and boolean queries easily mistaken for actions (`empty()`, `contains()`). Never on plain getters.
- No `using namespace` in headers. Avoid namespace-scope aliases in public headers unless they are part of the API.
- Trailing return types: `auto f() -> int`
- Non-owning inputs: `std::span`, `std::string_view`
- Explicit `#include` only — no transitive includes
- Do not introduce extra namespaces without a clear architectural need
- Simple getters/setters stay in headers
- Prefer constructor member initializer lists.
- Prefer `std::views::enumerate` over manual index loops when both index and value are needed
- When designing the API, adhere to STL conventions.
- `std::format` for every string built from parts; no `+` concatenation.
- Readable arithmetic over bit tricks: `bit / 64` and `bit % 64`, not `>> 6` and `& 63`.
- No speculative defence: the tool is used in-house. A single cheap check is the ceiling; ask "what is this defending against?" before adding one.
- Fail loudly on bad input (throw with a message naming contig, position, sample and value); never silently skip or coerce.
- Reply in Chinese; keep identifiers and technical terms as-is.

## Parentage data facts (don't re-derive)

- Parents are inbred lines with residual heterozygosity in chromosome blocks (real panel: 0.26% of GTs het, median parent 0.15%, max 3.3%; 40% of sites have some het parent). `Panel` carries a `het` bit plane; het sites contribute nothing to the integer ranking and are marginalised exactly in the float stage. Missing/haploid/multi-allelic in the panel is an error.
- Real panel: `data/f0.02/parents.len.bcf` (3.17M sites, 307 parents; from gpu:~/project/imputation_bench/panel/f0.02, concatenated and reheadered with contig lengths from the BAM). `sites.chrNN.tsv` there holds the per-site opposite-read fraction from clean samples.
- Samples are F1s; mother/father are symmetric, unordered pairs, diagonal = selfing.
- Real panel is ~3M sites; `data/panel_filled.vcf.gz` is a 47k-site fixture made by `bcftools +setGT ... -t . -n c:'0/0'` (its 202-site 080/111 difference is a fill artefact).
- Scoring is integer: 1-bit panel rows, AND/XOR + popcount over count bit-planes; no BLAS/gemm. Oracle: numpy float matrix products over the `bcftools query` genotype matrix agree on every count, score and LOD for both samples. On the 47k fixture the true pair of SH12277 shows hom_mismatch 0.038 (0.023 on sites both parents had called in `data/sites.parentage.vcf.gz`, 0.38 on `setGT`-filled sites); calibrate verdict thresholds on the real panel, not this fixture.
- Pileup runs on htslib's `bam_mplp_*` engine (what samtools mpileup uses) with samtools' read filters, so `samtools mpileup -B -q 20 -Q 20 --ff UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY` matches `vcfbox parentage` exactly per site, and `-x` matches `--keep-overlaps`. Non-proper pairs are dropped (samtools default). Only deliberate difference: MAPQ 255 is dropped (samtools keeps it); these BAMs have none.
- `data/SH12277.bam` (1.1 GB) + `data/SH17265.bam` (1.8 GB) pile up in ~39 s single-threaded against the 47k fixture (a hand-written CIGAR walk was 3x faster but not bit-identical to samtools; rejected).

## Architecture

- `CMakeLists.txt` builds `vcfbox_core` (everything under `src/` except `src/app/` and `main.cpp`) as a static library; the `vcfbox` binary and `tests/vcfbox_tests` both link it. New business `.cpp` files go into `vcfbox_core`, new command files into `vcfbox`.
- `src/app/` owns CLI11 wiring; never includes htslib.
- `src/combine/`, `src/parentage/` (and future `src/<name>/` modules) are pure business logic; never include CLI11.
- `src/app/subcommands.h` (Options structs + `add_*` / `run_*` decls) is the only contract between the two layers.
- `src/hts/hts_raii.h` provides `unique_ptr` deleters for `htsFile` / `bcf_hdr_t` / `bcf1_t` and a `Genotypes` RAII for the htslib-`malloc`'d GT buffer. Always use these.
- Business code throws `std::runtime_error`; `cli.cpp`'s `run_command` catches and converts to exit code. Don't catch in business logic.

### Adding a subcommand

1. Add `Options` struct + `add_*` / `run_*` decls to `src/app/subcommands.h`.
2. Create `src/app/<name>_command.cpp` mirroring `combine_command.cpp`.
3. Put logic in `src/<name>/` (no CLI deps).
4. Add the business `.cpp` to `vcfbox_core` and the command `.cpp` to `vcfbox` in `CMakeLists.txt`; tests under `tests/<name>/test_*.cpp` are picked up by glob.
5. Wire `add_*` + dispatch branch in `src/app/cli.cpp`.
