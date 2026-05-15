# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

```sh
pixi run build      # produces .build/vcfbox
```

htslib/cmake/ninja/compiler come from `.pixi/envs/default` — always go through pixi. Submodules `ext/CLI11` and `ext/eigen` must be initialized. No test framework; verify against fixtures in `data/`.

## CLI

Single binary with subcommands:

- `combine` — implemented; merges paired samples into pseudo-homozygous samples named `A~B`.
- `parentage build-matrix`, `parentage test` — stubs that throw `not implemented yet`.

## Architecture

- `src/app/` owns CLI11 wiring; never includes htslib.
- `src/combine/`, `src/parentage/` are pure business logic; never include CLI11.
- `src/app/subcommands.h` (Options structs + `add_*` / `run_*` decls) is the only contract between the two layers.
- `src/hts/hts_raii.h` provides `unique_ptr` deleters for `htsFile` / `bcf_hdr_t` / `bcf1_t` and a `Genotypes` RAII for the htslib-`malloc`'d GT buffer. Always use these.
- Business code throws `std::runtime_error`; `cli.cpp`'s `run_command` catches and converts to exit code. Don't catch in business logic.

### Adding a subcommand

1. Add `Options` struct + `add_*` / `run_*` decls to `src/app/subcommands.h`.
2. Create `src/app/<name>_command.cpp` mirroring `combine_command.cpp`.
3. Put logic in `src/<name>/` (no CLI deps).
4. Append both `.cpp` files to `add_executable(vcfbox …)` in `CMakeLists.txt`.
5. Wire `add_*` + dispatch branch in `src/app/cli.cpp`.
