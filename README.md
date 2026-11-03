# Integral knapsack and GPAP: exact verification

This repository contains the rational rule table and all programs needed to
reproduce the **integral** computer-assisted bound in *Improved Algorithms in the
Random-Order Model: Knapsack, GAP, and Beyond*, by Saksham Agrawal, Arindam Khan,
and Akash Pareek.

The computation checks the numerical hypotheses used for the paper's
**0.2564-competitive integral algorithm**, measured against the offline integral
optimum. It also checks the stronger finite coefficient **0.25647580** with the
larger augmented horizon specified below. Appendix B of the supplied manuscript
proves the one-resource construction and its finite-input analysis; Appendix C
extends it to integral GPAP using resource deletion.

The three integral checks described in Appendix E.1 are all included. This
repository contains no hardness verification, parameter-search archive, or
manuscript build files. Section references below refer to the supplied version
of the manuscript and may change when its appendices are renumbered.

## Quick start

1. Download and extract this repository, or clone it from GitHub.
2. Open a terminal in the directory containing `verify_all.py`.
3. Run:

   ```sh
   python -S verify_all.py
   ```

Use **Python 3.10 or newer**. If your system calls it `python3`, substitute
`python3` in the commands. On Windows, `py -3 -S verify_all.py` is an alternative
when Python is installed through the Python launcher.

There are **no packages to install**: all verification uses Python's standard
library. No optimizer, solver, network access, or random seed is required. The
compressed input is read directly; do not decompress or regenerate it.

**Do not use `-O`, `-OO`, or a nonzero `PYTHONOPTIMIZE`.** The separate spatial
audit uses assertions for mathematical checks. The root runner refuses optimized
mode rather than allowing these checks to disappear.

The programs perform hundreds of thousands of exact rational checks. Allow
several minutes, depending on your machine. The runner prints progress when each
check starts and finishes; a check can run silently for a while. Measured times
for the included run are in [the reference summary](reference_results/all_checks.json).

The successful ending is:

```text
ALL INTEGRAL CHECKS PASS: finite ratio > 0.2564.
Stronger finite check also passes: ratio > 0.25647580.
```

Success requires **exit code 0** and `"status": "PASS"` in the newly produced
`all_checks.json`. A previously supplied result file is not a fresh verification.

By default, output goes into `results/` beside `verify_all.py`. The directory must
be new or empty. To retain earlier results and run again:

```sh
python -S verify_all.py --output results-second-run
```

A relative `--output` path is resolved against the terminal's current directory.
The script locates all code and input data relative to its own location, so it
also works when invoked by its full path from another directory.

## What the runner does

Before executing mathematical code, `verify_all.py` checks all **10** entries in
`SHA256SUMS.json`: the nine preserved integral files and the new root runner.
It rejects missing files, changed bytes, and missing or extra manifest entries.
It separately checks the selected table's fixed SHA-256. The integrity check is
repeated before reporting success.

The runner then executes these three programs, in order:

| Stage | Entry point | What it establishes |
| --- | --- | --- |
| Full replay | `verification/integral/research/push027_cross_strict_replay.py` | Reconstructs the required index sets, verifies target tails, availability, crossing envelopes, exact packing maxima, quadratic lower functions, rational integral bounds, and finite-loss constants. |
| Spatial audit | `verification/integral/research/eplus1_audit_true_spatial.py` | Recomputes crossing gaps at actual boundaries and on open polynomial pieces, including one-sided limits at atom jumps. |
| Independent value audit | `verification/integral/verify_service_independent.py` | Recomputes all 38 endpoint value coefficients and the finite correction using individual companion subsets and positive Taylor/Beta integrals. It imports none of the research modules. |

The runner checks agreement on the selected table, the claimed exact constants,
the expected coverage counts, and the independent value bounds. Finally it
checks both finite guarantees using rational arithmetic.

Any failed subprocess or failed comparison makes the runner exit with a nonzero
status. A fresh successful summary is written only after every check succeeds.
`run_status.json` records `RUNNING`, then `PASS` or a caught `FAIL`; an interrupted
or forcibly terminated run must not be treated as successful.

## Mathematical context

For one resource, the paper studies an integral knapsack of capacity one. Items
arrive in uniformly random order. Each item has a nonnegative value and a size
in `[0,1]`; an accepted positive-size item is taken in full.

The online construction maintains an actual load and a virtual load that is at
least as large. A rule can accept an item only where it fits in the virtual
knapsack. Its two relevant quantities are its total acceptance probability and
its probability of moving the virtual load past a specified threshold.

A prefix LP chooses among a finite menu of rules. Its objective rewards
acceptance, while its threshold constraints limit the virtual-load drift.
Conditional on an unordered prefix of `k` items, each item is equally likely to
be current. This divides the prefix's total crossing cost by `k`. Upward padding
can then restore the prescribed virtual-load distribution. These online and
conditioning arguments are proved in the manuscript.

To compare the prefix LP with any offline integral packing, the analysis assigns
eligible rules to the observed members of that packing. The computer checks
that these assignments fit the crossing budget and yield a sufficiently large
coefficient for every distinguished item. The table uses:

| Parameter | Value | Meaning |
| --- | --- | --- |
| `K` | `6` | Size and load grid spacing `1/6`. |
| `h` | `2` | Small items have size at most `h/K = 1/3`; a packing has at most two larger items. |
| `factor` | `4` | Observed small-mass grid spacing `1/(K*factor) = 1/24`. |
| `M` | `16` | Number of time blocks. |
| `T` | `17/4` | Length of the continuous construction interval. |
| `denominator` | `100000000000` | Common denominator `Q` of the integer coefficient arrays. |

Target tails are cubic in time. Accepted masses, crossing envelopes, and lower
functions are quadratic in time. The analysis reduces the continuous
requirements to finitely many exact inequalities; it does not infer them by
sampling a grid. [Proof checks](docs/proof-checks.md) explains the reductions and
maps the five proof obligations to the programs.

## Expected coverage and bounds

The full replay must report:

| Family | Count |
| --- | ---: |
| Full large-item patterns | 7 |
| Comparison labels | 146 |
| Label/class rules | 290 |
| Availability checks | 226,688 |
| Crossing checks | 723,840 |
| Packing checks | 219,072 |
| Quadratic lower-function checks | 47,760 |
| Endpoint value coefficients | 38 |

The separate spatial audit performs **473,280** corner and quadratic checks,
**125,280** additional polynomial-identity checks, and records **1,914** jump
endpoint occurrences. Its minimum margin is zero. This is acceptable because
the required spatial inequalities are non-strict.

The exact limiting comparison coefficient is:

```text
gamma_star = 307770970018984529425873353861053278824139 /
             1200000005486015625000000000000000000000000
           = 0.256475807176628536517001723089915623... 
```

The finite-loss constant is:

```text
A_fin = 4718849996128173 / 17000000000
```

Here `A_fin` is a scalar from the finite-input analysis, distinct from the data
array named `A`. The finite comparison is
`gamma_star - A_fin/N - epsilon_samp`, with the following two checks:

| Minimum augmented horizon `N` | Global sampling-failure allowance | Verified lower coefficient, shown approximately | Strict target |
| --- | --- | --- | --- |
| `10^10` | `10^-8` | `0.25644803923547484138...` | `> 0.2564` |
| `10^14` | `10^-14` | `0.25647580440082442114...` | `> 0.25647580` |

For an input of length `n`, the paper uses `N >= max(n, N_min)` by inserting
zero-value dummy items. These two checks are algebraic substitutions into the
proved finite-loss bound; **the programs do not simulate `10^10` or `10^14`
rounds**. The sampling-failure allowance is global across all resources.

The value `gamma_star` is a limiting coefficient, not itself the stated
finite-input guarantee. The ratio `0.2564` has reciprocal factor
`2500/641 = 3.900156...`, approximately `3.9002`.

## Files and dependencies

All paths in the next table are relative to `verification/integral/`.

| File | Purpose |
| --- | --- |
| `research/runs/eplus1_search_expanded5_cert.json.gz` | The complete selected rational table: tails, accepted masses, crossing bounds, quadratic lower functions, reserves, and claimed constants. This is the only mathematical input data file. |
| `research/push027_cross_strict_replay.py` | Entry point for the full replay of the selected table. |
| `research/eplus1_audit_true_spatial.py` | Standalone spatial audit, including closed boundaries and one-sided open-segment quadratics. |
| `verify_service_independent.py` | Standalone independent value and finite-loss audit. Its default input is the selected table. |
| `research/exact_geometry.py` | Exact crossing expressions and the seven subdivided availability controls. The selected full replay imports its availability controls. |
| `research/push027_cross_geometry.py` | Crossing-envelope checks with subdivided diagonal segments for the selected construction. |
| `research/push027_strict_geometry.py` | Reconstructs labels and eligible pairs, computes strict-boundary packing maxima, and generates lower-function test cases. |
| `research/replay_exact.py` | Provides shared pattern enumeration, exact exponential bounds, signed interval sums, and integral weights. Its standalone replay targets an older schema; do not use that entry point for the selected table. |
| `research/full_geometry.py` | Geometry for an alternative size-affine-density representation. Imported by `replay_exact.py`, so retained even though that alternate representation is not selected here. |

The full replay imports four local modules: `exact_geometry`,
`push027_cross_geometry`, `push027_strict_geometry`, and `replay_exact`.
`replay_exact` also imports `full_geometry`. The other two entry points use only
the standard library. Keeping this complete import closure is essential.

Other repository files:

| File or directory | Purpose |
| --- | --- |
| `verify_all.py` | New integral-only orchestrator, integrity checks, exact cross-checks, finite substitutions, and result collection. |
| `SHA256SUMS.json` | SHA-256 digests of the nine preserved files and the new root runner. |
| `docs/proof-checks.md` | Mathematical role of each check, continuous-domain reductions, and finite-horizon constants. |
| `docs/table-format.md` | Complete array layout, index conventions, scaling, and meaning of the stored entries. |
| `PROVENANCE.md` | Source archive identity, preservation policy, packaging changes, and exclusions. |
| `reference_results/` | Fresh outputs and logs from the full run used to check this package. See the file table below. |
| `tests/test_runner.py` | Standard-library tests of failure handling: missing or changed files, manifest omissions, optimized mode, and output reuse. |
| `.github/workflows/verify.yml` | Runs the guard tests and all three integral checks on GitHub, then uploads the fresh results as a workflow artifact. |
| `.gitattributes` | Keeps Python/JSON text at LF line endings and treats the compressed table as binary, preserving hashes across checkouts. |
| `.gitignore` | Excludes generated local results and Python caches. |

[Table format](docs/table-format.md) also explains why stored metadata such as
`summary.status` and displayed decimal approximations are not accepted as proof.
Historical source filenames and output keys are retained to preserve the
mathematical programs byte for byte.

## Result files

Every successful run produces the following files in its output directory:

| File | Contents |
| --- | --- |
| `all_checks.json` | Combined status, interpreter/platform, table hash, coverage, exact constants, finite guarantees, and elapsed time. |
| `run_status.json` | Overall run state and UTC timestamps. |
| `integral_full.json` | Full replay result, counts, limiting bound, worst endpoint case, coefficient norm, finite-loss constant, and runtime. |
| `integral_spatial.json` | Spatial coverage, jump counts, exact minimum margin, method, table hash, and runtime. |
| `integral_service.json` | Independent value result, rational lower/upper intervals for all 38 endpoint cases, finite constants, and runtime. |
| `integral_full.stdout.log`, `integral_full.stderr.log` | Raw full-replay output and errors. |
| `integral_spatial.stdout.log`, `integral_spatial.stderr.log` | Raw spatial-audit output and errors. |
| `integral_service.stdout.log`, `integral_service.stderr.log` | Raw independent-value output and errors. |

`reference_results/` contains this same layout. These files are evidence of a
completed run, not input to a subsequent run. Time measurements and platform
metadata can differ; exact coefficients and coverage should agree.

## Running individual checks

Use the root runner for a complete verification. For inspection or debugging,
the three mathematical programs can be run separately from the repository root:

```sh
python -S verification/integral/research/push027_cross_strict_replay.py verification/integral/research/runs/eplus1_search_expanded5_cert.json.gz
python -S verification/integral/research/eplus1_audit_true_spatial.py verification/integral/research/runs/eplus1_search_expanded5_cert.json.gz spatial.json
python -S verification/integral/verify_service_independent.py verification/integral/research/runs/eplus1_search_expanded5_cert.json.gz --output service.json
```

The first prints one JSON result to stdout. The second writes `spatial.json`
and prints its result. The third writes `service.json`, including all endpoint
intervals, and prints a shorter summary. Individual invocations bypass the root
runner's manifest checks, agreement checks, and its combined two-horizon check.
Do not treat one individual `PASS` as a complete verification.

The guard tests are separate from the mathematical checks:

```sh
python -B -S -m unittest discover -s tests -v
```

They run quickly and do not replace the full replay. The included GitHub workflow
executes both this command and `python -B -S verify_all.py`.

## Troubleshooting

| Message or symptom | Action |
| --- | --- |
| Python executable not found | Install Python 3.10 or newer, then try `python3` or the Windows `py -3` launcher. |
| Optimized mode rejected | Remove `-O`/`-OO` and unset `PYTHONOPTIMIZE` or set it to `0`. |
| Missing file or local geometry module | Extract or clone the complete tree, retaining its directories. Do not copy just the entry points. |
| SHA-256 mismatch | Restore the original file from this bundle or the verified commit. Avoid editing or reformatting the preserved code or compressed data. |
| Manifest file-set error | Restore `SHA256SUMS.json`; omitting an entry is not a way to skip a check. |
| Output directory must be new or empty | Choose a fresh path with `--output results-another-run`. Previous outputs are deliberately preserved. |
| A check appears silent | Let it finish; exact rational calculations are CPU-intensive. Inspect the corresponding logs if it exits unsuccessfully. |
| A mathematical check fails | Keep the logs and nonzero exit status. Do not remove the failed inequality or change the constants to obtain `PASS`. |
