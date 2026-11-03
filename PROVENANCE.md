# Provenance and packaging

## Supplied source

The input was the user-supplied archive `gpap_arxiv.zip`, with SHA-256:

```text
683b5dc708fb20253f3027de678fe5a0918346c539be737b925efe4f081a790c
```

Its `gpap_arxiv/paper.pdf` was readable and used to identify the integral
reproduction instructions in Appendix E.1, the relevant scope in E.3, and the
supporting theory in Appendices B and C. Its SHA-256 is:

```text
b46f1c968c816d516ed3c4a9e9d72ba3d572eb5837ccbffe9a03cfd1cf76076c
```

No claim is made here about the history of earlier private drafts or code not
present in that archive.

## Preserved mathematical files

All nine files under `verification/integral/` are copied **byte for byte** from
the corresponding paths under `gpap_arxiv/` in the supplied ZIP. Their digests
were checked against the supplied `SHA256SUMS.json` before copying and retained
unchanged in this repository's manifest.

The preserved files comprise the compressed rational table, all three entry
points, and the five local geometry/integration support modules. The complete
dependency closure is included, including `full_geometry.py`, which is imported
by a shared module even though its alternate representation is not selected.

Original filenames, docstrings, comments, data keys, and historical usage text
inside those files remain unchanged. In particular, the standalone entry point
of `replay_exact.py` refers to an older data schema. The selected table must be
run through `push027_cross_strict_replay.py`, as explained in the README.

## New packaging files

- `verify_all.py` is a new integral-only runner. The supplied combined runner
  was used as a reference for its integral constants and subprocess invocation;
  it was not retained with hardness branches left active.
- `README.md` and `docs/` document execution, mathematical scope, arrays, and
  every file's role.
- `tests/test_runner.py` checks the new runner's integrity and failure handling.
- The GitHub workflow runs those tests and the full verification; it has not
  been executed on GitHub as part of preparing this local bundle.
- `.gitattributes` preserves LF text and binary gzip handling; `.gitignore`
  excludes generated output and caches.
- `reference_results/` contains a fresh complete run of the included integral
  checks, rather than copies of earlier supplied logs.

The new runner checks the same limiting bound and finite-loss constant as the
supplied combined runner, and substitutes the same two horizon/failure pairs.
It also rejects a weakened manifest, checks coverage and table identity, rejects
optimized mode, refuses to reuse nonempty output directories, and records run
state. No mathematical subcheck was deleted or weakened.

## Excluded material

Hardness data, hardness programs and logs, old combined output summaries,
discovery/search scripts and their archive, duplicate historical integral runs,
PDFs, LaTeX sources, and manuscript build files are excluded. They are not
dependencies of the selected integral verification. The uploaded source archive
itself is unchanged.

## Integrity and licensing

The complete integral pipeline passed locally on Python 3.12.14 on Linux. A
second complete run passed after extracting the generated ZIP to a different
directory and invoking the runner from outside the repository. The eight
failure-handling tests also passed in the original and extracted trees. The
included reference outputs are from the extracted-tree run. The workflow YAML
was parsed locally; its Python 3.10/3.12 GitHub jobs remain to run after upload.

`SHA256SUMS.json` covers the nine original mathematical files and the new root
runner, for 10 entries. It is an integrity record, not a digital signature or a
claim that any candidate table is optimal. Documentation, guard tests, workflow
configuration, and reference logs are not inputs to the mathematical replay.

No code license was present in the supplied files selected for this package.
No new license or ownership claim has been assigned during packaging.
