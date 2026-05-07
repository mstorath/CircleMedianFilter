# Changelog

All notable changes to the `pycirclemedianfilter` Python package are
documented here. This project follows
[Semantic Versioning](https://semver.org/spec/v2.0.0.html) and the
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) layout.

## [Unreleased]

## 0.1.7 — 2026-05-07 (maintenance)

A maintenance release fixing two real bugs and bringing CI / packaging
in line with the rest of the variational-methods family. No algorithm
changes; behaviour-preserving for any existing user code that already
followed the demos' `np.asfortranarray()` idiom.

### Fixed

- **CI verify step now actually verifies.** The post-publish step in
  `build_wheels.yml` ran `pip install … your-package-name` (literal
  placeholder, no relation to this package). `your-package-name`
  happens to exist as someone else's PyPI placeholder, so the step
  silently succeeded across every prior release including 0.1.6.
  Replaced with a real `pip install pycirclemedianfilter` followed by an
  import-and-call smoke check.
- **Default-C-order NumPy arrays now accepted.** The pybind11 bindings
  previously required `np.asfortranarray()` and threw
  `RuntimeError: Input array y must be Fortran contiguous` on default
  inputs. New `py::array_t<double, py::array::f_style | py::array::forcecast>`
  signature auto-converts on entry; both demos and any user code that
  was already passing C-order arrays now Just Work. Output is still
  Fortran-strided. The `np.asfortranarray()` calls in
  `demos_python/CMF_demo*.py` are now redundant but kept for continuity.

### Added

- **`tests_py/` directory with 12 automated tests.** Smoke (module API,
  C-order regression, F-order non-regression, R=1 identity, R-must-be-odd
  validation, quantized variant), invariants (constant-input idempotence,
  global phase shift commutativity, 2π periodicity, hand-checked 3×3
  cluster-vs-outlier median), and demo runners (both bundled demos
  exercised end-to-end under `MPLBACKEND=Agg`).
- **CI test step** runs `pytest tests_py/` against the freshly-built
  Linux wheel.
- **`CHANGELOG.md`** (this file).
- **`CITATION.cff`** at root with the IEEE TPAMI 2018 reference
  (`10.1109/TPAMI.2017.2692779`) and Storath + Weinmann ORCIDs.

### Changed

- **Packaging migrated to PEP 621.** Project metadata (name, version,
  authors, classifiers, keywords, deps, urls) lives in
  `pyproject.toml [project]`. `setup.py` is reduced to the C++ extension
  declaration. Single source of truth.
- **CI workflow modernised.**
  - `actions/checkout@v3` → `@v4`, `actions/setup-python@v4` → `@v5`.
  - Matrix dropped `macos-13` (Intel; family decision); added `cp313`.
  - Linux `cmake` apt install removed (this isn't a CMake build).
  - `strategy.fail-fast: false` so platform-specific failures don't
    cancel siblings.
  - **Publish job is now tag-gated.** Previous behaviour: every push to
    master ran `twine upload --skip-existing` using `PYPI_API_TOKEN`.
    New behaviour: publish runs only on `v*` tags via PyPI OIDC trusted
    publishing (no long-lived secret in the workflow).
- **`.gitignore`** rewritten from a 24-line list of specific old
  filenames into a glob-based layout. Adds `*.mex*` for future
  contributor-machine builds; existing committed Mex binaries
  (`filters/CMF_medfilt*.mex{maca64,maci64}`) explicitly retained via
  `!filters/...` re-includes.

### Security

- Long-lived `PYPI_API_TOKEN` is no longer used by the workflow. Per
  user decision the secret is **retained** in repo settings as a manual
  fallback (`twine upload` from a local terminal if OIDC ever breaks),
  but the published-by-OIDC posture is now the default.

### One-time user-side prerequisites for the next v* tag to publish

1. PyPI dashboard: register a Trusted Publisher at
   <https://pypi.org/manage/account/publishing/> with
   project=`pycirclemedianfilter`, owner=`mstorath`,
   repo=`CircleMedianFilter`, workflow=`build_wheels.yml`,
   environment=`pypi`.
2. GitHub: create a `pypi` environment in the repo's
   Settings → Environments. Optional: required-reviewer gate.

## 0.1.6 — 2025-02-18

First Python release; published as `pycirclemedianfilter` on PyPI. C++
algorithms (Storath & Weinmann TPAMI 2018) made available via pybind11
bindings; `medfilt_circ2d` and `medfilt_circ2d_quant` exposed.
