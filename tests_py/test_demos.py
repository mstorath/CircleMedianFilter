# Run the bundled demos end-to-end as regression tests. Catches breakage
# in the public-facing example code, which is the first thing a new user
# touches after `pip install`.

import importlib.util
import os
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
DEMOS_DIR = REPO_ROOT / "demos_python"


def _import_demo(filename):
    path = DEMOS_DIR / filename
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(autouse=True)
def _matplotlib_headless(monkeypatch):
    """Force matplotlib's non-interactive Agg backend so the demos render
    figures into memory instead of opening a GUI."""
    monkeypatch.setenv("MPLBACKEND", "Agg")
    yield


@pytest.mark.skipif(
    not (DEMOS_DIR / "data" / "CMF_imgInSAR.png").exists()
    and not (REPO_ROOT / "data" / "CMF_imgInSAR.png").exists(),
    reason="demo data file CMF_imgInSAR.png not present",
)
def test_demo_unquantized_runs():
    mod = _import_demo("CMF_demoUnquantized.py")
    # Demos chdir into demos_python/ relative paths; make sure cwd matches.
    cwd = os.getcwd()
    try:
        os.chdir(DEMOS_DIR)
        mod.main()
    finally:
        os.chdir(cwd)


@pytest.mark.skipif(
    not (DEMOS_DIR / "data" / "CMF_imgInSAR.png").exists()
    and not (REPO_ROOT / "data" / "CMF_imgInSAR.png").exists(),
    reason="demo data file CMF_imgInSAR.png not present",
)
def test_demo_quantized_runs():
    mod = _import_demo("CMF_demoQuantized.py")
    cwd = os.getcwd()
    try:
        os.chdir(DEMOS_DIR)
        mod.main()
    finally:
        os.chdir(cwd)
