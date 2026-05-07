# Basic smoke tests for the pycirclemedianfilter pybind11 module.
# Confirms the public API is callable, accepts default-C-order arrays
# (regression test for the v0.1.7 f_style|forcecast fix), and returns
# correctly-shaped output.

import numpy as np
import pytest

import pycirclemedianfilter as cmf


def test_module_exposes_expected_functions():
    public = {a for a in dir(cmf) if not a.startswith("_")}
    assert {"medfilt_circ2d", "medfilt_circ2d_quant"}.issubset(public)


def test_default_c_order_input_accepted(small_phase_image):
    """Regression test: prior to v0.1.7 a non-Fortran-contiguous input
    raised a confusing RuntimeError. The f_style|forcecast pybind11 flag
    auto-converts on entry so this should now just work."""
    assert small_phase_image.flags.c_contiguous
    out = cmf.medfilt_circ2d(small_phase_image, 1, 1)
    assert out.shape == small_phase_image.shape
    assert out.flags.f_contiguous  # output is always Fortran-strided


def test_fortran_order_input_still_works(small_phase_image):
    y_f = np.asfortranarray(small_phase_image)
    out_f = cmf.medfilt_circ2d(y_f, 1, 1)
    out_c = cmf.medfilt_circ2d(small_phase_image, 1, 1)
    np.testing.assert_array_equal(out_c, out_f)


def test_filter_size_one_is_identity(small_phase_image):
    """R=1, T=1 means a 1×1 window — no neighbours, output equals input.
    (R, T are full filter side lengths, must be odd; not radii.)"""
    out = cmf.medfilt_circ2d(small_phase_image, 1, 1)
    np.testing.assert_allclose(out, small_phase_image)


def test_filter_size_must_be_odd(small_phase_image):
    """R or T = 2, 4, … is rejected by the underlying C++."""
    with pytest.raises(ValueError, match="filter size must be odd"):
        cmf.medfilt_circ2d(small_phase_image, 2, 1)
    with pytest.raises(ValueError, match="filter size must be odd"):
        cmf.medfilt_circ2d(small_phase_image, 1, 2)


def test_quantized_variant(small_phase_image):
    v = np.linspace(-np.pi, np.pi, 16)
    out = cmf.medfilt_circ2d_quant(small_phase_image, 1, 1, v)
    assert out.shape == small_phase_image.shape
    # Output values should be drawn from the quantization set v.
    out_unique = np.unique(out)
    for val in out_unique:
        assert np.any(np.isclose(v, val, atol=1e-12)), (
            f"unquantized value {val} not in v"
        )
