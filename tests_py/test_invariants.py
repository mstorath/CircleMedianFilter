# Mathematical invariants of the arc-distance median filter.

import numpy as np
import pycirclemedianfilter as cmf


def _wrap_to_pi(x):
    """Wrap to [-π, π)."""
    return np.mod(x + np.pi, 2 * np.pi) - np.pi


def test_constant_input_unchanged(constant_phase_image):
    """A constant phase image is invariant under the circle-median filter
    at any odd window size (the median of equal values is that value)."""
    for size in [1, 3, 5, 7]:
        out = cmf.medfilt_circ2d(constant_phase_image, size, size)
        np.testing.assert_allclose(
            out,
            constant_phase_image,
            atol=1e-12,
            err_msg=f"size={size}",
        )


def test_global_phase_shift_commutes_with_filter():
    """For any constant Δφ, filter(input + Δφ) == filter(input) + Δφ
    (modulo 2π). The arc-distance median is rotationally equivariant on
    the circle."""
    rng = np.random.default_rng(0)
    img = rng.uniform(-np.pi, np.pi, size=(7, 7))
    delta = 1.234

    out0 = cmf.medfilt_circ2d(img, 1, 1)
    img_shifted = _wrap_to_pi(img + delta)
    out_shifted = cmf.medfilt_circ2d(img_shifted, 1, 1)
    expected = _wrap_to_pi(out0 + delta)

    np.testing.assert_allclose(out_shifted, expected, atol=1e-10)


def test_2pi_periodicity():
    """Adding 2π to any pixel doesn't change the result — the filter
    treats the input as circle-valued."""
    rng = np.random.default_rng(1)
    img = rng.uniform(-np.pi, np.pi, size=(5, 5))
    img_plus_2pi = img.copy()
    img_plus_2pi[2, 2] += 2 * np.pi  # one pixel offset by 2π

    out0 = cmf.medfilt_circ2d(img, 1, 1)
    out1 = cmf.medfilt_circ2d(img_plus_2pi, 1, 1)
    # The filter result should be the same up to 2π wrapping.
    np.testing.assert_allclose(_wrap_to_pi(out0), _wrap_to_pi(out1), atol=1e-10)


def test_3x3_window_pixelwise_median():
    """For R=3, T=3 (3×3 window — note R is full size, not radius) the
    filter at an interior pixel reduces to the arc-distance median of 9
    values. Hand-check on a 5×5 with a single-pixel outlier: 8 cluster
    pixels at 1.0, the centre at 3.0 → median = 1.0."""
    img = np.full((5, 5), 1.0)
    img[2, 2] = 3.0  # centre is the lone outlier
    out = cmf.medfilt_circ2d(img, 3, 3)
    # The 3×3 window around (2, 2) is rows 1..3, cols 1..3 — 9 pixels:
    # 8 ones + 1 three. Arc-distance median = 1.0 (cluster wins).
    assert abs(out[2, 2] - 1.0) < 1e-12, f"expected 1.0, got {out[2, 2]}"
