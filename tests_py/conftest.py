# Shared fixtures for the pycirclemedianfilter Python test suite.

import numpy as np
import pytest


@pytest.fixture
def small_phase_image():
    """Tiny 5×5 phase image with one obvious cluster — values in [-π, π]."""
    rng = np.random.default_rng(0)
    return rng.uniform(-np.pi, np.pi, size=(5, 5))


@pytest.fixture
def constant_phase_image():
    """8×8 image where every pixel equals 0.7 rad. The median filter must
    leave it unchanged regardless of filter radius."""
    return np.full((8, 8), 0.7)
