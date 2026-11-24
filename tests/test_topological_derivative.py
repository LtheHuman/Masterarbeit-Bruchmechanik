"""Grenzfälle des Polarisationstensors. FEM-Prüfung: snippets/05_topological_derivative.py."""
import numpy as np
import pytest

from bruchmechanik.topological_derivative import polarization_coefficients

E, NU = 30e3, 0.2
MU = E / (2 * (1 + NU))
LAM = E * NU / ((1 + NU) * (1 - 2 * NU))


def test_no_contrast_gives_zero():
    assert polarization_coefficients(1.0, LAM, MU) == (0.0, 0.0)


@pytest.mark.parametrize("eta", [1e-3, 1e-5])
def test_weak_contrast_limit(eta):
    # Störungsrechnung 1. Ordnung: Pi(eps)-Pi(0) = 1/2 (gamma-1) |B| sigma:eps  ->  P = 1/2 (1-gamma) II
    alpha, beta = polarization_coefficients(1 - eta, LAM, MU)
    assert np.isclose(alpha, 0.5 * eta, rtol=10 * eta)
    assert abs(beta) < 10 * eta**2


def test_softening_lowers_energy():
    # weicher Einschluss (gamma < 1): P positiv definit  ->  D_T Pi < 0
    alpha, beta = polarization_coefficients(1e-6, LAM, MU)
    assert alpha > 0 and alpha + 2 * beta > 0  # Eigenwerte auf Deviator- bzw. Kugelanteil (2D)
