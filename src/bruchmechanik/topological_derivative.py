"""Topologische Ableitung der Energie bei Nukleation eines kleinen kreisförmigen Einschlusses.

Xavier, Fancello, Farias, Van Goethem, Novotny, Gl. (3.3)-(3.5), 2D (ebener Verzerrungszustand,
Lamé-Konstanten lam, mu des 2D-Modells):

    D_T Pi(x) = -P_gamma sigma(u(x)) : eps(u(x))
    P_gamma   = 1/2 (1-gamma)/(1+gamma a2) [ (1+a2) II + 1/2 (a1-a2) (1-gamma)/(1+gamma a1) I(x)I ]
    a1 = (lam+mu)/mu,   a2 = (lam+3mu)/(lam+mu)

Normierung: Pi(eps) - Pi(0) = D_T Pi(x) * |B_eps| + o(eps^2), mit |B_eps| = pi eps^2.
gamma = rho0 im gesunden Material (Schädigung einfügen), gamma = 1/rho0 im geschädigten (Heilung).
Gilt für die Gesamtpotentialenergie (Lastvorgabe) und die gespeicherte Energie (Verschiebungsvorgabe).
"""
import ufl


def polarization_coefficients(gamma, lam, mu):
    """P_gamma = alpha * II + beta * I(x)I  ->  (alpha, beta)."""
    a1 = (lam + mu) / mu
    a2 = (lam + 3 * mu) / (lam + mu)
    c = 0.5 * (1 - gamma) / (1 + gamma * a2)
    alpha = c * (1 + a2)
    beta = c * 0.5 * (a1 - a2) * (1 - gamma) / (1 + gamma * a1)
    return alpha, beta


def dt_energy(sig, eps, gamma, lam, mu):
    """UFL-Ausdruck für D_T Pi = -(alpha sigma:eps + beta tr(sigma) tr(eps))."""
    alpha, beta = polarization_coefficients(gamma, lam, mu)
    return -(alpha * ufl.inner(sig, eps) + beta * ufl.tr(sig) * ufl.tr(eps))
