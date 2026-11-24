"""Vereinfachte Strain-Gradient-Elastizität mit Bewehrung, C0-Interior-Penalty, P2.

    1/2 int sigma:eps dx + mu l^2 int grad grad u : grad grad u dx + Stabenergie

C0-IP nach Engel et al. (2002). Erster Lauf mit l > 0: ca. 45 s JIT-Kompilierung.
"""
import ufl
from dolfinx import fem

from bruchmechanik.beam3d import (L, B, H, MU, beam_with_rebar, cantilever_bc_and_load, eps,
                     eval_point, rebar_form, sigma, solve)

msh, rt = beam_with_rebar(hmax=0.08)
V = fem.functionspace(msh, ("Lagrange", 2, (3,)))
u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
bc, rhs = cantilever_bc_and_load(msh, V)

n = ufl.FacetNormal(msh)("+")
h = ufl.avg(ufl.CellDiameter(msh))
Hess = lambda w: ufl.grad(ufl.grad(w))
jump_grad = lambda w: ufl.grad(w)("+") - ufl.grad(w)("-")


def a_strain_gradient(ell, alpha=20.0):
    c = fem.Constant(msh, MU * ell**2)
    return c * (ufl.inner(Hess(u), Hess(v)) * ufl.dx
                - ufl.inner(ufl.dot(ufl.avg(Hess(u)), n), jump_grad(v)) * ufl.dS
                - ufl.inner(ufl.dot(ufl.avg(Hess(v)), n), jump_grad(u)) * ufl.dS
                + alpha / h * ufl.inner(jump_grad(u), jump_grad(v)) * ufl.dS)


a_cl = ufl.inner(sigma(u), eps(v)) * ufl.dx + rebar_form(msh, rt, u, v)
for ell in [0.0, 0.02, 0.05]:
    a = a_cl + a_strain_gradient(ell) if ell > 0 else a_cl
    uh, dt = solve(a, rhs, [bc], prefix=f"sg{int(ell * 1000)}_")
    print(f"l = {ell:4.2f} m   w(L) = {eval_point(msh, uh, [L, B / 2, H / 2])[2]:.5e}   ({dt:.1f} s)")
