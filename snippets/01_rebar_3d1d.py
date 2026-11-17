"""3D-Beton + 1D-Bewehrung über Ridge-Integral (DOLFINx >= 0.10).

Prüft:
  1. Länge der Bewehrung über dr        -> exakt L
  2. Stabdehnung eines linearen Feldes  -> exakt 0.01
  3. Durchbiegung mit/ohne Bewehrung    -> Bewehrung versteift, Durchbiegung sinkt
Klassische Elastizität mit 1D-Einschluss ist nicht netzkonvergent (siehe Snippet 02).
"""
import numpy as np
import ufl
from dolfinx import fem

from _common import (L, B, H, beam_with_rebar, cantilever_bc_and_load, eval_point,
                     eps, rebar_form, sigma, solve, REBAR_TAG)

msh, rt = beam_with_rebar(hmax=0.08)
dr = ufl.Measure("dr", domain=msh, subdomain_data=rt)
print("Bewehrungslänge via dr:", fem.assemble_scalar(fem.form(1.0 * dr(REBAR_TAG))), "(soll", L, ")")

V = fem.functionspace(msh, ("Lagrange", 2, (3,)))
g = fem.Function(V)
g.interpolate(lambda x: np.vstack([0.01 * x[0] + 0.3 * x[2], 0 * x[0], 0 * x[0]]))
t = ufl.as_vector([1.0, 0.0, 0.0])
eps_t = fem.assemble_scalar(fem.form(ufl.dot(ufl.dot(ufl.grad(g), t), t) * dr(REBAR_TAG))) / L
print(f"mittlere Stabdehnung: {eps_t:.6f} (soll 0.010000)")

u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
bc, rhs = cantilever_bc_and_load(msh, V)
a_vol = ufl.inner(sigma(u), eps(v)) * ufl.dx
a_bar = rebar_form(msh, rt, u, v)
for name, a in [("ohne Bewehrung", a_vol), ("mit Bewehrung", a_vol + a_bar)]:
    uh, dt = solve(a, rhs, [bc], prefix=name.replace(" ", "_") + "_")
    w = eval_point(msh, uh, [L, B / 2, H / 2])[2]
    print(f"{name:15s} w(L) = {w:.5e}   ({dt:.1f} s, {V.dofmap.index_map.size_global * 3} DOF)")
