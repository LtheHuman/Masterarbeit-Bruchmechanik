"""Bausteine des Algorithmus von Xavier et al. in 3D, ein Schritt:

  1. Elastizität lösen (Schädigung als weiche Phase rho0*C, DG0-Feld)
  2. Topologische Ableitung D_T J = -P_gamma sigma : eps auswerten. Hier vereinfacht nur
     als Energiedichte-Indikator, der 3D-Polarisationstensor fehlt noch.
  3. Zellen mit stärkstem Indikator markieren, Netz lokal verfeinern (Plaza)
"""
import numpy as np
import ufl
from dolfinx import fem
from dolfinx import mesh as dmesh

from _common import beam_with_rebar, cantilever_bc_and_load, eps, sigma, solve

msh, rt = beam_with_rebar(hmax=0.1)
Q = fem.functionspace(msh, ("DG", 0))
rho = fem.Function(Q)
rho.x.array[:] = 1.0                                   # 1 = intakt, rho0 = geschädigt

V = fem.functionspace(msh, ("Lagrange", 1, (3,)))      # Xavier et al.: lineare Elemente
u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
bc, rhs = cantilever_bc_and_load(msh, V)
uh, _ = solve(rho * ufl.inner(sigma(u), eps(v)) * ufl.dx, rhs, [bc], prefix="td_")

# Indikator elementweise (Projektion auf DG0 = Interpolation des Ausdrucks im Schwerpunkt)
ind = fem.Function(Q)
ind.interpolate(fem.Expression(0.5 * ufl.inner(sigma(uh), eps(uh)), Q.element.interpolation_points))
thr = np.quantile(ind.x.array, 0.99)
cells = np.flatnonzero(ind.x.array >= thr).astype(np.int32)
print(f"markiert: {len(cells)} Zellen (oberstes 1 % der Energiedichte)")

msh.topology.create_connectivity(3, 1)
edges = dmesh.compute_incident_entities(msh.topology, cells, 3, 1)
res = dmesh.refine(msh, edges)
fine = res[0] if isinstance(res, tuple) else res
print("Zellen:", msh.topology.index_map(3).size_global, "->", fine.topology.index_map(3).size_global)
print("Offen: ridge_tags (Bewehrung) auf das verfeinerte Netz übertragen. refine() liefert",
      "die Parent-Maps, die Übertragung muss man selbst schreiben.")
