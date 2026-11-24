"""3D-Kragträger (Beton) mit eingebetteter Bewehrungslinie.

L x B x H, links eingespannt, Flächenlast auf der Stirnfläche x = L.
Bewehrung: Gerade y = B/2, z = COVER über die volle Länge, per gmsh fragment() eingebettet
und in DOLFINx als ridge_tags (Kodimension 2) verfügbar.
Einheiten: MN, m.
"""
import time

import gmsh
import numpy as np
import ufl
from dolfinx import fem, geometry
from dolfinx import mesh as dmesh
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import gmsh as dgmsh
from mpi4py import MPI

REBAR_TAG = 10

L, B, H, COVER = 2.0, 0.3, 0.5, 0.05
E_C, NU_C = 30e3, 0.2
EA_S = 200e3 * np.pi * 0.01**2       # Ø 20 mm

MU = E_C / (2 * (1 + NU_C))
LAM = E_C * NU_C / ((1 + NU_C) * (1 - 2 * NU_C))


def beam_with_rebar(hmax=0.08, comm=MPI.COMM_WORLD):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    occ = gmsh.model.occ
    vol = occ.addBox(0, 0, 0, L, B, H)
    line = occ.addLine(occ.addPoint(0, B / 2, COVER), occ.addPoint(L, B / 2, COVER))
    occ.fragment([(3, vol)], [(1, line)])
    occ.synchronize()
    vols = [t for _, t in gmsh.model.getEntities(3)]
    rebar = [t for _, t in gmsh.model.getEntities(1)
             if np.allclose(occ.getCenterOfMass(1, t), [L / 2, B / 2, COVER])
             and abs(occ.getMass(1, t) - L) < 1e-9]
    gmsh.model.addPhysicalGroup(3, vols, 1)
    gmsh.model.addPhysicalGroup(1, rebar, REBAR_TAG)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", hmax)
    gmsh.model.mesh.generate(3)
    md = dgmsh.model_to_mesh(gmsh.model, comm, 0, gdim=3)
    gmsh.finalize()
    return md.mesh, md.ridge_tags


def eps(w):
    return ufl.sym(ufl.grad(w))


def sigma(w):
    return 2 * MU * eps(w) + LAM * ufl.tr(eps(w)) * ufl.Identity(3)


def cantilever_bc_and_load(msh, V, q=1.0):
    fdim = msh.topology.dim - 1
    left = dmesh.locate_entities_boundary(msh, fdim, lambda x: np.isclose(x[0], 0.0))
    bc = fem.dirichletbc(np.zeros(3), fem.locate_dofs_topological(V, fdim, left), V)
    right = dmesh.locate_entities_boundary(msh, fdim, lambda x: np.isclose(x[0], L))
    ft = dmesh.meshtags(msh, fdim, right, np.full(len(right), 1, dtype=np.int32))
    v = ufl.TestFunction(V)
    load = ufl.dot(fem.Constant(msh, np.array([0.0, 0.0, -q])), v) * ufl.ds(subdomain_data=ft, subdomain_id=1)
    return bc, load


def rebar_form(msh, ridge_tags, u, v):
    """EA eps_t(u) eps_t(v) entlang der Bewehrung."""
    dr = ufl.Measure("dr", domain=msh, subdomain_data=ridge_tags)
    t = ufl.as_vector([1.0, 0.0, 0.0])
    eps_t = lambda w: ufl.dot(ufl.dot(ufl.grad(w), t), t)
    return EA_S * eps_t(u) * eps_t(v) * dr(REBAR_TAG)


def solve(a, rhs, bcs, prefix):
    t0 = time.perf_counter()
    prob = LinearProblem(a, rhs, bcs=bcs, petsc_options_prefix=prefix,
                         petsc_options={"ksp_type": "preonly", "pc_type": "lu",
                                        "pc_factor_mat_solver_type": "mumps"})
    uh = prob.solve()
    uh = uh[0] if isinstance(uh, tuple) else uh
    return uh, time.perf_counter() - t0


def eval_point(msh, f, point):
    pts = np.array([point], dtype=np.float64)
    tree = geometry.bb_tree(msh, msh.topology.dim)
    cells = geometry.compute_colliding_cells(msh, geometry.compute_collisions_points(tree, pts), pts)
    return f.eval(pts, cells.links(0)[:1])
