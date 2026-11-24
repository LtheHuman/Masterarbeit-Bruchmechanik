"""2D-Einfeldträger (ebener Verzerrungszustand, pro Meter Dicke) für die Bruchmechanik-Snippets.

Einheiten: MN, m (Spannung MPa = MJ/m^3, Energie pro Dicke MN = MJ/m).
Auflager als Streifen der Breite pad (links fest, rechts verschieblich), Gleichlast q oben.
Optional: kreisförmiger Einschluss (Zelltag INCL), sonst nur gesundes Material (Zelltag HEALTHY).
"""
from dataclasses import dataclass

import gmsh
import numpy as np
import ufl
from dolfinx import fem, geometry
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import gmsh as dgmsh
from mpi4py import MPI

HEALTHY, INCL = 1, 2
PAD_L, PAD_R, TOP = 1, 2, 3


@dataclass
class Beam2D:
    L: float = 2.0
    H: float = 0.5
    pad: float = 0.1
    E: float = 30e3
    nu: float = 0.2
    q: float = 0.3
    hmax: float = 0.02

    @property
    def mu(self):
        return self.E / (2 * (1 + self.nu))

    @property
    def lam(self):
        return self.E * self.nu / ((1 + self.nu) * (1 - 2 * self.nu))

    def eps(self, w):
        return ufl.sym(ufl.grad(w))

    def sigma(self, w, rho=1.0):
        e = self.eps(w)
        return rho * (2 * self.mu * e + self.lam * ufl.tr(e) * ufl.Identity(2))


def make_mesh(b: Beam2D, disc=None, refine_box=None, comm=MPI.COMM_WORLD):
    """Netz des Trägers.

    disc = (x0, y0, r):                  konform vernetzter Kreis (Zelltag INCL)
    refine_box = (x0, x1, y0, y1, h):     Elementgröße h im Rechteck (z.B. erwartete Nukleationszone)
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("beam2d")
    g = gmsh.model.geo
    coords = [(0, 0), (b.pad, 0), (b.L - b.pad, 0), (b.L, 0), (b.L, b.H), (0, b.H)]
    pts = [g.addPoint(x, y, 0) for x, y in coords]
    lines = [g.addLine(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
    outer = g.addCurveLoop(lines)
    arcs = []
    if disc is not None:
        x0, y0, r = disc
        c = g.addPoint(x0, y0, 0)
        rim = [g.addPoint(x0 + r * np.cos(t), y0 + r * np.sin(t), 0) for t in np.arange(4) * np.pi / 2]
        arcs = [g.addCircleArc(rim[i], c, rim[(i + 1) % 4]) for i in range(4)]
        loop = g.addCurveLoop(arcs)
        s_out = g.addPlaneSurface([outer, loop])
        s_in = g.addPlaneSurface([loop])
    else:
        s_out = g.addPlaneSurface([outer])
    g.synchronize()
    gmsh.model.addPhysicalGroup(2, [s_out], HEALTHY)
    if disc is not None:
        gmsh.model.addPhysicalGroup(2, [s_in], INCL)
    gmsh.model.addPhysicalGroup(1, [lines[0]], PAD_L)
    gmsh.model.addPhysicalGroup(1, [lines[2]], PAD_R)
    gmsh.model.addPhysicalGroup(1, [lines[4]], TOP)

    f = gmsh.model.mesh.field
    fields = []
    d_pad = f.add("Distance")
    f.setNumbers(d_pad, "PointsList", [pts[1], pts[2]])
    t_pad = f.add("Threshold")
    f.setNumber(t_pad, "InField", d_pad)
    for k, v in {"SizeMin": 2e-3, "SizeMax": b.hmax, "DistMin": 0.0, "DistMax": 0.1}.items():
        f.setNumber(t_pad, k, v)
    fields.append(t_pad)
    if disc is not None:
        r = disc[2]
        d_disc = f.add("Distance")
        f.setNumbers(d_disc, "CurvesList", arcs)
        f.setNumber(d_disc, "Sampling", 200)
        t_disc = f.add("Threshold")
        f.setNumber(t_disc, "InField", d_disc)
        for k, v in {"SizeMin": r / 12, "SizeMax": b.hmax, "DistMin": r / 2, "DistMax": 8 * r}.items():
            f.setNumber(t_disc, k, v)
        fields.append(t_disc)
    if refine_box is not None:
        x0, x1, y0, y1, h = refine_box
        box = f.add("Box")
        for k, v in {"VIn": h, "VOut": b.hmax, "XMin": x0, "XMax": x1, "YMin": y0, "YMax": y1,
                     "Thickness": 0.1}.items():
            f.setNumber(box, k, v)
        fields.append(box)
    fmin = f.add("Min")
    f.setNumbers(fmin, "FieldsList", fields)
    f.setAsBackgroundMesh(fmin)
    for opt in ["MeshSizeExtendFromBoundary", "MeshSizeFromPoints", "MeshSizeFromCurvature"]:
        gmsh.option.setNumber(f"Mesh.{opt}", 0)
    gmsh.model.mesh.generate(2)
    md = dgmsh.model_to_mesh(gmsh.model, comm, 0, gdim=2)
    gmsh.finalize()
    return md.mesh, md.cell_tags, md.facet_tags


def solve(b: Beam2D, msh, ft, rho=1.0, degree=2, prefix="beam_"):
    """Löst die Elastizität mit Steifigkeit rho*C. Rückgabe: uh, U (gespeichert), W (äußere Arbeit)."""
    V = fem.functionspace(msh, ("Lagrange", degree, (2,)))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    fdim = msh.topology.dim - 1
    msh.topology.create_connectivity(fdim, msh.topology.dim)
    bcs = [fem.dirichletbc(np.zeros(2), fem.locate_dofs_topological(V, fdim, ft.find(PAD_L)), V),
           fem.dirichletbc(fem.Constant(msh, 0.0),
                           fem.locate_dofs_topological(V.sub(1), fdim, ft.find(PAD_R)), V.sub(1))]
    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)
    t = fem.Constant(msh, np.array([0.0, -b.q]))
    prob = LinearProblem(ufl.inner(b.sigma(u, rho), b.eps(v)) * ufl.dx, ufl.dot(t, v) * ds(TOP),
                         bcs=bcs, petsc_options_prefix=prefix,
                         petsc_options={"ksp_type": "preonly", "pc_type": "lu",
                                        "pc_factor_mat_solver_type": "mumps"})
    uh = prob.solve()
    uh = uh[0] if isinstance(uh, tuple) else uh
    U = fem.assemble_scalar(fem.form(0.5 * ufl.inner(b.sigma(uh, rho), b.eps(uh)) * ufl.dx))
    W = fem.assemble_scalar(fem.form(ufl.dot(t, uh) * ds(TOP)))
    return uh, U, W


def rho_field(msh, ct, rho0):
    """DG0-Feld: 1 im gesunden Material, rho0 im Einschluss."""
    Q = fem.functionspace(msh, ("DG", 0))
    rho = fem.Function(Q)
    rho.x.array[:] = 1.0
    rho.x.array[ct.find(INCL)] = rho0
    return rho


def eval_at(msh, f, point):
    pts = np.array([[point[0], point[1], 0.0]])
    tree = geometry.bb_tree(msh, msh.topology.dim)
    cells = geometry.compute_colliding_cells(msh, geometry.compute_collisions_points(tree, pts), pts)
    return f.eval(pts, cells.links(0)[:1])


def cellwise(msh, expr):
    """Ausdruck je Zelle (DG0, im Schwerpunkt) als numpy-Array in Zellreihenfolge."""
    Q = fem.functionspace(msh, ("DG", 0))
    g = fem.Function(Q)
    g.interpolate(fem.Expression(expr, Q.element.interpolation_points))
    return g.x.array[Q.dofmap.list[:, 0]].copy()


def triangulation_data(msh, uh=None):
    """Knoten, Dreiecke (Zellreihenfolge) und optional Knotenverschiebungen für matplotlib."""
    V1 = fem.functionspace(msh, ("Lagrange", 1, (2,)))
    xy = V1.tabulate_dof_coordinates()[:, :2]
    disp = None
    if uh is not None:
        u1 = fem.Function(V1)
        u1.interpolate(uh)
        disp = u1.x.array.reshape(-1, 2)
    return xy, V1.dofmap.list, disp
