"""Griffith-Energie aus der Besprechung vom 17.11.2026, live visualisiert.

    E(a) = 1/2 int sigma(a):eps(a) dx + G_c a - W,      W = int f . u ds

Einfeldträger (2D, ebener Verzerrungszustand, pro Meter Dicke) unter Gleichlast oben,
Randriss der Länge a unten in Feldmitte. Für eine Folge von Risslängen wird jeweils
neu vernetzt und gelöst. Live angezeigt werden:
  oben links   Energiedichte auf dem verformten Träger (überhöht)
  oben rechts  Ausschnitt um die Rissspitze (Rissöffnung, Singularität)
  unten links  Energiebilanz relativ zu a = 0:  dPi, G_c a, dE = E(a) - E(0)
  unten rechts log-log: -dPi gegen G_c a. Wo -dPi die Gerade G_c a übersteigt, gilt E(a) < E(0).

Prüfungen:
  * Clapeyron: im Gleichgewicht ist W = 2 U  ->  Pi = U - W = -U
  * Skalierung: -dPi ~ a^2 für kleine a (die Normierung a^2 der Besprechung)
  * LEFM: -dPi(a) = int_0^a K(s)^2/E' ds mit K = sigma sqrt(pi a) F(a/H)
          (F für Randriss unter reiner Biegung, Tada/Brown-Srawley, gilt für a/H <= 0.6)

Der Riss ist eine Kerbe der Breite delta << a (kein Nullbreiten-Schlitz).

Aufruf:  python 04_griffith_energy.py            (Live-Fenster)
         python 04_griffith_energy.py --no-show  (nur Dateien schreiben)
Ausgabe: results/raw/04_griffith_energy/{runs.jsonl, final.png, verlauf.gif}
"""
import argparse
import json
import time
from pathlib import Path

import gmsh
import matplotlib
import numpy as np
import ufl
from dolfinx import fem, geometry
from dolfinx import mesh as dmesh
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import gmsh as dgmsh
from mpi4py import MPI

ap = argparse.ArgumentParser()
ap.add_argument("--no-show", action="store_true", help="kein Live-Fenster")
ap.add_argument("--q", type=float, default=0.3, help="Gleichlast oben [MN/m]")
ap.add_argument("--gc", type=float, default=1e-4, help="Bruchenergie G_c [MN/m] (1e-4 = 100 N/m)")
ap.add_argument("--n", type=int, default=14, help="Anzahl Risslängen")
args = ap.parse_args()
if args.no_show:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.tri import Triangulation  # noqa: E402

# ---- Parameter [MN, m] -------------------------------------------------------
L, H = 2.0, 0.5            # Träger
PAD = 0.1                  # Auflagerbreite links/rechts (vermeidet Punktlager-Singularität)
DELTA = 5e-4               # Kerbbreite
HMIN, HMAX = DELTA / 2, 0.025
E_MOD, NU = 30e3, 0.2      # Beton
E_PRIME = E_MOD / (1 - NU**2)
MU = E_MOD / (2 * (1 + NU))
LAM = E_MOD * NU / ((1 + NU) * (1 - 2 * NU))
Q, GC = args.q, args.gc
A_LIST = np.concatenate([[0.0], np.geomspace(0.01, 0.25, args.n)])
TAG_PAD_L, TAG_PAD_R, TAG_TOP, TAG_CRACK = 1, 2, 3, 4

OUT = Path(__file__).resolve().parents[1] / "results" / "raw" / "04_griffith_energy"
OUT.mkdir(parents=True, exist_ok=True)


# ---- Netz -------------------------------------------------------------------
def make_mesh(a):
    """Träger mit Randkerbe (Breite DELTA, Länge a) unten in Feldmitte; a = 0: intakt."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(f"beam_a{a:.4f}")
    g = gmsh.model.geo
    xm = L / 2
    if a > 0:
        bottom = [(xm - DELTA / 2, 0), (xm - DELTA / 2, a), (xm + DELTA / 2, a), (xm + DELTA / 2, 0)]
    else:
        bottom = [(xm, 0)]
    coords = [(0, 0), (PAD, 0), *bottom, (L - PAD, 0), (L, 0), (L, H), (0, H)]
    pts = [g.addPoint(x, y, 0) for x, y in coords]
    lines = [g.addLine(pts[i], pts[(i + 1) % len(pts)]) for i in range(len(pts))]
    g.addPlaneSurface([g.addCurveLoop(lines)])
    g.synchronize()

    nb = len(bottom)
    crack = lines[2:2 + nb - 1] if a > 0 else []
    gmsh.model.addPhysicalGroup(2, [1], 1)
    gmsh.model.addPhysicalGroup(1, [lines[0]], TAG_PAD_L)
    gmsh.model.addPhysicalGroup(1, [lines[nb + 2]], TAG_PAD_R)
    gmsh.model.addPhysicalGroup(1, [lines[nb + 4]], TAG_TOP)
    if crack:
        gmsh.model.addPhysicalGroup(1, crack, TAG_CRACK)

    # Netzfeinheit: fein an Kerbe/Rissspitze, moderat an den Auflagerrändern
    f = gmsh.model.mesh.field
    d_crack = f.add("Distance")
    if crack:
        f.setNumbers(d_crack, "CurvesList", crack)
        f.setNumber(d_crack, "Sampling", 400)
    else:
        f.setNumbers(d_crack, "PointsList", [pts[2]])
    t_crack = f.add("Threshold")
    f.setNumber(t_crack, "InField", d_crack)
    for k, v in {"SizeMin": HMIN, "SizeMax": HMAX, "DistMin": DELTA, "DistMax": 0.15}.items():
        f.setNumber(t_crack, k, v)
    d_pad = f.add("Distance")
    f.setNumbers(d_pad, "PointsList", [pts[1], pts[nb + 2]])
    t_pad = f.add("Threshold")
    f.setNumber(t_pad, "InField", d_pad)
    for k, v in {"SizeMin": 2e-3, "SizeMax": HMAX, "DistMin": 0.0, "DistMax": 0.1}.items():
        f.setNumber(t_pad, k, v)
    fmin = f.add("Min")
    f.setNumbers(fmin, "FieldsList", [t_crack, t_pad])
    f.setAsBackgroundMesh(fmin)
    for opt in ["MeshSizeExtendFromBoundary", "MeshSizeFromPoints", "MeshSizeFromCurvature"]:
        gmsh.option.setNumber(f"Mesh.{opt}", 0)
    gmsh.model.mesh.generate(2)
    md = dgmsh.model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, gdim=2)
    gmsh.finalize()
    return md.mesh, md.facet_tags


# ---- Mechanik ---------------------------------------------------------------
def eps(w):
    return ufl.sym(ufl.grad(w))


def sigma(w):
    return 2 * MU * eps(w) + LAM * ufl.tr(eps(w)) * ufl.Identity(2)


def solve_beam(a):
    msh, ft = make_mesh(a)
    V = fem.functionspace(msh, ("Lagrange", 2, (2,)))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    fdim = msh.topology.dim - 1
    msh.topology.create_connectivity(fdim, msh.topology.dim)
    dofs_l = fem.locate_dofs_topological(V, fdim, ft.find(TAG_PAD_L))
    dofs_r = fem.locate_dofs_topological(V.sub(1), fdim, ft.find(TAG_PAD_R))
    bcs = [fem.dirichletbc(np.zeros(2), dofs_l, V),                       # Festlager
           fem.dirichletbc(fem.Constant(msh, 0.0), dofs_r, V.sub(1))]     # Loslager
    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)
    traction = fem.Constant(msh, np.array([0.0, -Q]))
    prob = LinearProblem(ufl.inner(sigma(u), eps(v)) * ufl.dx, ufl.dot(traction, v) * ds(TAG_TOP),
                         bcs=bcs, petsc_options_prefix=f"griffith_{a:.4f}_",
                         petsc_options={"ksp_type": "preonly", "pc_type": "lu",
                                        "pc_factor_mat_solver_type": "mumps"})
    uh = prob.solve()
    uh = uh[0] if isinstance(uh, tuple) else uh
    U = fem.assemble_scalar(fem.form(0.5 * ufl.inner(sigma(uh), eps(uh)) * ufl.dx))
    W = fem.assemble_scalar(fem.form(ufl.dot(traction, uh) * ds(TAG_TOP)))
    return msh, uh, U, W, V.dofmap.index_map.size_global * 2


def sigma_xx_at(msh, uh, point):
    S = fem.functionspace(msh, ("DG", 1))
    s = fem.Function(S)
    s.interpolate(fem.Expression(sigma(uh)[0, 0], S.element.interpolation_points))
    pts = np.array([[*point, 0.0]])
    tree = geometry.bb_tree(msh, 2)
    cells = geometry.compute_colliding_cells(msh, geometry.compute_collisions_points(tree, pts), pts)
    return float(s.eval(pts, cells.links(0)[:1])[0])


def plot_data(msh, uh):
    """Triangulierung (P1-Knoten), Verschiebung an den Knoten, Energiedichte je Zelle."""
    V1 = fem.functionspace(msh, ("Lagrange", 1, (2,)))
    u1 = fem.Function(V1)
    u1.interpolate(uh)
    xy = V1.tabulate_dof_coordinates()[:, :2]
    tris = V1.dofmap.list
    disp = u1.x.array.reshape(-1, 2)
    Q0 = fem.functionspace(msh, ("DG", 0))
    psi = fem.Function(Q0)
    psi.interpolate(fem.Expression(0.5 * ufl.inner(sigma(uh), eps(uh)), Q0.element.interpolation_points))
    return xy, tris, disp, psi.x.array[Q0.dofmap.list[:, 0]]


# ---- LEFM-Referenz ------------------------------------------------------------
def F_bend(alpha):
    return 1.122 - 1.40 * alpha + 7.33 * alpha**2 - 13.08 * alpha**3 + 14.0 * alpha**4


def dPi_lefm(a, sig):
    s = np.linspace(0, a, 400)
    G = (sig * np.sqrt(np.pi * s) * F_bend(s / H))**2 / E_PRIME
    return -np.trapezoid(G, s)


# ---- Abbildung ----------------------------------------------------------------
plt.ion()
fig = plt.figure(figsize=(14, 8.5), layout="constrained")
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.1], width_ratios=[1.7, 1])
ax_beam, ax_tip = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
ax_en, ax_log = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
frames, rows = [], []
norm = scale = cbar = None


def draw_fields(a, xy, tris, disp, psi):
    global norm, scale, cbar
    if norm is None:  # aus dem intakten Zustand, dann fest für alle Bilder
        norm = LogNorm(vmin=psi.max() * 1e-3, vmax=psi.max() * 30)
        scale = 0.04 * H / np.abs(disp).max()
    tri = Triangulation(*(xy + scale * disp).T, tris)
    for ax in (ax_beam, ax_tip):
        ax.clear()
        pc = ax.tripcolor(tri, facecolors=np.clip(psi, norm.vmin, norm.vmax), norm=norm,
                          cmap="inferno", shading="flat", rasterized=True)
        ax.set_aspect("equal")
    ax_beam.set_title(f"Energiedichte ½σ:ε [MJ/m³], Verformung {scale:.0f}-fach, a = {a * 100:.1f} cm")
    ax_beam.set_xlim(-0.05, L + 0.05)
    ax_beam.set_ylim(-0.1, H + 0.05)
    for x0 in (0, L - PAD):
        ax_beam.add_patch(plt.Rectangle((x0, -0.04), PAD, 0.03, color="0.4"))
    ax_beam.annotate("", xy=(L / 2, H + 0.005), xytext=(L / 2, H + 0.045),
                     arrowprops=dict(arrowstyle="->", color="tab:blue"))
    ax_beam.text(L / 2 + 0.02, H + 0.02, f"q = {Q} MN/m", color="tab:blue", fontsize=9)
    w = max(0.06, 1.3 * a)
    ax_tip.set_xlim(L / 2 - w, L / 2 + w)
    ax_tip.set_ylim(-0.01, 1.6 * w)
    ax_tip.set_title("Ausschnitt Rissspitze")
    if cbar is None:
        cbar = fig.colorbar(pc, ax=[ax_beam, ax_tip], shrink=0.8, location="right")


def draw_energies(rows, sig0, final=False):
    a = np.array([r["a"] for r in rows])
    dPi = np.array([r["dPi"] for r in rows])
    gca_ = GC * a
    dE = dPi + gca_
    ax_en.clear()
    ax_en.axhline(0, color="0.5", lw=0.8)
    ax_en.plot(a * 100, gca_ * 1e6, "-o", color="tab:green", label="$G_c\\,a$  (Dissipation)")
    ax_en.plot(a * 100, dPi * 1e6, "-o", color="tab:blue", label="$\\Delta\\Pi = \\Pi(a)-\\Pi(0)$  (elastisch)")
    ax_en.plot(a * 100, dE * 1e6, "-o", color="tab:red", lw=2.2, label="$E(a)-E(0)$")
    if final:
        aa = np.linspace(1e-4, a.max(), 300)
        ax_en.plot(aa * 100, [(dPi_lefm(x, sig0) + GC * x) * 1e6 for x in aa], ":", color="tab:red",
                   label="$E(a)-E(0)$ LEFM")
        i = np.argmax(dE)
        ax_en.annotate(f"Energiebarriere\n$a^* \\approx$ {a[i] * 100:.1f} cm", (a[i] * 100, dE[i] * 1e6),
                       xytext=(75, 22), textcoords="offset points", ha="center", fontsize=9,
                       arrowprops=dict(arrowstyle="->"))
        k = np.flatnonzero((dE[1:] < 0) & (dE[:-1] >= 0) & (a[1:] > 0))
        if k.size:
            j = k[0]
            ac = a[j] + (a[j + 1] - a[j]) * dE[j] / (dE[j] - dE[j + 1])
            ax_en.axvspan(ac * 100, a.max() * 100, color="tab:red", alpha=0.08)
            ax_en.axvline(ac * 100, color="tab:red", lw=0.8, ls="--")
            ax_en.text(ac * 100, ax_en.get_ylim()[1] * 0.9, f" $a_c \\approx$ {ac * 100:.1f} cm\n E(a) < E(0) → Riss",
                       color="tab:red", va="top", fontsize=9)
    ax_en.set_xlabel("Risslänge a [cm]")
    ax_en.set_ylabel("Energie [J/m]  (pro Meter Dicke)")
    ax_en.set_title("Energiebilanz aus der Besprechung")
    ax_en.legend(loc="lower left", fontsize=9)

    ax_log.clear()
    m = a > 0
    if m.any():
        ax_log.loglog(a[m] * 100, -dPi[m] * 1e6, "o", color="tab:blue", label="$-\\Delta\\Pi$ FEM")
        aa = np.geomspace(a[m].min(), A_LIST.max(), 100)
        ax_log.loglog(aa * 100, [-dPi_lefm(x, sig0) * 1e6 for x in aa], "-", color="tab:blue", alpha=0.5,
                      label="$-\\Delta\\Pi$ LEFM")
        ax_log.loglog(aa * 100, GC * aa * 1e6, "-", color="tab:green", label="$G_c\\,a$  (Steigung 1)")
        ref = -dPi[m][0] * (aa / a[m][0])**2
        ax_log.loglog(aa * 100, ref * 1e6, "--", color="0.5", label="Steigung 2")
        ax_log.legend(fontsize=9)
    ax_log.set_xlabel("a [cm]")
    ax_log.set_title("Skalierung: E(a) < E(0) wo blau über grün")


def snapshot():
    fig.canvas.draw()
    frames.append(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy())
    if not args.no_show:
        plt.pause(0.05)


# ---- Hauptschleife ------------------------------------------------------------
print(f"{'a [cm]':>7} {'DOF':>7} {'W/2U-1':>9} {'dPi [J/m]':>11} {'dPi/a^2':>9} "
      f"{'G_c a [J/m]':>11} {'dE [J/m]':>10} {'dE/a^2':>9} {'FEM/LEFM':>9}  t[s]")
Pi0 = sig0 = None
for a in A_LIST:
    t0 = time.perf_counter()
    msh, uh, U, W, ndof = solve_beam(a)
    Pi = U - W
    if a == 0:
        Pi0 = Pi
        sig0 = sigma_xx_at(msh, uh, (L / 2, 1e-6))
    dPi = Pi - Pi0
    row = dict(a=a, ndof=ndof, U=U, W=W, Pi=Pi, dPi=dPi, GCa=GC * a, dE=dPi + GC * a,
               clapeyron=W / (2 * U) - 1, sigma0=sig0, q=Q, Gc=GC, delta=DELTA, hmin=HMIN,
               t=time.perf_counter() - t0)
    rows.append(row)
    r2 = f"{dPi / a**2 * 1e6:9.2f} " if a > 0 else f"{'-':>9} "
    e2 = f"{row['dE'] / a**2 * 1e6:9.2f}" if a > 0 else f"{'-':>9}"
    row["ratio_lefm"] = dPi / dPi_lefm(a, sig0) if a > 0 else None
    rl = f"{row['ratio_lefm']:9.3f}" if a > 0 else f"{'-':>9}"
    print(f"{a * 100:7.2f} {ndof:7d} {row['clapeyron']:9.1e} {dPi * 1e6:11.4f} {r2}"
          f"{GC * a * 1e6:11.4f} {row['dE'] * 1e6:10.4f} {e2} {rl}  {row['t']:.1f}")
    draw_fields(a, *plot_data(msh, uh))
    draw_energies(rows, sig0)
    snapshot()

# ---- Auswertung ---------------------------------------------------------------
a = np.array([r["a"] for r in rows])[1:]
dPi = np.array([r["dPi"] for r in rows])[1:]
small = a <= 0.05
p = np.polyfit(np.log(a[small]), np.log(-dPi[small]), 1)[0]
C_lefm = np.pi * (1.122 * sig0)**2 / (2 * E_PRIME)
print(f"\nsigma_xx(L/2, 0) intakt = {sig0:.3f} MPa")
print(f"Exponent -dPi ~ a^p (a <= 5 cm):   p = {p:.3f}   (Grenzwert 2; knapp darunter, weil F(a/H) mit a fällt)")
print(f"-dPi/a^2:  FEM bei a = {a[0] * 100:.0f} cm: {-dPi[0] / a[0]**2 * 1e6:.1f} J/m^3,  "
      f"LEFM-Grenzwert a -> 0: pi(1.122 sigma)^2/(2E') = {C_lefm * 1e6:.1f} J/m^3")
print(f"LEFM:  a* = G_c/(2C) = {GC / (2 * C_lefm) * 100:.1f} cm (Barriere),  "
      f"a_c = G_c/C = {GC / C_lefm * 100:.1f} cm (E(a) = E(0)),  jeweils für kleine a/H")
print("(E(a)-E(0))/a^2 = G_c/a - C  ->  +inf für a -> 0: Mit Dissipation ~ a gibt es keine Nukleation aus dem intakten Zustand.")

draw_energies(rows, sig0, final=True)
fig.suptitle(f"Griffith: E(a) = ½∫σ:ε dx + G_c a − W,   G_c = {GC * 1e6:.0f} N/m,   "
             f"Skalierung −ΔΠ ~ a^{p:.2f}", fontsize=13)
snapshot()
fig.savefig(OUT / "final.png", dpi=130)
with open(OUT / "runs.jsonl", "w") as fh:
    for r in rows:
        fh.write(json.dumps(r) + "\n")
from PIL import Image  # noqa: E402

imgs = [Image.fromarray(f) for f in frames]
imgs[0].save(OUT / "verlauf.gif", save_all=True, append_images=imgs[1:] + [imgs[-1]] * 4, duration=600, loop=0)
print(f"\nGeschrieben: {OUT}/final.png, verlauf.gif, runs.jsonl")
if not args.no_show:
    plt.ioff()
    plt.show()
