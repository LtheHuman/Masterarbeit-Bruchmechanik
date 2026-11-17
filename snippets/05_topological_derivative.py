"""Riss als Volumen: topologische Ableitung nach Xavier et al., geprüft und angewendet.

Teil 1, Prüfung der Formel (3.3)-(3.4) durch finite Störung:
  In einen gesunden Punkt x0 wird ein Kreis B_eps mit Steifigkeit rho0*C gesetzt.
      (Pi(eps) - Pi(0)) / (pi eps^2)  ->  D_T Pi(x0) = -P_gamma sigma:eps     für eps -> 0
  Pi(0) wird auf DEMSELBEN Netz gerechnet (Kreis vernetzt, aber rho = 1), damit sich der
  Diskretisierungsfehler heraushebt und nur die Störung bleibt.
  Drei Punkte: A Zugzone (Feldmitte unten), B Schubzone (halbe Höhe), C Druckzone (oben).

Teil 2, ein Schritt von Algorithmus 1 (Francfort-Marigo):
  F = Pi + kappa_delta |omega|,  kappa_delta = kappa_s / delta,  kappa_s = G_c.
  (Ein Schadensband der Breite delta und Länge a dissipiert kappa_delta*delta*a = G_c*a, also Griffith.)
  D_T F = D_T Pi + kappa_delta im gesunden Material. Nukleation, wo D_T F < 0.
  Da D_T Pi ~ q^2, folgt die Nukleationslast q_c = q sqrt(kappa_delta / max(-D_T Pi)).
  omega*:  {D_T F < 0},   omega^beta: die Zellen mit kleinstem D_T F bis |omega^beta| = pi delta^2/4.
  Wie im Paper werden die Auflagerränder (Spannungssingularität) ausgeblendet.

Aufruf:  python 05_topological_derivative.py [--no-show] [--delta 0.01]
Ausgabe: results/raw/05_topological_derivative/{final.png, verlauf.gif, runs.jsonl}
"""
import argparse
import json
import sys
import time
from pathlib import Path

import matplotlib
import numpy as np
import ufl
from dolfinx import fem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from bruchmechanik.beam2d import (Beam2D, INCL, cellwise, eval_at, make_mesh, rho_field,  # noqa: E402
                                  solve, triangulation_data)
from bruchmechanik.topological_derivative import dt_energy, polarization_coefficients  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--no-show", action="store_true")
ap.add_argument("--delta", type=float, default=0.01, help="Schadensbreite delta [m] (Paper: 0.01)")
ap.add_argument("--gc", type=float, default=1e-4, help="kappa_s = G_c [MN/m]")
args = ap.parse_args()
if args.no_show:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.tri import Triangulation  # noqa: E402

beam = Beam2D()
RHO0 = 1e-6                                   # Paper, Tabelle 1
KAPPA = args.gc / args.delta                  # kappa_delta [MN/m^2 = MJ/m^3]
POINTS = {"A Zug": (beam.L / 2, 0.04), "B Schub": (0.5, beam.H / 2), "C Druck": (beam.L / 2, beam.H - 0.04)}
COLORS = {"A Zug": "tab:red", "B Schub": "tab:purple", "C Druck": "tab:blue"}
EPS_LIST = np.geomspace(0.03, 0.0025, 7)
OUT = ROOT / "results" / "raw" / "05_topological_derivative"
OUT.mkdir(parents=True, exist_ok=True)
J = 1e6                                        # MJ -> J für die Anzeige

plt.ion()
fig = plt.figure(figsize=(14, 8.5), layout="constrained")
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.1], width_ratios=[1.7, 1])
ax_f, ax_z = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
ax_c, ax_e = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
frames = []


def snapshot():
    fig.canvas.draw()
    frames.append(np.asarray(fig.canvas.buffer_rgba())[..., :3].copy())
    if not args.no_show:
        plt.pause(0.05)


# ---- Intakter Träger: Feld der topologischen Ableitung -------------------------
# Zugzone fein (h = delta/5), damit omega^beta mit |omega^beta| = pi delta^2/4 auflösbar ist
msh0, ct0, ft0 = make_mesh(beam, refine_box=(beam.L / 2 - 0.25, beam.L / 2 + 0.25, 0.0, 0.06, args.delta / 5))
uh0, U0, W0 = solve(beam, msh0, ft0, prefix="intakt_")
dt_expr = dt_energy(beam.sigma(uh0), beam.eps(uh0), RHO0, beam.lam, beam.mu)
dtPi = cellwise(msh0, dt_expr)                              # [MJ/m^3], je Zelle
trs = cellwise(msh0, ufl.tr(beam.sigma(uh0)))
area = cellwise(msh0, ufl.CellVolume(msh0))
mid = cellwise(msh0, ufl.SpatialCoordinate(msh0)[0]), cellwise(msh0, ufl.SpatialCoordinate(msh0)[1])
near_support = (mid[1] < 0.15) & ((np.abs(mid[0] - beam.pad) < 0.15) | (np.abs(mid[0] - (beam.L - beam.pad)) < 0.15)
                                  | (mid[0] < beam.pad) | (mid[0] > beam.L - beam.pad))
xy0, tris0, _ = triangulation_data(msh0)
tri0 = Triangulation(*xy0.T, tris0)


def draw_field(highlight=None):
    ax_f.clear()
    val = np.clip(-dtPi * J, 1e-1, None)
    pc = ax_f.tripcolor(tri0, facecolors=val, norm=LogNorm(vmin=1, vmax=val[~near_support].max() * 3),
                        cmap="viridis", shading="flat", rasterized=True)
    ax_f.tripcolor(tri0, facecolors=np.where(near_support, 1.0, np.nan), cmap="Greys", vmin=0, vmax=2,
                   alpha=0.55, shading="flat")
    for name, (x, y) in POINTS.items():
        ax_f.plot(x, y, "o", ms=8, mfc="none", mec=COLORS[name], mew=2)
        ax_f.text(x + 0.03, y, name, color=COLORS[name], fontsize=9, va="center",
                  fontweight="bold" if name == highlight else "normal")
    ax_f.set_aspect("equal")
    ax_f.set_xlim(-0.03, beam.L + 0.03)
    ax_f.set_ylim(-0.03, beam.H + 0.03)
    ax_f.set_title(f"$-D_T\\Pi$ [J/m³] bei q = {beam.q} MN/m  (grau: Auflagerzonen ausgeblendet)")
    return pc


cbar = fig.colorbar(draw_field(), ax=ax_f, shrink=0.85)
snapshot()

# ---- Teil 1: finite Störung gegen Formel ------------------------------------------
rows, results = [], {n: [] for n in POINTS}
print(f"{'Punkt':8s} {'eps [mm]':>8} {'DOF':>7} {'dPi/(pi eps^2)':>15} {'Formel':>10} {'rel.Fehler':>10}  t[s]")
for name, x0 in POINTS.items():
    for e in EPS_LIST:
        t0 = time.perf_counter()
        msh, ct, ft = make_mesh(beam, disc=(*x0, e))
        u_ref, U_ref, W_ref = solve(beam, msh, ft, prefix="ref_")
        rho = rho_field(msh, ct, RHO0)
        u_inc, U_inc, W_inc = solve(beam, msh, ft, rho=rho, prefix="inc_")
        ratio = ((U_inc - W_inc) - (U_ref - W_ref)) / (np.pi * e**2)
        S = fem.functionspace(msh, ("DG", 2))
        dts = fem.Function(S)
        dts.interpolate(fem.Expression(dt_energy(beam.sigma(u_ref), beam.eps(u_ref), RHO0, beam.lam, beam.mu),
                                       S.element.interpolation_points))
        formula = float(eval_at(msh, dts, x0)[0])
        err = ratio / formula - 1
        results[name].append((e, ratio, formula))
        ndof = 2 * fem.functionspace(msh, ("Lagrange", 2)).dofmap.index_map.size_global
        rows.append(dict(part="check", point=name, x0=x0, eps=e, ndof=ndof, ratio=ratio, formula=formula,
                         rel_err=err, rho0=RHO0, t=time.perf_counter() - t0))
        print(f"{name:8s} {e * 1e3:8.2f} {ndof:7d} {ratio * J:15.3f} {formula * J:10.3f} {err:10.2e}  "
              f"{rows[-1]['t']:.1f}")

        # -- Bild: Ausschnitt um den Einschluss
        ax_z.clear()
        xy, tris, disp = triangulation_data(msh, u_inc)
        tri = Triangulation(*xy.T, tris)
        psi = cellwise(msh, 0.5 * ufl.inner(beam.sigma(u_inc, rho), beam.eps(u_inc))) * J
        ax_z.tripcolor(tri, facecolors=np.clip(psi, 1e-1, None), norm=LogNorm(), cmap="inferno",
                       shading="flat", rasterized=True)
        ax_z.triplot(tri, lw=0.15, color="w", alpha=0.35)
        ax_z.add_patch(plt.Circle(x0, e, fill=False, ec="c", lw=1.5))
        ax_z.set_xlim(x0[0] - 4 * e, x0[0] + 4 * e)
        ax_z.set_ylim(x0[1] - 3 * e, x0[1] + 3 * e)
        ax_z.set_aspect("equal")
        ax_z.set_title(f"{name}: Einschluss ε = {e * 1e3:.1f} mm, ρ₀ = {RHO0:g}\n"
                       "Energiedichte (Spannungen umfließen das Loch)", fontsize=10)
        # -- Bild: Konvergenz
        ax_c.clear()
        ax_e.clear()
        for n, res in results.items():
            if not res:
                continue
            ee, rr, ff = map(np.array, zip(*res))
            ax_c.semilogx(ee * 1e3, -rr * J, "o-", color=COLORS[n], label=f"{n}: FEM")
            ax_c.axhline(-ff[-1] * J, color=COLORS[n], ls="--", lw=1)
            ax_e.loglog(ee * 1e3, np.abs(rr / ff - 1), "o-", color=COLORS[n], label=n)
        ax_c.set_xlabel("ε [mm]")
        ax_c.set_ylabel("$-(\\Pi(\\varepsilon)-\\Pi(0))/(\\pi\\varepsilon^2)$ [J/m³]")
        ax_c.set_title("Finite Störung (Punkte) gegen Formel $-P_\\gamma\\sigma:\\varepsilon$ (gestrichelt)")
        ax_c.invert_xaxis()
        ax_c.legend(fontsize=9)
        ee = EPS_LIST * 1e3
        ax_e.loglog(ee, 0.4 * (ee / ee[0])**2, ":", color="0.5", label="Steigung 2")
        ax_e.axhspan(1e-4, 1e-2, color="0.85", alpha=0.4, lw=0)
        ax_e.text(ee[1], 1.3e-4, "grau: Auflösungsgrenze des Netzes", fontsize=8, color="0.4")
        ax_e.set_xlabel("ε [mm]")
        ax_e.set_ylabel("relativer Fehler")
        ax_e.set_title("Konvergenz gegen die Formel")
        ax_e.invert_xaxis()
        ax_e.legend(fontsize=9)
        draw_field(highlight=name)
        snapshot()

# ---- Teil 2: Nukleationslast, omega*, omega^beta -----------------------------------
ok = ~near_support
m = (-dtPi[ok]).max()
lam_c = np.sqrt(KAPPA / m)                     # D_T Pi skaliert mit q^2
q_c = beam.q * lam_c
lam_show = 1.15 * lam_c                        # Darstellung leicht über der Nukleationslast
dtF = lam_show**2 * dtPi + KAPPA
omega_star = ok & (dtF < 0)
order = np.argsort(np.where(omega_star, dtF, np.inf))
cum = np.cumsum(area[order])
n_beta = max(1, int(np.searchsorted(cum, np.pi * args.delta**2 / 4)))
omega_beta = np.zeros_like(omega_star)
omega_beta[order[:n_beta]] = True
omega_beta &= omega_star
alpha, beta = polarization_coefficients(RHO0, beam.lam, beam.mu)
ib = order[0]

draw_field()
ax_f.tricontour(Triangulation(mid[0], mid[1]), np.where(ok, -dtPi * lam_show**2 / KAPPA, 0), levels=[1.0],
                colors="w", linewidths=1.2)
ax_f.tripcolor(tri0, facecolors=np.where(omega_beta, 1.0, np.nan), cmap="autumn", shading="flat")
ax_f.annotate(f"erster Schaden $\\omega^\\beta$\n(tr σ {'> 0, Zug' if trs[ib] > 0 else '< 0, Druck!'})",
              (mid[0][ib], mid[1][ib]), xytext=(40, 60), textcoords="offset points", color="w", fontsize=9,
              arrowprops=dict(arrowstyle="->", color="w"))
ax_f.set_title(f"$-D_T\\Pi$ [J/m³] bei q = {beam.q} MN/m;  weiß: $\\omega^*$ bei {lam_show:.2f}·q,  "
               f"Nukleation ab $q_c$ = {q_c:.2f} MN/m", fontsize=10)

err_last = {n: r[-1][1] / r[-1][2] - 1 for n, r in results.items()}
fig.suptitle("Topologische Ableitung (Xavier et al., Gl. 3.3–3.4): Formel geprüft, erster Schritt von "
             f"Algorithmus 1   (δ = {args.delta * 100:.0f} cm, κ_δ = G_c/δ = {KAPPA * J / 1e3:.0f} kJ/m³)",
             fontsize=12)
snapshot()

print(f"\nPolarisationstensor (rho0 = {RHO0:g}, ebener Verzerrungszustand): alpha = {alpha:.4f}, beta = {beta:.4f}")
for n, e in err_last.items():
    rr = np.array(results[n])
    pre = slice(1, 5)   # 19.8 ... 5.7 mm: vor der Auflösungsgrenze
    p = np.polyfit(np.log(rr[pre, 0]), np.log(np.abs(rr[pre, 1] / rr[pre, 2] - 1)), 1)[0]
    best = np.abs(rr[:, 1] / rr[:, 2] - 1).min()
    print(f"  {n:8s}: Ordnung (19.8-5.7 mm) {p:.2f},  kleinster Fehler {best:.1e},  bei eps = 2.5 mm {e:+.1e}")
print(f"kappa_delta = G_c/delta = {KAPPA * J:.0f} J/m^3,  max(-D_T Pi) außerhalb der Auflager = {m * J:.1f} J/m^3")
print(f"Nukleationslast q_c = {q_c:.3f} MN/m  (= {lam_c:.2f} x q),  Ort ({mid[0][ib]:.3f}, {mid[1][ib]:.3f}) m, "
      f"tr(sigma) {'> 0 (Zug)' if trs[ib] > 0 else '< 0 (Druck)'}")
print(f"|omega^beta| = {area[omega_beta].sum() * 1e4:.3f} cm^2  (Ziel pi delta^2/4 = {np.pi * args.delta**2 / 4 * 1e4:.3f} cm^2), "
      f"{omega_beta.sum()} Zellen")
print(f"Max. -D_T Pi in den Auflagerzonen: {(-dtPi[near_support]).max() * J:.0f} J/m^3 "
      f"(-> ohne Ausblenden würde dort zuerst Schaden entstehen)")

rows.append(dict(part="nucleation", delta=args.delta, kappa=KAPPA, q=beam.q, q_c=q_c, x=mid[0][ib], y=mid[1][ib],
                 tr_sigma=trs[ib], omega_beta_area=float(area[omega_beta].sum()), alpha=alpha, beta=beta))
fig.savefig(OUT / "final.png", dpi=130)
with open(OUT / "runs.jsonl", "w") as fh:
    for r in rows:
        fh.write(json.dumps({k: (list(v) if isinstance(v, tuple) else v) for k, v in r.items()}) + "\n")
from PIL import Image  # noqa: E402

imgs = [Image.fromarray(f) for f in frames]
imgs[0].save(OUT / "verlauf.gif", save_all=True, append_images=imgs[1:] + [imgs[-1]] * 4, duration=500, loop=0)
print(f"\nGeschrieben: {OUT}/final.png, verlauf.gif, runs.jsonl")
if not args.no_show:
    plt.ioff()
    plt.show()
