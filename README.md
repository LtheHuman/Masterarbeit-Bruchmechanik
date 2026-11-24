# Masterarbeit Bruchmechanik

Rissentstehung über die topologische Ableitung (Riss als Teilvolumen) und Stahlbetonträger als
3D-Volumen mit 1D-Bewehrung. Betreuung: Prof. M. Gfrerer, TU Graz, Institut für Baumechanik.

Implementierung in FEniCSx (DOLFINx 0.11, PETSc/MUMPS). Ausschlaggebend waren die Laufzeiten für
Elastizität im Benchmark des Masterprojekts (`vorarbeiten/`) und die Ridge-Integrale (`ufl.dr`) für die
3D-1D-Kopplung.

## Stand

### Griffith-Energie (`snippets/04_griffith_energy.py`)

$E(a) = \tfrac12\int_\Omega \sigma:\varepsilon\,dx + G_c\,a - W$ für einen 2D-Einfeldträger
(ebener Verzerrungszustand, $q = 0{,}3$ MN/m, $G_c = 100$ N/m) mit Randkerbe der Länge $a$ in
Feldmitte, je Risslänge neu vernetzt.

![Griffith-Energie](docs/media/04_griffith_energy.gif)

- Clapeyron ($W = 2U$) auf $10^{-12}$ erfüllt; $-\Delta\Pi \sim a^2$, FEM/LEFM 0,97–1,02 für $a \le 9$ cm.
- Energiebarriere bei $a^* \approx 12$ cm, $E(a) < E(0)$ ab $a_c \approx 19{,}8$ cm.
- $(E(a)-E(0))/a^2 = G_c/a - C \to +\infty$ für $a\to 0$: Mit linienförmiger Dissipation keine Nukleation
  aus dem intakten Zustand.

### Topologische Ableitung nach Xavier et al. (`snippets/05_topological_derivative.py`)

Prüfung von $D_T\Pi = -P_\gamma\,\sigma:\varepsilon$ (Gl. 3.3–3.4) durch finite Störung mit weichem
Kreis ($\rho_0 = 10^{-6}$, Radius $\varepsilon\to0$) in Zug-, Schub- und Druckzone, danach ein Schritt
von Algorithmus 1 mit $\kappa_\delta = G_c/\delta$, $\delta = 1$ cm.

![Topologische Ableitung](docs/media/05_topological_derivative.gif)

- Relativer Fehler $\sim\varepsilon^2$ bis zur Auflösungsgrenze des Netzes, dort 0,2–0,3 % in Zug- und
  Druckzone (Lamé-Konstanten des ebenen Verzerrungszustands).
- Nukleation ab $q_c \approx 4{,}9\,q$ unten in Feldmitte; $|\omega^\beta| = \pi\delta^2/4$ aufgelöst.

### Machbarkeit 3D und 3D-1D (`snippets/01`–`03`)

Kragträger 2 × 0,3 × 0,5 m, Bewehrung Ø 20 mm, P2.

| | $w(L)$ [m] |
|---|---|
| ohne Bewehrung | −4,430e-3 |
| mit Bewehrung (`ufl.dr`) | −4,321e-3 |
| mit Bewehrung, Strain Gradient $\ell = 0{,}02$ m (C0-IP) | −4,214e-3 |
| mit Bewehrung, Strain Gradient $\ell = 0{,}05$ m (C0-IP) | −3,759e-3 |

Stablänge und Stabdehnung eines linearen Feldes über `dr` exakt. Lokale Verfeinerung in 3D (Plaza)
funktioniert; `ridge_tags` werden von `refine()` nicht übertragen.

## Offene Punkte

1. Zug und Druck: Die Druckzone hat fast dieselbe topologische Ableitung wie die Zugzone, das
   Francfort-Marigo-Funktional unterscheidet nicht. Energiesplit (spektral, volumetrisch-deviatorisch)?
2. An den Auflagerkanten ist $-D_T\Pi$ über 100-mal größer; derzeit wie im Paper ausgeblendet.
3. Das Maximum liegt am Rand ($y = 0$), die Formel gilt für innere Punkte. Polarisationstensor für
   Randpunkte.
4. Polarisationstensor in 3D.

## Nächste Schritte

1. Algorithmus 1 vollständig, Validierung an Xavier et al., Abschnitt 5.1 (`cases/01_mode1_2d`).
2. Verfügbare Methoden zur Rissentstehung am selben Testproblem vergleichen.

## Aufbau

```
src/bruchmechanik/      beam2d (2D-Einfeldträger), beam3d (3D-Kragträger mit Bewehrung),
                        topological_derivative (Polarisationstensor, D_T Pi)
snippets/               einzeln lauffähige Studien, Ausgabe nach results/raw/ (nicht versioniert)
cases/                  Modellprobleme (in Arbeit)
tests/                  pytest
scripts/                export_readme_media.py: GIFs nach docs/media/
organisation/           Besprechungen
literatur/              Literatur und literature.bib
vorarbeiten/            Masterprojekt, codeBaseStiffener (FreeFEM)
```

## Ausführen

```bash
mamba env create -f environment.yml
conda activate bruchmechanik
python snippets/00_check_env.py
python snippets/05_topological_derivative.py     # --no-show: ohne Fenster
python -m pytest
```
