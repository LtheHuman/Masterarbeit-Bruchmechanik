# Besprechung mit Prof. Gfrerer, 17.11.2026

Abschrift von `2026-11-17_Gfrerer_Besprechung.pdf` (2 Seiten).

## Seite 1: Griffith-Energie und topologische Ableitung

Balken unter Gleichlast, Riss der Länge $a$ am unteren Rand in Feldmitte.

$$E(a) = \tfrac12\int_\Omega \sigma(a):\varepsilon(a)\,dx + G_c\,a - W, \qquad W = \int_\Omega f\cdot v\,dx$$

- $a = 0$: intakt; $E(a) < E(0)$: Riss entsteht.
- Kriterium $\lim_{a\to 0} \dfrac{E(a)-E(0)}{a^2} < 0$ (topologische Ableitung, 2D).
- Frage: Was passiert bei Strain-Gradient-Theorie?
- Strain Gradient für die Kopplung von 3D-Volumen (Beton) und 1D-Element (Bewehrung).
- Riss als 3D-Volumen einführen.

## Seite 2: Skalierung und Modell

- Riemannsumme $\sum_{i=1}^{\ell/h} f(hi)\,h = \sum (hi)^p h = \sum h^{p+1} i^p$: Skalierung mit der
  Elementgröße $h$ bei singulärem Integranden ($p<0$). Daneben $2\pi r$, $\pi r^2$.
- Balken mit Kerbe der Größe $h$ am unteren Rand, Rissverlauf nach oben.
- $\sigma = C:\varepsilon + (\ldots):\nabla\varepsilon$, $\varepsilon = \tfrac12(\nabla u + \nabla u^T)$.
- 3D-Balken mit Bewehrungsstab (1D) knapp über dem unteren Rand, keilförmiger Riss als 3D-Teilvolumen.

## Anmerkungen

- Für $u \in H^1(\Omega)$, $\Omega\subset\mathbb{R}^3$, gibt es keine Spur auf einer Kurve; die Stabenergie
  $EA\,(\partial_s u_t)^2$ ist dann nicht definiert. Mit $u \in H^2(\Omega)$ liegt die Spur in $H^1$.
- Xavier et al. (2017): Francfort-Marigo-Funktional, Schädigung als weiche Phase $\rho_0 C$,
  Nukleation und Fortschritt allein über $D_T F < 0$ (Algorithmus 1), 2D, P1.
