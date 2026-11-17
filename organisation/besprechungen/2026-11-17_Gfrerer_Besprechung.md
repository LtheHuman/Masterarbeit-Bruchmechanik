# Besprechung mit Prof. Gfrerer, 17.11.2026

Scan: `2026-11-17_Gfrerer_Besprechung.pdf` (2 Seiten). Abschrift und Deutung von Claude,
bitte gegen den Scan prüfen.

## Seite 1: Griffith-Energie und topologische Ableitung

Balken unter Gleichlast, kleiner Riss der Länge `a` am unteren Rand in Feldmitte.

$$E(a) = \tfrac12\int_\Omega \sigma(a):\varepsilon(a)\,dx + G_c\,a - W, \qquad W = \int_\Omega f\cdot v\,dx$$

- `a = 0`: intakt.
- `E(a) < E(0)` → Riss entsteht.
- Kriterium: $\lim_{a\to 0} \dfrac{E(a)-E(0)}{a^2} < 0$. Das ist die topologische Ableitung.
  Die Skalierung $a^2$ entspricht dem 2D-Fall (Flächenmaß des Einschlusses, vgl. $\pi r^2$ auf Seite 2).
- Offene Frage: **„Was passiert bei Strain-Gradient-Theorie?“**
- Randnotiz: **Diese Theorie hilft bei der Kopplung von 3D-Element (Beton-Volumenkörper)
  mit 1D-Element (Bewehrung).**
- Randnotiz: **Riss als 3D-Volumen einführen.**

## Seite 2: Skalierung und Modell

- Riemannsumme $\sum_{i=1}^{\ell/h} f(hi)\,h = \sum (hi)^p h = \sum h^{p+1} i^p$: Wie skaliert
  eine Größe mit der Elementgröße `h`, wenn der Integrand singulär ist ($p<0$)?
  Notizen daneben: $2\pi r$, $\pi r^2$ (Umfang bzw. Fläche des Einschlusses).
- Balken mit Kerbe der Größe `h` am unteren Rand, Riss wandert nach oben (Rissverlauf skizziert).
- Konstitutivgesetz $\sigma = C:\varepsilon + (\ldots):\nabla\varepsilon$ (Strain-Gradient-Term,
  im Scan nicht ganz lesbar), $\varepsilon = \tfrac12(\nabla u + \nabla u^T)$.
- 3D-Balken mit Bewehrungsstab (blau, 1D) knapp über dem unteren Rand und keilförmigem Riss
  als 3D-Teilvolumen.

## Einordnung (Claude)

- Warum Strain Gradient für die 3D-1D-Kopplung hilft: Für $u \in H^1(\Omega)$ in 3D existiert
  keine Spur auf einer Kurve (Kodimension 2). Die Stabenergie $EA\,(\partial_s u_t)^2$ ist
  dann nicht definiert, ein Linienlager erzeugt eine log-singuläre Lösung, die bei
  Netzverfeinerung nicht konvergiert. Mit Strain Gradient ist $u \in H^2(\Omega)$, die Spur
  auf der Kurve liegt in $H^1$ und $\partial_s u$ in $L^2$. Die Stabenergie ist damit
  wohldefiniert. Das ist zu verifizieren: Spursatz für Kodimension 2.
- Paper `literatur/pdf/Topological_derivative_based_fracture_mo.pdf` (Xavier et al.):
  Francfort-Marigo-Funktional, Schädigung als weiche Phase $\rho_0 C$, Algorithmus 1
  (Nukleation über $D_T F < 0$, Netzverfeinerung an der Rissspitze). Bisher 2D, P1-Dreiecke.
  Die Arbeit erweitert das auf 3D und Strain Gradient.

## Offene Fragen an Prof. Gfrerer

1. Welches Strain-Gradient-Modell genau: Mindlin, vereinfacht nach Altan–Aifantis
   ($\ell^2\,\nabla\varepsilon$), oder nur volumetrisch?
2. Topologische Ableitung für Strain Gradient: Gibt es einen Polarisationstensor aus der
   Literatur oder ist die Herleitung Teil der Arbeit?
3. Bewehrung: nur axial (Fachwerkstab) oder auch Biegung? Verbund starr oder mit Schlupf?
4. Ist ein FEniCSx-basierter Workflow ok, oder nutzt die Gruppe FreeFEM + autofreefem als
   Standard?
