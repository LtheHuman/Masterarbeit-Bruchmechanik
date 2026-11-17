# Masterarbeit Bruchmechanik

Rissbildung über die **topologische Ableitung** (Riss als 3D-Teilvolumen) und
**Stahlbetonträger** als 3D-Volumenkörper mit 1D-Bewehrung, regularisiert über die
**Strain-Gradient-Theorie**. Betreuer: Prof. M. Gfrerer (TU Graz, Institut für Baumechanik).

Framework: **FEniCSx (DOLFINx 0.11)** mit PETSc/MUMPS. Die Begründung steht in `CONTEXT.md`.

## Ordnerstruktur

```
README.md            dieses Dokument
CONTEXT.md           Projektstand zum Wiedereinstieg (laufend aktualisieren)
environment.yml      conda-Umgebung (einzige Quelle für Versionen)
organisation/        Besprechungen (Scans + Abschrift), Zeitplan
literatur/           pdf/ + literature.bib (wird von thesis/ mitbenutzt)
vorarbeiten/         Masterprojekt (PDF) und codeBaseStiffener (FreeFEM), nur lesen
theory/              Herleitungen (Markdown/SymPy): Griffith-Energie, top. Ableitung, Strain Gradient
src/bruchmechanik/   gemeinsamer Python-Code (Formen, Material, Algorithmus, I/O)
meshes/              Gmsh-Generatoren (*.py); erzeugte *.msh werden nicht versioniert
snippets/            kleine Machbarkeitstests, einzeln lauffähig
cases/               ein Ordner je Modellproblem: case.yaml + run.py + README.md
  01_mode1_2d/         Kalibrierung gegen Xavier et al., Abschn. 5.1 (2D, Mode I)
  02_beam_crack_3d/    Riss als 3D-Teilvolumen im Biegebalken
  03_rc_beam_3d1d/     Beton 3D + Bewehrung 1D, Strain Gradient
tests/               pytest: Regressionssollwerte, FD-Check der topologischen Ableitung
results/raw/         Rohdaten je Lauf (JSONL), nicht versioniert
results/curated/     kuratierte CSV, versioniert, Quelle für alle Abbildungen
results/vtk/         ParaView-Ausgaben, nicht versioniert
scripts/             curate.py, make_figures.py, make_tables.py → thesis/figures, thesis/tables
thesis/              LaTeX, versioniert
```

## Quickstart

```bash
conda activate fenicsx-env          # vorhandene Umgebung, reicht für snippets/
python snippets/00_check_env.py
python snippets/01_rebar_3d1d.py
python snippets/02_strain_gradient.py
python snippets/03_local_refinement.py
```

Komplette Arbeitsumgebung: `mamba env create -f environment.yml && conda activate bruchmechanik`.

## Konventionen (Lehren aus dem Masterprojekt)

- **Git ab Tag 1, auch für `thesis/`.** Keine Ordnerkopien wie `alt/` oder `Version_DATUM/`.
- **Netze haben genau eine Quelle:** Generatorskripte in `meshes/`. Netze werden erzeugt, nicht kopiert.
- **Parameter stehen in `cases/*/case.yaml`,** nicht im Code. Jeder Lauf schreibt eine JSON-Zeile
  (Parameter, Versionen, Zeiten, Ergebnisgrößen) nach `results/raw/`.
- **Abbildungen und Tabellen werden nur aus `results/curated/` per Skript erzeugt.** Zahlen
  im Text gegen die kuratierten Daten prüfen.
- **Jeder neue Term bekommt einen Test:** Regressionssollwert oder Vergleich mit analytischer
  Lösung. Topologische Ableitungen per Finite-Differenzen-Störung prüfen (wie `checkGradient.edp`).
- **Linux-Pfade und Umgebungsdatei** statt PowerShell-Befehlen mit absoluten Windows-Pfaden.
- **Direkte Löser:** MUMPS. UMFPACK läuft über PETSc mit 32-Bit-Indizes und bricht bei
  ca. 3 GB ab (Masterprojekt, Abschn. 4.4.2). SuperLU nicht verwenden.
