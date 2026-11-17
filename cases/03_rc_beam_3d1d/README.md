# 03: Stahlbetonträger, Beton 3D + Bewehrung 1D, Strain Gradient

Beton als 3D-Volumen (Strain Gradient, C0-IP), Bewehrung als 1D-Fachwerkstab auf eingebetteten
Netzkanten (`ufl.dr`). Ausgangspunkt: `snippets/01_rebar_3d1d.py` und `02_strain_gradient.py`.
Kernfrage: Netzkonvergenz der Kopplung bei l = 0 (keine Konvergenz erwartet) gegenüber l > 0.
