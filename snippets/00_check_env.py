"""Prüft die Umgebung: Pakete, PETSc-Löser, Ridge-Integrale (ufl.dr)."""
import importlib

for name in ["dolfinx", "ufl", "basix", "petsc4py", "mpi4py", "gmsh", "numpy",
             "scipy", "sympy", "pandas", "matplotlib", "pyvista", "yaml", "h5py",
             "pytest", "adios4dolfinx"]:
    try:
        mod = importlib.import_module(name)
        print(f"  ok     {name:14s} {getattr(mod, '__version__', '')}")
    except ImportError:
        print(f"  FEHLT  {name}")

from petsc4py import PETSc  # noqa: E402

print("PETSc", ".".join(map(str, PETSc.Sys.getVersion())),
      "| Index:", PETSc.IntType.__name__, "| Skalar:", PETSc.ScalarType.__name__)
for pkg in ["mumps", "superlu_dist", "hypre"]:
    print(f"  {pkg:13s}", PETSc.Sys.hasExternalPackage(pkg))

import ufl  # noqa: E402

print("ufl.dr (Ridge-Integral) vorhanden:", hasattr(ufl, "dr"))
