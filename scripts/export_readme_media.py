"""Kopiert die Verlaufs-GIFs der Snippets verkleinert nach docs/media/ (für das README).

results/raw/ ist nicht versioniert. Nach einem neuen Lauf eines Snippets dieses Skript erneut ausführen:
    python scripts/export_readme_media.py
"""
from pathlib import Path

from PIL import Image, ImageSequence

ROOT = Path(__file__).resolve().parents[1]
WIDTH = 900
SOURCES = {
    "04_griffith_energy.gif": ROOT / "results/raw/04_griffith_energy/verlauf.gif",
    "05_topological_derivative.gif": ROOT / "results/raw/05_topological_derivative/verlauf.gif",
}

out_dir = ROOT / "docs" / "media"
out_dir.mkdir(parents=True, exist_ok=True)
for name, src in SOURCES.items():
    if not src.exists():
        print(f"fehlt: {src.relative_to(ROOT)} (Snippet zuerst ausführen)")
        continue
    im = Image.open(src)
    size = (WIDTH, round(im.height * WIDTH / im.width))
    frames = [f.convert("RGB").resize(size, Image.LANCZOS) for f in ImageSequence.Iterator(im)]
    durations = [f.info.get("duration", 500) for f in ImageSequence.Iterator(Image.open(src))]
    frames[0].save(out_dir / name, save_all=True, append_images=frames[1:], duration=durations,
                   loop=0, optimize=True)
    print(f"{name}: {len(frames)} Bilder, {(out_dir / name).stat().st_size / 1e6:.1f} MB")
