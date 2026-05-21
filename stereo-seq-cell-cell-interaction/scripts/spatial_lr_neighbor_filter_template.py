#!/usr/bin/env python3
"""Wrapper for the shared spatial LR neighbor filtering template.

The executable implementation lives in the spatial-communication skill because
that skill owns spatial constraints and interpretation. This wrapper keeps the
same article-derived template discoverable from the CCI skill.
"""

from __future__ import annotations

import runpy
from pathlib import Path


if __name__ == "__main__":
    script = (
        Path(__file__).resolve().parents[1].parent
        / "stereo-seq-spatial-communication"
        / "scripts"
        / "spatial_lr_neighbor_filter_template.py"
    )
    if not script.exists():
        raise SystemExit(f"Shared spatial LR template not found: {script}")
    runpy.run_path(str(script), run_name="__main__")
