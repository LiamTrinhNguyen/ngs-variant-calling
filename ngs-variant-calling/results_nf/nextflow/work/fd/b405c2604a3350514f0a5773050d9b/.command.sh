#!/bin/bash -ue
python3 - <<'PY'
import sys, shutil
from pathlib import Path
sys.path.insert(0, "/workspaces/ngs-variant-calling/ngs-variant-calling")
import run_interactive_pipeline as ngs
import dashboard
# The checks were staged in a sub-folder (their names can't clash with FASTP's); move them up
for f in Path("02_qc/checks").glob("*"):
    shutil.copy(f.resolve(), Path("02_qc") / f.name)
dashboard.build_dashboard(Path("."), min_qual=30, min_dp=15,
                          min_af=0.8, samples=["SRR1553607"], controls=[])
PY
