#!/bin/bash -ue
python3 - <<'PY'
import sys, shutil
from pathlib import Path
sys.path.insert(0, "/workspaces/ngs-variant-calling/ngs-variant-calling")
import run_interactive_pipeline as ngs
Path("logs").mkdir(exist_ok=True)
ngs.index_reference(Path("NC_002549.1.fasta"), Path("logs"))
PY
