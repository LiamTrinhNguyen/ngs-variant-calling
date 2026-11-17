#!/bin/bash -ue
python3 - <<'PY'
import sys, shutil
from pathlib import Path
sys.path.insert(0, "/workspaces/ngs-variant-calling/ngs-variant-calling")
import run_interactive_pipeline as ngs
ngs.run_r_stats(Path("SRR1553607.all_calls.tsv"), Path("NC_002549.1.fasta"), Path("."),
                "SRR1553607", 1000)
PY
