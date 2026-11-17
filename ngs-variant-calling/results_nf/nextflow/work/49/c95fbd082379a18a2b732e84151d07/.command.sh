#!/bin/bash -ue
python3 - <<'PY'
import sys, shutil
from pathlib import Path
sys.path.insert(0, "/workspaces/ngs-variant-calling/ngs-variant-calling")
import run_interactive_pipeline as ngs
ngs.screen_sample(Path("SRR1553607_R1.trimmed.fastq.gz"), Path("."), "SRR1553607", 2, expected_ref="NC_002549.1")
PY
