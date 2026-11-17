#!/bin/bash -ue
python3 - <<'PY'
import sys, shutil
from pathlib import Path
sys.path.insert(0, "/workspaces/ngs-variant-calling/ngs-variant-calling")
import run_interactive_pipeline as ngs
ngs.run_fastp(Path("SRR1553607_1.fastq.gz"), Path("SRR1553607_2.fastq.gz"), Path("."), "SRR1553607", 2)
PY
