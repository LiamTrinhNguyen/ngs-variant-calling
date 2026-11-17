#!/bin/bash -ue
python3 - <<'PY'
import sys, shutil
from pathlib import Path
sys.path.insert(0, "/workspaces/ngs-variant-calling/ngs-variant-calling")
import run_interactive_pipeline as ngs
ngs.align_and_sort(Path("NC_002549.1.fasta"), Path("SRR1553607_R1.trimmed.fastq.gz"), Path("SRR1553607_R2.trimmed.fastq.gz"),
                   Path("."), "SRR1553607", 2)
PY
