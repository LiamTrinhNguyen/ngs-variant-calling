#!/bin/bash -ue
python3 - <<'PY'
import sys, shutil
from pathlib import Path
sys.path.insert(0, "/workspaces/ngs-variant-calling/ngs-variant-calling")
import run_interactive_pipeline as ngs
ngs.filter_and_report(Path("SRR1553607.raw.vcf.gz"), Path("."), "SRR1553607",
                      min_qual=30, min_dp=15, min_af=0.8)
PY
