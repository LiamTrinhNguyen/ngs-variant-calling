#!/bin/bash -ue
python3 - <<'PY'
import sys, shutil
from pathlib import Path
sys.path.insert(0, "/workspaces/ngs-variant-calling/ngs-variant-calling")
import run_interactive_pipeline as ngs
tmp = Path("download"); tmp.mkdir()
_, fastqs = ngs.get_sra_metadata_and_fastq("SRR1553607", tmp)
if len(fastqs) != 2:
    sys.exit("SRR1553607 is not paired-end; this workflow needs paired reads")
shutil.move(str(fastqs[0]), "SRR1553607_1.fastq.gz")
shutil.move(str(fastqs[1]), "SRR1553607_2.fastq.gz")
shutil.rmtree(tmp)
PY
