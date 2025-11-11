# Configuration file

import os

REPO_PATH = "."
CODE_DIR = REPO_PATH + "/code"
SAMPLES = ['S01', 'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10',
           'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20', 'S21']
DEMULTIPLEX_MANIFEST = "annotation/demultiplex_manifest.csv"
RNA_REF = "annotation/snRNA_reference_v2.fasta"

DEMULTIPLEXER = "code/demultiplexer_v2.py"
TAIL_ANALYZER = "code/tail_analyzer_v2.py"

#RNAS = ['U1', 'U2', 'U3', 'U4', 'U5', 'U6']
#REPLICATES = ["R1", "R2", "R3"]
#PRIMER_IDS = [f"{RNA}-{REPLICATE}" for RNA in RNAS for REPLICATE in REPLICATES]






