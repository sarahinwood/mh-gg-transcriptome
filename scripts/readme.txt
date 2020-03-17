de novo transcriptome assembly using Trinity:

Snakefile now contains all major steps in transcriptome assembly and analysis pipeline, with other scripts for analysis and/or plotting in R.

To run must create a virtual environment, and install via pip3:
-pathlib2
-pandas
-os
-trinotate pipeline (see https://github.com/TomHarrop/trinotate_pipeline)

Must add below to bin folder:
-transrate
-trinotate folder from lab folder (also export path to it via:
					bin_dir="$(readlink -f bin/trinotate/bin)"
					export PATH="${bin_dir}:${PATH}"
					)