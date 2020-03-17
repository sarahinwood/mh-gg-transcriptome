#!/usr/bin/env bash

set -eu

bin/bbmap/filterbyname.sh \
in=output/trinity/Trinity.fasta \
include=t \
names=output/trinotate/viral/viral_transcript_ids.txt \
out=output/trinotate/viral/viral_transcripts.fasta