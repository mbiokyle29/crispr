#!/usr/bin/env python
"""
score_crispir_targets.py
author: Kyle McChesney

This script takes as input a fasta file, and optionally a GFF3 file.
It parses the fasta sequence, detects potential CRISPIR target sequences.

It then filters them, and scores the passing targets and outputs them.
It sets sane defaults for everything, but everything can be configured.
The scoring is based on the requirements, but each part can be weighted more 
or less given its multiplier argument
"""
import logging
from utils import CRISPIRArgumentParser, file_exists

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
log_formatter = logging.Formatter('CRISPR_SCORE|%(levelname)s|: %(message)s')

stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.INFO)
stream_handler.setFormatter(log_formatter)

# set it all up
log.addHandler(stream_handler)

def main():

    parser = CRISPIRArgumentParser(
        description=("Find & profile CRISPIR targets in a fasta file"),
    )

    args = parser.parse_args()
    
    # validate files
    # fasta must be given
    if not file_exists(args.fasta):
        log.error("fasta file %s does not exist -- Exiting", args.fasta)
        raise SystemExit

    if args.gff is not None and not file_exists(args.gff):
        log.error("gff file %s does not exist -- Exiting", args.gff)
        raise SystemExit

if __name__ == "__main__":
    main()
