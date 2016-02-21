#!/usr/bin/env python
"""
score_crispr_targets.py
author: Kyle McChesney

This script takes as input a fasta file, and optionally a GFF3 file.
It parses the fasta sequence, detects potential CRISPR target sequences.

It then filters them, and scores the passing targets and outputs them.
It sets sane defaults for everything, but everything can be configured.
The scoring is based on the requirements, but each part can be weighted more 
or less given its multiplier argument
"""
import logging
import re
from collections import defaultdict

from score.utils import CRISPRArgumentParser, file_exists, read_fasta_file
from score.bio import build_kmer_count, build_kmer_hamming, generate_targets, filter_target, score_target

def main():

    log = logging.getLogger(__name__)
    log.setLevel(logging.ERROR)
    log_formatter = logging.Formatter('CRISPR_SCORE|%(levelname)s|: %(message)s')
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.ERROR)
    stream_handler.setFormatter(log_formatter)
    log.addHandler(stream_handler)

    parser = CRISPRArgumentParser(
        description=("Find & profile CRISPR targets in a fasta file"),
    )

    args = parser.parse_args()
    
    if args.verbose:
        log.setLevel(logging.INFO)
        stream_handler.setLevel(logging.INFO)

    # validate files
    # fasta must be given
    if not file_exists(args.fasta):
        log.error("fasta file %s does not exist -- Exiting", args.fasta)
        raise SystemExit

    if args.gff is not None and not file_exists(args.gff):
        log.error("gff file %s does not exist -- Exiting", args.gff)
        raise SystemExit

    # set up data
    fasta = read_fasta_file(args.fasta)
    print fasta
    if args.uniqueness_method == "count":
        kmer_specta = build_kmer_count(fasta.seq)
    else:
        kmer_specta = build_kmer_hamming(fasta.seq)

    targets = []

    # yields a Crispir target seq instance
    for target in generate_targets(fasta):
        log.info("Target found: %s", target)
        if filter_target(target):
            
            target.score = score_target(target, kmer_specta)
            targets.append(target)

    for target in targets:
        print target.to_bed()


if __name__ == "__main__":
    main()
