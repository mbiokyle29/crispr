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

from score.utils import (CRISPRArgumentParser, file_exists,
                         read_fasta_file, read_gff_file, write_bed_file)
from score.bio import (annotate_gene_overlaps, build_kmer_count,
                       filter_target, generate_targets, score_target)


def main():

    # configure the logger here
    # default to error, --verbose flag will set to INFO
    log = logging.getLogger(__name__)
    log.setLevel(logging.ERROR)
    log_formatter = logging.Formatter(
        'CRISPR_SCORE|%(levelname)s|: %(message)s'
    )

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.ERROR)
    stream_handler.setFormatter(log_formatter)
    log.addHandler(stream_handler)

    # import the custom parser
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

    # GFF is optional
    if args.gff is not None and not file_exists(args.gff):
        log.error("gff file %s does not exist -- Exiting", args.gff)
        raise SystemExit

    # set up data
    fasta = read_fasta_file(args.fasta)

    # the kmer_spectra is a dictionary mapping kmers to a count score
    # the score = (1 / count) * multiplyer, so a kmer which appears once
    # has a base score of one, whereas a more common one (5) would have .2
    kmer_size = args.seed_end - args.seed_start + 1
    kmer_specta = build_kmer_count(fasta.seq, k=kmer_size)
    targets = []

    # yields a Crispir target seq instance
    for target in generate_targets(fasta, target_length=args.length):
        log.info("Target found: %s", target)

        # filter and score
        if filter_target(target, gc_low=args.gc_low, gc_high=args.gc_high,
                         homopolymer_length=args.homopolymer):

            # pass the scoring parameters from args Namespace as a dict
            target.score = score_target(target, kmer_specta, **args.__dict__)
            targets.append(target)

    log.info("Filtering/Scoring complete, %i targets passed", len(targets))

    if len(targets) == 0:
        log.warn("No targets passed filtering, exiting!")
        return

    # optionally check the gff file
    # overlaps are appened to the CrisprTarget genes[] data member
    if args.gff:
        log.info("Reading in GFF annotations")
        gff = read_gff_file(args.gff)
        log.info("Found %i genes from %s", len(gff), args.gff)
        if len(gff) > 0:
            targets = annotate_gene_overlaps(targets, gff)

    # write the output as a bed file
    write_bed_file(targets, fasta, args.output)

if __name__ == "__main__":
    main()
