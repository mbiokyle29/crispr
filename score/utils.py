"""
utils.py
author: Kyle McChesney

Utility functions and for scoring CRISPIR targets
"""
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from Bio import SeqIO

class CRISPRArgumentParser(ArgumentParser):

    def __init__(self, *args, **kwargs):

        super(CRISPRArgumentParser, self).__init__(*args, **kwargs)
        self.formatter_class = ArgumentDefaultsHelpFormatter
        self._optionals.title = 'Options'

        self.add_argument("--verbose", "-v", help="Set log to INFO",
                          action="store_true")

        # main args
        self.add_argument("--fasta", help="fasta file to search in",
                            required=True, metavar="fasta")
        self.add_argument("--length", default=20, metavar="len",
                            help="Target seq length (not PAM)")

        # filter config args
        self.add_argument("--gff", default="./data/codingTaskAnnotation.gff3",
                            help="GFF file for filtering (optional)",
                            metavar="gff")
        self.add_argument("--gc-low", help="GC percent bottom cutoff",
                            default=20, type=int, metavar="gc-low")
        self.add_argument("--gc-high", help="GC percent upper cutoff",
                            default=80, type=int, metavar="gc-high")
        self.add_argument("--homopolymer", default=5, type=int,
                            help="Homopolymer length cutoff (inclusive)",
                            metavar="homo-len")

        # score config args
        self.add_argument("--pam-gc-score", default=1, type=int,
                            help="Score added for a PAM starting with G/C",
                            metavar="pam-score")

        # Set a target GC content, defaults to 45%
        # and set a multiplier. We add |value-target| * target
        # to the score for a sequence
        self.add_argument("--gc-goal", default=0.45, type=float,
                            help="Target GC content percentage",
                            metavar="gc-goal")
        self.add_argument("--gc-multiplyer", default=1, type=int,
                            help="GC target multiplier: score * 1 / |val-target|",
                            metavar="gc-score")

        # seed region scoring
        # we can define the region and set a multiplyer
        self.add_argument("--seed-start", default=13, type=int,
                            help="1 based index for start of seed",
                            metavar="seed-start")
        self.add_argument("--seed-end", default=20, type=int,
                            help="1 based index for end of seed",
                            metavar="seed-end")
        self.add_argument("--uniqueness-method", choices=["count", "hamming"],
                            default="count", metavar="unique-meth",
                            help="Score uniqueness with counts or hamming dist")
        self.add_argument("--uniqueness-multiplyer", default=1, type=int,
                            help="Unique multiplier: mult * score",
                            metavar="unique-score")

def file_exists(file):
    fullpath = os.path.abspath(file)
    return os.path.exists(fullpath)

def read_fasta_file(fasta):

    # read it in
    # assume there is only one
    try:
        fasta_record = SeqIO.read(open(fasta), "fasta")
    except ValueError:
        log.error("Fasta file %s has more than one sequence -- Exiting", fasta)
        raise SystemExit

    return fasta_record