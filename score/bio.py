"""
bio.py
author: Kyle McChesney

Bioinformatics functions for crispr scoring
"""
import logging
from collections import defaultdict
from string import upper

from Bio.SeqRecord import SeqRecord

log = logging.getLogger("__main__")
DNA_ALPHABET = set(["A", "T", "C", "G"])
DNA_TO_BIN = {
    "A": "00",
    "T": "01",
    "C": "10",
    "G": "11"
}


class CrisprTarget(SeqRecord):

    def __init__(self, sequence, **kwargs):
        SeqRecord.__init__(self, sequence, **kwargs)

        sequence = upper(sequence)
        g_count = sequence.count('G')
        c_count = sequence.count('C')

        # make GC content as a percent
        content = sum([g_count, c_count]) / float(len(sequence))
        content = content * 100
        self.gc_content = round(content, 0)

        # set the start stop in parent seq
        self.start = int(kwargs.get('id', -1))
        self.stop = self.start + len(sequence) - 1

        # set up score
        self.score = -1
        self.genes = []

    def __str__(self):

        if self.score != -1:
            return "<{}: score: {} '{}' >".format(self.name, self.score,
                                                  self.seq)
        else:
            return "<{}: '{}' >".format(self.name, self.seq)

    def to_bed(self):
        source_seq = self.name.split("-")[0]
        start = self.id
        stop = len(self.seq) + int(start)

        pam = self.seq[-3:]
        target = self.seq[:-3]
        name = "target:\'{}\';pam:\'{}\'".format(target, pam)

        if len(self.genes) > 0:
            name += ";genes:[{}]".format(",".join(self.genes))

        return "\t".join([source_seq, str(start),
                          str(stop), name, str(self.score)])


def build_kmer_count(sequence, k=8):
    """
    sequence: a nuclecic acid sequence as string
    k: an integer value, the length of k-mers

    returns: a default dictionary keyed by k-mers
            with counts as values
    """
    spectra = defaultdict(int)

    if k < 1:
        raise ValueError("kmer size must be at least 1")

    start_idx = 0
    stop_idx = k

    while stop_idx <= (len(sequence)):

        # update the count
        kmer = sequence[start_idx:stop_idx]
        spectra[str(kmer.seq)] += 1

        # move pointer
        start_idx += 1
        stop_idx += 1

    for kmer in spectra:
        spectra[kmer] = (1 / float(spectra[kmer]))

    return spectra


def generate_targets(sequence, target_length=20, pam="GG"):

    seq_string = str(sequence.seq)
    pam_idx = seq_string.find(pam)

    while pam_idx != -1:

        if pam_idx > target_length:

            # yield the target
            start = pam_idx - target_length - 1
            stop = pam_idx + len(pam)  # add the length of the pam
            name = "{}-CRISPR-TARGET({})".format(sequence.name, start)
            yield CrisprTarget(seq_string[start:stop],
                               id=str(start), name=name)

        pam_idx = seq_string.find(pam, pam_idx+1)


def filter_target(target, gc_low=20, gc_high=80, homopolymer_length=5):

    # check GC
    if target.gc_content not in range(gc_low, gc_high):
        log.info("Target %s failed GC cutoff (%i%s)",
                 target, target.gc_content, "%")
        return False

    # check ATG
    if "ATG" in target.seq:
        log.info("Target %s failed ATG motif check", target)
        return False

    # check homopolymer
    homopolymers = [base * homopolymer_length for base in DNA_TO_BIN]
    for polymer in homopolymers:
        if polymer in target.seq:
            log.info("Target %s failed homopolymer check", target)
            return False

    log.info("Target %s passed filter", target)
    return True


def score_target(target, kmer_spectra, **scoring_params):

    score = 0
    log.info("Scoring: %s", target)

    # score the PAM region
    if target.seq[-3] in "GC":
        pam_gc_bonus = scoring_params.get('pam_gc_score', 1)
        log.info("Adding PAM region GC start bonus (%i)", pam_gc_bonus)
        score += pam_gc_bonus

    # GC scoring
    gc_target = scoring_params.get('gc_goal', 45)
    gc_multiplyer = scoring_params.get('gc_multiplyer', 1)
    gc_distance_measure = 1 / float(1 + abs(target.gc_content - gc_target))
    log.info("GC content distance measure: %.3f (value: %i, target: %i)",
             gc_distance_measure, target.gc_content, gc_target)
    score += gc_distance_measure * gc_multiplyer

    # seed region
    seed_start = scoring_params.get('seed_start', 13)
    seed_end = scoring_params.get('seed_end', 20)
    unique_multiplyer = scoring_params.get('uniqueness_multiplyer', 1)

    seed_sequence = target.seq[seed_start-1:seed_end]
    log.info("Seed: %s has a unique score of %f",
             seed_sequence, kmer_spectra[seed_sequence])
    score += kmer_spectra[seed_sequence] * unique_multiplyer

    return round(score, 2)


def annotate_gene_overlaps(targets, gff):

    for target in targets:
        for gene in gff:

            # note, GFF is 1 based
            gene_start, gene_stop = gff[gene]

            # BED is 0 so +1 to target bounds
            target_start = target.start + 1
            target_stop = target.stop + 1

            if target_start in range(gene_start, gene_stop):
                log.info("Target {} overlaps with gene: {}".format(target,
                                                                   gene))
                target.genes.append(gene)
            elif target_stop in range(gene_start, gene_stop):
                log.info("Target {} overlaps with gene: {}".format(target,
                                                                   gene))
                target.genes.append(gene)

    return targets
