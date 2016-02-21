"""
bio.py
author: Kyle McChesney

Bioinformatics functions for crispr scoring
"""
import logging
import itertools
import re
from collections import defaultdict
from string import upper

from Bio.SeqRecord import SeqRecord

log = logging.getLogger("__main__")
DNA_ALPHABET = set(["A","T","C","G"])
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
        self.gc_content = round(content,0)

        # set up score
        self.score = -1

    def __str__(self):
        return "<CrisprTarget ({}): '{}' >".format(self.id, self.seq)

    def to_bed(self):
        return -1

def build_kmer_count(sequence, k=8):
    """
    sequence: a nuclecic acid sequence as string
    k: an integer value, the length of k-mers

    returns: a default dictionary keyed by k-mers
            with counts as values
    """
    spectra = defaultdict(int)
    
    if k < 1:
        raise ValueError("kmer size must be at least 1" )

    start_idx = 0
    stop_idx = k

    while stop_idx <= (len(sequence)):

        # update the count
        kmer = sequence[start_idx:stop_idx]
        spectra[str(kmer)] += 1

        # move pointer
        start_idx += 1
        stop_idx += 1

    for kmer in spectra:
        spectra[kmer] = float(len(spectra)) / spectra[kmer] 

    return spectra

def build_kmer_hamming(sequence, k=8):

    spectra = set()
    if k < 1:
        raise ValueError("kmer size must be at least 1" )

    start_idx = 0
    stop_idx = k

    # build the set of kmers
    while stop_idx <= (len(sequence)):
        kmer = sequence[start_idx:stop_idx]
        spectra.add(str(kmer))

        # move pointer
        start_idx += 1
        stop_idx += 1

    hamming_spectra = {}
    size = len(spectra)
    for kmer in spectra:
        others = spectra - set(kmer)
        average_hamming_dist = sum([fast_hamming(kmer, x) for x in others]) / size
        hamming_spectra[kmer] = average_hamming_dist

    return hamming_spectra

def generate_targets(sequence, target_length=20, pam="GG"):

    seq_string = str(sequence.seq)
    pam_idx = seq_string.find(pam)

    while pam_idx != -1:

        if pam_idx > target_length:

            # yield the target
            start = pam_idx - target_length - 1
            stop = pam_idx + len(pam)  # add the length of the pam
            name = "{}-CRISPR-TARGET({})".format(sequence.name, start)
            yield CrisprTarget(seq_string[start:stop], id=str(start), name=name)

        pam_idx = seq_string.find(pam, pam_idx+1)

def filter_target(target, gc_low=20, gc_high=80, homopolymer_length=5):
    
    # check GC
    if not target.gc_content in range(gc_low, gc_high):
        log.info("Target %s failed GC cutoff (%i%s)\n", target, target.gc_content,"%")
        return False

    # check ATG
    if "ATG" in target.seq:
        log.info("Target %s failed ATG motif check\n", target)
        return False

    # check homopolymer
    regex = "".join(["(.)\1{", str(homopolymer_length), ",}"])
    if re.search(regex, target.seq) != None:
        log.info("Target %s failed homopolymer check\n", target)
        return False

    log.info("Target %s passed filter", target)
    return True

def score_target(target, kmer_spectra, **scoring_params):

    score = 0
    log.info("Scoring: %s", target)

    # score the PAM region
    if target.seq[-3] in "GC":
        pam_gc_bonus = scoring_params.get('pam_gc_score') or 1
        log.info("Adding PAM region GC start bonus (%i)", pam_gc_bonus)
        score += pam_gc_bonus

    # GC scoring
    gc_target = scoring_params.get('gc_goal') or 45
    gc_multiplyer = scoring_params.get('gc_multiplyer') or 1
    gc_distance_measure = 1 / float(1 + abs(target.gc_content - gc_target))
    log.info("GC content distance measure: %.3f (value: %i, target: %i)",
             gc_distance_measure, target.gc_content, gc_target)
    score += gc_distance_measure * gc_multiplyer

    # seed region
    seed_start = scoring_params.get('seed_start') or 13
    seed_end = scoring_params.get('seed_end') or 20
    unique_multiplyer = scoring_params.get('uniqueness_multiplyer') or 1

    seed_sequence = target.seq[seed_start-1:seed_end]
    log.info("Seed: %s has a unique score of %f",
             seed_sequence, kmer_spectra[seed_sequence])
    score += kmer_spectra[seed_sequence] * unique_multiplyer
    
    return score

# http://jhafranco.com/2012/02/12/hamming-distance/
def fast_hamming(left, right):

    if len(left) != len(right):
        log.error("length of %s not equal to %s cannot compute hamming distance",
            left, right)
        raise ValueError("Strings must be equal length for hamming distance")

    left = int(dna_to_bit_string(left), 2)
    right = int(dna_to_bit_string(right), 2)
    xor = right^left
    count = 0

    while xor:
        count += 1
        xor &= xor-1

    return count

def dna_to_bit_string(sequence):

    if not set(upper(sequence)) <= DNA_ALPHABET:
        log.error("Can only convert DNA (ATCG) to a bit string")
        raise ValueError("Sequence with elements not in ATCG")

    for char in DNA_TO_BIN:
        sequence = sequence.replace(char, DNA_TO_BIN[char])

    return sequence
