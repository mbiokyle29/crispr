"""
test_bio_functions.py
author: Kyle McChesney
"""
import unittest

from Bio.SeqUtils import GC as biopythonGC

from score.utils import read_fasta_file, read_gff_file
from score.bio import (annotate_gene_overlaps, CrisprTarget,
                       build_kmer_count, filter_target,
                       generate_targets, score_target)


class testBioFunctions(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.fasta_file = "./test/data/bio/test.fa"
        cls.gff_file = "./test/data/bio/test.gff3"
        cls.kmer_counts = "./test/data/bio/jellyfish_kmers.tsv"
        cls.num_targets = 7

    def setUp(self):
        self.fasta = read_fasta_file(self.fasta_file)
        self.gff = read_gff_file(self.gff_file)

    def test_kmer_building(self):

        # pull in kmer counts from jellyfish
        # as the True values
        valid_counts = {}
        with open(self.kmer_counts, "r") as fh:
            for line in fh:
                line = line.rstrip().split("\t")
                valid_counts[line[1]] = int(line[0])

        test_kmers = build_kmer_count(self.fasta)
        self.assertEqual(len(test_kmers), len(valid_counts))

        # retransform the data
        for kmer in test_kmers:
            count = 1 / test_kmers[kmer]
            self.assertEqual(count, valid_counts[kmer])

    def test_generate_targets(self):

        for target in generate_targets(self.fasta):
            self.assertEqual(len(target.seq), 23)
            self.assertTrue(target.seq.endswith("GG"))

    def test_filter_target(self):

        for target in generate_targets(self.fasta):

            # check if it is invalid
            valid = True
            homopolymers = [x*5 for x in ["A", "T", "C", "G"]]

            for polymer in homopolymers:
                if polymer in target:
                    valid = False

            if "ATG" in target:
                valid = False

            gc = biopythonGC(target.seq)
            if not 20 <= gc <= 80:
                valid = False

            status = filter_target(target)
            self.assertEqual(valid, status)

    def test_score_target(self):

        # fake this for now
        kmer_spectra = {'ATATATAT': 1, 'TATATATT': 1, 'GCGCGCGC': 1}

        # test PAM bonus
        score_params = {'pam_gc_score': 100,
                        'gc_multiplyer': 0,
                        'uniqueness_multiplyer': 0}

        pam_bonus = CrisprTarget("ATATATATATATATATATATCGG", id="1")
        no_pam_bonus = CrisprTarget("ATATATATATATATATATATTGG", id="1")

        self.assertEqual(score_target(pam_bonus, kmer_spectra,
                                      **score_params), 100)
        self.assertEqual(score_target(no_pam_bonus, kmer_spectra,
                                      **score_params), 0)

        # test GC stuff
        score_params['gc_goal'] = 0
        score_params['gc_multiplyer'] = 100
        no_gc = CrisprTarget("ATATATATATATATATATATATA", id="3")
        gc = CrisprTarget("GCGCGCGCGCGCGCGCGCGCGCG", id="4")
        no_gc_score = score_target(no_gc, kmer_spectra, **score_params)
        gc_score = score_target(gc, kmer_spectra, **score_params)
        self.assertTrue(no_gc_score < gc_score)

        # test uniqueness
        score_params['pam_gc_score'] = 0
        score_params['gc_multiplyer'] = 0
        score_params['uniqueness_multiplyer'] = 100
        kmer_spectra = {'AAAATTTT': 100, 'CCCCGGGG': 1}

        unique = CrisprTarget("TTATATATATATAAAATTTTGGG", id="5")
        not_unique = CrisprTarget("TTATATATATATCCCCGGGGGGG", id="6")
        unique_score = score_target(unique, kmer_spectra, **score_params)
        not_unique_score = score_target(not_unique, kmer_spectra,
                                        **score_params)
        self.assertTrue(unique_score > not_unique_score)

    def test_annotate_gene_overlaps(self):

        # the gene is at 25 - 50 in one index

        # from 0 --> 22
        left_of_gene = CrisprTarget("TGACTACGCCTTTCTCTAGAGGG", id="0")
        self.assertEqual(left_of_gene.start, 0)
        self.assertEqual(left_of_gene.stop, 22)

        # this should not get marked
        res = annotate_gene_overlaps([left_of_gene], self.gff)
        self.assertEqual(len(res[0].genes), 0)

        # from 2 --> 24 (0 index)
        # is 3 --> 25 (1 index)
        one_overlap = CrisprTarget("TGACTACGCCTTTCTCTAGAGGG", id="2")
        self.assertEqual(one_overlap.start, 2)
        self.assertEqual(one_overlap.stop, 24)
        res = annotate_gene_overlaps([one_overlap], self.gff)
        self.assertEqual(len(res[0].genes), 1)

        # marked
        all_in = CrisprTarget("TGACTACGCCTTTCTCTAGAGGG", id="24")
        self.assertEqual(all_in.start, 24)
        self.assertEqual(all_in.stop, 46)
        res = annotate_gene_overlaps([all_in], self.gff)
        self.assertEqual(len(res[0].genes), 1)

        # not marked
        right_of_gene = CrisprTarget("TGACTACGCCTTTCTCTAGAGGG", id="50")
        self.assertEqual(right_of_gene.start, 50)
        self.assertEqual(right_of_gene.stop, 72)
        res = annotate_gene_overlaps([right_of_gene], self.gff)
        self.assertEqual(len(res[0].genes), 0)
