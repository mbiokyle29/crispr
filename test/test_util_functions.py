"""
test_util_functions.py
author: Kyle McChesney
"""

from score.utils import file_exists, read_fasta_file, read_gff_file

import unittest


class testUtilFunctions(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.bad_file = "/home/bob/dna.fa"
        cls.invalid_fasta_file = "./test/data/util/invalid_fasta.fa"
        cls.valid_fasta_file = "./test/data/util/valid_fasta.fa"
        cls.invalid_gff_file = "./test/data/util/invalid_gff.gff3"
        cls.valid_gff_file = "./test/data/util/valid_gff.gff3"

    def test_file_doesnt_exist(self):
        self.assertFalse(file_exists(self.bad_file))

    def test_reading_valid_fasta(self):
        self.assertTrue(file_exists(self.valid_fasta_file))
        rec = read_fasta_file(self.valid_fasta_file)
        self.assertEqual(rec.id, "valid_fasta")

    def test_reading_invalid_fasta(self):
        self.assertTrue(file_exists(self.invalid_fasta_file))
        self.assertRaises(ValueError, read_fasta_file, self.invalid_fasta_file)

    def test_reading_valid_gff(self):
        self.assertTrue(file_exists(self.valid_gff_file))
        gff = read_gff_file(self.valid_gff_file)

        self.assertEqual(len(gff), 1)
        self.assertEqual(gff["ID=gene1;Name=GENE1"], (1000, 7000))

    def test_reading_invalid_gff(self):
        self.assertTrue(file_exists(self.invalid_gff_file))
        gff = read_gff_file(self.invalid_gff_file)

        # the gene line is invalid so its skipped
        self.assertEqual(len(gff), 0)
