# score_crispr_targets.py

Assume that the fasta file will have only one sequence (i.e. not a multi fasta)

>TEST_FASTA

# passing
TGACTACGCCTTTCTCTAGAGGG  # first place, GGG and 45%
CTTTACTTGACCAGATATACGGG  # second place, GGG and 40%
ACTCTTGAGTGTTAAAATCTTGG  # third place
TAATTCCTATATTTAAAATATGG  # fourth

```
jellyfish count -c 8 -s 1000 test.fasta
jellyfish dump mer_counts.js > counts.fa
```