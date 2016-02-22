# CRISPR Scoring

This program (`score_crispr_targets.py`) takes as input a fasta file (with a single sequence)
And outputs scored CRISPR targets, in a bed file format. In it's default running, it follows
the problem specification. However, there are many optional parameters that can be passed to the script, to allow for tweaking/etc.

#### Simple example:
This will run the script with the supplied data, and output to `crispr_targets.bed`
```bash
python score_crispr_targets.py --fasta data/codingTaskSequence.fa  --gff data/codingTaskAnnotation.gff3
```

#### More Complicated examples
```bash

# place more weight on uniqueness
python score_crispr_targets.py --fasta data/codingTaskSequence.fa  \
--gff data/codingTaskAnnotation.gff3 \
--uniqueness-multiplyer 4

# change the cutoff length of homopolyers
python score_crispr_targets.py --fasta data/codingTaskSequence.fa  \
--gff data/codingTaskAnnotation.gff3 \
--homopolymer 10

# relax the GC filter range
# change the goal GC level
# and place less scoring weight on GC
python score_crispr_targets.py --fasta data/codingTaskSequence.fa  \
--gff data/codingTaskAnnotation.gff3 \
--gc-low 15 \
--gc-high 95 \
--gc-goal 50 \
--gc-multiplyer .5
```


The `Makefile` provides three commands:
1. `lint` - Run flake8
2. `test` - Run `lint` as well as the unittest
3. `run`  - Run the script on the given data files, with parameters to the specification
