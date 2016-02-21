
lint:
	flake8 score_crispr_targets.py
	flake8 score/*.py
	flake8 test/*.py

test: lint
	nosetests

run:
	python score_crispr_targets.py --fasta data/codingTaskSequence.fa --output out.bed --gff data/codingTaskAnnotation.gff3
	@echo "CRISPR score calls are in out.bed"