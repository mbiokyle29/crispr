
lint:
	flake8 score_crispr_targets.py
	flake8 score/*.py

test: lint