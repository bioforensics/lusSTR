## help:      print this help message and exit
help: Makefile
	@echo ''
	@sed -n 's/^## //p' Makefile
	@echo ''

## test:      run the automated test suite
test:
	pytest lusSTR/tests/test_suite.py

## testcov:   run the automated test suite and print coverage information
testcov:
	pytest --cov=lusSTR lusSTR/tests/test_suite.py

## style:     check code style against PEP8
style:
	pycodestyle --max-line-length=99 lusSTR/*.py

## devenv:    configure a development environment
devenv:
	conda install pycodestyle pytest pytest-cov
	pip install -e .
	echo 'make style' > .git/hooks/pre-commit
	chmod 755 .git/hooks/pre-commit
