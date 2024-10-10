## help:      print this help message and exit
help: Makefile
	@echo ''
	@sed -n 's/^## //p' Makefile
	@echo ''

## test:      run the automated test suite and print coverage information
test:
	pytest --cov=lusSTR --doctest-modules lusSTR/tests/test_*.py

## style:     check code style
style:
	black --line-length=99 --check *.py lusSTR/cli/*.py lusSTR/scripts/*.py lusSTR/wrappers/*.py lusSTR/tests/test_*.py

## format:    auto-reformat code with Black
format:
	black --line-length=99 *.py lusSTR/cli/*.py lusSTR/scripts/*.py lusSTR/wrappers/*.py lusSTR/tests/test_*.py

## devenv:    configure a development environment
devenv:
	conda install black==22.6 pytest pytest-cov
	pip install -e .
	echo 'make style' > .git/hooks/pre-commit
	chmod 755 .git/hooks/pre-commit
