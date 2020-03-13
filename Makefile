test:
	pytest test_suite.py

testcov:
	pytest --cov=STR_annotation test_suite.py

style:
	pycodestyle --max-line-length=99 *.py
