test:
	pytest lusSTR/test_suite.py

testcov:
	pytest --cov=lusSTR lusSTR/test_suite.py

style:
	pycodestyle --max-line-length=99 *.py
