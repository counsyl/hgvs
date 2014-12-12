#! /usr/bin/make

VENV_DIR?=.virtualenvs/default
SCRIPTS=bin/hgvs

default:
	python setup.py check build

.PHONY: setup clean teardown venv lint test gitlint

setup: venv

venv: $(VENV_DIR)/bin/activate

$(VENV_DIR)/bin/activate: requirements-dev.txt
	test -d $(VENV_DIR) || virtualenv --python=python2.7 --system-site-packages $(VENV_DIR)
	. $(VENV_DIR)/bin/activate; pip install -r requirements-dev.txt
	touch $(VENV_DIR)/bin/activate

teardown:
	rm -rf $(VENV_DIR)/

clean:
	python setup.py clean
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg*/
	rm -rf __pycache__/
	rm -f MANIFEST
	rm -f nosetests.xml
	find pyhgvs -type f -name '*.pyc' -delete

lint: venv
	. $(VENV_DIR)/bin/activate; flake8 --jobs=auto pyhgvs/ $(SCRIPTS)


test: venv
	. $(VENV_DIR)/bin/activate; nosetests --verbosity=2 --with-xunit pyhgvs
