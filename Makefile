#! /usr/bin/make

PACKAGE_NAME=pyhgvs
SCRIPTS=bin/hgvs
TEST_OUTPUT?=nosetests.xml

VENV_DIR?=.venv
VENV_ACTIVATE=$(VENV_DIR)/bin/activate
WITH_VENV=. $(VENV_ACTIVATE);

default:
	python setup.py check build

.PHONY: setup clean teardown venv lint test gitlint

setup: venv

venv: $(VENV_ACTIVATE)

$(VENV_DIR)/bin/activate: requirements-dev.txt
	test -d $(VENV_DIR) || virtualenv --python=python2.7 --system-site-packages $(VENV_DIR)
	$(WITH_VENV) pip install -r requirements-dev.txt
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
	rm -f $(TEST_OUTPUT)
	find $(PACKAGE_NAME) -type f -name '*.pyc' -delete

lint: venv
	$(WITH_VENV) flake8 --jobs=auto $(PACKAGE_NAME)/ $(SCRIPTS)

test: venv
	$(WITH_VENV) nosetests --verbosity=2 --with-xunit --xunit-file=$(TEST_OUTPUT)

package:
	python setup.py sdist
