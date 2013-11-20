#!/bin/bash

TEST_PATH=$(pwd)
VIRT_ENV="$TEST_PATH/tmp_virt_env"
VIRT_ACTIVATE="$VIRT_ENV/bin/activate"
VIRT_PYTHON="$VIRT_ENV/bin/python"

echo "Creating virtual env..."
virtualenv $VIRT_ENV

echo "Activating virtual env..."
source $VIRT_ACTIVATE
pip install -r requirements-dev.txt

echo "Running tests..."
$VIRT_PYTHON setup.py nosetests
TESTS_STATUS=$?

echo "Removing virtualenv"
rm -r $VIRT_ENV

exit $TESTS_STATUS
