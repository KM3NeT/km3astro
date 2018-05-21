PKGNAME=km3astro
ALLNAMES = $(PKGNAME)

default: build

all: install

build: 
	@echo "No need to build anymore :)"

install: 
	pip install -U numpy
	pip install .

install-dev:
	pip install -U numpy
	pip install -e .

clean:
	python setup.py clean --all
	rm -f $(PKGNAME)/*.cpp
	rm -f $(PKGNAME)/*.c
	rm -f -r build/
	rm -f $(PKGNAME)/*.so

test: 
	py.test --junitxml=./reports/junit.xml $(PKGNAME)

test-cov:
	py.test --cov ./ --cov-report term-missing --cov-report xml:reports/coverage.xml --cov-report html:reports/coverage $(ALLNAMES)

test-loop: 
	py.test
	ptw --ext=.py,.pyx --ignore=doc

flake8: 
	py.test --flake8

pep8: flake8

docstyle: 
	py.test --docstyle

lint: 
	py.test --pylint

dependencies:
	pip install -U numpy
	pip install -Ur requirements.txt

.PHONY: all clean build install install-dev test test-km3modules test-nocov flake8 pep8 dependencies docstyle
