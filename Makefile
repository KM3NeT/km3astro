PKGNAME=km3astro

default: build

all: install

build: 
	python setup.py build_ext --inplace

build-trace:
	python setup.py build_ext --inplace --define CYTHON_TRACE

install: 
	pip install ".[full]"

install-dev: dev-dependencies
	pip install -e ".[full]"

clean:
	python setup.py clean --all
	rm -f $(PKGNAME)/*.cpp
	rm -f $(PKGNAME)/*.c
	rm -f -r build/
	rm -f $(PKGNAME)/*.so

test: build
	py.test --junitxml=./junit.xml || true

test-cov: build
	py.test --junitxml=./junit.xml \
		--cov ./ --cov-report term-missing --cov-report xml || true

test-loop: build
	# pip install -U pytest-watch
	py.test || true
	ptw --ext=.py,.pyx --beforerun "make build"

flake8: 
	py.test --flake8 || true

pep8: flake8

lint: 
	py.test --pylint || true

dependencies:
	pip install -Ur requirements.txt

dev-dependencies:
	pip install -Ur dev-requirements.txt

doc-dependencies:
	pip install -Ur doc-requirements.txt

.PHONY: all clean build install test test-nocov flake8 pep8 dependencies dev-dependencies doc-dependencies
