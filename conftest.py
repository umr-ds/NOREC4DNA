# -*- coding: utf-8 -*-
"""
Pytest configuration for NOREC4DNA tests.

This file ensures tests run correctly regardless of the current working directory
by:
1. Adding the NOREC4DNA directory to sys.path
2. Changing to the NOREC4DNA directory before running tests
3. Providing fixtures for test file paths
"""
import os
import sys
import pytest
from pathlib import Path

# Get the directory containing this conftest.py file (tests directory)
TESTS_DIR = Path(__file__).parent.absolute()

# Get the NOREC4DNA root directory (parent of tests)
NOREC4DNA_DIR = TESTS_DIR.parent.absolute()

# Store original working directory
ORIGINAL_CWD = os.getcwd()


def pytest_configure(config):
    """Called before test collection."""
    # Add NOREC4DNA directory to sys.path for imports
    if str(NOREC4DNA_DIR) not in sys.path:
        sys.path.insert(0, str(NOREC4DNA_DIR))
    
    # Change to NOREC4DNA directory so relative paths in tests work
    os.chdir(NOREC4DNA_DIR)


def pytest_unconfigure(config):
    """Called after all tests are done."""
    # Restore original working directory
    os.chdir(ORIGINAL_CWD)


@pytest.fixture(scope="session", autouse=True)
def setup_test_environment():
    """Set up test environment before all tests."""
    # Ensure we're in the right directory
    os.chdir(NOREC4DNA_DIR)
    yield
    # Cleanup after all tests
    os.chdir(ORIGINAL_CWD)


@pytest.fixture
def tests_dir():
    """Return the absolute path to the tests directory."""
    return TESTS_DIR


@pytest.fixture
def norec4dna_dir():
    """Return the absolute path to the NOREC4DNA directory."""
    return NOREC4DNA_DIR


@pytest.fixture
def test_file_path(filename):
    """
    Get absolute path to a test file.
    
    Usage:
        def test_something(test_file_path):
            filepath = test_file_path("logo.jpg")
            # filepath is now the absolute path to NOREC4DNA/tests/logo.jpg
    """
    return str(TESTS_DIR / filename)


@pytest.fixture
def cmp_file_path(filename):
    """
    Get absolute path to a comparison file in tests directory.
    
    Usage:
        def test_something(cmp_file_path):
            filepath = cmp_file_path("cmp_logo.jpg")
    """
    return str(TESTS_DIR / filename)
