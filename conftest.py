# -*- coding: utf-8 -*-
"""Pytest configuration for NOREC4DNA tests."""

import os
import sys
from collections.abc import Callable, Iterator
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).parent.resolve()
SRC_DIR = REPO_ROOT / "src"
TESTS_DIR = REPO_ROOT / "tests"
NOREC4DNA_DIR = REPO_ROOT

# Store original working directory
ORIGINAL_CWD = os.getcwd()


def pytest_configure(config: object) -> None:
    """Called before test collection."""
    if str(SRC_DIR) not in sys.path:
        sys.path.insert(0, str(SRC_DIR))
    if str(NOREC4DNA_DIR) not in sys.path:
        sys.path.insert(1, str(NOREC4DNA_DIR))

    # Change to NOREC4DNA directory so relative paths in tests work
    os.chdir(NOREC4DNA_DIR)


def pytest_unconfigure(config: object) -> None:
    """Called after all tests are done."""
    # Restore original working directory
    os.chdir(ORIGINAL_CWD)


@pytest.fixture(scope="session", autouse=True)
def setup_test_environment() -> Iterator[None]:
    """Set up test environment before all tests."""
    # Ensure we're in the right directory
    os.chdir(NOREC4DNA_DIR)
    yield
    # Cleanup after all tests
    os.chdir(ORIGINAL_CWD)


@pytest.fixture
def tests_dir() -> Path:
    """Return the absolute path to the tests directory."""
    return TESTS_DIR


@pytest.fixture
def norec4dna_dir() -> Path:
    """Return the absolute path to the NOREC4DNA directory."""
    return NOREC4DNA_DIR


@pytest.fixture
def test_file_path() -> Callable[[str], str]:
    """
    Get absolute path to a test file.

    Usage:
        def test_something(test_file_path):
            filepath = test_file_path("logo.jpg")
            # filepath is now the absolute path to NOREC4DNA/tests/logo.jpg
    """

    def _resolve(filename: str) -> str:
        return str(TESTS_DIR / filename)

    return _resolve


@pytest.fixture
def cmp_file_path() -> Callable[[str], str]:
    """
    Get absolute path to a comparison file in tests directory.

    Usage:
        def test_something(cmp_file_path):
            filepath = cmp_file_path("cmp_logo.jpg")
    """

    def _resolve(filename: str) -> str:
        return str(TESTS_DIR / filename)

    return _resolve
