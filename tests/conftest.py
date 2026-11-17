import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from inbrain_dementia_data_analyser.config import load_config  # noqa: E402


@pytest.fixture(scope="session")
def cfg():
    return load_config(ROOT / "configs")


@pytest.fixture
def config_copy(tmp_path):
    """A writable copy of the default configuration folder."""
    dst = tmp_path / "configs"
    shutil.copytree(ROOT / "configs", dst)
    return dst
