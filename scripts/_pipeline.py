#!/usr/bin/env python3
from __future__ import annotations

import os
from pathlib import Path


DEFAULT_PROJECT_ROOT = Path.home() / "projects" / "gv_project"

STRICT_CONSENSUS_CATEGORY = "rankpct_high_confidence"
RELAXED_CONSENSUS_CATEGORY = "rankpct_moderate_confidence"


def project_root() -> Path:
    return Path(os.environ.get("GV_PROJECT_ROOT", DEFAULT_PROJECT_ROOT))
